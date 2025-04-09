version 1.0

# Copyright (c) 2024 Sequencing Analysis Support Core - Leiden University Medical Center

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "tasks/sequali.wdl" as sequali 
import "tasks/samtools.wdl" as samtools
import "tasks/minimap2.wdl" as minimap2 
import "tasks/clair3.wdl" as clair3 
import "tasks/multiqc.wdl" as multiqc 
import "tasks/pbmm2.wdl" as pbmm2
import "tasks/chunked-scatter.wdl" as chunkedScatter
import "tasks/deepvariant.wdl" as deepvariant
import "tasks/picard.wdl" as picard
import "tasks/modkit.wdl" as modkit
import "tasks/vep.wdl" as vep
import "tasks/mosdepth.wdl" as mosdepth


struct SampleDataset {
    String readgroup_id
    File file
    String? lib_id
}

struct Sample {
    String id
    Array[SampleDataset]+ datasets
    File? clair3modelTar
    String? clair3builtinmodel
    String? deepvariantModelType
}


workflow LongReadVariantCalling {
    input {
        Array[Sample] samples
        File referenceFasta 
        File referenceFastaFai
        File? clair3modelTar
        String? clair3builtinmodel
        String clair3platform = "ont"
        String pbmm2Preset = "HIFI"
        String minimap2preset = "map-ont"
        String outputPrefix = "."
        String deepvariantModelType = "ONT_R104"


        File? vepCacheTar

        Boolean usePbmm2 = false
        Boolean runClair3 = true 
        Boolean runDeepVariant = false 
        Boolean runModKit = false

        Map[String, String] dockerImages = {
            "sequali": "quay.io/biocontainers/sequali:0.12.0--py312hf67a6ed_0",
            # Minimap 2.28 samtools 1.20
            "minimap2": "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0",
            "pbmm2": "quay.io/biocontainers/pbmm2:1.17.0--h9ee0642_0",
            "samtools": "quay.io/biocontainers/samtools:1.21--h96c455f_1",
            "mosdepth": "quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1",
            "clair3": "quay.io/biocontainers/clair3:1.0.11--py39hd649744_0",
            # Version 1.8.0 has a bug.
            # https://github.com/google/deepvariant/issues/912
            "deepvariant": "google/deepvariant:1.6.1",
            "vep": "quay.io/biocontainers/ensembl-vep:113.3--pl5321h2a3209d_0",
            "scatter-regions": "quay.io/biocontainers/chunked-scatter:1.0.0--py_0",
            "picard": "quay.io/biocontainers/picard:3.3.0--hdfd78af_0",
            "modkit": "quay.io/biocontainers/ont-modkit:0.4.3--hcdda2d0_0",
            "multiqc": "quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0",
        }
    }
    
    scatter (sample in samples) {
        String sampleDir = "~{outputPrefix}/~{sample.id}"

        scatter (dataset in sample.datasets) {
            String lib_id = select_first([dataset.lib_id, "lib1"])
            String readgroupID = "~{sample.id}-~{lib_id}-~{dataset.readgroup_id}"
            String libraryID = "~{sample.id}-~{lib_id}"

            call sequali.Sequali as sequaliTask {
                input: 
                    reads = dataset.file,
                    outDir = sampleDir,
                    dockerImage = dockerImages["sequali"],
            }

            String bamPrefix = if length(sample.datasets) == 1 then sample.id else readgroupID
            if (!usePbmm2) {
                call minimap2.Mapping as minimap2Mapping {
                    input:
                        presetOption = minimap2preset,
                        outputPrefix = "~{sampleDir}/~{bamPrefix}",
                        referenceFile = referenceFasta,
                        queryFile = dataset.file,
                        readgroup = "@RG\\tID:~{readgroupID}\\tLB:~{libraryID}\\tSM:~{sample.id}",
                        dockerImage = dockerImages["minimap2"],
                }
            }
            if (usePbmm2) {
                call pbmm2.Mapping as pacBioMapping {
                    input:
                        presetOption = pbmm2Preset,
                        sample = sample.id,
                        outputPrefix = "~{sampleDir}/~{bamPrefix}",
                        referenceMMI = referenceFasta,
                        queryFile = dataset.file,
                        unmapped = true,  # For archive purposes all reads should be in the BAM file.
                        dockerImage = dockerImages["pbmm2"],
                }
            }
            File sampleBamFiles = select_first([minimap2Mapping.bam, pacBioMapping.outputAlignmentFile])
            File sampleBamIndexes = select_first([minimap2Mapping.bamIndex, pacBioMapping.outputIndexFile])
        }

        if (length(sampleBamFiles) > 1) {
            call samtools.Merge as mergeBam {
                input:
                    bamFiles=sampleBamFiles,
                    outputBamPath="~{sampleDir}/~{sample.id}.bam",
                    dockerImage = dockerImages["samtools"]
            }
        }

        File bam = select_first([mergeBam.outputBam, sampleBamFiles[0]])
        File bamIndex = select_first([mergeBam.outputBamIndex, sampleBamIndexes[0]])

        call mosdepth.Mosdepth as mosdepthTask {
            input:
                bam = bam,
                bamIndex = bamIndex,
                prefix = "~{sampleDir}/~{sample.id}.bam",
                noPerBase = true, # Let's not waste time with this. 
                dockerImage = dockerImages["mosdepth"],
        }

        if (runClair3) {
            call clair3.Clair3 as clair3Task {
                input: 
                    outputPrefix = "~{sampleDir}/~{sample.id}.clair3",
                    bam = bam,
                    bamIndex = bamIndex,
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaFai,
                    modelTar = if defined(sample.clair3modelTar) then sample.clair3modelTar else clair3modelTar,
                    builtinModel = if defined(sample.clair3builtinmodel) then sample.clair3builtinmodel else clair3builtinmodel,
                    platform = clair3platform,
                    sampleName = sample.id,
                    dockerImage = dockerImages["clair3"],
            }

            if (defined(vepCacheTar)) {
                call vep.Vep as clair3Vep {
                    input: 
                        inputFile = clair3Task.vcf,
                        outputPath = "~{sampleDir}/~{sample.id}.clair3.vep.vcf.gz",
                        cacheTar = select_first([vepCacheTar]),
                        dockerImage = dockerImages["vep"],
                }
            }
        }

        if (runDeepVariant) {
            call chunkedScatter.ScatterRegions as scatterList {
                input:
                    inputFile = referenceFastaFai,
                    scatterSizeMillions = 100,
                    splitContigs = true,
                    dockerImage = dockerImages["scatter-regions"],
            }

            scatter (region in scatterList.scatters) {
                call deepvariant.RunDeepVariant as deepVariantTask {
                    input:
                        referenceFasta = referenceFasta,
                        referenceFastaIndex = referenceFastaFai,
                        inputBam = bam, 
                        inputBamIndex = bamIndex,
                        modelType = select_first([sample.deepvariantModelType, deepvariantModelType]),
                        outputVcf = "~{sample.id}.~{basename(region)}.vcf.gz",
                        regions = region,
                        dockerImage = dockerImages["deepvariant"],
                }
            }
            Array[File] deepVariantReports = flatten(deepVariantTask.outputVCFStatsReport)

            call picard.MergeVCFs as mergeDeepVariantVCFs {
                input:
                    inputVCFs = deepVariantTask.outputVCF,
                    inputVCFsIndexes = deepVariantTask.outputVCFIndex,
                    outputVcfPath = "~{sampleDir}/~{sample.id}.deepvariant.vcf.gz",
                    dockerImage = dockerImages["picard"],
            }

            if (defined(vepCacheTar)) {
                call vep.Vep as deepVariantVep {
                    input: 
                        inputFile = mergeDeepVariantVCFs.outputVcf,
                        outputPath = "~{sampleDir}/~{sample.id}.deepvariant.vep.vcf.gz",
                        cacheTar = select_first([vepCacheTar]),
                        dockerImage = dockerImages["vep"],
                }
            }
        }

        if (runModKit) {
            call modkit.Pileup as ModKitPileup {
                input: 
                    bam=bam, 
                    bamIndex=bamIndex, 
                    outputBed="~{sampleDir}/~{sample.id}.modkit.bed",
                    referenceFasta=referenceFasta,
                    referenceFastaFai=referenceFastaFai, 
                    logFilePath="~{sampleDir}/~{sample.id}.modkit.log",
                    dockerImage = dockerImages["modkit"],
            }
        }
    }

    call multiqc.MultiQC {
        input:
            reports = flatten([
                flatten(sequaliTask.json), 
                flatten(select_all(deepVariantReports)),
                select_all(clair3Vep.statsHtml),
                select_all(deepVariantVep.statsHtml),
                mosdepthTask.globalDist,
                mosdepthTask.summary,
                select_all(mosdepthTask.perBaseBed),
                select_all(mosdepthTask.regionsBed),
            ]),
            dataDir = false,
            dockerImage = dockerImages["multiqc"]
    }

    output {
        File multiqcReport = MultiQC.multiqcReport 
        Array[File] bamFiles = bam 
        Array[File] bamIndexes = bamIndex 
        Array[File] clair3VcfFiles = select_all(clair3Task.vcf) 
        Array[File] clair3VcfIndexes = select_all(clair3Task.vcfIndex) 
        Array[File] vepAnnotatedFiles = flatten([select_all(clair3Vep.outputFile), 
                                                 select_all(deepVariantVep.outputFile)]) 
        Array[File] deepVariantVcfFiles = select_all(mergeDeepVariantVCFs.outputVcf)
        Array[File] deepVariantVcfIndexes = select_all(mergeDeepVariantVCFs.outputVcfIndex)
        Array[File] sequaliReports = flatten(sequaliTask.html)
        Array[File] modKitBed = select_all(ModKitPileup.out)
        Array[File] modKitBedGraph = flatten(select_all(ModKitPileup.outFiles))
        Array[File] modKitLog = select_all(ModKitPileup.logFile)
        Array[File] vepHtmlReports = flatten([
                select_all(clair3Vep.statsHtml), 
                select_all(deepVariantVep.statsHtml),
            ])
        Array[File] mosdepthSummary = mosdepthTask.summary 
        Array[File] mosdepthGlobalDist = mosdepthTask.globalDist
        Array[File] mosdepthPerBaseBed = select_all(mosdepthTask.perBaseBed)
        Array[File] mosdepthRegionsBed = select_all(mosdepthTask.regionsBed)
    }

    parameter_meta {
        # input 
        dockerImages: {description: "specify which docker images should be used for running this pipeline",
                       category: "advanced" }
        samples: {description: "The samples with metadata and files.", category: "required"}
        referenceFasta: {description: "The reference FASTA file.", category: "required"}
        referenceFastaFai: {description: "The reference FASTA index file.", category: "required"}
        
        clair3modelTar: {description: "TAR file with clair3 model if no builtin model is used", category: "common"}
        clair3builtinmodel: {description: "String describing a builtin model if no TAR file is used.", category: "common"}
        clair3platform: {description: "String describing the clair3 platform", category: "common"}
        minimap2preset: {description: "Minimap2 preset string", category: "common"}
        vepCacheTar: {description: "A TAR file with a VEP cache, when given will cause VEP to run.", category: "common"}
        outputPrefix: {description: "Where to place the data.", category: "advanced"}
        deepvariantModelType: {description: "The DeepVariant model to use", category: "advanced"}

        runClair3: {description: "Whether to run clair3.", category: "common"} 
        runDeepVariant: {description: "Whether to run DeepVariant", category: "common"}
        runModKit: {description: "Whether to run ModKit", category: "common"}
        usePbmm2: {description: "Use pbmm2 instead of minimap2 for mapping.", category: "common"}
        pbmm2Preset: {description: "Pbmm2 preset for mapping reads.", category: "common"}

        # output
        multiqcReport: {description: "The MultiQC report."}
        bamFiles: {description: "The aligned BAM files generated by minimap2."}
        bamIndexes: {description: "The indexes for the aligned BAM files."}
        clair3VcfFiles: {description: "VCF files generated by clair3."}
        clair3VcfIndexes: {description: "VCF indexes for clair3 VCF files."}
        deepVariantVcfFiles: {description: "VCF files generated by DeepVariant."}
        deepVariantVcfIndexes: {description: "VCF indexes for DeepVariant VCF files."}
        sequaliReports: {description: "HTML reports generated by sequali."}
        modKitBed: {description: "BED file generated by modkit."}
        modLitBedGraph: {description: "All files generated by modkit if the --bedgraph option is used."}
        modKitLog: {description: "ModKit log file."}
        vepAnnotatedFiles: {description: "VCF file annotated by VEP."}
        vepHtmlReports: {description: "The VEP HTML reports."}
        modKitBedGraph: {description: "BedGraph output files for ModKit."}
        mosdepthSummary: {description: "Mosdepth summary files."}
        mosdepthGlobalDist: {description: "Mosdepth global dist files."}
        mosdepthPerBaseBed: {description: "Mosdepth per base bed files."}
        mosdepthRegionsBed: {description: "Mosdepth region bed files."}
    }
}

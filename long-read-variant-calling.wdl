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
import "tasks/chunked-scatter.wdl" as chunkedScatter
import "tasks/deepvariant.wdl" as deepvariant
import "tasks/picard.wdl" as picard
import "tasks/modkit.wdl" as modkit

task fileIsFastx {
    input {
        File file
    }
    command <<<
    python <<CODE
    import gzip
    with open("~{file}", 'rb') as f:
        begin = f.peek(1024)
        if begin[:2] == b"\x1f\x8b":
            with gzip.open(f, "rb") as gz: 
                begin = gz.read(1024)
    if begin[0] == "@" or begin[0] == ">":
        print("true")
        sys.exit(0)
    print("false")
    CODE
    >>>
    output {
        Boolean result = read_boolean(stdout())
    }

    runtime {
        # python:3.7-slim's sha256 digest. This image is based on debian buster.
        docker: "python@sha256:e0f6a4df17d5707637fa3557ab266f44dddc46ebfc82b0f1dbe725103961da4e"
    }
    
}

struct SampleDataset {
    String readgroup_id
    File file
    String? lib_id
}

struct Sample {
    String id
    Array[SampleDataset]+ datasets
}


workflow LongReadVariantCalling {
    input {
        Array[Sample] samples
        File referenceFasta 
        File referenceFastaFai
        File? clair3modelTar
        String? clair3builtinmodel
        String clair3platform
        String minimap2preset   
        String outputPrefix = "."
        String deepvariantModelType = "ONT_R104"

        Boolean runClair3 = true 
        Boolean runDeepVariant = false 
        Boolean runModKit = false
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
            }

            String bamPrefix = if length(sample.datasets) == 1 then sample.id else readgroupID
            call minimap2.Mapping as minimap2Mapping {
                input:
                    presetOption = minimap2preset,
                    outputPrefix = "~{sampleDir}/~{bamPrefix}",
                    referenceFile = referenceFasta,
                    queryFile = dataset.file,
                    readgroup = "@RG\\tID:~{readgroupID}\\tLB:~{libraryID}\\tSM:~{sample.id}",
            }
        }
        if (length(minimap2Mapping.bam) > 1) {
            call samtools.Merge as mergeBam {
                input:
                    bamFiles=minimap2Mapping.bam,
                    outputBamPath="~{sampleDir}/~{sample.id}.bam",
            }
        }

        File bam = select_first([mergeBam.outputBam, minimap2Mapping.bam[0]])
        File bamIndex = select_first([mergeBam.outputBamIndex, minimap2Mapping.bamIndex[0]])

        if (runClair3) {
            call clair3.Clair3 as clair3Task {
                input: 
                    outputPrefix = "~{sampleDir}/~{sample.id}.clair3",
                    bam = bam,
                    bamIndex = bamIndex,
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaFai,
                    modelTar = clair3modelTar,
                    builtinModel = clair3builtinmodel,
                    platform = clair3platform,
                    sampleName = sample.id,
            }
        }

        if (runDeepVariant) {
            call chunkedScatter.ScatterRegions as scatterList {
                input:
                    inputFile = referenceFastaFai,
                    scatterSizeMillions = 100,
                    splitContigs = true,
            }

            scatter (region in scatterList.scatters) {
                call deepvariant.RunDeepVariant as deepVariantTask {
                    input:
                        referenceFasta = referenceFasta,
                        referenceFastaIndex = referenceFastaFai,
                        inputBam = bam, 
                        inputBamIndex = bamIndex,
                        modelType = deepvariantModelType,
                        outputVcf = "~{sample.id}.~{basename(region)}.vcf.gz",
                        regions = region,
                }
            }
            Array[File] deepVariantReports = flatten(deepVariantTask.outputVCFStatsReport)

            call picard.MergeVCFs as mergeDeepVariantVCFs {
                input:
                    inputVCFs = deepVariantTask.outputVCF,
                    inputVCFsIndexes = deepVariantTask.outputVCFIndex,
                    outputVcfPath = "~{sampleDir}/~{sample.id}.deepvariant.vcf.gz",
            }

        }

        if (runModKit) {
            call modkit.Pileup as ModKitPileup {
                input: 
                    bam=bam, 
                    bamIndex=bamIndex, 
                    outputBed="~{sampleDir}/~{sample.id}.pseudoepic.bed",
                    referenceFasta=referenceFasta,
                    referenceFastaFai=referenceFastaFai, 
                    logFilePath="~{sampleDir}/~{sample.id}.modkit.log",
            }

        }

    }

    call multiqc.MultiQC {
        input:
            reports = flatten([
                flatten(sequaliTask.json), 
                flatten(select_all(deepVariantReports)),
            ]),
            dataDir = false,
    }

    output {
        File multiqcReport = MultiQC.multiqcReport 
        Array[File] bamFiles = bam 
        Array[File] bamIndexes = bamIndex 
        Array[File] clair3VcfFiles = select_all(clair3Task.vcf) 
        Array[File] clair3VcfIndexes = select_all(clair3Task.vcfIndex) 
        Array[File] deepVariantVcfFiles = select_all(mergeDeepVariantVCFs.outputVcf)
        Array[File] deepVariantVcfIndexes = select_all(mergeDeepVariantVCFs.outputVcfIndex)
        Array[File] sequaliReports = flatten(sequaliTask.html)
        Array[File] modKitBed = select_all(ModKitPileup.out) 
        Array[File] modKitLog = select_all(ModKitPileup.logFile)
    }
}

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
import "tasks/minimap2.wdl" as minimap2 
import "tasks/clair3.wdl" as clair3 
import "tasks/multiqc.wdl" as multiqc 

task fileIsFastx {
    input {
        File file
    }
    command <<<
    python <<CODE
    with open("~{file}", 'rb') as f: 
        begin = f.read(1024)
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

task BamToFastq {
    # A simpler task than in biowdl/tasks for this particular use case
    input {
        File inputBam 
        Int timeMinutes = 1 + ceil(size(inputBam, "G") * 2)
        String prefix = "sample"
        String dockerImage = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }

    command {
        samtools reset -u ~{inputBam} | samtools fastq | bgzip -l 1 > ~{prefix}.fastq.gz
    }
    output {
        File fastq = "~{prefix}.fastq.gz"
    }
    runtime {
        cpu: 2  # One for decompressing one for compressing
        memory: "2GiB"
        docker: dockerImage
        time_minutes: timeMinutes
    }

}
struct Sample {
    String id 
    File reads
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
    }
    
    scatter (sample in samples) {
        call sequali.Sequali as sequaliTask {
            input: 
                reads = sample.reads,
                outDir = outputPrefix + "/sequali/",
        }
        call fileIsFastx {
            input:
                file=sample.reads,
        }
        if (!fileIsFastx.result) {
            call BamToFastq {
                input:
                    inputBam = sample.reads,
                    prefix = sample.id,
            }
        }
        File reads = select_first([BamToFastq.fastq, sample.reads])
        call minimap2.Mapping as minimap2Mapping {
            input:
                presetOption = minimap2preset,
                outputPrefix = "~{outputPrefix}/bam/~{sample.id}",
                referenceFile = referenceFasta,
                queryFile = reads,
        }

        call clair3.Clair3 as clair3Task {
            input: 
                outputPrefix = "~{outputPrefix}/clair3/~{sample.id}",
                bam = minimap2Mapping.bam,
                bamIndex = minimap2Mapping.bamIndex,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                modelTar = clair3modelTar,
                builtinModel = clair3builtinmodel,
                platform = clair3platform,
        }
    }
    call multiqc.MultiQC {
        input:
            reports = sequaliTask.json,
            dataDir = false,
    }

    output {
        File multiqcReport = MultiQC.multiqcReport 
        Array[File] bamFiles = minimap2Mapping.bam 
        Array[File] bamIndexes = minimap2Mapping.bamIndex 
        Array[File] vcfFiles = clair3Task.vcf 
        Array[File] vcfIndexes = clair3Task.vcfIndex 
        Array[File] sequaliReports = sequaliTask.html
    }
}
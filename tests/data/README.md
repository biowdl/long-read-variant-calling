# Test file generation

The read IDs in [GM24385_1_mitochondrial_reads.txt](./GM24385_1_mitochondrial_reads.txt)
are from the GM@4385_1.fastq.gz file that can be found in the 
[GIAB data index](https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/sequence.index.AJtrio_UCSC_ONT_UL_Promethion_03312019.HG002).

Test data was generated with the following command:

```
samtools view -N GM24385_1_mitochondrial_reads.txt GM24385_1.fastq.gz | samtools fastq | gzip -9 -c > GM24385_1_mitochondrial_reads.fastq.gz
```

The mitochrondrial genome can be downloaded from:
https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=fasta 


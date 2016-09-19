# Scripts
TRIM_SRC=scripts/run_cutadapt.py
TRIM_EXE=python $(TRIM_SRC)

ASSEMBLY_SRC=scripts/assemble.py
ASSEMBLY_EXE=python $(ASSEMBLY_SRC)

#Variables
RAW_READS=$(wildcard test_fastq/fastq/*_001.fastq.gz)
RAW_READ1=$(wildcard test_fastq/fastq/*_R1_001.fastq.gz)
RAW_READ2=$(wildcard test_fastq/fastq/*_R2_001.fastq.gz)
RAW_QUALS=$(patsubst test_fastq/fastq/%.fastq.gz, test_fastq/QC/%_fastqc.html, $(RAW_READS))

TRIMMED_READS=$(patsubst test_fastq/fastq/%_001.fastq.gz, trimmed_reads/fastq/%_trimmed.fastq.gz, $(RAW_READS))
TRIMMED_QUALS=$(patsubst trimmed_reads/fastq/%.fastq.gz, trimmed_reads/QC/%_fastqc.html, $(TRIMMED_READS))

SPADES_ASSEMBLIES=$(patsubst test_fastq/fastq/%_R1_001.fastq.gz, assembly/spades_%/scaffolds.fasta, $(RAW_READ1))
IDBA_UD_ASSEMBLIES=$(patsubst test_fastq/fastq/%_R1_001.fastq.gz, assembly/idba_ud_%/scaffold.fa, $(RAW_READ1))
MEGAHIT_ASSEMBLIES=$(patsubst test_fastq/fastq/%_R1_001.fastq.gz, assembly/megahit_%/final.contigs.fa, $(RAW_READ1))
#Directories
RAW_DIR=raw_sequence
RAW_READ_DIR=$(RAW_DIR)/fastq
RAW_QUAL_DIR=$(RAW_DIR)/QC
TRIM_DIR=trimmed_reads
TRIM_READ_DIR=$(TRIM_DIR)/fastq
TRIM_QUAL_DIR=$(TRIM_DIR)/QC

ALIGN_DIR=alignment
FULL_IMS101_ALIGN_DIR=$(ALIGN_DIR)/full_dataset_v_IMS101


# Scripts
TRIM_SRC=scripts/run_cutadapt.py
TRIM_EXE=python $(TRIM_SRC)

ASSEMBLY_SRC=scripts/assemble.py
ASSEMBLY_EXE=python $(ASSEMBLY_SRC)

ALIGN_SRC=scripts/run_bowtie2.py
ALIGN_EXE=python $(ALIGN_SRC)

QUALIMAP_EXE=~/Downloads/qualimap_v2.2/qualimap


#Variables
RAW_READS=$(wildcard $(RAW_READ_DIR)/*_001.fastq.gz)
RAW_READ1=$(wildcard $(RAW_READ_DIR)/*_R1_001.fastq.gz)
RAW_READ2=$(wildcard $(RAW_READ_DIR)/*_R2_001.fastq.gz)
RAW_QUALS=$(patsubst $(RAW_READ_DIR)/%.fastq.gz, $(RAW_QUAL_DIR)/%_fastqc.html, $(RAW_READS))

TRIMMED_READS=$(patsubst $(RAW_READ_DIR)/%_001.fastq.gz, $(TRIM_READ_DIR)/%_trimmed.fastq.gz, $(RAW_READS))
TRIMMED_QUALS=$(patsubst $(TRIM_READ_DIR)/%.fastq.gz, $(TRIM_QUAL_DIR)/%_fastqc.html, $(TRIMMED_READS))

SPADES_ASSEMBLIES=$(patsubst $(RAW_READ_DIR)/%_R1_001.fastq.gz, assembly/spades_%/scaffolds.fasta, $(RAW_READ1))
IDBA_UD_ASSEMBLIES=$(patsubst $(RAW_READ_DIR)/%_R1_001.fastq.gz, assembly/idba_ud_%/scaffold.fa, $(RAW_READ1))
MEGAHIT_ASSEMBLIES=$(patsubst $(RAW_READ_DIR)/%_R1_001.fastq.gz, assembly/megahit_%/final.contigs.fa, $(RAW_READ1))
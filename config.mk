#Directories
RAW_DIR=raw_sequence
RAW_READ_DIR=$(RAW_DIR)/fastq
RAW_QUAL_DIR=$(RAW_DIR)/QC
TRIM_DIR=trimmed_reads
TRIM_READ_DIR=$(TRIM_DIR)/fastq
TRIM_QUAL_DIR=$(TRIM_DIR)/QC
FILTER_DIR=filtered_reads
FILTER_READ_DIR=$(FILTER_DIR)/fastq
FILTER_QUAL_DIR=$(FILTER_DIR)/QC

ALIGN_DIR=alignment
FULL_IMS101_ALIGN_DIR=$(ALIGN_DIR)/full_dataset_v_IMS101

REDUCTION_DIR=assembly_reduction
REDUNDANS_DIR=$(REDUCTION_DIR)/redundans

CDHIT_DIR=$(REDUCTION_DIR)/cdhit

# Scripts
TRIM_SRC=scripts/run_cutadapt.py
TRIM_EXE=python $(TRIM_SRC)

FILTER_SRC=scripts/filter_tagged_fastq.py
FILTER_EXE=python $(FILTER_SRC)

ASSEMBLY_SRC=scripts/assemble.py
ASSEMBLY_EXE=python $(ASSEMBLY_SRC)

ALIGN_SRC=scripts/run_bowtie2.py
ALIGN_EXE=python $(ALIGN_SRC)

QUALIMAP_EXE=~/Downloads/qualimap_v2.2/qualimap

QUAST_EXE=~/Downloads/quast-4.0/quast.py

CDHIT_EXE=~/Downloads/cd-hit-v4.6.6-2016-0711/cd-hit-est
CDHIT_CLUSTER_SRC=scripts/contig_clustering.py
CDHIT_CLUSTER_EXE=python $(CDHIT_CLUSTER_SRC)

VARSCAN_SRC=~/Downloads/VarScan.v2.3.9.jar
VARSCAN_EXE=java -jar $(VARSCAN_SRC)

#Variables
RAW_READS=$(wildcard $(RAW_READ_DIR)/*_001.fastq.gz)
RAW_READ1=$(wildcard $(RAW_READ_DIR)/*_R1_001.fastq.gz)
RAW_READ2=$(wildcard $(RAW_READ_DIR)/*_R2_001.fastq.gz)
RAW_QUALS=$(patsubst $(RAW_READ_DIR)/%.fastq.gz, $(RAW_QUAL_DIR)/%_fastqc.html, $(RAW_READS))

TRIMMED_READS=$(patsubst $(RAW_READ_DIR)/%_001.fastq.gz, $(TRIM_READ_DIR)/%_trimmed.fastq.gz, $(RAW_READS))
TRIMMED_QUALS=$(patsubst $(TRIM_READ_DIR)/%.fastq.gz, $(TRIM_QUAL_DIR)/%_fastqc.html, $(TRIMMED_READS))

FILTERED_READS=$(patsubst $(TRIM_READ_DIR)/%.fastq.gz, $(FILTER_READ_DIR)/%.filtered.fastq.gz, $(TRIMMED_READS))
FILTERED_QUALS=$(patsubst $(FILTER_READ_DIR)/%.fastq.gz, $(FILTER_QUAL_DIR)/%_fastqc.html, $(FILTERED_READS))

SPADES_ASSEMBLIES=$(patsubst $(RAW_READ_DIR)/%_R1_001.fastq.gz, assembly/spades_%/scaffolds.fasta, $(RAW_READ1))
IDBA_UD_ASSEMBLIES=$(patsubst $(RAW_READ_DIR)/%_R1_001.fastq.gz, assembly/idba_ud_%/scaffold.fa, $(RAW_READ1))
MEGAHIT_ASSEMBLIES=$(patsubst $(RAW_READ_DIR)/%_R1_001.fastq.gz, assembly/megahit_%/final.contigs.fa, $(RAW_READ1))
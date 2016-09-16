# Scripts
TRIM_SRC=scripts/run_cutadapt.py
TRIM_EXE=python $(TRIM_SRC)

#Variables
RAW_READS=$(wildcard testfastq/*_001.fastq.gz)
TRIMMED_READS=$(patsubst testfastq/%_001.fastq.gz, %_trimmed.fastq.gz, $(RAW_READS))
include config.mk

RAW_READ_DIR=test_fastq/fastq
TRIM_READ_DIR=trimmed_reads/fastq


.PHONY : all
all : quality_check trim quality_check assemble


# Check fastq quality metrics with fastqc
.PHONY : quality_check
quality_check : $(RAW_QUALS) $(TRIMMED_QUALS)

test_fastq/QC/%_fastqc.html : $(RAW_READ_DIR)/%_001.fastq.gz
	mkdir -p test_fastq/QC
	fastqc -o test_fastq/QC $(RAW_READ_DIR)/$*_001.fastq.gz

trimmed_reads/QC/%_fastqc.html : $(TRIM_READ_DIR)/%_trimmed.fastq.gz
	mkdir -p trimmed_reads/QC
	fastqc -o trimmed_reads/QC $(TRIM_READ_DIR)/$*_trimmed.fastq.gz



# Trim raw reads using run_cutadapt.py
.PHONY : trim
trim : $(TRIMMED_READS)

$(TRIM_READ_DIR)/%_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/%_R2_trimmed.fastq.gz trimmed_reads/%_trimlog.txt : $(RAW_READ_DIR)/%_R1_001.fastq.gz $(RAW_READ_DIR)/%_R2_001.fastq.gz $(TRIM_SRC)
	mkdir -p trimmed_reads
	mkdir -p $(TRIM_READ_DIR)
	$(TRIM_EXE) $(RAW_READ_DIR)/$*_R1_001.fastq.gz $(RAW_READ_DIR)/$*_R2_001.fastq.gz $(TRIM_READ_DIR)/$*_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/$*_R2_trimmed.fastq.gz trimmed_reads/$*_trimlog.txt



## Alignments
# Align reads vs IMS101

# Align reads vs T. theibautii H9-4

# Calculate alignment statistics



## Assembly of full dataset
# Assemble
.PHONY : assemble
assemble : $(SPADES_ASSEMBLIES) $(IDBA_UD_ASSEMBLIES) $(MEGAHIT_ASSEMBLIES)

assembly/spades_%/scaffolds.fasta : $(TRIM_READ_DIR)/%_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/%_R2_trimmed.fastq.gz $(ASSEMBLY_SRC)
	mkdir -p assembly
	$(ASSEMBLY_EXE) -m spades -f $(TRIM_READ_DIR)/$*_R1_trimmed.fastq.gz -r $(TRIM_READ_DIR)/$*_R2_trimmed.fastq.gz -o assembly -p $*

assembly/megahit_%/final.contigs.fa : $(TRIM_READ_DIR)/%_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/%_R2_trimmed.fastq.gz $(ASSEMBLY_SRC)
	mkdir -p assembly
	$(ASSEMBLY_EXE) -m megahit -f $(TRIM_READ_DIR)/$*_R1_trimmed.fastq.gz -r $(TRIM_READ_DIR)/$*_R2_trimmed.fastq.gz -o assembly -p $*

assembly/idba_ud_%/scaffold.fa : $(TRIM_READ_DIR)/%_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/%_R2_trimmed.fastq.gz $(ASSEMBLY_SRC)
	mkdir -p assembly
	$(ASSEMBLY_EXE) -m idba_ud -f $(TRIM_READ_DIR)/$*_R1_trimmed.fastq.gz -r $(TRIM_READ_DIR)/$*_R2_trimmed.fastq.gz -o assembly -p $*


# Calculate stats


## Extract Tricho reads by kmer-based binning
# Kmer counting

# Generate plots

# Do extraction


## Assembly of reduced dataset
# Assemble

# Calculate statistics


## Redundancy reduction


## Annotation


## Read alignment


## Variant detection


## Conservation statistics

.PHONY : clean
clean : 
	rm -rf $(TRIM_READ_DIR)
	
.PHONY : variables
variables:
	@echo RAW_READS: $(RAW_READS)
	@echo RAW_READ1: $(RAW_READ1)
	@echo RAW_READ2: $(RAW_READ2)
	@echo TRIMMED_READS: $(TRIMMED_READS)
	@echo SPADES_ASSEMBLIES: $(SPADES_ASSEMBLIES)
	@echo IDBA_UD_ASSEMBLIES: $(IDBA_UD_ASSEMBLIES)
	@echo MEGAHIT_ASSEMBLIES: $(MEGAHIT_ASSEMBLIES)
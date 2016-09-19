include config.mk


.PHONY : all
all : $(TRIMMED_READS)



# Trim raw reads using run_cutadapt.py
.PHONY : trim
trim : $(TRIMMED_READS)

%_R1_trimmed.fastq.gz %_R2_trimmed.fastq.gz %_trimlog.txt : test_fastq/fastq/%_R1_001.fastq.gz test_fastq/fastq/%_R2_001.fastq.gz $(TRIM_SRC)
	$(TRIM_EXE) test_fastq/fastq/$*_R1_001.fastq.gz testfastq/fastq/$*_R2_001.fastq.gz $*_R1_trimmed.fastq.gz $*_R2_trimmed.fastq.gz $*_trimlog.txt



## Alignments
# Align reads vs IMS101

# Align reads vs T. theibautii H9-4

# Calculate alignment statistics



## Assembly of full dataset
# Assemble
.PHONY : assemble
assemble : assembly/spades_test_reads/scaffolds.fasta assembly/idba_ud_*/scaffold.fa assembly/megahit_*/final.contigs.fa

assembly/spades_%/scaffolds.fasta : test_fastq/fastq/%_R1_001.fastq.gz test_fastq/fastq/%_R2_001.fastq.gz $(ASSEMBLY_SRC)
	mkdir -p assembly
	$(ASSEMBLY_EXE) -m spades -f test_fastq/fastq/$*_R1_001.fastq.gz -r test_fastq/fastq/$*_R2_001.fastq.gz -o assembly -p $*


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
	rm -f *_trimmed.fastq.gz
	
.PHONY : variables
variables:
	@echo RAW_READS: $(RAW_READS)
	@echo RAW_READ1: $(RAW_READ1)
	@echo RAW_READ2: $(RAW_READ2)
	@echo TRIMMED_READS: $(TRIMMED_READS)
	@echo SPADES_ASSEMBLIES: $(SPADES_ASSEMBLIES)
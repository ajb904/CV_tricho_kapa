include config.mk


.PHONY : all
all : $(TRIMMED_READS)



# Trim raw reads using run_cutadapt.py
.PHONY : trim
trim : $(TRIMMED_READS)

%_R1_trimmed.fastq.gz %_R2_trimmed.fastq.gz %_trimlog.txt : testfastq/%_R1_001.fastq.gz testfastq/%_R2_001.fastq.gz $(TRIM_SRC)
	$(TRIM_EXE) testfastq/$*_R1_001.fastq.gz testfastq/$*_R2_001.fastq.gz $*_R1_trimmed.fastq.gz $*_R2_trimmed.fastq.gz $*_trimlog.txt



## Alignments
# Align reads vs IMS101

# Align reads vs T. theibautii H9-4

# Calculate alignment statistics



## Assembly of full dataset
# Assemble

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
	@echo TRIMMED_READS: $(TRIMMED_READS)
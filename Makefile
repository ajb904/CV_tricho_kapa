include config.mk

.PHONY : all
all : quality_check_raw trim quality_check_trimmed


# Check fastq quality metrics with fastqc
.PHONY : quality_check_raw
quality_check_raw : $(RAW_QUALS)

$(RAW_QUAL_DIR)/%_fastqc.html : $(RAW_READ_DIR)/%.fastq.gz
	mkdir -p $(RAW_QUAL_DIR)
	fastqc -o $(RAW_QUAL_DIR) $<


.PHONY : quality_check_trimmed
quality_check_trimmed : $(TRIMMED_QUALS)

$(TRIM_QUAL_DIR)/%_fastqc.html : $(TRIM_READ_DIR)/%.fastq.gz
	mkdir -p $(TRIM_QUAL_DIR)
	fastqc -o $(TRIM_QUAL_DIR) $<



# Trim raw reads using run_cutadapt.py
.PHONY : trim
trim : $(TRIMMED_READS)

$(TRIM_READ_DIR)/%_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/%_R2_trimmed.fastq.gz trimmed_reads/%_trimlog.txt : $(RAW_READ_DIR)/%_R1_001.fastq.gz $(RAW_READ_DIR)/%_R2_001.fastq.gz $(TRIM_SRC)
	mkdir -p $(TRIM_DIR)
	mkdir -p $(TRIM_READ_DIR)
	$(TRIM_EXE) $(RAW_READ_DIR)/$*_R1_001.fastq.gz $(RAW_READ_DIR)/$*_R2_001.fastq.gz $(TRIM_READ_DIR)/$*_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/$*_R2_trimmed.fastq.gz $(TRIM_DIR)/$*_trimlog.txt



## Alignments

# Align reads vs IMS101
.PHONY : IMS101_alignment
IMS101_alignment : ref/IMS101.fasta $(FULL_IMS101_ALIGN_DIR)/aligned/Tn004_S1_L001_v_IMS101.bam  $(FULL_IMS101_ALIGN_DIR)/aligned/Tn019_S2_L001_v_IMS101.bam $(FULL_IMS101_ALIGN_DIR)/aligned/Tn004_S1_L001/qualimapReport.html $(FULL_IMS101_ALIGN_DIR)/aligned/Tn019_S2_L001/qualimapReport.html 

ref/IMS101.fasta : 
	mkdir -p ref
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000014265.1_ASM1426v1/GCF_000014265.1_ASM1426v1_genomic.fna.gz
	mv GCF_000014265.1_ASM1426v1_genomic.fna.gz ref/IMS101.fasta.gz
	gunzip ref/IMS101.fasta.gz

$(FULL_IMS101_ALIGN_DIR)/aligned/%_v_IMS101.bam : $(TRIM_READ_DIR)/%_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/%_R2_trimmed.fastq.gz $(ALIGN_SRC)
	mkdir -p $(ALIGN_DIR)
	mkdir -p $(FULL_IMS101_ALIGN_DIR)
	$(ALIGN_EXE) -f $(TRIM_READ_DIR)/$*_R1_trimmed.fastq.gz -r $(TRIM_READ_DIR)/$*_R2_trimmed.fastq.gz -i ref/IMS101 -a $(FULL_IMS101_ALIGN_DIR)/aligned -u $(FULL_IMS101_ALIGN_DIR)/unaligned
	
$(FULL_IMS101_ALIGN_DIR)/aligned/%/qualimapReport.html : $(FULL_IMS101_ALIGN_DIR)/aligned/%_v_IMS101.bam
	$(QUALIMAP_EXE) bamqc -bam $(FULL_IMS101_ALIGN_DIR)/aligned/$*_v_IMS101.bam -outdir $(FULL_IMS101_ALIGN_DIR)/aligned/$*_qualimap
	

# Align reads vs T. theibautii H9-4

# Calculate alignment statistics



## Assembly of full dataset
# Assemble
.PHONY : assemble
#assemble : $(SPADES_ASSEMBLIES) $(IDBA_UD_ASSEMBLIES) $(MEGAHIT_ASSEMBLIES)
assemble: $(MEGAHIT_ASSEMBLIES)

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
	rm -rf $(TRIM_DIR)
	rm -rf assembly
	rm -rf $(RAW_QUAL_DIR)


.PHONY : variables
variables:
	@echo RAW_READS: $(RAW_READS)
	@echo RAW_READ1: $(RAW_READ1)
	@echo RAW_READ2: $(RAW_READ2)
	@echo RAW_QUALS: $(RAW_QUALS)
	@echo TRIMMED_READS: $(TRIMMED_READS)
	@echo TRIMMED_QUALS: $(TRIMMED_QUALS)
	@echo SPADES_ASSEMBLIES: $(SPADES_ASSEMBLIES)
	@echo IDBA_UD_ASSEMBLIES: $(IDBA_UD_ASSEMBLIES)
	@echo MEGAHIT_ASSEMBLIES: $(MEGAHIT_ASSEMBLIES)
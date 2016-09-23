include config.mk

.PHONY : all
all : quality_check_raw trim quality_check_trimmed IMS101_alignment


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





.PHONY : filter_reads
filter_reads : $(FILTERED_READS) $(FILTERED_QUALS) ref/Pfluor.1.bt2
## There seems to be contamination (~25%) with P. fluorescens in both samples (this is not present in the D. tertiolecta samples that were processed at the same time). Remove contamination by alignment against P. fluorescens genome (strain LBUM636) and keep the read pairs that don't align. Continue with plan from here, but bear contamination in mind when drawing conclusions...

$(FILTER_READ_DIR)/%_R1_trimmed.filtered.fastq.gz $(FILTER_READ_DIR)/%_R2_trimmed.filtered.fastq.gz : $(TRIM_READ_DIR)/%_R1_trimmed.fastq.gz $(TRIM_READ_DIR)/%_R2_trimmed.fastq.gz ~/bin/fastq_screen.conf ref/Pfluor.1.bt2 $(FILTER_SRC)
	mkdir -p $(FILTER_DIR)
	mkdir -p $(FILTER_READ_DIR)
	$(FILTER_EXE) -1 $(TRIM_READ_DIR)/$*_R1_trimmed.fastq.gz -2 $(TRIM_READ_DIR)/$*_R2_trimmed.fastq.gz -o $(FILTER_READ_DIR)

$(FILTER_QUAL_DIR)/%_fastqc.html : $(FILTER_READ_DIR)/%.fastq.gz
	mkdir -p $(FILTER_QUAL_DIR)
	fastqc -o $(FILTER_QUAL_DIR) $<

ref/Pfluor.fasta : 
	mkdir -p ref
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_001612705.1_ASM161270v1/GCF_001612705.1_ASM161270v1_genomic.fna.gz
	mv GCF_001612705.1_ASM161270v1_genomic.fna.gz ref/Pfluor.fasta.gz
	gunzip ref/Pfluor.fasta.gz
	
ref/%.1.bt2 : ref/%.fasta
	bowtie2-build $< ref/$*
	

## Alignments

# Align reads vs IMS101
.PHONY : IMS101_alignment
IMS101_alignment : ref/IMS101.1.bt2 $(FULL_IMS101_ALIGN_DIR)/aligned/Tn004_S1_L001_v_IMS101.bam  $(FULL_IMS101_ALIGN_DIR)/aligned/Tn019_S2_L001_v_IMS101.bam $(FULL_IMS101_ALIGN_DIR)/aligned/Tn004_S1_L001/qualimapReport.html $(FULL_IMS101_ALIGN_DIR)/aligned/Tn019_S2_L001/qualimapReport.html 

ref/IMS101.fasta : 
	mkdir -p ref
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000014265.1_ASM1426v1/GCF_000014265.1_ASM1426v1_genomic.fna.gz
	mv GCF_000014265.1_ASM1426v1_genomic.fna.gz ref/IMS101.fasta.gz
	gunzip ref/IMS101.fasta.gz

$(FULL_IMS101_ALIGN_DIR)/aligned/%_v_IMS101.bam : $(FILTER_READ_DIR)/%_R1_trimmed.filtered.fastq.gz $(FILTER_READ_DIR)/%_R2_trimmed.filtered.fastq.gz $(ALIGN_SRC) ref/IMS101.1.bt2
	mkdir -p $(ALIGN_DIR)
	mkdir -p $(FULL_IMS101_ALIGN_DIR)
	$(ALIGN_EXE) -f $(FILTER_READ_DIR)/$*_R1_trimmed.filtered.fastq.gz -r $(FILTER_READ_DIR)/$*_R2_trimmed.filtered.fastq.gz -i ref/IMS101 -a $(FULL_IMS101_ALIGN_DIR)/aligned -u $(FULL_IMS101_ALIGN_DIR)/unaligned 2> $(FULL_IMS101_ALIGN_DIR)/aligned/$*_v_IMS101_align.log
	
$(FULL_IMS101_ALIGN_DIR)/aligned/%/qualimapReport.html : $(FULL_IMS101_ALIGN_DIR)/aligned/%_v_IMS101.bam
	$(QUALIMAP_EXE) bamqc -bam $(FULL_IMS101_ALIGN_DIR)/aligned/$*_v_IMS101.bam -outdir $(FULL_IMS101_ALIGN_DIR)/aligned/$*_qualimap
	

# Align reads vs T. theibautii H9-4

# Calculate alignment statistics







# 23 Sept 2016                                                                                                                                                        
# Plan: try a couple of methods to reduce the redundancy of the Tricho assembly made previously (in read_composition/filtered_assembly/0**_spades/) to see if we can get a more contiguous assembly.                                                                                                                                       
# First, try Redundans, using the new sequencing (filtered to remove most of the contaminant) as the paired-end read dataset (this will be more useful than using the original Nextera data as the insert size is larger and more consistent).

.PHONY : redundans
<<<<<<< HEAD
redundans : $(REDUCTION_DIR)/Tn004_S1_L001/redundans/scaffolds.filled.fa $(REDUCTION_DIR)/Tn019_S2_L001/redundans/scaffolds.filled.fa $(REDUCTION_DIR)/Tn004_S1_L001/quast_comparison/report.html $(REDUCTION_DIR)/Tn019_S2_L001/quast_comparison/report.html

$(REDUCTION_DIR)/%/redundans/scaffolds.filled.fa : $(REDUCTION_DIR)/%/spades_contigs.fasta $(FILTER_READ_DIR)/%_R1_trimmed.filtered.fastq.gz $(FILTER_READ_DIR)/%_R2_trimmed.filtered.fastq.gz
	redundans.py -v -t 8 -i $(FILTER_READ_DIR)/$*_R1_trimmed.filtered.fastq.gz $(FILTER_READ_DIR)/$*_R2_trimmed.filtered.fastq.gz -f $(REDUCTION_DIR)/$*/spades_contigs.fasta -o $(REDUCTION_DIR)/$*/redundans
=======
redundans : $(REDUNDANS_DIR)/Tn004_S1_L001/scaffolds.filled.fa $(REDUNDANS_DIR)/Tn019_S2_L001/scaffolds.filled.fa $(REDUNDANS_DIR)/Tn004_S1_L001/quast_comparison/report.html $(REDUNDANS_DIR)/Tn019_S2_L001/quast_comparison/report.html

$(REDUNDANS_DIR)/%/scaffolds.filled.fa : $(REDUNDANS_DIR)/%_spades_contigs.fasta $(FILTER_READ_DIR)/%_R1_trimmed.filtered.fastq.gz $(FILTER_READ_DIR)/%_R2_trimmed.filtered.fastq.gz
	redundans.py -v -t 8 -i $(FILTER_READ_DIR)/$*_R1_trimmed.filtered.fastq.gz $(FILTER_READ_DIR)/$*_R2_trimmed.filtered.fastq.gz -f $(REDUNDANS_DIR)/$*_spades_contigs.fasta -o $(REDUNDANS_DIR)/$*
>>>>>>> align_to_assembly

$(REDUNDANS_DIR)/Tn004_S1_L001_spades_contigs.fasta : ../CV_samples_nextera/read_composition/filtered_assembly/004_spades/contigs.fasta
	mkdir -p $(REDUCTION_DIR)
	mkdir -p $(REDUNDANS_DIR)
	cp $< $@

$(REDUNDANS_DIR)/Tn019_S2_L001_spades_contigs.fasta : ../CV_samples_nextera/read_composition/filtered_assembly/019_spades/contigs.fasta
	mkdir -p $(REDUCTION_DIR)
	mkdir -p $(REDUNDANS_DIR)
	cp $< $@

## Compare redundans results with old SPAdes assemblies (compare vs scaffolds for 004, not contigs)

$(REDUNDANS_DIR)/Tn004_S1_L001_spades_scaffolds.fasta : ../CV_samples_nextera/read_composition/filtered_assembly/004_spades/scaffolds.fasta
	mkdir -p $(REDUCTION_DIR)
	mkdir -p $(REDUNDANS_DIR)
	cp $< $@

$(REDUNDANS_DIR)/Tn019_S2_L001_spades_scaffolds.fasta : ../CV_samples_nextera/read_composition/filtered_assembly/019_spades/contigs.fasta
	mkdir -p $(REDUCTION_DIR)
	mkdir -p $(REDUNDANS_DIR)
	cp $< $@
	
$(REDUNDANS_DIR)/%/quast_comparison/report.html : $(REDUNDANS_DIR)/%_spades_scaffolds.fasta $(REDUNDANS_DIR)/%/scaffolds.filled.fa
	$(QUAST_EXE) -o $(REDUNDANS_DIR)/$*/quast_comparison -l SPAdes,Redundans -s $(REDUNDANS_DIR)/$*_spades_scaffolds.fasta $(REDUNDANS_DIR)/$*/scaffolds.filled.fa



### Align filtered reads against redundans assembly. Use scaffolds >= 1kbp to start with.
.PHONY : align_v_redundans
align_v_redundans : $(REDUNDANS_DIR)/Tn004_S1_L001/assembly1000.fasta $(REDUNDANS_DIR)/Tn004_S1_L001/assembly1000.1.bt2 $(REDUNDANS_DIR)/Tn004_S1_L001_v_assembly1000.bam $(REDUNDANS_DIR)/Tn019_S2_L001/assembly1000.fasta $(REDUNDANS_DIR)/Tn019_S2_L001/assembly1000.1.bt2 $(REDUNDANS_DIR)/Tn019_S2_L001_v_assembly1000.bam $(REDUNDANS_DIR)/Tn004_S1_L001_v_assembly1000_qualimap/qualimapReport.html $(REDUNDANS_DIR)/Tn019_S2_L001_v_assembly1000_qualimap/qualimapReport.html 

$(REDUNDANS_DIR)/%/assembly1000.fasta : $(REDUNDANS_DIR)/%/scaffolds.filled.fa scripts/filter_contigs_by_size.py
	python scripts/filter_contigs_by_size.py $(REDUNDANS_DIR)/$*/scaffolds.filled.fa 1000 $@
	
$(REDUNDANS_DIR)/%/assembly1000.1.bt2 : $(REDUNDANS_DIR)/%/assembly1000.fasta
	bowtie2-build $< $(REDUNDANS_DIR)/$*/assembly1000

$(REDUNDANS_DIR)/%_v_assembly1000.bam : $(FILTER_READ_DIR)/%_R1_trimmed.filtered.fastq.gz $(FILTER_READ_DIR)/%_R2_trimmed.filtered.fastq.gz $(ALIGN_SRC) $(REDUNDANS_DIR)/%/assembly1000.1.bt2
	$(ALIGN_EXE) -f $(FILTER_READ_DIR)/$*_R1_trimmed.filtered.fastq.gz -r $(FILTER_READ_DIR)/$*_R2_trimmed.filtered.fastq.gz -i $(REDUNDANS_DIR)/$*/assembly1000 -a $(REDUNDANS_DIR)/ -u $(REDUNDANS_DIR)/ 2> $(REDUNDANS_DIR)/$*_v_assembly1000_align.log
	
$(REDUNDANS_DIR)/%_qualimap/qualimapReport.html : $(REDUNDANS_DIR)/%_assembly1000.bam
	$(QUALIMAP_EXE) bamqc -bam $< -outdir $(REDUNDANS_DIR)/$*_qualimap



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
	@echo FILTERED_READS: $(FILTERED_READS)
	@echo FILTERED_QUALS: $(FILTERED_QUALS)
	@echo SPADES_ASSEMBLIES: $(SPADES_ASSEMBLIES)
	@echo IDBA_UD_ASSEMBLIES: $(IDBA_UD_ASSEMBLIES)
	@echo MEGAHIT_ASSEMBLIES: $(MEGAHIT_ASSEMBLIES)
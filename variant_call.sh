#!/bin/bash

## This is a shell script for variant calling using WES data
## Run this code preferably in the home "(~)" directory

## "PRE_REQUISITES"

### !!! install this pre-requisites by yourself, since you need to modify the .bashrc in some steps, and also some of them were initially installed in conda environment !!!

	# (optional) download and install "sratoolkit" (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)

	# install "fastqc"
		#sudo apt install fastqc

	# install "miniconda" (or, alternatives)
		#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
		#chmod +x Miniconda3-latest-Linux-x86_64.sh
		#bash Miniconda3-latest-Linux-x86_64.sh

	# install "samtools" (go inside conda environment first; preferably the "base")
		#conda install -c bioconda samtools
		#sudo apt install samtools

	# install "trimmomatic" (go insida conda environment first)
		#conda install -c bioconda trimmomatic
		#after installing trimmomatic in conda environment, if you want to run it outside conda, provide the file path to trimmomatic in the .bashrc

		# !!! a bug: when trimmomatric file path was provided in .bashrc, the samtools file path was getting disrupted in some way: need to fix !!!!

	# install "bwa"
		#sudo apt install bwa

	# download and install "gatk"
		#https://github.com/broadinstitute/gatk/releases
		#gatk pre-requisites-
			# "java 17 or greater"
				#sudo apt install openjdk-17-jdk
			# "python 2.6 or greater" (usually will be present)
		#gatk installation:
			#wget -c https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
			#unzip gatk-4.5.0.0.zip
			#provide the file path to gatk containing directory for user convenience

		# if you use gatk older version, it can be downloaded from bioconda
		# (must be run in a conda environment)
			#conda install -c bioconda gatk
			#gatk3 --help

	# (optional) install "picard" [in a conda environment] (only if you use gatk versions older than 4)
		#conda install -c bioconda picard

	# install "R" (for histogram generation in multi-qc analysis)
		#sudo apt install r-base

	# install "multiqc"
		#sudo apt install multiqc

	# download "functotator"
		#wget -c -P ~ https://storage.googleapis.com/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz
		#tar -xzvf ~/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz
		#cd ~/funcotator_dataSources.v1.8.hg38.20230908s
		#tar -zxf gnomAD_exome.tar.gz
		#tar -zxf gnomAD_genome.tar.gz


## "PROVIDE PERMISSION TO THIS SCRIPT"
#chmod 777 ~/variant_call.sh


## "MAKE DIRECTORIES FOR WELL ORGANIZATION OF DATA"

mkdir -p WES_project/breast_cancer/{reads,reference,data,results}
mkdir WES_project/breast_cancer/reads/{trimmed_reads,aligned_reads}

## "directories"

ref=~/WES_project/breast_cancer/reference
reads=~/WES_project/breast_cancer/reads
t_reads=~/WES_project/breast_cancer/reads/trimmed_reads
a_reads=~/WES_project/breast_cancer/reads/aligned_reads
known_sites=~/WES_project/breast_cancer/reference/Homo_sapiens_assembly38.dbsnp138.vcf.gz
data=~/WES_project/breast_cancer/data
results=~/WES_project/breast_cancer/results
func_path=~/funcotator_dataSources.v1.8.hg38.20230908s


## "Downloading and indexing Reference files"

wget -c -P ${ref} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ${ref}/hg38.fa.gz
samtools faidx ${ref}/hg38.fa

## "DOWNLOAD Sequence Reads (any of the two methods could be followed)"

### Method 1: first get the link from "SRA Explorer", then paste the link here-

wget -c -P ${reads} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/056/SRR25434456/SRR25434456_1.fastq.gz

wget -c -P ${reads} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR254/056/SRR25434456/SRR25434456_2.fastq.gz

### METHOD 2: alternatively the fastq files could be obtained using SRA toolkit of NCBI- (comment out the above two wget commands if you use this method)

#prefetch SRR25434456
#vdb-validate SRR25434456.sra
#fasterq-dump fasterq-dump SRR25434456.sra
#gzip *.fastq


## "QUALITY CHECK"

fastqc ${reads}/*fastq.gz


## "QUALITY CONTROL"
## !!!(this part depends on the quality of data, be sure to "modify" accordingly; for manual on trimmomatic: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) !!!

### The data was already in good shape as checked by fastqc. Some adapter sequences were trimmed at the 3` ends
trimmomatic PE -phred33 -threads 16 SRR25434456_1.fastq.gz SRR25434456_2.fastq.gz fwd_SRR25434456_p.fastq.gz fwd_SRR25434456_up.fastq.gz rev_SRR25434456_p.fastq.gz rev_SRR25434456_up.fastq.gz CROP:120 SLIDINGWINDOW:4:10
fastqc ${t_reads}/*fastq.gz


## "Reference Indexing and alignment"

bwa index ${ref}/hg38.fa
bwa mem -t 16 -R "@RG\tID:SRR25434456\tPL:ILLUMINA\tSM:SRR25434456" ${ref}/hg38.fa ${t_reads}/fwd_SRR25434456_p.fastq.gz ${t_reads}/rev_SRR25434456_p.fastq.gz > ${a_reads}/SRR25434456.paired.sam
	##use samtools to view the aligned file
		#samtools view SRR25434456.paired.sam


## "Duplicate marking and sorting"
samtools view -bS ${a_reads}/SRR25434456.paired.sam > ${a_reads}/SRR25434456.paired.bam
samtools sort ${a_reads}/SRR25434456.paired.bam -o ${a_reads}/SRR25434456.paired.sorted.bam
samtools rmdup ${a_reads}/SRR25434456.paired.sorted.bam ${a_reads}/SRR25434456.sorted.nodup.bam
samtools index ${a_reads}/SRR25434456.sorted.nodup.bam


## "Downloading known sites files for base quality score recalibration"
wget -c -P ${ref} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

## "sorting the known sites"
bgzip ${ref}/Homo_sapiens_assembly38.dbsnp138.vcf                 # ensure that you have bgzip with : sudo apt install tabix
tabix -p vcf ${ref}/Homo_sapiens_assembly38.dbsnp138.vcf.gz

## "Creating reference dictionary: used by gatk"
gatk CreateSequenceDictionary R=${ref}/hg38.fa O=${ref}/hg38.dict
## "or, (if you are using gatk versions older than 4, use picard)"
#picard CreateSequenceDictionary R=${ref}/hg38.fa O=${ref}/hg38.dict


## "Base quality recalibration"
### building the model:-
gatk BaseRecalibrator -I ${a_reads}/SRR25434456.sorted.nodup.bam -R ${ref}/hg38.fa -knownSites ${known_sites} -O ${data}/recal_data.table
### applying the model:-
gatk ApplyBQSR -I ${a_reads}/SRR25434456.sorted.nodup.bam -R ${ref}/hg38.fa --bqsr-recal-file ${data}/recal_data.table -O ${a_reads}/SRR25434456.sorted.nodup.bqsr.bam


## "Alignment Metrices and Insert Size Matrices for post alignment QC"
gatk CollectAlignmentSummaryMetrics R=${ref}/hg38.fa I=${a_reads}/SRR25434456.sorted.nodup.bqsr.bam O=${a_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${a_reads}/SRR25434456.sorted.nodup.bqsr.bam OUTPUT=${a_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${a_reads}/insert_size_histogram.pdf
multiqc -o ${a_reads} ${a_reads}/.


## "Variant calling"
gatk HaplotypeCaller -R ${ref}/hg38.fa -I ${a_reads}/SRR25434456.sorted.nodup.bqsr.bam -O ${results}/raw_variants.vcf

gatk SelectVariants -R ${ref}/hg38.fa -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref}/hg38.fa -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf


## "Variant filtering"
### SNP filtering:- (the filter scores are recommended by gatk on https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)
gatk VariantFiltration \
	-R ${ref}/hg38.fa \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \    # this showed a warning message: WARN  JexlEngine - ![0,9]: 'MQRankSum < -12.5;' undefined variable MQRankSum
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \  #this showed a warning message: WARN  JexlEngine - ![0,14]: 'ReadPosRankSum < -8.0;' undefined variable ReadPosRankSum
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"

### InDel filtering:-
gatk VariantFiltration \
	-R ${ref}/hg38.fa \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"



## "Variant selection based on the filter"
### SNP selection:-
gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis_ready_snps.vcf

### InDel selection:-
gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis_ready_indels.vcf

## "Variant selection considering the GT filters" (which were ignored by the SelectVariants)
cat ${results}/analysis_ready_snps.vcf | grep -v -E "DP_filter|GQ_filter" > ${results}/analysis_ready_snps_filteredGT.vcf
cat ${results}/analysis_ready_indels.vcf | grep -v -E "DP_filter|GQ_filter" > ${results}/analysis_ready_indels_filteredGT.vcf


## "Variant annotation" (for details of funcotator and its usage: https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator)

gatk Funcotator \
	--variant ${results}/analysis_ready_snps_filteredGT.vcf \
	--reference ${ref}/hg38.fa \
	--ref-version hg38 \
	--data-sources-path ${func_path} \
	--output ${results}/snps_functotated.vcf \
	--output-file-format VCF

gatk Funcotator \
	--variant ${results}/analysis_ready_indels_filteredGT.vcf \
	--reference ${ref}/hg38.fa \
	--ref-version hg38 \
	--data-sources-path ${func_path} \
	--output ${results}/indels_functotated.vcf \
	--output-file-format VCF


## "VCF file to table"

# for snps
gatk VariantsToTable \
	-V ${results}/snps_functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O ${results}/snps.table
# for indels
gatk VariantsToTable \
        -V ${results}/indels_functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/indels.table

## "Organizing the table using string manipulation"
# obtaining the header for the funcotator fields
grep "##INFO=<ID=FUNCOTATION" ${results}/snps_functotated.vcf | sed 's/.*Funcotation fields are: \(.*\)">/\1/' | sed 's/|/\t/g' > ${results}/curated_snps.txt
# appending the data from funcotator fields
cat ${results}/snps.table | cut -f 5 | sed 's/|/\t/g' >> ${results}/curated_snps.txt

# repeating the process for indels
grep "##INFO=<ID=FUNCOTATION" ${results}/indels_functotated.vcf | sed 's/.*Funcotation fields are: \(.*\)">/\1/' | sed 's/|/\t/g' > ${results}/curated_indels.txt
cat ${results}/indels.table | cut -f 5 | sed 's/|/\t/g' >> ${results}/curated_indels.txt

# WXS-workflow
Whole exome sequencing data analysis using GATK. Analysis is from fastq files to annotated SNPs.

### Source
- Project Source Data: Lee JH, Zhao XM, Yoon I, et al. Integrative analysis of mutational and transcriptional profiles reveals driver mutations of metastatic breast cancers. Cell Discovery. 2016;2:16025. DOI: 10.1038/celldisc.2016.25.
- SRA Project ID: PRJEB9083
- Metadata table: [SRA Run Selector table](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=8&WebEnv=MCID_6567e1fdc74a383354317088&o=acc_s%3Aa)

### Sample Details:
Total 78 individual with tumors contributed to the study.  Whole Exome is sequenced from each individual's tumor and Whole blood (as control).  RNA is siolated from the tumor cells.  Based on rsik of metastasis the samples were categorised as HRM (High Risk of Metastasis) and LRM (Low Risk of Metastasis). For further details on samples please refer to the [Metadata table](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=8&WebEnv=MCID_6567e1fdc74a383354317088&o=acc_s%3Aa). One individuals	`HRM_Tumor_WXS_6` is sequenced twice with run IDs ERR864175 and	ERR861009.  This makes the total 235.

Total libraries sequenced : 235

RNA sequnced libraries : 78

WXS samples : 157 (78 libraries from Tumor, 78 libraries from Whole blood and `HRM_Tumor_WXS_6` is sequenced twice)

The Scripts are located in a `scripts` directory with a `logs` directory to collect run logs. Scripts directory also have SRArunMetatable, here `WXS_Accession.txt`.  All the data and analysis directories are created outside scripts directory.


## Step1 Data Download

```
#!/bin/bash
#SBATCH -J WXS_Dwnld
#SBATCH -p general 
#SBATCH -q general 
#SBATCH -c 6
#SBATCH --mem=10G
#SBATCH -o ./logs/%x_%A_%a.out
#SBATCH --array=[0-156]%10

export TMPDIR="/scratch/vsingh/"

Accession="./WXS_Accession.txt"

# Column 1 of SRA Metadata table is the RUN id.
sampleID=($(awk '{print $1}' $Accession))

# Column 2 of SRA Metadata table is sample name.
samplename=($(awk '{print $2}' $Accession))

# Column 4 of SRA Metadata table describes assay type WXS or RNA.
AssayType=($(awk '{print $4}' $Accession))

# DIrectories to collect WXS and RNA data in seperate directories.
mkdir -p 01_rawdata_WXS 01_rawdata_RNA

# Setting variable based array job
SID=${sampleID[$SLURM_ARRAY_TASK_ID]}
SName=${samplename[$SLURM_ARRAY_TASK_ID]}
AT=${AssayType[$SLURM_ARRAY_TASK_ID]}

module load sratoolkit/3.0.2

fastq-dump --split-files --gzip ${SID}

# Renaming files and moving them to correct directory

if [ ${AT} == "WXS" ]
then
    mv ${SID}_1.fastq.gz 01_rawdata_WXS/${SNAME}_1.fastq.gz
    mv ${SID}_2.fastq.gz 01_rawdata_WXS/${SNAME}_2.fastq.gz
else
    mv ${SID}_1.fastq.gz 01_rawdata_RNA/${SNAME}_1.fastq.gz
    mv ${SID}_2.fastq.gz 01_rawdata_RNA/${SNAME}_2.fastq.gz 
fi

```

## Step2 Quality QC

```
#!/bin/bash
#SBATCH -J WXS01
#SBATCH -p general 
#SBATCH -q general 
#SBATCH -c 6
#SBATCH --mem=20G
#SBATCH -o ./logs/%x_%A_%a.out
#SBATCH --array=[0-156]%10

# Adapter file to be used in Trimommatic
ADAPTERFILE="/isg/shared/apps/Trimmomatic/0.39/adapters/NexteraPE-PE.fa"

export TMPDIR=/scratch/vsingh/

Accession="./WXS_Accession.txt"

#####################################################################
# Setting sample name to be processed
#####################################################################

samplename=($(awk '{print $2}' $Accession))
sampleID=${samplename[$SLURM_ARRAY_TASK_ID]}

#####################################################################
# Running FASTQC on Raw Reads
#####################################################################
mkdir -p ../01_fastqc

module load fastqc/0.11.7
fastqc -t 6 -o ../01_fastqc \
        ../01_rawdata_WXS/${sampleID}_1.fastq.gz \
        ../01_rawdata_WXS/${sampleID}_2.fastq.gz

#####################################################################
# Trimming Raw Reads
#####################################################################
mkdir -p ../02_trimreads/singles
module load Trimmomatic/0.39

java -jar $Trimmomatic PE -threads 6 \
        ../01_rawdata_WXS/${sampleID}_1.fastq.gz \
        ../01_rawdata_WXS/${sampleID}_2.fastq.gz \
        ../02_trimreads/trim_${sampleID}_R1.fastq.gz ../02_trimreads/singles/${sampleID}_singles.fastq.gz \
        ../02_trimreads/trim_${sampleID}_R2.fastq.gz ../02_trimreads/singles/${sampleID}_singles.fastq.gz \
        ILLUMINACLIP:${ADAPTERFILE}:2:30:10:5:true \
        SLIDINGWINDOW:4:20 MINLEN:45 HEADCROP:20


#####################################################################
# Running FASTQC on trimReads
#####################################################################
mkdir -p ../01_trimfastqc

module load fastqc

fastqc -t 4 -o ../01_trimfastqc \
    ../02_trimreads/trim_${sampleID}_R1.fastq.gz \
    ../02_trimreads/trim_${sampleID}_R1.fastq.gz


```


## Step3 Mapping of reads with BWA

```
#!/bin/bash
#SBATCH --job-name 02_bwaPicard
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --constraint="AVX512"
#SBATCH -p xeon
#SBATCH -q general
#SBATCH -o ./logs/%x_%A_%a 
#SBATCH --array=[0-156]%10

export TMPDIR=/scratch/vsingh/

Accession="./WXS_Accession.txt"

#####################################################################
# Setting sample name to be processed
#####################################################################

samplename=($(awk '{print $2}' $Accession))
INPUTFILE=${samplename[$SLURM_ARRAY_TASK_ID]}

# Directoty with trim reads
TRIMDIR=../02_trimreads/


#####################################################################
# Mapping reads
# NOTE
# BWA-mem2 failed on compute node with amd processors.
# Try it on Xeon processors.
# The reference genome fasta file is downloaded seperately and indexed with BWA.. 
# Refernce fasta file and BWA index location ../resources
#####################################################################

MAPDIR=../04_mapDNA
mkdir -p ${MAPDIR}
module load bwa-mem2/2.1
module load samtools

cd ../04_mapDNA

if [ -f "${INPUTFILE}.bam" ];then
        exit 0
else

        bwa-mem2 mem -t 10 \
                -R "@RG\tID:${INPUTFILE}\tPL:ILLUMINA\tSM:${INPUTFILE}" \
                ../resources/hs38 \
                $TRIMDIR/trim_${INPUTFILE}_R1.fastq.gz \
                $TRIMDIR/trim_${INPUTFILE}_R2.fastq.gz \
                | samtools view -bhS -@ 5 -o ${INPUTFILE}.bam -

        if [ $? -eq 0 ]; then
        echo "${INPUTFILE} mapping failed: Array JobID :${SLURM_ARRAY_TASK_ID}" >> FailedMappingJobs_${SLURM_ARRAY_JOB_ID}.out
        echo "${INPUTFILE}" >> FailedMappingAccession.txt
        else
        echo "${INPUTFILE}" >> SuccessfulMappinp.txt
        fi
fi

#####################################################################
# SORTing BAM file with Picard
#####################################################################
module load  picard/2.23.9
java -jar ${PICARD} SortSam \
        I=${INPUTFILE}.bam \
        O=${INPUTFILE}_sorted.bam \
        SORT_ORDER=coordinate

if [ $? -eq 0 ]; then
        rm ${INPUTFILE}.bam
else
        touch ${INPUTFILE}.bam_sorting_STEP1_FAILED
        exit 1
fi

```

## Step4 Marking Duplicates, Base Recalibration  and indelRealignment
```

#!/bin/bash
#SBATCH --job-name 02_bwaPicard
#SBATCH --cpus-per-task=15
#SBATCH --mem=40G
#SBATCH -p xeon
#SBATCH -q general
#SBATCH -o ./logs/%x_%A_%a 
#SBATCH --array=[0-156]%10


export TMPDIR=/scratch/vsingh/

Accession="./WXS_Accession.txt"

#####################################################################
# Setting sample name to be processed
#####################################################################

samplename=($(awk '{print $2}' $Accession))
INPUTFILE=${samplename[$SLURM_ARRAY_TASK_ID]}

#####################################################################
# Marking Duplicates
#####################################################################

cd ../04_mapDNA
module load  picard/2.23.9

java -jar ${PICARD} MarkDuplicates \
      --REMOVE_DUPLICATES=true \
      I=${INPUTFILE}_sorted.bam \
      O=${INPUTFILE}_duprem.bam \
      M=./stats/${INPUTFILE}_marked_dup_metrics.txt

if [ $? -eq 0 ]; then
        rm ${INPUTFILE}_sorted.bam
else
        touch ${INPUTFILE}_sorted.bam_MARKDUPLICATES_STEP2_FAILED
        exit 1
fi


java -jar ${PICARD} BuildBamIndex \
      I=${INPUTFILE}_duprem.bam


java -jar ${PICARD} BamIndexStats \
      I=${INPUTFILE}_duprem.bam \
      O=${INPUTFILE}

#####################################################################
# Base Recalibration
#####################################################################

mkdir -p ./BQSR
module load GATK/4.3.0.0

gatk BaseRecalibrator \
   -I ${INPUTFILE}_sorted_duprem.bam \
   -R ../resources/Homo_sapiens_assembly38.fasta \
   --known-sites ../Homo_sapiens_assembly38.dbsnp138.vcf \
   -O ./BQSR/${INPUTFILE}_recal_data.table

gatk ApplyBQSR \
   -R ../resources/Homo_sapiens_assembly38.fasta \
   -I ${INPUTFILE}_sorted_duprem.bam \
   --bqsr-recal-file ./BQSR/${INPUTFILE}_recal_data.table \
   -O ${INPUTFILE}_sorted_duprem_BQSR.bam

mkdir -p ./alignment_insert_metrics

java -Xmx60g -Djava.io.tmpdir=$TMPdir -jar picard.jar CollectAlignmentSummaryMetrics \
          R=../resources/Homo_sapiens_assembly38.fasta \
          I=${INPUTFILE}_duprem_BQSR.bam \
          O=./alignment_insert_metrics/${INPUTFILE}_Alignment_metrics.txt \
          -Xmx60g -Djava.io.tmpdir=$TMPdir

java -Xmx60g -Djava.io.tmpdir=$TMPdir -jar picard.jar CollectInsertSizeMetrics \
      I=${INPUTFILE}_duprem_BQSR.bam \
      O=./alignment_insert_metrics/${INPUTFILE}_insert_size_metrics.txt \
      H=./alignment_insert_metrics/${INPUTFILE}_insert_size_histogram.pdf \
      M=0.5 \
      TMP_DIR=$TMPdir

```

## Step5 Generating Panel of Normals (PON)

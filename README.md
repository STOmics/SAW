
# SAW: Stereo-seq Analysis Workflow
Workflow for analyzing Stereo-seq transcriptomic data. Stereo-seq Analysis Workflow (SAW) software suite is a set of pipelines bundled to map sequenced reads to their spatial location on the tissue section, quantify the corresponding gene expression levels and visually present spatial gene expression distribution.

**DockerHub Link: https://hub.docker.com/r/stomics/saw/tags**

##  Introduction
SAW processes the sequencing data of Stereo-seq to generate spatial gene expression matrices, and the users could take these files as the starting point to perform downstream analysis. SAW includes thirteen essential and suggest pipelines and auxiliary tools for supporting other handy functions.
![workflow.png](SAW_v5.1.3_workflow.jpg)

##  System Requirements
###   Hardware
Stereo-seq Analysis Workflow (SAW) should be run on a Linux system that meets the following requirements:
* 8-core Intel or AMD processor (24 cores recommended)
* 128GB RAM (256GB recommended)
* 1TB free disk space
* 64-bit CentOS/RedHat 7.8 or Ubuntu 20.04

###   Software
* Singularity: a container platform
* SAW in the Singularity Image File (SIF) format
* ImageQC version greater than v1.1.0

####   Quick installation of Singularity
```
## On Red Hat Enterprise Linux or CentOS install the following dependencies:
$ sudo yum update -y && \
     sudo yum groupinstall -y 'Development Tools' && \
     sudo yum install -y \
     openssl-devel \
     libuuid-devel \
     libseccomp-devel \
     wget \
     squashfs-tools \
     cryptsetup

## On Ubuntu or Debian install the following dependencies:
$ sudo apt-get update && sudo apt-get install -y \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin

## Install Go
$ export VERSION=1.14.12 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

$ echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

## Install singularity on CentOS without compile
$ yum install -y singularity
```
**For additional help or support, please visit https://sylabs.io/guides/3.8/admin-guide/installation.html**

####   Quick download SAW from DockerHub
Currently, the latest version of SAW is v5.1.3. You can download SAW by running the following command:
```
singularity build SAW_v5.1.3.sif docker://stomics/saw:05.1.3
```
_The bash script for v5.1.3 is coming very soon._



##   [Preparation](https://github.com/BGIResearch/SAW/tree/main/script/pre_buildIndexedRef)
###    Build index for reference genome
A genome index has to be constructed before performing data mapping. The index files are used as reference when aligning reads. You can prepare the indexed reference before run SAW as follow:
```
singularity exec <SAW_v5.1.3.sif> mapping --runMode genomeGenerate \
    --genomeDir reference/STAR_SJ100 \
    --genomeFastaFiles reference/genome.fa \
    --sjdbGTFfile reference/genes.gtf \
    --sjdbOverhang 99 \
    --runThreadN 12
```
**For more information, refer to "script/pre_buildIndexedRef"**

###    Get Stereo-seq Chip T mask file
- If you want to access mask file (.h5/.bin) for your own data, please contact BGI-FAS team.
- To access mask file for published paper, please go to [CNGBdb](https://db.cngb.org/) > [STOmicsDB](https://db.cngb.org/stomics) > [Collections](https://db.cngb.org/stomics/collections).



##  RUN
_The RUN examples and bash script for v5.1.3 are coming very soon._


### Usage
```
# for saw_v5.1.3
usage: sh <stereoRun.sh> -m maskFile -1 read1 -2 read2 -g indexedGenome -a annotationFile -o outDir -i image -t threads -s visualSif -c genomeSize
    -genomeSize : genome size
    -splitCount : count of splited stereochip mask file, usually 16 for SE+Q4 fq data and 1 for PE+Q40 fq data
    -maskFile : stereochip mask file
    -fq1 : fastq file path of read1, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_1.fq.gz,lane2_read_1.fq.gz
    -fq2 : fastq file path of read2, if there are more than one fastq file, please separate them with comma, not requested for 'SE+Q4' fastq data, e.g:lane1_read_2.fq.gz,lane2_read_2.fq.gz    
    -refIndex : reference genome indexed folder, please build IT before SAW analysis run
    -speciesName : specie of the sample
    -tissueType : tissue type of the sample
    -annotationsFile :  annotations file in gff or gtf format, the file must contain gene and exon annotations
    -outDir : output directory path
    -imageRecordFile : image file(*.ipr) generated by ImageQC software, not requested
    -imageCompressedFile : image file(*.tar.gz) generated by ImageQC software, not requested 
    -doCellBin : [Y/N]
    -threads : the number of threads to be used in running the pipeline
    -sif : the file format of the visual software

# 1GiB=1024M=10241024KB=10241024*1024B
# SAW version : v5.1.3
```
```
# for saw_beta_v4.1.0 & saw_v4.1.0
usage: sh <stereoRun.sh> -m maskFile -1 read1 -2 read2 -g indexedGenome -a annotationFile -o outDir -i image -t threads -s visualSif -c genomeSize
    -m stereochip mask file
    -1 fastq file path of read1, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_1.fq.gz,lane2_read_1.fq.gz
    -2 fastq file path of read2, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_2.fq.gz,lane2_read_2.fq.gz
    -g genome that has been indexed by star
    -a annotation file in gff or gtf format, the file must contain gene and exon annotation, and also the transript annotation
    -o output directory path
    -i image directory path, must contains SN*.tar.gz and SN*.json file generated by ImageQC software, not required
    -t thread number will be used to run this pipeline
    -s docker image that packed analysis softwares
    -c genome fasta file size (GiB, Gibibyte)

# 1GiB=1024M=10241024KB=10241024*1024B
# SAW version : v4.1.0
```

###   Example: Running the entire workflow
For SAW_v4, please use the [stereoRun_singleLane_v4.1.0.sh](https://github.com/BGIResearch/SAW/blob/main/script/stereoRun_singleLane_v4.1.0.sh) or [stereoRun_multiLane_v4.1.0.sh](https://github.com/BGIResearch/SAW/blob/main/script/stereoRun_multiLane_v4.0.0.sh) to run whole workflow.
For SAW_v5, please use the [stereoPipeline.sh](https://github.com/BGIResearch/SAW/blob/main/script/stereoPipeline.sh) to run whole workflow.

####    Run stereoPipeline.sh bash script
```
cd <个人工作目录>

ulimit -n 10240
ulimit -v 33170449147
NUMBA_CACHE_DIR=<Your Working Directory>

dataDir=<Your Working Directory>/rawData
outDir=<Your Working Directory>/resultresult

export SINGULARITY_BIND=$dataDir,$outDir
bash stereoPipeline.sh \
    -genomeSize 5 \
    -splitCount 1 \
    -maskFile $dataDir/mask/SN.h5 \
    -fq1 $dataDir/reads/lane1_read_1.fq.gz,...,$dataDir/reads/laneN_read_1.fq.gz  \
    -fq2 $dataDir/reads/lane1_read_2.fq.gz,...,$dataDir/reads/laneN_read_2.fq.gz \ # [option] when the sequenced data is in PE format
    -speciesName <speciesName> \
    -tissueType <tissueName> \
    -refIndex $dataDir/reference/STAR_SJ100 \
    -annotationFile $dataDir/reference/genes.gtf \
    -imageRecordFile $dataDir/SN/image_dir_path/<imageQC result>.ipr \ # [option] when tissue image was given
    -imageCompressedFile $dataDir/SN/image_dir_path/<imageQC result>.tar.gz \ # [option] when tissue image was given
    -sif $dataDir/SAW/SAW_version.sif \
    -threads 16 \
    -doCellBin Y \ # [option] when you want to do the cellBin part
    -outDir $outDir/result
```


####    Run stereoRun_singleLane_v4.1.0.sh bash script
If only one lane sequencing data was given, run the stereoRun_singleLane_v4.1.0.sh bash script as the following:
```
ulimit -n 10240 
dataDir=/Full/Path/Of/Input/File 
outDir=/Full/Path/Of/Output/File 
export SINGULARITY_BIND=$dataDir,$outDir
bash stereoRun_singleLane.sh \
    -m $dataDir/mask/SN.h5 \
    -1 $dataDir/reads/lane1_read_1.fq.gz \
    -2 $dataDir/reads/lane1_read_2.fq.gz \
    -g $dataDir/reference/STAR_SJ100 \
    -a $dataDir/reference/genes.gtf \
    -s $dataDir/SAW/SAW_version.sif \
    -c genome_size \ # (GiB, Gibibyte)
    -i $dataDir/SN/image_dir_path \ # [option] when tissue image was given
    -o $outDir/result
    
# 1GiB=1024M=10241024KB=10241024*1024B
# SAW version : v4.1.0
```
####    Run stereoRun_multiLane_v4.1.0.sh bash script
If more than one lane sequencing data was given, run the stereoRun_multiLane_v4.1.0.sh script as the following:
```
ulimit -n 10240 
dataDir=/Full/Path/Of/Input/File 
outDir=/Full/Path/Of/Output/File 
export SINGULARITY_BIND=$dataDir,$outDir
bash stereoRun_multiLane.sh \
    -m $dataDir/mask/SN.h5 \
    -1 $dataDir/reads/lane1_read_1.fq.gz,$dataDir/reads/lane2_read_1.fq.gz \
    -2 $dataDir/reads/lane1_read_2.fq.gz,$dataDir/reads/lane2_read_2.fq.gz \
    -g $dataDir/reference/STAR_SJ100 \
    -a $dataDir/reference/genes.gtf \
    -s $dataDir/SAW/SAW_version.sif \
    -c genome_size \ # (GiB, Gibibyte)
    -i $dataDir/SN/image_dir_path \ # [option] when tissue image was given
    -o $outDir/result
    

# 1GiB=1024M=10241024KB=10241024*1024B
# SAW version : v4.1.0
```
_The bash script for v5.1.3 is coming very soon._

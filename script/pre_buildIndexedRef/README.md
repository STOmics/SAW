# SAW : Stereo-seq Analysis Workflow
Building your reference before running SAW analysis.
##  Preparation : Indexing a reference genome
For more information visit https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

##  Run
###   Prepare input files
```
referenceDir=/Full/Path/Of/Reference/Folder/Path
refName=/Reference/Specie/Name
# make a new directory called reference specie name
mkdir -p $referenceDir/$refName
cd $referenceDir/$refName
mkdir genome genes

# download the genome fasta file, like "Genome sequence, primary assembly (GRCh38)"
cd genome
wget <ftp link for genome fasta file> #wget https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# download the annotation file that correspond to your reference genome (.gtf/.gff)
cd ../genes
wget <ftp link for annotation file in gtf/gff format> #wget https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz 
gunzip Homo_sapiens.GRCh38.93.gtf.gz
```
###   Run : build the genome index file
```
referenceDir=/Full/Path/Of/Reference/Folder/Path
export SINGULARITY_BIND=$referenceDir

singularity exec <SAW_version.sif> mapping \
    --runMode genomeGenerate \
    --genomeDir $reference/STAR_SJ100 \
    --genomeFastaFiles $referenceDir/genome/genome.fa \
    --sjdbGTFfile $referenceDir/genes/genes.gtf \
    --sjdbOverhang 99 \
    --runThreadN 12
Then you should get the mask file from our website through the slide number(SN)
```

###   Output files
Once the run is completed , you will have a STAR_SJ100 folder in the reference directory.

####    Description of each folder    
```
cd /data/dataManagement/reference
tree ./
.
└── specieName
    ├── STAR_SJ100
    │   ├── Genome
    │   ├── SA
    │   ├── SAindex
    │   ├── chrLength.txt
    │   ├── chrName.txt
    │   ├── chrNameLength.txt
    │   ├── chrStart.txt
    │   ├── exonGeTrInfo.tab
    │   ├── exonInfo.tab
    │   ├── geneInfo.tab
    │   ├── genomeParameters.txt
    │   ├── sjdbInfo.txt
    │   ├── sjdbList.fromGTF.out.tab
    │   ├── sjdbList.out.tab
    │   └── transcriptInfo.tab
    ├── genes
    │   └── genes.gtf
    └── genome
        └── genome.fa

4 directories, 17 files
```

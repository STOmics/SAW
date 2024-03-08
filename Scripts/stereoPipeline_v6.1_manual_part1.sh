#!/bin/bash
set -e

if [[ $# -lt 12 ]];then
    echo "usage: sh $0 -genomeSize -splitCount -maskFile -fq1 -fq2 -speciesName -tissueType -refIndex -annotationFile -imageRecordFile -imageCompressedFile -sif -threads -outDir
    -genomeSize : genome size
    -splitCount : count of splited stereochip mask file, usually 16 for SE+Q4 fq data and 1 for PE+Q40 fq data
    -maskFile : stereochip mask file
    -fq1 : fastq file path of read1, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_1.fq.gz,lane2_read_1.fq.gz
    -fq2 : fastq file path of read2, if there are more than one fastq file, please separate them with comma, not requested for 'SE+Q4' fastq data, e.g:lane1_read_2.fq.gz,lane2_read_2.fq.gz
    -speciesName : specie of the sample  
    -tissueType : tissue type of the sample	
	-refIndex : reference genome indexed folder, please build IT before SAW analysis run
    -annotationFile :  annotations file in gff or gtf format, the file must contain gene and exon annotations
    -imageRecordFile : image file(*.ipr) generated by ImageStudio software, not requested
    -imageCompressedFile : image file(*.tar.gz) generated by ImageStudio software, not requested
    -sif : the file format of the visual software
    -threads : the number of threads to be used in running the pipeline
    -outDir : output directory path
    "
    exit
fi

while [[ -n "$1" ]]
do
    case "$1" in
        -genomeSize) GSize="$2"
            shift ;;
        -splitCount) splitCnt="$2"
            shift ;;
        -maskFile) maskFile="$2"
            shift ;;
        -fq1) read1="$2"
            shift ;;
        -fq2) read2="$2"
            shift ;;
        -refIndex) GDir="$2"
            shift ;;
        -speciesName) refName="$2"
            shift ;;
        -tissueType) tissueType="$2"
            shift ;;
        -annotationFile) annoFile="$2"
            shift ;;
        -imageRecordFile) iprFile="$2"
            shift ;;
        -imageCompressedFile) imageTarFile="$2"
            shift ;;
        -outDir) outDir="$2"
            shift ;;
        -threads) threads="$2"
            shift ;;
        -sif) sif="$2"
            shift ;;
    esac
        shift
done


# Software Check
if [ `command -v singularity` ]
then
    singularityPath=`command -v singularity`
    echo `date` " singularity check: pass, and singularity path is ${singularityPath}"
else
    echo `date` " singularity check: singularity does not exits, please verify that you have installed singularity and exported it to your system PATH variable"
    exit
fi

if [[ -n $sif ]]
then
    echo `date` " singularity image file check: file exist and SIF path is ${sif}"
else
    echo `date` " singularity image file check: file does not exist, please double check your SIF file is in the current directory or the path given by the option -s is valid."
fi


if [[ ! -d $outDir ]];then
    mkdir -p $outDir
fi

# Get basic information
maskname=$(basename $maskFile)
SN=${maskname%%.*}

maskDIR=$(dirname $maskFile)
annoDIR=$(dirname $annoFile)
refDIR=$(dirname $GDir)

# Prepare output directories
result_00mapping=${outDir}/00.mapping
result_01merge=${outDir}/01.merge
result_02count=${outDir}/02.count
result_03register=${outDir}/03.register

arr_result=( $result_00mapping $result_01merge $result_02count $result_03register )
for each in "${arr_result[@]}";
do
    if [[ ! -d $each ]];then
        mkdir -p $each
    fi
done

## Run splitMask or CIDCount for preparation
# ulimit -c 100000000000
echo `date` "=> splitMask, compute CID count and predict the memory of mapping start......"
read1List=(`echo $read1 | tr ',' ' '`)
fqbases=()
starBams=()
bcStat=()
bcLogFinalOut=()
bcReadsCounts=()

if [[ ! -n "$read2" ]]; then
    fqType="SE"
    arr_result=( ${result_00mapping}/splitBin ${result_00mapping}/mergeList )
    for each in "${arr_result[@]}";do
        if [[ ! -d $each ]];then mkdir -p $each; fi
    done
    export SINGULARITY_BIND=$outDir,$maskDIR,$annoDIR,$refDIR
    /usr/bin/time -v singularity exec ${sif} splitMask \
        ${maskFile} ${result_00mapping}/splitBin $threads $splitCnt 2_25
    for ((i=1;i<=$splitCnt;i++)); do
        if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
        echo $read1 | sed 's/,/\n/g' | grep _$i.fq.gz > ${result_00mapping}/mergeList/$a.${SN}.Q4_SE.fq.list || continue
    done
else
    fqType="PE"
    read2List=(`echo $read2 | tr ',' ' '`)
    fqNumber=`echo ${#read1List[@]}`
    export SINGULARITY_BIND=$outDir,$maskDIR,$annoDIR,$refDIR
    /usr/bin/time -v singularity exec ${sif} CIDCount \
        -i ${maskFile} \
        -s ${refName} \
        -g ${GSize} > ${result_00mapping}/CIDCount
fi
echo "Your sequencing reads are in ${fqType} format."


## Run SAW mapping to perform CID mapping and STAR alignment.
echo `date` "=> CID mapping, adapter filtering and RNA alignment start......"
if [[ $fqType == 'PE' ]]; then
    for ((i=0;i<=`expr $(echo $fqNumber) - 1`;i++)); do
        fqname=$(basename ${read1List[i]})
        fqdir=$(dirname ${read1List[i]})
        fqbase=${fqname%%.*}
        fqbases[i]=$fqbase
        bcPara=${result_00mapping}/${fqbase}.bcPara
        barcodeReadsCount=${result_00mapping}/${fqbase}.barcodeReadsCount.txt
        echo  " ~~~ mapping - $fqname ~~~"
        echo "in=${maskFile}" > $bcPara
        echo "in1=${read1List[i]}" >> $bcPara
        echo "in2=${read2List[i]}" >> $bcPara
        echo "encodeRule=ACTG" >> $bcPara
        echo "action=4" >> $bcPara
        echo "barcodeReadsCount=${barcodeReadsCount}" >> $bcPara
        echo "platform=T10" >> $bcPara
        echo "out=${fqbase}" >> $bcPara
        echo "barcodeStart=0" >> $bcPara
        echo "barcodeLen=25" >> $bcPara
        echo "umiStart=25" >> $bcPara
        echo "umiLen=10" >> $bcPara
        echo "umiRead=1" >> $bcPara
        echo "mismatch=1" >> $bcPara
        echo "bcNum=`head -1 ${result_00mapping}/CIDCount`" >> $bcPara
        echo "polyAnum=15" >> $bcPara
        echo "mismatchInPolyA=2" >> $bcPara
        read1DIR=$(dirname ${read1List[i]})
        read2DIR=$(dirname ${read2List[i]})
        export SINGULARITY_BIND=$read1DIR,$read2DIR,$outDir,$maskDIR,$annoDIR,$refDIR
        /usr/bin/time -v singularity exec ${sif} mapping \
            --outSAMattributes spatial \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir ${GDir} \
            --runThreadN ${threads} \
            --outFileNamePrefix ${result_00mapping}/${fqbase}. \
            --sysShell /bin/bash \
            --stParaFile ${bcPara} \
            --readNameSeparator \" \" \
            --limitBAMsortRAM 63168332971 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize=480000000 \
            --outBAMsortingBinsN 50 > ${result_00mapping}/${fqbase}.run.log

        starBam=${result_00mapping}/${fqbase}.Aligned.sortedByCoord.out.bam
        starBams[i]=$starBam
        bcStat[i]=${result_00mapping}/${fqbase}.CIDMap.stat
        bcFinalOut[i]=${result_00mapping}/${fqbase}.Log.final.out
        bcReadsCounts[i]=$barcodeReadsCount
    done
elif [[ $fqType == 'SE' ]]; then
    for ((i=1;i<=$splitCnt;i++)); do
        if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
        /usr/bin/time -v singularity exec ${sif} CIDCount \
            -i $(ls ${result_00mapping}/splitBin/${a}.${SN}.barcodeToPos.bin) \
            -s ${refName} \
            -g ${GSize} > ${result_00mapping}/CIDCount
        fqbase=$a.${SN}.Q4_SE
        bcPara=${result_00mapping}/${fqbase}.bcPara
        barcodeReadsCount=${result_00mapping}/${fqbase}.barcodeReadsCount.txt
        echo  " ~~~ mapping SE FASTQ Data - $a ~~~"
        echo "in=$(ls ${result_00mapping}/splitBin/${a}.${SN}.barcodeToPos.bin)" > $bcPara
        read1List=${result_00mapping}/mergeList/$a.${SN}.Q4_SE.fq.list
        if [[ ! -s $read1List ]];then continue;else echo $read1List;fi
        echo "in1=${read1List}" >> $bcPara
        read1DIR=$(dirname $(cat $read1List)|tr '\n' ',')
        echo "encodeRule=ACTG" >> $bcPara
        echo "action=4" >> $bcPara
        echo "out=${fqbase}" >> $bcPara
        echo "barcodeReadsCount=${barcodeReadsCount}" >> $bcPara
        echo "platform=T10" >> $bcPara
        echo "barcodeStart=0" >> $bcPara
        echo "barcodeLen=24" >> $bcPara
        echo "umiStart=25" >> $bcPara
        echo "umiLen=10" >> $bcPara
        echo "umiRead=1" >> $bcPara
        echo "mismatch=1" >> $bcPara
        echo "bcNum=`head -1 ${result_00mapping}/CIDCount`" >> $bcPara
        echo "polyAnum=15" >> $bcPara
        echo "mismatchInPolyA=2" >> $bcPara
        export SINGULARITY_BIND=$read1DIR,$outDir,$maskDIR,$annoDIR,$refDIR
        /usr/bin/time -v singularity exec ${sif} mapping \
            --outSAMattributes spatial \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir ${GDir} \
            --runThreadN ${threads} \
            --outFileNamePrefix ${result_00mapping}/${fqbase}. \
            --sysShell /bin/bash \
            --stParaFile ${bcPara} \
            --readNameSeparator \" \" \
            --limitBAMsortRAM 63168332971 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize=280000000 \
            --outSAMmultNmax 1 \
            --outBAMsortingBinsN 50 > ${result_00mapping}/${fqbase}.run.log

        starBams[i]=${result_00mapping}/${fqbase}.Aligned.sortedByCoord.out.bam
        bcStat[i]=${result_00mapping}/${fqbase}.CIDMap.stat
        bcFinalOut[i]=${result_00mapping}/${fqbase}.Log.final.out
        bcReadsCounts[i]=$barcodeReadsCount
    done
fi

if [[ $(echo ${#bcReadsCounts[*]}) == '1' ]]; then
    bcReadsCountsStr=$bcReadsCounts
    starBamsStr=${starBams[0]}
    bcFinalOutStr=${bcFinalOut[0]}
    bcStatStr=${bcStat[0]}
else
    bcReadsCountsStr=$( IFS=','; echo "${bcReadsCounts[*]}" )
    starBamsStr=$( IFS=','; echo "${starBams[*]}" )
    bcFinalOutStr=$( IFS=','; echo "${bcFinalOut[*]}" )
    bcStatStr=$( IFS=','; echo "${bcStat[*]}" )
fi


## RUN SAW merge to merge barcodeReadsCount file
echo `date` "=> merge barcode reads count tables start......"
barcodeReadsCounts=${result_01merge}/${SN}.merge.barcodeReadsCount.txt
export SINGULARITY_BIND=$outDir,$maskDIR
if [[ $fqType == 'SE' ]] && [[ $(echo ${#bcReadsCounts[*]}) > '1' ]]; then
    /usr/bin/time -v singularity exec ${sif} merge \
        ${maskFile} \
        $bcReadsCountsStr \
        $barcodeReadsCounts
elif [[ $fqType == 'PE' ]]; then
    if [[ $(echo ${#bcReadsCounts[*]}) == '1' ]]
    then
        cp $bcReadsCountsStr $barcodeReadsCounts
    else
        /usr/bin/time -v singularity exec ${sif} merge \
            ${maskFile} \
            $bcReadsCountsStr \
            $barcodeReadsCounts
    fi
fi


## Run SAW count to do annotation, deduplication, and generate gene expression matrix
echo `date` "=> annotation, deduplication, and generate gene expression matrix start......"
geneExp=${result_02count}/${SN}.raw.gef
saturationFile=${result_02count}/${SN}_raw_barcode_gene_exp.txt
export SINGULARITY_BIND=$outDir,$annoDIR,$refDIR
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} count \
    -i ${starBamsStr} \
    -o ${result_02count}/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a ${annoFile} \
    -s ${result_02count}/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e ${geneExp} \
    --sat_file ${saturationFile} \
    --umi_on \
    --umi_len 10 \
    --save_lq \
    --save_dup \
    --sn ${SN} \
    -c ${threads} \
    -m 128
	
	
## Run cellcut for GEF completion
mkdir -p ${result_03register}/manual_register
export SINGULARITY_BIND=$outDir
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} cellCut bgef \
    -i ${geneExp} \
    -o ${result_03register}/manual_register/${SN}.gef \
    -O Transcriptomics

## Run ipr2img for TIFF
imgTarDIR=$(dirname $imageTarFile)
iprDIR=$(dirname $iprFile)
export SINGULARITY_BIND=$outDir,$imgTarDIR,$iprDIR
/usr/bin/time -v singularity exec ${sif} imageTools ipr2img \
            -i $imageTarFile \
            -c $iprFile \
            -r False \
            -o ${result_03register}/manual_register

## Run img2rpi for RPI 
registerTif=$(find ${result_03register}/manual_register -maxdepth 1 -name \*fov_stitched_transformed.tif)
if [[ -n $registerTif ]]
then
    regGroup=$(find ${result_03register}/manual_register -maxdepth 1 -name \*fov_stitched_transformed.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} +)
else
    registerTif=$(find ${result_03register}/manual_register -maxdepth 1 -name \*fov_stitched.tif)
    regGroup=$(find ${result_03register}/manual_register -maxdepth 1 -name \*fov_stitched.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched.tif;done' sh {} +)
fi
regTifStr=$(echo $registerTif | tr ' ' ',')
regGroupStr=$(echo $regGroup | sed 's/ \|$/\/Image,/g' | sed 's/.$//')
echo $regTifStr
echo $regGroupStr

/usr/bin/time -v singularity exec ${sif} imageTools img2rpi \
        -i ${regTifStr} \
        -g ${regGroupStr} \
        -b 1 10 50 100 \
        -o ${result_03register}/manual_register/fov_stitched.rpi


echo `date` " Now you got GEF and RPI for manual registration. "

#!/bin/bash
set -e

if [[ $# -lt 12 ]];then
    echo "usage: sh $0 -genomeSize -splitCount -maskFile -rnafq1 -rnafq2 -adtfq1 -adtfq2 -proteinList -refIndex -annotationFile -speciesName -tissueType -rRNAremove -imageRecordFile -imageCompressedFile -sif -threads -outDir
    -genomeSize : genome size
    -splitCount : count of splited stereochip mask file, usually 16 for Q4 fq data and 1 for Q40 fq data
    -maskFile : stereochip mask file
    -rnaFq1 : RNA fastq file path of read1, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_1.fq.gz,lane2_read_1.fq.gz
    -rnaFq2 : RNA fastq file path of read2, if there are more than one fastq file, please separate them with comma, not requested for Q4 fastq data, e.g:lane1_read_2.fq.gz,lane2_read_2.fq.gz
    -adtFq1 : ADT fastq file path of read1, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_1.fq.gz,lane2_read_1.fq.gz
    -adtFq2 : ADT fastq file path of read2, if there are more than one fastq file, please separate them with comma, not requested for Q4 fastq data, e.g:lane1_read_2.fq.gz,lane2_read_2.fq.gz
    -proteinList : protein list file which contain protein sequences and names.
	  -refIndex : reference genome indexed folder, please build IT before SAW analysis run
    -annotationFile :  annotations file in gff or gtf format, the file must contain gene and exon annotations
    -speciesName : specie of the sample
    -tissueType : tissue type of the sample
    -rRNAremove : [Y/N]
    -pidStart : PID start position. 21 for sequencing transcriptome and ADT libraries together. 0 for sequencing separately.
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
        -splitCount) splitCnt="$2"
            shift ;;
        -maskFile) maskFile="$2"
            shift ;;
        -rnaFq1) rnaRead1="$2"
            shift ;;
        -rnaFq2) rnaRead2="$2"
            shift ;;
        -adtFq1) adtRead1="$2"
            shift ;;
        -adtFq2) adtRead2="$2"
            shift ;;
        -proteinList) proteinList="$2"
            shift ;;
        -speciesName) refName="$2"
            shift ;;
        -tissueType) tissueType="$2"
            shift ;;
        -refIndex) GDir="$2"
            shift ;;
        -annotationFile) annoFile="$2"
            shift ;;
        -rRNAremove) rRNAremove="$2"
            shift ;;
        -imageRecordFile) iprFile="$2"
            shift ;;
        -imageCompressedFile) imageTarFile="$2"
            shift ;;
        -pidStart) pidStart="$2"
            shift ;;
        -sif) sif="$2"
            shift ;;
        -threads) threads="$2"
            shift ;;
        -outDir) outDir="$2"
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


# Get basic information
maskname=$(basename $maskFile)
SN=${maskname%%.*}

proteinListDIR=$(dirname $proteinList)
maskDIR=$(dirname $maskFile)
annoDIR=$(dirname $annoFile)
refDIR=$(dirname $GDir)

# Prepare output directories
if [[ ! -d $outDir ]];then
    mkdir -p $outDir
fi

arr_result=(${outDir}/00T.mapping ${outDir}/00P.mapping ${outDir}/01T.merge ${outDir}/02T.count ${outDir}/03.calibration  ${outDir}/04.register )
for each in "${arr_result[@]}";
do
    if [[ ! -d $each ]];then
        mkdir -p $each
    fi
done


# Run SAW splitMask or CIDCount for preparation
echo `date` "=> splitMask, compute CID count and predict the memory of mapping start......"
rnaRead1List=(`echo $rnaRead1 | tr ',' ' '`)
fqbases=()
starBams=()
bcStat=()
bcLogFinalOut=()
bcReadsCounts=()
GSize=(`du -sh --block-size=G ${GDir}/Genome | cut -f 1`)

if [[ ! -n "$rnaRead2" ]]; then
    fqType="Q4"
    arr_result=( ${outDir}/00T.mapping/splitBin ${outDir}/00T.mapping/mergeList )
    for each in "${arr_result[@]}";do
        if [[ ! -d $each ]];then mkdir -p $each; fi
    done
    export SINGULARITY_BIND=$outDir,$maskDIR,$annoDIR,$refDIR
    /usr/bin/time -v singularity exec ${sif} splitMask \
        ${maskFile} ${outDir}/00T.mapping/splitBin $threads $splitCnt 2_25
    for ((i=1;i<=$splitCnt;i++)); do
        if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
        echo $rnaRead1 | sed 's/,/\n/g' | grep _$i.fq.gz > ${outDir}/00T.mapping/mergeList/$a.${SN}.Q4.fq.list
    done
else
    fqType="Q40"
    rnaRead2List=(`echo $rnaRead2 | tr ',' ' '`)
    fqNumber=`echo ${#rnaRead1List[@]}`
    export SINGULARITY_BIND=$outDir,$maskDIR,$annoDIR,$refDIR
    /usr/bin/time -v singularity exec ${sif} CIDCount \
        -i ${maskFile} \
        -s ${refName} \
        -g ${GSize} > ${outDir}/00T.mapping/CIDCount
fi
echo "Your rna sequencing reads are in ${fqType} format."

## Run SAW mapping to perform CID mapping and STAR alignment
echo `date` "=> CID mapping, adapter filtering and RNA alignment start......"
if [[ $fqType == 'Q40' ]]; then
    for ((i=0;i<=`expr $(echo $fqNumber) - 1`;i++)); do
        fqname=$(basename ${rnaRead1List[i]})
        fqdir=$(dirname ${rnaRead1List[i]})
        fqbase=${fqname%%.*}
        fqbases[i]=$fqbase
        bcPara=${outDir}/00T.mapping/${fqbase}.bcPara
        barcodeReadsCount=${outDir}/00T.mapping/${fqbase}.barcodeReadsCount.txt
        echo  " ~~~ mapping - $fqname ~~~"
        echo "in=${maskFile}" > $bcPara
        echo "in1=${rnaRead1List[i]}" >> $bcPara
        echo "in2=${rnaRead2List[i]}" >> $bcPara
        echo "barcodeReadsCount=${barcodeReadsCount}" >> $bcPara
        echo "barcodeStart=0" >> $bcPara
        echo "barcodeLen=25" >> $bcPara
        echo "umiStart=25" >> $bcPara
        echo "umiLen=10" >> $bcPara
        echo "mismatch=1" >> $bcPara
        echo "bcNum=`head -1 ${outDir}/00T.mapping/CIDCount`" >> $bcPara
        echo "polyAnum=15" >> $bcPara
        echo "mismatchInPolyA=2" >> $bcPara
        if [[ $rRNAremove == "Y" ]]; then
            echo "rRNAremove" >> $bcPara
        fi
        rnaRead1DIR=$(dirname ${rnaRead1List[i]})
        rnaRead2DIR=$(dirname ${rnaRead2List[i]})
        export SINGULARITY_BIND=$rnaRead1DIR,$rnaRead2DIR,$outDir,$maskDIR,$annoDIR,$refDIR
        /usr/bin/time -v singularity exec ${sif} mapping \
            --outSAMattributes spatial \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir ${GDir} \
            --runThreadN ${threads} \
            --outFileNamePrefix ${outDir}/00T.mapping/${fqbase}. \
            --sysShell /bin/bash \
            --stParaFile ${bcPara} \
            --readNameSeparator \" \" \
            --limitBAMsortRAM 63168332971 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize=280000000 \
            --outBAMsortingBinsN 50 \
            --outSAMmultNmax 1 \
            > ${outDir}/00T.mapping/${fqbase}.run.log

        starBam=${outDir}/00T.mapping/${fqbase}.Aligned.sortedByCoord.out.bam
        starBams[i]=$starBam
        bcStat[i]=${outDir}/00T.mapping/${fqbase}.CIDMap.stat
        bcFinalOut[i]=${outDir}/00T.mapping/${fqbase}.Log.final.out
        bcReadsCounts[i]=$barcodeReadsCount
    done
elif [[ $fqType == 'Q4' ]]; then
    for ((i=1;i<=$splitCnt;i++)); do
        if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
        /usr/bin/time -v singularity exec ${sif} CIDCount \
            -i $(ls ${outDir}/00T.mapping/splitBin/${a}.${SN}.barcodeToPos.bin) \
            -s ${refName} \
            -g ${GSize} > ${outDir}/00T.mapping/CIDCount
        fqbase=$a.${SN}.Q4
        bcPara=${outDir}/00T.mapping/${fqbase}.bcPara
        barcodeReadsCount=${outDir}/00T.mapping/${fqbase}.barcodeReadsCount.txt
        echo  " ~~~ mapping Q4 FASTQ Data - $a ~~~"
        echo "in=$(ls ${outDir}/00T.mapping/splitBin/${a}.${SN}.barcodeToPos.bin)" > $bcPara
        rnaRead1List=${outDir}/00T.mapping/mergeList/$a.${SN}.Q4.fq.list
        echo "in1=${rnaRead1List}" >> $bcPara
        rnaRead1DIR=$(dirname $(cat $rnaRead1List)|tr '\n' ',')
        echo "barcodeReadsCount=${barcodeReadsCount}" >> $bcPara
        echo "barcodeStart=0" >> $bcPara
        echo "barcodeLen=24" >> $bcPara
        echo "umiStart=25" >> $bcPara
        echo "umiLen=10" >> $bcPara
        echo "mismatch=1" >> $bcPara
        echo "bcNum=`head -1 ${outDir}/00T.mapping/CIDCount`" >> $bcPara
        echo "polyAnum=15" >> $bcPara
        echo "mismatchInPolyA=2" >> $bcPara
        if [[ $rRNAremove == "Y" ]]; then
            echo "rRNAremove" >> $bcPara
        fi
        export SINGULARITY_BIND=$rnaRead1DIR,$outDir,$maskDIR,$annoDIR,$refDIR
        /usr/bin/time -v singularity exec ${sif} mapping \
            --outSAMattributes spatial \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir ${GDir} \
            --runThreadN ${threads} \
            --outFileNamePrefix ${outDir}/00T.mapping/${fqbase}. \
            --sysShell /bin/bash \
            --stParaFile ${bcPara} \
            --readNameSeparator \" \" \
            --limitBAMsortRAM 63168332971 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize=280000000 \
            --outBAMsortingBinsN 50 \
            --outSAMmultNmax 1 \
            > ${outDir}/00T.mapping/${fqbase}.run.log

        starBams[i]=${outDir}/00T.mapping/${fqbase}.Aligned.sortedByCoord.out.bam
        bcStat[i]=${outDir}/00T.mapping/${fqbase}.CIDMap.stat
        bcFinalOut[i]=${outDir}/00T.mapping/${fqbase}.Log.final.out
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


# Run SAW merge to integrate barcodeReadsCount file
echo `date` "=> merge barcode reads count tables start......"
barcodeReadsCounts=${outDir}/01T.merge/${SN}.merge.barcodeReadsCount.txt
export SINGULARITY_BIND=$outDir,$maskDIR
if [[ $fqType == 'Q4' ]] && [[ $(echo ${#bcReadsCounts[*]}) > '1' ]]; then
    echo 'Q4'
    /usr/bin/time -v singularity exec ${sif} merge \
        ${maskFile} \
        $bcReadsCountsStr \
        $barcodeReadsCounts
elif [[ $fqType == 'Q40' ]]; then
    echo 'Q40'
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

# Run SAW count to do annotation, deduplication, and generate gene expression matrix
echo `date` "=> annotation, deduplication, and generate gene expression matrix start......"
saturationFile=${outDir}/02T.count/${SN}_raw_barcode_gene_exp.txt
export SINGULARITY_BIND=$outDir,$annoDIR,$refDIR
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} count \
    -i ${starBamsStr} \
    -o ${outDir}/02T.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a ${annoFile} \
    -s ${outDir}/02T.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e ${outDir}/02T.count/${SN}.raw.gef \
    --sat_file ${saturationFile} \
    --umi_on \
    --umi_len 10 \
    --save_lq \
    --save_dup \
    --sn ${SN} \
    -c ${threads}


# Run SAW mapping-SP to perform CID mapping, MID filtration and PID alignment
echo `date` "=> CID mapping and PID alignment start......"
echo `date` "=> Please make sure your protein list only records actual input proteins"
echo $adtRead1|sed 's/,/\n/g'>${outDir}/00P.mapping/adtFq1.list
adtRead1List=${outDir}/00P.mapping/adtFq1.list
adtRead1DIR=$(dirname $(cat $adtRead1List)|tr '\n' ',')

if [[ ! -n "$adtRead2"  ]]; then
    adtFqType="Q4"
    echo "Your sequencing reads are in ${adtFqType} format."
    find ${outDir}/00T.mapping/splitBin/*${SN}.barcodeToPos.bin>${outDir}/00P.mapping/splitMask.txt
    export HDF5_USE_FILE_LOCKING=FALSE
    export SINGULARITY_BIND=$adtRead1DIR,$outDir,$maskDIR,${proteinListDIR}
    /usr/bin/time -v singularity exec ${sif} mapping-SP \
        --thread ${threads} \
        --cidStart 0 --cidLen 24 --cidMismatch 1 \
        --midStart 25 --midLen 10 \
        --pidStart ${pidStart} --pidLen 15 --pidMismatch 1 \
        --mask ${outDir}/00P.mapping/splitMask.txt \
        --fastq ${outDir}/00P.mapping/adtFq1.list \
        --reference ${proteinList} \
        --output ${outDir}/00P.mapping \
        --sn ${SN}
else
    adtFqType="Q40"
    echo "Your sequencing reads are in ${adtFqType} format."
    echo "~~~ mapping Q40 ${adtRead1} ~~~"
    echo $adtRead2|sed 's/,/\n/g'>${outDir}/00P.mapping/adtFq2.list
    adtRead2List=${outDir}/00P.mapping/adtFq2.list
    adtRead2DIR=$(dirname $(cat $adtRead2List)|tr '\n' ',')
    export HDF5_USE_FILE_LOCKING=FALSE
    export SINGULARITY_BIND=$adtRead1DIR,$adtRead2DIR,$outDir,$maskDIR,${proteinListDIR}
    /usr/bin/time -v singularity exec ${sif} mapping-SP \
        --thread ${threads} \
        --cidStart 0 --cidLen 25 --cidMismatch 1 \
        --midStart 25 --midLen 10 \
        --pidStart ${pidStart} --pidLen 15 --pidMismatch 1 \
        --mask ${maskFile} \
        --fastq ${outDir}/00P.mapping/adtFq1.list \
        --fastq2 ${outDir}/00P.mapping/adtFq2.list \
        --reference ${proteinList} \
        --output ${outDir}/00P.mapping \
        --sn ${SN}
fi
gzip -f ${outDir}/00P.mapping/${SN}.protein.raw.gem

# Run SAW calibration to do multi-GEF alignment
echo `date` "=> calibration start......."
export HDF5_USE_FILE_LOCKING=FALSE
export SINGULARITY_BIND=$outDir
/usr/bin/time -v singularity exec ${sif} calibration \
    -i ${outDir}/00P.mapping/${SN}.protein.raw.gef,${outDir}/02T.count/${SN}.raw.gef \
    -o ${outDir}/03.calibration/${SN}.protein.calibrated.raw.gef,${outDir}/03.calibration/${SN}.calibrated.raw.gef \
    -O Proteomics,Transcriptomics

## Convert calibrated GEF to GEM
export HDF5_USE_FILE_LOCKING=FALSE
export SINGULARITY_BIND=$outDir
/usr/bin/time -v singularity exec ${sif} cellCut view \
    -s ${SN} \
    -i ${outDir}/03.calibration/${SN}.calibrated.raw.gef \
    -o ${outDir}/03.calibration/${SN}.calibrated.gem
gzip -f ${outDir}/03.calibration/${SN}.calibrated.gem

export HDF5_USE_FILE_LOCKING=FALSE
export SINGULARITY_BIND=$outDir
/usr/bin/time -v singularity exec ${sif} cellCut view \
    -s ${SN} \
    -i ${outDir}/03.calibration/${SN}.protein.calibrated.raw.gef \
    -o ${outDir}/03.calibration/${SN}.protein.calibrated.gem
gzip -f ${outDir}/03.calibration/${SN}.protein.calibrated.gem

# Run cellcut for GEF completion
mkdir -p ${outDir}/04.register/manual_register
export SINGULARITY_BIND=$outDir
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} cellCut bgef \
    -i ${outDir}/03.calibration/${SN}.calibrated.raw.gef \
    -o ${outDir}/04.register/manual_register/${SN}.gef \
    -O Transcriptomics


# Run ipr2img for TIFF
imgTarDIR=$(dirname $imageTarFile)
iprDIR=$(dirname $iprFile)
export SINGULARITY_BIND=$outDir,$imgTarDIR,$iprDIR
/usr/bin/time -v singularity exec ${sif} imageTools ipr2img \
            -i $imageTarFile \
            -c $iprFile \
            -r False \
            -o ${outDir}/04.register/manual_register


# Run img2rpi for RPI

registerTif=$(find ${outDir}/04.register/manual_register -maxdepth 1 -name \*fov_stitched_transformed.tif)
if [[ -n $registerTif ]]
then
    regGroup=$(find ${outDir}/04.register/manual_register -maxdepth 1 -name \*fov_stitched_transformed.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} +)
else
    registerTif=$(find ${outDir}/04.register/manual_register -maxdepth 1 -name \*fov_stitched.tif)
    regGroup=$(find ${outDir}/04.register/manual_register -maxdepth 1 -name \*fov_stitched.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched.tif;done' sh {} +)
fi
regTifStr=$(echo $registerTif | tr ' ' ',')
regGroupStr=$(echo $regGroup | sed 's/ \|$/\/Image,/g' | sed 's/.$//')
echo $regTifStr
echo $regGroupStr

/usr/bin/time -v singularity exec ${sif} imageTools img2rpi \
        -i ${regTifStr} \
        -g ${regGroupStr} \
        -b 1 10 50 100 \
        -o ${outDir}/04.register/manual_register/fov_stitched.rpi


echo `date` " Now you got GEF and RPI files for manual registration. "

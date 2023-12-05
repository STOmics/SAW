#!/bin/bash
set -e

if [[ $# -lt 12 ]];then
    echo "usage: sh $0 -splitCount -maskFile -fq1 -fq2 -refIndex -genomeFile -speciesName -tissueType -annotationFile -outDir -imageRecordFile -imageCompressedFile -doCellBin -rRNARemove -threads -sif
    -splitCount : count of splited stereochip mask file, usually 16 for SE+Q4 fq data and 1 for PE+Q40 fq data
    -maskFile : stereochip mask file
    -fq1 : fastq file path of read1, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_1.fq.gz,lane2_read_1.fq.gz
    -fq2 : fastq file path of read2, if there are more than one fastq file, please separate them with comma, not requested for 'Q4' fastq data, e.g:lane1_read_2.fq.gz,lane2_read_2.fq.gz
    -refIndex : reference genome indexed folder, please build index before SAW analysis run
    -speciesName : specie of the sample
    -tissueType : tissue type of the sample
    -annotationFile :  annotations file in gff or gtf format, the file must contain gene and exon annotations
    -outDir : output directory path
    -imageRecordFile : image file(*.ipr) generated by ImageQC/ImageStudio software, not requested
    -imageCompressedFile : image file(*.tar.gz) generated by ImageQC/ImageStudio software, not requested  
    -doCellBin : [Y/N]
    -rRNAremove : [Y/N]
    -threads : the number of threads to be used in running the pipeline
    -sif : the file format of the visual software
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
        -outDir) outDir="$2"
            shift ;;
        -imageRecordFile) iprFile="$2"
            shift ;;
        -imageCompressedFile) imageTarFile="$2"
            shift ;;
        -doCellBin) doCell="$2"
            shift ;;
        -rRNAremove) rRNAremove="$2"
            shift ;;
        -threads) threads="$2"
            shift ;;
        -sif) sif="$2"
            shift ;;
        # -preRegDir) preRegDir="$2"
        #     shift ;;
    esac
        shift
done


# Software check
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
if [[ $iprFile ]] && [[ $imageTarFile ]];then
    iprDIR=$(dirname $iprFile)
    imgTarDIR=$(dirname $imageTarFile)
fi

# Prepare output directories
arr_result=( ${outDir}/00.mapping ${outDir}/01.merge ${outDir}/02.count ${outDir}/04.tissuecut ${outDir}/05.spatialcluster ${outDir}/06.saturation ${outDir}/07.report ${outDir}/visualization)
for each in "${arr_result[@]}";
do
    if [[ ! -d $each ]];then
        mkdir -p $each
    fi
done

if [[ $iprFile ]] && [[ $imageTarFile ]];then
    mkdir -p ${outDir}/03.register
fi

if [[ $doCell == "Y" ]]; then
    mkdir -p ${outDir}/041.cellcut
    mkdir -p ${outDir}/051.cellcluster
fi


# Run SAW splitMask or CIDCount for preparation

echo `date` "=> splitMask, compute CID count and predict the memory of mapping start......"
read1List=(`echo $read1 | tr ',' ' '`)
fqbases=()
starBams=()
bcStat=()
bcLogFinalOut=()
bcReadsCounts=()
GSize=(`du -sh --block-size=G ${GDir}/Genome | cut -f 1`)

if [[ ! -n "$read2" ]]; then
    fqType="Q4"
    arr_result=( ${outDir}/00.mapping/splitBin ${outDir}/00.mapping/mergeList )
    for each in "${arr_result[@]}";do
        if [[ ! -d $each ]];then mkdir -p $each; fi
    done
    export SINGULARITY_BIND=$outDir,$maskDIR,$annoDIR,$refDIR
    /usr/bin/time -v singularity exec ${sif} splitMask \
        ${maskFile} ${outDir}/00.mapping/splitBin $threads $splitCnt 2_25
    for ((i=1;i<=$splitCnt;i++)); do
        if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
        echo $read1 | sed 's/,/\n/g' | grep _$i.fq.gz > ${outDir}/00.mapping/mergeList/$a.${SN}.Q4.fq.list
    done
else
    fqType="Q40"
    read2List=(`echo $read2 | tr ',' ' '`)
    fqNumber=`echo ${#read1List[@]}`
    export SINGULARITY_BIND=$outDir,$maskDIR,$annoDIR,$refDIR
    /usr/bin/time -v singularity exec ${sif} CIDCount \
        -i ${maskFile} \
        -s ${refName} \
        -g ${GSize} > ${outDir}/00.mapping/CIDCount
fi
echo "Your sequencing reads are in ${fqType} format."


# Run SAW mapping to perform CID mapping and STAR alignment
echo `date` "=> CID mapping, adapter filtering and RNA alignment start......"
if [[ $fqType == 'Q40' ]]; then
    for ((i=0;i<=`expr $(echo $fqNumber) - 1`;i++)); do
        fqname=$(basename ${read1List[i]})
        fqdir=$(dirname ${read1List[i]})
        fqbase=${fqname%%.*}
        fqbases[i]=$fqbase
        bcPara=${outDir}/00.mapping/${fqbase}.bcPara
        barcodeReadsCount=${outDir}/00.mapping/${fqbase}.barcodeReadsCount.txt
        echo  " ~~~ mapping - $fqname ~~~"
        echo "in=${maskFile}" > $bcPara
        echo "in1=${read1List[i]}" >> $bcPara
        echo "in2=${read2List[i]}" >> $bcPara
        echo "barcodeReadsCount=${barcodeReadsCount}" >> $bcPara
        echo "barcodeStart=0" >> $bcPara
        echo "barcodeLen=25" >> $bcPara
        echo "umiStart=25" >> $bcPara
        echo "umiLen=10" >> $bcPara
        echo "mismatch=1" >> $bcPara
        echo "bcNum=`head -1 ${outDir}/00.mapping/CIDCount`" >> $bcPara
        echo "polyAnum=15" >> $bcPara
        echo "mismatchInPolyA=2" >> $bcPara
        if [[ $rRNAremove == "Y" ]]; then
            echo "rRNAremove" >> $bcPara
        fi
        read1DIR=$(dirname ${read1List[i]})
        read2DIR=$(dirname ${read2List[i]})
        export SINGULARITY_BIND=$read1DIR,$read2DIR,$outDir,$maskDIR,$annoDIR,$refDIR
        /usr/bin/time -v singularity exec ${sif} mapping \
            --outSAMattributes spatial \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir ${GDir} \
            --runThreadN ${threads} \
            --outFileNamePrefix ${outDir}/00.mapping/${fqbase}. \
            --sysShell /bin/bash \
            --stParaFile ${bcPara} \
            --readNameSeparator \" \" \
            --limitBAMsortRAM 63168332971 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize=280000000 \
            --outBAMsortingBinsN 50 \
            --outSAMmultNmax 1 \
            > ${outDir}/00.mapping/${fqbase}.run.log

        starBam=${outDir}/00.mapping/${fqbase}.Aligned.sortedByCoord.out.bam
        starBams[i]=$starBam
        bcStat[i]=${outDir}/00.mapping/${fqbase}.CIDMap.stat
        bcFinalOut[i]=${outDir}/00.mapping/${fqbase}.Log.final.out
        bcReadsCounts[i]=$barcodeReadsCount
    done
elif [[ $fqType == 'Q4' ]]; then
    for ((i=1;i<=$splitCnt;i++)); do
        if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
        /usr/bin/time -v singularity exec ${sif} CIDCount \
            -i $(ls ${outDir}/00.mapping/splitBin/${a}.${SN}.barcodeToPos.bin) \
            -s ${refName} \
            -g ${GSize} > ${outDir}/00.mapping/CIDCount
        fqbase=$a.${SN}.Q4
        bcPara=${outDir}/00.mapping/${fqbase}.bcPara
        barcodeReadsCount=${outDir}/00.mapping/${fqbase}.barcodeReadsCount.txt
        echo  " ~~~ mapping Q4 FASTQ Data - $a ~~~"
        echo "in=$(ls ${outDir}/00.mapping/splitBin/${a}.${SN}.barcodeToPos.bin)" > $bcPara
        read1List=${outDir}/00.mapping/mergeList/$a.${SN}.Q4.fq.list
        echo "in1=${read1List}" >> $bcPara
        read1DIR=$(dirname $(cat $read1List)|tr '\n' ',')
        echo "barcodeReadsCount=${barcodeReadsCount}" >> $bcPara
        echo "barcodeStart=0" >> $bcPara
        echo "barcodeLen=24" >> $bcPara
        echo "umiStart=25" >> $bcPara
        echo "umiLen=10" >> $bcPara
        echo "mismatch=1" >> $bcPara
        echo "bcNum=`head -1 ${outDir}/00.mapping/CIDCount`" >> $bcPara
        echo "polyAnum=15" >> $bcPara
        echo "mismatchInPolyA=2" >> $bcPara
        if [[ $rRNAremove == "Y" ]]; then
            echo "rRNAremove" >> $bcPara
        fi
        export SINGULARITY_BIND=$read1DIR,$outDir,$maskDIR,$annoDIR,$refDIR
        /usr/bin/time -v singularity exec ${sif} mapping \
            --outSAMattributes spatial \
            --outSAMtype BAM SortedByCoordinate \
            --genomeDir ${GDir} \
            --runThreadN ${threads} \
            --outFileNamePrefix ${outDir}/00.mapping/${fqbase}. \
            --sysShell /bin/bash \
            --stParaFile ${bcPara} \
            --readNameSeparator \" \" \
            --limitBAMsortRAM 63168332971 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize=280000000 \
            --outBAMsortingBinsN 50 \
            --outSAMmultNmax 1 \
            > ${outDir}/00.mapping/${fqbase}.run.log

        starBams[i]=${outDir}/00.mapping/${fqbase}.Aligned.sortedByCoord.out.bam
        bcStat[i]=${outDir}/00.mapping/${fqbase}.CIDMap.stat
        bcFinalOut[i]=${outDir}/00.mapping/${fqbase}.Log.final.out
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
barcodeReadsCounts=${outDir}/01.merge/${SN}.merge.barcodeReadsCount.txt
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


# Run SAW count to perform annotation, deduplication, and generate gene expression matrix
echo `date` "=> annotation, deduplication, and generate gene expression matrix start......"
export SINGULARITY_BIND=$outDir,$annoDIR,$refDIR
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} count \
    -i ${starBamsStr} \
    -o ${outDir}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a ${annoFile} \
    -s ${outDir}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e ${outDir}/02.count/${SN}.raw.gef \
    --sat_file ${outDir}/02.count/${SN}_raw_barcode_gene_exp.txt \
    --umi_on \
    --umi_len 10 \
    --save_lq \
    --save_dup \
    --sn ${SN} \
    -c ${threads}

# Run SAW register to stitch microscope tile images into a panoramic image, perform tissue and cell (optional, depends on -doCellBin) segmentation, and register the panoramic image and the segmented images with the gene expression matrix.
if [[ -f $imageTarFile ]] && [[ -f $iprFile ]]  && [[ $doCell == "Y" ]]; then
    ## Run register (stitch, tissue segmentation, cell segmentation) + imageTools
    echo `date` "=> image processing and registration start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    imgTarDIR=$(dirname $imageTarFile)
    iprDIR=$(dirname $iprFile)
    export SINGULARITY_BIND=$outDir,$imgTarDIR,$iprDIR

    /usr/bin/time -v singularity exec ${sif} register \
        -i ${imageTarFile} \
        -c ${iprFile} \
        -v ${outDir}/02.count/${SN}.raw.gef \
        -w True --core ${threads} \
        -o ${outDir}/03.register
    
    out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)    
    /usr/bin/time -v singularity exec ${sif} imageTools ipr2img \
        -i ${imageTarFile} \
        -c ${out_iprFile} \
        -d tissue cell \
        -r True \
        -o ${outDir}/03.register
    
    registerTif=$(find ${outDir}/03.register -maxdepth 1 -name \*fov_stitched_transformed.tif)
    regTifStr=$(echo $registerTif | tr ' ' ',')
    regGroup=$(find ${outDir}/03.register -maxdepth 1 -name \*fov_stitched_transformed.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} +)
    regGroupStr=$(echo $regGroup | sed 's/ \|$/\/Image,/g' | sed 's/.$//')
    echo $regTifStr
    echo $regGroupStr
    /usr/bin/time -v singularity exec ${sif} imageTools img2rpi \
        -i ${regTifStr} \
        -g ${regGroupStr} \
        -b 1 10 50 100 \
        -o ${outDir}/03.register/fov_stitched_transformed.rpi
elif [[ -f $imageTarFile ]] && [[ -f $iprFile ]]  && [[ $doCell == "N" ]]; then
    ## Run register (stitch, tissue segmentation) + imageTools
    echo `date` "=> image processing and registration start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    imgTarDIR=$(dirname $imageTarFile)
    iprDIR=$(dirname $iprFile)

    export SINGULARITY_BIND=$outDir,$imgTarDIR,$iprDIR
    /usr/bin/time -v singularity exec ${sif} register \
        -i ${imageTarFile} \
        -c ${iprFile} \
        -v ${outDir}/02.count/${SN}.raw.gef \
        -w False --core ${threads} \
        -o ${outDir}/03.register
    out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)
    /usr/bin/time -v singularity exec ${sif} imageTools ipr2img \
        -i ${imageTarFile} \
        -c ${out_iprFile} \
        -d tissue \
        -r True \
        -o ${outDir}/03.register
    registerTif=$(find ${outDir}/03.register -maxdepth 1 -name \*fov_stitched_transformed.tif)
    regTifStr=$(echo $registerTif | tr ' ' ',')
    regGroup=$(find ${outDir}/03.register -maxdepth 1 -name \*fov_stitched_transformed.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} +)
    regGroupStr=$(echo $regGroup | sed 's/ \|$/\/Image,/g' | sed 's/.$//')
    echo $regTifStr
    echo $regGroupStr
    /usr/bin/time -v singularity exec ${sif} imageTools img2rpi \
        -i ${regTifStr} \
        -g ${regGroupStr} \
        -b 1 10 50 100 \
        -o ${outDir}/03.register/fov_stitched_transformed.rpi
fi


# Run SAW tissueCut
echo `date` "=> tissueCut start......."
if [[ -f $imageTarFile ]] && [[ -f $iprFile ]]; then
    ## Run tissueCut to get the spatial gene expression profile of the tissue-covered region
    nucleusLayer=$(find ${outDir}/03.register -maxdepth 1 -name \*fov_stitched*.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched*.tif;done' sh {} + | grep -v IF | awk -F_ '{print$1}')
    tissueMaskFile=$(find ${outDir}/03.register -maxdepth 1 -name *_${SN}_tissue_cut.tif)
    export HDF5_USE_FILE_LOCKING=FALSE
    export SINGULARITY_BIND=$outDir
    /usr/bin/time -v singularity exec ${sif} tissueCut \
        -i ${outDir}/02.count/${SN}.raw.gef \
        --dnbfile ${barcodeReadsCounts} \
        -s ${tissueMaskFile} \
        --sn ${SN} \
        -O Transcriptomics \
        -d -t 8 \
        -o ${outDir}/04.tissuecut

    export HDF5_USE_FILE_LOCKING=FALSE
    /usr/bin/time -v singularity exec ${sif} cellCut bgef \
        -i ${outDir}/04.tissuecut/${SN}.tissue.gef \
        -o ${outDir}/04.tissuecut/${SN}.${nucleusLayer}.gef \
        -O Transcriptomics \
        -b 1,5,10,20,50,100,150,200

    # output labeled GEF [optional]
    out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)
    for i in `singularity exec ${sif} h5dump -n $out_iprFile | grep 'Labeling/' | grep -v 'Labeling/.*/Canvas'|grep -v 'Labeling/.*/Element'|awk '{print$2}'`;
    do
    label=`basename $i`
    mkdir -p ${outDir}/04.tissuecut/tissuecut_${label}
    labelmask=$(find ${outDir}/03.register -name \*${label}_tissue_cut.tif)
    echo $labelmask
    export HDF5_USE_FILE_LOCKING=FALSE
    export SINGULARITY_BIND=$outDir
    /usr/bin/time -v singularity exec ${sif} tissueCut \
        -l $label \
        -i ${outDir}/02.count/${SN}.raw.gef \
        --dnbfile ${barcodeReadsCounts} \
        -s ${labelmask} \
        --sn ${SN} \
        -O Transcriptomics \
        -d -t 8 \
        -o ${outDir}/04.tissuecut/tissuecut_${label}
    done
else
    ## Run tissueCut based on the gene expression matrix directly
    export SINGULARITY_BIND=$outDir,$annoDIR,$refDIR
    echo `date` "=> there is no image, tissueCut based on the gene expression matrix start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    /usr/bin/time -v singularity exec ${sif} tissueCut \
        -i ${outDir}/02.count/${SN}.raw.gef \
        --dnbfile ${barcodeReadsCounts} \
        --sn ${SN} \
        -O Transcriptomics \
        -d -t 8 \
        -o ${outDir}/04.tissuecut
fi

# Complete raw GEF to visual GEF
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} cellCut bgef \
    -i ${outDir}/02.count/${SN}.raw.gef \
    -o ${outDir}/02.count/${SN}.gef \
    -O Transcriptomics \
    -b 1,5,10,20,50,100,150,200

# Complete raw labeled tissue GEF to visual GEF 
out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)
for i in `singularity exec ${sif} h5dump -n $out_iprFile | grep 'Labeling/' | grep -v 'Labeling/.*/Canvas'|grep -v 'Labeling/.*/Element'|awk '{print$2}'`;
do
label=`basename $i`
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} cellCut bgef \
    -i ${outDir}/04.tissuecut/tissuecut_${label}/${SN}.${label}.raw.label.gef \
    -o ${outDir}/04.tissuecut/tissuecut_${label}/${SN}.${label}.label.gef \
    -O Transcriptomics
done

## Convert GEF to GEM [optional]
# /usr/bin/time -v singularity exec ${sif} cellCut view \
#     -s ${SN} \
#     -i ${outDir}/02.count/${SN}.gef \
#     -o ${outDir}/02.count/${SN}.gem
# gzip ${outDir}/02.count/${SN}.gem
# /usr/bin/time -v singularity exec ${sif} cellCut view \
#     -s ${SN} \
#     -i ${outDir}/04.tissuecut/${SN}.tissue.gef \
#     -o ${outDir}/04.tissuecut/${SN}.tissue.gem
# gzip ${outDir}/04.tissuecut/${SN}.tissue.gem


# Run SAW spatialCluster
binSize=200
resolution=1.0
export SINGULARITY_BIND=$outDir
echo `date` "=> spatialCluster start......."
export HDF5_USE_FILE_LOCKING=FALSE
mkdir -p ${outDir}/tmp
export NUMBA_CACHE_DIR=${outDir}/tmp
export MPLCONFIGDIR=${outDir}/tmp

/usr/bin/time -v singularity exec ${sif} spatialCluster \
    -i ${outDir}/04.tissuecut/${SN}.tissue.gef \
    -s ${binSize} \
	-r ${resolution} \
	-o ${outDir}/05.spatialcluster/${SN}.bin${binSize}_${resolution}.spatial.cluster.h5ad

## output labeled spatial H5AD [optional]
# out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)
# for i in `singularity exec ${sif} h5dump -n $out_iprFile | grep 'Labeling/' | grep -v 'Labeling/.*/Canvas'|grep -v 'Labeling/.*/Element'|awk '{print$2}'`;
# do
# label=`basename $i`
# echo $labelmask
# export HDF5_USE_FILE_LOCKING=FALSE
# /usr/bin/time -v singularity exec ${sif} spatialCluster \
#     -i ${outDir}/04.tissuecut/tissuecut_${label}/${SN}.${label}.label.gef \
#     -s ${binSize} \
#     -r ${resolution} \
#     -o ${outDir}/05.spatialcluster/${SN}.${label}.bin${binSize}_${resolution}.spatial.cluster.h5ad
# done

# Run SAW cellCut, cellCorrect, cellCluster and cellChunk
export SINGULARITY_BIND=$outDir
if [[ $doCell == 'Y' ]]; then
    echo `date` "=> cellCut start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    nucleusLayer=$(find ${outDir}/03.register -maxdepth 1 -name \*fov_stitched*.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched*.tif;done' sh {} + | grep -v IF | awk -F_ '{print$1}')
    nucleusMask=$(find ${outDir}/03.register -maxdepth 1 -name ${nucleusLayer}\*_mask.tif)
    /usr/bin/time -v singularity exec ${sif} cellCut cgef \
        -i ${outDir}/02.count/${SN}.raw.gef \
        -m ${nucleusMask} \
        -o ${outDir}/041.cellcut/${SN}.cellbin.gef

    echo `date` "=> cellCorrect start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    /usr/bin/time -v singularity exec ${sif} cellCorrect \
            -i ${outDir}/02.count/${SN}.raw.gef \
            -m ${nucleusMask} \
            -d 10 \
            -o ${outDir}/041.cellcut

    ## Write the cellCorrect mask image into SN.rpi
    cellCorrectMask=$(find ${outDir}/041.cellcut -maxdepth 1 -name ${nucleusLayer}\*_mask_edm_dis\*.tif)
    /usr/bin/time -v singularity exec ${sif} imageTools img2rpi \
        -i ${cellCorrectMask} \
        -g ${nucleusLayer}/CellMask_adjusted \
        -b 2 10 50 100 150 \
        -o ${outDir}/03.register/${SN}.rpi

    ### Convert cellbin GEF to GEM
    # /usr/bin/time -v singularity exec ${sif} cellCut view \
    #     -i ${outDir}/041.cellcut/${SN}.cellbin.gef \
    #     -d ${outDir}/02.count/${SN}.gef \
    #     -o ${outDir}/041.cellcut/${SN}.cellbin.gem \
    #     -s ${SN}

    echo `date` "=> cellCluster start......."
    mkdir -p ${outDir}/tmp
    export NUMBA_CACHE_DIR=${outDir}/tmp
    export MPLCONFIGDIR=${outDir}/tmp
    
    /usr/bin/time -v singularity exec ${sif} cellCluster \
        -i ${outDir}/041.cellcut/${SN}.adjusted.cellbin.gef \
        -o ${outDir}/051.cellcluster/${SN}.adjusted.cell.cluster.h5ad

    echo `date` "=> cellChunk start......."
    # Write rendering data into cellbin.gef
    /usr/bin/time -v singularity exec ${sif} cellChunk \
        -i ${outDir}/041.cellcut/${SN}.adjusted.cellbin.gef \
        -o ${outDir}/041.cellcut/
fi


# Run SAW saturation
export SINGULARITY_BIND=$outDir
echo `date` "=> saturation start ......"
export HDF5_USE_FILE_LOCKING=FALSE
bcStatStr=$(find ${outDir}/00.mapping -name \*stat | tr '\n' ',' | sed 's/.$//')
/usr/bin/time -v singularity exec ${sif} saturation \
    -i ${outDir}/02.count/${SN}_raw_barcode_gene_exp.txt \
    --tissue ${outDir}/04.tissuecut/${SN}.tissue.gef \
    -o ${outDir}/06.saturation \
    --bcstat ${bcStatStr} \
    --summary ${outDir}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat
## plot_200x200_saturation.png also named as ${SN}.saturation.bin200.png in some case.


# Run SAW report to generate HTML file
echo `date` "=> report generation start......"
export HDF5_USE_FILE_LOCKING=FALSE
export SINGULARITY_BIND=$outDir
out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)

bcStatStr=$(find ${outDir}/00.mapping -name \*stat | tr '\n' ',' | sed 's/.$//')
bcFinalOutStr=$(find ${outDir}/00.mapping -name \*.final.out | tr '\n' ',' | sed 's/.$//')
pipever=$(basename ${sif} .sif)

if [[ -n ${out_iprFile} ]] && [[ -e ${out_iprFile} ]] && [[ $doCell == 'Y' ]]; then
    /usr/bin/time -v singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${outDir}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${outDir}/04.tissuecut/tissuecut.stat \
        -n ${outDir}/02.count/${SN}.gef \
        -i ${outDir}/03.register/${SN}.rpi \
        -d ${outDir}/05.spatialcluster/${SN}.bin${binSize}_${resolution}.spatial.cluster.h5ad \
        -t ${outDir}/06.saturation/plot_200x200_saturation.png \
        -b ${outDir}/04.tissuecut/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${outDir}/04.tissuecut/tissue_fig/violin_200x200_MID_gene.png \
        -c ${outDir}/04.tissuecut/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        --bin20Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_20x20_MID_gene_counts.png \
        --bin20violin ${outDir}/04.tissuecut/tissue_fig/violin_20x20_MID_gene.png \
        --bin20MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_20x20_MID_gene_DNB.png \
        --bin50Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${outDir}/04.tissuecut/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${outDir}/04.tissuecut/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${outDir}/04.tissuecut/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --cellBinGef ${outDir}/041.cellcut/${SN}.adjusted.cellbin.gef \
        --cellCluster ${outDir}/051.cellcluster/${SN}.adjusted.cell.cluster.h5ad \
        --iprFile ${out_iprFile} \
        --pipelineVersion ${pipever} \
        -s ${SN} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${outDir}/07.report
elif [[ -n ${out_iprFile} ]] && [[ -e ${out_iprFile} ]] && [[ $doCell == 'N' ]]; then
    /usr/bin/time -v singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${outDir}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${outDir}/04.tissuecut/tissuecut.stat \
        -n ${outDir}/02.count/${SN}.gef \
        -i ${outDir}/03.register/${SN}.rpi \
        -d ${outDir}/05.spatialcluster/${SN}.bin${binSize}_${resolution}.spatial.cluster.h5ad \
        -t ${outDir}/06.saturation/plot_200x200_saturation.png \
        -b ${outDir}/04.tissuecut/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${outDir}/04.tissuecut/tissue_fig/violin_200x200_MID_gene.png \
        -c ${outDir}/04.tissuecut/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        --bin20Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_20x20_MID_gene_counts.png \
        --bin20violin ${outDir}/04.tissuecut/tissue_fig/violin_20x20_MID_gene.png \
        --bin20MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_20x20_MID_gene_DNB.png \
        --bin50Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${outDir}/04.tissuecut/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${outDir}/04.tissuecut/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${outDir}/04.tissuecut/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --iprFile ${out_iprFile} \
        --pipelineVersion ${pipever} \
        -s ${SN} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${outDir}/07.report
else
    /usr/bin/time -v singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${outDir}/02.count/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${outDir}/04.tissuecut/tissuecut.stat \
        -n ${outDir}/02.count/${SN}.gef \
        -d ${outDir}/05.spatialcluster/${SN}.bin${binSize}_${resolution}.spatial.cluster.h5ad \
        -t ${outDir}/06.saturation/plot_200x200_saturation.png \
        -b ${outDir}/04.tissuecut/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${outDir}/04.tissuecut/tissue_fig/violin_200x200_MID_gene.png \
        -c ${outDir}/04.tissuecut/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        --bin20Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_20x20_MID_gene_counts.png \
        --bin20violin ${outDir}/04.tissuecut/tissue_fig/violin_20x20_MID_gene.png \
        --bin20MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_20x20_MID_gene_DNB.png \
        --bin50Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${outDir}/04.tissuecut/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${outDir}/04.tissuecut/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${outDir}/04.tissuecut/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${outDir}/04.tissuecut/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${outDir}/04.tissuecut/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --pipelineVersion ${pipever} \
        -s ${SN} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${outDir}/07.report
fi



# Organize files for visualization required by StereoMap

## 02.count
if [[ -f ${outDir}/02.count/${SN}.gef ]]
then
    ln -s ${outDir}/02.count/${SN}.gef ${outDir}/visualization/${SN}.gef
fi

## 03.register
if [[ -f ${outDir}/03.register/fov_stitched_transformed.rpi ]] || [[ -f ${outDir}/03.register/${SN}.rpi ]] || [[ -f out_iprFile ]]
then
    ln -s ${outDir}/03.register/fov_stitched_transformed.rpi ${outDir}/visualization/fov_stitched_transformed.rpi
    ln -s ${outDir}/03.register/${SN}.rpi ${outDir}/visualization/${SN}.rpi
    ln -s $out_iprFile ${outDir}/visualization/$(basename "$out_iprFile")
fi

# 04.tissuecut
out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)
for i in `singularity exec ${sif} h5dump -n $out_iprFile | grep 'Labeling/' | grep -v 'Labeling/.*/Canvas'|grep -v 'Labeling/.*/Element'|awk '{print$2}'`;
do
label=`basename $i`
    if [[ -f ${outDir}/04.tissuecut/tissuecut_${label}/${SN}.${label}.label.gef ]]
    then
        ln -s ${outDir}/04.tissuecut/tissuecut_${label}/${SN}.${label}.label.gef ${outDir}/visualization/${SN}.${label}.label.gef
    fi
done

## 041.cellcut
if [[ -f ${outDir}/041.cellcut/${SN}.adjusted.cellbin.gef ]]
then
    ln -s ${outDir}/041.cellcut/${SN}.adjusted.cellbin.gef ${outDir}/visualization/${SN}.adjusted.cellbin.gef
fi

## 05.spatialcluster
if [[ -f ${outDir}/05.spatialcluster/${SN}.bin${binSize}_${resolution}.spatial.cluster.h5ad ]]
then
    ln -s ${outDir}/05.spatialcluster/${SN}.bin${binSize}_${resolution}.spatial.cluster.h5ad ${outDir}/visualization/${SN}.bin${binSize}_${resolution}.spatial.cluster.h5ad
fi

#out_iprFile=$(find ${outDir}/03.register -maxdepth 1 -name \*.ipr | head -1)
#for i in `singularity exec ${sif} h5dump -n $out_iprFile | grep 'Labeling/' | grep -v 'Labeling/.*/Canvas'|grep -v 'Labeling/.*/Element'|awk '{print$2}'`;
#do
#label=`basename $i`
#    if [[ -f ${outDir}/05.spatialcluster/${SN}.${label}.bin${binSize}_${resolution}.spatial.cluster.h5ad ]]
#    then
#        ln -s ${outDir}/05.spatialcluster/${SN}.${label}.bin${binSize}_${resolution}.spatial.cluster.h5ad ${outDir}/visualization/${SN}.${label}.bin${binSize}_${resolution}.spatial.cluster.h5ad
#    fi
#done

echo `date` " All done! "

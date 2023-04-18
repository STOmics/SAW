#!/bin/bash
set -e

if [[ $# -lt 12 ]];then
    echo "usage: sh $0 -genomeSize -splitCount -maskFile -fq1 -fq2 -refIndex -genomeFile -speciesName -tissueType -annotationFile -outDir -imageRecordFile -imageCompressedFile -doCellBin -threads -sif
    -genomeSize : genome size
    -splitCount : count of splited stereochip mask file, usually 16 for SE+Q4 fq data and 1 for PE+Q40 fq data
    -maskFile : stereochip mask file
    -fq1 : fastq file path of read1, if there are more than one fastq file, please separate them with comma, e.g:lane1_read_1.fq.gz,lane2_read_1.fq.gz
    -fq2 : fastq file path of read2, if there are more than one fastq file, please separate them with comma, not requested for 'SE+Q4' fastq data, e.g:lane1_read_2.fq.gz,lane2_read_2.fq.gz
    -refIndex : reference genome indexed folder, please build IT before SAW analysis run
    -speciesName : specie of the sample
    -tissueType : tissue type of the sample
    -annotationFile :  annotations file in gff or gtf format, the file must contain gene and exon annotations
    -outDir : output directory path
    -imageRecordFile : image file(*.ipr) generated by ImageQC software, not requested
    -imageCompressedFile : image file(*.tar.gz) generated by ImageQC software, not requested
    -doCellBin : [Y/N]
    -threads : the number of threads to be used in running the pipeline
    -sif : the file format of the visual software
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
        -outDir) outDir="$2"
            shift ;;
        -imageRecordFile) iprFile="$2"
            shift ;;
        -imageCompressedFile) imageTarFile="$2"
            shift ;;
        -doCellBin) doCell="$2"
            shift ;;
        -threads) threads="$2"
            shift ;;
        -sif) sif="$2"
            shift ;;
    esac
        shift
done


#software check
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
result_00mapping=${outDir}/00.mapping
result_01merge=${outDir}/01.merge
result_02count=${outDir}/02.count
result_04tissuecut=${outDir}/04.tissuecut
result_05spatialcluster=${outDir}/05.spatialcluster
result_06saturation=${outDir}/06.saturation
result_07report=${outDir}/07.report
result_visualization=${outDir}/visualization
arr_result=( $result_00mapping $result_01merge $result_02count $result_04tissuecut $result_05spatialcluster $result_06saturation $result_07report $result_visualization)
for each in "${arr_result[@]}";
do
    if [[ ! -d $each ]];then
        mkdir -p $each
    fi
done

if [[ $iprFile ]] && [[ $imageTarFile ]];then
    result_03register=${outDir}/03.register
    mkdir -p $result_03register
fi

if [[ $doCell == "Y" ]]; then
    result_041cellcut=${outDir}/041.cellcut
    result_051cellcluster=${outDir}/051.cellcluster
    mkdir -p $result_041cellcut
    mkdir -p $result_051cellcluster
fi


#Run splitMask or CIDCount for preparation
#  ulimit -c 100000000000
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
        echo $read1 | sed 's/,/\n/g' | grep _$i.fq.gz > ${result_00mapping}/mergeList/$a.${SN}.Q4_SE.fq.list
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
            --limitIObufferSize=280000000 \
            --outBAMsortingBinsN 50 \
            > ${result_00mapping}/${fqbase}_barcodeMap.stat

        starBam=${result_00mapping}/${fqbase}.Aligned.sortedByCoord.out.bam
        starBams[i]=$starBam
        bcStat[i]=${result_00mapping}/${fqbase}_barcodeMap.stat
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
        export SINGULARITY_BIND=$read1DIR$outDir,$maskDIR,$annoDIR,$refDIR
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
            --outBAMsortingBinsN 50 \
            > ${result_00mapping}/${fqbase}_barcodeMap.stat

        starBams[i]=${result_00mapping}/${fqbase}.Aligned.sortedByCoord.out.bam
        bcStat[i]=${result_00mapping}/${fqbase}_barcodeMap.stat
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


## Run SAW register to stitch microscope tile images to a panoramic image, perform tissue and cell (optional, depends on -doCellBin) segmentation, and register the panoramic image and the segmentated images with the gene expression matrix.
if [[ -f $imageTarFile ]] && [[ -f $iprFile ]]  && [[ $doCell == "Y" ]]; then
    # Run SAW register (stitch, tissue segmentation, cell segmentation) + SAW imageTools
    echo `date` "=> image processing and registration start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    imgTarDIR=$(dirname $imageTarFile)
    iprDIR=$(dirname $iprFile)
    export SINGULARITY_BIND=$outDir,$imgTarDIR,$iprDIR
    /usr/bin/time -v singularity exec ${sif} register \
            -i ${imageTarFile} \
            -c ${iprFile} \
            -v ${geneExp} \
            -o ${result_03register}
    out_iprFile=$(find ${result_03register} -maxdepth 1 -name *.ipr | head -1)
    /usr/bin/time -v singularity exec ${sif} imageTools ipr2img \
            -i ${imageTarFile} \
            -c ${out_iprFile} \
            -d tissue cell \
            -r True \
            -o ${result_03register}
    registerTif=$(find ${result_03register} -maxdepth 1 -name *fov_stitched_transformed.tif)
    regTifStr=$(echo $registerTif | tr ' ' ',')
    regGroup=$(find ${result_03register} -maxdepth 1 -name '*fov_stitched_transformed.tif' -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} +)
    regGroupStr=$(echo $regGroup | sed 's/ \|$/\/Image,/g' | sed 's/.$//')
    echo $regTifStr
    echo $regGroupStr
    /usr/bin/time -v singularity exec ${sif} imageTools img2rpi \
        -i ${regTifStr} \
        -g ${regGroupStr} \
        -b 1 10 50 100 \
        -o ${result_03register}/fov_stitched_transformed.rpi
elif [[ -f $imageTarFile ]] && [[ -f $iprFile ]]  && [[ $doCell == "N" ]]; then
    # Run SAW rapidRegister (stitch, tissue segmentation) + SAW imageTools
    echo `date` "=> image processing and registration start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    imgTarDIR=$(dirname $imageTarFile)
    iprDIR=$(dirname $iprFile)

    export SINGULARITY_BIND=$outDir,$imgTarDIR,$iprDIR
    /usr/bin/time -v singularity exec ${sif} rapidRegister \
            -i ${imageTarFile} \
            -c ${iprFile} \
            -v ${geneExp} \
            -o ${result_03register}
    out_iprFile=$(find ${result_03register} -maxdepth 1 -name *.ipr | head -1)
    /usr/bin/time -v singularity exec ${sif} imageTools ipr2img \
            -i ${imageTarFile} \
            -c ${out_iprFile} \
            -d tissue \
            -r True \
            -o ${result_03register}
    registerTif=$(find ${result_03register} -maxdepth 1 -name *fov_stitched_transformed.tif)
    regTifStr=$(echo $registerTif | tr ' ' ',')
    regGroup=$(find ${result_03register} -maxdepth 1 -name '*fov_stitched_transformed.tif' -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} +)
    regGroupStr=$(echo $regGroup | sed 's/ \|$/\/Image,/g' | sed 's/.$//')
    echo $regTifStr
    echo $regGroupStr
    /usr/bin/time -v singularity exec ${sif} imageTools img2rpi \
        -i ${regTifStr} \
        -g ${regGroupStr} \
        -b 1 10 50 100 \
        -o ${result_03register}/fov_stitched_transformed.rpi
fi


# Run SAW tissueCut
if [[ -f $imageTarFile ]] && [[ -f $iprFile ]]; then
    # Run tissueCut to get the spatial gene expression profile of the tissue-covered region
    nucleusLayer=$(find ${result_03register} -maxdepth 1 -name *fov_stitched_transformed.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} + | grep -v IF)
    tissueMaskFile=$(find ${result_03register} -maxdepth 1 -name ${nucleusLayer}*_tissue_cut.tif)
    echo `date` "=> tissueCut start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    export SINGULARITY_BIND=$outDir
    /usr/bin/time -v singularity exec ${sif} tissueCut \
        -i ${geneExp} \
        --dnbfile ${barcodeReadsCounts} \
        -s ${tissueMaskFile} \
        --sn ${SN} \
        -O Transcriptomics \
        -d -t 8 \
        -o ${result_04tissuecut}
else
    # Run SAW tissueCut based on the gene expression matrix directly
    export SINGULARITY_BIND=$outDir,$annoDIR,$refDIR
    echo `date` "=> there is no image, tissueCut based on the gene expression matrix start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    /usr/bin/time -v singularity exec ${sif} tissueCut \
        -i ${geneExp} \
        --dnbfile ${barcodeReadsCounts} \
        --sn ${SN} \
        -O Transcriptomics \
        -d -t 8 \
        -o ${result_04tissuecut}
fi


## Complete raw GEF to visual GEF
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} cellCut bgef \
    -i ${geneExp} \
    -o ${result_04tissuecut}/${SN}.gef \
    -O Transcriptomics
## Convert GEF to GEM [optional]
/usr/bin/time -v singularity exec ${sif} cellCut view \
    -s ${SN} \
    -i ${result_04tissuecut}/${SN}.gef \
    -o ${result_04tissuecut}/${SN}.gem
gzip ${result_04tissuecut}/${SN}.gem
/usr/bin/time -v singularity exec ${sif} cellCut view \
    -s ${SN} \
    -i ${result_04tissuecut}/${SN}.tissue.gef \
    -o ${result_04tissuecut}/${SN}.tissue.gem
gzip ${result_04tissuecut}/${SN}.tissue.gem


# Run spatialCluster
binSize=200
export SINGULARITY_BIND=$outDir
echo `date` "=> spatialCluster start......."
export HDF5_USE_FILE_LOCKING=FALSE
mkdir -p ${outDir}/tmp
export NUMBA_CACHE_DIR=${outDir}/tmp
export MPLCONFIGDIR=${outDir}/tmp

/usr/bin/time -v singularity exec ${sif} spatialCluster \
    -i ${result_04tissuecut}/${SN}.tissue.gef \
    -o ${result_05spatialcluster}/${SN}.spatial.cluster.h5ad \
    -s $binSize


# Run cellCut and cellCluster
export SINGULARITY_BIND=$outDir
if [[ $doCell == 'Y' ]]; then
    echo `date` "=> cellCut start......."
    export HDF5_USE_FILE_LOCKING=FALSE
    nucleusLayer=$(find ${result_03register} -maxdepth 1 -name *fov_stitched_transformed.tif -exec sh -c 'for f do basename -- "$f" _fov_stitched_transformed.tif;done' sh {} + | grep -v IF)
    nucleusMask=$(find ${result_03register} -maxdepth 1 -name ${nucleusLayer}*_mask.tif)
    /usr/bin/time -v singularity exec ${sif} cellCut cgef \
        -i ${geneExp} \
        -m ${nucleusMask} \
        -o ${result_041cellcut}/${SN}.cellbin.gef
    ## convert cellbin GEF to GEM
    # /usr/bin/time -v singularity exec ${sif} cellCut view \
    #     -i ${result_041cellcut}/${SN}.cellbin.gef \
    #     -d ${result_04tissuecut}/${SN}.gef \
    #     -o ${result_041cellcut}/${SN}.cellbin.gem \
    #     -s ${SN}
    echo `date` "=> cellCluster start......."
    mkdir -p ${outDir}/tmp
    export NUMBA_CACHE_DIR=${outDir}/tmp
    export MPLCONFIGDIR=${outDir}/tmp
    /usr/bin/time -v singularity exec ${sif} cellCluster \
        -i ${result_041cellcut}/${SN}.cellbin.gef \
        -o ${result_051cellcluster}/${SN}.cell.cluster.h5ad
fi


# Run saturation
export SINGULARITY_BIND=$outDir
echo `date` "=> saturation start ......"
export HDF5_USE_FILE_LOCKING=FALSE
/usr/bin/time -v singularity exec ${sif} saturation \
    -i ${saturationFile} \
    --tissue ${result_04tissuecut}/${SN}.tissue.gef \
    -o ${result_06saturation} \
    --bcstat ${bcStatStr} \
    --summary ${result_02count}/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat
## plot_200x200_saturation.png also named as ${SN}.saturation.bin200.png in some case.


# Run report to generate HTML report
echo `date` "=> report generation start......"
export HDF5_USE_FILE_LOCKING=FALSE
export SINGULARITY_BIND=$outDir
out_iprFile=$(find ${result_03register} -maxdepth 1 -name *.ipr | head -1)

if [[ -n ${out_iprFile} ]] && [[ -e ${out_iprFile} ]] && [[ $doCell == 'Y' ]]; then
    /usr/bin/time -v singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${result_02count}/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${result_04tissuecut}/tissuecut.stat \
        -n ${result_04tissuecut}/${SN}.gef \
        -i ${result_03register}/${SN}.rpi \
        -d ${result_05spatialcluster}/${SN}.spatial.cluster.h5ad \
        -t ${result_06saturation}/plot_200x200_saturation.png \
        -b ${result_04tissuecut}/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${result_04tissuecut}/tissue_fig/violin_200x200_MID_gene.png \
        -c ${result_04tissuecut}/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        --bin20Saturation ${result_04tissuecut}/tissue_fig/scatter_20x20_MID_gene_counts.png \
        --bin20violin ${result_04tissuecut}/tissue_fig/violin_20x20_MID_gene.png \
        --bin20MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_20x20_MID_gene_DNB.png \
        --bin50Saturation ${result_04tissuecut}/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${result_04tissuecut}/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${result_04tissuecut}/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${result_04tissuecut}/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${result_04tissuecut}/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${result_04tissuecut}/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --cellBinGef ${result_041cellcut}/${SN}.cellbin.gef \
        --cellCluster ${result_051cellcluster}/${SN}.cell.cluster.h5ad \
        --iprFile ${out_iprFile} \
        --pipelineVersion saw_v6.0.0 \
        -s ${SN} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${result_07report}
elif [[ -n ${out_iprFile} ]] && [[ -e ${out_iprFile} ]] && [[ $doCell == 'N' ]]; then
    /usr/bin/time -v singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${result_02count}/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${result_04tissuecut}/tissuecut.stat \
        -n ${result_04tissuecut}/${SN}.gef \
        -i ${result_03register}/${SN}.rpi \
        -d ${result_05spatialcluster}/${SN}.spatial.cluster.h5ad \
        -t ${result_06saturation}/plot_200x200_saturation.png \
        -b ${result_04tissuecut}/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${result_04tissuecut}/tissue_fig/violin_200x200_MID_gene.png \
        -c ${result_04tissuecut}/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        --bin20Saturation ${result_04tissuecut}/tissue_fig/scatter_20x20_MID_gene_counts.png \
        --bin20violin ${result_04tissuecut}/tissue_fig/violin_20x20_MID_gene.png \
        --bin20MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_20x20_MID_gene_DNB.png \
        --bin50Saturation ${result_04tissuecut}/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${result_04tissuecut}/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${result_04tissuecut}/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${result_04tissuecut}/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${result_04tissuecut}/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${result_04tissuecut}/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --iprFile ${out_iprFile} \
        --pipelineVersion saw_v6.0.0 \
        -s ${SN} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${result_07report}
else
    /usr/bin/time -v singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${result_02count}/${SN}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${result_04tissuecut}/tissuecut.stat \
        -n ${result_04tissuecut}/${SN}.gef \
        -d ${result_05spatialcluster}/${SN}.spatial.cluster.h5ad \
        -t ${result_06saturation}/plot_200x200_saturation.png \
        -b ${result_04tissuecut}/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${result_04tissuecut}/tissue_fig/violin_200x200_MID_gene.png \
        -c ${result_04tissuecut}/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        --bin20Saturation ${result_04tissuecut}/tissue_fig/scatter_20x20_MID_gene_counts.png \
        --bin20violin ${result_04tissuecut}/tissue_fig/violin_20x20_MID_gene.png \
        --bin20MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_20x20_MID_gene_DNB.png \
        --bin50Saturation ${result_04tissuecut}/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${result_04tissuecut}/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${result_04tissuecut}/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${result_04tissuecut}/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${result_04tissuecut}/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${result_04tissuecut}/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --pipelineVersion saw_v6.0.0 \
        -s ${SN} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${result_07report}
fi

## Organize files required by StereoMap
# 03.register
printf $result_03register/fov_stitched_transformed.rpi;
if [[ -f $result_03register/fov_stitched_transformed.rpi ]] || [[ -f $result_03register/${SN}.rpi ]] || [[ -f out_iprFile ]]
then
    ln -s $result_03register/fov_stitched_transformed.rpi $result_visualization/fov_stitched_transformed.rpi
    ln -s $result_03register/${SN}.rpi $result_visualization/${SN}.rpi
    ln -s $out_iprFile $result_visualization/$(basename "$out_iprFile")
fi

# 04.tissuecut
if [[ -f $result_04tissuecut/${SN}.gef ]]
then
    ln -s $result_04tissuecut/${SN}.gef $result_visualization/${SN}.gef
fi

# 041.cellcut
if [[ -f $result_041cellcut/${SN}.cellbin.gef ]]
then
    ln -s $result_041cellcut/${SN}.cellbin.gef $result_visualization/${SN}.cellbin.gef
fi

echo `date` " All done! "

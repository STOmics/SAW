#####saw_v5.1.3#####
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
        -annotationFile) AFile="$2"
            shift ;;
        -outDir) outDir="$2"
            shift ;;
        -imageRecordFile) imageRF="$2"
            shift ;;
        -imageCompressedFile) imageCF="$2"
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
    echo `date` " singularity check: pass and path is ${singularityPath}"
else
    echo `date` " singularity check: singularity does not exits, please verify that you have installed singularity and exported it to system PATH variable"
    exit
fi

if [[ -n $sif ]]
then
    echo `date` "  visualSif check: file exists and path is ${sif}"
else
    echo `date` "  visualSif check: file does not exists, please verify that your visualSif file in the current directory or the path given by the option -s is correct."
fi


if [[ ! -d $outDir ]];then
    mkdir -p $outDir
fi

#basic information get
maskname=$(basename $maskFile)
SNid=${maskname%%.*}

maskDIR=$(dirname $maskFile)
annoDIR=$(dirname $AFile)
refDIR=$(dirname $GDir)
if [[ $imageRF ]] && [[ $imageCF ]];then 
    imgRFDIR=$(dirname $imageRF)
    imgCFDIR=$(dirname $imageCF)
fi
#create result path
##result path
result_00mapping=${outDir}/00.mapping
result_01merge=${outDir}/01.merge
result_02count=${outDir}/02.count
result_03register=${outDir}/03.register
result_04tissuecut=${outDir}/04.tissuecut
result_041cellcut=${outDir}/041.cellcut
result_05spatialcluster=${outDir}/05.spatialcluster
result_051cellcluster=${outDir}/051.cellcluster
result_06saturation=${outDir}/06.saturation
result_07report=${outDir}/07.report
result_sn=${outDir}/${SNid}_result
arr_result=( $result_00mapping $result_01merge $result_02count $result_03register $result_04tissuecut $result_041cellcut $result_05spatialcluster $result_051cellcluster $result_06saturation $result_07report $result_sn)
for each in "${arr_result[@]}";
do
    if [[ ! -d $each ]];then
        mkdir -p $each
    fi
done


### STEP1: preparation step => calculate the memory of reference  
##  ulimit -c 100000000000


#barcode mapping and star alignment.
echo `date` "=> barcode mapping, adapter filter and RNA alignment start......"
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
    singularity exec ${sif} splitMask \
        ${maskFile} ${result_00mapping}/splitBin $threads $splitCnt 2_25 
    for ((i=1;i<=$splitCnt;i++)); do
        if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
        echo $read1 | sed 's/,/\n/g' | grep _$i.fq.gz > ${result_00mapping}/mergeList/$a.${SNid}.Q4_SE.fq.list
    done    
else
    fqType="PE"
    read2List=(`echo $read2 | tr ',' ' '`)
    fqNumber=`echo ${#read1List[@]}`
    export SINGULARITY_BIND=$outDir,$maskDIR,$annoDIR,$refDIR
    singularity exec ${sif} CIDCount \
        -i ${maskFile} \
        -s ${refName} \
        -g ${GSize} > ${result_00mapping}/CIDCount
fi
echo "your fastq is in ${fqType} data format"

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
        read1DIR=$(dirname ${read1List[i]})
        read2DIR=$(dirname ${read2List[i]})
        export SINGULARITY_BIND=$read1DIR,$read2DIR,$outDir,$maskDIR,$annoDIR,$refDIR
        singularity exec ${sif} mapping \
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
        singularity exec ${sif} CIDCount \
            -i $(ls ${result_00mapping}/splitBin/${a}.${SNid}.barcodeToPos.bin) \
            -s ${refName} \
            -g ${GSize} > ${result_00mapping}/CIDCount
        fqbase=$a.${SNid}.Q4_SE
        bcPara=${result_00mapping}/${fqbase}.bcPara
        barcodeReadsCount=${result_00mapping}/${fqbase}.barcodeReadsCount.txt
        echo  " ~~~ mapping SE fastqData - $a ~~~"
        echo "in=$(ls ${result_00mapping}/splitBin/${a}.${SNid}.barcodeToPos.bin)" > $bcPara
        read1List=${result_00mapping}/mergeList/$a.${SNid}.Q4_SE.fq.list
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
        echo "useF14" >> $bcPara
        echo "polyAnum=15" >> $bcPara
        echo "mismatchInPolyA=2" >> $bcPara
        echo "bcNum=`head -1 ${result_00mapping}/CIDCount`" >> $bcPara
        export SINGULARITY_BIND=$read1DIR$outDir,$maskDIR,$annoDIR,$refDIR
        singularity exec ${sif} mapping \
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

if [[ $(echo ${#bcReadsCounts[*]}) == '1' ]];then
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

#merge barcode reads count file
echo `date` "=> merge barcode reads count tables start......"
barcodeReadsCounts=${result_01merge}/${SNid}.merge.barcodeReadsCount.txt
export SINGULARITY_BIND=$outDir,$maskDIR
if [[ $fqType == 'PE' ]] && [[ $(echo ${#bcReadsCounts[*]}) == '1' ]]; then
    cp $bcReadsCountsStr $barcodeReadsCounts
else
    singularity exec ${sif} merge \
        ${maskFile} \
        $bcReadsCountsStr \
        $barcodeReadsCounts
fi

#annotation and deduplication
echo `date` "=> annotation and deduplication start......"
geneExp=${result_02count}/${SNid}.raw.gef
saturationFile=${result_02count}/${SNid}_raw_barcode_gene_exp.txt
export SINGULARITY_BIND=$outDir,$annoDIR,$refDIR
singularity exec ${sif} count \
    -i ${starBamsStr} \
    -o ${result_02count}/${SNid}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
    -a ${AFile} \
    -s ${result_02count}/${SNid}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
    -e ${geneExp} \
    --sat_file ${saturationFile} \
    --umi_on \
    --save_lq \
    --save_dup \
    --sn ${SNid} \
    -c ${threads} \
    -m 128 

if [[ -f $imageCF ]] && [[ -f $imageRF ]];then
    #firstly do registration, then cut the gene expression matrix based on the repistration result
    echo `date` "=> registration and tissueCut start......."
    imageDirCF=$(dirname $imageCF)
    imageDirRF=$(dirname $imageRF)
    export SINGULARITY_BIND=$outDir,$imageDirCF,$imageDirRF
    singularity exec ${sif} register \
        -i ${imageCF} \
        -c ${imageRF} \
        -v ${geneExp} \
        -o ${result_03register} 
    imageRF_out=$(find ${result_03register} -maxdepth 1 -name *.ipr | head -1)
    singularity exec ${sif} ipr2img \
        -i ${imageCF} \
        -c ${imageRF_out} \
        -d tissue cell \
        -r True \
        -o ${result_03register} 

    registTiff=$(find ${result_03register} -maxdepth 1 -name *_regist.tif | head -1)


    ## organize your outputs (optional)
    mkdir -p $result_sn/IMAGE
    if [[ ! -f $result_sn/IMAGE/${SNid}.transform.attrs.json ]]
    then
            ln ${result_03register}/attrs.json $result_sn/IMAGE/${SNid}.transform.attrs.json
    fi

    if [[ ! -f $result_sn/IMAGE/${SNid}.transform.thumbnail.png ]]
    then
            ln ${result_03register}/transform_thumb.png $result_sn/IMAGE/${SNid}.transform.thumbnail.png
    fi


    echo `date` "=>   tissuecut start......."
    export SINGULARITY_BIND=$outDir
    singularity exec ${sif} tissueCut \
        --dnbfile ${barcodeReadsCounts} \
        -i ${geneExp} \
        -o ${result_04tissuecut} \
        -s ${result_03register} \
        -t tissue \
        --platform T10 \
        --snId ${SNid} \
        --omics=Transcriptomics \
        --develop 

    singularity exec ${sif} cellCut \
        view -s ${SNid} -i ${result_04tissuecut}/${SNid}.gef -o ${result_04tissuecut}/${SNid}.gem
    gzip ${result_04tissuecut}/${SNid}.gem

    singularity exec ${sif} cellCut \
        view -s ${SNid} -i ${result_04tissuecut}/${SNid}.tissue.gef -o ${result_04tissuecut}/${SNid}.tissue.gem
    gzip ${result_04tissuecut}/${SNid}.tissue.gem

    ## organize your outputs (optional)
    if [[ ! -f $result_sn/${SNid}.ssDNA.rpi ]] || [[ ! -f $result_sn/${SNid}.gef ]] || [[ ! -f $result_sn/${SNid}.thumbnail.png ]] || [[ ! -f $result_sn/${SNid}.tissue.gef ]]
    then
         ln ${result_04tissuecut}/${SNid}.gef $result_sn/${SNid}.gef
         ln ${result_04tissuecut}/${SNid}.gem.gz $result_sn/${SNid}.gem.gz
         ln ${result_04tissuecut}/tissue_fig/${SNid}.ssDNA.rpi $result_sn/${SNid}.ssDNA.rpi
         ln ${result_04tissuecut}/dnb_merge/bin200.png $result_sn/${SNid}.thumbnail.png
         ln ${result_04tissuecut}/${SNid}.tissue.gef $result_sn/${SNid}.tissue.gef
         ln ${result_04tissuecut}/${SNid}.tissue.gem.gz $result_sn/${SNid}.tissue.gem.gz
    fi

else
    #cut the gene expression matrix directly
    export SINGULARITY_BIND=$outDir,$annoDIR,$refDIR
    echo `date` " there is no image, tissueCut start......."
    singularity exec ${sif} tissueCut \
        --dnbfile ${barcodeReadsCounts} \
        -i ${geneExp} \
        -o ${result_04tissuecut} \
        -t tissue \
        --platform T10 \
        --snId ${SNid} \
        --omics=Transcriptomics \
        --develop
fi

## organize your outputs (optional)
if [[ ! -f $result_sn/${SNid}.gef ]] || [[ ! -f $result_sn/${SNid}.thumbnail.png ]] || [[ ! -f $result_sn/${SNid}.tissue.gef ]]
then
        ln ${result_04tissuecut}/${SNid}.gef $result_sn/${SNid}.gef
        ln ${result_04tissuecut}/dnb_merge/bin200.png $result_sn/${SNid}.thumbnail.png
        ln ${result_04tissuecut}/${SNid}.tissue.gef $result_sn/${SNid}.tissue.gef
fi

#spatialCluster
binSize=200
export SINGULARITY_BIND=$outDir
echo `date` "=>   spatialCluster start......."
singularity exec ${sif} spatialCluster \
    -i ${result_04tissuecut}/${SNid}.tissue.gef \
    -o ${result_05spatialcluster}/${SNid}.spatial.cluster.h5ad \
    -s $binSize 

## organize your outputs (optional)
ln ${result_05spatialcluster}/${SNid}.spatial.cluster.h5ad $result_sn/${SNid}.spatial.cluster.h5ad

#cellCut
export SINGULARITY_BIND=$outDir
if [[ $doCell == 'Y' ]]; then
    echo `date` "=>   cellCut start......."
    singularity exec ${sif} cellCut cgef \
        -i ${geneExp} \
        -m ${result_03register}/${SNid}*_mask.tif \
        -o ${result_041cellcut}/${SNid}.cellbin.gef
    singularity exec ${sif} cellCut view \
        -i ${result_041cellcut}/${SNid}.cellbin.gef \
        -d ${result_04tissuecut}/${SNid}.gef \
        -o ${result_041cellcut}/${SNid}.cellbin.gem \
        -s ${SNid}
    echo `date` "=>   cellCluster start......."
    singularity exec ${sif} cellCluster \
        -i ${result_041cellcut}/${SNid}.cellbin.gef \
        -o ${result_051cellcluster}/${SNid}.cell.cluster.h5ad
    ## organize your outputs (optional)
    ln ${result_051cellcluster}/${SNid}.cell.cluster.h5ad $result_sn/${SNid}.cell.cluster.h5ad
fi



#saturation
export SINGULARITY_BIND=$outDir
echo `date` "=> saturation start ......"
singularity exec ${sif} saturation \
    -i ${saturationFile} \
    --tissue ${result_04tissuecut}/${SNid}.tissue.gef \
    -o ${result_06saturation} \
    --bcstat ${bcStatStr} \
    --summary ${result_02count}/${SNid}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat 
## organize your outputs (optional)
if [[ ! -f $result_sn/${SNid}.saturation.bin200.png ]]
then
    ln ${result_06saturation}/plot_200x200_saturation.png $result_sn/${SNid}.saturation.bin200.png
fi


#generate report file in json format
echo `date` "=> report generation start......"
echo "${result_04tissuecut}/tissue_fig/${SNid}.ssDNA.rpi    ${result_04tissuecut}/tissue_fig/${SNid}.ssDNA.rpi"
export SINGULARITY_BIND=$outDir


packageFile=/opt/saw_v5.1.3_software/pipeline/report/plotly_package.txt
iprFile=`ls ${result_03register}/${SNid}*.ipr`
export SINGULARITY_BIND=$outDir

if [[ -n ${iprFile} ]] && [[ -e ${iprFile} ]] && [[ $doCell == 'Y' ]]
then
  echo 333
    singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${result_02count}/${SNid}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${result_04tissuecut}/tissuecut.stat \
        -n ${result_04tissuecut}/${SNid}.gef \
        -i ${result_04tissuecut}/tissue_fig/${SNid}.ssDNA.rpi \
        -d ${result_05spatialcluster}/${SNid}.spatial.cluster.h5ad \
        -t ${result_06saturation}/plot_200x200_saturation.png \
        -b ${result_04tissuecut}/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${result_04tissuecut}/tissue_fig/violin_200x200_MID_gene.png \
        -c ${result_04tissuecut}/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        -p ${packageFile} \
        --bin1Saturation ${result_06saturation}/plot_1x1_saturation.png \
        --bin50Saturation ${result_04tissuecut}/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${result_04tissuecut}/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${result_04tissuecut}/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${result_04tissuecut}/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${result_04tissuecut}/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${result_04tissuecut}/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --cellBinGef ${result_041cellcut}/${SNid}.cellbin.gef \
        --cellCluster ${result_051cellcluster}/${SNid}.cell.cluster.h5ad \
        --iprFile ${iprFile} \
        --pipelineVersion saw_v5.1.3 \
        -s ${SNid} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${result_07report} \
        --logo /opt/saw_v5.1.3_software/pipeline/report/logo.png \
        -r standard_version

elif [[ -n ${iprFile} ]] && [[ -e ${iprFile} ]]
then
    singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${result_02count}/${SNid}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${result_04tissuecut}/tissuecut.stat \
        -n ${result_04tissuecut}/${SNid}.gef \
        -i ${result_04tissuecut}/tissue_fig/${SNid}.ssDNA.rpi \
        -d ${result_05spatialcluster}/${SNid}.spatial.cluster.h5ad \
        -t ${result_06saturation}/plot_200x200_saturation.png \
        -b ${result_04tissuecut}/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${result_04tissuecut}/tissue_fig/violin_200x200_MID_gene.png \
        -c ${result_04tissuecut}/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        -p ${packageFile} \
        --bin1Saturation ${result_06saturation}/plot_1x1_saturation.png \
        --bin50Saturation ${result_04tissuecut}/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${result_04tissuecut}/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${result_04tissuecut}/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${result_04tissuecut}/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${result_04tissuecut}/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${result_04tissuecut}/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --iprFile ${iprFile} \
        --pipelineVersion saw_v5.1.3 \
        -s ${SNid} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${result_07report} \
        --logo /opt/saw_v5.1.3_software/pipeline/report/logo.png \
        -r standard_version
else
    singularity exec ${sif} report \
        -m ${bcStatStr} \
        -a ${bcFinalOutStr} \
        -g ${result_02count}/${SNid}.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
        -l ${result_04tissuecut}/tissuecut.stat \
        -n ${result_04tissuecut}/${SNid}.gef \
        -d ${result_05spatialcluster}/${SNid}.spatial.cluster.h5ad \
        -t ${result_06saturation}/plot_200x200_saturation.png \
        -b ${result_04tissuecut}/tissue_fig/scatter_200x200_MID_gene_counts.png \
        -v ${result_04tissuecut}/tissue_fig/violin_200x200_MID_gene.png \
        -c ${result_04tissuecut}/tissue_fig/statistic_200x200_MID_gene_DNB.png \
        -p ${packageFile} \
        --bin1Saturation ${result_06saturation}/plot_1x1_saturation.png \
        --bin50Saturation ${result_04tissuecut}/tissue_fig/scatter_50x50_MID_gene_counts.png \
        --bin50violin ${result_04tissuecut}/tissue_fig/violin_50x50_MID_gene.png \
        --bin50MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_50x50_MID_gene_DNB.png \
        --bin100Saturation ${result_04tissuecut}/tissue_fig/scatter_100x100_MID_gene_counts.png \
        --bin100violin ${result_04tissuecut}/tissue_fig/violin_100x100_MID_gene.png \
        --bin100MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_100x100_MID_gene_DNB.png \
        --bin150Saturation ${result_04tissuecut}/tissue_fig/scatter_150x150_MID_gene_counts.png \
        --bin150violin ${result_04tissuecut}/tissue_fig/violin_150x150_MID_gene.png \
        --bin150MIDGeneDNB ${result_04tissuecut}/tissue_fig/statistic_150x150_MID_gene_DNB.png \
        --pipelineVersion saw_v5.1.3 \
        -s ${SNid} \
        --species ${refName} \
        --tissue ${tissueType} \
        --reference ${refName} \
        -o ${result_07report} \
        --logo /opt/saw_v5.1.3_software/pipeline/report/logo.png \
        -r standard_version

fi
## organize your outputs (optional)
if [[ ! -f $result_sn/${SNid}.statistics.json ]] || [[ ! -f $result_sn/${SNid}.report.html ]]
then
    ln ${result_07report}/${SNid}.statistics.json $result_sn/${SNid}.statistics.json
    ln ${result_07report}/${SNid}.report.html $result_sn/${SNid}.report.html
fi

echo `date` " all finish "

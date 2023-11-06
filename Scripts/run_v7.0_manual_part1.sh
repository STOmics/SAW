# The first step is to get gene expression matrix (SN.gef) and visualizing file (fov_stitched*.rpi) for manual registration.

ulimit -n 10240
ulimit -v 33170449147

sif=/path/to/sif

dataDir=/path/to/data
imgDir=/path/to/image/files
outDir=/path/to/output

export SINGULARITY_BIND=${dataDir},${outDir}
export PATH=/path/to/singularity/bin:$PATH
bash stereoPipeline_v7.0.0_manual_part1.sh \
    -splitCount <split_num> \
    -maskFile ${dataDir}/path/to/chipMask \
    -fq1 ${dataDir}/path/to/read1 \
    -fq2 ${dataDir}/path/to/read2 \
    -speciesName <speciesName> \
    -tissueType <tissueName> \
    -imageRecordFile ${imgDir}/<SN_date_time_version>.ipr \
    -imageCompressedFile ${imgDir}/<SN_date_time_version>.tar.gz \
    -refIndex  ${dataDir}/path/to/reference\
    -annotationFile ${dataDir}/path/to/gtf/or/gff \
    -rRNAremove <Y or N> \  # whether to remove rRNA
    -sif ${sif} \
    -threads <threads_num> \
    -outDir ${outDir}
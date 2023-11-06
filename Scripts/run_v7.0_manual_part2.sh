# The second step is to pass manually processed information (.json) to IPR (.ipr) file, and then output subsequent analysis results.

ulimit -n 10240
ulimit -v 33170449147

sif=/path/to/sif

dataDir=/path/to/data	# notice that dataDir in step2 should be the outDir from step1, to get basic information like statistics
iprDir=/path/to/ipr  # when QCPass, changed the directory to processed ipr from register module; when QCFail, set the directory to original ipr
tarDir=/path/to/tar.gz
outDir=/path/to/output
registJsonDir=/path/to/<date_time>.regist.json


export SINGULARITY_BIND=${dataDir},${iprDir},${outDir},${registJsonDir}
export PATH=/path/to/singularity/bin:$PATH
bash stereoPipeline_v7.0.0_manual_part2.sh \
    -SN <SN> \
    -dataDir ${dataDir} \
    -registJson ${registJsonDir} \
    -speciesName <speciesName> \
    -tissueType <tissueName> \
    -imageRecordFile ${iprDir} \
    -imageCompressedFile ${tarDir} \
    -sif ${sif} \
	-doCellBin <Y or N> \
    -threads <threads_num> \
    -outDir ${outDir}
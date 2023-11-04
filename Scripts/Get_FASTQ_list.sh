FASTQDirList=$1
outDir=$2
SN=$3
splitCnt=16

mkdir -p ${outDir}/00.mapping/mergeList
for ((i=1;i<=$splitCnt;i++)); do
    if [[ $(echo ${#i}) == '1' ]];then a=0$i; else a=$i;fi
    while IFS= read -r line
    do
      ls $line/* | grep _$i.fq.gz >> ${outDir}/00.mapping/mergeList/$a.${SN}.Q4.fq.list
    done < "$FASTQDirList"
done
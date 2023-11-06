# TEST DATA
##  Demo data basic information

|  | SN | ChipSize | Reference | Tissue Type |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| Demo 1 | SS200000135TL_D1 | 1\*1 | mouse | brain |
| Demo 2 | SS200000154TR_F5 | 1\*1 | mouse | tongue |
| Demo 3 | SS200000464BL_C4 | 1\*1 | mouse | heart |
| Demo 4 | SS200000059_NC | 1\*1 | mouse | testis |
| Demo 5 | D02070C3D3 | 1\*2 | mouse | embryo |
| Demo 6 | FP200009107_E414 | 0.5\*0.5 | mouse | olfactory bulb |
| Demo 7 | C02533C1 | 1\*1 | mouse | kidney |


##  Download data from Raysync
**Link: http://116.6.21.110:8090/share/21bb9df9-e6c5-47c5-9aa8-29f2d23a6df4**

![demo_data.png](demo_data.png)

| Demo Data Directory | SAW Version |
| ----------- | ----------- |
| SS200000135TL_D1_v4_brain | <= V4.1.0  |
| SS200000135TL_D1_v5_brain | \>= V5.1.3 |
| SS200000154TR_F5_v5_tongue | \>= V5.1.3 |
| SS200000464BL_C4_v5_heart | \>= V5.1.3 |
| SS200000059_NC_v6_testis | \>= V6.0.0 |
| D02070C3D3_v6.1_embryo | \>= V6.1.0 |
| FP200009107_E414_v6.1_olfactory_bulb | \>= V6.1.0 |
| C02533C1_v7.0_kidney | \>= V7.0.0 |

## Raw Data Directory Structure
Here we take `SS200000135TL_D1_v5` as an example. 
```
# v1
cd /data/dataManagement/reference
tree ./
.
└── specieName
    ├── STAR_SJ100
    │   ├── Genome
    │   ├── SA
    │   ├── SAindex
    │   ├── chrLength.txt
    │   ├── chrName.txt
    │   ├── chrNameLength.txt
    │   ├── chrStart.txt
    │   ├── exonGeTrInfo.tab
    │   ├── exonInfo.tab
    │   ├── geneInfo.tab
    │   ├── genomeParameters.txt
    │   ├── sjdbInfo.txt
    │   ├── sjdbList.fromGTF.out.tab
    │   ├── sjdbList.out.tab
    │   └── transcriptInfo.tab
    ├── genes
    │   └── genes.gtf
    └── genome
        └── genome.fa

4 directories, 17 files


# v2
cd /data/dataManagement/reference
tree ./
.
└── specieName
    ├── STAR_SJ100
    |   ├── Genome 
    |   ├── SA
    |   ├── SAindex
    |   ├── SAindexAux
    |   ├── chrLength.txt
    |   ├── chrName.txt
    |   ├── chrNameLength.txt
    |   ├── chrStart.txt
    |   ├── exonGeTrInfo.tab
    |   ├── exonInfo.tab
    |   ├── FMindex
    |   ├── geneInfo.tab
    |   ├── genomeParameters.txt
    |   ├── sjdbInfo.txt
    |   ├── sjdbList.fromGTF.out.tab
    |   ├── sjdbList.out.tab
    |   └── transcriptInfo.tab
    ├── genes
    │   └── genes.gtf
    └── genome
        └── genome.fa

4 directories, 19 files
```

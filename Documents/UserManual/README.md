# Important Notifications
This notification page is a place to record and respond the user-reported issues/bugs and inform other users who may be affected by the same bug. Meanwhile, we inform users of important improvement and version changes here. We apologize for any inconvenience caused you.

## SAW ST v6.0.1 & v5.5.4 (Modified Apr. 17, 2023)

1. Bug fix: fixed the bug that the execution path of `checkGTF` function was invalid.


## SAW ST v5.5.3 (Modified Mar. 21, 2023)

1. Improvement description: upgraded `mapping` pipeline module, which performs polyA filtering after CID mapping, for improved statistical output in `report`. SAW V5 has an embedded polyA filtration function during `mapping`. In V5.5.2, polyA reads identification and filtration were done prior to CID mapping (an embedded function within `mapping` module) and these filtrated reads were recorded as 'Invalid CID Reads' which resulted in significantly low 'Valid CID Reads' outputs in certain cases. Hence, we adjusted the sequential order of polyA filtration step to be performed after CID mapping.
  - Modified the algorithm to optimize ordering capacity without compromising ordering stability.
  - PolyA filtering includes two main steps: 1) indentify and trim qualified reads according to length limitation, 2) remove trimmed reads of which length is less than 30.
  - Depending on the polyA ratio of the sequenced data, this adjustment will affect the final output expression matrix when much high polyA ratio were present. Because where to perform polyA filtering would have a slight influence on adapter and DNB filtering of `mapping` pipeline. 


## SAW ST v5.1.3 (Modified Dec. 2, 2022)
1. Issue description: SAW tissueCut pipeline failed and the error message shows `xxx line 1: 324233 Segmentation fault (core dumped). xxx`. (reported on: Oct. 24, 2022)
  - R: In SAW ST v5.1.3, the data type used for recording gene numbers is uint16 which can store a maximum of 65535 genes (using the ID index). However, for some exceptions where the gene number can exceed the upper limit, it will cause memory corruption and kill the program. 
  - Fixing this issue will be done in an upcoming release.
  - Update: This issue has been fixed in SAW V5.4.0, please update to the latest version if your sample may have more than 65535 genes. (released in: Dec. 2022)

2. Issue description: There is a typo in the example for `-p` and `--logo` in 5.1.3 Manual for report pipeline. (reported on: Oct. 21)
  - R: Yes, thanks for the comment. Sorry for the typo. The correct path should be `-p /opt/saw_v5.1.3_software/pipeline/report/plotly_package.txt` and `--logo /opt/saw_v5.1.3_software/pipeline/report/logo.png`.
  - This wil be fixed in the next manual release. 
  - Update: This has been fixed in manual A3.1 and A4.(fixed on: Dec. 2, 2022)

3. Issue description: Manual typo: `--bin1Saturation` is not show in the script exmaple, but it is a required parameter. (reported on: Oct. 25, 2022)
  - R: `--bin1Saturation` is actually an optional parameter, this is a typo in the manual.
  - Update: This has been fixed in manual A3.1 and A4.(fixed on: Dec. 2, 2022)

4. Issue description: Exon group in .tissue.gef is missing. (reported on: Nov. 8, 2022)
  - R: In SAW v5.1.3, the `tissueCut` output GEF didn't process the exon group. This issue will be fixed in the upcoming release of SAW. A circuitous way to get the exon group in v5.1.3 is to 1) convert .raw.gef and .tissue.gef to GEM format, 2) extract the rows in the raw.gem that has the same coordinate (x, y) with tissue.gem. In this way, you can get a tissue-coverage GEM table that has the ExonCount column. You may then convert this new GEM back into GEF format using gefpy or geftools.
  - Update: This issue has been fixed in SAW V5.4.0, please update to the latest version and rerun `tissueCut`. (released in: Dec. 2022)

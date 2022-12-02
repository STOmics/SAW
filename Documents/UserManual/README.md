# Important Notifications
This notification page is a place to record and respond the user-reported issues/bugs and inform other users who may be affected by the same bug. We apologize for any inconvenience caused you.


## SAW ST v5.1.3 (Modified Oct. 25, 2022)
1. Issue description: SAW tissueCut pipeline failed and the error message shows `xxx line 1: 324233 Segmentation fault (core dumped). xxx`. (reported on: Oct. 24, 2022)
  - R: In SAW ST v5.1.3, the data type used for recording gene numbers is uint16 which can store a maximum of 65535 genes (using the ID index). However, for some exceptions where the gene number can exceed the upper limit, it will cause memory corruption and kill the program. 
  - Fixing this issue will be done in an upcoming release.

2. Issue description: There is a typo in the example for `-p` and `--logo` in 5.1.3 Manual for report pipeline. (reported on: Oct. 21)
  - R: Yes, thanks for the comment. Sorry for the typo. The correct path should be `-p /opt/saw_v5.1.3_software/pipeline/report/plotly_package.txt` and `--logo /opt/saw_v5.1.3_software/pipeline/report/logo.png`.
  - This wil be fixed in the next manual release. 

3. Issue description: Manual typo: `--bin1Saturation` is not show in the script exmaple, but it is a required parameter. (reported on: Oct. 25, 2022)
  - R: `--bin1Saturation` is actually an optional parameter, this is a typo in the manual.

4. Issue description: Exon group in .tissue.gef is missing. (reported on: Nov. 8, 2022)
  - R: In SAW v5.1.3, the `tissueCut` output GEF didn't process the exon group. This issue will be fixed in the upcoming release of SAW. A circuitous way to get the exon group in v5.1.3 is to 1) convert .raw.gef and .tissue.gef to GEM format, 2) extract the rows in the raw.gem that has the same coordinate (x, y) with tissue.gem. In this way, you can get a tissue-coverage GEM table that has the ExonCount column. You may then convert this new GEM back into GEF format using gefpy or geftools.

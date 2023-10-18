# Important Notifications
This notification page is a place to record and respond the user-reported issues/bugs and inform other users who may be affected by the same bug. Meanwhile, we inform users of important improvement and version changes here. We apologize for any inconvenience caused you.

## SAW ST v6.1.3 (Modified Oct, 2023)
For `imageTools` module:    
1. Bug fix: fixed the problem that manually registered images were incorrectly cropped after scale operation.

For `report` module:    
1. Bug fix: fixed the missing percentage value of floating text box for 'Sequencing Saturation', in HTML diagram.

## SAW ST v6.1.2 (Modified Sep, 2023)
For `manualRegigster` module:   
1. Bug fix: fixed the offset calculation error of whcih the scale direction was incorrect.

## SAW ST v6.1.1 (Modified Aug. 11, 2023)
For `manualRegigster` module:      
1. Bug fix: fixed the issue that obtaining information from manual registration file `(.json)` incorrectly resulted in wrong width and height of images.

For `tissueCut` module:  
1. Bug fix: corrected the statistical method of `Reads_under_tissue` item in `tissueCut`, which now represents "Number of reads with position prior to filtration under tissue coverage".

For `report` module:    
1. Bug fix:  
compatible with the `.ipr` file which is failed with image quality control (QC) and has no record of heat map matrix information needed for html report.

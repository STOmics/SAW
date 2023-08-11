# SAW: Manual Processing Tutorial

**STOmics offline software version recommend**: ImageStudio V2.1 + SAW V6.1 + StereoMap V2.1

## Glossary

**Pre-registration**: Automatic stitching, tissue segmentation and cell segmentation based on Image. This function belongs to SAW `register`. 

**Registration**:  Image Pre-registration + Registration between Expression Matrix & Image. This function also belongs to SAW `register`. 


## QC-Success
### Scenarios

We have listed common scenarios. You can find the corresponding help according to the last module you used  in ImageStudio. 

For example, if you have finished manual processing in ImageStudio at the step of tissue segmentation, please check "After ImageStudio-Tissue/Cell Segmentation" for next.

	Access SAW first time:
		A. ImageStudio: QC -> SAW: register
		B. ImageStudio: QC + Stitching -> SAW: register
		C. ImageStudio: QC + Stitching + Tissue/Cell Segmentation -> SAW: register
	
	Not the first time to access SAW
		D. ImageStudio: QC -> SAW: Pre-registration/Registration -> ImageStudio: Stitching + Tissue/Cell Segmentation -> SAW: register
		E. ImageStudio: QC -> SAW: Pre-registration/Registration -> ImageStudio: Tissue/Cell Segmentation -> SAW: imageTools ipr2img

- ### After ImageStudio-QC

  - #### Pre-registration

    - Input Files: QC TARGZ + QC IPR

    ````bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## run pre-registration by SAW register which performs automatic stitching + tissue segmentation + cell segmentation.
    ## Module 'register' can be replaced by 'rapidRegister' if no need to automatic cell segmentation.
    singularity exec ${sif} register \
        -i <QC TARGZ> \
        -c <QC IPR> \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <QC IPR> replace to the path (directory required) of IPR obtained from ImageStudio-QC.
    # -o <outDir> replace to the result directory
    
    ## 'imageTools ipr2img' extracts images from IPR
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <register output IPR in outDir> \
        -d tissue cell \
        -r False \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <register output IPR in outDir> replace to the path (directory required) of IPR from the register output directory. 
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'False' to output pre-registered images.
    # -o <outDir> replace to the result directory which can be same as the Module 'register' output directory
    ````

    ```bash
    ## Use Cases
    singularity exec ${sif} register \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.ipr \
    	-o /path/to/test/01.qc/preregistration
    
    singularity exec ${sif} imageTools ipr2img \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/test/01.qc/preregistration/<ipr> \
    	-d tissue cell \
    	-r False \
    	-o /path/to/test/01.qc/preregistration
    ```

  - #### Registration

    - Input Files: QC TARGZ + QC IPR + GEF (raw )

    ```bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## run registration by SAW register which performs automatic stitching + tissue segmentation + cell segmentation + registration. 
    ## Module 'register' can be replaced by 'rapidRegister' if no need of automatic cell segmentation.
    singularity exec ${sif} register \
        -i <QC TARGZ> \
        -c <QC IPR> \
        -v <SAW raw GEF> \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <QC IPR> replace to the path (directory required) of IPR obtained from ImageStudio-QC.
    # -v <SAW raw GEF> replace to the path (directory required) of SN.raw.gef from SAW count directory.
    # -o <outDir> replace to the result directory
    
    ## 'imageTools ipr2img' extracts images from IPR
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <register output IPR in outDir> \
        -d tissue cell \
        -r True \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <register output IPR in outDir> replace to the path (directory required) of IPR from the register output directory. 
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'True' to output registered images.
    # -o <outDir> replace to the result directory which can be same as the Module 'register' output directory.
    ```

    ```bash
    ## Use Cases
    singularity exec ${sif} register \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.ipr \
    	-v /path/to/SAW_OUTPUT/02.count/SS200000059_NC.raw.gef \
    	-o /path/to/test/01.qc/registration
    
    singularity exec ${sif} imageTools ipr2img \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/test/01.qc/registration/<ipr> \
    	-d tissue cell \
    	-r True \
    	-o /path/to/test/01.qc/registration
    ```

- ### After ImageStudio-Image Stitching

  - #### Pre-registration

    - Input Files: QC TARGZ + Manual Stitching IPR

    ```bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## run pre-registration by SAW 'register' which performs automatic tissue segmentation + cell segmentation based on manual stitching.
    ## Module 'register' can be replaced by 'rapidRegister' if no need of automatic cell segmentation.
    singularity exec ${sif} register \
        -i <QC TARGZ> \
        -c <manual stitch IPR> \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual stitch IPR> replace to the path (directory required) of IPR obtained from ImageStudio-Image Stitching.
    # -o <outDir> replace to the result directory
    
    ## 'imageTools ipr2img' extracts images from IPR
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <register output IPR in outDir> \
        -d tissue cell \
        -r False \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <register output IPR in outDir> replace to the path (directory required) of IPR from the register output directory.
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'False' to output pre-registered images.
    # -o <outDir> replace to the result directory which can be same as the Module 'register' output directory
    ```

    ```Bash
    ## Use Cases
    singularity exec ${sif} register \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
        -c <manual stitch IPR> \
        -o /path/to/test/02.stitch/preregistration
    
    singularity exec ${sif} imageTools ipr2img \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/test/02.stitch/preregistration/<ipr> \
    	-d tissue cell \
    	-r False \
    	-o /path/to/test/02.stitch/preregistration
    ```

  - #### Registration

    - Input Files: QC TARGZ + Manual Stitching IPR + GEF (raw)

    ```bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## run registration by SAW register which performs automatic tissue segmentation + cell segmentation + registration based on manual stitching.
    ## Module 'register' can be replaced by 'rapidRegister' if no need of automatic cell segmentation.
    singularity exec ${sif} register \
        -i <QC TARGZ> \
        -c <manual stitch IPR> \
        -v <SAW raw GEF> \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual stitch IPR> replace to the path (directory required) of IPR obtained from ImageStudio-Image Stitching.
    # -v <SAW raw GEF> replace to the path (directory required) of SN.raw.gef from SAW count directory.
    # -o <outDir> replace to the result directory
    
    ## 'imageTools ipr2img' can convert images from IPR
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <register output IPR in outDir> \
        -d tissue cell \
        -r True \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <register output IPR in outDir> replace to the path (directory required) of IPR from the register output directory. 
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'True'  to output registered images.
    # -o <outDir> replace to the result directory which can be same as the Module 'register' output directory
    ```

    ```bash
    ## Use Cases
    singularity exec ${sif} register \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c <manual stitch IPR> \
    	-v /path/to/SAW_OUTPUT/02.count/SS200000059_NC.raw.gef \
    	-o /path/to/test/02.stitch/registration
    
    singularity exec ${sif} imageTools ipr2img \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/test/02.stitch/registration/<ipr> \
    	-d tissue cell \
    	-r True \
    	-o /path/to/test/02.stitch/registration
    ```

- ### After ImageStudio-Tissue/Cell Segmentation

  - #### Pre-registration again

    - Input Files: QC TARGZ+ Manual Segmentation IPR
    - Applicable scenes:
      - ImageStudio: QC -> SAW: Pre-registration/Registration -> **ImageStudio: Tissue/Cell Segmentation** -> SAW: imageTools ipr2img
      - ImageStudio: QC + Stitching -> SAW: Pre-registration/Registration -> **ImageStudio: Tissue/Cell Segmentation** -> SAW: imageTools ipr2img

    ````bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Because the manual segmentation info has already written into the IPR file, and no registration is required, you can skip the part of running the SAW register and directly use imageTools ipr2img to output the image
    ## 'imageTools ipr2img' can convert images from IPR
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <manual seg IPR> \
        -d tissue cell \
        -r False \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual seg IPR> replace to the path (directory required) of IPR from ImageStudio-Tissue/Cell segmentation.
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'False'  to output pre-registered images.
    # -o <outDir> replace to the result directory which can be same/not same as the pre-regisration output directory. It's recommended that put <manual seg IPR> in the <outDir>
    ````

    ```Bash
    ## Use Cases
    singularity exec ${sif} imageTools ipr2img \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \ 
    	-c /path/to/test/03.tissueseg/<ipr> \
    	-d tissue cell \
    	-r False \
    	-o /path/to/test/03.tissueseg/preregistration
    ```

  - #### Registration

    We list two different inputs (pre-registration/registration directory or QC TARGZ) in parameter `-i` in SAW `register` .

    

    - Input Files: **Pre-registration/Registration Directory** + Manual Segmentation IPR + GEF(raw)
      - Applicable scenes:
        - ImageStudio: QC -> SAW: Pre-registration/Registration -> **ImageStudio: Tissue/Cell Segmentation** -> SAW: register
        - ImageStudio: QC + Stitching -> SAW: Pre-registration/Registration -> **ImageStudio: Tissue/Cell Segmentation** -> SAW: register

    ```bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Module 'register' writes the registration information in the manual segmentation IPR file. In essence, it performs registration based on the pre-registration stitching results;
    ## Run 'register' to register the pre-registration stitching & manual-segmented image with the matrix
    singularity exec ${sif} register \
        -i <preregistration directory path> \
        -c <manual seg IPR> \
        -v <SAW raw GEF> \
        -o <outDir>
    # -i <preregistration directory path> replace to the prer-registration directory. Also accept the registration directorty.
    # -c <manual seg IPR> replace to the path (directory required) of IPR obtained from ImageStudio-Tissue/Cell Segmentation.
    # -v <SAW raw GEF> replace to the path (directory required) of SN.raw.gef from SAW count directory.
    # -o <outDir> replace to the result directory which stores registered manually-segmented IPR
    
    ## 'imageTools ipr2img' can convert images from IPR
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <register output IPR in outDir> \
        -d tissue cell \
        -r True \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <register output IPR in outDir> replace to the path (directory required) of IPR from the register output directory. 
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'True' to output registered images. You can also fill in 'False' if need to output pre-registered images
    # -o <outDir> replace to the result directory which can be same as the Module 'register' output directory
    ```
  
    ```bash
    ## Use Cases
    singularity exec ${sif} register \
    	-i /path/to/test/02.stitch/registration \
    	-c <manual seg IPR> \
    	-v </path/to/SS200000059_NC.raw.gef> \
    	-o /path/to/test/03.tissueseg/registration
    
    singularity exec ${sif} imageTools ipr2img \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/test/03.tissueseg/registration/<ipr> \
    	-d tissue cell \
    	-r True \
    	-o /path/to/test/03.tissueseg/registration
    ```
  
    
  
    - Input Files: **QC TARGZ** + Manual Segmentation IPR + GEF(raw)
      - Applicable scenes:
        - ImageStudio: QC + Stitching + **Tissue/Cell Segmentation** -> SAW: register
  
  - ```bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Run 'register' to generate the stitching & manual-segmentation result, and align the image with the expression matrix
    singularity exec ${sif} register \
        -i <QC TARGZ> \
        -c <manual seg IPR> \
        -v <SAW raw GEF> \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual seg IPR> replace to the path (directory required) of IPR obtained from ImageStudio-Tissue/Cell Segmentation.
    # -v <SAW raw GEF> replace to the path (directory required) of SN.raw.gef from SAW count directory.
    # -o <outDir> replace to the result directory which stores registered manually-segmented IPR
    
    ## 'imageTools ipr2img' can convert images from IPR
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <register output IPR in outDir> \
        -d tissue cell \
        -r True \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <register output IPR in outDir> replace to the path (directory required) of IPR from the register output directory. 
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'True' to output registered images. You can also fill in 'False' if output pre-registered images.
    # -o <outDir> replace to the result directory which can be same as the Module 'register' output directory
    ```
  
    ````bash
    ## Use Cases
    singularity exec ${sif} register \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c <manual seg IPR> \
    	-v </path/to/SS200000059_NC.raw.gef> \
    	-o /path/to/test/03.tissueseg/registration
    
    singularity exec ${sif} imageTools ipr2img \
    	-i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
    	-c /path/to/test/03.tissueseg/registration/<ipr> \
    	-d tissue cell \
    	-r True \
    	-o /path/to/test/03.tissueseg/registration
    ````


​    

- ### After Manual Immunofluorescence-threshold-segmentation (ImageStudio-Tissue/Cell Segmentation) 

  - #### Extract IF(Immunofluorescence) Expression Matrix + Visualization

    - Input Files: Manual-threshold-segmented IPR + GEF(raw) + SN.merge.barcodeReadsCount.txt + {stainType}_SN_tissue_cut.tif (from ipr2img)
    - Applicable scenes: Registration has been done

    ```Bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Because there is no need for registration, you can skip the part of running 'register' and directly use imageTools ipr2img to output images
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <manual seg IPR> \
        -d tissue cell \
        -r True \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual seg IPR> replace to the pathof manual-threshold-segmentated IPR.
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r Fill in 'True'. Here is the result after registration. 
    # -o <outDir> is replaced by the directory path to save the output results.It is recommended to create an empty directory and put <manual seg IPR> in it.
    
    # run 'tissueCut' to extract the manual-IF-threshold-segmentated expression file
    singularity exec ${sif} tissueCut \
        --dnbfile <SAW barcodeReadsCount TXT> \
        -i <SAW raw GEF>  \
        -o <outDir> \
        -s <IF Tissue Mask> \
        --sn <SN> -O Transcriptomics -d
    # --dnbfile <SAW barcodeReadsCount TXT> replace the path (directory required) of SN.merge.barcodeReadsCount.txt from SAW merge
    # -i <SAW raw GEF> replace to the path (directory required) of SN.raw.gef from SAW count directory.
    # -o <outDir> output directory
    # -s <IF Tissue Mask> replace to the path (directory required) of IF tissue mask image (from imageTools ipr2img)
    # --sn chip serial number
    # -O omics name
    # -d Generate statistical pictures
    
    # run cellCut bgef to complete GEF for visualization
    singularity exec ${sif} cellCut bgef \
        -i <outDir>/<SN Tissue GEF> \
        -o <outDir>/<complete SN Tissue GEF> \
        -b 1,10,20,50,100,200,500 \
        -O Transcriptomics
    # -i <outDir>/<SN Tissue GEF> replace to the path (directory required) of tissue.gef from 'tissuecut'
    # -o <outDir>/<complete SN Tissue GEF> replace to the path (directory required) of the complete gef
    # -b bin size
    # -O omics Name
    
    ```

    ```Bash
    ## Use Cases
    mkdir -p /path/to/test/03.tissueseg/IF_threshold/AKAP3_IF
    
    singularity exec ${sif} imageTools ipr2img \
        -i /path/to/<dataDir>/QC/SS200000059_NC_SC_20230309_094034_2.0.0.tar.gz \
        -c <manual seg IPR> \
        -d tissue cell\
        -r True \
        -o /ldfssz1/ST_BI/USER/st_stereotools/test/03.tissueseg/IF_threshold/AKAP3_IF
    
    singularity exec ${sif} tissueCut \
        --dnbfile /path/to/01.merge/SS200000059_NC.merge.barcodeReadsCount.txt \
        -i /path/to/03.count/SS200000059_NC.raw.gef \
        -o /path/to/test/03.tissueseg/IF_threshold/AKAP3_IF \
        -s /path/to/test/03.tissueseg/IF_gef_3_c/AKAP3_IF/AKAP3_IF_SS200000059_NC_tissue_cut.tif \
        --sn SS200000059_NC -O Transcriptomics -d
    
    singularity exec ${sif} cellCut bgef \
        -i /path/to/test/03.tissueseg/IF_threshold/AKAP3_IF/SS200000059_NC.tissue.gef \
        -o /path/to/test/03.tissueseg/IF_threshold/AKAP3_IF/SS200000059_NC.AKAP3_IF.gef \
        -b 1,10,20,50,100,200,500 \
        -O Transcriptomics
    ```

- ### View Results

  - The output TIFF image can be opened in imageJ, and the output RPI file can be viewed in StereoMap

  - Generate `fov_stitched_transformed.rpi` for manual registration (StereoMap). Other RPI can also be generated with `imageTools img2rpi`

    ```Bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    # imageTools img2rpi converts the output TIFF image into an image pyramid RPI file readable by StereoMap
    singularity exec ${sif} imageTools img2rpi \
            -i <*_fov_stitched_transformed.tif> \
            -g <stainType>/Image \
            -b 1 10 50 100 \
            -o <fov_stitched_transformed.rpi>
    # -i <*_fov_stitched_transformed.tif> replace to the file path (directory required) of the TIFF image that needs to be made into an image pyramid and provided by StereoMap
    # -g <stainType> Indicates the group name of the input TIFF image in RPI (such as DAPI)
    # -b The binsize of each layer in the image pyramid, and multiple binsizes are separated by spaces. The minimum binsize used for manual registration is 1 in order to provide more accurate images for manual registration, and the minimum binsize of SN.rpi, which is usually 2 to reduce file size
    # -o <fov_stitched_transformed.rpi> replace to the file path (directory required) of the output RPI
    ```


## QC-Fail

### Demo of SAW manual pipeline

- ####  ImageStudio: QC + Stitch + Tissue Segmentation (Case 1)

   - In this case, all three requirements of IPR from ImageStudio must be met before running SAW manual pipeline:

     - ImageStudio-QC has been done

     - ImageStudio-Image Stitching has been done

     - ImageStudio-Tissue Segmentation has been done

   - Run `stereoPipeline_v6.1.0_manual_part1.sh`  to obtain Expression Matrix and RPI for manual registration

     - [stereoPipeline_v6.1.0_manual_part1.sh](https://github.com/STOmics/SAW/blob/main/Scripts/stereoPipeline_v6.1.0_manual_part1.sh)

       ```Bash
       ## run_v6.1.0_manual_part1.sh
       ulimit -n 10240
       ulimit -v 33170449147
       
       NUMBA_CACHE_DIR=<absolute path to output directory>/tmp
       export NUMBA_CACHE_DIR=<absolute path to output directory>/tmp
       MPLCONFIGDIR=<absolute path to output directory>/tmp
       export MPLCONFIGDIR=<absolute path to output directory>/tmp
       
       sif=SAW_<version>.sif
       
       dataDir=<absolute path to data directory>
       outDir=<absolute path to output directory>
       QCDir=<absolute path to manual-processed IPR directory> ## replace to the directory of manual-processed IPR
       
       ## no cellbin
       export SINGULARITY_BIND=${dataDir},${outDir},${QCDir},${NUMBA_CACHE_DIR},${MPLCONFIGDIR}
       export PATH=/share/app/singularity/3.8.1/bin:$PATH
       bash stereoPipeline_v6.1.0_manual1.sh \
           -genomeSize 2.7 \
           -splitCount 1 \
           -maskFile ${dataDir}/mask/<SN>.barcodeToPos.h5 \
           -fq1 ${dataDir}/reads/<lane>_read_1.fq.gz \
           -fq2 ${dataDir}/reads/<lane>_read_2.fq.gz \
           -speciesName <speciesName> \
           -tissueType <tissueName> \
           -imageRecordFile ${QCDir}/<SN_date_time_version>.ipr \
           -imageCompressedFile ${QCDir}/<SN_date_time_version>.tar.gz \
           -refIndex /path/to/reference/STAR_SJ100 \
           -annotationFile /path/to/reference/genes.gtf \
           -sif ${sif} \
           -threads 24 \
           -outDir ${outDir}
       ```
   
- #### StereoMap: Manual Registration

  - Download following files and input these two files in StereoMap:

    - The complete expression matrix `<part1_outdir>/03.register/manual_register/SN.gef`
    - The stitched panoramic image `<part1_outdir>/03.register/manual_register/fov_stitched.rpi`

  - Click the Manual Registration Logo to do manual registration

  - Click `Submit` to save the registration result

  - Upload `.../Regist/<date_time>.regist.json` to SAW analysis directory

  - Run `stereoPipeline_v6.1.0_manual_part2.sh` to obtain tissue segmentation + cluster analysis + report

    - [stereoPipeline_v6.1.0_manual_part2.sh](https://github.com/STOmics/SAW/blob/main/Scripts/stereoPipeline_v6.1.0_manual_part2.sh)
    
    ```Bash
    ## run_v6.1.0_manual_part2.sh
    ulimit -n 10240
    ulimit -v 33170449147
    
    NUMBA_CACHE_DIR=<absolute path to output directory>/tmp
    export NUMBA_CACHE_DIR=<absolute path to output directory>/tmp
    MPLCONFIGDIR=<absolute path to output directory>/tmp
    export MPLCONFIGDIR=<absolute path to output directory>/tmp
    
    sif=SAW_<version>.sif
    
    dataDir=<absolute path to data directory>
    outDir=<absolute path to output directory> ## replace to part1_outdir
    QCDir=<absolute path to manual-processed IPR directory> ## replace to the directory of manual-processed IPR
    registJsonDir=<absolute path to regist.json directory> ## replace to the directory of regist.json directory
    
    ## no cellbin
    export SINGULARITY_BIND=${dataDir},${outDir},${QCDir},${NUMBA_CACHE_DIR},${registJsonDir},${MPLCONFIGDIR}
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    bash stereoPipeline_v6.1.0_manual2.sh \
        -SN <SN> \
        -countDir ${outDir} \
        -registJson ${registJsonDir}/<date_time>.regist.json \
        -speciesName <speciesName> \
        -tissueType <tissueName> \
        -imageRecordFile ${QCDir}/<SN_date_time_version>.ipr \
        -imageCompressedFile ${QCDir}/<SN_date_time_version>.tar.gz \
        -doCellBin N \
        -sif ${sif} \
        -threads 24 \
        -outDir ${outDir}
    ```
  
- #### View Results

  - The RPI file can be viewed in StereoMap. The output TIFF image can be opened in imageJ. 

  - Input the following files in StereoMap and view the results.

    ```Bash
    ## RPI
    SN.rpi
    ## Expression Matrix
    SN.gef
    ## bin200 spatial cluster and UMAP
    SN.spatial.cluster.h5ad
    ```



### Other cases in ImageStudio  

- #### ImageStudio: QC + Stitch

  - Input Files: QC TARGZ + Manual Stitching IPR

  - Modify `stereoPipeline_v6.1.0_manual_part1.sh` corresponding line

    ```bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Run imageTools ipr2img for stitched panoramic TIFF
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <manual stitch IPR> \
        -r False \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual stitch IPR> replace to the path (directory required) of IPR from ImageStudio-Image Stitching.
    # -r False No registration has been done. The output TIFF image will not be reversed, rotated and translated according to the registration parameters
    # -o <outDir> The output directory
    
    ## Run imageTools img2rpi for RPI                        
    singularity exec ${sif} imageTools img2rpi \
        -i <*_fov_stitched.tif> \
        -g <stainType>/Image \
        -b 1 10 50 100 \
        -o <fov_stitched.rpi>
    # -i <*_fov_stitched_transformed.tif> replace to the file path (directory required) of the TIFF image that needs to be made into an image pyramid
    # -g <stainType> Indicates the group name of the input TIFF image in RPI (such as DAPI)
    # -b The binsize of each layer in the image pyramid, and multiple binsizes are separated by spaces. The minimum binsize used for manual registration is 1 in order to provide more accurate images for manual registration, and the minimum binsize of SN.rpi, which is usually 2 to reduce file size
    # -o <fov_stitched_transformed.rpi> replace to the file path (directory required) of the output RPI
    ```

    ```bash
    #Use Cases
    singularity exec ${sif} imageTools ipr2img \
        -i /path/to/qc+stitch/${SN}_B6_SC_20230516_174736_2.0.0.tar.gz \
        -c /path/to/qc+stitch/${SN}_B6_SC_20230516_183205_2.0.0.ipr \
        -r False \
        -o /path/to/QC_Fail/result/QC_Stitch
                             
    singularity exec ${sif} imageTools img2rpi \
        -i /path/to/QC_Fail/result/QC_Stitch/ssDNA_fov_stitched.tif \
        -g ssDNA/Image \
        -b 1 10 50 100 \
        -o /path/to/QC_Fail/result/QC_Stitch/fov_stitched.rpi
    ```

  - Complete the manual registration in StereoMap and upload `regist.json` to the analysis path

  - Confirm the cellbin paramter in `run_v6.1.0_manual_part2.sh `  is  `-doCellBin N`

  - Modify `stereoPipeline_v6.1.0_manual_part2.sh`: Due to the **lack** of tissue segmentation and cell segmentation results, 

    - use SAW `tissueCut` **without image**, 
    - do **not** run `cellCut`  and  `cellCluster`,
    - and check `report ` input parameters

- #### ImageStudio: QC + Stitch + Tissue Segmentation+Cell Segmentation

  - Input Files: QC TARGZ + Manual Segmentation IPR

  - Modify `stereoPipeline_v6.1.0_manual_part1.sh` corresponding line

    ```bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Run imageTools ipr2img for stitched panoramic TIFF, tissue & cell mask TIFF
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <manual seg IPR> \
        -d tissue cell \
        -r False \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual stitch IPR> replace to the path (directory required) of IPR from last operated module in ImageStudio.
    # -d 'tissue' and 'cell' indicate the tissue/cell segmentation TIFF image separately. If two kinds of images need to be output, separate them with spaces.
    # -r False No registration has been done. The output TIFF image will not be reversed, rotated and translated according to the registration parameters
    # -o <outDir> The output directory
    
    ## Run imageTools img2rpi for RPI                   
    singularity exec ${sif} imageTools img2rpi \
        -i <*_fov_stitched.tif> \
        -g <stainType>/Image \
        -b 1 10 50 100 \
        -o <fov_stitched.rpi>
    # -i <*_fov_stitched_transformed.tif> replace to the file path (directory required) of the TIFF image that needs to be made into an image pyramid
    # -g <stainType> Indicates the group name of the input TIFF image in RPI (such as DAPI)
    # -b The binsize of each layer in the image pyramid, and multiple binsizes are separated by spaces. The minimum binsize used for manual registration is 1 in order to provide more accurate images for manual registration, and the minimum binsize of SN.rpi, which is usually 2 to reduce file size
    # -o <fov_stitched_transformed.rpi> replace to the file path (directory required) of the output RPI
    ```

    ````bash
    #Use Cases
    singularity exec ${sif} imageTools ipr2img \
        -i /path/to/qc+stitch/${SN}_B6_SC_20230516_174736_2.0.0.tar.gz \
        -c /path/to/qc+stitch/${SN}_B6_SC_20230516_183205_2.0.0.ipr \
        -r False \
        -o /path/to/QC_Fail/result/QC_Stitch
                             
    singularity exec ${sif} imageTools img2rpi \
        -i /path/to/QC_Fail/result/QC_Stitch/ssDNA_fov_stitched.tif \
        -g ssDNA/Image \
        -b 1 10 50 100 \
        -o /path/to/QC_Fail/result/QC_Stitch/fov_stitched.rpi
    ````

  - Complete  the manual registration in StereoMap and upload `regist.json` to the analysis path

  - Confirm the cellbin paramter in `run_v6.1.0_manual_part2.sh ` is  `-doCellBin Y`

  - Modify `stereoPipeline_v6.1.0_manual_part2.sh`: Due to the **generation** of tissue segmentation and cell segmentation results, 

    - use SAW `tissueCut` **with image**, 
    - run `cellCut`  and  `cellCluster`,
    - and check `report ` input parameters

- #### ImageStudio: QC + Tissue Segmentation

  - Input Files: QC TARGZ + Manual Tissue Segmentation IPR

  - Modify `stereoPipeline_v6.1.0_manual_part1.sh` corresponding line

    ```Bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Run imageTools ipr2img for stitched panoramic TIFF and tissue mask TIFF
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <manual seg IPR> \
        -d tissue \
        -r False \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual stitch IPR> replace to the path (directory required) of IPR from ImageStudio-Tissue Segmentation.
    # -d tissue need to output the tissue segmentation TIFF image
    # -r False No registration has been done. The output TIFF image will not be reversed, rotated and translated according to the registration parameters
    # -o <outDir> The output directory
    
    ## Run imageTools img2rpi for RPI             
    singularity exec ${sif} imageTools img2rpi \
        -i <*_fov_stitched.tif> \
        -g <stainType>/Image \
        -b 1 10 50 100 \
        -o <fov_stitched.rpi>
    # -i <*_fov_stitched_transformed.tif> replace to the file path (directory required) of the TIFF image that needs to be made into an image pyramid
    # -g <stainType> Indicates the group name of the input TIFF image in RPI (such as DAPI)
    # -b The binsize of each layer in the image pyramid, and multiple binsizes are separated by spaces. The minimum binsize used for manual registration is 1 in order to provide more accurate images for manual registration, and the minimum binsize of SN.rpi, which is usually 2 to reduce file size
    # -o <fov_stitched_transformed.rpi> replace to the file path (directory required) of the output RPI
    ```

    ```Bash
    # Use Cases
    singularity exec ${sif} imageTools ipr2img \
        -i /path/to/qc+Tissue/${SN}_B6_SC_20230516_175051_2.0.0.tar.gz \
        -c /path/to/qc+Tissue/${SN}_B6_SC_20230516_185938_2.0.0.ipr \
        -d tissue \
        -r False \
        -o /path/to/QC_Fail/result/QC_Tissue/
                 
    singularity exec ${sif} imageTools img2rpi \
        -i /path/to/QC_Fail/result/QC_Tissue/DAPI_fov_stitched.tif \
        -g DAPI/Image \
        -b 1 10 50 100 \
        -o /path/to/QC_Fail/result/QC_Tissue/fov_stitched.rpi
    ```

  - Complete manual registration in StereoMap and upload `regist.json` to the analysis path

  - Confirm that the cellbin paramter in `run_v6.1.0_manual_part2.sh`  is `-doCellBin N`

  - Modify `stereoPipeline_v6.1.0_manual_part2.sh`: Due to the **generation** of tissue segmentation results, 

    - use SAW `tissueCut` **with image**, 
    - do **not** run `cellCut`  and  `cellCluster`,
    - and check `report ` input parameters

- #### ImageStudio-QC + Cell Segmentation

  - Input Files: QC TARGZ + Manual Cell Segmentation IPR

  - Modify `stereoPipeline_v6.1.0_manual_part1.sh` corresponding line

    ```Bash
    ## copy
    export HDF5_USE_FILE_LOCKING=FALSE
    export PATH=/share/app/singularity/3.8.1/bin:$PATH
    sif=SAW_<version>.sif
    
    ## replace with real paths
    dataDir=<absolute path to data directory>  ## raw data directory, such as FASTQ,TARGZ,IPR
    outDir=<absolute path to output directory>  ## output directory
    export SINGULARITY_BIND=${dataDir},${outDir},"<absolute path to other directories>" ## bind paths
    
    ## Run imageTools ipr2img for stitched panoramic TIFF and cell mask TIFF
    singularity exec ${sif} imageTools ipr2img \
        -i <QC TARGZ> \
        -c <manual seg IPR> \
        -d cell \
        -r False \
        -o <outDir>
    # -i <QC TARGZ> replace to the path (directory required) of TARGZ obtained from ImageStudio-QC.
    # -c <manual stitch IPR> replace to the path (directory required) of IPR from ImageStudio-Cell Segmentation.
    # -d cell need to output the cell segmentation TIFF image
    # -r False No registration has been done. The output TIFF image will not be reversed, rotated and translated according to the registration parameters
    # -o <outDir> The output directory
    
    ## Run imageTools img2rpi for RPI
    singularity exec ${sif} imageTools img2rpi \
        -i <*_fov_stitched.tif> \
        -g <stainType>/Image \
        -b 1 10 50 100 \
        -o <fov_stitched.rpi>
    # -i <*_fov_stitched_transformed.tif> replace to the file path (directory required) of the TIFF image that needs to be made into an image pyramid
    # -g <stainType> Indicates the group name of the input TIFF image in RPI (such as DAPI)
    # -b The binsize of each layer in the image pyramid, and multiple binsizes are separated by spaces. The minimum binsize used for manual registration is 1 in order to provide more accurate images for manual registration, and the minimum binsize of SN.rpi, which is usually 2 to reduce file size
    # -o <fov_stitched_transformed.rpi> replace to the file path (directory required) of the output RPI
    ```

  - Complete manual registration in StereoMap and upload `regist.json` to the analysis path

  - Confirm that the cellbin paramter in `run_v6.1.0_manual_part2.sh `  is  `-doCellBin Y`

  - Modify `stereoPipeline_v6.1.0_manual_part2.sh`: Due to the **generation** of cell segmentation results, 

    - use SAW `tissueCut` **without image**, 
    - run `cellCut`  and  `cellCluster`,
    - and check `report ` input parameters.

​	




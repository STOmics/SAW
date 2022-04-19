path_data=/home/ubuntu/Test_Data/SS200000003BR_B3
ulimit -n 10240
bash stereo_mul_4.1.0.sh \
        -m $path_data/mask/SS200000003BR_B3.barcodeToPos.h5 \
        -1 $path_data/reads/V350044321_L01_read_1.fq.gz,$path_data/reads/E100026545_L01_cp_trim_read_1.fq.gz \
        -2 $path_data/reads/V350044321_L01_read_2.fq.gz,$path_data/reads/E100026545_L01_cp_trim_read_2.fq.gz \
        -g $path_data/reference/STAR_SJ100 \
        -a $path_data/reference/genes.gtf \
	-c 5 \
        -i $path_data/image/SS200000003BR_B3 \
	-s /home/ubuntu/TEST/SAW_v4/image_sif/SAW_v4.1.0.sif \
        -o /home/ubuntu/TEST/SAW_v4/result/SS200000003BR_B3_trim_t2

#!/bin/bash


#BSUB -L /bin/bash
#BSUB -o transcript4.txt
#BSUB -e err4.txt
#BSUB -J trans4_test
#BSUB -u tunaiya@gmail.com
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10194304

cd /scratch/cluster/monthly/tina/
mkdir cavefish/
cd cavefish/

wget ftp://ftp.ensembl.org/pub/release-92/fasta/astyanax_mexicanus/cdna/Astyanax_mexicanus.AstMex102.cdna.all.fa.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_11_CTTGTA_L006_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_11_CTTGTA_L006_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_12_GTGGCC_L008_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_12_GTGGCC_L008_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_13_ACAGTG_L003_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_13_ACAGTG_L003_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_2_ATCACG_L005_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_2_ATCACG_L005_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_3_ACTTGA_L005_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_3_ACTTGA_L005_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_4_TAGCTT_L005_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_4_TAGCTT_L005_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_5_GGCTAC_L005_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_5_GGCTAC_L005_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_6_TAGCTT_L008_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_6_TAGCTT_L008_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_7_GATCAG_L005_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_7_GATCAG_L005_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_8_ACTGAT_L005_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_8_ACTGAT_L005_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_9_GCCAAT_L006_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/DnD2XmsxUR/raw/T_Am_9_GCCAAT_L006_R2.fastq.gz


wait
echo "Downloading done..."



module add UHTS/Analysis/trimmomatic/0.36

trimmomatic PE -threads 32 T_Am_11_CTTGTA_L006_R1.fastq.gz T_Am_11_CTTGTA_L006_R2.fastq.gz T_Am_11_CTTGTA_L006_R1_trim.fastq.gz unpaired_1.fq T_Am_11_CTTGTA_L006_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_11_CTTGTA_L006_R1_trim.fastq.gz T_Am_11_CTTGTA_L006_R1.fastq.gz
mv T_Am_11_CTTGTA_L006_R2_trim.fastq.gz T_Am_11_CTTGTA_L006_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_12_GTGGCC_L008_R1.fastq.gz T_Am_12_GTGGCC_L008_R2.fastq.gz T_Am_12_GTGGCC_L008_R1_trim.fastq.gz unpaired_1.fq T_Am_12_GTGGCC_L008_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_12_GTGGCC_L008_R1_trim.fastq.gz T_Am_12_GTGGCC_L008_R1.fastq.gz
mv T_Am_12_GTGGCC_L008_R2_trim.fastq.gz T_Am_12_GTGGCC_L008_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_13_ACAGTG_L003_R1.fastq.gz T_Am_13_ACAGTG_L003_R2.fastq.gz T_Am_13_ACAGTG_L003_R1_trim.fastq.gz unpaired_1.fq T_Am_13_ACAGTG_L003_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_13_ACAGTG_L003_R1_trim.fastq.gz T_Am_13_ACAGTG_L003_R1.fastq.gz
mv T_Am_13_ACAGTG_L003_R2_trim.fastq.gz T_Am_13_ACAGTG_L003_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_2_ATCACG_L005_R1.fastq.gz T_Am_2_ATCACG_L005_R2.fastq.gz T_Am_2_ATCACG_L005_R1_trim.fastq.gz unpaired_1.fq T_Am_2_ATCACG_L005_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_2_ATCACG_L005_R1_trim.fastq.gz T_Am_2_ATCACG_L005_R1.fastq.gz
mv T_Am_2_ATCACG_L005_R2_trim.fastq.gz T_Am_2_ATCACG_L005_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_3_ACTTGA_L005_R1.fastq.gz T_Am_3_ACTTGA_L005_R2.fastq.gz T_Am_3_ACTTGA_L005_R1_trim.fastq.gz unpaired_1.fq T_Am_3_ACTTGA_L005_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_3_ACTTGA_L005_R1_trim.fastq.gz T_Am_3_ACTTGA_L005_R1.fastq.gz
mv T_Am_3_ACTTGA_L005_R2_trim.fastq.gz T_Am_3_ACTTGA_L005_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_4_TAGCTT_L005_R1.fastq.gz T_Am_4_TAGCTT_L005_R2.fastq.gz T_Am_4_TAGCTT_L005_R1_trim.fastq.gz unpaired_1.fq T_Am_4_TAGCTT_L005_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_4_TAGCTT_L005_R1_trim.fastq.gz T_Am_4_TAGCTT_L005_R1.fastq.gz
mv T_Am_4_TAGCTT_L005_R2_trim.fastq.gz T_Am_4_TAGCTT_L005_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_5_GGCTAC_L005_R1.fastq.gz T_Am_5_GGCTAC_L005_R2.fastq.gz T_Am_5_GGCTAC_L005_R1_trim.fastq.gz unpaired_1.fq T_Am_5_GGCTAC_L005_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_5_GGCTAC_L005_R1_trim.fastq.gz T_Am_5_GGCTAC_L005_R1.fastq.gz
mv T_Am_5_GGCTAC_L005_R2_trim.fastq.gz T_Am_5_GGCTAC_L005_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_6_TAGCTT_L008_R1.fastq.gz T_Am_6_TAGCTT_L008_R2.fastq.gz T_Am_6_TAGCTT_L008_R1_trim.fastq.gz unpaired_1.fq T_Am_6_TAGCTT_L008_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_6_TAGCTT_L008_R1_trim.fastq.gz T_Am_6_TAGCTT_L008_R1.fastq.gz
mv T_Am_6_TAGCTT_L008_R2_trim.fastq.gz T_Am_6_TAGCTT_L008_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_7_GATCAG_L005_R1.fastq.gz T_Am_7_GATCAG_L005_R2.fastq.gz T_Am_7_GATCAG_L005_R1_trim.fastq.gz unpaired_1.fq T_Am_7_GATCAG_L005_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_7_GATCAG_L005_R1_trim.fastq.gz T_Am_7_GATCAG_L005_R1.fastq.gz
mv T_Am_7_GATCAG_L005_R2_trim.fastq.gz T_Am_7_GATCAG_L005_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_8_ACTGAT_L005_R1.fastq.gz T_Am_8_ACTGAT_L005_R2.fastq.gz T_Am_8_ACTGAT_L005_R1_trim.fastq.gz unpaired_1.fq T_Am_8_ACTGAT_L005_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_8_ACTGAT_L005_R1_trim.fastq.gz T_Am_8_ACTGAT_L005_R1.fastq.gz
mv T_Am_8_ACTGAT_L005_R2_trim.fastq.gz T_Am_8_ACTGAT_L005_R2.fastq.gz

wait

trimmomatic PE -threads 32 T_Am_9_GCCAAT_L006_R1.fastq.gz T_Am_9_GCCAAT_L006_R2.fastq.gz T_Am_9_GCCAAT_L006_R1_trim.fastq.gz unpaired_1.fq T_Am_9_GCCAAT_L006_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv T_Am_9_GCCAAT_L006_R1_trim.fastq.gz T_Am_9_GCCAAT_L006_R1.fastq.gz
mv T_Am_9_GCCAAT_L006_R2_trim.fastq.gz T_Am_9_GCCAAT_L006_R2.fastq.gz

wait

module add UHTS/Analysis/kallisto/0.43.1

kallisto index -i cavefish_index Astyanax_mexicanus.AstMex102.cdna.all.fa.gz
wait
echo "kallisto index is build for this organism"

kallisto quant -i cavefish_index -o Embryo T_Am_11_CTTGTA_L006_R1.fastq.gz T_Am_11_CTTGTA_L006_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o NIovary T_Am_12_GTGGCC_L008_R1.fastq.gz T_Am_12_GTGGCC_L008_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o NItestis T_Am_13_ACAGTG_L003_R1.fastq.gz T_Am_13_ACAGTG_L003_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Brain T_Am_2_ATCACG_L005_R1.fastq.gz T_Am_2_ATCACG_L005_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Gills T_Am_3_ACTTGA_L005_R1.fastq.gz T_Am_3_ACTTGA_L005_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Heart T_Am_4_TAGCTT_L005_R1.fastq.gz T_Am_4_TAGCTT_L005_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Muscle T_Am_5_GGCTAC_L005_R1.fastq.gz T_Am_5_GGCTAC_L005_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Liver T_Am_6_TAGCTT_L008_R1.fastq.gz T_Am_6_TAGCTT_L008_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Headkidney T_Am_7_GATCAG_L005_R1.fastq.gz T_Am_7_GATCAG_L005_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Bones T_Am_8_ACTGAT_L005_R1.fastq.gz T_Am_8_ACTGAT_L005_R2.fastq.gz
wait
kallisto quant -i cavefish_index -o Intestine T_Am_9_GCCAAT_L006_R1.fastq.gz T_Am_9_GCCAAT_L006_R2.fastq.gz
wait
echo "Finished abundance work...."


rm *.gz

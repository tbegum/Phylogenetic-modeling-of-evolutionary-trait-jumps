#!/bin/bash


#BSUB -L /bin/bash
#BSUB -o transcript2.txt
#BSUB -e err2.txt
#BSUB -J trans2_test
#BSUB -u tunaiya@gmail.com
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10194304

cd /scratch/cluster/monthly/tina/
mkdir medaka/
cd medaka/

wget ftp://ftp.ensembl.org/pub/release-92/fasta/oryzias_latipes/cdna/Oryzias_latipes.MEDAKA1.cdna.all.fa.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_11_CTTGTA_L002_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_11_CTTGTA_L002_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_12_ATCACG_L003_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_12_ATCACG_L003_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_13_CGATGT_L003_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_13_CGATGT_L003_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_2_TGACCA_L001_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_2_TGACCA_L001_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_3_GCCAAT_L001_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_3_GCCAAT_L001_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_4_TTAGGC_L001_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_4_TTAGGC_L001_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_5_CAGATC_L002_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_5_CAGATC_L002_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_6_ACTTGA_L002_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_6_ACTTGA_L002_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_7_GATCAG_L002_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_7_GATCAG_L002_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_8_TAGCTT_L002_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_8_TAGCTT_L002_R2.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_9_GGCTAC_L002_R1.fastq.gz
wget http://phylofish.sigenae.org/ngspipelines/data/fF78IfGgrl/raw/L_Ol_9_GGCTAC_L002_R2.fastq.gz


wait
echo "Downloading done..."



module add UHTS/Analysis/trimmomatic/0.36

trimmomatic PE -threads 32 L_Ol_11_CTTGTA_L002_R1.fastq.gz L_Ol_11_CTTGTA_L002_R2.fastq.gz L_Ol_11_CTTGTA_L002_R1_trim.fastq.gz unpaired_1.fq L_Ol_11_CTTGTA_L002_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_11_CTTGTA_L002_R1_trim.fastq.gz L_Ol_11_CTTGTA_L002_R1.fastq.gz
mv L_Ol_11_CTTGTA_L002_R2_trim.fastq.gz L_Ol_11_CTTGTA_L002_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_12_ATCACG_L003_R1.fastq.gz L_Ol_12_ATCACG_L003_R2.fastq.gz L_Ol_12_ATCACG_L003_R1_trim.fastq.gz unpaired_1.fq L_Ol_12_ATCACG_L003_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_12_ATCACG_L003_R1_trim.fastq.gz L_Ol_12_ATCACG_L003_R1.fastq.gz
mv L_Ol_12_ATCACG_L003_R2_trim.fastq.gz L_Ol_12_ATCACG_L003_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_13_CGATGT_L003_R1.fastq.gz L_Ol_13_CGATGT_L003_R2.fastq.gz L_Ol_13_CGATGT_L003_R1_trim.fastq.gz unpaired_1.fq L_Ol_13_CGATGT_L003_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_13_CGATGT_L003_R1_trim.fastq.gz L_Ol_13_CGATGT_L003_R1.fastq.gz
mv L_Ol_13_CGATGT_L003_R2_trim.fastq.gz L_Ol_13_CGATGT_L003_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_2_TGACCA_L001_R1.fastq.gz L_Ol_2_TGACCA_L001_R2.fastq.gz L_Ol_2_TGACCA_L001_R1_trim.fastq.gz unpaired_1.fq L_Ol_2_TGACCA_L001_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_2_TGACCA_L001_R1_trim.fastq.gz L_Ol_2_TGACCA_L001_R1.fastq.gz
mv L_Ol_2_TGACCA_L001_R2_trim.fastq.gz L_Ol_2_TGACCA_L001_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_3_GCCAAT_L001_R1.fastq.gz L_Ol_3_GCCAAT_L001_R2.fastq.gz L_Ol_3_GCCAAT_L001_R1_trim.fastq.gz unpaired_1.fq L_Ol_3_GCCAAT_L001_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_3_GCCAAT_L001_R1_trim.fastq.gz L_Ol_3_GCCAAT_L001_R1.fastq.gz
mv L_Ol_3_GCCAAT_L001_R2_trim.fastq.gz L_Ol_3_GCCAAT_L001_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_4_TTAGGC_L001_R1.fastq.gz L_Ol_4_TTAGGC_L001_R2.fastq.gz L_Ol_4_TTAGGC_L001_R1_trim.fastq.gz unpaired_1.fq L_Ol_4_TTAGGC_L001_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_4_TTAGGC_L001_R1_trim.fastq.gz L_Ol_4_TTAGGC_L001_R1.fastq.gz
mv L_Ol_4_TTAGGC_L001_R2_trim.fastq.gz L_Ol_4_TTAGGC_L001_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_5_CAGATC_L002_R1.fastq.gz L_Ol_5_CAGATC_L002_R2.fastq.gz L_Ol_5_CAGATC_L002_R1_trim.fastq.gz unpaired_1.fq L_Ol_5_CAGATC_L002_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_5_CAGATC_L002_R1_trim.fastq.gz L_Ol_5_CAGATC_L002_R1.fastq.gz
mv L_Ol_5_CAGATC_L002_R2_trim.fastq.gz L_Ol_5_CAGATC_L002_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_6_ACTTGA_L002_R1.fastq.gz L_Ol_6_ACTTGA_L002_R2.fastq.gz L_Ol_6_ACTTGA_L002_R1_trim.fastq.gz unpaired_1.fq L_Ol_6_ACTTGA_L002_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_6_ACTTGA_L002_R1_trim.fastq.gz L_Ol_6_ACTTGA_L002_R1.fastq.gz
mv L_Ol_6_ACTTGA_L002_R2_trim.fastq.gz L_Ol_6_ACTTGA_L002_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_7_GATCAG_L002_R1.fastq.gz L_Ol_7_GATCAG_L002_R2.fastq.gz L_Ol_7_GATCAG_L002_R1_trim.fastq.gz unpaired_1.fq L_Ol_7_GATCAG_L002_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_7_GATCAG_L002_R1_trim.fastq.gz L_Ol_7_GATCAG_L002_R1.fastq.gz
mv L_Ol_7_GATCAG_L002_R2_trim.fastq.gz L_Ol_7_GATCAG_L002_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_8_TAGCTT_L002_R1.fastq.gz L_Ol_8_TAGCTT_L002_R2.fastq.gz L_Ol_8_TAGCTT_L002_R1_trim.fastq.gz unpaired_1.fq L_Ol_8_TAGCTT_L002_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_8_TAGCTT_L002_R1_trim.fastq.gz L_Ol_8_TAGCTT_L002_R1.fastq.gz
mv L_Ol_8_TAGCTT_L002_R2_trim.fastq.gz L_Ol_8_TAGCTT_L002_R2.fastq.gz

wait

trimmomatic PE -threads 32 L_Ol_9_GGCTAC_L002_R1.fastq.gz L_Ol_9_GGCTAC_L002_R2.fastq.gz L_Ol_9_GGCTAC_L002_R1_trim.fastq.gz unpaired_1.fq L_Ol_9_GGCTAC_L002_R2_trim.fastq.gz unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25


wait
rm unpaired_1.fq
rm unpaired_2.fq
mv L_Ol_9_GGCTAC_L002_R1_trim.fastq.gz L_Ol_9_GGCTAC_L002_R1.fastq.gz
mv L_Ol_9_GGCTAC_L002_R2_trim.fastq.gz L_Ol_9_GGCTAC_L002_R2.fastq.gz

wait

module add UHTS/Analysis/kallisto/0.43.1

kallisto index -i medaka_index Oryzias_latipes.MEDAKA1.cdna.all.fa.gz
wait
echo "kallisto index is build for this organism"

kallisto quant -i medaka_index -o Embryo L_Ol_11_CTTGTA_L002_R1.fastq.gz L_Ol_11_CTTGTA_L002_R2.fastq.gz
wait
kallisto quant -i medaka_index -o NIovary L_Ol_12_ATCACG_L003_R1.fastq.gz L_Ol_12_ATCACG_L003_R2.fastq.gz
wait
kallisto quant -i medaka_index -o NItestis L_Ol_13_CGATGT_L003_R1.fastq.gz L_Ol_13_CGATGT_L003_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Brain L_Ol_2_TGACCA_L001_R1.fastq.gz L_Ol_2_TGACCA_L001_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Gills L_Ol_3_GCCAAT_L001_R1.fastq.gz L_Ol_3_GCCAAT_L001_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Heart L_Ol_4_TTAGGC_L001_R1.fastq.gz L_Ol_4_TTAGGC_L001_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Muscle L_Ol_5_CAGATC_L002_R1.fastq.gz L_Ol_5_CAGATC_L002_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Liver L_Ol_6_ACTTGA_L002_R1.fastq.gz L_Ol_6_ACTTGA_L002_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Headkidney L_Ol_7_GATCAG_L002_R1.fastq.gz L_Ol_7_GATCAG_L002_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Bones L_Ol_8_TAGCTT_L002_R1.fastq.gz L_Ol_8_TAGCTT_L002_R2.fastq.gz
wait
kallisto quant -i medaka_index -o Intestine L_Ol_9_GGCTAC_L002_R1.fastq.gz L_Ol_9_GGCTAC_L002_R2.fastq.gz
wait
echo "Finished abundance work...."


rm *.gz

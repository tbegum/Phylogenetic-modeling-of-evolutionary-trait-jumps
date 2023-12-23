#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o transcript.txt
#BSUB -e err.txt
#BSUB -J trans_test
#BSUB -u tunaiya@gmail.com
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10194304

cd /scratch/cluster/monthly/tina/

mkdir tilapia
cd tilapia/
mkdir Testis/
mkdir Liver/
mkdir Brain/
mkdir Heart/
mkdir Muscle/
mkdir Kidney/
mkdir Ovary/
mkdir Embryo/

cd Testis/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391690/SRR391690_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391690/SRR391690_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391695/SRR391695_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391695/SRR391695_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391701/SRR391701_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391701/SRR391701_2.fastq.gz
wait

cd ../Liver/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391688/SRR391688_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391688/SRR391688_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391698/SRR391698_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391698/SRR391698_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391708/SRR391708_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391708/SRR391708_2.fastq.gz
wait
 
cd ../Brain/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391697/SRR391697_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391697/SRR391697_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391699/SRR391699_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391699/SRR391699_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391709/SRR391709_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391709/SRR391709_2.fastq.gz
wait

cd ../Heart/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391681/SRR391681_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391681/SRR391681_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391686/SRR391686_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391686/SRR391686_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391696/SRR391696_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391696/SRR391696_2.fastq.gz
wait

cd ../Muscle/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391702/SRR391702_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391702/SRR391702_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391704/SRR391704_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391704/SRR391704_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391706/SRR391706_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391706/SRR391706_2.fastq.gz
wait

cd ../Kidney/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391680/SRR391680_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391680/SRR391680_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391684/SRR391684_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391684/SRR391684_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391689/SRR391689_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391689/SRR391689_2.fastq.gz
wait 

cd ../Ovary/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391687/SRR391687_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391687/SRR391687_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391691/SRR391691_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391691/SRR391691_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391693/SRR391693_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391693/SRR391693_2.fastq.gz
wait

cd ../Embryo/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391705/SRR391705_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391705/SRR391705_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391707/SRR391707_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391707/SRR391707_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391711/SRR391711_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR391/SRR391711/SRR391711_2.fastq.gz
wait
   


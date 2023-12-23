## Script to download pike tissue data

#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o transcript1.txt
#BSUB -e err1.txt
#BSUB -J trans1_test
#BSUB -u tunaiya@gmail.com
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10194304

cd /scratch/cluster/weekly/tina/

mkdir northern_pike
cd northern_pike/
mkdir Testis/
mkdir Liver/
mkdir Brain/
mkdir Bones/
mkdir Gills/
mkdir Heart/
mkdir Muscle/
mkdir Kidney/
mkdir Intestine/
mkdir Ovary/
mkdir Embryo/

cd Testis/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/001/SRR1533661/SRR1533661_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/001/SRR1533661/SRR1533661_2.fastq.gz
wait

cd ../Liver/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/006/SRR1533656/SRR1533656_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/006/SRR1533656/SRR1533656_2.fastq.gz
wait
 
cd ../Brain/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/002/SRR1533652/SRR1533652_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/002/SRR1533652/SRR1533652_2.fastq.gz
wait

cd ../Bones/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/008/SRR1533658/SRR1533658_1.fastq.gz 
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/008/SRR1533658/SRR1533658_2.fastq.gz
wait

cd ../Gills/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/003/SRR1533653/SRR1533653_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/003/SRR1533653/SRR1533653_2.fastq.gz
wait

cd ../Heart/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/004/SRR1533654/SRR1533654_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/004/SRR1533654/SRR1533654_2.fastq.gz
wait

cd ../Muscle/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR1533655/SRR1533655_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/005/SRR1533655/SRR1533655_2.fastq.gz
wait

cd ../Kidney/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/007/SRR1533657/SRR1533657_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/007/SRR1533657/SRR1533657_2.fastq.gz
wait 

cd ../Intestine/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/009/SRR1533659/SRR1533659_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/009/SRR1533659/SRR1533659_2.fastq.gz
wait

cd ../Ovary/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/001/SRR1533651/SRR1533651_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/001/SRR1533651/SRR1533651_2.fastq.gz
wait

cd ../Embryo/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/000/SRR1533660/SRR1533660_1.fastq.gz
wait

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR153/000/SRR1533660/SRR1533660_2.fastq.gz




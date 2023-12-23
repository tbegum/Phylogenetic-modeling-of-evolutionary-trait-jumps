#!/usr/bin/env perl -w
## This program downloads the raw fastq (forward and backward) and (denovo assembly) transcript fasta file from ENSEMBL for each organism
## The raw fastq files are contaminated with adapter sequences, so we aim to filter it 
## and process it to get the abundance files against each tissue
use warnings;
use strict;

##Declaring variables
my$i;
my$outfile;
my$urlfile;
my$libraryfile;
my$index_file;
my$forward_read_file;
my%library;
my@fdline;

## Declearing the organism data we need from Phylofish
my@organism = ("medaka","cavefish","northern_pike");
 
## Making specific .sh files for each organism 
for($i=0;$i<@organism;$i++)
{
	chomp ($organism[$i]);
	my$file=$organism[$i];	
	$outfile="$organism[$i]_1.sh";
	$urlfile="$organism[$i]_url.txt";
	$libraryfile="$organism[$i]_library.txt";
	$index_file="$organism[$i]_index";
	$forward_read_file="$organism[$i]_forward.txt";

	##Printing the path of directories inside the .sh file
	## where the files should be downloaded
	open(FH,">$outfile");
	print FH "#!/bin/bash","\n\n\n";
	print FH '#BSUB -L /bin/bash',"\n",'#BSUB -o transcript',$i,'.txt',"\n",'#BSUB -e err',$i,'.txt',"\n",'#BSUB -J trans',$i,'_test',"\n",'#BSUB -u tunaiya@gmail.com',"\n",'#BSUB -R "rusage[mem=10000]"',"\n",'#BSUB -M 10194304',"\n\n";
	print FH 'cd /scratch/cluster/monthly/tina/',"\n";
	print FH "mkdir $organism[$i]/","\n";
	print FH "cd $organism[$i]/","\n\n"; 

	
	##Opening file containing the URL list for the organism
	open(FX1,"$urlfile");
	my@urls=<FX1>;

	##Printing the URLs to download the files
	foreach(@urls)
	{
		print FH "wget $_";
		
	}

	##To check the process is done
	print FH "\n\nwait\n";
	print FH 'echo "Downloading done..."',"\n\n";

	##Now we want to make library for all files so that the abundance file generates inside the specific tissue
	##So we need to open and read the library files
	open(FX2,"$libraryfile");
	my@tissuedata=<FX2>;

	##Now we need to add the trimmomatic module on cluster
	print FH "\n\nmodule add UHTS/Analysis/trimmomatic/0.36\n\n";
	
	## Next we need filter the raw reads from adapter sequence contamination
	## So we need to make a list of forward file
	foreach(@urls)
        {
		chomp $_;
		if($_=~/_R1.fastq.gz/)
		{
			my@read_file=split('/',$_);
        		my$readline=$read_file[-1];
			$readline=~s/_R1.fastq.gz//g;
			my$read1="$readline"."_R1.fastq.gz";
			my$read2="$readline"."_R2.fastq.gz";
			my$out1="$readline"."_R1_trim.fastq.gz";
			my$out2="$readline"."_R2_trim.fastq.gz";
			print FH "trimmomatic PE -threads 32 $read1 $read2 $out1 unpaired_1.fq $out2 unpaired_2.fq ILLUMINACLIP:/scratch/cluster/monthly/tina/adapter.fa:2:30:10 SLIDINGWINDOW:4:20  HEADCROP:12 MINLEN:25\n\n";
			print FH "\nwait\n";
			print FH "rm unpaired_1.fq\n";
			print FH "rm unpaired_2.fq\n";
			print FH "mv $out1 $read1\n";
			print FH "mv $out2 $read2\n";  
			print FH "\nwait\n\n";
		}
	}

	##Now adding module of Kallisto to work
	print FH "module add UHTS/Analysis/kallisto/0.43.1\n\n";

	## Now matching the cDNA file to make index file
	my$kline=$urls[0];
	my@kkallisto_indexfile=split('/',$kline);
	my$kallisto_indexfile=$kkallisto_indexfile[-1];
 
	##Making index file for each organism
	print FH "kallisto index -i $index_file $kallisto_indexfile\nwait\n";
	print FH 'echo "kallisto index is build for this organism"',"\n\n";   

	##Using the index file to make the abundance data against each tissue
	foreach(@tissuedata)
        {
                my@line=split("\t",$_);
		chop($line[0]);
		$line[3]=~s/\s//g;
		
		
		##Making a hash file for easy indexing when the "merge" is not there
		if($line[0]!~/merge/)
		{	
			$library{$line[0]}=$line[3];
			my$R1="$line[0]_R1.fastq.gz";
			my$R2="$line[0]_R2.fastq.gz";
                	print FH "kallisto quant -i $index_file -o $library{$line[0]} $R1 $R2";
			print FH "\nwait\n";
		}
		
		##for "merge" files the format is a bit different 
                ##Making a hash file for easy indexing when the "merge" is  there
                if($line[0]=~/merge/)
                {
			## Extracting tissue name
                        my$tissue_merge_1="$line[3]_1";
			my$tissue_merge_2="$line[3]_2";

			## Modifying the library name for proper taking data
			$line[0]=~s/merge_1/merge/g;
                        my$R11="$line[0]_1_R1.fastq.gz";
                        my$R12="$line[0]_1_R2.fastq.gz";
			my$R21="$line[0]_2_R1.fastq.gz";
                        my$R22="$line[0]_2_R2.fastq.gz";
	
                        print FH "kallisto quant -i $index_file -o $tissue_merge_1 $R11 $R12";
                        print FH "\nwait\n";
                        print FH "kallisto quant -i $index_file -o $tissue_merge_2 $R21 $R22";
                        print FH "\nwait\n";
                }
		
        }
	print FH 'echo "Finished abundance work...."',"\n\n\n";

	##Finally removing the .gz files
	print FH 'rm *.gz',"\n"; 

}

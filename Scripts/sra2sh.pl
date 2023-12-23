#!/usr/bin/env perl
use warnings;

# This program is written to generate .sh file for automatic download and unzip of data

open(FH,"SRA.txt");
@file=<FH>;
close FH;

open(FX,">geo_extrach.sh");
print FX "#!/usr/bin/bash\n\n";
print FX "mkdir dog\n";
print FX "mkdir rabbit\n";
print FX "mkdir ferret\n";


foreach $line(@file)
{
	chomp $line; ##To remove new line character at the end of a line
	@column=split("\t",$line);
        #print "$column[5]\t$column[8]\t$column[10]\n";
	if($column[8]=~/skeletal muscle/)
	{
		$column[8]="muscle";
	}
        if($column[5]=~/Canis lupus familiaris/)
	{
		$species="dog";
		print FX "\nmkdir $species/$column[8]\n";
		my$loc="$species/$column[8]/";
		print FX "cd $loc\n";
		if($column[8]=~/muscle/){$column[8]=skMuscle;}
                print FX "wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/",$column[10],"/suppl/",$column[10],"%5F",$column[8],"%2E",$species,"%2Egenes%2Eresults%2Etxt%2Egz";
		print FX "\nwait\n";
		print FX "gunzip *.gz\n";
		print FX "wait\n";
		print FX "cd ../../\n\n";
	}
	if($column[5]=~/Oryctolagus cuniculus/)
        {
		$species="rabbit";
                print FX "\nmkdir $species/$column[8]\n";
		my$loc="$species/$column[8]/";
                print FX "cd $loc\n";
		if($column[8]=~/muscle/){$column[8]=skMuscle;}
                print FX "wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/",$column[10],"/suppl/",$column[10],"%5F",$column[8],"%2E",$species,"%2Egenes%2Eresults%2Etxt%2Egz";		    
		print FX "\nwait\n";
                print FX "gunzip *.gz\n";
                print FX "wait\n";
		print FX "cd ../../\n\n";
        }
	if($column[5]=~/Mustela putorius furo/)
        {
		$species="ferret";
                print FX "\nmkdir $species/$column[8]\n";
                my$loc="$species/$column[8]/";
                print FX "cd $loc\n";
		if($column[8]=~/muscle/){$column[8]=skMuscle;}
                print FX "wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/",$column[10],"/suppl/",$column[10],"%5F",$column[8],"%2E",$species,"%2Egenes%2Eresults%2Etxt%2Egz";
                print FX "\nwait\n";
                print FX "gunzip *.gz\n";
                print FX "wait\n";
		print FX "cd ../../\n\n";
        }
 
}
	




#!/usr/bin/bash

mkdir dog
mkdir rabbit
mkdir ferret


mkdir dog/brain
cd dog/brain/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829087/suppl/GSM2829087%5Fbrain%2Edog%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir dog/liver
cd dog/liver/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829101/suppl/GSM2829101%5Fliver%2Edog%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir dog/muscle
cd dog/muscle/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829111/suppl/GSM2829111%5FskMuscle%2Edog%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir dog/testis
cd dog/testis/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829116/suppl/GSM2829116%5Ftestis%2Edog%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir dog/heart
cd dog/heart/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829091/suppl/GSM2829091%5Fheart%2Edog%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir dog/kidney
cd dog/kidney/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829096/suppl/GSM2829096%5Fkidney%2Edog%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir dog/lung
cd dog/lung/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829106/suppl/GSM2829106%5Flung%2Edog%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../

mkdir ferret/brain
cd ferret/brain/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829088/suppl/GSM2829088%5Fbrain%2Eferret%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir ferret/heart
cd ferret/heart/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829092/suppl/GSM2829092%5Fheart%2Eferret%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir ferret/kidney
cd ferret/kidney/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829097/suppl/GSM2829097%5Fkidney%2Eferret%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir ferret/lung
cd ferret/lung/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829107/suppl/GSM2829107%5Flung%2Eferret%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir ferret/liver
cd ferret/liver/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829102/suppl/GSM2829102%5Fliver%2Eferret%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir ferret/muscle
cd ferret/muscle/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829112/suppl/GSM2829112%5FskMuscle%2Eferret%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir ferret/testis
cd ferret/testis/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829117/suppl/GSM2829117%5Ftestis%2Eferret%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir rabbit/heart
cd rabbit/heart/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829090/suppl/GSM2829090%5Fheart%2Erabbit%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir rabbit/kidney
cd rabbit/kidney/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829095/suppl/GSM2829095%5Fkidney%2Erabbit%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir rabbit/liver
cd rabbit/liver/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829100/suppl/GSM2829100%5Fliver%2Erabbit%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir rabbit/lung
cd rabbit/lung/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829105/suppl/GSM2829105%5Flung%2Erabbit%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir rabbit/brain
cd rabbit/brain/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829086/suppl/GSM2829086%5Fbrain%2Erabbit%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir rabbit/muscle
cd rabbit/muscle/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829110/suppl/GSM2829110%5FskMuscle%2Erabbit%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


mkdir rabbit/testis
cd rabbit/testis/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2829nnn/GSM2829115/suppl/GSM2829115%5Ftestis%2Erabbit%2Egenes%2Eresults%2Etxt%2Egz
wait
gunzip *.gz
wait
cd ../../


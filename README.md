# Scaladbg
Parallel assembly algorithm building on top of IDBA. 

#Instructions

Get idba from either https://github.com/loneknightpy/idba or http://hku-idba.googlecode.com/files/idba-1.1.1.tar.gz

To build the ScalaDBG_SP version:

    cp scaladbg_sp_ud.cpp idba/src/release/idba_ud.cpp
    cd idba/
    # Now change build.sh: Change the line with ./configure to: ./configure CXX=mpic++ CXXFLAGS=-lrt
    ./build.sh   
    #./configure CXX=mpic++ CXXFLAGS=-lrt

To build the ScalaDBG_PP version:

    cp scaladbg_pp_ud.cpp idba/src/release/idba_ud.cpp
    cd idba/
    # Now change build.sh: Change the line with  ./configure to: ./configure CXX=mpic++ CXXFLAGS=-lrt
    ./build.sh   
    ./configure CXX=mpic++ CXXFLAGS=-lrt

The binary to use is idba_ud

Datasets were retrieved from the following:
Ecoli(real): http://spades.bioinf.spbau.ru/spades_test_datasets/ecoli_sc/
Ecoli(reference): Escherichia coli str. K-12 substr. MG1655: http://www.ncbi.nlm.nih.gov/nuccore/49175990?report=fasta
Simulated Ecoli: Used WGSIM and used the following command to generate the reads: ./wgsim -e .01 -1 75 -2 75 -d 250 -N 11909298 reference_ecoli.fa read1.fq read2.fq
Staphylococcus aureus(real): http://spades.bioinf.spbau.ru/spades_test_datasets/saureus/
Staphylococcus aureus(reference): http://www.ncbi.nlm.nih.gov/nuccore/?term=Staphylococcus+aureus+subsp.+aureus+USA300_FPR3757
SAR324(real): http://spades.bioinf.spbau.ru/spades_test_datasets/SAR324/reads/
Rhodobacter sphaeroides(real): http://gage.cbcb.umd.edu/data/Rhodobacter_sphaeroides/
Rhodobacter sphaeroides(reference): http://gage.cbcb.umd.edu/data/Rhodobacter_sphaeroides/
 

CAMI 

https://s3-eu-west-1.amazonaws.com/cami-data-eu/CAMI_low/RL_S001__insert_270.fq.gz


All .fastq files were converted using the following:

    ./idba/bin/fq2fa —-merge read1.fastq read2.fastq output.fa

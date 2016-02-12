#!/bin/bash

# install dependencies and libraries
sudo apt-get update && sudo apt-get -y upgrade && sudo apt-get -y install autoconf automake bison build-essential default-jdk default-jre expat fastqc fastx-toolkit  g++ gcc git libboost-all-dev libbz2-dev libncurses5-dev libpcre++-dev libpcre3-dev make parallel python-dev python-setuptools trimmomatic unzip wget zlib1g-dev

# install bioawk
sudo mkdir /sw 
sudo chown ubuntu /sw
chmod 775 /sw
cd /sw
git clone https://github.com/lh3/bioawk.git
cd bioawk
make

# install pullseq
cd /sw
git clone https://github.com/bcthomas/pullseq.git
cd pullseq
./bootstrap
./configure
make
sudo make install

# install trimmomatic
cd /sw && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip
unzip Trimmomatic-0.35.zip
echo 'trimmomatic=/sw/Trimmomatic-0.35/trimmomatic-0.35.jar' >> ~/.bashrc

#install samtools
cd /sw
wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
tar -xvjf samtools-1.2.tar.bz2
cd samtools-1.2
make

#install bowtie
cd /sw
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-src.zip
unzip bowtie-1.1.2-src.zip
cd bowtie-1.1.2
make

#set PATH variables
echo 'PATH=~/tnseq/scripts:$PATH' >> ~/.bashrc
echo 'PATH=/sw/bioawk:$PATH' >> ~/.bashrc
echo 'PATH=/sw/samtools-1.2:$PATH' >> ~/.bashrc
echo 'PATH=/sw/bowtie-1.1.2:$PATH' >> ~/.bashrc

#install python modules
sudo easy_install pip setuptools
sudo pip install --upgrade pip setuptools
sudo -H pip install pyopenssl ndg-httpsclient pyasn1
sudo -H pip install biopython



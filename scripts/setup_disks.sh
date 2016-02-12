#!/bin/bash

mkdir ~/data
sudo mkfs.ext4 -E nodiscard /dev/xvdc
sudo mount /dev/xvdc ~/data
sudo chown ubuntu ~/data
cd ~/data
wget http://dib-training.ucdavis.edu.s3.amazonaws.com/2016-bodega/tnseq_reads.fastq.gz
gunzip tnseq_reads.fastq.gz && md5sum tnseq_reads.fastq.gz

sudo umount /dev/xvdb
mkdir ~/analysis
sudo mount /dev/xvdb ~/analysis
sudo chown ubuntu ~/analysis

exit 0

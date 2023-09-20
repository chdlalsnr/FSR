#!/bin/bash
export WORK=/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/Current/5_Full_Dataset_Analysis

rm -rf $WORK/TASK/OUTPUT_$1
mkdir $WORK/TASK/OUTPUT_$1
cp $WORK/generateCSVfile.py $WORK/dataset.json $WORK/TASK/OUTPUT_$1
cd $WORK/TASK/OUTPUT_$1
python generateCSVfile.py $1

cd $WORK

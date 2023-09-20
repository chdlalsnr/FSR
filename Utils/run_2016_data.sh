#!/bin/bash
export WORK=/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Full_Dataset_Analysis

rm -rf $WORK/TASK_2016_data/OUTPUT_$1
mkdir $WORK/TASK_2016_data/OUTPUT_$1
cp $WORK/generateCSVfile_2016_data.py $WORK/TASK_2016_data/OUTPUT_$1
cd $WORK/TASK_2016_data/OUTPUT_$1
ln -s $WORK/dataset_2016_data .
python generateCSVfile_2016_data.py $1

cd $WORK

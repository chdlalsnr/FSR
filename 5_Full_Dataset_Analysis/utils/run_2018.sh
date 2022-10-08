#!/bin/bash
export WORK=/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/Present/5_Full_Dataset_Analysis

rm -rf $WORK/TASK_2018/OUTPUT_$1
mkdir $WORK/TASK_2018/OUTPUT_$1
cp $WORK/generateCSVfile_2018.py $WORK/TASK_2018/OUTPUT_$1
cd $WORK/TASK_2018/OUTPUT_$1
ln -s $WORK/dataset_2018 .
python generateCSVfile_2018.py $1

cd $WORK

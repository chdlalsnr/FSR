#!/bin/bash
export WORK=/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Make_Control_Plots

rm -rf $WORK/TASK_2017/OUTPUT_$1
mkdir $WORK/TASK_2017/OUTPUT_$1
cp $WORK/generateCSVfile_2017.py $WORK/TASK_2017/OUTPUT_$1
cd $WORK/TASK_2017/OUTPUT_$1
ln -s $WORK/dataset_2017 .
python generateCSVfile_2017.py $1

cd $WORK

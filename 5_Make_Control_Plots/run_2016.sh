#!/bin/bash
export WORK=/cms/ldap_home/chdlalsnr/HH/CMSSW_11_0_0/src/HHAnalysis/HH/python/5_Make_Control_Plots

rm -rf $WORK/TASK_2016/OUTPUT_$1
mkdir $WORK/TASK_2016/OUTPUT_$1
cp $WORK/generateCSVfile_2016.py $WORK/TASK_2016/OUTPUT_$1
cd $WORK/TASK_2016/OUTPUT_$1
ln -s $WORK/dataset_2016 .
python generateCSVfile_2016.py $1

cd $WORK

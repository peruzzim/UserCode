#!/bin/bash
cd /shome/peruzzi/shape_studies || exit 1
source /swshare/ROOT/thisroot.sh
if [ "$2" != "EEEE" ]
then 
root -q -b -l template_studies_2d_variablebinning.C+O\(\"$1\",\"$2\",$3,\"templateshapeMCpromptdrivenEB\"\)
fi
if [ "$2" != "EBEB" ]
then
root -q -b -l template_studies_2d_variablebinning.C+O\(\"$1\",\"$2\",$3,\"templateshapeMCpromptdrivenEE\"\)
fi



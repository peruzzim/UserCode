#!/bin/bash
cd /shome/peruzzi/shape_studies || exit 1
source /swshare/ROOT/thisroot.sh
root -q -b -l template_studies_2d_variablebinning.C+O\(\"outphoton_allmc_sigsig.root\",\"outphoton_allmc_1prcone1fgen.root\",\"outphoton_allmc_2fgen.root\",\"outphoton_allmc_standard.root\",\"$1\",\"$2\",$3,\"templateshapeMCpromptdriven\"\)




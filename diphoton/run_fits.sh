#!/bin/bash
cd /shome/peruzzi/shape_studies || exit 1
source /swshare/ROOT/thisroot.sh
root -q -b -l template_studies_2d_variablebinning.C+O\(\"outphoton_data_sigsig.root\",\"outphoton_data_sigbkg.root\",\"outphoton_data_bkgbkg.root\",\"outphoton_data_standard.root\",\"$1\",\"$2\",$3,\"$4\"\)

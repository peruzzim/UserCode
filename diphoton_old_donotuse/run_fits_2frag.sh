#!/bin/bash
cd /shome/peruzzi/shape_studies || exit 1
source /swshare/ROOT/thisroot.sh
root -q -b -l template_studies_2d_variablebinning.C+O\(\"$1\",\"$2\",$3,\"templateshape2frag\"\)



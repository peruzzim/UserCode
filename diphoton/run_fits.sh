#!/bin/bash
cd /Users/peruzzi/temp/UserCode/diphoton || exit 1
source /Users/peruzzi/root/root_5_34_08/bin/thisroot.sh
root -l template_studies_2d_variablebinning.C+O\(\"$1\",\"$2\",$3,\"$4\"\)




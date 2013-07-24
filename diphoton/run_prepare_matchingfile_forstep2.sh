#!/bin/bash
WORKDIR=/scratch/`whoami`/sgejob-${JOB_ID}-$1
source /swshare/ROOT/thisroot.sh
mkdir -p ${WORKDIR}
cd ${WORKDIR} || exit 1
root -q -b -l /shome/peruzzi/shape_studies/prepare_matchingfile_forstep2.C+O\(\"/scratch/peruzzi/matchingfile.root\",\"/scratch/peruzzi/Photon_Run2011AB_16Jan2012_v1_AOD.root\",$1\)
cd ${WORKDIR} || exit 1
mv matchingtree_*root /shome/peruzzi/shape_studies/matchingtrees/
rm -rf ${WORKDIR}

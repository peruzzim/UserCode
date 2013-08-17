#!/bin/bash
REMDIR=gg_minitree_020616_data2011_newtemplates_ago17
T3DCAPHOME=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/peruzzi/
WORKDIR=/scratch/`whoami`/sgejob-${JOB_ID}-$1
source /swshare/ROOT/thisroot.sh
mkdir -p ${WORKDIR}
cd ${WORKDIR} || exit 1
root -q -b -l /shome/peruzzi/shape_studies/prepare_matchingfile_forstep2.C+O\(\"/scratch/peruzzi/matchingfile_data_ago17.root\",\"/scratch/peruzzi/Photon_Run2011AB_16Jan2012_v1_AOD_ago17.root\",$1\)
cd ${WORKDIR} || exit 1
for i in matchingtree_*root; do lcg-cp -b -D srmv2 file:///`pwd`/${i} ${T3DCAPHOME}/${REMDIR}/matchingtrees/${i}; done
rm -rf ${WORKDIR}

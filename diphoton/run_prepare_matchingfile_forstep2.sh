#!/bin/bash
WORKDIR=/scratch/`whoami`/sgejob-${JOB_ID}-$1
source /swshare/ROOT/thisroot.sh
mkdir -p ${WORKDIR}
cd ${WORKDIR} || exit 1
#root -q -b -l /shome/peruzzi/shape_studies/prepare_matchingfile_forstep2.C+O\(\"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/peruzzi/gg_minitree_020616_data2011_newtemplates_jul23/matchingfile.root\",\"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/peruzzi/gg_minitree_020616_data2011_newtemplates_jul23/Photon_Run2011AB_16Jan2012_v1_AOD.root\",$1\)
#root -q -b -l /shome/peruzzi/shape_studies/prepare_matchingfile_forstep2.C+O\(\"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/peruzzi/gg_minitree_020616_data2011_newtemplates_jul23/matchingfile.root\",\"/shome/peruzzi/shape_studies/Photon_Run2011AB_16Jan2012_v1_AOD.root\",$1\)
root -q -b -l /shome/peruzzi/shape_studies/prepare_matchingfile_forstep2.C+O\(\"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/peruzzi/gg_minitree_020616_data2011_newtemplates_jul23/matchingfile.root\",\"/scratch/peruzzi/Photon_Run2011AB_16Jan2012_v1_AOD.root\",$1\)
cd ${WORKDIR} || exit 1
mv matchingtree_*root /shome/peruzzi/shape_studies/matchingtrees/
rm -rf ${WORKDIR}

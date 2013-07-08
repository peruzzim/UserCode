cd /shome/peruzzi/shape_studies || exit 1
for i in invmass; do for j in EBEB EBEE EEEE; do for k in `seq 0 14`; do qsub -o /dev/null -e /dev/null  -q short.q run_fits.sh $i $j $k; done ; done; done;
for i in diphotonpt; do for j in EBEB EBEE EEEE; do for k in `seq 0 19`; do qsub -o /dev/null -e /dev/null  -q short.q run_fits.sh $i $j $k; done ; done; done;
for i in costhetastar; do for j in EBEB EBEE EEEE; do for k in `seq 0 6`; do qsub -o /dev/null -e /dev/null  -q short.q run_fits.sh $i $j $k; done ; done; done;
for i in dphi; do for j in EBEB EBEE EEEE; do for k in `seq 0 12`; do qsub -o /dev/null -e /dev/null  -q short.q run_fits.sh $i $j $k; done ; done; done;
for i in dR; do for j in EBEB EBEE EEEE; do for k in `seq 0 6`; do qsub -o /dev/null -e /dev/null  -q short.q run_fits.sh $i $j $k; done ; done; done;

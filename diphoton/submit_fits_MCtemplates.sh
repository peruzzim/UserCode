cd /shome/peruzzi/shape_studies || exit 1

#for i in invmass; do for j in EBEB EBEE EEEE; do qsub -o /dev/null -e /dev/null -q long.q run_fits_MCtrue.sh $i $j 26; done ; done; 
for i in invmass; do for j in EBEB EBEE EEEE; do qsub -o /dev/null -e /dev/null -q long.q run_fits_MCpromptdriven.sh $i $j 26; done ; done; 
for i in invmass; do for j in EBEB EBEE EEEE; do qsub -o /dev/null -e /dev/null -q long.q run_fits_MCfakedriven.sh $i $j 26; done ; done; 


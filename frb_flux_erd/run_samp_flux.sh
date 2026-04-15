export OMP_NUM_THREADS=2
#fgtype=ALG_YMW16
for fgtype in ALG_NE2001 ALG_YMW16
do
    echo "mpiexec -n 4 ./nest_samp_flux.py -fc chime_frb_cat.txt -fs chime_tel_svy.txt -g ${fgtype} -upper 1 -o flux_${fgtype}_upper"
    mpiexec -n 4 ./nest_samp_flux.py -fc chime_frb_cat.txt -fs chime_tel_svy.txt -g ${fgtype} -upper 1 -o flux_${fgtype}_upper 
    #echo "mpiexec -n 4 ./nest_samp_flux.py -fc chime_frb_cat.txt -fs chime_tel_svy.txt -g ${fgtype} -upper 1 halo 1 -o flux_${fgtype}_upper_halo"
    #mpiexec -n 4 ./nest_samp_flux.py -fc chime_frb_cat.txt -fs chime_tel_svy.txt -g ${fgtype} -upper 1 -halo 1 -o flux_${fgtype}_upper_halo
done

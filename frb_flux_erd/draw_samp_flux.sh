#for fgtype in ETG_NE2001_upper ETG_YMW16_upper LTG_NE2001_upper LTG_YMW16_upper ALG_NE2001_upper ALG_YMW16_upper
#for fgtype in ETG_NE2001_upper_halo ETG_YMW16_upper_halo LTG_NE2001_upper_halo LTG_YMW16_upper_halo ALG_NE2001_upper_halo ALG_YMW16_upper_halo

for fgtype in ALG_NE2001 ALG_YMW16
do
    prefix="flux_${fgtype}_upper"
    echo "New plot of ${prefix} has been finished."
    ./pltpost.py -f ./nest_out/samp/${prefix} -o ./plots/samp/${prefix}.eps -title "Real Sample (Flux)" -up 1 -bo 1
done

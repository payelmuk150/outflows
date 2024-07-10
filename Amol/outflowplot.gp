reset
#set term epslatex standalone color colortext dashed 'cmr,12' header "\\usepackage{amsmath,stackengine}\n\\usepackage{grffile}"
set terminal pdfcairo enhanced color font 'cmr,12'
set grid

set linetype 1 lc rgb "dark-violet" lw 2 pt 0
set linetype 2 lc rgb "dark-green"  lw 2 pt 7
set linetype 3 lc rgb "#00bbbb"     lw 2 pt 6 pi -1 #instead of "cyan" I chose a 73.3...%-intensity, yet still fully saturated cyan
set linetype 4 lc rgb "dark-red"    lw 2 pt 5 pi -1
set linetype 5 lc rgb "blue"        lw 2 pt 8
set linetype 6 lc rgb "goldenrod"   lw 2 pt 5
set linetype 7 lc rgb "magenta"     lw 2 pt 11
set linetype 8 lc rgb "dark-orange" lw 2
set linetype 9 lc rgb "#d55e00"     lw 2 #vermillion
set linetype 10 lc rgb "#cc79a7"    lw 2 #reddish purple

set key noenhanced
set logscale xy
set format xy "10^{%T}"


if(1) { 

strprog = "s12.75"
strsup = strprog."_Bollig_3D_lum_mix_rns17"
strsub = strprog."_Bollig_3D_lum"

file_new = "Nucleo_Amol_".strsup.".txt"
file_sub = "Nucleo_Amol_".strsub.".txt"

print file_new


set xlabel "Time (s)"
set ylabel "Temperature (MeV)"
set output "comparison_Tvst_".strsup.".pdf"
plot file_new using 1:3 w l title strsup 


set xlabel "Time (s)"
set ylabel "Distance (km)"
set output "comparison_rvst_".strsup.".pdf"
plot file_new using 1:($2/1.0e5) w l title strsup, file_sub using 1:($2/1.0e5) w l title strsub


set xlabel "Distance (km)"
set ylabel "Temperature (MeV)"
set output "comparison_Tvsr_".strsup.".pdf"
plot file_new using (($2)/(1.0e5)):3 w l title strsup, file_sub using (($2)/(1.0e5)):3 w l title strsub


set xrange [*:1e6]
set yrange [*:*]
set xlabel "Distance (km)"
set ylabel "Outflow velocity (cm/s)"
set output "vvsr_expansion_".strsup.".pdf"
plot file_new u ($2/1.0e5):5 w l title "Outflow velocity"


file_vtest = "vin_tests_".strsup.".txt"

#unset logscale y
#set format y "%g"

set xrange [*:1e4]
set yrange [1e6:1e10]
set xlabel "Distance (km)"
set ylabel "Outflow velocity and sound speed (cm/s)"
set output "vvsr_".strsup.".pdf"
plot file_vtest u ($1/1.0e5):4 w l title "Outflow velocity" , file_vtest u ($1/1.0e5):5 w l title "Sound speed"

set xlabel "Distance (km)"
set xrange [*:*]
set yrange [*:*]
unset ylabel
set output "Tf_error_".strsup.".pdf"
plot file_vtest u ($1/1.0e5):6 w l title "Tf_error = GM/(3v_s^2r)"

}




##################################### For testing only ####################################


if(0) { # enable if testing

file_test = "testout_noterm.txt"
file_ref = "Bollig_3D_rs12000/Alex_checks_Amol_s12.75_Bollig_3D_lum_mix.txt"

set logscale x

set ylabel "Velocity (cm/s)"
set xlabel "Distance (km)"

set xrange [*:1e6]

set output "test_vvsr_noterm.pdf"
plot file_test u ($2/1.0e5):5 w l title "transonic", file_ref u ($2/1.0e5):5 w l title "critical"

}


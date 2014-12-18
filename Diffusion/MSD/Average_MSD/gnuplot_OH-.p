set ytics font "Times-Roman, 20"
set xtics font "Times-Roman, 20"
set ylabel "MSD (A^2/ps)" font "TimesNewRoman, 20"
set xlabel "time (ps)" font "TimesNewRoman, 20"
set key font "TimesNewRoman, 15"
set key center top
plot "OH-_PBE.msd" u 1:2 w l title "OH- PBE" , \
     "OH-_PBE_vdW.msd" u 1:2 w l title "OH- PBE vdW" , \
     "OH-_PBE0_vdW.msd" u 1:2 w l title "OH- PBE0 vdW", \
     "H+_PBE.msd" u 1:2 w l title "H+ PBE" , \
     "H+_PBE_vdW.msd" u 1:2 w l title "H+ PBE vdW" , \
     "H+_PBE0_vdW.msd" u 1:2 w l title "H+ PBE0 vdW"

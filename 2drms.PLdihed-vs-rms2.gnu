set size square
set pm3d map corners2color c1
set xlabel "Clusters from RMS :2-61@P"
set ylabel "Clusters from dihedrals :13-20"
set yrange [   1.000:   6.000]
set xrange [   1.000:   6.000]
set cbrange [   0.0:    16.0]
splot "-" with pm3d title "2drms.dihed-vs-rms2.gnu"
   1.000    1.000    3.400
   1.000    2.000   12.133
   1.000    3.000    9.159
   1.000    4.000    8.852
   1.000    5.000   16.455
   1.000    6.000 0

   2.000    1.000    9.886
   2.000    2.000    2.590
   2.000    3.000   12.867
   2.000    4.000   11.080
   2.000    5.000    9.576
   2.000    6.000 0

   3.000    1.000   10.173
   3.000    2.000   11.486
   3.000    3.000    4.014
   3.000    4.000    3.357
   3.000    5.000   14.973
   3.000    6.000 0

   4.000    1.000    9.043
   4.000    2.000   12.124
   4.000    3.000    6.148
   4.000    4.000    7.331
   4.000    5.000   14.813
   4.000    6.000 0

   5.000    1.000   13.605
   5.000    2.000    6.090
   5.000    3.000   14.799
   5.000    4.000   13.402
   5.000    5.000    3.072
   5.000    6.000 0

   6.000    1.000 0
   6.000    2.000 0
   6.000    3.000 0
   6.000    4.000 0
   6.000    5.000 0
   6.000    6.000 0

end
pause -1

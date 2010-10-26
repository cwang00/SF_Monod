reset
#set origin 0,0
#set size 0.9, 0.9
#set terminal postscript eps enhanced dashed lw 2 "Helvetica" 16 
set terminal X11
#set size 2, 2
#set origin 0,0
#set multiplot layout 3,2 columnsfirst scale 2,2
set style data linespoints
#set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0
#set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0
#set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0
set view map
#set tic scale 0
set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set xrange [ 0.500000 : 50.50000 ] noreverse nowriteback
set yrange [ 0.500000 : 100.50000 ] noreverse nowriteback
set zrange [ * : * ] noreverse nowriteback  # (currently [0.00000:5.00000] )
set cblabel "Conc mg/l" 
#set cbrange [ 0.00000 : 5.00000 ] noreverse nowriteback
#set palette rgbformulae -7, 2, -7
set title "Ca, t = 1"
splot './more_particles/conc.1.00000001.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 1, Press enter to continue "

set title "Ca, t = 2"
splot './more_particles/conc.1.00000002.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 2,Press enter to continue"

set title "Ca, t = 3"
splot './more_particles/conc.1.00000003.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 3, Press enter to continue"

set title "Ca, t = 4"
splot './more_particles/conc.1.00000004.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 4, Press enter to continue"

set title "Ca, t = 5"
splot './more_particles/conc.1.00000005.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 5, Press enter to continue"

set title "Ca, t = 6"
splot './more_particles/conc.1.00000006.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 6, Press enter to continue"

set title "Ca, t = 7"
splot './more_particles/conc.1.00000007.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 7, Press enter to continue"

set title "Ca, t = 8"
splot './more_particles/conc.1.00000008.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 8, Press enter to continue"

set title "Ca, t = 9"
splot './more_particles/conc.1.00000009.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 9, Press enter to continue"

set title "Ca, t = 10"
splot './more_particles/conc.1.00000010.gnuplot' using (1+$1):($2+1):3 index 67 matrix with image
pause -1 "Ca, t = 10, Press enter to continue"
#set output "conc.eps"
#set output "conc.eps"
#replot

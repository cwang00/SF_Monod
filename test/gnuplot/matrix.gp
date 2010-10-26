reset
set style data linespoints
set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0
set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0
set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0
set view map
set tic scale 0
set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set xrange [ 0.500000 : 3.50000 ] noreverse nowriteback
set yrange [ 0.500000 : 3.50000 ] noreverse nowriteback
set zrange [ * : * ] noreverse nowriteback  # (currently [0.00000:5.00000] )
set cblabel "Conc mg/l" 
#set cbrange [ 0.00000 : 5.00000 ] noreverse nowriteback
#set palette rgbformulae -7, 2, -7
splot 'martix.txt' using (1+$1):($2+1):3 index 2 matrix with image

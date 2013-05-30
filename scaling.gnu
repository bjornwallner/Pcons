set terminal postscript color
set output 'scaling.ps'
set xlabel 'threads'
set ylabel 'time (s)'
set pointsize 2
plot 'scaling.dat' u 1:2 w lp pt 5
set yrange [0:5] 
plot 'scaling.dat' u 1:2 w lp pt 5
set xrange [0:150]
plot 'scaling.dat' u 1:2 w lp pt 5
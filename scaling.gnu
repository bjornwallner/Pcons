set terminal postscript color
set output 'scaling.ps'
set xlabel 'threads'
set ylabel 'time (s)'
plot 'scaling.dat' u 1:2 w lp 
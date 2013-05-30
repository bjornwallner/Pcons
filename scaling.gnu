set terminal postscript color solid
set output 'scaling.ps'
set xlabel 'threads'
set ylabel 'time (s)'
set pointsize 1.5
set xrange [0:20]
plot 'scaling.dat.no.schedule' u 1:2 w lp pt 5 t 'no schedule','scaling.dat' u 1:2 w lp pt 6 t 'schedule(dynamic)'
set xrange [*:*]
plot 'scaling.dat.no.schedule' u 1:2 w lp pt 5 t 'no schedule','scaling.dat' u 1:2 w lp pt 6 t 'schedule(dynamic)'
set yrange [0:5] 
plot 'scaling.dat.no.schedule' u 1:2 w lp pt 5 t 'no schedule','scaling.dat' u 1:2 w lp pt 6 t 'schedule(dynamic)'
set xrange [0:150]
plot 'scaling.dat.no.schedule' u 1:2 w lp pt 5 t 'no schedule','scaling.dat' u 1:2 w lp pt 6 t 'schedule(dynamic)'

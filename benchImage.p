# Gnuplot script

set terminal png
set output "avril.png"

set title "Parallel N-body implementation for time complexity N²"

set key right center

set xlabel "Number of threads"
set ylabel "Time in s"

set xrange [0:16]
set yrange [0:12]

plot "benchSeqNN.dat" u 1:2 w linespoints title "Sequential N²", \
     "benchParNN.dat" u 1:2 w linespoints title "Parallel N²"

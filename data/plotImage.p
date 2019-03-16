# Gnuplot script

set terminal png
set output "avril.png"

set title "N-body implementation for time complexity N²"

set key right center

set xlabel "Number of threads"
set ylabel "Time in s"

set xrange [0:16]
set yrange [0:70]

plot "bench120.dat" u 1:2 w linespoints title "Parallel N² 120", \
     "bench180.dat" u 1:2 w linespoints title "Parallel N² 180", \
     "bench240.dat" u 1:2 w linespoints title "Parallel N² 240", \
     "sbench120.dat" u 1:2 w linespoints title "Sequential N² 120", \
     "sbench180.dat" u 1:2 w linespoints title "Sequential N² 180", \
     "sbench240.dat" u 1:2 w linespoints title "Sequential N² 240"

# Gnuplot script

set terminal png
set output "movement.png"

set xrange [0:20 < * < 1000]
set yrange [0:20 < * < 1000]

plot("data.dat") with points palette

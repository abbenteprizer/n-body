# Gnuplot script

set terminal png
set output "movementPar.png"
#set output "images/movement.png"

set xrange [-1000:1000 < * < 1000]
set yrange [-1000:1000 < * < 1000]

plot("dataPar.dat") with points palette
#do for [t=0:10] {
#		outfile = sprintf('images/movement_%03.0f.png',t)
#		set output outfile
#		dataline = sprintf('<(sed -n '%d, %dp' data.dat)', t, t)

#		plot dataline with points palette
#}

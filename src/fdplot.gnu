unset key
set grid
set terminal png
set output "test.png"
set xrange[0:1]
set yrange[-2:2]

filename(n) = sprintf("data/sol%d.dat",i)
plot for [i=0:5] filename(i) using 1:2 with lines

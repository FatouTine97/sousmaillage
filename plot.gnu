set terminal png
set output 'sous-maillage.png'
set title "TFD - Courbe 3"
set xlabel "t"
set ylabel "Eres"
plot 'ez_s.dat' using 1:2 with linespoints title "H"

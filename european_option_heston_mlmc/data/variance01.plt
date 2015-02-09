set termoption enhanced
set termoption dashed
set xtics 0, 1, 4
set xlabel 'l'
set ylabel 'log_M variance'
set xrange [0:4]
plot "variance01.txt" using 1:(log($2)/log(4)) with lp pt 1 lt 3 lw 1 lc 1title "EM"
replot "variance01.txt" using 1:(log($3)/log(4)) with lp pt 2 lt 3 lw 1 lc 1title "NV"
replot "variance01.txt" using 1:(log($4)/log(4)) with lp pt 3 lt 3 lw 1 lc 1title "NN"
set terminal eps dashed enhanced monochrome
set output "variance01.eps"
replot
set terminal aqua
replot


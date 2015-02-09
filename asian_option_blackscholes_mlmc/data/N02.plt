set termoption enhanced
set termoption dashed
set xtics 0, 1, 6
set logscale y
set format y "10^{%L}"
set xlabel 'l'
set ylabel 'N_l'
set xrange [0:6]
set yrange [1:10000000000.0]
plot "N02.txt" using 1:($2) with lp pt 1 lt 3 lw 1 lc 1 title "{/Symbol e}=0.00005"
replot "N02.txt" using 1:($3) with lp pt 2 lt 3 lw 1 lc 1 title "{/Symbol e}=0.0001"
replot "N02.txt" using 1:($4) with lp pt 3 lt 3 lw 1 lc 1 title "{/Symbol e}=0.0002"
replot "N02.txt" using 1:($5) with lp pt 4 lt 3 lw 1 lc 1 title "{/Symbol e}=0.0005"
replot "N02.txt" using 1:($6) with lp pt 5 lt 3 lw 1 lc 1 title "{/Symbol e}=0.001"
set terminal eps dashed enhanced monochrome
set output "N02.eps"
replot
set terminal aqua
replot


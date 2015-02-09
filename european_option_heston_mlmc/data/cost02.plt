set termoption enhanced
set termoption dashed
set xlabel '{/Symbol e}'
set ylabel '{/Symbol e}^2 Ccost'
plot "cost02.txt" using 1:($2) with lp pt 1 lt 3 lw 1 lc 1 title "EM"
replot "cost02.txt" using 1:($3) with lp pt 2 lt 3 lw 1 lc 1 title "NV"
replot "cost02.txt" using 1:($4) with lp pt 3 lt 3 lw 1 lc 1 title "NN"
set terminal eps dashed enhanced monochrome
#set terminal pngcairo dashed enhanced monochrome
set output "cost02.eps"
replot
set terminal aqua
replot


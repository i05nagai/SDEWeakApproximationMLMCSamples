set termoption enhanced
set termoption dashed
set xlabel '{/Symbol e}'
set ylabel 'error[%]'
plot "P01.txt" using 1:(100*(1.099826125805e-01-$2)/1.099826125805e-01) with lp pt 1 lt 3 lw 1 lc 1 title "EM"
replot "P01.txt" using 1:(100*(1.099826125805e-01-$3)/1.099826125805e-01) with lp pt 2 lt 3 lw 1 lc 1 title "NV"
replot "P01.txt" using 1:(100*(1.099826125805e-01-$4)/1.099826125805e-01) with lp pt 3 lt 3 lw 1 lc 1 title "NN"
set terminal eps dashed enhanced monochrome
#set terminal pngcairo dashed enhanced monochrome
set output "P01.eps"
replot
set terminal aqua
replot


set terminal pngcairo enhanced font "times new roman,35" size 600,600 transparent

set size square

set key bottom right samplen 0.5 font "times new roman,30"

# set output 'Plot_Egel_R.png'
set xtic 1 offset 0,0.25
set xr [0:5]
set mxtics 2
set xl "{/:Italic E}_{gel} (kPa)" offset 0,0.750
set yl "{/:Italic R}_{lum}/{/:Italic R}^{ini}_{cell}"

set arrow 1 from 3.5,0.4 to 3.5,1.8 nohead lc rgb 'gray20' dt 2 lw 2

# plot 'E_gel.dat'  u 1:2  w l lc "red" lw 3  t ''

unset arrow 1

# set output 'Plot_Egel_Acell.png'
set ytic 0.2 offset 0.25,0
set yl "{/:Italic A}_{cell}/{/:Italic A}^{ini}_{cell}"
# plot 'E_gel.dat'  u 1:3  w l lc "red" lw 3  t ''

# set output 'Plot_Egel_Lcleft.png'
set ytic 0.1 offset 0.25,0
set yl "{/:Italic L}_{cleft}/{/:Italic R}^{ini}_{cell}"
# plot 'E_gel.dat'  u 1:6  w l lc "red" lw 3  t ''

# set output 'Plot_Egel_P.png'
set ytic 100 offset 0.25,0
set yr [0.0:]
set yl "{/:Italic P}-{/:Italic P}_{ext} (Pa)"
# plot 'E_gel.dat'  u 1:4  w l lc "orange" lw 3  t 'cell','E_gel.dat'  u 1:5  w l lc "dark-magenta" lw 3  t 'lumen'

set xl "time (h)" offset 0,1.0

set xtic 12 offset 0,0.25
set mxtics 6
set mytics 2
set ytic offset 0.25,0

set ytic 0.5 offset 0.25,0
set yr [0:2]
unset xr
unset yr

set output 'Plot1.png'
#set yl "c_{cell}/c_{ext}"
set xr [0:28]
set yl "{/:Italic R}_{lum}/{/:Italic R}^{ini}_{cell}"  offset 1,0
plot 'Data_c1.dat'  u ($1/60/60):7  w l lc "red" lw 3  t 'cell 1',\
'Data_c2.dat'  u ($1/60/60):7  w l lc "blue" lw 3 dt "-"  t 'cell 2'

#set yr [1:1.03]
set output 'Plot2.png'
set key top right
 set ytic 0.05 offset 0.25,0
set yr [0.0:0.2]
set yl "{/:Italic c}-{/:Italic c}_{ext} (mM)"
plot 'Data_c1.dat'  u ($1/60/60):(($2-1)*300)  w l lc "orange" lw 3  t 'cell',\
'Data_c1.dat'  u ($1/60/60):(($3-1)*300)  w l lc "dark-magenta" lw 3  t 'lumen',\
'Data_c2.dat'  u ($1/60/60):(($2-1)*300)  w l lc "orange" lw 3 dt "-" t '',\
'Data_c2.dat'  u ($1/60/60):(($3-1)*300)  w l lc "dark-magenta" lw 3 dt "-" t ''

unset yr
set ytic 0.2 offset 0.25,0
set output 'Plot4.png'
set key top right
set yl "{/:Italic A}_{cell}/{/:Italic A}^{ini}_{cell}"
plot 'Data_c1.dat'  u ($1/60/60):4  w l lc "red" lw 3  t 'cell1 ',\
'Data_c2.dat'  u ($1/60/60):4  w l lc "blue" lw 3 dt "-"  t 'cell 2'


set output 'Plot5.png'

set yl "{/:Italic L}_{cleft}/{/:Italic R}^{ini}_{cell}"
plot 'Data_c1.dat'  u ($1/60/60):8  w l lc "red" lw 3  t 'cell 1',\
'Data_c2.dat'  u ($1/60/60):8  w l lc "blue" lw 3 dt "-"  t 'cell 2'


set output 'Plot6.png'
set key bottom right
set ytic 100 offset 0.25,0
set yr [0.0:]
set yl "{/:Italic P}-{/:Italic P}_{ext} (Pa)"
plot 'Data_c1.dat'  u ($1/60/60):(($9-1)*101325)  w l lc "orange" lw 3  t 'cell',\
'Data_c1.dat'  u ($1/60/60):(($10-1)*101325)  w l lc "dark-magenta" lw 3  t 'lumen',\
'Data_c2.dat'  u ($1/60/60):(($9-1)*101325)  w l lc "orange" lw 3 dt "-" t '',\
'Data_c2.dat'  u ($1/60/60):(($10-1)*101325)  w l lc "dark-magenta" lw 3 dt "-" t ''

set output 'Plot7a.png'
unset yr
unset ytics
set ytics
set yl "∆{/:Italic P}_{apical} (kPa)"
plot 'Data_c1.dat'  u ($1/60/60):(($9-$10)*101.325)  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):(($9-$10)*101.325)  w l lc "blue" lw 3 dt "-"  t 'c2'
set output 'Plot7b.png'
set yl "∆{/:Italic P}_{basal} (kPa)"
plot 'Data_c1.dat'  u ($1/60/60):(($9-1.0)*101.325)  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):(($9-1.0)*101.325)  w l lc "blue" lw 3 dt "-"  t 'c2'

set output 'Plot8.png'
set yl "A_{lum}/A_{cell}^{ini}"
plot 'Data_c1.dat'  u ($1/60/60):6  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):6  w l lc "blue" lw 3 dt "-"  t 'c2'

set output 'Plot9.png'
set yl "L_{a}/L_{a}^{ini}"
plot 'Data_c1.dat'  u ($1/60/60):11  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):11  w l lc "blue" lw 3 dt "-"  t 'c2'

set output 'Plot10.png'
unset xr
unset yr
set yr
# reset ytic
set ytic 0.1
set yl "L_{b}/L_{b}^{ini}"
plot 'Data_c1.dat'  u ($1/60/60):12  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):12  w l lc "blue" lw 3 dt "-"  t 'c2'


#set yr [-0.01:0.01]

set output 'Plot11.png'
set yl "dW_{cell}"
plot 'Data_c1.dat'  u ($1/60/60):14  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):14  w l lc "blue" lw 3 dt "-"  t 'c2'

set output 'Plot12.png'
set yl "dW_{lum}"
plot 'Data_c1.dat'  u ($1/60/60):15  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):15  w l lc "blue" lw 3 dt "-"  t 'c2'

set output 'Plot_Rb.png'
unset yr
set ytics 0.1
set yl "R_{b}"
plot 'Data_c1.dat'  u ($1/60/60):25  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):25  w l lc "blue" lw 3 dt "-"  t 'c2'

unset yr
set ytics
set ytics (0)
#set yr [-7*10**-13:7*10**12]
set key top right maxrows 2 width -4.5
set output 'Plot_Jw.png'
set yl ""
plot 'Data_c1.dat'  u ($1/60/60):16  w l lc "red" lw 3 t '{/:Italic J}^{ w}_{basal}',\
'Data_c1.dat'  u ($1/60/60):17  w l lc "blue" lw 3  t '{/:Italic J}^{ w}_{apical}',\
'Data_c1.dat'  u ($1/60/60):27  w l lc "dark-magenta" lw 3  t '{/:Italic J}^{ w}_{cleft}',\
'Data_c1.dat'  u ($1/60/60):18  w l lc "web-green" lw 3  t '{/:Italic J}^{ w}_{leak}',\
'Data_c1.dat'  u ($1/60/60):(0*$16)  w l lc "black" lw 0.5  t '',\
'Data_c2.dat'  u ($1/60/60):16  w l lc "orange-red" lw 3 dt "_" t '',\
'Data_c2.dat'  u ($1/60/60):17  w l lc "royalblue" lw 3 dt "_" t '',\
'Data_c2.dat'  u ($1/60/60):27  w l lc "magenta" lw 3 dt "_" t '',\
'Data_c1.dat'  u ($1/60/60):18  w l lc "web-green" lw 3 dt "_" t ''

set output 'Plot_Ji.png'
set key top center maxrows 2 width -4
set yr [0:6*10**-12]
set yl ""
plot 'Data_c1.dat'  u ($1/60/60):22  w l lc "red" lw 3  t '{/:Italic J}^{ i}_{basal}',\
'Data_c1.dat'  u ($1/60/60):23  w l lc "blue" lw 3  t '{/:Italic J}^{ i}_{apical}',\
'Data_c1.dat'  u ($1/60/60):28  w l lc "dark-magenta" lw 3  t '{/:Italic J}^{ i}_{cleft}',\
'Data_c1.dat'  u ($1/60/60):24  w l lc "web-green" lw 3  t '{/:Italic J}^{ i}_{leak}',\
'Data_c1.dat'  u ($1/60/60):(0*$16)  w l lc "black" lw 0.5  t '',\
'Data_c2.dat'  u ($1/60/60):22  w l lc "orange-red" lw 3 dt "_" t '',\
'Data_c2.dat'  u ($1/60/60):23  w l lc "royalblue" lw 3 dt "_"  t '',\
'Data_c2.dat'  u ($1/60/60):28  w l lc "magenta" lw 3 dt "_" t '',\
'Data_c2.dat'  u ($1/60/60):24  w l lc "web-green" lw 3 dt "_"  t ''

set key bottom right

set yr
set ytics

set ytics 5

set output 'Plot16.png'
set yl "Sectoral angle"
plot 'Data_c1.dat'  u ($1/60/60):($21*180/pi)  w l lc "red" lw 3  t 'c1','Data_c2.dat'  u ($1/60/60):($21*180/pi)  w l lc "blue" lw 3 dt "-"  t 'c2'

r=10;

set output 'Plot17.png'
set grid
set logscale xy
set key bottom right
set xr [1:200]
set yr [1:1000]
set xtic 10
set xl "Lumen radius (μm)"
set yl "Number of cells"
set ytic 10
set mytics 10
set mxtics 10
s=100;
plot 'Experiment.dat'  u 2:4 w p pt 7 lc "skyblue" t '',\
'< tail -n 1 Data_c1.dat'  u ($7*r):(4*pi*(($7*r)**2))/(pi*($26/2*r)**2)  w p pt 4 lc "red" lw 3  t 'c1',\
'< tail -n 1 Data_c2.dat'  u ($7*r):(4*pi*(($7*r)**2))/(pi*($26/2*r)**2)  w p pt 4 lc "blue" lw 3   t 'c2'

set output 'Plot18.png'
set grid
set logscale xy
set key top left
set xr [0.1:100]
set yr [0.1:1000]
set xtic 10
set xl "Lumen radius (μm)"
set yl "Cell apical area (μm^2)"
set ytic 10
set mytics 10
plot 'Experiment.dat'  u 2:5 w p pt 7 lc "skyblue" t '', '< tail -n 1 Data_c1.dat'  u ($7*r):(pi*($26/2*r)**2)  w p pt 4 lc "red" lw 3  t '','< tail -n 1 Data_c2.dat'  u ($7*r):(pi*($26/2*r)**2)  w p pt 4 lc "blue" dt "-" lw 3  t ''

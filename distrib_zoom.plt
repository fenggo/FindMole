set terminal png truecolor enhanced size 1280,960 font arial 30 linewidth 4 #transparent
set out 'distrib_b.png'


set xlabel 'Time step(fs)'
set ylabel 'Molecular structures populations'
#set xtics ("{/Symbol G}" 0.0, "F"  0.5, "Q" 1.0, "Z'" 1.5, "{/Symbol G}" 2.0)

set xr [0:500]
set yr [0:3]

#set key at 3.0,8.2

#set arrow 1 from 0.7071,0  to 0.7071,8.5   nohead lt 9
#set arrow 2 from 1.06066,0 to 1.06066,8.5  nohead lt 9
#set arrow 3 from 1.41421,0 to 1.41421,8.5  nohead lt 9
#set arrow 4 from 2.63896,0 to 2.63896,8.5  nohead lt 9

#pt 7 red points
#plot 'hbond_986.923.txt'      u 1:2 w p pt 7 ps 2 t 'CH3-NO2'  
#plot 'hbr-gj.txt'  using 1:2 with  points  pointswidth 2 title '  Cocrystal'
#plot 'energy.txt'      using 1:2 with lp lw 2 ps 4 title 'ReaxFF'

plot  'distribution.txt'      using 1:2 with  points  ps 2   title 'CH3NO2',\
        'distribution.txt'       using 1:3 with  points  ps 2   title '           H',\
        'distribution.txt'       using 1:4 with  points  ps 2   title '         HO',\
        'distribution.txt'       using 1:5 with  points  ps 2   title '      H2O',\
        'distribution.txt'       using 1:7 with  points  ps 2   title '       NO',\
        'distribution.txt'       using 1:8 with  points  ps 2   title '      NO2',\
        'distribution.txt'      using 1:11 with  points  ps 2  title '        CO',\
        'distribution.txt'      using 1:12 with  points  ps 2  title '       CO2',\
        'distribution.txt'      using 1:13 with  points  ps 2  title '     CH2O',\
        'distribution.txt'      using 1:17 with  points  ps 2  title '        CH3',\
        'distribution.txt'      using 1:18 with  points  ps 2  title '        CH4'






!Set the output to a png file
set terminal png size 500,500

!The file we'll write to
set output 'marea_graf.png'

!The graphic title
!set title 

!X axis
set xlabel "Tiempo (s)"

!Y axis
set ylabel "Altura de la marea (m)"

!plot the graphic
plot "datos.dat"

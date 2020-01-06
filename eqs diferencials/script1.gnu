# Format i nom de la imatge
#set term png enhanced
#set output "P7-1920-fig1.png"

# Permet escriure lletres gregues i altres mogudes
set encoding utf8

# Mostra els eixos
#set xzeroaxis
#set yzeroaxis

# Títol del gràfic
set title "{/Symbol q}(t)"

# Rang dels eixos
#set xrange[-1.00e-10:1.00e10]
#set yrange[-1.00e-10:1.00e10]

# Títols dels eixos
set xlabel "Temps, t(s)"
set ylabel "Angle, {/Symbol q}(rad)"

# Canvia els nombres dels eixos per nombres personalitzats
#set ytics("1x10^-^1^0" 1.00e-10,"1x10^-^0^5" 1.00e-05,"1x10^0" 1.00e+00,"1x10^5" 1.00e+05,"1x10^1^0" 1.00e+10)
#set xtics("1x10^-^3" 1.00e-03,"1x10^-^2" 1.00e-02,"1x10^-^1" 1.00e-01,"1x10^0" 1.00e+00,"1x10^1" 1.00e+01)

# Format dels nombres dels eixos
set format y '%.2f'
set format x '%.2f'

# Escala dels eixos logarítmica
#set logscale y
#set logscale x

# Posició de la llegenda
set key top left

# Plot 
plot "aux.dat" using 1:2 with points t "Euler"
pause -1

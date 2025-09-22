#Plots dsigdptdy

set terminal postscript eps color

set xlabel font 'Times-Roman,22'
set ylabel font 'Times-Roman,22'
set label font 'Helvetica,18'
set key box
set key font 'Helvetica,18'
set key sample 2 
set key spacing 1.5
set tics scale 1.5
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"

dist_Y   = '../Output/dsig_dy.dat'
dist_pT  = '../Output/dsig_dpt.dat'
dist_m   = '../Output/dsig_dm.dat'
#dist_Y_T   = '../Output/dsig_dy_T.dat'
#dist_pT_T  = '../Output/dsig_dpt_T.dat'
#dist_Y_L   = '../Output/dsig_dy_L.dat'
#dist_pT_L  = '../Output/dsig_dpt_L.dat'
ypTdist_2p25 = '../Output/dsig_dydpt_y2p25.dat'
ypTdist_2p75 = '../Output/dsig_dydpt_y2p75.dat'
ypTdist_3p25 = '../Output/dsig_dydpt_y3p25.dat'
ypTdist_3p75 = '../Output/dsig_dydpt_y3p75.dat'
ypTdist_4p25 = '../Output/dsig_dydpt_y4p25.dat'

dist_Y_grid  = '../Output/y_dist_grid.dat'

DATA_Ydist  = '../DATA/Ydist_dataset.csv'
DATA_pTdist = '../DATA/pTdist_dataset.csv'
DATA_ypTdist_2p25 = '../DATA/YpTdist_2p25_dataset.csv'
DATA_ypTdist_2p75 = '../DATA/YpTdist_2p75_dataset.csv'
DATA_ypTdist_3p25 = '../DATA/YpTdist_3p25_dataset.csv'
DATA_ypTdist_3p75 = '../DATA/YpTdist_3p75_dataset.csv'
DATA_ypTdist_4p25 = '../DATA/YpTdist_4p25_dataset.csv'

set style increment default
set samples 50, 50

#set title 'LHCb Collaboration - R. Aaij et al., JHEP 07 (2022), 026' 


set logscale x
set logscale y

set output 'pTdist.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dp_T [pb/GeV]' 
set xrange[0:200]
plot dist_pT u 1:2 w l lc rgb 'red' lw 3 t 'IP-SAT' , \
     DATA_pTdist     using 1:4:2:5 with xyerrorbars pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV'
#     dist_pT_T u 1:2 w l lc rgb 'blue' lw 3  t 'T' , \
#     dist_pT_L u 1:2 w l lc rgb 'green' lw 3 t 'L' , \
unset labeldist_Y_alt = '../Output/y_dist_grid.dat'
unset logscale y 
unset xrange

set output 'YpTdist_2p25.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '2 < Y < 2.5' at graph 0.85, graph 0.75

plot ypTdist_2p25 u 1:2 w l lc rgb 'blue' lw 3  t 'IP-SAT', \
     DATA_ypTdist_2p25     using 1:4:2:5 with xyerrorbars pt 7 lc rgb 'black' t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label

set output 'YpTdist_2p75.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '2.5 < Y < 3' at graph 0.85, graph 0.75

plot ypTdist_2p75 u 1:2 w l lc rgb 'blue' lw 3  t 'IP-SAT', \
     DATA_ypTdist_2p75     using 1:4:2:5 with xyerrorbars  pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV' 
unset label



set output 'YpTdist_3p25.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '3 < Y < 3.5' at graph 0.85, graph 0.75

plot ypTdist_3p25 u 1:2 w l lc rgb 'blue' lw 3  t 'IP-SAT', \
     DATA_ypTdist_3p25     using 1:4:2:5 with xyerrorbars  pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label



set output 'YpTdist_3p75.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '3.5 < Y < 4' at graph 0.85, graph 0.75

plot ypTdist_3p75 u 1:2 w l lc rgb 'blue' lw 3  t 'IP-SAT', \
     DATA_ypTdist_3p75    using 1:4:2:5 with xyerrorbars pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label



set output 'YpTdist_4p25.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '4 < Y < 4.5' at graph 0.85, graph 0.75

plot  ypTdist_4p25 u 1:2 w l lc rgb 'blue' lw 3  t 'IP-SAT', \
     DATA_ypTdist_4p25     using 1:4:2:5 with xyerrorbars pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label



     

unset logscale x 
set output 'Ydist.eps'
set xlabel 'Y' 
set ylabel 'd{/Symbol s}/dY [pb]' 
set xrange[2:4.5]
set yrange[0:350]
#set label 'Only transverse' at graph 0.75, graph 0.7
#set label 'PT= 5' at graph 0.8, graph 0.8
plot dist_Y u 1:2 w l lc rgb 'red' lw 3 t 'IP-SAT' , \
     DATA_Ydist u 1:4:2:5 w xyerrorbars pt 7 lc rgb 'black'   t 'LHCb data {/Symbol=\326}s = 13 TeV'#,  \
     #dist_Y_grid u 1:2 w l dt 2 lc rgb 'blue' lw 3 t 'grid'
      #dist_Y_T u 1:2 w l lc rgb 'blue' lw 3  t 'T' , \
     #dist_Y_L u 1:2 w l lc rgb 'green' lw 3  t 'L' , \
unset label



     


     
     





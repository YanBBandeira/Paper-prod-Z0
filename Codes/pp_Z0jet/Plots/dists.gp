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

gbw_dist_Y   = '../Output/gbw_dsig_dy.dat'
gbw_dist_pT  = '../Output/gbw_dsig_dpt.dat'
gbw_dist_m   = '../Output/gbw_dsig_dm.dat'
gbw_ypTdist_2p25 = '../Output/gbw_dsig_dydpt_y2p25.dat'
gbw_ypTdist_2p75 = '../Output/gbw_dsig_dydpt_y2p75.dat'
gbw_ypTdist_3p25 = '../Output/gbw_dsig_dydpt_y3p25.dat'
gbw_ypTdist_3p75 = '../Output/gbw_dsig_dydpt_y3p75.dat'
gbw_ypTdist_4p25 = '../Output/gbw_dsig_dydpt_y4p25.dat'


kslinear_dist_Y   = '../Output/kslinear_dsig_dy.dat'
kslinear_dist_pT  = '../Output/kslinear_dsig_dpt.dat'
kslinear_dist_m   = '../Output/kslinear_dsig_dm.dat'
kslinear_ypTdist_2p25 = '../Output/kslinear_dsig_dydpt_y2p25.dat'
kslinear_ypTdist_2p75 = '../Output/kslinear_dsig_dydpt_y2p75.dat'
kslinear_ypTdist_3p25 = '../Output/kslinear_dsig_dydpt_y3p25.dat'
kslinear_ypTdist_3p75 = '../Output/kslinear_dsig_dydpt_y3p75.dat'
kslinear_ypTdist_4p25 = '../Output/kslinear_dsig_dydpt_y4p25.dat'

ksnonlinear_dist_Y   = '../Output/ksnonlinear_dsig_dy.dat'
ksnonlinear_dist_pT  = '../Output/ksnonlinear_dsig_dpt.dat'
ksnonlinear_dist_m   = '../Output/ksnonlinear_dsig_dm.dat'
ksnonlinear_ypTdist_2p25 = '../Output/ksnonlinear_dsig_dydpt_y2p25.dat'
ksnonlinear_ypTdist_2p75 = '../Output/ksnonlinear_dsig_dydpt_y2p75.dat'
ksnonlinear_ypTdist_3p25 = '../Output/ksnonlinear_dsig_dydpt_y3p25.dat'
ksnonlinear_ypTdist_3p75 = '../Output/ksnonlinear_dsig_dydpt_y3p75.dat'
ksnonlinear_ypTdist_4p25 = '../Output/ksnonlinear_dsig_dydpt_y4p25.dat'

ccfmset_dist_Y   = '../Output/ccfmset_dsig_dy.dat'
ccfmset_dist_pT  = '../Output/ccfmset_dsig_dpt.dat'
ccfmset_dist_m   = '../Output/ccfmset_dsig_dm.dat'
ccfmset_ypTdist_2p25 = '../Output/ccfmset_dsig_dydpt_y2p25.dat'
ccfmset_ypTdist_2p75 = '../Output/ccfmset_dsig_dydpt_y2p75.dat'
ccfmset_ypTdist_3p25 = '../Output/ccfmset_dsig_dydpt_y3p25.dat'
ccfmset_ypTdist_3p75 = '../Output/ccfmset_dsig_dydpt_y3p75.dat'
ccfmset_ypTdist_4p25 = '../Output/ccfmset_dsig_dydpt_y4p25.dat'



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

plot gbw_dist_pT u 1:2 w l lc rgb 'red' lw 3 t 'GBW' , \
     kslinear_dist_pT u 1:2 w l lc rgb 'blue' lw 3  t 'kslinear' , \
     ksnonlinear_dist_pT u 1:2 w l lc rgb 'green' lw 3 t 'ksnonlinear' , \
     ccfmset_dist_pT u 1:2 w l lc rgb 'purple' lw 3 t 'ccfmset' , \
     DATA_pTdist using 1:4:($1+$2):($1+$3):($4-$5):($4+$6) with xyerrorbars pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV' 
#unset labeldist_Y_alt = '../Output/y_dist_grid.dat'
#unset logscale y 
unset xrange

set output 'YpTdist_2p25.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '2 < Y < 2.5' at graph 0.2, graph 0.25

plot gbw_ypTdist_2p25 u 1:2 w l lc rgb 'red' lw 3  t 'GBW', \
     kslinear_ypTdist_2p25 u 1:2 w l lc rgb 'red' lw 3  t 'kslinear' , \
     ksnonlinear_ypTdist_2p25 u 1:2 w l lc rgb 'green' lw 3 t 'ksnonlinear' , \
     ccfmset_ypTdist_2p25 u 1:2 w l lc rgb 'purple' lw 3 t 'ccfmset' , \
     DATA_ypTdist_2p25     using 1:4:2:5 with xyerrorbars pt 7 lc rgb 'black' t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label

set output 'YpTdist_2p75.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '2.5 < Y < 3' at graph 0.2, graph 0.25

plot gbw_ypTdist_2p75 u 1:2 w l lc rgb 'red' lw 3  t 'GBW', \
     kslinear_ypTdist_2p75 u 1:2 w l lc rgb 'blue' lw 3  t 'kslinear', \
     ksnonlinear_ypTdist_2p75 u 1:2 w l lc rgb 'green' lw 3  t 'ksnonlinear', \
     ccfmset_ypTdist_2p75 u 1:2 w l lc rgb 'purple' lw 3  t 'ccfmset', \
     DATA_ypTdist_2p75     using 1:4:2:5 with xyerrorbars  pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV' 
unset label



set output 'YpTdist_3p25.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '3 < Y < 3.5' at graph 0.2, graph 0.25

plot gbw_ypTdist_3p25 u 1:2 w l lc rgb 'red' lw 3  t 'GBW', \
     kslinear_ypTdist_3p25 u 1:2 w l lc rgb 'blue' lw 3  t 'kslinear', \
     ksnonlinear_ypTdist_3p25 u 1:2 w l lc rgb 'green' lw 3  t 'ksnonlinear', \
     ccfmset_ypTdist_3p25 u 1:2 w l lc rgb 'purple' lw 3  t 'ccfmset', \
     DATA_ypTdist_3p25     using 1:4:2:5 with xyerrorbars  pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label



set output 'YpTdist_3p75.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '3.5 < Y < 4' at graph 0.2, graph 0.25

plot gbw_ypTdist_3p75 u 1:2 w l lc rgb 'red' lw 3  t 'GBW', \
     kslinear_ypTdist_3p75 u 1:2 w l lc rgb 'blue' lw 3  t 'kslinear', \
     ksnonlinear_ypTdist_3p75 u 1:2 w l lc rgb 'green' lw 3  t 'ksnonlinear', \
     ccfmset_ypTdist_3p75 u 1:2 w l lc rgb 'purple' lw 3  t 'ccfmset', \
     DATA_ypTdist_3p75    using 1:4:2:5 with xyerrorbars pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label



set output 'YpTdist_4p25.eps'
set xlabel 'p_T [GeV]' 
set ylabel 'd{/Symbol s}/dYdp_T [pb/GeV]' 
set label '4 < Y < 4.5' at graph 0.2, graph 0.25

plot gbw_ypTdist_4p25 u 1:2 w l lc rgb 'red' lw 3  t 'GBW', \
     kslinear_ypTdist_4p25 u 1:2 w l lc rgb 'blue' lw 3  t 'kslinear', \
     ksnonlinear_ypTdist_4p25 u 1:2 w l lc rgb 'green' lw 3  t 'ksnonlinear', \
     ccfmset_ypTdist_4p25 u 1:2 w l lc rgb 'purple' lw 3  t 'ccfmset', \
     DATA_ypTdist_4p25     using 1:4:2:5 with xyerrorbars pt 7 lc rgb 'black'  t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label


unset logscale x 
unset logscale y
set output 'Ydist.eps'
set xlabel 'Y' 
set ylabel 'd{/Symbol s}/dY [pb]' 
set xrange[2:4.5]
set yrange[0:350]
#set label 'Only transverse' at graph 0.75, graph 0.7
#set label 'PT= 5' at graph 0.8, graph 0.8
plot gbw_dist_Y u 1:2 w l lc rgb 'red' lw 3 t 'GBW' , \
     kslinear_dist_Y u 1:2 w l lc rgb 'blue' lw 3  t 'kslinear' , \
     ksnonlinear_dist_Y u 1:2 w l lc rgb 'green' lw 3 t 'ksnonlinear' , \
     ccfmset_dist_Y u 1:2 w l lc rgb 'purple' lw 3 t 'ccfmset' , \
     DATA_Ydist u 1:4:2:5 w xyerrorbars pt 7 lc rgb 'black'   t 'LHCb data {/Symbol=\326}s = 13 TeV'
unset label

     





     


     
     





#1. Plot g(r): Al-Cl
set terminal png
set termopt enhanced
set output 'figures/gr_Al_Cl.png'
set title 'g(r): Al-Cl'
set xlabel 'r (nm)' font ",15"
set ylabel 'g(r)' font ",15"
unset key
p [0:4] "rdf_ALCL3_MRAT_13-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "rdf_ALCL3_MRAT_14-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1:4:1",\
  "rdf_ALCL3_MRAT_17-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.7:1"

#2. Plot g(r): Al-Al
set terminal png
set termopt enhanced
set output 'figures/gr_Al_Al.png'
set title 'g(r): Al-Al'
set xlabel 'r (nm)' font ",15"
set ylabel 'g(r)' font ",15"
unset key
p [0:4.5] "rdf_ALCL3_MRAT_13-pos-1.xyz" u 1:3 w l lw 3 lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "rdf_ALCL3_MRAT_14-pos-1.xyz" u 1:3 w l lw 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1:4:1",\
  "rdf_ALCL3_MRAT_17-pos-1.xyz" u 1:3 w l lw 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.7:1"


# Plot Cl-neighbors
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/nCl_neigh.png'
set auto x
set title 'Cl Neighbors around Al'
set xlabel 'n_{Cl}' font ",15"
set ylabel 'f(n_{Cl})' font ",15"
unset key
p [0:4.8] "catanneigh_ALCL3_MRAT_13-pos-1.xyz" u 3:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "catanneigh_ALCL3_MRAT_14-pos-1.xyz" u 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "catanneigh_ALCL3_MRAT_14-pos-1.xyz" u 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.7:1"

# Plot Al-neighbors
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/nAl_neigh.png'
set auto x
set title 'Al Neighbors around Cl'
set xlabel 'n_{Al}' font ",15"
set ylabel 'f(n_{Al})' font ",15"
unset key
p [0:4.8] "catanneigh_ALCL3_MRAT_13-pos-1.xyz" u 5:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "catanneigh_ALCL3_MRAT_14-pos-1.xyz" u 5 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "catanneigh_ALCL3_MRAT_14-pos-1.xyz" u 5 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.7:1"

# Plot Al-Cl clusters
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/nclust.png'
set auto x
set title 'Al-Cl Clusters'
set xlabel 's_i' font ",15"
set ylabel 'f(s_{i})' font ",15"
unset key
p [0.01:12] "clust_ALCL3_MRAT_13-pos-1.xyz" u 3:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "clust_ALCL3_MRAT_14-pos-1.xyz" u 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "clust_ALCL3_MRAT_14-pos-1.xyz" u 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.7:1"

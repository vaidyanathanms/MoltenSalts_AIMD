#1. Plot g(r): Al-Al
set terminal png
set termopt enhanced
set output 'figures/AcONH2/gr_Al_Al.png'
set title 'g(r): Al-Al'
set xlabel 'r (nm)' font ",15"
set ylabel 'g(r)' font ",15"
unset key
p [0:4] "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/rdf_traj_nvtprod_14-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1:4:1",\
  "allresults_AcONH2/results_15/rdf_traj_nvtprod_15-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/rdf_traj_nvtprod_17-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"

#2. Plot g(r): Al-Cl
set terminal png
set termopt enhanced
set output 'figures/AcONH2/gr_Al_Cl.png'
set title 'g(r): Al-Cl'
set xlabel 'r (Angstroms)' font ",15"
set ylabel 'g(r)' font ",15"
unset key
p [0:4] "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:3 w l lw 3 lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/rdf_traj_nvtprod_14-pos-1.xyz" u 1:3 w l lw 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1:4:1",\
  "allresults_AcONH2/results_15/rdf_traj_nvtprod_15-pos-1.xyz" u 1:3 w l lw 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/rdf_traj_nvtprod_17-pos-1.xyz" u 1:3 w l lw 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"


#3. Plot g(r): Al-N
set terminal png
set termopt enhanced
set output 'figures/AcONH2/gr_Al_N.png'
set title 'g(r): Al-N'
set xlabel 'r (Angstroms)' font ",15"
set ylabel 'g(r)' font ",15"
unset key
p [0:4] "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:4 w l lw 3 lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/rdf_traj_nvtprod_14-pos-1.xyz" u 1:4 w l lw 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1:4:1",\
  "allresults_AcONH2/results_15/rdf_traj_nvtprod_15-pos-1.xyz" u 1:4 w l lw 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/rdf_traj_nvtprod_17-pos-1.xyz" u 1:4 w l lw 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"

#4. Plot g(r): Al-O
set terminal png
set termopt enhanced
set output 'figures/AcONH2/gr_Al_O.png'
set title 'g(r): Al-O'
set xlabel 'r (Angstroms)' font ",15"
set ylabel 'g(r)' font ",15"
unset key
p [0:4] "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:5 w l lw 3 lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/rdf_traj_nvtprod_14-pos-1.xyz" u 1:5 w l lw 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1:4:1",\
  "allresults_AcONH2/results_15/rdf_traj_nvtprod_15-pos-1.xyz" u 1:5 w l lw 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/rdf_traj_nvtprod_17-pos-1.xyz" u 1:2 w l lw 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"

#5. Plot g(r): Cl-C
set terminal png
set termopt enhanced
set output 'figures/AcONH2/gr_Cl_C.png'
set title 'g(r): Cl-C'
set xlabel 'r (Angstroms)' font ",15"
set ylabel 'g(r)' font ",15"
unset key
p [0:4] "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:6 w l lw 3 lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/rdf_traj_nvtprod_14-pos-1.xyz" u 1:6 w l lw 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1:4:1",\
  "allresults_AcONH2/results_15/rdf_traj_nvtprod_15-pos-1.xyz" u 1:6 w l lw 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/rdf_traj_nvtprod_17-pos-1.xyz" u 1:6 w l lw 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"


#6. Plot g(r): all_for_one_case
set terminal png
set termopt enhanced
set output 'figures/AcONH2/gr_Al_all_13.png'
set title 'g(r): All Al neighbors'
set key
set xlabel 'r (Angstroms)' font ",15"
set ylabel 'g(r)' font ",15"
p [0:4] "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:2 w l lw 2 lc rgb 'red'  title "Al - Al",\
  "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:3 w l lw 2 lc rgb 'green'  title "Al - Cl",\
  "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:4 w l lw 2 lc rgb 'blue'   title "Al - N",\
  "allresults_AcONH2/results_13/rdf_traj_nvtprod-pos-1.xyz" u 1:5 w l lw 2 lc rgb 'black'  title "Al - O"

#7. Plot Cl neighbors around Al
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/AcONH2/nCl_Al_neigh.png'
set auto x
set title 'Cl Neighbors around Al'
set xlabel 'n_{Cl}' font ",15"
set ylabel 'f(n_{Cl})' font ",15"
unset key
p [-0.8:5] "allresults_AcONH2/results_13/catanneigh_12_traj_nvtprod-pos-1.xyz" u 3:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/catanneigh_12_traj_nvtprod_14-pos-1.xyz" u 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "allresults_AcONH2/results_15/catanneigh_12_traj_nvtprod_15-pos-1.xyz" u 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/catanneigh_12_traj_nvtprod_17-pos-1.xyz" u 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"

#8. Plot N neighbors around Al
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/AcONH2/nN_Al_neigh.png'
set auto x
set title 'N Neighbors around Al'
set xlabel 'n_{N}' font ",15"
set ylabel 'f(n_{N})' font ",15"
unset key
p [-0.8:5] "allresults_AcONH2/results_13/catanneigh_14_traj_nvtprod-pos-1.xyz" u 3:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/catanneigh_14_traj_nvtprod_14-pos-1.xyz" u 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "allresults_AcONH2/results_15/catanneigh_14_traj_nvtprod_15-pos-1.xyz" u 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/catanneigh_14_traj_nvtprod_17-pos-1.xyz" u 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"

#9. Plot O neighbors around Al
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/AcONH2/nO_Al_neigh.png'
set auto x
set title 'O Neighbors around Al'
set xlabel 'n_{O}' font ",15"
set ylabel 'f(n_{O})' font ",15"
unset key
p [-0.8:5] "allresults_AcONH2/results_13/catanneigh_16_traj_nvtprod-pos-1.xyz" u 3:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/catanneigh_16_traj_nvtprod_14-pos-1.xyz" u 3 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "allresults_AcONH2/results_15/catanneigh_16_traj_nvtprod_15-pos-1.xyz" u 3 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/catanneigh_16_traj_nvtprod_17-pos-1.xyz" u 3 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"


#10. Plot Al-Cl clusters
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/AcONH2/AlCl_clust.png'
set auto x
set title 'Al-Cl Clusters'
set xlabel 's_i' font ",15"
set ylabel 'f(s_{i})' font ",15"
unset key
p [-0.8:15] "allresults_AcONH2/results_13/clust_12_traj_nvtprod-pos-1.xyz" u 2:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/clust_12_traj_nvtprod_14-pos-1.xyz" u 2 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "allresults_AcONH2/results_15/clust_12_traj_nvtprod_15-pos-1.xyz" u 2 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/clust_12_traj_nvtprod_17-pos-1.xyz" u 2 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"


#11. Plot Al-N clusters
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/AcONH2/AlN_clust.png'
set auto x
set title 'Al-N Clusters'
set xlabel 's_i' font ",15"
set ylabel 'f(s_{i})' font ",15"
unset key
p [-0.8:5] "allresults_AcONH2/results_13/clust_14_traj_nvtprod-pos-1.xyz" u 2:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/clust_14_traj_nvtprod_14-pos-1.xyz" u 2 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "allresults_AcONH2/results_15/clust_14_traj_nvtprod_15-pos-1.xyz" u 2 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/clust_14_traj_nvtprod_17-pos-1.xyz" u 2 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"

# Plot Al-O clusters
set terminal png
set termopt enhanced
set style data histogram
set style histogram clustered gap 1
set style fill solid
set boxwidth 0.9
set output 'figures/AcONH2/AlO_clust.png'
set auto x
set title 'Al-O Clusters'
set xlabel 's_i' font ",15"
set ylabel 'f(s_{i})' font ",15"
unset key
p [-1:8] "allresults_AcONH2/results_13/clust_16_traj_nvtprod-pos-1.xyz" u 2:xtic(1) lc rgb 'red'   title "AlCl_3:AcONH_2 = 1.3:1",\
  "allresults_AcONH2/results_14/clust_16_traj_nvtprod_14-pos-1.xyz" u 2 lc rgb 'green' title "AlCl_3:AcONH_2 = 1.4:1",\
  "allresults_AcONH2/results_15/clust_16_traj_nvtprod_15-pos-1.xyz" u 2 lc rgb 'blue'  title "AlCl_3:AcONH_2 = 1.5:1",\
  "allresults_AcONH2/results_17/clust_16_traj_nvtprod_17-pos-1.xyz" u 2 lc rgb 'black' title "AlCl_3:AcONH_2 = 1.7:1"

import os, sys
from runnning_plots import variables

input_root_file = '../GF_HH_Benchmark222.root'

for key in variables:
  #print variables[key]
  #python make_aQGC_plots.py -n ../GF_HH_Benchmark222.root --grid -v gen_leading_photon_Pt gen_Subleading_photon_Pt --leg "Leading" "Sub Leading" --xmin 0 --xmax 500 --xlabel "Photon p_{T}" -o photon_pt.png --legPos "tr" --nbin 51
  print 'python make_aQGC_plots.py -n '+input_root_file + ' --grid  -v "' + key + '" --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + variables[key][1] + '" -o ' + key + '.png  --legPos "tr"  --nbin ' + str(variables[key][3])
  toPlot = 'python make_aQGC_plots.py -n '+input_root_file + ' --grid  -v "' + key + '" --leg "' + variables[key][1] +'" --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + variables[key][1] + '" -o ' + key + '.png  --legPos "tr"  --nbin ' + str(variables[key][3])
  os.system(toPlot)

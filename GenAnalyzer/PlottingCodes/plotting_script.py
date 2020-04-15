import os, sys
from runnning_plots import variables
import runnning_plots2

input_root_file = '../GF_HH_Benchmark222.root'

singleVariable = False
multipleVariable = True

if (singleVariable):
  for key in variables:
    #print variables[key]
    #python make_aQGC_plots.py -n ../GF_HH_Benchmark222.root --grid -v gen_leading_photon_Pt gen_Subleading_photon_Pt --leg "Leading" "Sub Leading" --xmin 0 --xmax 500 --xlabel "Photon p_{T}" -o photon_pt.png --legPos "tr" --nbin 51
      print 'python make_aQGC_plots.py -n '+input_root_file + ' --grid  -v "' + key + '" --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + variables[key][1] + '" -o ' + key + '.png  --legPos "tr"  --nbin ' + str(variables[key][3])
      toPlot = 'python make_aQGC_plots.py -n '+input_root_file + ' --grid  -v "' + key + '" --leg "' + variables[key][1] +'" --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + variables[key][1] + '" -o ' + key + '.png  --legPos "tr"  --nbin ' + str(variables[key][3])
      os.system(toPlot)

if (multipleVariable):
  for var_arr in runnning_plots2.variables:
    varToAdd = ""
    legendToAdd = ""
    key = []
    for local_var in var_arr:
      key = local_var
      varToAdd += '"' + local_var + '" '
      legendToAdd += '"' + variables[key][1] + '" '
    #print varToAdd
    toPlot = 'python make_aQGC_plots.py -n '+input_root_file + ' --grid  -v ' + varToAdd + ' --leg ' + legendToAdd +' --xmin ' + str(variables[key][4]) + ' --xmax ' + str(variables[key][5]) + ' --xlabel "' + (variables[key][1]).replace('Leading','').replace('SubLeading','').replace('Sub','') + '" -o ' + key + '_multi.png  --legPos "tr"  --nbin ' + str(variables[key][3])
    print "\n\n","*"*51
    print toPlot
    os.system(toPlot)

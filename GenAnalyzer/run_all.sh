inputFiles=( GF_HH_Benchmark1 GF_HH_Benchmark2 GF_HH_Benchmark3 GF_HH_Benchmark4 GF_HH_Benchmark5 GF_HH_Benchmark6 GF_HH_Benchmark7 GF_HH_Benchmark8 GF_HH_Benchmark9 GF_HH_Benchmark10 GF_HH_Benchmark11 GF_HH_Benchmark12 GF_HH_BenchmarkSM )
#inputFiles=( GF_HH_Benchmark1 )

nFiles=${#inputFiles[@]}

for ((ifile=0; ifile<$nFiles; ifile++))
  do
  echo "==================================================================="
  sed -i "s/inputTxtFile=.*/inputTxtFile=\"${inputFiles[${ifile}]}\"/g" python/ConfFile_cfg.py
  cat python/ConfFile_cfg.py
  echo "==================================================================="
  cmsRun python/ConfFile_cfg.py
done

# inputFiles=( Radion_hh_narrow_M1000_LHEBqrk Radion_hh_narrow_M1500_LHEBqrk Radion_hh_narrow_M260_LHEBqrk Radion_hh_narrow_M270_LHEBqrk Radion_hh_narrow_M500_LHEBqrk )
inputFiles=( GF_HH_10_slc6_LHEBqrk GF_HH_11_slc6_LHEBqrk GF_HH_12_slc6_LHEBqrk GF_HH_1_slc6_LHEBqrk GF_HH_2_slc6_LHEBqrk GF_HH_3_slc6_LHEBqrk GF_HH_4_slc6_LHEBqrk GF_HH_5_slc6_LHEBqrk GF_HH_6_slc6_LHEBqrk GF_HH_7_slc6_LHEBqrk GF_HH_8_slc6_LHEBqrk GF_HH_9_slc6_LHEBqrk GF_HH_SM_slc6_LHEBqrk )
nFiles=${#inputFiles[@]}

for ((ifile=0; ifile<$nFiles; ifile++))
  do
  echo "==================================================================="
  sed -i "s/inputTxtFile=.*/inputTxtFile=\"${inputFiles[${ifile}]}\"/g" python/ConfFile_cfg.py
  cat python/ConfFile_cfg.py
  echo "==================================================================="
  cmsRun python/ConfFile_cfg.py | tee ${inputFiles[${ifile}]}.log
done

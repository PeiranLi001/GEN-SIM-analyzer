# this script file is useful when we need to run configuration file for more than 255 files
# Just we have to replace the output of this file with the lines in configuration file where we need to put outputs
# i.e. we need to replace following part 
# #####	process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
# ##### 'file:/eos/uscms/store/user/rasharma/merged_patTuple_500_8T/WW_TuneZ2star.root'
# #####  ))
# in configuration with the output of this file

# here you need to enter the path of files to list
# remember that the path should be end with / otherwise script will show you an error

#cd /pnfs/cms/WAX/11	# I put this line because I don't need this part in my list of files
echo "Enter the path of input files to list:"
set pth=$<
set count=1
echo "myfilelist = cms.untracked.vstring()"
echo "myfilelist.extend( ["
foreach f (${pth}*)
	@ count=${count} + 1
	@ y=${count} % 250
	if ( $y == 0 )	then
		echo "\t'file:$f'" 
		echo "])"
		echo "myfilelist.extend( ["
	else
		echo "\t'file:$f'," 
	endif
end
echo "])"
echo 'process.source = cms.Source("PoolSource", fileNames = myfilelist)'
#cd -

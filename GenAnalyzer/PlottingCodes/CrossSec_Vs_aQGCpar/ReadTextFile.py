f = open("Test.txt","r")
#print f.readlines()
aQGC_Val_Initializer = "FT0_12p5"
FileCategory = ["fs0", "fs1", "ft0", "ft1", "ft2", "fm0", "fm1", "fm6", "fm7"]
FileContent = f.readlines()
fs0 = [row for row in FileContent if 'fs0' in row]
fs1 = [row for row in FileContent if 'fs1' in row]
ft0 = [row for row in FileContent if 'ft0' in row]
ft1 = [row for row in FileContent if 'ft1' in row]
ft2 = [row for row in FileContent if 'ft2' in row]
fm0 = [row for row in FileContent if 'fm0' in row]
fm1 = [row for row in FileContent if 'fm1' in row]
fm6 = [row for row in FileContent if 'fm6' in row]
fm7 = [row for row in FileContent if 'fm7' in row]

f.close()



count = 0
for f in [fs0, fs1, ft0, ft1, ft2, fm0, fm1, fm6, fm7]:
	text_file = open(aQGC_Val_Initializer+"_"+FileCategory[count]+".txt", "w")
	for i in f:
		temp =  i.split("\t")
		temp[0] = temp[0].replace(FileCategory[count]+"_","")
		temp[0] = temp[0].replace("m","-")
		temp[0] = temp[0].replace("p",".")
		temp[2] = temp[2].replace("\n","")
		#print temp[0],"\t",temp[1],"\t",temp[2]
		text_file.write("%s \t %s \t %s \n" % (temp[0], temp[1], temp[2]))
	text_file.close()
	count += 1

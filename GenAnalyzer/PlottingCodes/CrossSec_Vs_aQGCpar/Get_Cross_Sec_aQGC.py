#!/usr/bin/python

# First argument = Input path 
import sys
import os
from glob import glob
from ROOT import *


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 2:
	print "please enter path of input directory as one of argumetns of command"
	exit(0)

gROOT.SetBatch(True)

inputFolder = sys.argv[1]

inputfolderSplit = sys.argv[1].split("/")
print inputfolderSplit
scaleFactor = 1.0


#files = glob(inputFolder+"/*.root")
#files = glob(inputFolder+"/*NPle1_mjj200*.root")
files = glob(inputFolder+"/*NPle1_RunCardChanged_InitilaizeAllaQGCPar*.root")
chain = TChain("otree")
for f in files:
	print "Reading file: ",f
	chain.Add(f)
entries = chain.GetEntries()
print "Total number of events: ",entries


if inputfolderSplit[-1] == "":
	TreeName = inputfolderSplit[-2]
else:	
	TreeName = inputfolderSplit[-1]

A = [0]*718
SumA = [0]*718
AvgA = [0]*718
DevA = [0]*718
Verbose = 0

nEvents = 0
for i in xrange(entries):
	if i%10000==0:
		print "nEntries = ",i
	#if i>3:
	#	break
	chain.GetEntry(i)

	count = 0
	for j in xrange(len(chain.LHEWeightIDs)):
		if j<446:
			continue
		else:
			SumA[count] += chain.LHEWeights[j]/chain.LHEWeights[0]
			#if count == 0:
			#	print "\tSumA = ", chain.LHEWeights[j]/chain.LHEWeights[0]
			#	print "\tSumA[count ] = ",SumA[count]
			count += 1
	nEvents += 1

# Calculate St. Dev.

for i in xrange(entries):
	if i%10000==0:
		print "nEntries = ",i
	#if i>3:
	#	break
	chain.GetEntry(i)
	count = 0
	for j in xrange(len(chain.LHEWeightIDs)):
		if j<446:
			continue
		else:
			AvgA = SumA[count]/entries
			#if count < 3:
			#	print SumA[count]/entries
			DevA[count] += ( chain.LHEWeights[j]/chain.LHEWeights[0] -AvgA)*(chain.LHEWeights[j]/chain.LHEWeights[0] -AvgA)
			#if count < 3:
			#	print  ( chain.LHEWeights[j]/chain.LHEWeights[0] -AvgA)*(chain.LHEWeights[j]/chain.LHEWeights[0] -AvgA)
			count += 1


for j in xrange(len(DevA)):
	print chain.LHEWeightIDs[446+j],"\t",SumA[j]*17.11," +/- ",DevA[j]/len(DevA)

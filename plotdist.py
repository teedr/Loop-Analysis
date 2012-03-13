#!/usr/bin/python

import sys
import re
import operator
import numpy as np
import json
import start_api as st
import numpy as np
import matplotlib.pyplot as plt

def sort_lengths(SequenceList):
	sorted_sequence_dict = {}
	for sequence in SequenceList:
		length = len(sequence)
		if length not in sorted_sequence_dict:
			sorted_sequence_dict[length] = {}
			sorted_sequence_dict[length]['sequences'] = []
		sorted_sequence_dict[length]['sequences'].append(sequence)
	return sorted_sequence_dict

def initPositionInfo(MasterDict):
	for key in MasterDict:
		MasterDict[key]['positions'] = {}
		for position in range(0, key):
			MasterDict[key]['positions'][position] = {}
	return MasterDict

def generatePositionInfo(MasterDict, SequenceList):
	for sequence in SequenceList:
		length = len(sequence)
		for idx, residue in enumerate(sequence):
			if residue not in MasterDict[length]['positions'][idx]:
				MasterDict[length]['positions'][idx][residue] = 1
			else:
				MasterDict[length]['positions'][idx][residue] += 1
	return MasterDict
def elementwise_addition(arr1, arr2):
	if len(arr1) != len(arr2):
		print "error in elementwise_addition: arrays must be same length"
	sumArray = [0]*len(arr1)
	for idx, element in enumerate(arr1):
		sumArray[idx] = arr1[idx] + arr2[idx]
	return sumArray

def color(key):
	if key == 'A':
		color = '#A50B5E'
	elif key == 'C':
		color = '#C72C48'
	elif key == 'E':
		color = '#C19A6B'
	elif key == 'D':
		color = '#00009C'
	elif key == 'G':
		color = '#F8DE7E'
	elif key == 'F':
		color = '#E1A95F'
	elif key == 'I':
		color = '#85BB65'
	elif key == 'H':
		color = '#FF9933'
	elif key == 'K':
		color = '#BDDA57'
	elif key == 'M':
		color = '#004B49'
	elif key == 'L':
		color = '#704241'
	elif key == 'N':
		color = '#B94E48'
	elif key == 'Q':
		color = '#6C541E'
	elif key == 'P':
		color = '#86608E'
	elif key == 'S':
		color = '#D73B3E'
	elif key == 'R':
		color = '#29AB87'
	elif key == 'T':
		color = '#4D5D53'
	elif key == 'W':
		color = '#F64A8A'
	elif key == 'V':
		color = '#FADA5E'
	elif key == 'Y':
		color = '#00A86B'
	else:
		color = 'r'

	return color

def plot(masterDict):
	#for N in masterDict:
		N = 15
		dataDict = {}
		ind = np.arange(N)
		width = 0.5
		plotbars = []
		plotbarsnames = []

		for residue in xrange(0,N,1):
			for key in masterDict[N]['positions'][residue]:
				if key not in dataDict:
					dataDict[key] = np.array([0.]*N)
				dataDict[key][residue] = masterDict[N]['positions'][residue][key]
		#totals = [0]*N
		prevKey = ''
		bottomTotal = [0]*N
		for key in dataDict:
			#print 'previous key:', prevKey
			if prevKey == '':
				plotbars.append(plt.bar(ind, 100.*(dataDict[key]/len(masterDict[N]['sequences'])), width, color=color(key)))
				plotbarsnames.append(key)
			else:
				bottomTotal = elementwise_addition(bottomTotal, 100.*(dataDict[prevKey]/len(masterDict[N]['sequences'])))
				plotbars.append(plt.bar(ind, 100.*(dataDict[key]/len(masterDict[N]['sequences'])), width, color=color(key), bottom=bottomTotal))
				plotbarsnames.append(key)

			prevKey = key
		
		plt.title('Positional AA Distributions for loop length: %s' % str(N))
		plt.ylabel('Percent of Amino Acids (%)')
		plt.xlabel('Residue')
		xArray = [str(x) for x in range(0,N)]
		plt.xticks(ind+width/2.,xArray)
		plt.yticks(np.arange(0,101,10))
		plt.ylim([0,103])
		plt.figlegend(plotbars[::-1], plotbarsnames[::-1], loc='right')
		#plt.show()
		plt.savefig('media/CDRH3_length_%s' % str(N), dpi=200)



def main():
	args = sys.argv[1:]

	#Check for proper input arguments- args[0] is loop name
	if not args or not args[0]:
		print 'usage: required filename'
		print args
		sys.exit(1)

	parsedFile = open('parsed/%scontact.txt' % args[0], 'rU')
	loopSequenceList = []
	if parsedFile:
		for line in parsedFile:
			loopSequenceList.append(line[:-1])
	print 'Filtered Number of Sequences: ', len(loopSequenceList)

	masterDict = initPositionInfo(sort_lengths(loopSequenceList))
	masterDict = generatePositionInfo(masterDict, loopSequenceList)
	print masterDict.keys()
	print masterDict[4].keys()
	print masterDict[4]['positions'].keys()
	print masterDict[4]['positions'][1]
	print '\n\n\n\n\n'
	#print json.dumps(masterDict, sort_keys=True, indent=4)
	
	#print masterDict[12]['positions'][0]
	plot(masterDict)
	


if __name__ == '__main__':
	main()
import sys
import re
import operator
import numpy as np
import matplotlib.pyplot as plt

def parseFile(filename):
	file = open(filename ,  'rU')
	list = []
	for line in file:
		match = re.search(r'<p class=\'aasequence\'> (\w+) </p>', line)
		if match:
			list.append(match.group(1))
	file.close()
	print 'Raw Number of Sequences: ', len(list)
	return list

def addToDict(item, gendict):
	if item not in gendict:
		gendict[item] = 1
	else:
		gendict[item] += 1
	return gendict

def sortDict(dictToSort, sortIndex):
	sortedList = sorted(dictToSort.iteritems(), key=operator.itemgetter(sortIndex))
	sortedList.reverse()
	return sortedList

def lengthFilter(min,unfiltList):
	filteredList = []
	for element in unfiltList:
		if len(element) > min and len(element) < 20:
			filteredList.append(element)
	return filteredList

def makeFig(plotList, xIndex, yIndex, title):
	N = len(plotList)
	ind = np.arange(N)
	width = 0.25

	fig = plt.figure()
	ax = fig.add_subplot(111)

	yValues = []
	xValues = []
	totalPoints = 0
	for datapoint in plotList:
		yValues.append(datapoint[yIndex])
		xValues.append(datapoint[xIndex])
		totalPoints += datapoint[yIndex]
	yValues = [ float(point)/float(totalPoints) for point in yValues ]
	rects = ax.bar(ind, yValues, width)
	ax.set_ylabel('Percent abundance')
	ax.set_xlabel('Motif')
	ax.set_title(title)
	ax.set_xticks(ind+width)
	ax.set_xticklabels(xValues)
	#plt.rc('xtick', labelsize=3)

def filterForSpecifiedLength(length, seqs):
	filteredList = []
	for sequence in seqs:
		if len(sequence) is length:
			filteredList.append(sequence)
	return filteredList

def redunancyFilter(unfiltList):
	filteredList = []
	for element in unfiltList:
		if element not in filteredList:
			filteredList.append(element)
	return filteredList

def track(length, seqs):
	merDict = {}
	for sequence in seqs:
		for idx, residue in enumerate(sequence):
			if idx + (length-1) < len(sequence):
				merDict = addToDict(sequence[idx:idx+length],merDict)
	#return sortDict(merDict,1)[0:100]
	return merDict

def lengthTrack(seqs):
	lengthDict = {}
	lengthList = [ len(sequence) for sequence in seqs ]
	for length in lengthList:
		lengthDict = addToDict(length, lengthDict)
	return lengthDict.items()

def calcAAdist(args):
	if not args or not args[0] or not args[1]:
		print 'usage: required filename and mer length'
		print args
		sys.exit(1)
	try:
		parsedFile = open('./parsed/%s' % args[0], 'rU')
	except:
		parsedFile = False
		print 'File not parsed yet'
	
	loopSeqList = []
	if parsedFile:
		for line in parsedFile:
			loopSeqList.append(line[:-1])
	else:
		print 'Parsing & Filtering', args[0]
		loopSeqList = redunancyFilter(lengthFilter(3, parseFile(args[0])))
		newParsedFile = open('./parsed/%s' % args[0], 'w')
		for item in loopSeqList:
			newParsedFile.write("%s\n" % item)

	print 'Filtered Number of Sequences: ', len(loopSeqList)
	lengthDist = lengthTrack(loopSeqList)
	lengthInput = raw_input("Enter length of loop for AAdist (type 'all' for all lengths): ")
	if lengthInput != 'all':
		loopSeqList = filterForSpecifiedLength(int(lengthInput), loopSeqList)

	AAdist = track(int(args[1]), loopSeqList)
	return AAdist
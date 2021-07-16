# Description: Find the ratios of different exons in BrainVar data over time
# Usage: python3 olego_brainvar.py chr2:165308703-165308794 chr2:165309165-165309256 chr2:165309352-165309443 chr2:165310323-165310414 
# Author: Stephan Sanders

import os
import sys
import subprocess
import re

# Quick way of printing variable name and variable, useful when debugging
def printv(x):
	import inspect
	frame = inspect.currentframe().f_back
	s = inspect.getframeinfo(frame).code_context[0]
	r = re.search(r"\((.*)\)", s).group(1)
	print("{} = {}".format(r,x))

# Get the region list
regionList = sys.argv.copy()
regionList.pop(0)
regionStarts = []
regionEnds = []
for region in regionList:
	regionPos = region.split(':')[1]
	regionStarts.append(int(regionPos.split('-')[0]))
	regionEnds.append(int(regionPos.split('-')[1]))

# Location of key files
homePath = os.path.expanduser('~')
sampPcdFileName = f'/data/sample_pcd.txt'
bamDir = f'/data/'
outputFile = f'{bamDir}SCN2A_readCount.txt'
out = open(outputFile, 'w')
headOut = '\t'.join(regionList)
out.write(f'sample\tage_pcd\t{headOut}\n')

# Get the age for each sample (postconceptual days)
sampPcdFile = open(sampPcdFileName, 'r')
headLine = sampPcdFile.readline()

sampPcdDict = {}
sampSexDict = {}
for line in sampPcdFile:
	item = line.strip().split('\t') # HSB272	42.98	6.14	PCW	1	M
	sampPcdDict[item[0]] = item[1]
	sampSexDict[item[0]] = item[5]
	bamFileName = f'{bamDir}{item[0]}.sorted.bam'
	countReg = []
	countReg.append(item[0])
	countReg.append(item[1])
	
	# Get the unique reads for the region
	bamDict = {}
	for region in regionList:
		countReg.append(0)
		command = f'samtools view {bamFileName} region | wc -l'
		result = subprocess.run(['samtools', 'view', bamFileName, region], stdout=subprocess.PIPE)
		bamData = result.stdout.decode('utf-8').split('\n')
		for bamLine in bamData:
			bamDict[bamLine] = 5

	# For each read work out the regions with data
	for bamLine in bamDict:
		bamField = bamLine.strip().split('\t')
		if len(bamField) > 6:
			startPos = int(bamField[3])
			cigar = bamField[5]
			# printv(startPos)
			# printv(cigar)
			end = 0
			exonStarts = [startPos]
			exonEnds = []

			# Find the start and ends of the exons
			while end == 0:
				z = re.search(r'^(\d+)(M|N)', cigar)
				if z:
					num = int(z.group(1))
					let = z.group(2)
					if let == 'M':
						# match
						exonEnds.append(exonStarts[-1]+num)
						cigar = re.sub(r'^(\d+)(M|N)', '', cigar)
					elif let == 'N':
						# skip
						exonStarts.append(exonEnds[-1]+num)
						cigar = re.sub(r'^(\d+)(M|N)', '', cigar)
				elif cigar == '':
					end = 1
				else:
					print(f'ERROR: unexpexted CIGAR string: {cigar}')
					exonStarts = []
					exonEnds = []
					end = 1

			# printv(exonStarts)
			# printv(exonEnds)

			# For each region assess overlap with each exon
			for i in range(len(regionList)):
				startReg = regionStarts[i]
				endReg = regionEnds[i]

				# for each exon
				for j in range(len(exonStarts)):
					startExon = exonStarts[j]
					endExon = exonEnds[j]
					if startReg <= endExon and endReg >= startExon:
						countReg[i+2] = countReg[i+2] + 1

	
	countOut = '\t'.join(map(str, countReg))
	print(countOut)
	out.write(f'{countOut}\n')





	# print(item[0])

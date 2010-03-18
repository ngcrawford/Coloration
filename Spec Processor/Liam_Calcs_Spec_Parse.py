#!/usr/bin/env python
# encoding: utf-8
"""
specparse.py

Created by Nicholas Crawford on 2009-09-24.
Copyright (c) 2009 Boston Univeristy. All rights reserved.
"""

import sys
import os
import glob
import numpy

def getFilenames(path2dir):
	filenames = []
	for infile in glob.glob( os.path.join(path2dir, '*.txt') ):
		filenames.append(infile)
	return filenames

def parseFile(filename,min_reflct,max_reflct):
	""" Read in ocean optics datafile (with headers) and return array of reflectance measurments
		The user can provide min and max reflectance values (e.g. 300-700)
	"""
	fin = open(filename,'r')
	save_sum = 0
	save_count = 0
	save_reflectance = []
	for count, line in enumerate(fin):
		line_parts = line.strip().split()
		if float(line_parts[0]) >= max_reflct: continue
		if float(line_parts[0]) >= min_reflct:
			save_sum += float(line_parts[1])
			save_reflectance.append(float(line_parts[1]))
			save_count += 1
	mean = save_sum/save_count	
	header =  os.path.split(filename)[-1]
	save_reflectance = [header] + save_reflectance
	save_reflectance = numpy.array(save_reflectance)
	return save_reflectance

def calcColorMeasurments(data_array,min_reflct,max_reflct):
	
	results = []
	for index, row in enumerate(data_array):
		x = numpy.asfarray(data_array[index][1:])
		title = data_array[0][0]
		U = ((sum(x[71:286])/216)*75)/((sum(x[71:1200])/1130)*375)
		B = ((sum(x[287:506])/220)*75)/((sum(x[71:1200])/1130)*375)
		G = ((sum(x[507:731])/225)*75)/((sum(x[71:1200])/1130)*375)
		Y = ((sum(x[732:961])/230)*75)/((sum(x[71:1200])/1130)*375)
		R = ((sum(x[962:1200])/239)*75)/((sum(x[71:1200])/1130)*375)
		MU = G-U
		MS = Y-B
		LM = R-G
		C = numpy.sqrt(pow(LM,2)+pow(MS,2)+pow(MU,2))
		C1 = numpy.sqrt(pow(LM,2)+pow(MS,2))
		H1 = numpy.degrees(numpy.arccos(LM/C1))
		Qt = numpy.mean(x[71:1200])/100
	
		results.append([title, U, B, G, Y, R, MU, LM, C, C1, H1, Qt])
	
	results = numpy.array(results)
	print results[1]
	
	
def main():
	path = '/Users/nick/Desktop/Projects/Coloration/Spec Processor/testfiles_no_headers/'
	#path, min_reflct, max_reflct = sys.argv[1:]
	filenames = getFilenames(path)
	min_reflct, max_reflct = (400, 700)
	data_set = []
	for filename in filenames:
		data = parseFile(filename,float(min_reflct),float(max_reflct))
		data_set.append(data)
	data_set = numpy.array(data_set)
	
	calcColorMeasurments(data_set,float(min_reflct),float(max_reflct))

if __name__ == '__main__':
	main()


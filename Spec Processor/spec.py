#!/usr/bin/env python
# encoding: utf-8
"""
specparse.py

Created by Nicholas Crawford on 2009-09-24.
Copyright (c) 2009 Boston Univeristy. All rights reserved.
"""
import os
import sys
import glob
import numpy
import argparse
import itertools

def get_args():
    parser = argparse.ArgumentParser(prog='Spec Parser', description='Convert directory of Spec Files to CSV')
    parser.add_argument('-i','--in-dir', help='the input directory containing the spec files.')
    parser.add_argument('-o','--out-file', help='a csv file to contain the merged specs.')
    parser.add_argument('--header', action='store_true', help='setting this flag will skip headers.')
    parser.add_argument('--min-nm', type=int, default=400, help='lowest nm to include. Default is 400 nm.')
    parser.add_argument('--max-nm', type=int, default=700, help='highest nm to include. Default is 700 nm.')
    parser.add_argument('-v','-verbose', action='store_true', help='write verbose output (non functional)')
    parser.add_argument('--version', action='version', version='%(prog)s beta', help='prints version.')
    args = parser.parse_args()
    
    # CHECK ARGUEMENTS FOR ERRORS
    if os.path.exists(os.path.abspath(args.in_dir)) != True:
        print 'Input directory does not exit.'
        sys.exit() 
    
    if os.path.exists(os.path.abspath(args.out_file)) == True:
        print "\n\t\t\tWARNING: OVERWRITING EXISTING OUTPUT AT %s\n" % (args.out_file)
    
    return args


def getFilenames(path2dir):
    filenames = []
    for infile in glob.glob( os.path.join(path2dir, '*.txt') ):
        filenames.append(infile)
    return filenames

def parseFile(filename, min_reflct, max_reflct, header):
    """ Read in ocean optics datafile (with headers) and return array of reflectance measurments
        The user can provide min and max reflectance values (e.g. 300-700)
    """
    min_reflct = float(min_reflct)
    max_reflct = float(max_reflct)
    in_data_flag = False
    
    fin = open(filename,'r')
    save_reflectance = []
    if header == True:
        for count, line in enumerate(fin):
            if "End" in line: break
                
            if in_data_flag == True:
                line_parts = line.strip().split()
                if float(line_parts[0]) >= max_reflct: 
                    continue
                if float(line_parts[0]) >= min_reflct:
                    save_reflectance.append(float(line_parts[1]))
                    
            if "Begin" in line: 
                in_data_flag = True
    else:
        for count, line in enumerate(fin):
            line_parts = line.strip().split()
            if float(line_parts[0]) >= max_reflct: continue
            if float(line_parts[0]) >= min_reflct:
                save_sum += float(line_parts[1])
                save_reflectance.append(float(line_parts[1]))
                save_count += 1
                
    basename = os.path.splitext(os.path.basename(filename))[0]      
    save_reflectance = numpy.array(save_reflectance)
    return (save_reflectance, basename)

def calcColorMeasurments(data_array,min_reflct,max_reflct):
    
    results = []
    for index, row in enumerate(data_array):
        x = numpy.asfarray(data_array[index])
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
    
        results.append([U, B, G, Y, R, MU, LM, C, C1, H1, Qt])
    
    print "U\tB\tG\tY\tR\tMU\tLM\tC\tC1\tH1\tQt"
    print "%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f" % \
    (U, B, G, Y, R, MU, LM, C, C1, H1, Qt)
        
    results = numpy.array(results)
    return results
    
def saveCSV(data_set, header_list, fout):
    fout = open(fout, 'w')
    
    s = ','.join(itertools.chain(header_list)) + '\n'
    fout.write(s)
    data_set = numpy.transpose(data_set)
    numpy.savetxt(fout, data_set, delimiter=',', fmt='%1.4f')   # X is an array
    pass
        
def main():
    args = get_args()
    filenames = getFilenames(args.in_dir)
    col_headers = []
    data_set = []
    base_dir_name = os.path.split(args.in_dir)[-1]
    
    header_list = []
    for filename in filenames:
      data = parseFile(filename, args.min_nm, args.max_nm, args.header)
      header_list.append(data[1])
      data_set.append(data[0])
    
    data_set = numpy.array(data_set)
    saveCSV(data_set, header_list, args.out_file)
    calcColorMeasurments(data_set,float(args.min_nm),float(args.max_nm))

if __name__ == '__main__':
    main()


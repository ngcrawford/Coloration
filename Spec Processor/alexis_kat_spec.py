#!/usr/bin/env python
# encoding: utf-8
"""
Script for Alexix and Kat to parse spec files.
"""

import csv
import sys
from utilities import spectra
from optparse import OptionParser

def getOptsAndArgs():
    parser = OptionParser()
                
    parser.add_option("-b","--beginning-nm", type = 'int', dest = "beginning_nm", default = 300,
                    help = "First nm to begin reporting output. Default = 300 nm")
                
    parser.add_option("-e","--ending-nm", type = 'int', dest = "ending_nm", default = 700,
                    help = "Last nm to report output. Default = 700 nm")
                
    parser.add_option("-s", "--step-nm", type = 'float', dest = "step_nm",default = 1,
                    help = "Number of nanometers between reported measurements. Default = 1 nm")
                
    options, args = parser.parse_args()
    return (options, args)
    
def checkOptsAndArgs(options, args):
    
    if len(args) == 0:
        print """
        This script requires arguements to run. The general format is:
        
        python alexis_kat_spec.py infile outfile [arguements]
        
        To get help type: python alexis_kat_spec.py -h
        """
        sys.exit()
    
    if len(args) == 1:
        path = args[0]
        outfile = False

    if len(args) == 2:
        path = args[0]
        outfile = args[1]
    
    return outfile
    
def processSpectra(path,options):
    y = spectra.Spectra()
    y.parseSpecFile(path, parsefname = False)
    y.interpolate(options.beginning_nm,options.ending_nm,options.step_nm)
    nm_reflect = zip(y.nanometers,y.reflectance)
    return nm_reflect

def createOutput(outfile,nm_reflect):
    if outfile != False:
        csvWriter = csv.writer(open(outfile, 'w'), delimiter='\t', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        csvWriter.writerow(['nm','reflectance'])
        for row in nm_reflect:
            csvWriter.writerow(row)
    else:
        print '\nnm\treflectance'
        for row in nm_reflect:
            print str(row[0]) + '\t' + str(row[1])
            
def main():
    options, args = getOptsAndArgs()
    outfile = checkOptsAndArgs(options, args)
    nm_reflect = processSpectra(args[0], options)
    createOutput(outfile, nm_reflect)    
    pass

if __name__ == '__main__':
    main()

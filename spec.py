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
import argparse
import itertools
from pylab import *
import numpy
from scipy import interpolate

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(prog='Spec.py', description='Convert directory of Spec Files to CSV, interpolate nanometers, and smooth and plot data.')
    parser.add_argument('-i','--input-dir', help='The input directory containing the spec files.')
    parser.add_argument('-o','--output-file', help='A csv file to contain the merged specs suitable for opening in excel.')
    parser.add_argument('--header', action='store_true', help='Setting this flag will skip headers.')
    parser.add_argument('--min-nm', type=int, default=300, help='Lowest nm to include. Default is 400 nm.')
    parser.add_argument('--max-nm', type=int, default=700, help='Highest nm to include. Default is 700 nm.')
    parser.add_argument('--intrp', type=float, default=1.0, help='Interpolate nm increments. Default is 1 nm.')
    parser.add_argument('-s', '--smooth', action='store_true', help='Add smoothing function. Default is a 100 nm hanning window.')
    parser.add_argument('--window-type', type=complex, choices=['flat','hanning','hamming','bartlett','blackman'])
    parser.add_argument('--window-length', type=int, default=100, help='Window size for smoothing. Longer is more aggressive. Default is 100.')
    parser.add_argument('-p','--plot', action='store_true', help='Produce interactive plots with matplotlib.')
    parser.add_argument('-v','-verbose', action='store_true', help='Write verbose output (non functional).')
    parser.add_argument('--version', action='version', version='%(prog)s beta', help='Print version.')
    args = parser.parse_args()
    
    # CHECK ARGUEMENTS FOR ERRORS
    if os.path.exists(os.path.abspath(args.in_dir)) != True:
        print 'Input directory does not exit.'
        sys.exit() 
    
    if os.path.exists(os.path.abspath(args.out_file)) == True:
        print "\n\t\t\tWARNING: Overwriting existing output at %s\n" % (args.out_file)
    
    if args.window_type == None: args.window_type = 'hanning'
    return args

def getFilenames(path2dir):
    """Parse file names in directory"""
    txt = glob.glob(os.path.join(path2dir, '*.txt'))
    b = glob.glob(os.path.join(path2dir, '*.b'))
    trans = glob.glob(os.path.join(path2dir, '*.transmission'))
    filenames = trans + txt + b
    return filenames

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also: 

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def calcColorMeasurments(data_array):
    """Print Liam's color measurments to SDOUT"""        
    Macedonia_values = []
    Endler_values = []
    data_array = data_array.transpose()
    for author in ['Macedonia', 'Endler']:
        # create slice indices
        if author == 'Macedonia': Qt = (data_array[:,0] >= 325) & (data_array[:,0] <= 700)
        else: Qt = (data_array[:,0] >= 400) & (data_array[:,0] <= 700)
        U = (data_array[:,0] >= 325) & (data_array[:,0] <= 400)
        B = (data_array[:,0] > 400) & (data_array[:,0] <= 475)
        G = (data_array[:,0] > 475) & (data_array[:,0] <= 550)
        Y = (data_array[:,0] > 550) & (data_array[:,0] <= 625)
        R = (data_array[:,0] > 625) & (data_array[:,0] <= 700) 
        
        # do basic calculations
        B = data_array[:,1:].compress(B,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        G = data_array[:,1:].compress(G,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        Y = data_array[:,1:].compress(Y,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        R = data_array[:,1:].compress(R,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        if author == 'Macedonia': U = data_array[:,1:].compress(U,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        else: U = numpy.zeros(data_array[:,1:].shape[-1])
        Qt = data_array[:,1:].compress(Qt,0).sum(0)
        
        # calculate MU, MS, LM, C, and H
        if author == 'Macedonia': MU = G-U
        else: MU = numpy.zeros(data_array[:,1:].shape[-1])
        MS = Y-B
        LM = R-G
        if author == 'Macedonia': C = numpy.sqrt(pow(LM,2)+pow(MS,2)+pow(MU,2))
        else: C = numpy.sqrt(pow(LM,2)+pow(MS,2))
        H = numpy.degrees(numpy.arccos(LM/C))
        if author == 'Macedonia': 
            Macedonia_values = [U,B,G,Y,R,Qt,MU,MS,LM,C,H]
        else: 
            Endler_values = [U,B,G,Y,R,Qt,MU,MS,LM,C,H]
            
    results = numpy.array([Macedonia_values, Endler_values])    
    return results
        
def printCSV(data_set, column_names, row_names):
    """Print files as table"""
    s = "Value,"+','.join(itertools.chain(column_names[1:]))
    print s
    for count, line in enumerate(data_set):
        print row_names[count] + ',' + ','.join(["%.3f" % f for f in line])
    
    # data_set = numpy.transpose(data_set)
    # numpy.savetxt(fout, data_set, delimiter=',', fmt='%1.4f')   # X is an array

def saveCSV(data_set, column_names, fout):
    """Save files as table"""
    fout = open(fout,'w')
    s = 'nanometers,' + ','.join(itertools.chain(column_names[1:])) + '\n'
    fout.write(s)
    data_set = numpy.transpose(data_set)
    numpy.savetxt(fout, data_set, delimiter=',', fmt='%1.4f')   # X is an array

    
def parseFile(filename, min_reflct, max_reflct, header, intrp):
    """ Read in ocean optics datafile (with headers) and return array of reflectance measurments
        The user can provide min and max reflectance values (e.g., 300-700)
    """
    min_reflct = float(min_reflct)
    max_reflct = float(max_reflct)
    in_data_flag = False

    fin = open(filename,'r')
    reflectances = []
    nanometers = []
    if header == True:
        for count, line in enumerate(fin):
            if "End" in line: break

            if in_data_flag == True:
                line_parts = line.strip().split()
                if float(line_parts[0]) >= max_reflct: break
                if float(line_parts[0]) >= min_reflct:
                    nanometers.append(float(line_parts[0]))
                    reflectances.append(float(line_parts[1]))

            if "Begin" in line: 
                in_data_flag = True
    else:
        for count, line in enumerate(fin):
            line_parts = line.strip().split()
            if float(line_parts[0]) > max_reflct: break
            if float(line_parts[0]) >= min_reflct:
                nanometers.append(float(line_parts[0]))
                reflectances.append(float(line_parts[1]))

    basename = os.path.basename(filename)     
    reflectances = numpy.array(reflectances)
    nanometers = numpy.array(nanometers)

    # INTERPOLATE VALUES TO 1 NM INCREMENTS
    tck = interpolate.splrep(nanometers,reflectances,s=0)
    nanometers = numpy.arange(min_reflct,max_reflct,intrp)
    reflectances = interpolate.splev(nanometers,tck,der=0)
    return (numpy.array(reflectances), numpy.array(nanometers), basename)

def plotMean(data_set):
    mean = data_set[1:].mean(axis=0)
    x = data_set.transpose()[:,0]
    mean = data_set[1:].mean(axis=0)
    var = data_set[1:].var(axis=0)
    upper_var = mean + var
    lower_var = mean - var
    xlabel('Nanometers')
    ylabel('Reflectance')
    maxy = mean.max() + 3
    ylim(0, maxy)
    xlim(data_set.transpose()[:,0].min(), data_set.transpose()[:,0].max())
    fill_between(x, upper_var, lower_var, alpha=0.15, color='k')
    plot(x,mean,'k')

def plotThumbs(data_set, header_list):
    data_set = data_set.transpose()
    numb_cols = data_set.shape[1]
    cols = int(sqrt(numb_cols))
    rows = cols + 1
    x = data_set[:,0]
    plt.figure()
    counter = 0
    for r in arange(0,rows):
        for c in arange(0,cols):
            if counter == data_set.shape[1]-1: break
            ax = plt.subplot2grid((rows,cols),(r,c))
            y = data_set[:,counter+1]
            ax.annotate(header_list[counter], xy=(.5, .5),  xycoords='axes fraction',
                            horizontalalignment='center', verticalalignment='center')
            plt.plot(x,y)
            counter += 1
       
def main():

    args = get_args()
    filenames = getFilenames(args.in_dir)
    base_dir_name = os.path.split(args.in_dir)[-1]
    
    # SETUP DATASET
    col_headers = []
    data_set = []
    header_list = []
    for count, filename in enumerate(filenames):
      reflectances, nm, header = parseFile(filename, args.min_nm, args.max_nm, args.header, args.intrp)
      if args.smooth:
          reflectances = smooth(reflectances, args.window_length, args.window_type,)
      header_list.append(header)
      if count == 0:
          data_set.append(nm)
      data_set.append(reflectances)
    data_set = numpy.array(data_set)
    
    # DO SHIT WITH THE DATA!!!!
    if args.plot == True: 
        plotMean(data_set) 
        plotThumbs(data_set,header_list)
        plt.show()
    saveCSV(data_set, header_list, args.out_file)
    macedonia, endler = calcColorMeasurments(data_set)
    row_names = ['U (325-400nm)', 'B (401-475nm)', 'G (476-550nm)', 'Y (551-625)',\
                 'R (626-700)', 'Qt','MU', 'MS', 'LM', 'C', 'H']
    print 'Macedonia Values' 
    printCSV(macedonia, header_list, row_names)
    print '\n' +'Endler Values'
    printCSV(endler, header_list, row_names)
    return data_set

data_set = main()

if __name__ == '__main__':
    try: z = main()
    except KeyboardInterrupt: sys.exit(1) # makes clean control-C exit 




#!/usr/bin/env python
# encoding: utf-8
"""
specparse.py

Created by Nicholas Crawford on 2009-09-24.
Copyright (c) 2009 Boston Univeristy. All rights reserved.

Example: 

python spec.py -i test_data/testfiles_with_headers/ \
-o test_data/out.test  \
--header   \
--plot \
-s \
--window-length 25

"""

import os
import sys
import glob
import argparse
import itertools
import jellyfish
from pylab import *
from scipy import interpolate
from Coloration import Coloration

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser(prog='Spec.py', 
        description='Convert directory of Spec Files to CSV, interpolate nanometers, and smooth and plot data.')
    
    parser.add_argument('-i','--input-dir', 
        help='The input directory containing the spec files.')
    
    parser.add_argument('-o','--output-file', 
        help='A csv file to contain the merged specs suitable for opening in excel.')
    
    parser.add_argument('--header', action='store_true', 
        help='Setting this flag will skip headers.')
    
    parser.add_argument('--min-nm', type=int, default=300, 
        help='Lowest nm to include. Default is 300 nm.')
    
    parser.add_argument('--max-nm', type=int, default=700, 
        help='Highest nm to include. Default is 700 nm.')
    
    parser.add_argument('--intrp', type=float, default=1.0, 
        help='Interpolate nm increments. Default is 1 nm.')
    
    parser.add_argument('-s', '--smooth', action='store_true', 
        help='Add smoothing function. Default is a 100 nm hanning window.')
    
    parser.add_argument('--window-type', type=complex, choices=['flat','hanning','hamming','bartlett','blackman'],
        help= 'Define window smoothing type.')

    parser.add_argument('--window-length', type=int, default=100, 
        help='Window size for smoothing. Longer is more aggressive. Default is 100.')

    parser.add_argument('-p','--plot', action='store_true', 
        help='Produce interactive plots with matplotlib.')

    parser.add_argument('-v','-verbose', action='store_true', 
        help='Write verbose output (non functional).')

    parser.add_argument('--version', action='version', version='%(prog)s beta', 
        help='Print version.')
    
    parser.add_argument("--DDV",action='store_true', 
        help="process files based on dewlap/dorsal/ventral in filename")
    
    args = parser.parse_args()
    
    # CHECK ARGUEMENTS FOR ERRORS
    if os.path.exists(os.path.abspath(args.input_dir)) != True:
        print 'Input directory does not exit.'
        sys.exit() 
    
    if os.path.exists(os.path.abspath(args.output_file)) == True:
        print "\n\t\t\tWARNING: Overwriting existing output at %s\n" % (args.output_file)
    
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

    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def calcColorMeasurments(data_array):
    """Print Liam's color measurments to SDOUT"""        
    Macedonia_values = []
    Endler_values = []
    data_array = data_array.transpose()
    for author in ['Macedonia', 'Endler']:
        # CREATE SLICE INDICES
        if author == 'Macedonia': Qt = (data_array[:,0] >= 325) & (data_array[:,0] <= 700)
        else: Qt = (data_array[:,0] >= 400) & (data_array[:,0] <= 700)
        U = (data_array[:,0] >= 325) & (data_array[:,0] < 400)
        B = (data_array[:,0] >= 400) & (data_array[:,0] < 475)
        G = (data_array[:,0] >= 475) & (data_array[:,0] < 550)
        Y = (data_array[:,0] >= 550) & (data_array[:,0] < 625)
        R = (data_array[:,0] >= 625) & (data_array[:,0] <= 700) 
        
        # DO BASIC CALCULATIONS
        B = data_array[:,1:].compress(B,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        G = data_array[:,1:].compress(G,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        Y = data_array[:,1:].compress(Y,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        R = data_array[:,1:].compress(R,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        if author == 'Macedonia': U = data_array[:,1:].compress(U,0).sum(0) / data_array[:,1:].compress(Qt,0).sum(0)
        else: U = np.zeros(data_array[:,1:].shape[-1])
        Qt = data_array[:,1:].compress(Qt,0).sum(0)
        
        # CALCULATE MU, MS, LM, C, AND H
        if author == 'Macedonia': MU = G-U
        else: MU = np.zeros(data_array[:,1:].shape[-1])
        MS = Y-B
        LM = R-G
        if author == 'Macedonia': C = np.sqrt(pow(LM,2)+pow(MS,2)+pow(MU,2))
        else: C = np.sqrt(pow(LM,2)+pow(MS,2))
        H = np.degrees(np.arccos(LM/C))
        if author == 'Macedonia': 
            Macedonia_values = [U,B,G,Y,R,Qt,MU,MS,LM,C,H]
        else: 
            Endler_values = [U,B,G,Y,R,Qt,MU,MS,LM,C,H]
            
    results = np.array([Macedonia_values, Endler_values])    
    return results
        
def printCSV(data_set, column_names, row_names):
    """Print files as table"""
    
    s = "Value,"+','.join(itertools.chain(column_names))
    print s
    for count, line in enumerate(data_set):
        print row_names[count] + ',' + ','.join(["%.3f" % f for f in line])
    
    # data_set = np.transpose(data_set)
    # np.savetxt(fout, data_set, delimiter=',', fmt='%1.4f')   # X is an array

def saveCSV(data_set, column_names, fout):
    """Save files as table"""
    fout = open(fout,'w')
    s = 'nanometers,' + ','.join(itertools.chain(column_names)) + '\n'
    fout.write(s)
    data_set = np.transpose(data_set)
    np.savetxt(fout, data_set, delimiter=',', fmt='%1.4f')   # X is an array

    
def parseFile(filename, min_reflct, max_reflct, header, intrp):
    """ Read in ocean optics datafile (with headers) and return array of reflectance measurments
        The user can provide min and max reflectance values (e.g., 300-700)
    """
    min_reflct = float(min_reflct)
    max_reflct = float(max_reflct) + 1.0
        
    in_data_flag = False

    fin = open(filename,'r')
    reflectances = []
    nanometers = []
    
    if header == True:
        for count, line in enumerate(fin):

            if "End" in line: break

            if in_data_flag == True:
                line_parts = line.strip().split()
                if float(line_parts[0]) >= min_reflct:
                    nanometers.append(float(line_parts[0]))
                    reflectances.append(float(line_parts[1]))
                if float(line_parts[0]) > max_reflct: break
                
            if "Begin" in line: 
                in_data_flag = True
    
    else:
        for count, line in enumerate(fin):
            line_parts = line.strip().split()
            if float(line_parts[0]) >= min_reflct:
                nanometers.append(float(line_parts[0]))
                reflectances.append(float(line_parts[1]))
            if float(line_parts[0]) > max_reflct: break
            
    basename = os.path.basename(filename)     
    reflectances = np.array(reflectances)
    nanometers = np.array(nanometers)
    
    # INTERPOLATE VALUES TO 1 NM INCREMENTS
    tck = interpolate.splrep(nanometers,reflectances,xb=min_reflct,s=0)
    nanometers = np.arange(min_reflct,max_reflct,intrp)
    reflectances = interpolate.splev(nanometers,tck,der=0)
    return (np.array(reflectances), np.array(nanometers), basename)

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
            ax.set_ylim(0, 50)
            plt.plot(x,y)
            counter += 1


def process_dewlap_dorsal_ventral():

    args = get_args()
    filenames = getFilenames(args.input_dir)
    base_dir_name = os.path.split(args.input_dir)[-1]
    
    # SETUP DATASET
    col_headers = []
    data_set = []
    header_list = []
    
    # Sort files by tissue allowing for missspellings
    organized_by_tissue = {'dewlap':[], 'dorsal':[], 'ventral':[]}
    for count, filename in enumerate(filenames):

        fname = os.path.split(filename)[-1]
        
        for tissue in organized_by_tissue.keys():
           distances = sort([jellyfish.hamming_distance(tissue, item) for item in fname.split("_")])
           if distances[0] <= 2:
                organized_by_tissue[tissue].append(filename)

    
    print organized_by_tissue



    #     reflectances, nm, header = parseFile(filename, args.min_nm, args.max_nm, args.header, args.intrp)
    #     if args.smooth:
    #         reflectances = smooth(reflectances, args.window_length, args.window_type,)            
    #     if count == 0:
    #         data_set.append(nm)
    #     header_list.append(header)
    #     data_set.append(reflectances)
    # data_set = np.array(data_set)
    # # DO #$%^ WITH THE DATA!!!!
    # if args.plot == True: 
    #     plotMean(data_set) 
    #     plotThumbs(data)
       
def main():

    args = get_args()
    
    if args.DDV == True: 
        process_dewlap_dorsal_ventral()
    
    else:
        filenames = getFilenames(args.input_dir)
        base_dir_name = os.path.split(args.input_dir)[-1]
        
        # SETUP DATASET
        col_headers = []
        data_set = []
        header_list = []
        for count, filename in enumerate(filenames):
            reflectances, nm, header = parseFile(filename, args.min_nm, args.max_nm, args.header, args.intrp)
            if args.smooth:
                reflectances = smooth(reflectances, args.window_length, args.window_type,)            
            if count == 0:
                data_set.append(nm)
            header_list.append(header)
            data_set.append(reflectances)
        data_set = np.array(data_set)
        
        # DO #$%^ WITH THE DATA!!!!
        if args.plot == True: 
            plotMean(data_set) 
            plotThumbs(data_set,header_list)
            plt.show()
            plt.savefig('test.png')

        saveCSV(data_set, header_list, args.output_file)
        macedonia, endler = calcColorMeasurments(data_set)
        row_names = ['U (325-399nm)', 'B (40-474nm)', 'G (475-549nm)', 'Y (550-624)',\
                     'R (625-700)', 'Qt','MU', 'MS', 'LM', 'C', 'H']
        
        print 'Macedonia Values' 
        printCSV(macedonia, header_list, row_names)
        
        print '\n' +'Endler Values'
        printCSV(endler, header_list, row_names)
        return data_set


if __name__ == '__main__':

    try: main()
    except KeyboardInterrupt: sys.exit(1) # makes clean control-C exit 




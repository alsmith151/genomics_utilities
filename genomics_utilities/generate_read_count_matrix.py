#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 10:30:06 2019

@author: asmith
"""

import os
import sys
import subprocess as sub
import pandas as pd
import argparse
import itertools
import glob

from pybedtools import BedTool

p = argparse.ArgumentParser()
p.add_argument('peaks', help='Bed file to process')
p.add_argument('--data_loc', help='Directories containing files to process', nargs='+')
p.add_argument('--file_type', help='File extension to process', default='bw', nargs='+')
p.add_argument('--peak_length', 
               help='Determines if peaks are to be standardised to a specific length',
               default=None)
p.add_argument('-o', '--output', help='Counts matrix file name')
p.add_argument('--exclude', help='List of words to exclude if they occur within a file name',
                   default=None, nargs='+')
    
args = p.parse_args()


def standardise_peak_width(peaks, output='bt', length=1000):
    
    '''Uses the midpoint of each interval to standardise each interval to a
       predefined length. Output is returned as a dataframe or a Bedtool object'''
    
    length = int(length)
    print(f'Standardising peaks to {length/1000}kb')
    
    if isinstance(peaks, BedTool):
        peaks = peaks.to_dataframe()
    elif isinstance(peaks, str):
        peaks = BedTool(peaks).to_dataframe()
    
       
    #Generate a new dataframe and determine start and end coordinates based on the midpoints
    df_coords = pd.DataFrame()
    df_coords['chrom'] = peaks['chrom']
    df_coords['midpoint'] = (peaks['start'] + (peaks['end'] - peaks['start']) / 2).astype(int)
    df_coords['start'] = (df_coords['midpoint'] - length / 2).astype(int)
    df_coords['end'] = (df_coords['midpoint'] + length / 2).astype(int)
    df_coords['name'] = peaks['name']
    
    #Drop unwanted columns
    df_coords = df_coords.drop(columns=['midpoint'])
    
    if output == 'df':
        return df_coords
    else:
        return BedTool.from_dataframe(df_coords)

def get_readcounts_bigwig(peaks, data_files, output_raw, output_matrix=None, n_cores=12):
    
    '''Function wraps deeptools multiBigwigSummary to produce readcounts'''
    print(f'Running deeptools multiBigwigSummary with {n_cores} cores')
    
    if isinstance(peaks, pd.DataFrame):
        peaks = BedTool.from_dataframe(peaks[['chrom', 'start', 'end', 'name']])
    
    if not output_matrix:
        tmp_matrix = 'tmp.npz'  
        output_matrix = tmp_matrix
    
    cmd = ['multiBigwigSummary',
                'BED-file',
                '-b',  *data_files,
                '--BED', peaks.fn,
                '-o',  output_matrix,
                '--binSize', args.peak_length,
                '--outRawCounts', output_raw,
                '-p', str(n_cores),
                '-v']
    #print(cmd)
    sub.run(cmd)
    
    if os.path.exists(tmp_matrix):
        os.remove(tmp_matrix)
    
    return output_raw


def format_readcount_matrix(df, peaks):
    
    '''Formats the output of multiBigwigSummary removing unwanted characters
       and also merging peak ids if these are present in the original bed file'''
    
    if isinstance(peaks, BedTool):
        peaks = peaks.to_dataframe()
    
    #Correct column headings
    df.columns = (df.columns.str.replace('#', '')
                            .str.replace("'", ''))
    #Fill NaN
    df = df.fillna(0)
    
    #If peak IDs exist then merge these
    if 'name' in peaks.columns:
        df = (df.merge(peaks[['chrom', 'start', 'end', 'name']],
                       left_on=['chr', 'start', 'end'], 
                       right_on=['chrom', 'start', 'end'])
                    .drop(columns=['chrom', 'chr', 'start', 'end'])
                    .set_index('name')
                )
    
    return df

def main():

    
    if args.peak_length:
        peaks = standardise_peak_width(peaks=args.peaks, length=args.peak_length)
    else:
        peaks = BedTool(args.peaks)
    
    data_files = []
    for ext in args.file_type:
        if len(args.data_loc) > 1:
            for loc in args.data_loc:
                fp = os.path.join(loc, f'*.{ext}')
                fnames = glob.glob(fp)         
                if fnames:
                    data_files.extend(fnames)
        
        else:
            fp = os.path.join(args.data_loc[0], f'*.{ext}')
            fnames = glob.glob(fp)         
            if fnames:
                data_files.extend(fnames)
    
    
    if args.exclude:
        data_files = [fn for fn in data_files if not any(y.lower() in fn.lower() for y in args.exclude)]
        
    read_counts_fn = get_readcounts_bigwig(peaks=peaks, 
                                           data_files=data_files,
                                           output_raw=args.output)
    
    read_counts_df = pd.read_csv(read_counts_fn, sep='\t')
    read_counts_df = format_readcount_matrix(read_counts_df, peaks)
    read_counts_df.to_csv(read_counts_fn)
    
    
if __name__ == '__main__':
    main()
    
    




    
    
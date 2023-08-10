'''
This code was very generously shared by Sebastian Gomez
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import os
import sys
from astropy.nddata import CCDData
from astropy.wcs import WCS

date = '2021.0522'
#date = float(sys.argv[1])
data_directory = '/data/kepler/Archive/rawdata/kepcam/'

name_contains = 'AT_'

def flatten_flat(flat, bias, edge_buffer = 0.2):
    # Flatten the flats, minus bias and normalized to average
    flat_debias = flat - bias
    x_shape, y_shape = flat_debias.shape
    flat_median = np.average(flat_debias[int(x_shape*edge_buffer):int(x_shape*(1 - edge_buffer)),int(y_shape*edge_buffer):int(y_shape*(1 - edge_buffer))])
    final_flat = flat_debias / flat_median
    return final_flat

def correct_data(science, bias_data, g_flat_flat, r_flat_flat, i_flat_flat, z_flat_flat, date, extension = 2):
    # Get Data
    science_file = fits.open(science)
    science_data = science_file[extension].data

    # Get Filter
    science_head = science_file[0].header + science_file[extension].header
    color    = science_head['FILTER']

    # Bias and FLat correct
    if color == 'g': 
        reduced_data = (science_data - bias_data) / g_flat_flat
    if color == 'r': 
        reduced_data = (science_data - bias_data) / r_flat_flat
    if color == 'i': 
        reduced_data = (science_data - bias_data) / i_flat_flat
    if color == 'z': 
        reduced_data = (science_data - bias_data) / z_flat_flat

    # Make directory if it doesn't exist
    print('hi')
    if len(glob.glob('reduced_'+date)) == 0: os.system(f'mkdir reduced_{date}')

    # Format Output Name    
    output_name = f'reduced_{date}/' + science.split('/')[-1][:-5] + '_BiasFlat' + science.split('/')[-2] + '.fits'

    # Save reduced data
    for i in range(4):
        try:
            science_head.remove('COMMENT')
        except:
            pass

    ccddata = CCDData(reduced_data, header=science_head, wcs = WCS(science_file[extension].header),unit='adu')
    ccddata.write(output_name, overwrite = True)

def reduce_data(our_data, data_directory, date, extension = 2):
    # Check which filters were used
    filters = np.array([fits.getval(file, 'FILTER') for file in our_data])
    filters_used = np.unique(filters)

    # Gather calibrations
    print('Gethering Bias and Flats ...')
    bias   = glob.glob(f'{data_directory}{date}/*.BIAS.fits' )
    flat_g = glob.glob(f'{data_directory}{date}/*.FLATg.fits')
    flat_r = glob.glob(f'{data_directory}{date}/*.FLATr.fits')
    flat_i = glob.glob(f'{data_directory}{date}/*.FLATi.fits')
    flat_z = glob.glob(f'{data_directory}{date}/*.FLATz.fits')

    # Reduce bias/flats
    print('Reducing Bias and Flats ...')
    bias_data   = np.median([fits.getdata(bias_file, ext = extension) for bias_file in bias  ], axis = 0)
    g_flat_data = np.median([fits.getdata(flat_file, ext = extension) for flat_file in flat_g], axis = 0)
    r_flat_data = np.median([fits.getdata(flat_file, ext = extension) for flat_file in flat_r], axis = 0)
    i_flat_data = np.median([fits.getdata(flat_file, ext = extension) for flat_file in flat_i], axis = 0)
    z_flat_data = np.median([fits.getdata(flat_file, ext = extension) for flat_file in flat_z], axis = 0)

    g_flat_flat = flatten_flat(g_flat_data, bias_data)
    r_flat_flat = flatten_flat(r_flat_data, bias_data)
    i_flat_flat = flatten_flat(i_flat_data, bias_data)
    z_flat_flat = flatten_flat(z_flat_data, bias_data)

    for science in our_data:
        print('Reducing Data ...', science)
        correct_data(science, bias_data, g_flat_flat, r_flat_flat, i_flat_flat, z_flat_flat, date)

def lookup_data(data_directory, date, name_contains):
    # Check if night exists
    data_exists = False
    night = glob.glob(f'{data_directory}{date}')
    our_data = []
    if len(night) > 0:
        # Check if night has "AT_" data
        our_data = glob.glob(f'{data_directory}{date}/*{name_contains}*')
        if len(our_data) > 0:
            data_exists = True
            print(f'Relevant {name_contains} data found in {date}.')
        else:
            print(f'No data with {name_contains} found in {date}.')
    else:
        print(f'Night {date} does not exist.')
    return data_exists, our_data
data_exists, our_data = lookup_data(data_directory, date, name_contains)
if data_exists:
    reduce_data(our_data, data_directory, date)

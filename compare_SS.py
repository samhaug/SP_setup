#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : compare_SS.py
Purpose : compare time difference and amplitude difference of SS
Creation Date : 02-08-2017
Last Modified : Wed 02 Aug 2017 05:14:05 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from os import listdir
import obspy
import seispy
import scipy

def main():
    std = setup_data('/home/samhaug/work1/SP_sims/20140504_crust/st_R.pk')
    sts = setup_data('/home/samhaug/work1/SP_sims/20140504_D/st_R.pk')
    m1_list,m2_list = M_diff(std,sts)
    return std,sts,m1_list,m2_list

def M_diff(std,sts):
    m1_list = []
    m2_list = []
    for idx,tr in enumerate(std):
        m1_n = scipy.signal.correlate(std[idx].data,sts[idx].data,mode='same')
        m1_d = scipy.signal.correlate(sts[idx].data,sts[idx].data,mode='same')
        m2_n = scipy.signal.correlate(std[idx].data,std[idx].data,mode='same')
        m2_d = scipy.signal.correlate(std[idx].data,sts[idx].data,mode='same')
        m1_list.append(np.sum(m1_n)/np.sum(m1_d))
        m2_list.append(np.sum(m2_n)/np.sum(m2_d))
    return m1_list,m2_list

def setup_data(path_to_data):
    st = obspy.read(path_to_data)
    st.interpolate(10)
    for tr in st:
        tr.data = scipy.signal.hilbert(tr.data).imag
    for idx,tr in enumerate(st):
        st[idx] = seispy.data.phase_window(tr,phase=['SS'],window=(-50,50))
    return st

def setup_sim(path_to_sim):
    st = obspy.read(path_to_sim)
    st.interpolate(10)
    for tr in st:
        tr.data = scipy.signal.hilbert(tr.data).imag
    for idx,tr in enumerate(st):
        st[idx] = seispy.data.phase_window(tr,phase=['SS'],window=(-50,50))
    return st

std,sts,m1_list,m2_list = main()


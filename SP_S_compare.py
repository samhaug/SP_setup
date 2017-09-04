#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : SP_S_compare.py
Purpose : Determine if P branch of SP can be explained with 1D model
Creation Date : 10-08-2017
Last Modified : Mon 04 Sep 2017 04:32:04 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from scipy.signal import hilbert
from os import listdir
import obspy
import seispy
from scipy import correlate
from obspy.taup import TauPyModel
model = TauPyModel(model='prem')

def main():
    gc_bin = 3
    az_bin = 5
    #st = read_waveforms('/home/samhaug/work1/SP_data/2016-08-04-mww63-volcano-islands-japan-region/st_Z.pk')
    #st = read_waveforms('/home/samhaug/work1/SP_data/20170224_172844.a/processed/st_Z.pk')
    st = read_waveforms('/home/samhaug/work1/SP_data/2014-05-04-mw66-south-of-fiji-islands/st_Z.pk')
    #sts = read_waveforms('/home/samhaug/work1/SP_sims/20160804_iasp91/st_Z.pk')
    #sts = read_waveforms('/home/samhaug/work1/SP_sims/20170224_iasp91/st_Z.pk')
    sts = read_waveforms('/home/samhaug/work1/SP_sims/20140504_iasp91/st_Z.pk')
    bounce_list = find_correspond(st)
    S_list,SP_list,s_S_list,s_SP_list = organize_traces(st,sts,bounce_list,
                                                        gc_bin)
    S_map_list,SP_map_list = align_on_SP(S_list,SP_list,s_S_list,
                                         s_SP_list,az_bin,plot=True)
    write_file(S_map_list,SP_map_list,gc_bin,az_bin,st)

def write_file(S_map_list,SP_map_list,gc_bin,az_bin,st):
    name = '{}_{}_{}_{}az_{}gc.dat'.format(str(st[0].stats.starttime.year),
                                       str(st[0].stats.starttime.month),
                                       str(st[0].stats.starttime.day),
                                       str(az_bin),
                                       str(gc_bin))
    with open(name,'w') as f:
        f.write('SP_la, SP_lo, SP_time, SP_az, SP_gcarc, S_la, S_lo, S_time, S_az, S_gcarc \n')
        for ii in range(len(S_map_list)):
            f.write('{} {} {} {} {} {} {} {} {} {}\n'.format(SP_map_list[ii][0],
                                               SP_map_list[ii][1],
                                               SP_map_list[ii][2],
                                               SP_map_list[ii][3],
                                               SP_map_list[ii][4],
                                               S_map_list[ii][0],
                                               S_map_list[ii][1],
                                               S_map_list[ii][2],
                                               S_map_list[ii][3],
                                               S_map_list[ii][4],
                                               ))

def align_on_SP(S_list,SP_list,s_S_list,s_SP_list,az_bin,**kwargs):
    plot = kwargs.get('plot',False)
    SP_map_list = []
    S_map_list = []
    for idx,ii in enumerate(S_list):
        j = min(len(S_list[idx]),len(SP_list[idx]))
        for jj in range(j):
            if abs(S_list[idx][jj].stats.sac['az']-
                   SP_list[idx][jj].stats.sac['az']) > az_bin:
                continue
            else:
                samp = S_list[idx][jj].stats.sampling_rate
                d_s = seispy.data.phase_window(S_list[idx][jj],
                                             phase=['S'],window=(-100,100)).data
                s_s = seispy.data.phase_window(s_S_list[idx][jj],
                                             phase=['S'],window=(-100,100)).data
                S_time = (len(d_s)/2.-np.argmax(correlate(s_s,d_s,mode='same')))/samp
                #SP_list[idx][jj].data = hilbert(SP_list[idx][jj].data).imag
                SP_list[idx][jj].data = SP_list[idx][jj].data
                #s_SP_list[idx][jj].data = hilbert(s_SP_list[idx][jj].data).imag
                s_SP_list[idx][jj].data = s_SP_list[idx][jj].data
                d_sp = seispy.data.phase_window(SP_list[idx][jj],
                                             phase=['SP'],window=(-100,100)).data
                s_sp = seispy.data.phase_window(s_SP_list[idx][jj],
                                             phase=['SP'],window=(-100,100)).data
                SP_time = (len(d_sp)/2.-np.argmax(correlate(s_sp,d_sp,mode='same')))/samp
                print 'SP_az: ',SP_list[idx][jj].stats.sac['az']
                print 'SP_gc: ',SP_list[idx][jj].stats.sac['gcarc']
                print 'S_az: ',S_list[idx][jj].stats.sac['az']
                print 'S_gc: ',S_list[idx][jj].stats.sac['gcarc']
                print SP_time,S_time
                SP_map_list.append((SP_list[idx][jj].stats.sac['stla'],
                                    SP_list[idx][jj].stats.sac['stlo'],
                                    SP_time,
                                    SP_list[idx][jj].stats.sac['az'],
                                    SP_list[idx][jj].stats.sac['gcarc'],
                                    ))
                S_map_list.append((S_list[idx][jj].stats.sac['stla'],
                                    S_list[idx][jj].stats.sac['stlo'],
                                    S_time,
                                    S_list[idx][jj].stats.sac['az'],
                                    S_list[idx][jj].stats.sac['gcarc'],
                                    ))
                if plot == True:
                    fig,ax = plt.subplots(2,1)
                    ax[0].plot(d_s,color='k')
                    ax[0].plot(s_s,color='r')
                    ax[1].plot(d_sp,color='k')
                    ax[1].plot(s_sp,color='r')
                    plt.show()

    return S_map_list,SP_map_list

def organize_traces(st,sts,bounce_list,gc_bin):
    S_list = []
    SP_list = []
    s_S_list = []
    s_SP_list = []
    for ii in bounce_list:
        bounce_S = []
        s_bounce_S = []
        bounce_SP = []
        s_bounce_SP = []
        for idx,tr in enumerate(st):
            gcarc = tr.stats.sac['gcarc']
            if np.abs(gcarc-ii[0]) < gc_bin:
                bounce_S.append(st[idx])
                s_bounce_S.append(sts[idx])
            if np.abs(gcarc-ii[1]) < gc_bin:
                bounce_SP.append(st[idx])
                s_bounce_SP.append(sts[idx])
        S_list.append(bounce_S)
        SP_list.append(bounce_SP)
        s_S_list.append(s_bounce_S)
        s_SP_list.append(s_bounce_SP)
    return S_list, SP_list, s_S_list, s_SP_list

def read_waveforms(path_to_pk):
    st = obspy.read(path_to_pk)
    st.interpolate(10)
    st.filter('lowpass',freq=1./22,zerophase=True).normalize()
    st.sort(['station'])
    return st

def find_correspond(st):
    range_list = []
    for tr in st:
        range_list.append(tr.stats.sac['gcarc'])
    minrange = np.min(range_list)
    maxrange = np.max(range_list)

    bounce_list = []
    for ii in np.arange(95,120,1):
        a = model.get_pierce_points(source_depth_in_km=st[0].stats.sac['evdp'],
                  distance_in_degree=ii,phase_list=['SP'])[0]
        b = np.array([list(ii) for ii in a.pierce])
        i = np.where(b[:,-1]==0)[0]
        bounce_list.append(np.degrees(b[i,-2]))
    return np.array(bounce_list)

main()





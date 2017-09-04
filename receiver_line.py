#!/home/samhaug/anaconda2/bin/python

'''
==============================================================================

File Name : receiver_line.py
Purpose : make line of stations from some source at some azimuth
Creation Date : 20-07-2017
Last Modified : Thu 20 Jul 2017 01:05:04 PM EDT
Created By : Samuel M. Haugland

==============================================================================
'''

import numpy as np
from matplotlib import pyplot as plt
from subprocess import call
from os import listdir
import h5py
import obspy
import seispy
import geopy
from sys import argv

def main():
    st = obspy.read('/home/samhaug/work1/SP_sims/PREM_5s/st_Z.pk')
    stout = obspy.core.stream.Stream()
    start = geopy.Point(st[0].stats.sac['evla'],st[0].stats.sac['evlo'])
    for idx,ii in enumerate(range(180)):
        tr = obspy.core.trace.Trace()
        tr.stats.sac = {}
        d = geopy.distance.VincentyDistance(kilometers=111.195*ii)
        e = d.destination(point=start,bearing=45)
        tr.stats.sac['evla'] = st[0].stats.sac['evla']
        tr.stats.sac['evlo'] = st[0].stats.sac['evlo']
        tr.stats.sac['stla'] = e.latitude
        tr.stats.sac['stlo'] = e.longitude
        tr.stats.station = 'FUCK'
        tr.stats.network = 'II'
        stout.append(tr)
    seispy.mapplot.plot(stout)
    seispy.convert.gemini_stations(stout)

main()

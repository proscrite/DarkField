#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import sys
import linecache
import re
r_wire = 0.125

dat = sys.argv[1]

def check_nelectrons(line, col, col_step):
    try:
        nelectrons = int(linecache.getline(dat, line)[col : col + col_step])
        return nelectrons
    except ValueError:
        return False

col_step = 10
while(True):
    nelectrons = check_nelectrons(27, 25, col_step)
    if nelectrons != False:
        break
    else:
        col_step -= 1
if dat == None:
    dat = "av_out.dat"

# Read file as bytes and get the number of lines printed for each electron avalanche
with open(dat, 'r') as filebytes:
    b_lines = filebytes.read()
    nLines = [26]
# Find float pattern attributed to 'avalanche complete' and append line number to array
    for match in re.finditer('... avalanche complete', b_lines):
        st = match.end()
        nLines.append(b_lines.count('\n', 0, st))

nBufferLines = 4 #Number of lines between each avalanche dataset

print(nLines)
for i in range(nelectrons):
    print('Loading electron trace nr: {}, with Nlines: {}-{}'.format(i, nLines[i]+4, nLines[i+1]))
    ELdat = np.loadtxt(dat, skiprows=nLines[i] + nBufferLines, max_rows=nLines[i+1]-nLines[i]-nBufferLines) #Option max_rows available from np version 1.16
    #nLines[i+1] += nLines[i]

    X_av = ELdat[:,0]*10
    Y_av = ELdat[:,1]*10
    Z_av = ELdat[:,2]*10
    T_av = ELdat[:,3]
    type_av = ELdat[:,4]
    lvl_av = ELdat[:,5]
    fig = plt.figure(1)
    fig.set_figheight(5.0)
    fig.set_figwidth(8.0)
    ax = fig.add_subplot(111)
    ax.plot(-Z_av[type_av == 4],Y_av[type_av == 4],'*',c='blue',markersize=1.5)
    ax.plot(-Z_av[type_av == -1],Y_av[type_av == -1],'--',c='gray',alpha=0.2)
    ax.set_xlim(-0.5,1.5)
    ax.set_ylim(-0.5, 1.0)
    ax.set_xlabel('z (mm)')
    ax.set_ylabel('y (mm)')
    ax.axvspan(-r_wire, r_wire, color='gray')
    plt.text(-0.35, 0.42, "Buffer region", size=14, rotation=90.)
    plt.text(0.8, 0.2, "Active region", size=14, rotation=0.)
    plt.text(0.9, -0.4, "{} excitations".format(len(Z_av[type_av == 4])), size=14, rotation=0.,color='blue')

file_out = sys.argv[2]
if file_out == None: file_out = 'EL_wire'
plt.savefig(file_out+".pdf", bbox_inches='tight')

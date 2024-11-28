# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 20:52:46 2024

@author: DELL
"""

"""readme:
Output the distribution of individual survival. 
In the main code (main_dgds.py) for outputting heatmap data,
set the growth rate at the time of individual extinction and death to -1. 
Therefore, the code classifies areas with a growth rate of -1 and represents 
them as white regions, while areas with a growth rate other than -1 are represented 
as light blue regions."""

import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
n=10
M,SAM=[100,10]
c2=1

dg=[0,2]
ds=[0,1]

num_dup=20  #个数
grid_num=30
dg_range=np.array([dg[0],dg[1]])
grid_dg=np.linspace(dg_range[0],dg_range[1],num=grid_num,endpoint=True)
dg_list=grid_dg

ds_range=np.array([ds[0],ds[1]])
grid_ds=np.linspace(ds_range[0],ds_range[1],num=grid_num,endpoint=True)
ds_list=grid_ds
#print(dg_list)
result=[[[] for i in range(grid_num)] for i in range(grid_num)]
all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
np.shape(all_result[0][0][0])

for dg_num in range(0,grid_num):
    for ds_num in range(0,grid_num):
        max_dup=[[] for i in range(num_dup)]
        for dup in range(0,num_dup):
            max_dup[dup]=np.loadtxt("data/dgs_b1_c1/%s_%s_%s.txt"%(dup,dg_num,ds_num))
            all_result[dup][dg_num][ds_num]=max_dup[dup]


        grate_list=[]
        for dup in range(0,num_dup):
            grate_list.append(max_dup[dup][n,0])


        result[dg_num][ds_num]=max_dup[np.argmax(grate_list)]



all_result_classify1=[[[0 for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
all_result_grate=[[[0 for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
            
            
for dup in range(0,num_dup):
    for dg_num in range(0,grid_num):
        for ds_num in range(0,grid_num):
            item =all_result[dup][dg_num][ds_num]
            all_result_grate[dup][ds_num][dg_num]=item[n][0]
            item1=item[:n,:]   
            if all_result_grate[dup][ds_num][dg_num]==-1:
                t=1
            else:
                for i in range(1,n+1):
                    if all(x==1 for x in item1[:,5]):#ND
                        #all_result_classify[dup][ds_num][dg_num]=1  
                        all_result_classify1[dup][ds_num][dg_num]=4
                # ISD
                    elif all(x==1 for x in item1[-i:,0]) and item1[-1:,5]!=1:#原代码
    
                        all_result_classify1[dup][ds_num][dg_num]=4
                # IGD
                    elif all(x==1 for x in item1[-i:,5]) and item1[-1:,0]!=1:#原代
                        #all_result_classify[dup][ds_num][dg_num]=3
                        all_result_classify1[dup][ds_num][dg_num]=4
                # IGSD
                    elif all(x==1 for x in item1[-i:,0])and all(x==1 for x in item1[-i:,5]):#原代码
                        #all_result_classify[dup][ds_num][dg_num]=3 
                        all_result_classify1[dup][ds_num][dg_num]=4
                    else:
                        #all_result_classify1[dup][ds_num][dg_num]=2
                        all_result_classify1[dup][ds_num][dg_num]=4
                

all_classify_narry1=np.array([[np.array(i) for i in item] for item in all_result_classify1])
all_grate_narry=np.array([[np.array(i) for i in item] for item in all_result_grate])

class_percent_1=[[[] for i in range(grid_num)] for i in range(grid_num)]


for dg_num in range(0,grid_num):
    for ds_num in range(0,grid_num):
        dup_num=np.count_nonzero(~np.isnan(all_classify_narry1[:,dg_num,ds_num]))
        class_percent_1[dg_num][ds_num]=np.count_nonzero(all_classify_narry1[:,dg_num,ds_num]==4)/dup_num #no

class_percent_1_narry=np.array([np.array(i) for i in class_percent_1])



#------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
plt.rc('font',family='Times New Roman')
fig, ax = plt.subplots(1, 1, gridspec_kw={'width_ratios': [1]},figsize=(3, 3))

color=['#ffffff','#f2f3f8','#e5e8f0','#d8dce9','#cbd0e2','#bfc5da','#b2b9d3','#a5adcc']#浅黄色

cmap= LinearSegmentedColormap.from_list('mycmap', color)

norm1 = plt.Normalize(0, 1)
im = ax.imshow(class_percent_1_narry, interpolation=None, origin='lower',cmap=cmap, norm=norm1)


ax.set_ylabel(r'Soma-like death rate, $d_s$',fontsize=15)
ax.set_xlabel(r'Germ-like death rate, $d_g$',fontsize=15)
#ax.set_xlabel(r'Fraction of germ-like, $f_g$',fontsize=13)

#ax.tick_params(direction='inout', length=3, width=1, colors='k')
test0=np.linspace(0,grid_num-1,5,endpoint=True)
ax.set_xticks( test0, minor=False)
test=np.linspace(dg[0],dg[1],5,endpoint=True)
x_label=["%.2f"% i for i in test]
ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax.set_yticks( test0, minor=False)
y_test=np.linspace(ds[0],ds[1],5,endpoint=True)
y_label=["%.2f"% i for i in y_test]
ax.yaxis.set_major_formatter(mpl.ticker.FixedFormatter(y_label))

cbar_ax = fig.add_axes([0.91, 0.3, 0.02, 0.51])
cbar=fig.colorbar(im, cax=cbar_ax, orientation="vertical",norm=norm1, boundaries=None)
cbar.ax.set_ylabel('Persentage of survival', rotation=90,fontsize=15)
cbar.ax.tick_params(labelsize=8,length=2,direction='in')
cbar.outline.set_visible(False)
plt.show()


fig.savefig('./figure/fig2A.pdf' ,bbox_inches='tight')   # save figures

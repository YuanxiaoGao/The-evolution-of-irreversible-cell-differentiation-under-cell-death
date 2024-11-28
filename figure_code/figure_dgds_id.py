# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:44:15 2024

@author: DELL
"""

import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
n=10
M,SAM=[100,10]
c2=1

dg=[0,0.6]
ds=[0,0.6]



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
            max_dup[dup]=np.loadtxt("D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/dgs_b1_c1_small/%s_%s_%s.txt"%(dup,dg_num,ds_num))
            all_result[dup][dg_num][ds_num]=max_dup[dup]


        grate_list=[]
        for dup in range(0,num_dup):
            grate_list.append(max_dup[dup][n,0])


        result[dg_num][ds_num]=max_dup[np.argmax(grate_list)]


all_result_classify=[[[0 for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
all_result_grate=[[[0 for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
            
            
for dup in range(0,num_dup):
    for dg_num in range(0,grid_num):
        for ds_num in range(0,grid_num):
            item =all_result[dup][dg_num][ds_num]
            all_result_grate[dup][ds_num][dg_num]=item[n][0]
            item1=item[:n,:]   
            if all_result_grate[dup][ds_num][dg_num]==-1:
                all_result_classify[dup][ds_num][dg_num]=3
            else:
                for i in range(1,n+1):
                    if all(x==1 for x in item1[:,5]):#ND
                        all_result_classify[dup][ds_num][dg_num]=1 
                # ISD
                    elif all(x==1 for x in item1[-i:,0]) and item1[-1:,5]!=1:#原代码
                    #elif all(x==1 for x in item1[7:9,0]) and item1[-1:,5]!=1:
                        all_result_classify[dup][ds_num][dg_num]=2
                # IGD
                    elif all(x==1 for x in item1[-i:,5]) and item1[-1:,0]!=1:#原代
                    #elif all(x==1 for x in item1[7:9,5]) and item1[-1:,0]!=1:
                        all_result_classify[dup][ds_num][dg_num]=2
                # IGSD
                    elif all(x==1 for x in item1[-i:,0]) and all(x==1 for x in item1[-i:,5]):#原代码
                    #elif all(x==1 for x in item1[7:9,0]) and all(x==1 for x in item1[7:9,5]):#原代码
                        all_result_classify[dup][ds_num][dg_num]=2 
                       
all_classify_narry=np.array([[np.array(i) for i in item] for item in all_result_classify])
all_grate_narry=np.array([[np.array(i) for i in item] for item in all_result_grate])


class_percent_nd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_rd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_did=[[[] for i in range(grid_num)] for i in range(grid_num)]



for dg_num in range(0,grid_num):
    for ds_num in range(0,grid_num):
        dup_num=np.count_nonzero(~np.isnan(all_classify_narry[:,dg_num,ds_num]))
        class_percent_nd[dg_num][ds_num]=np.count_nonzero(all_classify_narry[:,dg_num,ds_num]==1)/dup_num #nd
        class_percent_rd[dg_num][ds_num]=np.count_nonzero(all_classify_narry[:,dg_num,ds_num]==0)/dup_num #rd
        class_percent_did[dg_num][ds_num]=np.count_nonzero(all_classify_narry[:,dg_num,ds_num]==2)/dup_num #id
        
class_percent_nd_narry=np.array([np.array(i) for i in class_percent_nd])
class_percent_rd_narry=np.array([np.array(i) for i in class_percent_rd])
class_percent_did_narry=np.array([np.array(i) for i in class_percent_did])



#------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
plt.rc('font',family='Times New Roman')
fig, ax = plt.subplots(1, 1, gridspec_kw={'width_ratios': [1]},figsize=(3, 3))

color=['#ffffff','#fbf7ed','#f6efdb','#f2e8c9','#eee0b7','#e9d8a6','#e5d094','#e1c982','#dcc170','#d8b95e']#浅黄色

cmap= LinearSegmentedColormap.from_list('mycmap', color)
ax.set_title(r'$ID$',fontsize=14)
norm1 = plt.Normalize(0, 1)
im = ax.imshow(class_percent_did_narry, interpolation=None, origin='lower',cmap=cmap, norm=norm1)


ax.set_ylabel(r'Soma-like death rate, $d_s$',fontsize=13)
ax.set_xlabel(r'Germ-like death rate, $d_g$',fontsize=13)

# xy ticks
# artifical x and y ticks
test0=np.linspace(0,grid_num-1,4,endpoint=True)
ax.set_xticks( test0, minor=False)
test=np.linspace(dg[0],dg[1],4,endpoint=True)
x_label=["%.1f"% i for i in test]
ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax.set_yticks( test0, minor=False)
y_test=np.linspace(ds[0],ds[1],4,endpoint=True)
y_label=["%.1f"% i for i in y_test]
ax.yaxis.set_major_formatter(mpl.ticker.FixedFormatter(y_label))

#ax.tick_params(direction='inout', length=3, width=1, colors='k')
# fig 2 x
#ax.set_xticks( test0, minor=False)
#ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

#ax.set_yticks([] )


cbar_ax2 = fig.add_axes([0.91, 0.3, 0.02, 0.51])
cbar2=fig.colorbar(im, cax=cbar_ax2, orientation="vertical",norm=norm1, boundaries=None)
cbar2.ax.set_ylabel('Persentage of stategy', rotation=90,fontsize=13)
cbar2.ax.tick_params(labelsize=8,length=2,direction='in')
cbar2.outline.set_visible(False)

plt.show()
#fig.savefig('D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/second_fig/app_b0_c1.3.pdf' ,bbox_inches='tight')   # save figures
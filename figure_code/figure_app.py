# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 18:27:25 2024

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



num_dup=10#个数


grid_num=10
dg_range=np.array([dg[0],dg[1]])
grid_dg=np.linspace(dg_range[0],dg_range[1],num=grid_num,endpoint=True)
dg_list=grid_dg

ds_range=np.array([ds[0],ds[1]])
grid_ds=np.linspace(ds_range[0],ds_range[1],num=grid_num,endpoint=True)
ds_list=grid_ds
#print(dg_list[8],ds_list[17])

result=[[[] for i in range(grid_num)] for i in range(grid_num)]
all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
np.shape(all_result[0][0][0])

for dg_num in range(0,grid_num):
    for ds_num in range(0,grid_num):
        max_dup=[[] for i in range(num_dup)]
        for dup in range(0,num_dup):
            max_dup[dup]=np.loadtxt("D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/fl_dgs_b0_c0/%s_%s_%s.txt"%(dup,dg_num,ds_num))
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

fig, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1,1,1]}, figsize=(9, 3.2),)

color0=['#ffffff','#eef8f0','#def0e1','#cde9d2','#bce1c3','#acdab4','#9bd2a5','#8acb96','#7ac387','#69bc78']#浅绿色
color1=['#ffffff','#f0f5fc','#e2eaf9','#d3e0f5','#c4d6f2','#b6cbef','#a7c1ec','#98b7e8','#8aace5','#7ba2e2']#浅蓝色	
color2=['#ffffff','#fbf7ed','#f6efdb','#f2e8c9','#eee0b7','#e9d8a6','#e5d094','#e1c982','#dcc170','#d8b95e']#浅黄色

cmap0 = LinearSegmentedColormap.from_list('mycmap', color0)
cmap1 = LinearSegmentedColormap.from_list('mycmap', color1)
cmap2 = LinearSegmentedColormap.from_list('mycmap', color2)
norm1 = plt.Normalize(0, 1)


im0 = ax[0].imshow(class_percent_nd_narry, interpolation=None, origin='lower',cmap=cmap0, norm=norm1)
im1 = ax[1].imshow(class_percent_rd_narry, interpolation=None, origin='lower',cmap=cmap1, norm=norm1)
im2 = ax[2].imshow(class_percent_did_narry, interpolation=None, origin='lower',cmap=cmap2, norm=norm1)
# title
ax[0].set_title(r'$ND$',fontsize=14)
ax[1].set_title(r'$RD$',fontsize=14)
ax[2].set_title(r'$ID$',fontsize=14)


# x and y label
ax[0].set_xlabel(r'Germ-like death rate, $d_g$',fontsize=13)
ax[0].set_ylabel(r'Soma-like death rate, $d_s$',fontsize=13)
ax[1].set_xlabel(r'Germ-like death rate, $d_g$',fontsize=13)
#ax[2].set_ylabel(r'Soma-like death rate, $d_s$',fontsize=15)
ax[2].set_xlabel(r'Germ-like death rate, $d_g$',fontsize=13)

# xy ticks
ax[0].tick_params(direction='inout', length=3, width=1, colors='k')
ax[1].tick_params(direction='inout', length=3, width=1, colors='k')
ax[2].tick_params(direction='inout', length=3, width=1, colors='k')

# artifical x and y ticks
test0=np.linspace(0,grid_num-1,5,endpoint=True)
ax[0].set_xticks( test0, minor=False)
test=np.linspace(dg_range[0],dg_range[1],5,endpoint=True)
x_label=["%.1f"% i for i in test]
ax[0].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax[0].set_yticks( test0, minor=False)
y_test=np.linspace(ds_range[0],ds_range[1],5,endpoint=True)
y_label=["%.1f"% i for i in y_test]
ax[0].yaxis.set_major_formatter(mpl.ticker.FixedFormatter(y_label))

# fig 2 x
ax[1].set_xticks( test0, minor=False)
ax[1].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))
ax[2].set_xticks( test0, minor=False)
ax[2].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax[1].set_yticks([] )
ax[2].set_yticks([] )

# color bar
cbar_ax0 = fig.add_axes([0.36, 0.3, 0.008, 0.5])
cbar0=fig.colorbar(im0, cax=cbar_ax0, orientation="vertical",norm=norm1, boundaries=None)
cbar0.ax.tick_params(labelsize=8,length=2,direction='in')
cbar0.outline.set_visible(False)

cbar_ax1 = fig.add_axes([0.635, 0.3, 0.008, 0.5])
cbar1=fig.colorbar(im1, cax=cbar_ax1, orientation="vertical",norm=norm1, boundaries=None)
cbar1.ax.tick_params(labelsize=8,length=2,direction='in')
cbar1.outline.set_visible(False)

cbar_ax2 = fig.add_axes([0.91, 0.3, 0.008, 0.5])
cbar2=fig.colorbar(im2, cax=cbar_ax2, orientation="vertical",norm=norm1, boundaries=None)
cbar2.ax.set_ylabel('Persentage of stategy', rotation=90,fontsize=13)
cbar2.ax.tick_params(labelsize=8,length=2,direction='in')
cbar2.outline.set_visible(False)

plt.show()
#fig.savefig('D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/second_fig/app_b0_c0.pdf' ,bbox_inches='tight')
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:01:49 2024

@author: DELL
"""


import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
n=10#分裂次数
M,SAM=[100,10]#M：100个dp0   SAM：每个dp0随机生成10个策略
c2=1#相当于现在的alpha

dg=[0,0.6]#g细胞死亡率
ds=[0,0.6]#s细胞死亡率



num_dup=20  #重复次数


grid_num=30#30*30的矩阵
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
                all_result_classify[dup][ds_num][dg_num]=1
            else:
                for i in range(1,n+1):
                    if all(x==1 for x in item1[:,5]):#ND
                        all_result_classify[dup][ds_num][dg_num]=1 
                # ISD
                    elif all(x==1 for x in item1[-i:,0]) and item1[-1:,5]!=1:
                        all_result_classify[dup][ds_num][dg_num]=2
                # IGD
                    elif all(x==1 for x in item1[-i:,5]) and item1[-1:,0]!=1:
                        all_result_classify[dup][ds_num][dg_num]=3
                # IGSD
                    elif all(x==1 for x in item1[-i:,0]) and all(x==1 for x in item1[-i:,5]):
                        all_result_classify[dup][ds_num][dg_num]=4 
                       
all_classify_narry=np.array([[np.array(i) for i in item] for item in all_result_classify])
all_grate_narry=np.array([[np.array(i) for i in item] for item in all_result_grate])


class_percent_isd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_igd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_igsd=[[[] for i in range(grid_num)] for i in range(grid_num)]



for dg_num in range(0,grid_num):
    for ds_num in range(0,grid_num):
        dup_num=np.count_nonzero(~np.isnan(all_classify_narry[:,dg_num,ds_num]))
        class_percent_isd[dg_num][ds_num]=np.count_nonzero(all_classify_narry[:,dg_num,ds_num]==2)/dup_num #nd
        class_percent_igd[dg_num][ds_num]=np.count_nonzero(all_classify_narry[:,dg_num,ds_num]==3)/dup_num #rd
        class_percent_igsd[dg_num][ds_num]=np.count_nonzero(all_classify_narry[:,dg_num,ds_num]==4)/dup_num #id
        
class_percent_isd_narry=np.array([np.array(i) for i in class_percent_isd])
class_percent_igd_narry=np.array([np.array(i) for i in class_percent_igd])
class_percent_igsd_narry=np.array([np.array(i) for i in class_percent_igsd])
print(class_percent_igsd_narry)


#------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

fig, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1,1,1]},figsize=(9, 3.2))

color0=['#ffffff','#eef8f5','#ddf1eb','#ccebe1','#bbe4d7','#aaddcd','#99d6c3','#88d0b9','#77c9af','#66c2a5']#浅绿色
color1=['#ffffff','#fff2ee','#fee6dc','#fed9cb','#feccb9','#fdc0a8','#fdb396','#fda685','#fc9a73','#fc8d62']#浅蓝色	
color2=['#ffffff','#f2f4f9','#e6eaf3','#d9dfee','#ccd5e8','#c0cae2','#b3c0dc','#a6b5d7','#9aabd1','#8da0cb']#浅黄色

cmap0 = LinearSegmentedColormap.from_list('mycmap', color0)
cmap1 = LinearSegmentedColormap.from_list('mycmap', color1)
cmap2 = LinearSegmentedColormap.from_list('mycmap', color2)

norm1 = plt.Normalize(0,0.05)
norm2 = plt.Normalize(0,0.9)
norm3 = plt.Normalize(0,0.1)

im0 = ax[0].imshow(class_percent_isd_narry, interpolation=None, origin='lower',
                   cmap=cmap0, norm=norm1)
im1 = ax[1].imshow(class_percent_igd_narry, interpolation=None, origin='lower',
                   cmap=cmap1, norm=norm2)
im2 = ax[2].imshow(class_percent_igsd_narry, interpolation=None, origin='lower',
                   cmap=cmap2, norm=norm3)
# title
ax[0].set_title(r'$ISD$',fontsize=15)
ax[1].set_title(r'$IGD$',fontsize=15)
ax[2].set_title(r'$IGSD$',fontsize=15)


# x and y label
ax[0].set_xlabel(r'Germ-like death rate, $d_g$',fontsize=15)
ax[0].set_ylabel(r'Soma-like death rate, $d_s$',fontsize=15)
ax[1].set_xlabel(r'Germ-like death rate, $d_g$',fontsize=15)
ax[2].set_xlabel(r'Germ-like death rate, $d_g$',fontsize=15)

# xy ticks
ax[0].tick_params(direction='inout', length=3, width=1, colors='k')
ax[1].tick_params(direction='inout', length=3, width=1, colors='k')
ax[2].tick_params(direction='inout', length=3, width=1, colors='k')

# artifical x and y ticks
test0=np.linspace(0,grid_num-1,4,endpoint=True)
ax[0].set_xticks( test0, minor=False)
test=np.linspace(dg_range[0],dg_range[1],4,endpoint=True)
x_label=["%.1f"% i for i in test]
ax[0].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax[0].set_yticks( test0, minor=False)
y_test=np.linspace(ds_range[0],ds_range[1],4,endpoint=True)
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
cbar_ax0 = fig.add_axes([0.3555, 0.3, 0.008, 0.5])
cbar0=fig.colorbar(im0, cax=cbar_ax0, orientation="vertical",norm=norm1, boundaries=None)
cbar0.ax.tick_params(labelsize=8,length=2,direction='in')
cbar0.outline.set_visible(False)

cbar_ax1 = fig.add_axes([0.635, 0.3, 0.008, 0.5])
cbar1=fig.colorbar(im1, cax=cbar_ax1, orientation="vertical",norm=norm1, boundaries=None)
cbar1.ax.tick_params(labelsize=8,length=2,direction='in')
cbar1.outline.set_visible(False)

cbar_ax2 = fig.add_axes([0.91, 0.3, 0.008, 0.5])
cbar2=fig.colorbar(im2, cax=cbar_ax2, orientation="vertical",norm=norm1, boundaries=None)
cbar2.ax.set_ylabel('Persentage of stategy', rotation=90,fontsize=15)
cbar2.ax.tick_params(labelsize=8,length=2,direction='in')
cbar2.outline.set_visible(False)

plt.show()
#fig.savefig('D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/fig/fig2_7.pdf' ,bbox_inches='tight')   # save figures
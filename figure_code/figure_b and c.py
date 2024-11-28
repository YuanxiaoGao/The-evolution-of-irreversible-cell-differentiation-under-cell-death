# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 14:13:13 2023

@author: 赵雪艳
"""
#([[s_ss,s_gs,s_gg,g_ss,g_gs,g_gg]])
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
n=10
M,SAM=[100,10]
alpha=1

b=[-2,0]
c=[-2,0]
bx=[0.01,1]
cx=[0.01,1]
b_expo_rangex=np.array([bx[0],bx[1]])
c_expo_rangex=np.array([cx[0],cx[1]])

num_dup=20


grid_num=30
b_expo_range=np.array([b[0],b[1]])
grid_b=np.linspace(b_expo_range[0],b_expo_range[1],num=grid_num,endpoint=True)
b_list=10**grid_b

c_expo_range=np.array([c[0],c[1]])
grid_c=np.linspace(c_expo_range[0],c_expo_range[1],num=grid_num,endpoint=True)
c_list=10**grid_c

result=[[[] for i in range(grid_num)] for i in range(grid_num)]
all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
np.shape(all_result[0][0][0])

for b_num in range(0,grid_num):
    for c_num in range(0,grid_num):
        max_dup=[[] for i in range(num_dup)]
        for dup in range(0,num_dup):
            max_dup[dup]=np.loadtxt("D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/bc_dg0.1_ds0/%s_%s_%s.txt"%(dup,b_num,c_num))
            all_result[dup][b_num][c_num]=max_dup[dup]

        
        grate_list=[]
        for dup in range(0,num_dup):
            grate_list.append(max_dup[dup][n,0])

        
        result[b_num][c_num]=max_dup[np.argmax(grate_list)]


all_result_classify=[[[0 for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
all_result_grate=[[[0 for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
       

for dup in range(0,num_dup):
    for b_num in range(0,grid_num):
        for c_num in range(0,grid_num):
            
            item=all_result[dup][b_num][c_num]
            all_result_grate[dup][c_num][b_num]=item[n][0]
            
            item1=item[:n,:]
            if all_result_grate[dup][c_num][b_num]==-1:
                all_result_classify[dup][c_num][b_num]=3
            else:
                for i in range(1,n+1):
                    if all(x==1 for x in item1[:,5]):#ND
                        all_result_classify[dup][c_num][b_num]=1 
                # ISD
                    elif all(x==1 for x in item1[-i:,0]) and item1[-1:,5]!=1:#原代码
                    #elif all(x==1 for x in item1[7:9,0]) and item1[-1:,5]!=1:
                        all_result_classify[dup][c_num][b_num]=2
                # IGD
                    elif all(x==1 for x in item1[-i:,5]) and item1[-1:,0]!=1:#原代
                    #elif all(x==1 for x in item1[7:9,5]) and item1[-1:,0]!=1:
                        all_result_classify[dup][c_num][b_num]=2
                # IGSD
                    elif all(x==1 for x in item1[-i:,0]) and all(x==1 for x in item1[-i:,5]):#原代码
                    #elif all(x==1 for x in item1[7:9,0]) and all(x==1 for x in item1[7:9,5]):#原代码
                        all_result_classify[dup][c_num][b_num]=2 


all_classify_narry=np.array([[np.array(i) for i in item] for item in all_result_classify])
all_grate_narry=np.array([[np.array(i) for i in item] for item in all_result_grate])


class_percent_nd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_rd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_did=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_no=[[[] for i in range(grid_num)] for i in range(grid_num)]

for b_num in range(0,grid_num):
    for c_num in range(0,grid_num):
        dup_num=np.count_nonzero(~np.isnan(all_classify_narry[:,b_num,c_num]))
        class_percent_nd[b_num][c_num]=np.count_nonzero(all_classify_narry[:,b_num,c_num]==1)/dup_num #ND fraction
        class_percent_rd[b_num][c_num]=np.count_nonzero(all_classify_narry[:,b_num,c_num]==0)/dup_num #rd
        class_percent_did[b_num][c_num]=np.count_nonzero(all_classify_narry[:,b_num,c_num]==2)/dup_num 
        class_percent_no[b_num][c_num]=np.count_nonzero(all_classify_narry[:,b_num,c_num]==3)/dup_num #no
        
class_percent_nd_narry=np.array([np.array(i) for i in class_percent_nd])
class_percent_rd_narry=np.array([np.array(i) for i in class_percent_rd])
class_percent_did_narry=np.array([np.array(i) for i in class_percent_did])
class_percent_no_narry=np.array([np.array(i) for i in class_percent_no])

#------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

fig, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1,1,1]},figsize=(9,3.2))


color0=['#ffffff','#eef8f0','#def0e1','#cde9d2','#bce1c3','#acdab4','#9bd2a5','#8acb96','#7ac387','#69bc78']#浅绿色
color1=['#ffffff','#f0f5fc','#e2eaf9','#d3e0f5','#c4d6f2','#b6cbef','#a7c1ec','#98b7e8','#8aace5','#7ba2e2']#浅蓝色	
color2=['#ffffff','#fbf7ed','#f6efdb','#f2e8c9','#eee0b7','#e9d8a6','#e5d094','#e1c982','#dcc170','#d8b95e']#浅黄色

cmap0 = LinearSegmentedColormap.from_list('mycmap', color0)
cmap1 = LinearSegmentedColormap.from_list('mycmap', color1)
cmap2 = LinearSegmentedColormap.from_list('mycmap', color2)
norm1 = plt.Normalize(0, 1)


#----figure------------------
im0 = ax[0].imshow(class_percent_nd_narry, interpolation=None, origin='lower',cmap=cmap0, norm=norm1)
im1 = ax[1].imshow(class_percent_rd_narry, interpolation=None, origin='lower',cmap=cmap1, norm=norm1)
im2 = ax[2].imshow(class_percent_did_narry, interpolation=None, origin='lower',cmap=cmap2, norm=norm1)
# title
ax[0].set_title(r'$ND$',fontsize=14)
ax[1].set_title(r'$RD$',fontsize=14)
ax[2].set_title(r'$ID$',fontsize=14)

# x and y label
ax[0].set_xlabel(r'Cell differentiation benefit, $b$',fontsize=12)
ax[0].set_ylabel(r'Cell differentiation cost, $c$',fontsize=12)
ax[1].set_xlabel(r'Cell differentiation benefit, $b$',fontsize=12)
ax[2].set_xlabel(r'Cell differentiation benefit, $b$',fontsize=12)

# xy ticks
ax[0].tick_params(direction='inout', length=3, width=1, colors='k')
ax[1].tick_params(direction='inout', length=3, width=1, colors='k')
ax[2].tick_params(direction='inout', length=3, width=1, colors='k')



# artifical x and y ticks
test0=np.linspace(0,grid_num-1,4,endpoint=True)
ax[0].set_xticks( test0, minor=False)
test=np.linspace(b_expo_rangex[0],b_expo_rangex[1],4,endpoint=True)
#x_label=[r'$10^{'+str(int(i))+'}$' for i in test]
x_label=["%.2f"% i for i in test]
ax[0].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax[0].set_yticks( test0, minor=False)
y_test=np.linspace(c_expo_rangex[0],c_expo_rangex[1],4,endpoint=True)
#y_label=[r'$10^{'+str(int(i))+'}$' for i in y_test]
y_label=["%.2f"% i for i in y_test]
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

cbar_ax1 = fig.add_axes([0.63, 0.3, 0.008, 0.5])
cbar1=fig.colorbar(im1, cax=cbar_ax1, orientation="vertical",norm=norm1, boundaries=None)
cbar1.ax.tick_params(labelsize=8,length=2,direction='in')
cbar1.outline.set_visible(False)

cbar_ax2 = fig.add_axes([0.91, 0.3, 0.008, 0.5])
cbar2=fig.colorbar(im2, cax=cbar_ax2, orientation="vertical",norm=norm1, boundaries=None)
cbar2.ax.set_ylabel('Persentege of stategy', rotation=90,fontsize=13)
cbar2.ax.tick_params(labelsize=8,length=2,direction='in')
cbar2.outline.set_visible(False)

plt.show()

fig.savefig('D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/second_fig/fig3_6(dg=0.1  ds=0).pdf' ,bbox_inches='tight')   # save figures
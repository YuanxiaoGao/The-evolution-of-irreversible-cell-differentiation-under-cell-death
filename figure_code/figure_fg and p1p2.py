# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 17:12:54 2024

@author: 赵雪艳
"""

import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
import matplotlib.pyplot as plt

dg_l=[0,0.6]
ds_l=[0,0.6]

b=1
c=1
alpha=1
n=10
M,SAM=[100,10]
num_dup=20   #个数

grid_num=30
dg_range=np.array([dg_l[0],dg_l[1]])
grid_dg=np.linspace(dg_range[0],dg_range[1],num=grid_num,endpoint=True)
dg_list=grid_dg

ds_range=np.array([ds_l[0],ds_l[1]])
grid_ds=np.linspace(ds_range[0],ds_range[1],num=grid_num,endpoint=True)
ds_list=grid_ds

dg=dg_list[29]
ds=ds_list[0]

result=[[[] for i in range(grid_num)] for i in range(grid_num)]
all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
np.shape(all_result[0][0][0])

dg_num=29#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ds_num=0#   0,29 ; 20 25 ;29 0
max_dup=[[] for i in range(num_dup)]
for dup in range(0,num_dup):
    max_dup[dup]=np.loadtxt("D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/dgs_b1_c1_small/%s_%s_%s.txt"%(dup,dg_num,ds_num))
    all_result[dup][dg_num][ds_num]=max_dup[dup]

#----------分类策略--------------

l_ND=[]
l_ID=[]
l_RD=[]
#顺序([[s_ss,s_gs,s_gg,g_ss,g_gs,g_gg]])
for i in range(num_dup):
    item=all_result[i][29][0]#20个策略!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    item1=item[:n,:]#策略矩阵
    #ND
    if all(x==1 for x in item1[:9,5]):#如果前九次g_gg=1
        if item[9,5]==1:#如果最后一次也g_gg=1
            l_ND.append(item)
        else:
            #RD
            l_RD.append(item)
    else:
        #ISD
        if item1[-1,0]==1 and item1[-1,5]!=1:
            l_ID.append(item)
        #IGD
        elif item1[-1,0]!=1 and item1[-1,5]==1:
            l_ID.append(item)
        #IGSD
        elif item1[-1,0]==1 and item1[-1,5]==1:
            l_ID.append(item)
        #RD
        else:
            l_RD.append(item)
print(len(l_RD))#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def Fra(n,dg,ds,dp,b,c,alpha):
    p1=[0]
    p2=[0]
    fgsx=[np.array([1,0])]
    fgs=[np.array([1,0])]
    t=0
    f_g=[1]
    f_s=[0]
    time=[0]
    for j in range(n):
        p_1=dp[j][3]+0.5*dp[j][4]
        p1.append(p_1)
        
        p_2=dp[j][2]+0.5*dp[j][1]
        p2.append(p_2)
    
    for i in range(n):
        fgx=fgs[i][0]*(1-p1[i+1])+fgs[i][1]*p2[i+1]
        fsx=fgs[i][1]*(1-p2[i+1])+fgs[i][0]*p1[i+1]
        fgsx.append(np.array([fgx,fsx]))#每次分裂后的gs细胞比例
        
        cost=1+c*((fgsx[i][0]*p1[i+1])+(alpha*fgsx[i][1]*p2[i+1]))
        bene=1+b*fgsx[i][1]
        time.append(cost/bene)
        t_i=cost/bene
        t=t+t_i
        
        fg=(fgsx[i+1][0]*(1-(dg*time[i+1])))/(fgsx[i+1][0]*(1-(dg*time[i+1]))+fgsx[i+1][1]*(1-(ds*time[i+1])))
        fs=(fgsx[i+1][1]*(1-(ds*time[i+1])))/(fgsx[i+1][0]*(1-(dg*time[i+1]))+fgsx[i+1][1]*(1-(ds*time[i+1])))
             
        if 0<=fg<=1 and 0<=fs<=1:
            fg=fg
            fs=fs
        else:
            if fg<0:
                fg=0
                fs=1
            else:
                fg=1
                fs=0
        f_g.append(fg)
        f_s.append(fs)
        fgs.append(np.array([fg,fs]))#每次分裂后加上死亡的gs细胞比例       
    return[p1,p2,f_g,f_s]

#所有策略的p11,p22,fgg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
p11_rd=[]
p22_rd=[]
f_gg_rd=[]
f_ss_rd=[]
for i in range(len(l_RD)):
    p_11=Fra(n,dg,ds,l_RD[i],b,c,alpha)[0]
    p_22=Fra(n,dg,ds,l_RD[i],b,c,alpha)[1]
    fgg=Fra(n,dg,ds,l_RD[i],b,c,alpha)[2]
    fss=Fra(n,dg,ds,l_RD[i],b,c,alpha)[3]
    p11_rd.append(p_11)
    p22_rd.append(p_22)
    f_gg_rd.append(fgg)
    f_ss_rd.append(fss)
avg_p1=np.mean(p11_rd,axis=0)
avg_p2=np.mean(p22_rd,axis=0)
avg_fg=np.mean(f_gg_rd,axis=0)
avg_fs=np.mean(f_ss_rd,axis=0)


plt.rc('font',family='Times New Roman')
fig,ax1=plt.subplots()
x=np.linspace(0,11,11)
ls_w =0.3
for i in range(len(p11_rd)):
    ax1.plot(x,p11_rd[i],linewidth=ls_w,c='#4898D4',linestyle=':')
    ax1.plot(x,p22_rd[i],linewidth=ls_w,c='#65D838',linestyle='--')
    ax1.plot(x,f_gg_rd[i],linewidth=ls_w,c='#DF662E',linestyle='-.')
    ax1.plot(x,f_ss_rd[i],linewidth=ls_w,c='#ff42b3',linestyle='-.')


ax1.plot(x,avg_p1,c='#4898D4',linestyle=':',linewidth=3,label=r'$g_{g \rightarrow s}$')
ax1.plot(x,avg_p2,c='#65D838',linestyle='--',linewidth=3,label=r"$s_{s \rightarrow g}$")
ax1.plot(x,avg_fg,c='#DF662E',linestyle='-.',label=r"$f_g$")
ax1.plot(x,avg_fs,c='#ff42b3',linestyle='-',label=r"$f_s$")
#  ,,,
ax1.scatter(x,avg_p1,c='#4898D4',marker='*',linewidth=3)
ax1.scatter(x,avg_p2,c='#65D838',marker='*',linewidth=3)
ax1.scatter(x,avg_fg,c='#DF662E')
ax1.scatter(x,avg_fs,c='#ff42b3')

#ax1.set_ylabel(r'$ g_{g \rightarrow s}$',fontsize=15)
#ax1.set_xlabel(r'$ s_{s \rightarrow g}$',fontsize=15)
#ax1.legend(loc='best')

ax1.set_title(r'$RD$',fontsize=15)

ax1.set_xlim(-0.15,11.15)
#ax1.set_ylim(-0.05,1.05)
#y_values=
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
x_ticks=ax1.set_xticks([0,2.2,4.4,6.6,8.8,11])#位置
labels=ax1.set_xticklabels(['0','2','4','6','8','10'],fontsize=13)
y_ticks=ax1.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])#位置
labels1=ax1.set_yticklabels(['0','0.2','0.4','0.6','0.8','1.0'],fontsize=13)
#ax1.yticks(size=12)
ax1.set_xlabel(r'Cell division times, $i$',fontsize=15)



plt.show()

#fig.savefig('D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/second_fig/fig2_4.pdf' ,bbox_inches='tight')   # save figures 
    
    
    

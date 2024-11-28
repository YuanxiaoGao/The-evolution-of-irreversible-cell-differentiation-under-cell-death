# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:35:01 2024

@author: DELL
"""
#import sys
import numpy as np
#np.set_printoptions(threshold=sys.maxsize)
#np.get_printoptions()['linewidth']
#np.set_printoptions(linewidth=2000)
n=10
alpha=1


b=1
c=0.5

dg_fanwei=[0,0.6]
ds_fanwei=[0,0.6]

num_dup=10
grid_num=13

dg_range=np.array([dg_fanwei[0],dg_fanwei[1]])
grid_dg=np.linspace(dg_range[0],dg_range[1],num=grid_num,endpoint=True)
dg_list=grid_dg

ds_range=np.array([ds_fanwei[0],ds_fanwei[1]])
grid_ds=np.linspace(ds_range[0],ds_range[1],num=grid_num,endpoint=True)
ds_list=grid_ds

def Fration(n,dp,dg,ds,b,c,alpha):
    p12=[np.array([0,0])]
    fgsx=[np.array([1,0])]
    fgs=[np.array([1,0])]
    t=0
    time=[0]
   
    for j in range(n):
        p1=dp[j][3]+0.5*dp[j][4]
        p2=dp[j][2]+0.5*dp[j][1]
        p12.append(np.array([p1,p2]))
   
    for i in range(n):
        fgx=fgs[i][0]*(1-p12[i+1][0])+fgs[i][1]*p12[i+1][1]
        fsx=fgs[i][1]*(1-p12[i+1][1])+fgs[i][0]*p12[i+1][0]
        fgsx.append(np.array([fgx,fsx]))#每次分裂后的gs细胞比例
        
        cost=1+c*((fgsx[i][0]*p12[i+1][0])+(alpha*fgsx[i][1]*p12[i+1][1]))
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
        fgs.append(np.array([fg,fs]))#每次分裂后加上死亡的gs细胞比例   
        
    return[p12,fgs,time,t,fgsx]

    
def Num(n,dp,dg,ds):
    N=1
    fgsx=Fration(n,dp,dg,ds,b,c,alpha)[4]      
    time=Fration(n,dp,dg,ds,b,c,alpha)[2]
    for i in range(n): 
        num=1-(time[i+1]*((fgsx[i+1][0]*dg)+(fgsx[i+1][1]*ds))) 
        if num>0:
            N=N*num
        else:
            N=0
    return (2**n)*N



result=[[[] for i in range(grid_num)] for i in range(grid_num)]
all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
np.shape(all_result[0][0][0])


item_juzhen=[]
l_ID=[] 
l_NOT=[]
l_si=[]
l_size=[] 
for dg_num in range(0,grid_num):
    dg=dg_list[dg_num]
    for ds_num in range(0,grid_num):
        ds=ds_list[ds_num]
        max_dup=[[] for i in range(num_dup)]
        for dup in range(0,num_dup):
            max_dup[dup]=np.loadtxt("D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/qipao_04/%s_%s_%s.txt"%(dup,dg_num,ds_num))
            all_result[dup][dg_num][ds_num]=max_dup[dup]
            
            item=all_result[dup][dg_num][ds_num]#num_dup个策略
            item_juzhen.append(item[:n,:])#num_dup个纯策略
        #print(item_juzhen)
        for i in range(len(item_juzhen)):#存储所有的增长率为l_ID
            item1=item_juzhen[i]
            if all(x==1 for x in item1[:9,5]):#如果前九次g_gg=1
                if item1[-1,0]==1 and item1[-1,5]!=1:
                    #ISD
                    l_ID.append(item1)
                else:
                    l_NOT.append(item1)
            else:
                #ISD
                if item1[-1,0]==1 and item1[-1,5]!=1:
                    l_ID.append(item1)
                #IGD
                elif item1[-1,0]!=1 and item1[-1,5]==1:
                    l_ID.append(item1)
                #IGSD
                elif item1[-1,0]==1 and item1[-1,5]==1:
                    l_ID.append(item1)
                #非ID
                else:
                    l_NOT.append(item1)
        if len(l_ID)==0:
            size=0
            l_size.append(size)
        else:
            for j in range(len(l_ID)):
                id_sta=l_ID[j]
                #s=Num(n,id_sta,dg,ds)
                s=Fration(n,id_sta,dg,ds,b,c,alpha)[1][n][0]
                l_si.append(s)
            size=np.mean(l_si)*100#求比例×100
            siz=round(size)
            l_size.append(siz)
        item_juzhen.clear()
        l_ID.clear()
        l_si.clear()
print(l_size)          
            
            
            
            
 


            
        

        
            

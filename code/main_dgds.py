# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 15:06:46 2024

@author: DELL
"""
"""We set the cell differentiation benefits b and differentiation costs c as 
constant values (b and c), and the mortality rates of g cells and s cells 
(d_g and d_s) as the circulation volume. We repeated the calculation of k 
optimal strategies for each cycle, and finally obtain k * b * c optimal strategies's data."""

import numpy as np
import random

n=10
portion=[80,10,10]
M=100
SAM=10
alpha=1

b=0
c=1.3

for k in range(10):
    num_duplicate=k
    for i in range(10):
        dg_num=i
        for j in range(10):
            ds_num=j
            grid_num=10
            dg_range=np.array([0,0.6])
            grid_dg=np.linspace(dg_range[0],dg_range[1],num=grid_num,endpoint=True)
            dg_list=grid_dg
            dg=dg_list[dg_num]
            
            ds_range=np.array([0,0.6])
            grid_ds=np.linspace(ds_range[0],ds_range[1],num=grid_num,endpoint=True)
            ds_list=grid_ds
            ds=ds_list[ds_num]
            
            def Fration(n,dp,dg,ds,b,c,alpha):
                p12=[np.array([0,0])]
                fgsx=[np.array([1,0])]
                fgs=[np.array([1,0])]
                t=0
                time=[0]
                
                for j in range(n):
                    p1=dp[j][1][0]+0.5*dp[j][1][1]
                    p2=dp[j][0][2]+0.5*dp[j][0][1]
                    p12.append(np.array([p1,p2]))
                    
                for i in range(n):
                    fgx=fgs[i][0]*(1-p12[i+1][0])+fgs[i][1]*p12[i+1][1]
                    fsx=fgs[i][1]*(1-p12[i+1][1])+fgs[i][0]*p12[i+1][0]
                    fgsx.append(np.array([fgx,fsx]))
                    
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
                    fgs.append(np.array([fg,fs])) 
                    
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


            def Growth_num(n,dp,c,alpha,b):
                fgs=Fration(n,dp,dg,ds,b,c,alpha)[1]
                NUM=Num(n,dp,dg,ds)
                if NUM==0:
                    lamb=-1
                else:
                    dno=np.log(NUM*fgs[n][0])
                    numo=Fration(n,dp,dg,ds,b,c,alpha)[3]
                    lambdaa=np.nan_to_num((dno)/(numo),nan=-1)
                    if lambdaa<=0:
                        lamb=-1
                    else:
                        lamb=lambdaa
                return lamb

            
            def Dp_space():
                size=11
                q=np.linspace(0.0,1.0,size)
                for i in range(size):
                    q[i]=round(q[i],2)
                initial_para=[]
                for i in range(size):
                    for j in range(size):
                        for k in range(size):
                            for l in range(0,size):
                                if round(1-q[i]-q[j],2)>=0 and round(1-q[k]-q[l],2)>=0:
                                    initial_para.append(np.array([[q[i],1-q[i]-q[j],q[j]],[q[k],1-q[k]-q[l],q[l]]]))
                return initial_para
            n_dp=4356

            
            Dp_pool=[[] for i in range(3)]  # RD; ID; ND
            for stra in Dp_space():
                if stra[1][2]==1:
                    Dp_pool[2].append(stra)       # ND
                elif stra[0][0]==stra[1][2]==1:
                    Dp_pool[1].append(stra)       # IGSD
                elif stra[0][0]==1 and stra[1][2]<1:
                    Dp_pool[1].append(stra)       # ISD
                elif stra[0][0]<1 and stra[1][2]==1:
                    Dp_pool[1].append(stra)       # IGD
                else:
                    Dp_pool[0].append(stra)#RD
            

            def Neighbour_dps(dp):
                delta_dp=np.array([[1,-1,0],[1,0,-1],[-1,1,0],[0,1,-1],[-1,0,1],[0,-1,1]])
                dp_set=[dp]
                for i in range(2):
                    for j in range(6):
                        if i==0:
                            new0=dp[i]+0.1*delta_dp[j]
                            new1=dp[1]
                            new=np.array([new0,new1])
                        else:
                            new1=dp[i]+0.1*delta_dp[j]
                            new0=dp[0]
                            new=np.array([new0,new1])
                        if np.all(new>=0):
                            dp_set.append(new)
                return dp_set
            
            
            def Dynamic_strategy(dp0,n):
                dpn=[dp0]
                initial_dp=dp0
                for i in range(n-1):
                    neighbour_dps=np.array(Neighbour_dps(initial_dp))
                    random_num=random.randint(0,len(neighbour_dps)-1)  
                    random_neighbour=neighbour_dps[random_num]
                    dpn.append(random_neighbour)
                    initial_dp=random_neighbour
                return dpn
            
            def Dynamic_stra_set(dp0,SAM,n):
                dy_sta_set=[]
                for i in range(SAM):
                    one_sta=Dynamic_strategy(dp0,n)
                    one_sta=np.array(one_sta)
                    dy_sta_set.append(one_sta)
                return dy_sta_set
            
            
            def Dp_lambda(dp0,c,alpha,b,SAM,n):
                growth_table=[]
                dynamic_strategy=Dynamic_stra_set(dp0,SAM,n) 
    
                for i in range(len(dynamic_strategy)):
                    dp=dynamic_strategy[i]
                    lambda0=Growth_num(n,dp,c,alpha,b)
                    growth_table.append(np.array([i,lambda0]))
               
                array_data=np.array(growth_table)  
                growth_max=np.nanmax(array_data[:,1])   #
                index_growth_max=np.where(array_data==growth_max)  
                dp_max=dynamic_strategy[index_growth_max[0][0]]
    
                return [dp_max,growth_max]
            
            
            
            
            def Dp_lambda_set(neighbour_dp_squeue,grate,c,alpha,b,n):
                growth_table=[np.array([0,grate])]
                for i in range(1,len(neighbour_dp_squeue)):
                    dp=neighbour_dp_squeue[i]
                    lambda0=Growth_num(n,dp,c,alpha,b)
                    growth_table.append(np.array([i,lambda0]))
                growth_tab=np.array(growth_table)
                
                growth_max=np.nanmax(growth_tab[:,1])
                index_growth_max=np.where(growth_tab==growth_max)
                dp_max=neighbour_dp_squeue[index_growth_max[0][0]]
                return [dp_max,growth_max]
            
            
            
            def Optimal_one_step(ddp0,grate,c,alpha,b,n):
                neighbour_dps=np.array(Neighbour_dps(ddp0[-1]))
                neighbour_dy_set=[ddp0]
                for i in range(len(neighbour_dps)):
                    new=neighbour_dps[i]
                    new_dyna=    np.append(ddp0[1:],[new],axis=0)
                    neighbour_dy_set.append(new_dyna)
                optimal=Dp_lambda_set(neighbour_dy_set,grate,c,alpha,b,n)
                return optimal
            
            
            def Optimal_one_step_l(ddp0,grate,c,alpha,b,n):
                neighbour_dps=np.array(Neighbour_dps(ddp0[0]))
                neighbour_dy_set=[ddp0]
                for i in range(len(neighbour_dps)):
                    new=neighbour_dps[i]
                    new_dyna=    np.append([new],ddp0[:-1],axis=0)
                    neighbour_dy_set.append(new_dyna)
                optimal=Dp_lambda_set(neighbour_dy_set,grate,c,alpha,b,n)
                return optimal
            
            
            def Local_Optimal(ddp0,grate,c,alpha,b,n):
                dyna_stra_list=[ddp0]
                dyna_rate_list=[grate]
                one_step=Optimal_one_step(ddp0,grate,c,alpha,b,n)
                one_step_l=Optimal_one_step_l(ddp0,grate,c,alpha,b,n)
                
                if one_step[1]>one_step_l[1]:
                    dyna_stra_list.append(one_step[0])
                    dyna_rate_list.append(one_step[1])
                    while abs(dyna_rate_list[-1]-dyna_rate_list[-2])>10**(-8):
                        new_optimal_dp=dyna_stra_list[-1]
                        grate=dyna_rate_list[-1]
                        new_optimal=Optimal_one_step(new_optimal_dp,grate,c,alpha,b,n)
                        dyna_stra_list.append(new_optimal[0])
                        dyna_rate_list.append(new_optimal[1])
                else:
                    dyna_stra_list.append(one_step_l[0])
                    dyna_rate_list.append(one_step_l[1])
                    while abs(dyna_rate_list[-1]-dyna_rate_list[-2])>10**(-8):
                        new_optimal_dp=dyna_stra_list[-1]
                        grate=dyna_rate_list[-1]
                        new_optimal=Optimal_one_step_l(new_optimal_dp,grate,c,alpha,b,n)
                        dyna_stra_list.append(new_optimal[0])
                        dyna_rate_list.append(new_optimal[1])
                if len(dyna_rate_list)==2:
                    return [[ddp0],[grate]]
                else:
                    return [dyna_stra_list,dyna_rate_list]
            
            def Golbal_initial_ds(M,c,alpha,b,n,SAM):
                initial_dps=[]
                for i in range(3):  
                    M_1=portion[i]  
                    sta=random.choices(Dp_pool[i],k=M_1)
                    initial_dps.extend(sta)
                    
                neighbour_optimal=[]
                for dp0 in initial_dps:
                    optimal_dp_grate=Dp_lambda(dp0,c,alpha,b,SAM,n)
                    neighbour_optimal.append(optimal_dp_grate)
                neighbour_optimal=np.array(neighbour_optimal,dtype=object)
                max_grate=np.nanmax(neighbour_optimal[:,1])
                op_index=np.where(neighbour_optimal[:,1]==max_grate)
                if  np.shape(op_index[0])[0]==1 or n==1:
                    return neighbour_optimal[op_index[0][0]]
                else:
                    multiple_ddy0=[]
                    for index in op_index[0]:
                        ddp_flatten=np.array([neighbour_optimal[index][0][i].flatten()    for i in range(n)])
                        multiple_ddy0.append(ddp_flatten)
                    multiple_ddy0=np.array(multiple_ddy0)
                   
                    if np.all(multiple_ddy0[:,:,-1].flatten()==1):
                        return neighbour_optimal[op_index[0][0]]
                    else:
                        result_list=[]
                        for index in op_index[0]:
                            result_list.append(neighbour_optimal[index][0])
                            
                        result_list.append(max_grate)
                        return np.array(result_list,dtype=object)
            
            def Golbal_ddp(M,c,alpha,b,SAM,n):
                local_ddp_grate=Golbal_initial_ds(M,c,alpha,b,n,SAM)
                grate=local_ddp_grate[-1]
                num=np.shape(local_ddp_grate)[0]
               
                if num==2:
                    ddp0=local_ddp_grate[0]
                    history_dynamic=Local_Optimal(ddp0,grate,c,alpha,b,n)
                    
                    return [history_dynamic[0][-1],history_dynamic[1][-1]]
                else:
                    ddp_list=[]
                    grate_list=[]
                    for i in range(num-1):
                        ddp0=local_ddp_grate[i]
                        history_dynamic=Local_Optimal(ddp0,grate,c,alpha,b,n)
                        ddp_list.append(history_dynamic[0][-1])
                        grate_list.append(history_dynamic[1][-1])
                    grate_list=np.array(grate_list)
                    max_grate=np.nanmax(grate_list)
                    op_index=np.where(max_grate==max_grate)
                        
                    if np.shape(op_index[0])[0]==1:
                        return [ddp_list[op_index[0][0]],max_grate]
                    else:
                        ddp_global_list=[]
                        for item in op_index[0]:
                            ddp_global_list.append(ddp_list[item])
                        ddp_global_list.append(max_grate)
                        return ddp_global_list
            
            
            def Output_global_ddp(M,c,alpha,b,SAM,n):
                data_ddp_grate=Golbal_ddp(M,c,alpha,b,SAM,n)
                num_global_ddp=len(data_ddp_grate)
                result=[]
               
                if num_global_ddp==2:
                    for i in data_ddp_grate[0]:
                        result.append(i.flatten())
                    result.append(np.array([data_ddp_grate[-1],0,0,0,0,0]))
                    result=np.array(result)
                else:
                    for i in range(num_global_ddp-1):
                        ddp_g=data_ddp_grate[i]
                        for j in ddp_g:
                            result.append(j.flatten())
                    result.append(np.array([data_ddp_grate[-1],0,0,0,0,0]))
                    result=np.array(result)
                return np.array(result)
    
    
            result=Output_global_ddp(M,c,alpha,b,SAM,n)
            with open('fl_dgs_b0_c1.3/%s_%s_%s.txt'%(num_duplicate,dg_num,ds_num),'wb') as f:
                np.savetxt(f,result,fmt='%1.8f')

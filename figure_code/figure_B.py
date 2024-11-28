# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 23:32:55 2024

@author: DELL
"""




import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

c1='w' #blue
c2='#6CAAD1' #any given color
n=100

color_list=[]
for x0 in range(n+1):
    color_list.append(colorFader(c1,c2,x0/n))

color=np.array(color_list)
len(color)


import pandas as pd
import matplotlib.pyplot as plt
plt, ax=plt.subplots(1, 1, figsize=(5,5))
d = {
'Germ-like death rate, $d_g$':[0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0,
      0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
      0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
      0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,
      0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,
      0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,
      0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,
      0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,0.35,
      0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,
      0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45,
      0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
      0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,
      0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6], 
'Soma-like death rate, $d_s$':[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
      0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6],
'size':[0, 0, 0, 647, 544, 464, 451, 448, 419, 339, 0, 0, 0, 
        0, 0, 0, 490, 416, 378, 317, 314, 311, 310, 0, 0, 0, 
        0, 0, 0, 396, 324, 261, 227, 225, 211, 208, 202, 0, 0, 
        0, 0, 0, 0, 0, 208, 162, 145, 126, 125, 122, 110, 0, 
        0, 0, 0, 0, 221, 171, 131, 100, 88, 71, 70, 69, 65, 
        0, 0, 0, 0, 0, 0, 110, 0, 65, 51, 43, 40, 39, 
        0, 0, 0, 0, 0, 0, 87, 70, 0, 38, 29, 25, 23, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20, 17, 13, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 8, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 8, 4, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 
        0, 0, 216, 0, 0, 0, 0, 0, 0, 22, 0, 0, 0],
'density':[0, 0, 0, 80, 81, 85, 92, 96, 98, 98, 0, 0, 0, 
           0, 0, 0, 72, 75, 82, 83, 91, 93, 98, 0, 0, 0, 
           0, 0, 0, 46, 54, 69, 82, 89, 92, 98, 99, 0, 0, 
           0, 0, 0, 0, 0, 72, 75, 86, 88, 96, 97, 99, 0, 
           0, 0, 0, 0, 0, 61, 74, 76, 85, 86, 94, 98, 99, 
           0, 0, 0, 0, 0, 0, 56, 59, 68, 82, 90, 92, 99, 
           0, 0, 0, 0, 0, 0, 50, 54, 60, 74, 83, 92, 94, 
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 83, 87, 92, 
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 86, 88, 
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 54, 81, 86, 
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 72, 
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67, 
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}

df = pd.DataFrame(data=d)
x1 = df['Germ-like death rate, $d_g$'] 
y1 = df['Soma-like death rate, $d_s$'] 
s = df['size']
c = df['density']

ax.set_title(r'Fraction of germ-like,$f_g$',fontsize=14)
ax.scatter(x1, y1, s=s, c=color[c])
ax.set_xlabel('Germ-like death rate, $d_g$',fontsize=14)
ax.set_ylabel('Soma-like death rate, $d_s$',fontsize=14)

plt.show()
plt.savefig('D:/Users/DELL/Desktop/fenhua/workspace/script/jingsi/second_fig/fig4_3.pdf' ,bbox_inches='tight')   # save figures
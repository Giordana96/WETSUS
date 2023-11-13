# -*- coding: utf-8 -*-
"""
start with using pyplot to make graphs with pandas dataframe
although pandas also defines a plot command I want to write one myself to
get experience with OO programming 
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def make_test_dataframe(max):
    y=[]
    for i in range(max+1):
        y.concat([i,(2*i/max)**2,np.sin(i/max*2*np.pi),i/2.0])
    yd = pd.DataFrame(y,columns = ['i', 'isq','sin(i)','i/2'])     
    return yd

def plot_XY(title,data,x_name,y_name, size=(4,4), width=3.0, filesave = 0):
    #plots a number of y variables versus one x
    # x, y are indicated by the column name of the dataframe
    # title is a list with upper graph title, y label and unit of x lable if any
    # several options can be changed: thickness, size, filesave (1  is saving ouder title) ,  etc. 
    fig,ax = plt.subplots(figsize=size)
    x = data[x_name]        
    for y_n in y_name:
        y = data[y_n]
        ax.plot(x,y,linewidth=width)
    ax.set_title(title[0])
    if title[2]!='':
        x_name = x_name+'  ['+title[2]+']'
    ax.set_xlabel(x_name)
    ax.set_ylabel(title[1])
    ax.legend(loc = 'best')
    if filesave: fig.savefig(title[0]+'.tiff')
    return

def plot_all(title,data,x_name,width=3.0, filesave = 0, n_col = 2):
    #plot all cols in the datframe as funcrion of selected col
    # each col gets single graph
    n_axes = len(data.columns)-1
    n_row = (n_axes // n_col)
    if n_axes>n_row*n_col:
        n_row = n_row+1
    fig,axes = plt.subplots(nrows=n_row,ncols=n_col,sharex = True,figsize=(n_col*6,n_row*4))
#    print(n_axes,'   ', n_row,'   ',n_col)
    axes_list = data.columns.tolist()
    axes_list.remove(x_name)
#    print(axes_list)
    for i in range(n_axes):
#        row_i= (i// n_col)
#        col_i = i - row_i*n_col
#        print(i)
        plt.subplot(n_row,n_col,i+1)
        plt.plot(data[x_name],data[axes_list[i]],linewidth=width)       
        plt.xlabel(title[1])
#        print(axes_list[i])
        plt.legend((axes_list[i],), loc ='best')
    if filesave: fig.savefig(title[0]+'.tiff')    
    return

def main(choice):
    test = make_test_dataframe(100)
    if choice ==1:
        plot_XY(['single plot ','y','kg/m3'],test,'i',['sin(i)'])
    if choice==2:
        plot_XY(['double plot','yy',''],test,'i',['sin(i)','isq'], size=(7,3), filesave=1)
    if choice==3:
        plot_all(['all plot','yy'],test,'i')
 

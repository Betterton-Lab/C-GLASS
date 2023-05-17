# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:36:07 2022

@author: Daniel Steckhahn
"""
#Import
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from PIL import Image, ImageFilter, ImageEnhance

#Set red and green colors
bottom = cm.get_cmap('Reds', 128)
newcolors = np.vstack(
                       bottom(np.linspace(.7, .9, 128)))
redcmp = ListedColormap(newcolors)

bottom = cm.get_cmap('hsv', 128)
newcolors = np.vstack(
                       bottom(np.linspace(.3, .4, 128)))
greencmp = ListedColormap(newcolors)

#Plot a microtubule 
def plot_tube(y, start, stop, cmp=redcmp, points_per_length=200):
    start=start+40
    stop=stop+40
    length=int(stop-start)
    total_points=length*points_per_length
    x1 = np.linspace(start,stop,total_points)
    y1 = np.linspace(start,stop,total_points)
    x1=[random.gauss(0, 6) for i in x1]
    y1=[random.gauss(i, 6) for i in y1]
    plt.hist2d(x1, y1, bins=(1000, 1000), cmap=cmp,cmin = .9,range=np.array([(-30, 30), (-60,140)]))
    plt.gca().invert_yaxis()

#Plot a cut7
def plot_cut7(x_data,y_data, points_per_cut7=600,cmp=greencmp):
    x_points=[]
    y_points=[]

    for x in x_data:
        for i in range(points_per_cut7):
            x_points.append(random.gauss(x+80, 4))
    for y in y_data:
        for i in range(points_per_cut7):
            y_points.append(random.gauss(y, 4))            

    plt.hist2d(y_points, x_points, bins=(1000, 1000), cmap=cmp,cmin = .1,range=np.array([(-30, 30), (-60,140)]))                      
    plt.gca().invert_yaxis()

#Add background noise
def add_noise(points, cmp):    
    x1 = np.linspace(0,1, points)
    y1 = np.linspace(0,1, points)
    x1 = [-30+60*random.uniform(0, 1) for x in x1] 
    y1 = [-60+210*random.uniform(0, 1) for y in y1]
    plt.hist2d(x1, y1, bins=(1000, 1000), cmap=cmp,cmin = .9,range=np.array([(-30, 30), (-60,140)]))
    plt.gca().invert_yaxis()

#Ask for inputs

#Ask for simualtion name
base=input("What is the base name of the simulation: ")

#Ask for frames
frames=input("What frames would you like to use for images? Enter seperated by spaces: ")
frame_list = frames.split()
frame_array=[int(i) for i in frame_list]

#Ask for microtubules in protrusion
pro_microtubules=input("How many microtubules are in protrusion?: ")
pro_microtubules=int(pro_microtubules)

#Cycle through frames that images are being created for
for f in frame_array:
    print('Making flourescent image for frame', f)
    fig, ax = plt.subplots()

    #Set background color
    color_num=.14
    ax.set_facecolor((color_num, color_num, color_num))
    fig.patch.set_facecolor((color_num, color_num, color_num))

    #Adjust shape of plot
    plt.subplots_adjust(bottom=0, top=.9, left=0, right=.1)

    #Turn off ticks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off

    #Turn off axis border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    #add noise
    add_noise(220000, redcmp)
    add_noise(80000, greencmp)

    print('Warning: microtubules set for locations in protrusion simulations, if locations change microtubule fluorescence wont be accurate')

    if pro_microtubules>3:
        print("Warning: fluorescent image only set up for 3 or less protrusion microtubules")

    #plot the microtubules
    frame=f
    loc=-(frame-505)/572*80
    #1
    plot_tube(0, loc,loc+80)
    if pro_microtubules==1:	
       loc=0
    #2
    plot_tube(0, loc,loc+80)    
    if pro_microtubules==2:	
       loc=0
    #3
    plot_tube(0, loc,loc+80)
    #4
    plot_tube(.6, 0,35)
    #5
    plot_tube(-.6, 0,20)
    #6
    plot_tube(-.6, 0,10)
    #7
    plot_tube(-1.2, 0,22)
    #8
    plot_tube(1.2, 0,50)
    #9
    plot_tube(0, 0,60)
    #10
    plot_tube(0, 0,30)
    #11
    plot_tube(0, 0,70)
    #12
    plot_tube(0, 0,25)    
    #13
    plot_tube(.6, 0,65)
    #14
    plot_tube(.6, 70,80)
    #15
    plot_tube(-.6, 53,80)
    #16
    plot_tube(-.6, 30,80)
    #17
    plot_tube(-.6, 75,80)
    #18
    plot_tube(-.6, 20,80)
    #19
    plot_tube(0,45,80)
    #20
    plot_tube(0, 65,80)
    #21
    plot_tube(0, 55,80)   
    #22
    plot_tube(0, 48,80)    
    #23
    plot_tube(-.6, 10,80)
    #24
    plot_tube(0, 60,80)
    #25
    plot_tube(-.6, 15,80)
    #26
    plot_tube(0, 50,80)

    #Load and plot cut7 locations
    stringx='./data/'+base+'_x_'+str(frame)+'.npy'
    stringy='./data/'+base+'_y_'+str(frame)+'.npy'
    try:
        x=np.load(stringx)
    except FileNotFoundError:
        print(base+'_x_'+str(frame)+'.npy'+' not in the data folder')
    try:
        y=np.load(stringy)
    except FileNotFoundError:
        print(base+'_y_'+str(frame)+'.npy'+' not in the data folder')
    plot_cut7(x,y)

    #Saving as an png, then opening it as a png instead of a plot for image processing
    unprocessed_image='unprocessed_image.png'
    plt.savefig(unprocessed_image,bbox_inches='tight')
    OriImage = Image.open(unprocessed_image)

    #Blur the image
    blurImage = OriImage.filter(ImageFilter.GaussianBlur(3.5))
    blurImage = blurImage.filter(ImageFilter.GaussianBlur(1.8))
    enhancer = ImageEnhance.Brightness(blurImage)
    
    #Make the image brighter
    enhancer = ImageEnhance.Brightness(blurImage)
    factor = 1.9
    brighter = enhancer.enhance(factor)

    #Save final image
    processed_image='./data/processed_'+base+'_'+str(frame)+'.png'
    brighter.save(processed_image)
    print('Fluorescent image made for frame', str(frame), '. Image is in file "data".')

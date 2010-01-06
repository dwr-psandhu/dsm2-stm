import pylab as pl
import numpy as np
import sys

for arg in sys.argv:
    print arg

#x = np.linspace(xmin,xmax,n)

xmin = 0
xmax = 100.
ymin = 0
ymax = 6.
nx = 0
nplot = 0 

datafile      = open('../../../stm/build/vs2008/example_single_channel/temp.dat', 'r')
datafile_ref  = open('../../../stm/build/vs2008/example_single_channel/temp_ref.dat', 'r')
datalines     = datafile.readlines()
datalines_ref = datafile_ref.readlines()
dataraw=[]
data=[]
i=0

# todo: rewrite below code into a function

while i< len(datalines):

    if "zone" in datalines[i]:
        nx = int(datalines[i].split()[3])
        nplot += 1 
        
        # read x y 
        dataraw.append(line.split() for line in datalines[i+1:i+nx+1])
        
        # fast forward
        i += nx

    else:
        i += 1
        
x=[]
y=[]

for set in dataraw:
    # transpose to [(x1, x2,... xn), (y1, y2, ... yn)]
    datatranspose =  zip(*set)
    
    x.append(np.array(map(float,datatranspose[0]) ))
    y.append(np.array(map(float,datatranspose[1]) ))
    



v = [xmin,xmax,ymin,ymax]

def plot():

    pl.figure(1)
    pl.xlabel('x (m)')
    pl.ylabel('concentration')    
    pl.plot(x[0],y[0],'ro',x[0],y[0],'g')
    pl.axis(v)

    pl.savefig("test1.pdf", dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='pdf',
        transparent=False)

    pl.figure(2)
    pl.xlabel('x (m)')
    pl.ylabel('concentration')
    pl.plot(x[1],y[1],'ro',x[1],y[1],'g')
    pl.axis(v)

    pl.savefig("test2.pdf", dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='pdf',
        transparent=False)

plot()
pl.show()


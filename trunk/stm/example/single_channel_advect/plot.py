import pylab as pl
import numpy as np
import sys
import re

# todo: read these file names from a input file or system arguments.
datafile      = open('../../../stm/build/vs2008/example_single_channel/uniform_rectangular.dat', 'r')
datafile_ref  = open('../../../stm/build/vs2008/example_single_channel/uniform_rectangular_ref.dat', 'r')
#for arg in sys.argv:
#    print arg

xmin = 0
xmax = 100.
ymin = 0
ymax = 6.
v = [xmin,xmax,ymin,ymax] #axis range for plotting
nx = 0
nplot = 0 
x=[]
y=[]

def parse2xy(infile):
    
    lines = infile.readlines()
    i=0
    nzone = 0
    data = []
    dataTransposed =[]
    xc=[]
    yc=[]
    while i< len(lines):
    
        if "zone" in lines[i]:
    
            p = re.compile('i.*\= *(\d+)[, ].+\"(.+)\"')            
            ps = p.search(lines[i])
            
            nsize = int(ps.group(1))
            nzone += 1 
            
            # read x y
            data=[]
            
            # list of [['x0','y0'],['x1','y1'],...]
            for line in lines[i+1:i+nsize+1]:
            
                data.append(line.split())
                        
            # transpose to [('x1','x2',...),('y1','y2',...)]
            dataTransposed = zip(*data)
            
            xc.append(np.array(map(float,dataTransposed[0])))
            yc.append(np.array(map(float,dataTransposed[1])))
            
            # fast forward to the end of this zone
            i += nsize
    
        else:
            i += 1

    return xc,yc,nsize,nzone





def plot(figureID,x,y,xr,yr,axisRange,outFile):

    pl.figure(figureID)
    pl.xlabel('x (m)')
    pl.ylabel('concentration')    
    pl.plot(x,y,'ro',xr,yr,'g')
    pl.axis(axisRange)

    pl.savefig(outFile+".pdf", dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='pdf',
        transparent=False)



x,y,nx,nplot = parse2xy(datafile)
xr,yr,nx,nplot = parse2xy(datafile_ref)

plot(1,x[3],y[3],xr[3],yr[3],v,'3')
plot(2,x[6],y[6],xr[6],yr[6],v,'6')
pl.show()


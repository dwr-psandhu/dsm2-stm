import pylab as pl
import numpy as np


#x = np.linspace(xmin,xmax,n)

datafile= open('D:/delta/trunk/stm/build/vs2008/report/uniform_rectangular_at_itime_  3.txt', 'r')

data = [line.split() for line in datafile]

# transpose to [(x1, x2,... xn), (y1, y2, ... yn)]
data = zip(*data)

x = np.array(map(float, data[0]))
y = np.array(map(float, data[1]))


def plot_fit():

    pl.xlabel('x (m)')
    pl.ylabel('concentration')
    pl.plot(x,y,'ro')

plot_fit()
pl.show()
import pylab as pl
import numpy as np

#x = np.linspace(xmin,xmax,n)

datafile   = open('D:/delta/trunk/stm/build/vs2008/example_single_channel/uniform_rectangular_at_itime_ 20.txt', 'r')
datafile0  = open('D:/delta/trunk/stm/build/vs2008/example_single_channel/uniform_rectangular_at_itime_0.txt', 'r')

data  = [line.split() for line in datafile ]
data0 = [line.split() for line in datafile0]

# transpose to [(x1, x2,... xn), (y1, y2, ... yn)]
data = zip(*data)
data0 = zip(*data0)

x = np.array(map(float, data[0]))
y = np.array(map(float, data[1]))

x0 = np.array(map(float, data0[0]))
y0 = np.array(map(float, data0[1]))

n=10000

x_exact = np.linspace(float(x[0]), float(x[-1]),n)
y_exact = np.zeros(n).tolist()

for i in range(n):
    if (x_exact[i] >= 10. and x_exact[i] <= 40.):
        y_exact[i] = 5.

y_exact =np.array(y_exact)

def plot():

    pl.figure(1)
    pl.xlabel('x (m)')
    pl.ylabel('concentration')    
    pl.plot(x0,y0,'ro',x_exact,y_exact,'g')
    #pl.axis('scaled')

    pl.savefig("test.pdf", dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='pdf',
        transparent=False)

    pl.figure(2)
    pl.xlabel('x (m)')
    pl.ylabel('concentration')
    pl.plot(x,y,'ro',x_exact,y_exact,'g')


plot()
pl.show()


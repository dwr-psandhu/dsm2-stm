import pylab as pl
import numpy as np

#x = np.linspace(xmin,xmax,n)

datafile= open('D:/delta/trunk/stm/build/vs2008/report/uniform_rectangular_at_itime_  3.txt', 'r')


data = [line.split() for line in datafile]

# transpose to [(x1, x2,... xn), (y1, y2, ... yn)]
data = zip(*data)

x = np.array(map(float, data[0]))
y = np.array(map(float, data[1]))

n=10000

x2 = np.linspace(float(x[0]), float(x[-1]),n)
y2 = np.zeros(n)
y2.tolist()


for i in range (n):
    if (x2[i] >= 20. and x2[i] <= 70.):
        y2[i] = 5.

y2 =np.array(y2)

def plot_fit():

    pl.xlabel('x (m)')
    pl.ylabel('concentration')
    pl.plot(x,y,'ro',x2,y2,'g')

plot_fit()
pl.show()
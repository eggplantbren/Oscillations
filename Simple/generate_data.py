from pylab import *

seed(123)

t = linspace(0., 100., 1001)
y = sin(2.*pi*t/30.) + sin(2.*pi*t/2. + 1.)
Y = y + randn(y.size)

data = empty((t.size, 2))
data[:,0], data[:,1] = t, Y

savetxt('data.txt', data)

plot(t, y, 'r')
plot(t, Y, 'bo')
show()


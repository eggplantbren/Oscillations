import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

rng.seed(0)

x_min, x_max = 0., 2000.
num_points = 100001
x = np.linspace(0., 2000., num_points)
y = np.zeros(num_points)

num = 100
for i in xrange(0, num):
	# Peak, width, center
	A  = -np.log(rng.rand())
	w  = 3*rng.rand()
	xc = x_min + (x_max - x_min)*rng.rand()

	y += A/(1. + ((x - xc)/w)**2)

y += 0.2 # Constant noise level
y *= -np.log(rng.rand(y.size))

data = np.empty((y.size, 2))
data[:,0], data[:,1] = x, y
np.savetxt('simulated_data.txt', data)

plt.plot(x, y)
plt.show()


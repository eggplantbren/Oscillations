import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

samples = np.atleast_2d(np.loadtxt('posterior_sample.txt'))
data = np.loadtxt('simulated_data.txt')
maxNumComponents = (samples.shape[1] - 3)//3

for i in xrange(0, samples.shape[0]):
	model = samples[i, 3:]

	A = model[0:maxNumComponents]
	f = model[maxNumComponents:2*maxNumComponents]
	w = model[2*maxNumComponents:]
	which = A > 0.
	A, f, w = A[which], f[which], w[which]

	y = np.zeros(data.shape[0]) + 0.2
	for j in xrange(0, A.size):
		y += A[j]/(1. + ((data[:,0] - f[j])/w[j])**2)

	plt.figure(figsize=(12, 6))
	plt.clf()
	plt.hold(True)
	plt.plot(data[:,0], data[:,1], 'b')
	plt.plot(data[:,0], y, 'r')
	plt.show()


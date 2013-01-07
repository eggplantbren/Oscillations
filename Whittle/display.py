import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

samples = np.atleast_2d(np.loadtxt('posterior_sample.txt'))
data = np.loadtxt('simulated_data.txt')
maxNumComponents = (samples.shape[1] - 3)//3

plt.ion()
plt.figure(figsize=(12, 6))

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

	plt.hold(False)
	plt.plot(data[:,0], data[:,1], 'b', alpha=0.3, label='Data')
	plt.hold(True)
	plt.plot(data[:,0], y, 'r', label='Model')
	plt.xlabel('Frequency ($\\mu$Hz)', fontsize=16)
	plt.ylabel('Power', fontsize=16)
	plt.title('Model ' + str(i+1))
	plt.legend()
	plt.draw()

plt.ioff()
plt.show()

#	plt.figure(2)
#	plt.hist(1. - np.exp(-data[:,1]/y), 100)
#	plt.show()


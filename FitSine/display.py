# Copyright (c) 2009, 2010, 2011, 2012 Brendon J. Brewer.
#
# This file is part of DNest3.
#
# DNest3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DNest3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DNest3. If not, see <http://www.gnu.org/licenses/>.

"""
Load models and display them
"""

import numpy as np
import matplotlib.pyplot as plt

maxNumComponents = 20
data = np.atleast_2d(np.loadtxt('fake_data.txt'))
sample = np.atleast_2d(np.loadtxt('posterior_sample.txt'))

start = 0

numComponents = sample[:, start].astype('int')
start += 1

muAmplitudes = sample[:, start]
start += 1

sigmaBoost = sample[:, start]
start += 1

dof = sample[:, start]
start += 1

staleness = sample[:, start]
start += 1

amplitudes = sample[:, start:start+maxNumComponents]
start += maxNumComponents

frequencies = sample[:, start:start+maxNumComponents]
start += maxNumComponents

phases = sample[:, start:start+maxNumComponents]
start += maxNumComponents

plt.ion()
for i in xrange(0, sample.shape[0]):
	
	mockData = np.zeros(data.shape[0])
	for j in xrange(0, numComponents[i]):
		mockData += amplitudes[i, j]*np.sin(2*np.pi*frequencies[i, j]*data[:,0] + phases[i, j])

	plt.subplot(2,1,1)
	plt.hold(False)
	plt.plot(data[:,0], data[:,1], 'b.')
	plt.hold(True)
	plt.plot(data[:,0], mockData, 'r')
	plt.axis([-1., 101., -15., 15.])
	plt.xlabel('Time')
	plt.ylabel('y')
	plt.title('Model %i, %i components.'%(i+1, numComponents[i]))

	plt.subplot(2,1,2)
	plt.hold(False)
	num = numComponents[0:(i+1)]
	plt.hist(num, bins=np.arange(0, num.max()+2), align='left', rwidth=0.3)
	plt.xlim([-0.5, num.max() + 0.5])
	plt.xlabel('Number of Components')
	plt.ylabel('Probability')
	plt.draw()

plt.ioff()
plt.show()


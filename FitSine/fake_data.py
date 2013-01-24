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
generate fake data
"""

import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt

if __name__ == '__main__':
	rng.seed(0)

	t = np.sort(100.*rng.randn(100))
	y = np.zeros(t.shape)

	# Number of components
	N = 10

	# Frequencies
	f = np.exp(np.log(0.1) + rng.randn(N))

	# Amplitudes
	A = -np.log(rng.rand(N))

	# Phases
	phi = 2*np.pi*rng.rand(N)
	
	# Add components to the fake data
	for i in xrange(0, N):
		y += A[i]*np.sin(2*np.pi*f[i]*t + phi[i])

	# Errorbars
	sig = 0.3*np.ones(t.size)

	# Add noise
	y += sig*rng.randn(t.size)

	np.savetxt('fake_data.txt', np.vstack([t, y, sig]).T)
	plt.errorbar(t, y, yerr=sig, fmt='b.')
	plt.show()


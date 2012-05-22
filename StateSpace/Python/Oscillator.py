import numpy as np
import numpy.random as rng

class Oscillator:
	"""
	A point in phase space for an oscillator.
	"""
	def __init__(self, state, omega=1., tau=1., beta=1.):
		"""
		Constructor: takes initial yition
		and velocity as argument. Sets the time to zero
		"""
		self.state = state
		self.omega, self.tau, self.beta = omega, tau, beta
		self.time = 0.

	def deriv(self, time, state, dt):
		"""
		Compute the derivatives from the given state
		(not necessarily self.state, you need to pass that in
		if that's what you want!)
		"""
		a = np.empty(2)
		a[0] = state[1]
		a[1] = -self.omega**2*state[0] - state[1]/self.tau\
				+ self.beta*rng.randn()/np.sqrt(dt)
		return a

	def update(self, dt):
		"""
		Take a step using RK4
		"""
		f1 = self.deriv(self.time, self.state, dt)
		f2 = self.deriv(self.time + 0.5*dt, self.state + 0.5*dt*f1, dt)
		f3 = self.deriv(self.time + 0.5*dt, self.state + 0.5*dt*f2, dt)
		f4 = self.deriv(self.time + dt, self.state + dt*f3, dt)
		self.state += dt/6.*(f1 + 2*f2 + 2*f3 + f4)
		self.time += dt


if __name__ == '__main__':
	import matplotlib.pyplot as plt

	# Initial conditions
	oscillator = Oscillator(np.array([0., 0.]))

	# Timestep
	dt = 0.01

	steps = 100000   # Take this many steps
	skip = 100	 # Store and plot results this often
	keep = np.empty((steps/skip, 3)) # Store results in here
					 # Columns: time, yition, vocity

	plt.ion() # Turn "interactive mode" for plotting on, so plots
		  # can update without the user having to close the window
	plt.hold(False) # Clear the plot every time we plot something new
	# Main loop
	for i in xrange(0, steps):
		# Saving and plotting
		if i%skip == 0:
			index = i/skip
			# Save state to keep array
			keep[index, :] = \
				np.array([oscillator.time, oscillator.state[0],\
				oscillator.state[1]])
			# Plot yition vs time
			plt.plot(keep[0:(index+1), 0], keep[0:(index+1), 1], 'b')
			plt.xlabel('Time')
			plt.ylabel('y')
			plt.title('Stdev = %.03f'%keep[0:(index+1), 1].std())
			plt.draw() # Refresh the plot

		# Update the oscillator
		oscillator.update(dt)

	# At end of run, leave the last plot showing until the user closes it
	plt.ioff()
	plt.show()


#Iterate chisq_chisq.py over a large set of different variables and return all the output in csv files
#Finds SNR values and Chisq values for matched filtering at different chirp masses and eccentricities
import time

start = time.time()

import pycbc.noise
import pycbc.psd
import pycbc.filter
import pycbc.waveform
import pycbc.vetoes
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

#functions to generate mass1 and mass2 from chirp mass and symmetric mass ratio

def func_mass1(chirp,smr):
	x = chirp/(smr**(0.6))
	y = (chirp**2)/(smr**(0.2))
	if smr ==0.25:
		return x/2
	else:
		mass1 = (x/2) + (((x**2)/4)-y)**0.5
		return(mass1)
	
def func_mass2(chirp,smr):
	x = chirp/(smr**(0.6))
	y = (chirp**2)/(smr**(0.2))
	if smr ==0.25:
		return x/2
	else:
		mass2 = y/((x/2) + (((x**2)/4)-y)**0.5)
		return(mass2)

# Generate some noise with an advanced ligo psd
flow = 15.0
delta_f = 1.0 / 16
flen = int(2048 / delta_f) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

#Number of variable values to iterate over
nchirp = 50
neccentricity = 50
teststeps = nchirp * neccentricity

#Arrays to store the output
chirps = np.zeros(teststeps)
eccentricities = np.zeros(teststeps)
smrs = np.zeros(teststeps)
chisqs = np.zeros(teststeps)
counter = 0

for a in range(nchirp):
	for b in range(neccentricity):
		#input parameters to vary
		chirp = 1+0.28*a
		eccentricity = 0.006*b
		smr = 0.2
		snr_ideal = 12 #the SNR that the input would be if matched perfectly

		mass1 = func_mass1(chirp,smr)
		mass2 = func_mass2(chirp,smr)


		#Find snr and fix it to a chosen value by changing the distance parameter
		max_snr = 0.0
		distlow = 0.0
		disthigh = 0.0
		distance = 500.0
		threshold = 0.01
		# Use a circular waveform as a matched filter, only generate this once
		Tp, Tc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
									mass1=mass1, mass2=mass2, eccentricity = 0.0,
									f_lower=flow, delta_f=delta_f)

		#find the distance which makes the snr very close to snr_ideal
		#This is slow but since I am no longer using this function it is not worth changing
		while(abs(max_snr-snr_ideal) > threshold):
	
			# Generate an eccentric waveform as a signal
			Sp, Sc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
										mass1=mass1, mass2=mass2, eccentricity = eccentricity, distance = distance, 
										f_lower=flow, delta_f=delta_f)
	
			Tp.resize(len(psd))
			Sp.resize(len(psd))
			snr = pycbc.filter.matched_filter(Tp, Sp, psd=psd,
                                      low_frequency_cutoff=flow)

			#find the max snr and it's index
			max_snr = max(abs(snr))
			max_index = np.argmax(abs(snr))
			if abs(max_snr-snr_ideal) <threshold:
				break
			if max_snr > snr_ideal:
				distlow = distance
				if disthigh ==0:
					distance = distance*2
				else:
					distance = (disthigh+distlow)/2
			else:
				disthigh = distance
				distance = (disthigh+distlow)/2

		#Find the chisq value for this signal
		num_bins = 16
		Tp.resize(len(psd))
		Sp.resize(len(psd))

		chisq = pycbc.vetoes.power_chisq(Tp, Sp, num_bins, psd=psd,
											low_frequency_cutoff=flow)

		#create a reduced chisq value
		chisq /= (num_bins * 2) - 2
		chisq_value = chisq[max_index]
		
		#Print the values to show progress
		
		print str(chirp) + ' ' + str(eccentricity) + ' ' +  str(max_snr) + ' ' + str(chisq_value)
		chirps[counter] = chirp
		eccentricities[counter] = eccentricity
		smrs[counter] = smr
		chisqs[counter] = chisq_value
		counter+=1
np.savetxt('chirps.csv',chirps,delimiter =',')
np.savetxt('eccentricities.csv',eccentricities,delimiter =',')
np.savetxt('smrs.csv',smrs,delimiter =',')
np.savetxt('chisqs.csv',chisqs,delimiter =',')


end = time.time()

print 'time elapsed is ' + str(end - start) + 'seconds'

#Finds the modified SNR for an eccentric signal by increasing the power of the n=2 mode by a theoretical factor
#This has very little effect, it is the phase which causes the greater discrepancy (see phase.py)
import time

start = time.time()

import pycbc.noise
import pycbc.psd
import pycbc.filter
import pycbc.waveform
import pycbc.vetoes
import scipy.special as sp
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
#Physical Constants
G = 1.3275E+20 #in (m)^3 (Solar mass)^-1 (s)^-2

#functions to find mass1 and mass2 from chirp mass and smr

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

#generate a non-eccentric template
chirp = 2
smr = 0.25
mass1 = func_mass1(chirp,smr)
mass2 = func_mass2(chirp,smr)
m= mass1+mass2
mu = (mass1*mass2)/m
distance = 400
e_init = 0.2
Tp, Tc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
					mass1=mass1, mass2=mass2, eccentricity = 0.0, distance = distance,
					f_lower=flow, delta_f=delta_f)
Tp.resize(len(psd))					


#Find the positions of the bins used in the chisq test

num_bins = 16
bins = pycbc.vetoes.chisq.power_chisq_bins(Tp,num_bins = num_bins,psd = psd, low_frequency_cutoff = flow, high_frequency_cutoff = 2000.0)

#Find the positions of the centres of the bins

bincentres = np.zeros(num_bins)
for i in range(len(bincentres)):
	bincentres[i] = (bins[i]+bins[i+1])/2

#Find the eccentricity at the frequencies specified in bincentres

def gfunc(e):
	return ((e**(12.0/19.0))/(1-e**2))*((1+(121.0/304.0)*(e**2))**(870.0/2299.0))

def gdiff(e,a):
	num = 12.0*a*(1.0+((73.0/24.0)*e**2)+((37.0/96.0)*e**4))
	den = 19.0*e*(1-(e**2))*(1+((121.0/304.0)*e**2))
	return num/den

#Use Newton's method to find eccentricity for a given frequency(via a, the semi major axis)

flow_orbit = flow/2.0 #since lowest GW frequency is double orbital frequency
a_init = (G*m/((2*np.pi*flow_orbit)**2))**(1.0/3.0)
if e_init > 0: #To avoid divide by 0 errors
	c0 = a_init / gfunc(e_init) #This is a factor which fixes the initial values of e and a
eccs = np.zeros(num_bins)
for i in range(num_bins):
	j = bincentres[i]
	w = (delta_f*j)
	a_target = (G*m/((2*np.pi*w/2)**2))**(1.0/3.0)
	n=1
	eguess = e_init
	if e_init > 0:
		while n < 10: #Perform Newton's method for 10 iterations
			n+=1
			a = c0*gfunc(eguess)
			eguess = eguess - (a-a_target)/gdiff(eguess,a)
			if eguess<0:
				eguess = -eguess #this prevents the method breaking when e becomes negative,
	eccs[i] = eguess
	#print str(w) + '    ' + str(eguess) + '   ' + str(a_target) + '    ' + str(a)

#For each frequency in bincentres, find the corresponding power in n=2 band relative to e=0 case

powerfactors = np.zeros(num_bins)

for i in range(num_bins):
	orbitfreq = bincentres[i]/2.0
	e = eccs[i]
	a = (G*m/((2*np.pi*orbitfreq)**2))**(1.0/3.0)
	b = a*((1.0-e**2)**0.5)
	A2 = ((a**2)/2.0)*(sp.jv(0,2.0*e)-sp.jv(4,2.0*e)-2.0*e*(sp.jv(1,2.0*e)-sp.jv(3,2.0*e)))
	B2 = ((b**2)/2.0)*(sp.jv(4,2.0*e)-sp.jv(0,2.0*e))
	C2 = ((a*b)/2.0)*(sp.jv(0,2.0*e)+sp.jv(4,2.0*e)-e*(sp.jv(1,2.0*e)+sp.jv(3,2.0*e)))
	g = (A2**2 + B2**2 + 3*(C2**2) - A2*B2)
	power_e0 = 1.5*(a**4)
	power = g
	powerfactor = power/power_e0
	powerfactors[i] = powerfactor
	#print str(e) + '    ' + str(powerfactor)
	
	
#Find the new SNR by breaking the source and template into sections and weighting by the factors which have been calculated

#Generate a source with the correct eccentric form

Sp, Sc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
					mass1=mass1, mass2=mass2, eccentricity = e_init, distance = distance,
					f_lower=flow, delta_f=delta_f)
Sp.resize(len(psd))

#Normalise the source such that SNR=12 if matched perfectly with itself

snr = pycbc.filter.matched_filter(Sp, Sp, psd=psd,
                                      low_frequency_cutoff=flow, high_frequency_cutoff=2000.0)
max_snr = max(abs(snr))

distance = distance * max_snr/12.0


#Redefine Sp and Tp at this distance

Tp, Tc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
					mass1=mass1, mass2=mass2, eccentricity = 0.0, distance = distance,
					f_lower=flow, delta_f=delta_f)
Tp.resize(len(psd))		

Sp, Sc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
					mass1=mass1, mass2=mass2, eccentricity = e_init, distance = distance,
					f_lower=flow, delta_f=delta_f)
Sp.resize(len(psd))

#TEST

snr = pycbc.filter.matched_filter(Tp, Sp, psd=psd,
                                      low_frequency_cutoff=flow, high_frequency_cutoff = 2000.0)
print max(abs(snr))
#END TEST


snrs = np.zeros(num_bins)
for i in range(num_bins):
	low_frequency_cutoff = (delta_f*bins[i])
	high_frequency_cutoff = (delta_f*bins[i+1])
	snr = pycbc.filter.matched_filter(Tp, Sp, psd=psd,
                                      low_frequency_cutoff=low_frequency_cutoff, 
									  high_frequency_cutoff = high_frequency_cutoff,)
	snrs[i] = (max(abs(snr))/(num_bins**0.5))/powerfactors[i] #num_bins**0.5 is needed to weight the output so total SNR is independent of num_bins
	print str(low_frequency_cutoff)+'   '+str(high_frequency_cutoff)+'    '+str(snrs[i])
	
total = sum(snrs)
print total

#Generate a chi^2 value using these values for the SNR

chisq = 0
for i in range(num_bins):
	chisq += num_bins*((total/num_bins)-snrs[i])**2

reduced_chisq = chisq/(2.0*num_bins-2)
print reduced_chisq
end = time.time()

print 'time elapsed is ' + str(end - start) + 'seconds'
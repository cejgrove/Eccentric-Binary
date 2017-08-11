#Investigate possible uses of SNR time series to find eccentricity
import pycbc.noise
import pycbc.psd
import pycbc.filter
import pycbc.waveform
import matplotlib
#matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
# Generate some noise with an advanced ligo psd
flow = 15.0
delta_f = 1.0 / 16
flen = int(2048 / delta_f) + 1
psd = pycbc.psd.aLIGOZeroDetHighPower(flen, delta_f, flow)

#Initiate the signals with different eccentricities

Tp, Tc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
                             mass1=5, mass2=5, eccentricity = 0.0,
                             f_lower=flow, delta_f=delta_f, distance = 400)

Sp, Sc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
                             mass1=5, mass2=5, eccentricity = 0.1,
                             f_lower=flow, delta_f=delta_f, distance = 400)

Rp, Rc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
                             mass1=5, mass2=5, eccentricity = 0.2,
                             f_lower=flow, delta_f=delta_f, distance = 400)
Qp, Qc = pycbc.waveform.get_fd_waveform(approximant="EccentricFD",
                             mass1=5, mass2=5, eccentricity = 0.3,
                             f_lower=flow, delta_f=delta_f, distance = 400)

Tp.resize(len(psd))									  
Sp.resize(len(psd))
Rp.resize(len(psd))
Qp.resize(len(psd))

#Create and add noise to the signals
noise = pycbc.noise.frequency_noise_from_psd(psd, seed=127)
Tpnew = Tp #+ noise	#remove these comments to add the noise						  
Spnew = Sp #+ noise
Rpnew = Rp #+ noise
Qpnew = Qp #+ noise

#Create the SNR time series	
	
snr0 = pycbc.filter.matched_filter(Tp, Tpnew, psd=psd,
                                      low_frequency_cutoff=flow)
snr1 = pycbc.filter.matched_filter(Tp, Spnew, psd=psd,
                                      low_frequency_cutoff=flow)
snr2 = pycbc.filter.matched_filter(Tp, Rpnew, psd=psd,
                                      low_frequency_cutoff=flow)
snr3 = pycbc.filter.matched_filter(Tp, Qpnew, psd=psd,
                                      low_frequency_cutoff=flow)



									  



plt.plot(snr0.sample_times[0:800], abs(snr0)[0:800],label = 'e = 0.0')
plt.plot(snr1.sample_times[0:800], abs(snr1)[0:800],label = 'e = 0.1')
plt.plot(snr2.sample_times[0:800], abs(snr2)[0:800],label = 'e = 0.2')
plt.plot(snr3.sample_times[0:800], abs(snr3)[0:800],label = 'e = 0.3')
plt.legend()
plt.ylabel('signal-to-noise ratio')
plt.xlabel('time (s)')
#plt.savefig("snrtime.png")


#np.savetxt('times.csv',snr0.sample_times,delimiter = ',')
#np.savetxt('snr0.csv',abs(snr0),delimiter =',')
#np.savetxt('snr1.csv',abs(snr1),delimiter =',')
#np.savetxt('snr2.csv',abs(snr2),delimiter =',')
#np.savetxt('snr3.csv',abs(snr3),delimiter =',')



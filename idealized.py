'''
Idealized experiments.
'''

import numpy as np
import matplotlib.pyplot as plt


# Define the base idealized functions
# stuff for the common first part of the spectra
f = [1e-5, 1e-4]; fs = np.logspace(np.log10(f[0]), np.log10(f[1])) # x endpoints
uf = [2e-2, 3e-4] # y endpoints
m = (np.log10(uf[0])-np.log10(uf[1]))/(np.log10(f[0])-np.log10(f[1])) # slope
c = uf[1]/(f[1]**m) 
ufs = c*fs**m

# stuff for the separate lines
f = [1e-4, 1e-3]; fss = np.logspace(np.log10(f[0]), np.log10(f[1])) # x endpoints
uf = [3e-4, 1e-4] # y endpoints
m = (np.log10(uf[0])-np.log10(uf[1]))/(np.log10(f[0])-np.log10(f[1])) # slope
c = uf[1]/(f[1]**m) 
ufshigh = c*fss**m

uf = [3e-4, 4e-5] # y endpoints
m = (np.log10(uf[0])-np.log10(uf[1]))/(np.log10(f[0])-np.log10(f[1])) # slope
c = uf[1]/(f[1]**m) 
ufsmed = c*fss**m

uf = [3e-4, 2e-5] # y endpoints
m = (np.log10(uf[0])-np.log10(uf[1]))/(np.log10(f[0])-np.log10(f[1])) # slope
c = uf[1]/(f[1]**m) 
ufslow = c*fss**m

# combine lines
f = np.concatenate((fs, fss))
ufh = np.concatenate((ufs, ufshigh))
ufm = np.concatenate((ufs, ufsmed))
ufl = np.concatenate((ufs, ufslow))

# add noise
noise = np.random.uniform(-.1, .1, f.size)
ufh = ufh + ufh*noise
ufm = ufm + ufm*noise
ufl = ufl + ufl*noise

# plot the lines
plt.figure()
plt.loglog(f, ufh, 'k-', f, ufm, 'r-', f, ufl, 'b-', lw=3)
# plt.loglog(fs, ufs, '-', color='0.5', lw=3)
# plt.loglog(fss, ufshigh, 'k-', fss, ufsmed, 'r-', fss, ufslow, 'b-', lw=3)
plt.legend(('high', 'med', 'low'))
plt.xlabel('frequency')
plt.ylabel('|u(f)|')
plt.show()

# Convert from frequency to time space, using code from time_res.ipynb
Ts = 60 # sampling interval: 60 seconds
Fs = 1./Ts # sampling rate
n = f.size # length of signal
k = np.arange(n)
T = n/Fs
frq = k/T # two sides frequency range
frq = frq[range(n/2)] # one side frequency range

# double frequency spectrum
Y2 = np.concatenate((Y[::-1],Y))

y = np.fft.ifft(ufh)
plt.figure()
plt.plot(y.real)
plt.plot(y.imag)
plt.show()




## Test case - do a tidal test or similar ##
n = np.zeros((400,))
ufh =  # a spike of 10 at 0.5 day frequency
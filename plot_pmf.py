#plot that sucker
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

print "I look for pmf.dat and plot it."

plt.figure(figsize=(8,6), dpi=90)
data = np.genfromtxt('pmf.dat')
#entropy correction
ent_corr = 2 * 2.494 * np.log(data[:,0] / data[np.argmin(data[:,1]),0])
data[:,1] += ent_corr
#get 0-line
# 0.5 < x < 1
zero_point = np.mean(data[np.where(data[:,0] > 0.4)[0],1][np.where(data[np.where(data[:,0] > 0.4)[0], 0] < 1)[0]])
data[:,1] -= zero_point
#plot
plt.plot(data[:,0], data[:,1])
ax = plt.gca()
ax.set_xlim([0.2,1.0])
ymax = ax.get_ylim()[0]
ax.set_ylim(ymax, -ymax)
plt.xlabel('r [nm]')
plt.ylabel(r'$\Delta$ A [kJ/mol]')
plt.savefig('pmf.png')
np.savetxt('pmf_processed.txt', data)

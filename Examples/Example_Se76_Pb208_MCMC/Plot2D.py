import numpy as np
import matplotlib.pyplot as plt
import sys

x = []
y = []

print sys.argv[1]

par1 = int(sys.argv[1])
par2 = int(sys.argv[2])

argcounter = 0
for filename in sys.argv:
	if argcounter > 2:
		with open(filename) as f:
			counter = 0
			for line in f:	
				if counter > 5000:
					x.append(float(line.split()[par1]))
					y.append(float(line.split()[par2]))
				counter += 1
	argcounter += 1

x = np.asarray(x)
y = np.asarray(y)
		
print x, y
	
fig = plt.figure()
plt.subplots_adjust(wspace=0,hspace=0)
axs = []
axs.append(plt.subplot2grid((6,6),(0,2),colspan=4,rowspan=4))
axs.append(plt.subplot2grid((6,6),(4,2),colspan=4,rowspan=2))
axs.append(plt.subplot2grid((6,6),(0,0),colspan=2,rowspan=4))

axs[0].hist2d(y,x,bins=(40,40))
axs[0].get_xaxis().set_ticks([])
axs[0].get_yaxis().set_ticks([])
axs[1].hist(y,bins=(40))
axs[1].yaxis.tick_right()
axs[2].hist(x,bins=(40),orientation='horizontal')
axs[1].set_xlabel(r"$<{2^+_1}|E2|{2^+_1}>$ [eb]")
axs[2].set_ylabel(r"$<{0^+_1}|E2|{2^+_1}>$ [eb]")

plt.show()

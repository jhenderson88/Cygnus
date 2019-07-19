import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.figure import figaspect

datafile = sys.argv[1]
fitfile = sys.argv[2]

dInit = []
dFina = []
dExpt = []
dData = []
dErro = []
exptnum = 0
with open(datafile) as file:
	for line in file:
		if line.find('EXPT') != -1 and line.find("!") == -1:
			exptnum += 1
		if line.find('EXPT') == -1 and line.find("!") == -1 and line.find('END'):
			dInit.append(float(line.split()[0]))
			dFina.append(float(line.split()[1]))
			dData.append(float(line.split()[2]))
			dErro.append(float(line.split()[3]))
			dExpt.append(exptnum)

fInit = []
fFina = []
fExpt = []
fData = []
exptnum = 0
with open(fitfile) as file:
	for line in file:
		if line.find("Init") != -1:
			exptnum += 1
		else:
			fInit.append(float(line.split()[0]))
			fFina.append(float(line.split()[1]))
			fData.append(float(line.split()[4]))
			fExpt.append(exptnum)

fig,ax = plt.subplots(len(fInit)/exptnum)
fig.subplots_adjust(hspace=0,wspace=0)

fData = np.ma.array(fData)
dData = np.ma.array(dData)
fExpt = np.ma.array(fExpt)
dExpt = np.ma.array(dExpt)

for a in range(len(fInit)/exptnum):

	tmp_dData = []
	tmp_fData = []
	tmp_dErro = []
	for i in range(len(fData)):
		if fInit[i] == fInit[a] and fFina[i] == fFina[a]:
			tmp_fData.append(fData[i])	
	for i in range(len(dData)):
		if dInit[i] == fInit[a] and dFina[i] == fFina[a]:
			tmp_dData.append(dData[i])	
			tmp_dErro.append(dErro[i])
	tmp_dData = np.asarray(tmp_dData)
	tmp_dErro = np.asarray(tmp_dErro)
	tmp_fData = np.asarray(tmp_fData)
	print a,'\n'
	print tmp_fData
	print tmp_dData
	print tmp_dErro
	print np.arange(len(tmp_fData))
	print '\n'
	ax[a].errorbar(np.arange(len(tmp_fData)),tmp_dData,yerr=tmp_dErro,fmt='ko',label='Expt: '+str(fInit[a])+r'$\rightarrow$'+str(fFina[a]))	
	ax[a].plot(np.arange(len(tmp_fData)),tmp_fData,'k--',label='Fitted: '+str(fInit[a])+r'$\rightarrow$'+str(fFina[a]))
	ax[a].legend()

ax[len(fInit)/exptnum/2].set_ylabel('Counts')
ax[len(fInit)/exptnum-1].set_xlabel('Experiment number')

plt.show()

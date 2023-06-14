try: np
except: import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
from matplotlib import cm

def PlotMoulin(MoulNum,filename):

	PosX=np.floor((95./4.5**2.)*((1-np.random.rand(MoulNum))*4.5)**2.)+5
	PosY=np.floor(18.*np.random.rand(MoulNum))+1

	print PosX
	print PosY

	PosX=[x for (y,x) in sorted(zip(PosY,PosX))]
	PosY=sorted(PosY)

	for i,val in enumerate(PosX):
		print PosX[i],PosY[i]

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax1.scatter(PosX,PosY)

	TotalInput=90.0
	WaterInput=TotalInput/float(MoulNum)

	with	open('../../input_functions/source/'+filename,'w') as moulinfile:
		Writer = csv.writer(moulinfile)
		for i,posx in enumerate(PosX):
			Writer.writerow([i,posx,PosY[i],WaterInput])
	
	moulinfile.close()


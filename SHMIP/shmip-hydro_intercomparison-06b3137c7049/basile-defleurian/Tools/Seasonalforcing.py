try: np
except: import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def ForcingPlot():

	Xcoords={x*1.0e3 for x in range(0,105,5)}
	thickness=np.zeros(np.shape(Xcoords))
	time={t*(3600*24)for t in range(0,365,1)}
	print	Xcoords
	Melt=np.zeros((np.shape(Xcoords)[0],np.shape(time)[0]))

	for i,x in enumerate(Xcoords):
		thickness[i]=1.0+(6.0*pow(x+5000.0,0.5)-6.0*pow(5000.0,0.5))
		
	rm= 2.15e3 #peak rate at 500m elevation m s-1
	rs= 5.0e-3#lapse rate m s-1 m-1 
	sm= 500.#peak position m
	tspr= 135.*(3600*24)#spring date in days
	taut= 244.*(3600*24)#fall date in days
	Dt= 50.*(3600*24)


	for i,elevation in enumerate(thickness):
		for date in time:
			Melt[i][date]=max(0,(rm+rs*sm)*((0.5*np.tanh(date-tspr)/Dt)-(0.5*np.tanh(date-taut)/Dt))-rs*elevation)

	fig = plt.figure()
	ax = fig.add_subplot(111)

	ax.plot(time,Melt[0][:])
	ax.plot(time,Melt[10][:])
	ax.plot(time,Melt[20][:])

	plt.show()

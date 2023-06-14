try: np
except: import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

def TroughGeom():

	Xcoords=np.arange(0,6050,50)
	Ycoords=np.arange(0,1050,50)
	surface=np.zeros((np.shape(Xcoords)[0],np.shape(Ycoords)[0]))
	bed=np.zeros((np.shape(Xcoords)[0],np.shape(Ycoords)[0]))
	thickness=np.zeros((np.shape(Xcoords)[0],np.shape(Ycoords)[0]))
	surfX=np.zeros((np.shape(Xcoords)[0],np.shape(Ycoords)[0]))
	surfY=np.zeros((np.shape(Xcoords)[0],np.shape(Ycoords)[0]))

	print np.shape(bed)


	Ysurf=-((Ycoords-500)/100.)**2
	Ybed=((Ycoords-500)/100.)**4
	Xsurf=1.0+(9.0*pow(Xcoords+50.0,0.5)-9.0*pow(50.0,0.5))
	Xbed=(-9.0*pow(6000-Xcoords,0.5)+9.0*pow(6000.0,0.5))

	
	for i,x in enumerate(Xcoords):
		for j,y in enumerate(Ycoords):
			surface[i][j]=1.0-((y-500)/100.)**2+(9.0*pow(x+50.0,0.5)-9.0*pow(50.0,0.5))
			bed[i][j]=((y-500)/100.)**4-9.0*(pow(6000-x,0.5)-pow(6000.0,0.5))
			surfX[i][j]=x
			surfY[i][j]=y
			thickness[i][j]=max(0,(surface[i][j]-bed[i][j]))

	fig = plt.figure()
	ax1 = fig.add_subplot(121,projection='3d')
	ax1.plot_wireframe(surfX,surfY,bed)
	ax1.plot_wireframe(surfX,surfY,(bed+thickness),color='r')
	ax2 = fig.add_subplot(122,)
	surf=ax2.contour(Xcoords,Ycoords,thickness.T)
	ax2.axis('scaled')

	fig.colorbar(surf)

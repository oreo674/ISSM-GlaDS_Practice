try: np
except: import numpy as np

# {{{def ComputeCenter(Elts,x,y):
def ComputeCenter(Elts,x,y):

	numberofnodes    = len(x);
	nodesperelements = np.shape(Elts)[0]
	numberofelements = np.shape(Elts)[1]

	#center have 2 columns, X,Y
	center = np.zeros((2,numberofelements))

	for j in range(0,nodesperelements):
		center[0,:]+=x[Elts[j,:]]
		center[1,:]+=y[Elts[j,:]]

	center=center/nodesperelements

	return center
# }}}
# {{{def ComputeMean(Elts,x,y,data):
def ComputeMean(Elts,x,y,data):

	numberofnodes    = len(x);
	numberofelements = np.shape(Elts)[1]
	nodesperelements = np.shape(Elts)[0]
	if numberofelements<nodesperelements:
		Elts = np.transpose(Elts)
		numberofelements = np.shape(Elts)[1]
		nodesperelements = np.shape(Elts)[0]

	weights = np.zeros(np.shape(data));
	areas   = np.zeros(numberofelements)
	if np.min(Elts)==1.:
		Elts    = Elts.astype(int) -1
	elif np.min(Elts)==0.:
		Elts    = Elts.astype(int)
	else:
		print 'Element tab has a bad format'
		raise
	#Computation of element area
	x1=x[Elts[0,:]]
	x2=x[Elts[1,:]]
	x3=x[Elts[2,:]]
	y1=y[Elts[0,:]]
	y2=y[Elts[1,:]]
	y3=y[Elts[2,:]]
	areas=(0.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)))

	#Computing weights associated to each node
	for i,element in enumerate(Elts):
		weights[element]+=areas[i]/nodesperelements
	WeightedDat = weights*data
	mean=sum(WeightedDat)/sum(weights)
	return mean
# }}}
# {{{def GetWeights(Elts,x,y):
def GetWeights(Elts,x,y):
	numberofnodes    = len(x);
	numberofelements = np.shape(Elts)[1]
	nodesperelements = np.shape(Elts)[0]
	if numberofelements<nodesperelements:
		Elts = np.transpose(Elts)
		numberofelements = np.shape(Elts)[1]
		nodesperelements = np.shape(Elts)[0]

	weights = np.zeros(np.shape(x));
	areas   = np.zeros(numberofelements)
	if np.min(Elts)==1.:
		Elts    = Elts.astype(int) -1
	elif np.min(Elts)==0.:
		Elts    = Elts.astype(int)
	else:
		print 'Element tab has a bad format'
		raise
	# Computation of element area
	x1=x[Elts[0,:]]
	x2=x[Elts[1,:]]
	x3=x[Elts[2,:]]
	y1=y[Elts[0,:]]
	y2=y[Elts[1,:]]
	y3=y[Elts[2,:]]
	areas=(0.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)))

	# Computing weights associated to each node
	for i,element in enumerate(Elts):
		weights[element]+=areas[i]/nodesperelements

	return weights
# }}}
# {{{def GetGridWeight(grid_x,grid_y)
def GetGridWeight(grid_x,grid_y):

		#grid weights to take the sides into account, just for rectangle
		Xstep=grid_x[1,0]-grid_x[0,0]
		Ystep=grid_y[0,1]-grid_y[0,0]

		Weights_x=np.zeros(np.shape(grid_x))
		Weights_x[np.where(grid_x+(Xstep*0.5)<(np.max(grid_x)))]=(Xstep*0.5)
		Weights_x[np.where(grid_x-(Xstep*0.5)>(np.min(grid_x)))]=Weights_x[np.where(grid_x-(Xstep*0.5)>(np.min(grid_x)))]+(Xstep*0.5)

		Weights_y=np.zeros(np.shape(grid_y))
		Weights_y[np.where(grid_y+(Ystep*0.5)<(np.max(grid_y)))]=(Ystep*0.5)
		Weights_y[np.where(grid_y-(Ystep*0.5)>(np.min(grid_y)))]=Weights_y[np.where(grid_y-(Ystep*0.5)>(np.min(grid_y)))]+(Ystep*0.5)
		Weights=Weights_x*Weights_y
		return Weights
# }}}

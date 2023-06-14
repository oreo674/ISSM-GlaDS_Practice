from netCDF4 import Dataset
import matplotlib.pyplot as plt
import csv
import numpy as np

def TsFit(order):
	rho_i=910.0
	rho_w=1.0e3
	g=9.8

	with open('../../../parameters/tuning_A3') as csvfile:
		spam_read=csv.reader(csvfile)
		pos_x=[0]
		eff_press=[0]
		for r,row in enumerate(spam_read):
			if r>0:
				pos_x.append(row[0])
				eff_press.append(row[1])

	pos_x=np.asarray(pos_x,'float')
	eff_press=np.asarray(eff_press,'float')
	thick=1.0+6.0*(pow(pos_x+5000.0,0.5)-pow(5000.0,0.5))

	head=((rho_i*g*thick)-eff_press)/(rho_w*g)
	
	grad_h=[(head[p+1]-head[p])/(pos_x[p+1]-pos_x[p]) for p,pos in enumerate(pos_x[:-1])]
	grad_h=np.abs(np.asarray(grad_h))

	grad_pos=[(pos_x[p+1]+pos_x[p])/2. for p,pos in enumerate(pos_x[:-1])]
	grad_pos=np.abs(np.asarray(grad_pos))

	wat_input=[5.79e-9*(100.0e3-x) for x in grad_pos]
	wat_input=np.asarray(wat_input)
	
	ts_fit=wat_input/grad_h
	section=np.where(pos_x[np.where(pos_x<80000.)]>10000.)
	z = np.polyfit(pos_x[section], ts_fit[section], order)
	p = np.poly1d(z)
	
	fiter=p(pos_x)
	fiter[np.where(pos_x<13920.)]=1.025-np.tanh((pos_x[np.where(pos_x<13920.)]+2000.)/4000.)

	# fig = plt.figure(tight_layout=True)
	# ax = fig.add_subplot(311)

	# ax.plot(pos_x,fiter,'r')
	# ax.plot(pos_x,p(pos_x),'k')
	# ax.plot(pos_x,1.025-np.tanh((pos_x+2000.)/4000.),'g+')
	# ax.plot(grad_pos,ts_fit,'bx')
	# ax2=fig.add_subplot(312)
	# ax2.plot(pos_x,thick,'b')
	# ax2.plot(pos_x,head,'r')
	# ax3=fig.add_subplot(313)
	# ax3.plot(pos_x,eff_press,'b')
	# plt.show(block=False)

	return p

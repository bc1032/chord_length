from mpl_toolkits.mplot3d import Axes3D as ax  # noqa: F401 unused import
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.stats import kde
from scipy import stats as st
import numpy as np
import pandas as pd
import pylab as pl
import itertools
from itertools import islice
import os
from scipy import stats
from mayavi import mlab
import matplotlib.colors as colors
import datetime


mkbinary = 'on'
numx,numy,numz=46,46,140
numbins=numx*numz*numy
numframes=210
cutoff = 0.005
fz=open("chordlength0.005.dat", "w")
fz.write("#frame, lambda, error\n")
if mkbinary == 'on':
	alldata = np.loadtxt('3dimagechordlength.dat')
	for frame in range(0,numframes):
		print(datetime.datetime.now().time())
		print(frame)
		with open('binarisedz.dat', 'w') as output_file:
			print(numbins)
			for nbin in range((frame*numbins), (frame+1)*numbins):
				if ( alldata[nbin,3] >= cutoff ):
					output_file.write(str(alldata[nbin,0]) + "	" + str(alldata[nbin,1]) + "	" + str(alldata[nbin,2]) + "	"  + str(1) +"\n")
					#break
				else:
					output_file.write(str(alldata[nbin,0]) + "	" + str(alldata[nbin,1]) + "	" + str(alldata[nbin,2]) + "	"  + str(0) +"\n")
		# 		#break
		# with open('binarisedy.dat', 'w') as output_file:
		# 	for k in range(0,numz):
		# 		for i in range(0,numx):
		# 			for j in range(0,numy):
		# 				coeff = i*numy*numz + j*numz + k
		# 				if ( data[coeff,3] >= cutoff ):
		# 					output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(1) +"\n")
		# 				else:
		# 					output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(0) +"\n")
		#
		# with open('binarisedx.dat', 'w') as output_file:
		# 	for k in range(0,numz):
		# 		for j in range(0,numy):
		# 			for i in range(0,numx):
		# 				coeff = i*numx*numz +j*numz+ k
		# 				if ( data[coeff,3] >= cutoff ):
		# 					output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(1) +"\n")
		# 				else:
		# 					output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(0) +"\n")

		data = np.loadtxt('binarisedz.dat')
		print(data.shape)
		#print(np.max(data[:,3]), np.min(data[:,3]))
		#print(data)
		#domains
		x = data[:,0]#(data[8,:]+data[9,:])/2. # [0.1, 5]
		y = data[:,1]#(data[6,:]+data[7,:])/2.            # [6, 9]
		z = data[:,2]#(data[4,:]+data[5,:])/2.           # [-1, 1]
		c = data[:,3]

		count = 0
		chordlengthz,chordlengthx,chordlengthy=[],[],[]
		for k in range(0,numbins):
			if ( k % (numz-1) == 0 ) and count >= 1:
				chordlengthz.append(count/2.)
				count = 0
			elif c[k] == 1:
				count+=1
			elif c[k] == 0 and count >= 1:
				chordlengthz.append(count/2.)
				count = 0
			elif c[k] == 0:
				#chordlengthz.append(count)
				count = 0
			#elif ( k % (numz-1) == 0 ) and count >= 1:
			#	chordlengthz.append(count/2.)
			#	count = 0
					#count = 0
		chrdln = np.sum(np.square(chordlengthz))/np.sum(chordlengthz)
		fz.write("%d	%0.5f\n" % (frame, chrdln ))
#		print(len(chordlengthz))
		print(chrdln, np.sum(np.square(chordlengthz)), np.sum(chordlengthz))
		print(np.square(chordlengthz))
#		print(np.std(chordlengthz))
fz.close()
## 	with open('binarisedy.dat', 'w') as output_file:
# 		for k in range(0,numz):
# 			for i in range(0,numx):
# 				for j in range(0,numy):
# 					coeff = i*numy*numz + j*numz + k
# 					if ( data[coeff,3] >= cutoff ):
# 						output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(1) +"\n")
# 					else:
# 						output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(0) +"\n")
#
# with open('binarisedx.dat', 'w') as output_file:
# 	for k in range(0,numz):
# 		for j in range(0,numy):
# 			for i in range(0,numx):
# 				coeff = i*numx*numz +j*numz+ k
# 				if ( data[coeff,3] >= cutoff ):
# 					output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(1) +"\n")
# 				else:
# 					output_file.write(str(data[coeff,0]) + "	" + str(data[coeff,1]) + "	" + str(data[coeff,2]) + "	"  + str(0) +"\n")

#
# data = np.loadtxt('binarisedz.dat')
# print(data.shape)
# print(np.max(data[:,3]), np.min(data[:,3]))
# #domains
# x = data[:,0]#(data[8,:]+data[9,:])/2. # [0.1, 5]
# y = data[:,1]#(data[6,:]+data[7,:])/2.            # [6, 9]
# z = data[:,2]#(data[4,:]+data[5,:])/2.           # [-1, 1]
# c = data[:,3]
#
# count = 0
# chordlengthz,chordlengthx,chordlengthy=[],[],[]
# for k in range(0,numbins):
# 	if c[k] == 1:
# 		count+=1
# 	elif c[k] == 0 and count >= 1:
# 		chordlengthz.append(count/2.)
# 		count = 0
# 	elif c[k] == 0:
# 		#chordlengthz.append(count)
# 		count = 0
# 	elif ( k % (numz-1) == 0 ) and count >= 1:
# 		chordlengthz.append(count/2.)
# 		count = 0
# 			#count = 0
# loopnum,loopk,loopi = 0,0,0
# data = np.loadtxt('binarisedy.dat')
# print(data.shape)
# print(np.max(data[:,3]), np.min(data[:,3]))
# #domains
# x = data[:,0]#(data[8,:]+data[9,:])/2. # [0.1, 5]
# y = data[:,1]#(data[6,:]+data[7,:])/2.            # [6, 9]
# z = data[:,2]#(data[4,:]+data[5,:])/2.           # [-1, 1]
# c = data[:,3]
#
# count = 0
#
# for j in range(0,numbins):
# 	coeff = j*numx*numz + (i+k)
# 	if c[j] == 1:
# 		count+=1
# 		loopnum += 1
# 	elif c[j] == 0 and count >= 1:
# 		chordlengthy.append(count/2.)
# 		count = 0
# 		loopnum += 1
# 	elif c[j] == 1:
# 		#chordlengthz.append(count)
# 		count = 0
# 		loopnum += 1
# 	elif ( j % (numy-1) == 0 ) and count >= 1:
# 		chordlengthy.append(count/2.)
# 		count = 0
#
# data = np.loadtxt('binarisedx.dat')
# print(data.shape)
# print(np.max(data[:,3]), np.min(data[:,3]))
# #domains
# x = data[:,0]#(data[8,:]+data[9,:])/2. # [0.1, 5]
# y = data[:,1]#(data[6,:]+data[7,:])/2.            # [6, 9]
# z = data[:,2]#(data[4,:]+data[5,:])/2.           # [-1, 1]
# c = data[:,3]
# count=0
# # for k in range(0,140):
# # 	for j in range(0,46):
# # 		for i in range(0,46):
# for i in range(0,numbins):
# 	if c[i] == 1:
# 		count+=1
# 	#	loopnum += 1
# 	elif c[i] == 0 and count >= 1:
# 		chordlengthx.append(count/2.)
# 		count = 0
# 	#	loopnum += 1
# 	elif c[i] == 1:
# 		#chordlengthz.append(count)
# 		count = 0
# 	#	loopnum += 1
# 	elif ( i % (numx-1) == 0 ) and count >= 1:
# 		chordlengthx.append(count/2.)
# 		count = 0
#
# #chordlength = np.mean(c[0:70])
# #print(chordlengthx)
# print(len(chordlengthz),len(chordlengthy),len(chordlengthx))
# print(np.mean(chordlengthz), np.mean(chordlengthy), np.mean(chordlengthx))
# print(np.std(chordlengthz), np.std(chordlengthy), np.std(chordlengthx))

#print(c[0:10])


# bounds=[0,.1,.2,.3]
# cmap = plt.colors.ListedColormap(['b','g','y','r'])
# norm = plt.colors.BoundaryNorm(bounds, cmap.N)
#im=ax.imshow(data[None], aspect='auto',cmap=cmap, norm=norm)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# img = ax.scatter(x, y, z, c=c,cmap='binary') #,norm=norm)#, cmap='binary')plt.spring()
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# # plt.colorbar(boundaries=np.linspace(0,1,5))
# #fig.colorbar(img,ticks=range(0,1,1), vmin=0, vmax=0.3)
# fig.colorbar(img,ticks = np.linspace(0, 1, 1, endpoint=True))#, boundaries=np.linspace(0,1,4))
#
# plt.savefig('/home/bc1032/Desktop/Work/Gels/polydisperse/morse20/rerun/test/binarydensity.pdf')
# plt.show()

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

mk3dimage = 'off'
mkchordfile = 'on'
rddmp = 1
numpart=13999
#framestep = 10000
numxyz = 16848
cut = 80
skplns = 9
#numframes = int((numxyz-1)/2)
browntime=step = 2.44034302759135
rangex, rangey, rangez = 22.8,22.8,70.0
xwallmin, ywallmin, zwallmin = -11.4, -11.4, -35.0
nbinsx=46
nbinsy=46
nbinsz = 140

frames = [0]
counter = 0
bad_words = ['13999']
j = 0
state = 0
i = 0

if mkchordfile == 'on':
	with open('xyzgel.xyz') as oldfile, open('cutinter.dat', 'w') as newfile:
		for line in oldfile:
			if any(bad_word in line for bad_word in bad_words) and (j % cut == 0) and len(line.split())<=2:# and (j != 0):
				j+=1
				state = 1
			elif any(bad_word in line for bad_word in bad_words) and (j % cut != 0) and len(line.split())<=2:# and (j != 0):# and (i % 13139 == 0)
				state = 0
				j+=1
			if state == 1:
				newfile.write(line)
				counter += 1
				#print(counter)

	#
	numframes = int(numxyz/cut)
	print(numframes)
	#i = numpart*numframes + (numframes*2)
	i=0
	with open('cutinter.dat') as oldfile, open('chordxyzinter.dat', 'w') as newfile:
		for line in oldfile:#.readlines()[numpart*numxyz:numpart*(numxyz + 1) + skplns]:
			if (len(line.split()) == 4):# and ( numpart*numframes + (numframes*2) <= i < (numpart*(numframes + 1) + (numframes*2) )):
				#	for i in range(numxyz, numxyz+numpart):
				newfile.write(line)
				i+=1
			elif i == numpart*(numxyz + 1):
				i=0
				#i+=1
	#			print(i)
				break

	os.remove('cutinter.dat')

	i=0
	pos = []
	data = np.loadtxt('chordxyzinter.dat')
	dim=0
	for dim in range(1,4):
		for i in range (0,numpart*numframes):
			pos.append(data[i,dim])
		#dim += 3
	posarray=np.array(pos)
	pos = posarray.reshape(3,numpart*numframes)
	npart = 1

	i,j,k, = 0,0,0

################NEW STUFF###################################
ilow, iup, imid = [],[],[]
jlow, jup, jmid = [],[],[]
klow, kup, kmid = [],[],[]
#expx = np.exp(abs)

for i in range(0, nbinsx):
	il, iu = (xwallmin + i*(rangex/nbinsx)), (xwallmin + (i+1)*(rangex/nbinsx))
	imid.append((il+iu)/2.0)

for j in range(0,nbinsy):
	jl, ju = (ywallmin + j*(rangey/nbinsy)), (ywallmin + (j+1)*(rangey/nbinsy))
	jmid.append((jl+ju)/2.0)

for k in range(0,nbinsz):
	kl,ku = (zwallmin+(k*(rangez/nbinsz))), (zwallmin + (k+1)*(rangez/nbinsz))
	klow.append(kl)
	kup.append(ku)
	kmid.append((kl+ku)/2.0)
#print(len(kmid))
locdens=[]
im,jm,km = 0,0,0
loopcount=1
#for coords in range(0, len(imid)*len(jmid)*len(kmid)):
with open("3dimagechordlength.dat", 'w') as output_file:
	for frame in range(0,numframes):
		print(datetime.datetime.now().time())
		print(frame)
		for im in range(0, len(imid)):
	#		print(datetime.datetime.now().time())
		#	print(im,jm)
			for jm in range(0, len(jmid)):
		#		print(im,jm)
				#print("testy")
				for km in range(0, len(kmid)):
					fxexponent,fyexponent,fzexponent = 0,0,0
					fcoorddist = 0
					funcr = []
					#print(len(funcr))
					#print("testz")
					#for npart in range((frame*numpart), (frame+1)*numpart):
					fxexponent = pos[0,(frame*numpart):(frame+1)*numpart] - imid[im]
					fyexponent = pos[1,(frame*numpart):(frame+1)*numpart] - jmid[jm] #abs(pos[1,npart] - jmid[jm])
					fzexponent = pos[2,(frame*numpart):(frame+1)*numpart] - kmid[km] #bs(pos[2,npart] - kmid[km])
					fcoorddist = (fxexponent**2 + fyexponent**2 + fzexponent**2)
						#print(fcoorddist)
					density = np.sum(np.exp(-16*fcoorddist))
						#print(np.exp(-16*fcoorddist))
					#print(np.sum(funcr))
					#locdens.append(density)
					#print(density)
					#locdens = np.sum(funcr)
					#output_file.write(str(imid[im]) + "	" + str(jmid[jm]) + "	" + str(kmid[km]) + "	"  + str(density[(im*jm)+loopcount*km])  +  str(im) + "	" + str(jm) + "	" + str(km) + "	"  + "\n")
					output_file.write(str(imid[im]) + "	" + str(jmid[jm]) + "	" + str(kmid[km]) + "	"  + str(density)  +  str(im) + "	" + str(jm) + "	" + str(km) + "	"  + "\n")

				#	print(locdens[(im*jm)+loopcount*km])
				loopcount += 1


if mk3dimage == 'on':
	with open("3dimagechordlength.dat", 'w') as output_file:
		for i in range(0, nbinsx):
			ilow, iup  = (xwallmin + i*(rangex/nbinsx)), (xwallmin + (i+1)*(rangex/nbinsx))
			imid = (ilow+iup)/2.0

			for j in range(0,nbinsy):
				jlow, jup = (ywallmin + j*(rangey/nbinsy)), (ywallmin + (j+1)*(rangey/nbinsy))
				jmid = (jlow+jup)/2.0

				for k in range(0,nbinsz):
					klow, kup  = (zwallmin+(k*(rangez/nbinsz))), (zwallmin + (k+1)*(rangez/nbinsz))
					kmid = (ilow+iup)/2.0
					for npart in range(0, len(pos[0])):
					#	if ( klow <= any(pos[2,part])  < kup ) and ( jlow <= pos[1,part]  < jup) and ( ilow <= pos[0,part]  < iup):
						if ( klow <= pos[2,npart]  < kup ) and ( jlow <= pos[1,npart]  < jup) and ( ilow <= pos[0,npart]  < iup):
							output_file.write(str(i) + "	" + str(j) + "	" + str(k) + "	"  + str(1) + "	" + str(klow) + "	" + str(kup) + "	" + str(jlow) + "	"  + "	" + str(jup) + "	" + str(ilow) + "	" + str(iup) + "	"  +"\n")
							break
							#npart = 13999
						elif npart < numpart:
							npart += 1
						else:
							output_file.write(str(i) + "	" + str(j) + "	" + str(k) + "	"  + str(0) + "	" + str(klow) + "	" + str(kup) + "	" + str(jlow) + "	"  + "	" + str(jup) + "	" + str(ilow) + "	" + str(iup) + "	"  +"\n")
							break
			#print(i)

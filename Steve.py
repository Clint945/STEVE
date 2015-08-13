#!/usr/bin/env python


import numpy as np
import os.path
import pickle

from findtools.find_files import (find_files, Match)
from matplotlib import pyplot as plt
from pylab import *



 					  #DEFINE ANALYSIS PARAMETERS#
#====================================================================#

#Define Misc parameters, Do Not Change Unless Needed.
numfolders = 1
numdatas = 1
foldercounter = -1
densityloopcounter = 0

#Define how simulation data is saved for file read-in.
#Example  '/bz__001.dat'
DataExt = '*.dat'				#'.dat'
VarLength = 4					#'bz__'
DumpLength = 3					#'001'
SaveLength = VarLength + DumpLength + len(DataExt) - 1

#Define some simulation parameters.
rspot = 5.0              				#E-6 m
reta = 10.0								#E-7 m
divangle = 0.7           				#Rads
loglambda = 5.0							#Collison Parameter
Etaconv = 0.3            				#Laser to Electron conversion efficiency
Lambda = 1.053                          #Laser Wavelength (E-6 m)
Ilaser = 2.5E19			 				#Wcm-2
Plaser = Ilaser*np.pi*(rspot/1E4)**2	#W

#Define the size of simulation. [in microns]
nx=int(150)
ny=int(200)
nz=int(200)
nn = nx*ny*nz


#Define simulation total length in picoseconds and number of timesteps.
totaltime = float(2.0)
timesteps = int(4)

#Define where slices for analysis are to be made. [in microns]
slice = list([1,30,75,120])
endonslice = slice[2]

#Define EEDF timestamps to be taken and resolution of binsize.
ts = list([150])        		#Picoseconds
res = int(5)  					#Kev  (eV for background EEDF)


#Define if data is to be pickled and what variables are to be processed.
selectiveprocessing = True
savepickledata = True
processlist = ('bz__','nf__','eta_')

#Plotting and single lineouts.
savefig_topdown = True
savefig_endon = True
plotlineout = False
ColdDensity = False			

#Diagnostics and required variables for each diagnostic. 
plotBdata = True					#Requires ('bz__')
plotlineBdata = False				#Requires ('bz__')
plotGammadata = False				#Requires ('bz__','tb__')  
useBfieldtheory = False				#Requires ('bz__','ez__','tb__')
plotEEDF = False			#BROKEN #Requires ('nf__','ez__','eta_')
plotEEDFbackground = False			#Requires ('nf__','tb__')
plotbackwalldensity = True			#Requires ('nf__')
plotbackwalltimehistory = False		#Requires ('nf__')


	   #MAKE INITIAL LISTS#
#================================#

#Generate names for graph legends.
slices = list()
for i in range(0,len(slice)):
	slices.append('X='+str(slice[i])+' Zeph B')
for i in range(0,len(slice)):
	slices.append('X='+str(slice[i])+' Theory B')
#endfor

#Create time array for plotting.
timelist = list()
for i in range (0,timesteps):
	timelist.append((i+1)*(totaltime/timesteps)) 
#endfor

#Calculate global fast electron temperature.
Tf = 8.19E-14*( (1.0+(Ilaser*Lambda**2)/1.38E18)**(0.5)-1.0)
Tf511 = Tf/(5.11E5*1.602E-19)





 					   #INITIATE ARRAYS AND LISTS#
#====================================================================#



#Create lists for basic processing
rawdata = list()
data = list()
dir = list()

#Create lists to keep the analysed data
timearray = list()
dirreplace = list()
graphnames = list()

maxbzarray = list()
meanbzarray = list()
slicebzarray = list()
newtheoryBzarray = list()
deltabz = list()
bzratio = list()

maxtbarray = list()
meantbarray = list()
slicetbarray = list()
tbarray = list()

maxezarray = list()
meanezarray = list()
sliceezarray = list()
ezarray = list()

maxetaarray = list()
meanetaarray = list()
sliceetaarray = list()
etaarray = list()

meannfarray = list()
backwalldensityarray = list()
backwallFWHMarray = list()
nfdensityfolderlist = list()
slicenfarray = list()
nfarray = list()

tfarray = list()

meantiarray = list()


oldgammaarray = list()
oldgammaplot1 = list()
oldgammaplot2 = list()
oldgammaplot3 = list()
oldgammaplot4 = list()

newgammaarray = list()
newgammaplot1 = list()
newgammaplot2 = list()
newgammaplot3 = list()
newgammaplot4 = list()

newtheorygammaarray = list()
newtheorygammaplot1 = list()
newtheorygammaplot2 = list()
newtheorygammaplot3 = list()
newtheorygammaplot4 = list()





                          #GETTING DIRECTORIES#
#====================================================================#

#/media/sjd549/Transcend/Scott/Assignments/MSc/MScResearch/1E20baseintensity/RawMovies/   {savename = savename[:(len(savename)-3)] + 'png'}

#./          {savename = savename[2:(len(savename)-3)] + 'png'}


#Find all files ending in dir recursively from current directory
sh_files_pattern = Match(filetype='f', name=DataExt)
found_files = find_files(path='./', match=sh_files_pattern)


#Organize the files into an array and sort them alphabetically
for found_file in found_files:
	dir.append(found_file)
#endfor
dir.sort()

#Remove various unwanted files from folders for plotting convinence.
for i in range(0,len(dir)):
	temp = dir[i]
	if temp[(len(dir[i])-10):] not in ('eta_00.dat','Z___00.dat','ni__00.dat'):
		if temp[(len(dir[i])-8):] not in ('zinp.dat','ninp.dat'):
			if temp[(len(dir[i])-12):] != 'totalnrg.dat':	
				dirreplace.append(temp)
			#endif
		#endif
	#endif
#endfor
dir = dirreplace

#Calculate the number of seperate simulations involved for plotting
for i in range(0, len(dir)-1):
	tempname1 = dir[i]
	tempname2 = dir[i+1]
	if tempname1[:(len(tempname1)-SaveLength)] != tempname2[:(len(tempname2)-SaveLength)]:
		numfolders += 1
	#endif
#endfor

#Calculate how many seperate bits of information are being analysed
for i in range(0, len(dir)-1):
	tempname1 = dir[i]
	tempname2 = dir[i+1]
	if tempname1[:(len(tempname1)-(SaveLength-VarLength))] != tempname2[:(len(tempname2)-(SaveLength-VarLength))]:
		numdatas += 1
	#endif
#endfor





                 #WELCOME TEXT AND ANALYSIS INFORMATION#
#====================================================================#
print''
print'---------------------------------------------------------------'
print'       _______.___________. ___________    ____  _______       '
print'      /       |           ||   ____\   \  /   / |   ____|      '
print'     |   (----`---|  |----`|  |__   \   \/   /  |  |__         '
print'      \   \       |  |     |   __|   \      /   |   __|        '
print'  .----)   |      |  |     |  |____   \    /    |  |____       '
print'  |_______/       |__|     |_______|   \__/     |_______|      '
print'                                                         v1.0.5'
print'---------------------------------------------------------------'
print ''
print '+=======================================================================+'
print 'Analysing',numfolders,'datasets containing',numdatas,'variables and',len(dir),'individual entries.'
print '+=======================================================================+'
print ''
print 'The following diagnostics will be used:'
print '---------------------------------------'
if savefig_topdown or savefig_endon == True:
	print'# Image Processing'
if plotlineout == True:
	print'# Plotting Lineouts'
if plotBdata == True:
	print'# Magnetic Evolution'
if ColdDensity == True:
	print'# Plotting Cold Ion Density'
if plotbackwalldensity == True:
	print'# Back-Wall Fast Electron Density'
if plotbackwalltimehistory == True:
	print'# Plotting Back-Wall nfast Time History'
if plotlineBdata == True:
	print'# Magnetic Gradients'
if plotGammadata and useBfieldtheory == True:
	print'# Collimation Thresholds inc Theory B-fields'
elif plotGammadata == True:
	print'# Collimation Thresholds'
if plotEEDF or plotEEDFbackground == True:
	print'# Electron Energy Dist Function'
print '---------------------------------------'





                         #DEFINE FUNCTIONS#
#====================================================================#


#Works out what figure has been requested, arranges data accordingly and saves.
def savefigs(array,savename,endonslice):
	if any ([savefig_topdown, savefig_endon]) == True:
		fig = plt.figure(figsize=(10,10),facecolor='w') 
		variable = savename[len(savename)-SaveLength:len(savename)-DumpExtLen]
		dataset = graphname1
	#endif
	if array.min() < 0 and savefig_topdown == True:
		plt.imshow((array[100,:,:]))
		plt.title('Midplane slice of '+dataset+' '+variable)
		plt.colorbar().set_label(ylabelselector(savename))
		plt.ylabel('distance [microns]')
		plt.xlabel('distance [microns]')

#		plt.show()
		plt.savefig(savename)
		plt.clf()
	elif savefig_topdown == True:
		plt.imshow(np.log(array[100,:,:]))
		plt.title('Midplane slice of '+dataset+' '+variable)
		plt.colorbar().set_label(ylabelselector(savename))
		plt.ylabel('distance [microns]')
		plt.xlabel('distance [microns]')

#		plt.show()
		plt.savefig(savename)
		plt.clf()
	#endif
	if array.min() < 0 and savefig_endon == True:
		plt.imshow((array[:,:,endonslice]))
		plt.title('Midplane slice of '+dataset+' '+variable+' at x='+str(endonslice))
		plt.colorbar().set_label(ylabelselector(savename))
		plt.ylabel('distance [microns]')
		plt.xlabel('distance [microns]')

		ext = (len(savename)-4)
		savename = savename[:ext] + 'b' + savename[ext:]
#		plt.show()
		plt.savefig(savename[:ext]+'_X='+str(endonslice)+'.png')
		plt.clf() 
	elif savefig_endon == True:
		plt.imshow(np.log(array[:,:,endonslice]))
		plt.title('Midplane slice of '+dataset+' '+variable+' at x='+str(endonslice))
		plt.colorbar().set_label(ylabelselector(savename))
		plt.ylabel('distance [microns]')
		plt.xlabel('distance [microns]')

		ext = (len(savename)-4)
		savename = savename[:ext] + 'b' + savename[ext:]
#		plt.show()
		plt.savefig(savename[:ext]+'_X='+str(endonslice)+'.png')
		plt.clf()

	#Plot a time history of the density profile if required.
	if plotbackwalltimehistory == True and len(nfdensityfolderlist) != 0:
		folderID = nfdensityfolderlist[densityloopcounter-1] 
		plt.plot(array[100,:,(nx-1)]/1.0E23)
		plt.title('Back-Wall nfast Density for '+graphnames[folderID]+' @'+str(timelist[(densityloopcounter-1)%(timesteps)])+'Ps')
		plt.legend(graphnames,loc=2)
		plt.ylabel('nfast [1E23 m-3]')
		plt.xlabel('distance [microns]')
		plt.grid(True)

		plt.savefig(graphnames[folderID]+'_nfast_walldensity'+str((densityloopcounter-1)%(timesteps))+'.png')
#		plt.show()
		plt.clf() 
	#endif
#enddef



#Selects an appropriate label based on the variable type.
def ylabelselector(savename):
	if savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('bz__','by__','bx__'):
		ylabel = ('B-field strength [T]')
	elif savename[len(savename)-SaveLength-1:len(savename)-DumpExtLen-1] in ('bz__','by__','bx__'):
		ylabel = ('B-field strength [T]')
	#endif
	if savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('ez__','ey__','ex__'):
		ylabel = ('E-field strength [V]')
	elif savename[len(savename)-SaveLength-1:len(savename)-DumpExtLen-1] in ('ez__','ey__','ex__'):
		ylabel = ('E-field strength [V]')
	#endif
	if savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('nf__','ni__'):
		ylabel = ('Number Density [1E23 m-3]')
	elif savename[len(savename)-SaveLength-1:len(savename)-DumpExtLen-1] in ('nf__','ni__'):
		ylabel = ('Number Density [1E23 m-3]')
	#endif
	if savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('tb__','ti__'):
		ylabel = ('Temperature [eV]')
	elif savename[len(savename)-SaveLength-1:len(savename)-DumpExtLen-1] in ('tb__','ti__'):
		ylabel = ('Temperature [eV]')
	#endif
	if savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('eta_'):
		ylabel = ('Resistivity [Ohm m]')
	elif savename[len(savename)-SaveLength-1:len(savename)-DumpExtLen-1] in ('eta_'):
		ylabel = ('Resistivity [Ohm m]')
	#endif
	return ylabel
#enddef




                 #UNPACKING DIR AND SETTING FILENAMES#
#====================================================================#

DumpExtLen = (len(DataExt)-1+DumpLength)

#Selectively process only a limited number of variables
if selectiveprocessing == True:
	dirreplace = list()
	for i in range(0,len(dir)):
		dirseperate = dir[i]
		#Collect new style Zeph save files
		if dirseperate[(len(dirseperate)-SaveLength):(len(dirseperate)-DumpExtLen)] in processlist:
			dirreplace.append(dirseperate)
		#Collect original style Zeph save files
		elif dirseperate[(len(dirseperate)-(SaveLength-1)):(len(dirseperate)-(DumpExtLen-1))] in processlist:
			dirreplace.append(dirseperate)
		#endif
	#endfor
	dir = dirreplace
	print ''
	print 'Analysing only',processlist
	print 'processing',len(dir),'total entries.'
	print ''
#endif


#Create backwalldensity storage arrays depending on number of datas processed.
if plotbackwalldensity or plotbackwalltimehistory == True:
	#Create correctly sized list for the data
	for i in range(0,numfolders):
		backwalldensityarray.append(list())
		for k in range(0,ny):
			backwalldensityarray[i].append(0)
		#endfor
	#endfor

	#Create an array to effectivly partition backwalldensityarray between folders.
	if type(processlist) == str:
		lenprocesslist = 1
	else:
		lenprocesslist = len(processlist)
	#endif
	for j in range(0,numfolders):
		for k in range(0,((len(dir)/lenprocesslist)/numfolders)):
			nfdensityfolderlist.append(j)
		#endfor
	#endfor
#endif

 

print '----------------------'
print 'Beginning Data Read-in'
print '----------------------'

#Start main data sorting and storage loop.
for i in range(0,len(dir)):
	#Refresh old lists to prevent memory issues.
	rawdata = list()
	data = list()
	array = list()

	#Prepare the savename for the figure as saved in this section.
	savename = dir[i]
	savename = savename[2:(len(savename)-len(DataExt)+2)] + 'png'
	
	#Prepare a list of names for autotitles in the graph section.
	graphname1 = dir[i]
	graphname1 = graphname1[2:len(graphname1)-12]
	try:
		graphname2 = dir[i+1]
		graphname2 = graphname2[2:len(graphname2)-12]
	except:
		This_needs_to_be_here_for_some_reason=1
	#endtry
	if graphname1 != graphname2 or len(graphnames) == 0:
		foldercounter += 1
		try:
			graphnames.append(graphname2)
		except:
			graphnames.append(graphname1)
		#endtry
	#endif

	#Load previously analysed data or pickle new data
	try:
		with open(savename[:(len(savename)-len(DataExt)+2)]+'p', "rb" ) as f:
			array = pickle.load(f)
		#endwith
	except:
		#extract data and remove trailing \n from numbers.
		rawdata = open(dir[i]).readlines()
		for j in range(0,nn):
			rawdata[j] = str.strip(rawdata[j])
			try:
				data.append(float(rawdata[j]))
			except:
				data.append(0)
			#endtry
		#endfor

		#reshape array for plotting and pickle array
		array = np.reshape(data,(nz,ny,nx))
		if savepickledata == True:
			with open(savename[:(len(savename)-len(DataExt)+2)]+'p', "wb" ) as f:
				pickle.dump(array, f)
			#endwith
		#endif
	#endtry
	

    				   #SINGLE LINE DATA ANALYSIS#
#====================================================================#




	#Plot a single linecut at requested slice.
	if plotlineout == True:
		plt.plot(array[100,:,endonslice])
		variable = savename[len(savename)-SaveLength:len(savename)-DumpExtLen]
		dataset = graphname1
		plt.title('Line plot for '+dataset+' '+variable+' at x='+str(endonslice))
		plt.ylabel(ylabelselector(savename))
		plt.xlabel('Simulation Width [microns]')
		plt.grid(True)
		plt.show()
	#endif

	#Compare |B| against |gradB| across a requested slice.
	if savename[len(savename)-SaveLength:len(savename)-DumpExtLen] == 'bz__' and plotlineBdata == True:
		for j in range(0,ny-1):
			deltabz.append(abs(array[100,j+1,endonslice])-abs(array[100,j,endonslice]))
		#endfor
		deltabz.append(0)
		for j in range(1,(ny)):
			if abs(deltabz[j-1])/abs(array[100,j,endonslice]) > 10:
				bzratio.append(0)
			else:
				bzratio.append(abs(deltabz[j-1])/abs(array[100,j,endonslice]))
			#endif
		#endfor
		bzratio.append(0)

		l=len(deltabz)

		plt.subplot(2,1,1)
		plt.plot(range(0,ny),array[100,:,endonslice])
		plt.plot(range(0,ny),deltabz[-200:])
		plt.legend(['B','Grad B'])
		plt.ylabel('B-field strength [T]')
		#plt.xlabel('Simulation Width [microns]')
		plt.grid(True)
		
		plt.subplot(2,1,2)
		plt.plot(range(0,ny),bzratio[-200:])
		plt.ylabel('GradB to B ratio')
		plt.xlabel('Simulation Width [microns]')
		plt.grid(True)

#		plt.show()
		plt.savefig(graphnames[foldercounter]+'_GradBratios'+str(i+1)+'.png')
		plt.clf()	
		#endfor
	#endif




			#DATA SORTING AND SAVING FOR FURTHER ANALYSIS#
#====================================================================#




	#Seperate data types into specific arrays for analysis.
	if savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('bz__','/bz_'):
		maxbzarray.append(array[:,:,:].max())
		meanbzarray.append(np.mean(abs(array)))

		#Gamma Data Collectors
		slicebzarray.append(array[100,:,slice[0]].max())
		slicebzarray.append(array[100,:,slice[1]].max())
		slicebzarray.append(array[100,:,slice[2]].max())
		slicebzarray.append(array[100,:,slice[3]].max())
		
		savefigs(array,savename,endonslice)


	elif savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('ez__','/ez_'):
		maxezarray.append(array.max())
		meanezarray.append(np.mean(abs(array)))

		#Gamma Data Collectors
		sliceezarray.append(array[100,:,slice[0]].max())
		sliceezarray.append(array[100,:,slice[1]].max())
		sliceezarray.append(array[100,:,slice[2]].max())
		sliceezarray.append(array[100,:,slice[3]].max())
		
		savefigs(array,savename,endonslice)
		
		#EEDF timestamp selector and data collection
		tempname = savename[len(savename)-5:len(savename)-8:-1]
		tempname = tempname[::-1]
		if plotEEDF == True and int(tempname) in ts:
			array = np.reshape(array,(1,1,nn)).tolist()
			ezarray.append(array[0][0][:])
		#endif


	elif savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('eta_','/eta'):
		maxetaarray.append(array.max())
		meanetaarray.append(np.mean(abs(array)))

		#Gamma Data Collectors
		sliceetaarray.append(array[100,:,slice[0]].min())
		sliceetaarray.append(array[100,:,slice[1]].min())
		sliceetaarray.append(array[100,:,slice[2]].min())
		sliceetaarray.append(array[100,:,slice[3]].min())
		
		savefigs(array,savename,endonslice)

		#EEDF timestamp selector and data collection
		tempname = savename[len(savename)-5:len(savename)-8:-1]
		tempname = tempname[::-1]
		if plotEEDF == True and int(tempname) in ts:
			array = np.reshape(array,(1,1,nn)).tolist()
			etaarray.append(array[0][0][:])
		#endif


	elif savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('nf__','/nf_'):
		meannfarray.append(np.mean(abs(array)))

		#fast density arbitary slice collectors
		slicenfarray.append(array[100,:,slice[0]].mean())
		slicenfarray.append(array[100,:,slice[1]].mean())
		slicenfarray.append(array[100,:,slice[2]].mean())
		slicenfarray.append(array[100,:,slice[3]].mean())
		
		#Back-wall density collector
		if plotbackwalldensity == True:
			folderID = nfdensityfolderlist[densityloopcounter]
			for k in range(0,ny):
				#Save back wall density profile for further analysis and plotting.
				backwalldensityarray[folderID][k] += array[100,k,(nx-1)]/(1.0E23*((len(dir)/lenprocesslist)/numfolders))

			#Calculate the FWHM to plot the speed of collimation.
			FWHM = 0
			Max = max(backwalldensityarray[folderID])
			for q in range(0,ny):
				if backwalldensityarray[folderID][q] > Max/2:
					FWHM += 1
				#endif
			#endfor
			backwallFWHMarray.append(FWHM)
		#endif
		densityloopcounter += 1

		savefigs(array,savename,endonslice)

		#EEDF timestamp selector and data collection
		tempname = savename[len(savename)-5:len(savename)-8:-1]
		tempname = tempname[::-1]
		if plotEEDF or plotEEDFbackground == True and int(tempname) in ts:
			array = np.reshape(array,(1,1,nn)).tolist()
			nfarray.append(array[0][0][:])
		#endif


	elif savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('tb__','/tb_'):
		maxtbarray.append(array.max())
		meantbarray.append(np.mean(abs(array)))

		#Gamma Data Collectors
		slicetbarray.append(array[100,:,slice[0]].max())
		slicetbarray.append(array[100,:,slice[1]].max())
		slicetbarray.append(array[100,:,slice[2]].max())
		slicetbarray.append(array[100,:,slice[3]].max())

		savefigs(array,savename,endonslice)

		#EEDF timestamp selector and data collection
		tempname = savename[len(savename)-len(DataExt):len(savename)-DumpExtLen-1:-1]
		tempname = tempname[::-1]
		if plotEEDFbackground == True and int(tempname) in ts:
			array = np.reshape(array,(1,1,nn)).tolist()
			tbarray.append(array[0][0][:])
		#endif


	elif savename[len(savename)-SaveLength:len(savename)-DumpExtLen] in ('ti__','/ti_'):
		meantiarray.append(np.mean(abs(array)))
		
		savefigs(array,savename,endonslice)


	else:
		savefigs(array,savename,endonslice)
	#endif


	



	#Percentage complete printout.
	oldpercentage = int(((float(i)/len(dir))*100.0))
	newpercentage = int(((float(i+1)/len(dir))*100.0))
	if oldpercentage != newpercentage:
		print newpercentage,'%'
	#endif
#endfor

#Alert user that process has ended or continue with selected diagnostics.
if any([plotBdata, plotlineBdata, plotGammadata, useBfieldtheory, plotEEDF, plotEEDFbackground, ColdDensity, plotbackwalldensity, plotbackwalltimehistory]) == True:
	print '----------------------------------------'
	print 'Data Readin Complete, Starting Analysis.'
	print '----------------------------------------'
else:
	print '------------------'
	print 'Analysis Complete.'
	print '------------------'
#endif




						#COLD-DENSITY ANALYSIS#
#====================================================================#


dir2 = list()
Zcounter = 0
if ColdDensity == True:

	#Find files again to get the ni__00 files that were ignored previously
	sh_files_pattern = Match(filetype='f', name=DataExt)
	found_files = find_files(path='./', match=sh_files_pattern)
	for found_file in found_files:
		dir2.append(found_file)
	#endfor
	dir2.sort()

	#For each file attempt to match ni__00 name
	for i in range(0,len(dir2)):
		temp = dir2[i]
		if temp[(len(dir2[i])-10):] == 'ni__00.dat':
			data = list()
			if Zcounter == 0:
				Z = 13.0
			else:
				Z = 2.66
			Zcounter += 1
			#endif

			#Open data once found and read out to an array as usual.
			rawdata = open(dir2[i]).readlines()
			for j in range(0,nn):
				rawdata[j] = str.strip(rawdata[j])
				try:
					data.append(float(rawdata[j]))
				except:
					data.append(0)
				#endtry
			#endfor

			#arrange and print with basic colourbar and max cold density.
			array = np.reshape(data,(nz,ny,nx))
			plt.imshow((array[100,:,:]))
			plt.colorbar().set_label(Z*array.max())
			plt.show()
			plt.clf
		#endif
	#endfor


	print '-------------------------------'
	print 'Cold-Density Analysis Complete.'
	print '-------------------------------'
#endif





						#BACKWALL DENSITY ANALYSIS#
#====================================================================#

if plotbackwalldensity == True:

	#Plot individual graph of each folder for single use.
	for i in range(0,numfolders):
		fig = plt.figure(figsize=(10,10),facecolor='w') 

		plt.plot(backwalldensityarray[i])
		plt.title('Time Integrated Back-Wall nfast Density for '+str(graphnames[i])+' after '+str(totaltime)+'Ps')
		plt.ylabel('nfast [1E23 m-3]')
		plt.xlabel('distance [microns]')
		plt.grid(True)

		plt.savefig(graphnames[i]+'_nfast_walldensity.png')
#		plt.show()
		plt.clf()
	#endfor

	#Print one graph for comparison of all folders of data.
	if numfolders not in ([0,1]):
		fig = plt.figure(figsize=(10,10),facecolor='w') 
		for i in range(0,len(backwalldensityarray)):

			plt.plot(backwalldensityarray[i])
			plt.title('Time Integrated Comparison of Back-Wall nfast Density for \n '+str(graphnames)+' @'+str(totaltime)+'Ps')
			plt.legend(graphnames,loc=2)
			plt.ylabel('nfast [1E23 m-3]')
			plt.xlabel('distance [microns]')
			plt.grid(True)
		#endfor

		plt.savefig(graphnames[i]+'_nfast_walldensitycomparison.png')
#		plt.show()
		plt.clf()
	#endif

	#Plot the FWHM time history to see speed of collimation.
	l=len(backwallFWHMarray)
	fig = plt.figure(figsize=(10,10),facecolor='w')

	for k in range(1,numfolders+1):
		plt.plot(backwallFWHMarray[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.title('Time Integrated Back-Wall Density Profile FWHM '+str(graphnames[i])+' after '+str(totaltime)+'Ps')
		plt.legend(graphnames,loc=2)
		plt.ylabel('FWHM [microns]')
		plt.xlabel('time [Ps]')
		plt.grid(True)
	#endfor

	plt.savefig(graphnames[i]+'_nfast_walldensityFWHM.png')
#	plt.show()
	plt.clf()


	print '-------------------------------'
	print 'Fast-Density Analysis Complete.'
	print '-------------------------------'
#endif


                           #B-FIELD ANALYSIS#
#====================================================================#


l=len(maxbzarray)

if plotBdata == True:

	fig = plt.figure(figsize=(15,10),facecolor='w') 
	plt.subplot(2,1,1)
	for k in range(1,numfolders+1):
		plt.plot(timelist,maxbzarray[((k-1)*l)/numfolders:(k*l)/numfolders])
	#endfor
	plt.title('Maximum B-field evolution for '+str(graphnames)+' cases')
	plt.legend(graphnames,loc=2)
	plt.ylabel('B-field strength [T]')
	#plt.xlabel('Time [ps]')
	plt.grid(True)
	
	plt.subplot(2,1,2)
	for k in range(1,numfolders+1):
		plt.plot(timelist,meanbzarray[((k-1)*l)/numfolders:(k*l)/numfolders])
	#endfor
	plt.title('Mean B-field evolution for '+str(graphnames)+' cases')
	plt.legend(graphnames,loc=2)
	plt.ylabel('B-field strength [T]')
	plt.xlabel('Time [ps]')
	plt.grid(True)

	plt.savefig(graphnames[k-1]+'_Bfields.png')
#	plt.show()
	plt.clf()


	print '--------------------------'
	print 'B-Field Analysis Complete.'
	print '--------------------------'
#endif




                            #GAMMA ANALYSIS#
#====================================================================#



if plotGammadata == True:

#	#TEMPERATURE THEORY IS WRONG HERE, DO NOT USE.
#	#Calculate temperature from eta and ez arrays using my method.
#	Tf511array = list()
#	rbeam = 16
#	for i in range(len(sliceezarray)):
#		eta = float(sliceetaarray[i])
#		ez = float(sliceezarray[i])
#		
#		#endif
#		print abs((Etaconv*Plaser*eta)/(np.pi*((rbeam/1E6)**2)*ez)/511000)
#		tfarray.append(abs((Etaconv*Plaser*eta)/(np.pi*((rspot/1E6)**2)*ez)))	
#	#endfor


	#Calculate BellGamma (old Gamma) from simulated B-field data.
	for i in range(0,(len(slicebzarray))):
		temp=(1/(np.sqrt(Tf511)))*1/(np.sqrt(2+(Tf511))) 
		oldgammaarray.append(.0591*rspot*(1/(divangle**2))*temp*(slicebzarray[i]/100))
	#endfor 

	#Calculate ScottGamma (new Gamma) from simulated B-field data.
	for i in range(0,(len(oldgammaarray))):
		radii= 1 + rspot/(reta*0.1)
		newgammaarray.append(oldgammaarray[i]*radii)
	#endfor

	for i in range(0,len(oldgammaarray),4):
		oldgammaplot1.append(oldgammaarray[i])
		oldgammaplot2.append(oldgammaarray[i+1])
		oldgammaplot3.append(oldgammaarray[i+2])
		oldgammaplot4.append(oldgammaarray[i+3])
	#endfor

	for i in range(0,len(newgammaarray),4):
		newgammaplot1.append(newgammaarray[i])
		newgammaplot2.append(newgammaarray[i+1])
		newgammaplot3.append(newgammaarray[i+2])
		newgammaplot4.append(newgammaarray[i+3])
	#endfor




	if useBfieldtheory == True:
		
		#Create a timearray to create B-fields.
		for j in range(0,numfolders):
			for i in range(0,len(sliceezarray)/numfolders):
				timearray.append(timelist[int(i/4)]*1E-12)
			#endfor
		#endfor

		#Generate theoretical prediction for B to be compared against simulation B.
		for i in range(0,len(sliceezarray)):
			newtheoryBzarray.append(abs(sliceezarray[i])*timearray[i]*((1/(rspot*1E-6))+(1/(reta*1E-7))))
		#endfor

		#Calculate ScottGamma (new Gamma) from theoretical B-field data.
		for i in range(0,(len(slicebzarray))):
			temp=(1/(np.sqrt(Tf511)))*1/(np.sqrt(2+(Tf511)))
			radii=(rspot+((0.5)*(reta/10)))*(1/(divangle**2))
			newtheorygammaarray.append(.0591*radii*temp*(newtheoryBzarray[i]))
		#endfor

		for i in range(0,len(newgammaarray),4):
			newtheorygammaplot1.append(newtheorygammaarray[i])
			newtheorygammaplot2.append(newtheorygammaarray[i+1])
			newtheorygammaplot3.append(newtheorygammaarray[i+2])
			newtheorygammaplot4.append(newtheorygammaarray[i+3])
		#endfor
	#endif
	


        #PLOTTING GAMMA#
#================================#


	l=len(oldgammaplot1)


	#Plot oldgamma and newgamma for each slice over one all folders of data
	for k in range(1,numfolders+1):
		fig = plt.figure(figsize=(15,7),facecolor='w')
		plt.plot(timelist,oldgammaplot1[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.plot(timelist,oldgammaplot2[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.plot(timelist,oldgammaplot3[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.plot(timelist,oldgammaplot4[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.title('Temporal plot of maximum Bell_Gamma at given slices for '+graphnames[k-1])
		plt.legend(slices[0:len(slice)], loc=2)
		plt.ylabel('Gamma []')
		plt.xlabel('Time [ps]')
		plt.grid(True)

		plt.savefig(graphnames[k-1]+'_BellGammas.png')
#		plt.show()
		plt.clf()

		fig = plt.figure(figsize=(15,7),facecolor='w') 
		plt.plot(timelist,newgammaplot1[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.plot(timelist,newgammaplot2[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.plot(timelist,newgammaplot3[((k-1)*l)/numfolders:(k*l)/numfolders])
		plt.plot(timelist,newgammaplot4[((k-1)*l)/numfolders:(k*l)/numfolders])
			
		plt.title('Temporal plot of maximum New_Gamma at given slices for '+graphnames[k-1])
		plt.legend(slices, loc=2)
		plt.ylabel('Gamma []')
		plt.xlabel('Time [ps]')
		plt.grid(True)
	
		plt.savefig(graphnames[k-1]+'_NewGammas.png')
#		plt.show()
		plt.clf()
	#endfor

	
	#Plot Gamma comparisons at each slice over all folders of data
	if useBfieldtheory == True:
		for k in range(1,numfolders+1):
			fig = plt.figure(figsize=(15,10),facecolor='w') 

			plt.subplot(2,2,1)
			plt.plot(timelist,oldgammaplot1[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.plot(timelist,newtheorygammaplot1[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.title(str(slices[0]))
			plt.legend(['Oldgamma','Newgamma'], loc=2)
			plt.ylabel('Gamma []')
			#plt.xlabel('Time [ps]')
			plt.grid(True)

			plt.subplot(2,2,2)
			plt.plot(timelist,oldgammaplot2[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.plot(timelist,newtheorygammaplot2[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.title(str(slices[1]))
			plt.legend(['Oldgamma','Newgamma'], loc=2)
			#plt.ylabel('Gamma []')
			#plt.xlabel('Time [ps]')
			plt.grid(True)

			plt.subplot(2,2,3)
			plt.plot(timelist,oldgammaplot3[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.plot(timelist,newtheorygammaplot3[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.title(str(slices[2]))
			plt.legend(['Oldgamma','Newgamma'], loc=2)
			plt.ylabel('Gamma []')
			plt.xlabel('Time [ps]')
			plt.grid(True)

			plt.subplot(2,2,4)
			plt.plot(timelist,oldgammaplot4[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.plot(timelist,newtheorygammaplot4[((k-1)*l)/numfolders:(k*l)/numfolders])
			plt.title(str(slices[3]))
			plt.legend(['Oldgamma','Newgamma'], loc=2)
			#plt.ylabel('Gamma []')
			plt.xlabel('Time [ps]')
			plt.grid(True)

			plt.suptitle('Comparative plots of both Bell and New Gamma values for '+graphnames[k-1], size=16)

			plt.savefig(graphnames[k-1]+'_BellvsTheory_Gammas.png')
#			plt.show()
			plt.clf()
		#endfor
	#endif


	print '------------------------------'
	print 'Collimation Analysis Complete.'
	print '------------------------------'
#endif




							#EEDF ANALYSIS#
#====================================================================#


EEDFeVarray = list()
EEDFarray = list()
tfarray = list()



if plotEEDF == True:

#	#TEMPERATURE THEORY IS WRONG HERE, DO NOT USE.
#	#Calculate temperature from eta and ez arrays using my method.
#	for i in range(len(ezarray)):
#		eta = float(etaarray[i])
#		ez = float(ezarray[i])
#		if ez == 0.0:
#			ez = meanezarray[0]
#		#endif
#		print ez
#		print eta
#		print abs((Etaconv*Plaser*eta)/(np.pi*((rspot/1E6)**2)*ez)/1E6)
#		tfarray.append(abs((Etaconv*Plaser*eta)/(np.pi*((rspot/1E6)**2)*ez)))	
#	#endfor

	#Calculate the number of temperature 'bins' for given resolution.
#	maxtemp = max(tfarray)/1E12
#	print maxtemp
#	print len(range(0,int(maxtemp),(res*1000)))
#	EEDFeVarray.append(range(0,int(maxtemp),(res*1000)))
#	EEDFeVarray = EEDFeVarray[0]

	print'# Non-background EEDF analysis is disabled.'

#endif



		#EEDF BACKGROUND#
#================================#


if plotEEDFbackground == True:

	for k in range(len(nfarray)):

		EEDFeVarray = list()
		EEDFarray = list()

		#Calculate the number of temperature 'bins' for given resolution.
		maxtemp = max(tbarray[k])
		EEDFeVarray.append(range(0,int(maxtemp),res))
		EEDFeVarray = EEDFeVarray[0]

		#Create array of correct length to hold number of particles at given temp.
		for i in range(len(EEDFeVarray)):
			EEDFarray.append(0)
		#endfor
	
		#Go over all elements in the electron population and organize into bins
		for i in range(0,len(nfarray[k])):
			for j in range(0,len(EEDFeVarray)):
				if tbarray[k][i] >= EEDFeVarray[j]:
					EEDFarray[j] += nfarray[k][i]/nn
					j = len(EEDFeVarray)
				#endif		
			#endfor
		#endfor


		#Individual Plotting
		plt.plot(EEDFeVarray[1:],EEDFarray[1:])
		plt.title(graphnames[k-1]+'_T='+str(ts[int(k/len(graphnames))])+'_EEDFbackground')
		plt.ylabel('Electron number density [m-3]')
		plt.xlabel('Electron Energy [eV]')
		plt.grid(True)

		plt.savefig(graphnames[k-1]+'_T='+str(ts[int(k/len(graphnames))])+'_EEDF.png')
#		plt.show()
		plt.clf()


		#Percentage Printout
		percentage1 = int( (float(k)/len(nfarray))*100 )
		percentage2 = int( ((k+1.0)/len(nfarray))*100 )
		if percentage1 != percentage2:
			print percentage2,'%'
		#endif
	#endfor

	print '-----------------------'
	print 'EEDF Analysis Complete.'
	print '-----------------------'
#endif	




							   #CODE DUMP#
#====================================================================#


#Random Stuff Goes here!



#Location Collector for Temperature remover.
#Not Currently needed as temperature extraction is erronious.
#		print np.where(array == slicebzarray[(4*i)-3])
#		print np.where(array == slicebzarray[(4*i)-2])
#		print np.where(array == slicebzarray[(4*i)-1])
#		print np.where(array == slicebzarray[(4*i)])



#OLD METHOD OF CALCULATING GAMMA USING TBARRAY, MAY BE USEFUL LATER.
#	#Calculate BellGamma (old Gamma) from simulated B-field data.
#	for i in range(0,(len(slicebzarray))):
#		temp=(1/(np.sqrt(slicetbarray[i]/511)))*1/(np.sqrt(2+(slicetbarray[i]/511))) 
#		oldgammaarray.append(.0591*rspot*(1/(divangle**2))*temp*(slicebzarray[i]/100))
#	#endfor 
#
#	#Calculate ScottGamma (new Gamma) from simulated B-field data.
#	for i in range(0,(len(slicebzarray))):
#		temp=(1/(np.sqrt(slicetbarray[i]/511)))*1/(np.sqrt(2+(slicetbarray[i]/511)))
#		radii=(rspot+((0.5)*(reta/10)))*(1/(divangle**2))
#		newgammaarray.append(.0591*radii*temp*(slicebzarray[i]/100))
#	#endfor

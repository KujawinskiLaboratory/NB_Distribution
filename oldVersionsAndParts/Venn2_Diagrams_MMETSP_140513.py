#!/usr/bin/env python

'''
Created on November 9, 2013

VennDiagram of Species composition

@author: harrietalexander
'''

import sys
import glob
import os
import numpy as np
import itertools
import re
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_pdf as PdfPages
import matplotlib as mpl

from random import shuffle
os.chdir('/Users/harrietalexander/Analysis/NB_Distribution/')

ManClus=csv.reader(open('MMETSP_140513_Cluster.tab'),delimiter='\t') #Manual clusters

AllClus=csv.reader(open('MMETSP_HigherOrder.tab'),delimiter='\t')
ClusCount=csv.reader(open('SummedSpecies.tab'), delimiter='\t')
mpl.rcParams['pdf.fonttype'] = 42

#LOAD CLUSTERING
MC_hash={}
for line in ManClus:
	Clus=line[-1].strip()
	MMETSP=line[0].strip()
	if Clus in MC_hash:
		MC_hash[Clus].append(MMETSP)
	else:
		MC_hash[Clus]=[MMETSP]
All_hash=[]
hash={}
#LOAD MMETSP DATA
MMETSP_Hash={}
for line in AllClus:
	key=line[0].strip()
	Class=line[1].strip()
	Order=line[2].strip()
	Family=line[3].strip()
	Genus=line[4].strip()
	MMETSP_Hash[key]=[Class,Order,Family,Genus]
#LOUAD COUNTS
Count_hash={}
for line in ClusCount:
	c=0
	numList=[]
	for i in line:
		i=i.strip()
		if c==0:
			Clus=i
			c+=1
		else:
			numList.append(int(i))
	Count_hash[Clus]=numList
newHash={}
for key in MC_hash:
	gID=MC_hash[key][0]
	newHash[key]=MMETSP_Hash[gID]

#CREAT HASH OF SPECIES FOR EACH CATAGORY
Chash={}
Ohash={}
Fhash={}
Ghash={}
Ahash=[Chash,Ohash,Fhash,Ghash]
for key in newHash:
	MM=newHash[key]
	for h,m in zip(Ahash,MM):
		if m in h:
			h[m].append(key)
		else:
			h[m]=[key]
#CREAT HASH OF COUNTS FOR EACH CATAGORY
Ccount={}
Ocount={}
Fcount={}
Gcount={}
Acount=[Ccount,Ocount,Fcount,Gcount]

for m,c in zip(Ahash, Acount):
	for key in m:
		print key
		for i in m[key]:
			print i
			if i in Count_hash:
				if key in c:
					c[key]=[x+y for x,y in zip(c[key], Count_hash[i])]
				else:
					c[key]=Count_hash[i]
					
# #Plot for E7
outdir='/Users/harrietalexander/Dropbox/NB_Paper/'
nns=['Family','Class','Order','Genus']
cc=1
for Count_hash,nn in zip(Acount,nns):
	nums=[]
	for i in numList:
		nums.append([])
	Name_list=[]
	for key in Count_hash:
		Name_list.append(key)
		for i in range(len(numList)):
			nums[i].append(Count_hash[key][i])		
	fig=plt.figure(cc)
	fig.suptitle(nn)
	names=['S1', 'S2', 'S3', 'S4','S5', '+N', '-N', '+P', '-P','C']
	axs=[]
	for x in range(len(names)):
		ax1=fig.add_subplot(len(names),1,x+1,aspect='equal')
		axs.append(ax1)
	
	##Create color map
	cmap=plt.cm.gist_rainbow
	colors=cmap(np.linspace(0.,1.,len(nums[1])))
	colors_shuffle=[]
	index_shuf=range(len(colors))
	shuffle(index_shuf)
	
	for i in index_shuf:
		colors_shuffle.append(colors[i])
	
	both=zip(*nums)
	tboth=zip(both, Name_list)
	sboth = sorted(tboth)
	Name_list = [p[-1] for p in sboth]
	tboth = [p[0] for p in sboth]
	nums=[list(t) for t in zip(*tboth)]
	for J,ax,t in zip(nums,axs,names):	
		
		
		
		slices=J

		pie_wedge_collection=ax.pie(slices, colors=colors)
		ax.set_title(t)
		for pie_wedge in pie_wedge_collection[0]:
			pie_wedge.set_edgecolor('none')
	ll=ax.legend(Name_list, loc=4, ncol=4,fontsize='10')
	cc+=1
	plt.savefig(outdir+"NB_"+nn+".pdf")
plt.show()




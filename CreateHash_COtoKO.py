#!/usr/bin/env python
import pandas as pd
import urllib2
from bs4 import BeautifulSoup
import re
import matplotlib.pyplot as plt
import cPickle as cpk
import sys

def getROData(RO):
    #get the KO/CO data associated with a reaction
    RO_httpstr='http://www.genome.jp/dbget-bin/www_bget?rn:'
    ROsite=urllib2.urlopen(RO_httpstr+RO)
    Rsoup=BeautifulSoup(ROsite)
    Rtable=Rsoup.table
    KOlist=Rtable.findAll(text=re.compile('K[0-9]{,5}$'))
    KOlist=[str(item) for item in KOlist]
    COlist=Rtable.findAll(text=re.compile('^C[0-9]{,5}$'))
    COlist=[str(item) for item in COlist]

    return KOlist, COlist
    
def getKOfromCO(COid):
    #Function to find and locate pertenant information about the CO (compound) ids
    #Get the website that we want by adding this string to the input
    CO_httpstr='http://www.genome.jp/dbget-bin/www_bget?cpd:'
    CO_httpfull=CO_httpstr+COid
    COsite=urllib2.urlopen(CO_httpfull)
#     Use beautiful soup to parse the website
    Csoup = BeautifulSoup(COsite)
    #Identify the table on the website
    COtable=Csoup.table
    #Find all of the reaction numbers
    if COtable==None:
        ReactionDict=None
        KO_Set=None 
        CO_Set=None
        ModuleNums=None
        pass
    else: 
        ReactionNums=COtable.findAll(text=re.compile('^R[0-9]{,6}$'))
        ReactionNums=[str(item) for item in ReactionNums]

        ModuleNums=COtable.findAll(text=re.compile('^M[0-9]{,6}$'))
        ModuleNums=[str(item) for item in ModuleNums]

        #Open reaction and identify the KO numbers
        ReactionDict={}
        
        for RO in ReactionNums:
            KOlist,COlist=getROData(RO)
            rx={}
            rx['KO']=KOlist
            rx['CO']=COlist
            ReactionDict[RO]=rx
        KO=[v['KO'] for v in ReactionDict.itervalues()]
        CO=[v['CO'] for v in ReactionDict.itervalues()]
        KO_Set=list(set([item for sublist in KO for item in sublist]))
        CO_Set=list(set([item for sublist in CO for item in sublist]))
    return ReactionDict, KO_Set, CO_Set,  ModuleNums




if __name__=='__main__': 
    #load in all of the CO data that Krista was able to identify.
    InputFile=sys.argv[1]
    if len(sys.argv)>2:
        indexCol=sys.argv[2]
    else: 
        indexCol='cNumber'
    CO_RawData=pd.read_csv(InputFile, index_col=indexCol)
    CompoundHash={}
    for i,CO in enumerate(CO_RawData.index):
        if i%100==0:
            print i, CO
            cpk.dump(CompoundHash, open('running_Script.pickle', 'w'))    
        ReactionDict, KO_Set, CO_Set, ModuleNums=getKOfromCO(CO)
        if ReactionDict==None:
            pass
        else:
            outHash={}
            outHash['Reaction']=ReactionDict
            outHash['Related KO']=KO_Set
            outHash['Related CO']=CO_Set
            outHash['Modules']=ModuleNums
            outHash['Name']=CO
            CompoundHash[CO]=outHash
    cpk.dump(CompoundHash, open(sys.argv[1]+'.pickle', 'w'))    
    

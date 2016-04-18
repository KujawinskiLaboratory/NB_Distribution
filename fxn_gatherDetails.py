#make one py file with all the plotting functions and information gathering
#KLongnecker, 13 April 2016, updated 4/15/2016
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import palettable as pal
import glob

from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from IPython.display import Image, HTML
   
def gatherDetails(makeNclusters,trimPath,forRelatedness,folderName,CO_fromMATLAB,KO_Norm2Mean,Insitu_TPM_DIA,Insitu_TPM_DIN,Insitu_TPM_Oth):
    colLabel = ['nCpds','nGenes'] #starting with this is easiest - makes one list, no need to flatten

    for item in range(makeNclusters):
        colLabel.append('Km' + str(item) + '_cpd')
        colLabel.append('Km' + str(item) + '_gene')

    gatherCounts = pd.DataFrame(0, index = trimPath, columns = colLabel)

    #setup the strings to match first
    rnString = re.compile('(?:[rn:R])(\d+)$') #will return R00190
    cpdString = re.compile('(?:[cpd:C])(\d+)$') #will return C00190

    size = 20 #turns out I can increase the size of the compounds in the plots

    for kn in range(makeNclusters):
        fullSet = set(forRelatedness.KEGG)
        oneK = forRelatedness[forRelatedness.kmeans == kn] #get gene & transcript information for one Kmeans group
        getKm = 'Km' + str(kn)

        #check if the directories exist, one for pathway files
        directoryPDF = folderName + str(kn) + '/pathway_files'
        if not os.path.exists(directoryPDF):
            os.makedirs(directoryPDF)

        #check if the directories exist, one for reaction files
        directoryPNG = folderName + str(kn) + '/reaction_files'
        if not os.path.exists(directoryPNG):
            os.makedirs(directoryPNG) 
            
        #check if the directories exist, one for species 
        directorySpecies = folderName + str(kn) + '/species_files'
        if not os.path.exists(directorySpecies):
            os.makedirs(directorySpecies) 
        
        for item in trimPath: #searching within one pathway at a time
            plotPathway = [] #gather up yes/no and will only plot if have linked genes/mtabs    
            genes = getKfrom_ko(item)
            compounds = getCfrom_ko(item)
            gatherCounts.loc[item,'nCpds'] = len(compounds)
            gatherCounts.loc[item,'nGenes'] = len(genes)     
            #have to track genes and compounds differently for the biopython plotting later on 
            setG = set(genes)
            setC = set(compounds)
            setB = set(oneK.KEGG)
            intGenes = setG.intersection(setB)
            intCompounds = setC.intersection(setB)
            gatherCounts.loc[item,(getKm + '_gene')] = len(intGenes)
            gatherCounts.loc[item,(getKm + '_cpd')] = len(intCompounds)
            for gen in intGenes: #go through each gene...one at a time
                rnList = kegg_link('reaction',gen).read() #get the list of reactions for that gene
                #can have cases where there is a gene and no reaction (K02906 for example). This returns rnList = '\n'
                #since this is not actually empty...need a few way to filter those out
                test = '\n'
                if test != rnList:
                    for line in rnList.rstrip().split('\n'):
                        countCpd = []
                        countGene = []
                        m = rnString.search(line) #get the reaction number
                        cpdList = kegg_link('cpd',m.group(0)).read() #now go get the compounds for that reaction
                        #can have no compounds in a reaction (only glycans, begin with G, nothing I have matched)
                        if len(cpdList) > 1: #will be true if cpdList includes compounds
                            for line2 in cpdList.rstrip().split('\n'):
                                m2 = cpdString.search(line2).group(0)
                                #now that I have a compound, check if it is in intCompounds
                                if m2 in intCompounds:
                                    countCpd.append(m2) 
                                    countGene.append(gen)
                                    plotPathway.append('yes')
                        ##Now, plot the PNG files (one for each reaction within a pathway)
                        if len(countCpd) > 0:
                            dayList = ['S1','S2','S3','S4','S5']
                            kData = pd.DataFrame(columns = dayList)
                            for k in set(countGene):
                                kData = kData.append(oneK.ix[k,dayList])
                            cData = pd.DataFrame(columns = dayList)
                            for co in set(countCpd):
                                #convert CO to RI, can have multiple options
                                j = findRInumber(oneK,co)
                                cData = cData.append(oneK.loc[j,dayList])
                            fig,ax = plt.subplots(1)
                            cData.T.plot(color = 'k',ax=ax)
                            kData.T.plot(color = 'r',ax=ax)
                            handles, labels = ax.get_legend_handles_labels()
                            #convert the RI numbers to COnumbers for the figure
                            for ia, a in enumerate(labels):
                                #add compound/gene name to the legend
                                if a[0]== 'R':
                                    tLabel = convertRItoCO(CO_fromMATLAB,a)
                                    fn = kegg_list(tLabel).read()                          
                                    labels[ia] = fn
                                elif a[0] == 'K':
                                    fn = kegg_list(a).read()
                                    labels[ia] = fn
                            ax.legend(handles, labels, bbox_to_anchor = ([-1, 0.5]))
                            fig.suptitle('pathway ' + item + ', Kmeans grp ' + str(kn))
                            pngName = 'pathway' + item + '_' + m.group(0) + '.png'
                            fig.savefig(directoryPNG + '/' + pngName, bbox_inches = 'tight')
                            pngName = None #empty it in case that is where I am having issues
                            plt.close()
            if len(plotPathway)>0:
                ## plot the pathway map for this pathway, get details from KEGG for plotting
                useColors = pal.colorbrewer.qualitative.Set1_4.hex_colors
                useColors.insert(0,'#f7f7f7') ## insert white at beginning
                # order of colors: white, red, blue,green,purple
                sd = 0 #not in dataset
                sk = 1 #in K means group and pathway
                sa = 2 #in pathway, in any K means (for genes, bc overlap in numbers)
                sn = 3 #in pathway, not in K means group (compounds only)               
                su = 4 #unconnected gene or compound
                line1 = useColors[sd] + ', not in dataset' + '\n'
                line2 = useColors[sk] + ', in K means group and pathway' + '\n'
                line3 = useColors[sa] + ', #in pathway, in any K means (for genes, bc overlap in numbers)' +'\n'
                line4 = useColors[sn] +  ', #in pathway, not in K means group (compounds only)' + '\n'               
                line5 = useColors[su] + ', #unconnected gene or compound' + '\n'
                file = open("readme_colorsInPathways.txt", "w")
                file.write(line1 + line2 + line3 + line4 + line5)
                file.close()
                
                pathway = KGML_parser.read(kegg_get(item, "kgml"))
                for element in pathway.orthologs:
                    #print element.name
                    for graphic in element.graphics:
                        tg = element.name[3:9] #skip over the 'ko:'
                        if (tg in intGenes):
                            #in the pathway AND in the set for this particular K means group
                            graphic.bgcolor = useColors[sk] #
                            
                            #if this is something in the pathway, plot up the species for the K number
                            if tg in Insitu_TPM_DIA.index.tolist():
                                Dk=Insitu_TPM_DIA.loc[tg]
                            else: 
                                Dk = 0/Insitu_TPM_DIA.iloc[0] #make an empty frame
                            if tg in Insitu_TPM_DIN.index.tolist():
                                Nk=Insitu_TPM_DIN.loc[tg]
                            else:
                                Nk = 0/Insitu_TPM_DIN.iloc[0]
                            if tg in Insitu_TPM_Oth.index.tolist():
                                Ok=Insitu_TPM_Oth.loc[tg]
                            else:
                                Ok = 0/Insitu_TPM_Oth.iloc[0]
                            fig,ax=plt.subplots(1)
                            ax.stackplot(range(5), Dk, Nk, Ok, colors=pal.colorbrewer.qualitative.Set3_6_r.hex_colors, lw=0)
                            ax.set_xticks(range(5))
                            ax.set_xticklabels([1,2,3,4,5])
                            ax.set_ylabel('In situ TPM')
                            plt.title(tg + ', lt orange=diatoms, blue=dinos, dk orange=other')
                            fig.savefig(directorySpecies + '/' + tg + '_species.png',bbox_inches='tight')
                            plt.close()
                        elif (tg in fullSet) and (tg in genes) and (tg not in intGenes):
                            #in the pathway AND in the set of genes from RI, allow any Kmeans group for genes
                            graphic.bgcolor = useColors[sa] #
                        elif (tg not in fullSet) and (tg in genes) and (tg not in KO_Norm2Mean.index.tolist()):
                            #in the pathway, but *not* in anything from the RI samples
                            graphic.bgcolor = useColors[sd] #
                        elif (tg not in fullSet) and (tg in genes) and (tg in KO_Norm2Mean.index.tolist()): 
                            #an unconnected gene in the RI data
                            graphic.bgcolor = useColors[su] #
                # Change the colours of compounds (mostly same as genes
                for element in pathway.compounds:
                    for graphic in element.graphics:
                        tc = element.name[4:10] #skip over the 'cpd:'
                        if (tc in intCompounds):
                            #in the pathway AND in the set for this particular K means group
                            graphic.bgcolor = useColors[sk] #
                            graphic.width = size
                            graphic.height = size
                        elif (tc in fullSet) and (tc in compounds) and (tc not in intCompounds):
                            #in the pathway AND in the set of compounds from RI, but *not* in this Kmeans group
                            graphic.bgcolor = useColors[sn] #
                            graphic.width = size
                            graphic.height = size
                        elif (tc not in fullSet) and (tc in compounds) and (tc not in CO_fromMATLAB.cNumber.values):
                            #in the pathway, but *not* in anything from the RI samples
                            graphic.bgcolor = useColors[sd] #  
                        elif (tc not in fullSet) and (tc in compounds) and (tc in CO_fromMATLAB.cNumber.values): #seems like a hack
                            #unconnected compound in the RI data
                            graphic.bgcolor = useColors[su] #
                            graphic.width = size
                            graphic.height = size
                canvas = KGMLCanvas(pathway, import_imagemap=True)
                pdfName = 'mapWithColors_' + str(item) + '.pdf'
                canvas.draw(directoryPDF + '/' + pdfName)
                pdfName = None #empty it in case that is where I am having issues
    #stick the pathway information into gatherCounts before I export...
    #want to export gatherCounts, with the added pathway name as a new column
    gatherCounts['pathwayInfo'] = ''
    gatherCounts['pathwayGroup_A'] = ''
    gatherCounts['pathwayGroup_B'] = ''
    gatherCounts['pathwayGroup_C'] = ''
    #go read in the file from KEGG
    D = glob.glob('br08901.keg') #from http://www.genome.jp/kegg-bin/get_htext?br08901.keg; 3/15/2016
    allBRITE=[]
    for idx,nof in enumerate(D):
        allBRITE = ReadBRITEfile(nof) 

    #put the pathway name and group into the data frame before exporting it
    for item in gatherCounts.index:
        #if this error appears: IndexError: index 0 is out of bounds for axis 0 with size 0
        #KEGG has updated a pathway, but not the BRITE file (see below for work around)
        pathstr = kegg_list(item).read()
        #this next line splits the string at the '\t', then keeps the piece at index = 1, and strips off the '\n'
        gatherCounts.loc[item,('pathwayInfo')] = pathstr.split('\t')[1].rstrip()
        t = allBRITE.loc[allBRITE['map']==item[2:]]  
        #put in a check to see if t.empty ...will be empty if KEGG updated pathway and not BRITE file
        if t.empty is False: 
            gatherCounts.set_value(item,'pathwayGroup_A',t['A'].values[0])
            gatherCounts.set_value(item,'pathwayGroup_B',t['B'].values[0])
            gatherCounts.set_value(item,'pathwayGroup_C',t['C'].values[0])
    
    return gatherCounts

#set up a function to get the list of K orthologues for a given pathway (must be defined as ko00140 NOT map00140)
def getKfrom_ko(ko_id):
    pathway_file = kegg_get(ko_id).read()  # query and read the pathway
    K_list = []

    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section
        if current_section == "ORTHOLOGY":
            K_identifiers = line[12:].split("; ")
            t = K_identifiers[0]
            K_id = t[0:6]

            if not K_id in K_list:
                K_list.append(K_id)
    return K_list

#set up a function to get the list of compounds for a given pathway (must be defined as ko00140 NOT map00140)
def getCfrom_ko(ko_id):
    pathway_file = kegg_get(ko_id).read()  # query and read the pathway
    compound_list = []

    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section
        if current_section == "COMPOUND":
            compound_identifiers = line[12:].split("; ")
            t = compound_identifiers[0]
            compound_id = t[0:6]

            if not compound_id in compound_list:
                compound_list.append(compound_id)
    return compound_list

def findRInumber(dataIn,KEGGin):
    #find possible RI numbers for a given KEGG number. 
    dataOut = []
    for i,KEGG in enumerate(dataIn['KEGG']):
        if KEGG == KEGGin:
            t = dataIn.index[i]
            dataOut.append(t)
    return dataOut

def convertRItoCO(dataIn,RIin):
    #do the reverse, given an RInumber find the cNumber
    dataOut = dataIn.loc[RIin].loc['cNumber']
    return dataOut

# A bit of code that will help us display the PDF output
def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

# A bit of helper code to shorten long text
def head(text, lines=10):
    """ Print the first lines lines of the passed text.
    """
    print '\n'.join(text.split('\n')[:lines] + ['[...]'])


#organize pathways into the groups defined in the BRITE file
def ReadBRITEfile(briteFile):
    forBrite = pd.DataFrame(columns = ['map','A','B','C','wholeThing'])
    # set up the expressions to match each level in the BRITE hierarchy
    
    textA = re.compile(r'(^A<b>)(.+)(</b>)\s*(.*)$')
    textB = re.compile(r'(^B)\s*(.*)$')
    textC = re.compile(r'(\d+)\s*(.*)$')
    #this relies on the fact that the rows are in order: A, with B subheadings, then C subheadings
    setA = []
    idxA = []

    setB = []
    setC = []

    with open(briteFile) as f:
        for idx,line in enumerate(f):
            if line[0] is not '#': #skip over the comments
                mA = textA.search(line) 
                mB = textB.search(line) 
                mC = textC.search(line) 
                if mA:
                    setA = mA.group(2)
                    #house cleaning (probably c)
                    idxA = idx
                    forBrite.loc[idx,'A'] = setA
                    forBrite.loc[idx,'wholeThing'] = line #using this as a double check for now
                    #forBrite.loc[idx,'map'] = mC.group(1)
                elif mB:
                    setB = mB.group(2)
                    forBrite.loc[idx,'A'] = setA
                    forBrite.loc[idx,'B'] = setB
                    forBrite.loc[idx,'wholeThing'] = line
                    #forBrite.loc[idx,'map'] = mC.group(1)
                elif mC:
                    #Tracer()()
                    setC = mC.group(2)
                    forBrite.loc[idx,'A'] = setA
                    forBrite.loc[idx,'B'] = setB
                    forBrite.loc[idx,'C'] = setC
                    forBrite.loc[idx,'wholeThing'] = line
                    forBrite.loc[idx,'map'] = mC.group(1)

        return forBrite

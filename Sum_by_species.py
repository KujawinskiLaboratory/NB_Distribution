'''
Created on Sept 11, 2013

Sum_by_species.py

@author: harrietalexander
'''

import sys
import numpy
import re
import csv
if len(sys.argv)<2:
	print "Usage: python FastaHeader_to_Seq_CommandLine.py HTSeq.tab"
	sys.exit(1)
fileIn=sys.argv[1]

handle = csv.reader(open(fileIn), delimiter="\t")

Spe_hash={}
for record in handle: 
	m=re.search('^.*?(?=[0-9])',record[0].strip())
	id=m.group(0)
	if re.search('CCMP', id):
		id=record[0][0:8]
	if re.search('RCC',id):
		id=record[0][0:7] 
	num=record[1:]
	num = [int(i) for i in num]
	if id in Spe_hash:
		num_=Spe_hash[id]
		new=[x+y for x,y in zip(num, num_)]
		Spe_hash[id]=new
	else:
		Spe_hash[id]=num
test = file("SummedSpecies.tab", "w")
for key in Spe_hash:
	test.write(key)
	test.write("\t")
	test.write("\t".join(str(x) for x in Spe_hash[key]))
	test.write('\n')
test.close()
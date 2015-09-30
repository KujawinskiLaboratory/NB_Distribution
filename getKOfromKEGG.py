#want something that I can use to repeatedly get the KO information ...

import urllib2
from bs4 import BeautifulSoup
 
def getKO(KO):
	KO_httpstr='http://www.genome.jp/dbget-bin/www_bget?ko:'
	KOsite=urllib2.urlopen(KO_httpstr+KO)
	Ksoup=BeautifulSoup(KOsite)
	results = Ksoup.find_all("td",{"class":"td40"})
	for i, cell in enumerate(results):
		if i==1:
			t = cell.get_text(strip = True)
			labels[ia] = a + ', ' + str(t)
	
	return labels


                            #for cell in results:
                             #   t = results.getText(strip=True)
                              #  labels[ia] = m + ', ' + str(t)

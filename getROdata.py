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
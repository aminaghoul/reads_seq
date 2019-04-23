import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import beta
from scipy import stats
import random
import time
import re

import sys

"""
monFichier = open("../Documents/MINION_FASTQ/SRR6472704.fastq", "r")
monFichier1 = open("../Documents/MINION_FASTQ/SRR7548026.fastq", "r")
monFichier2 = open("../Documents/MINION_FASTQ/SRR7989207.fastq", "r")
monFichier3 = open("../Documents/MINION_FASTQ/SRR8335319.fastq", "r")
monFichier4 = open("../Documents/MINION_FASTQ/SRR8467877.fastq", "r")
monFichier5 = open("../Documents/MINION_FASTQ/SRR8494915.fastq", "r")
monFichier6 = open("../Documents/MINION_FASTQ/SRR8494916.fastq", "r")
monFichier7 = open("../Documents/MINION_FASTQ/SRR8494917.fastq", "r")
monFichier8 = open("../Documents/MINION_FASTQ/SRR8494918.fastq", "r")
monFichier9 = open("../Documents/MINION_FASTQ/SRR8494919.fastq", "r")
monFichier10 = open("../Documents/MINION_FASTQ/SRR8494921.fastq", "r")
monFichier11 = open("../Documents/MINION_FASTQ/SRR8494924.fastq", "r")
monFichier12 = open("../Documents/MINION_FASTQ/SRR8494938.fastq", "r")
monFichier13 = open("../Documents/MINION_FASTQ/SRR8494942.fastq", "r")
monFichier14 = open("../Documents/MINION_FASTQ/SRR8538952.fastq", "r")
monFichier15 = open("../Documents/MINION_FASTQ/SRR8560598.fastq", "r")
monFichier16 = open("../Documents/MINION_FASTQ/SRR8588895.fastq", "r")
monFichier17 = open("../Documents/MINION_FASTQ/SRR8691535.fastq", "r")
monFichier18 = open("../Documents/MINION_FASTQ/SRR8848723.fastq", "r")

Fichiers = [monFichier18]

id = []
nbr = []
tab = []

def b(x,g,n):
    return (1/g)*n*(1-x)**(n-1)




for fichier in Fichiers:
    for line in fichier:
        fields = line.strip().split()
        for word in fields:
            if word.startswith("length="):
                id.append(word.strip("length="))
    id = list(map(int, id))
    tmp = 0
    for i in id:
        if i != tmp:
            tab.append(i)
        tmp = i
    nbr.append(len(tab))

    plt.hist(tab, range=(0, 80000))

    plt.xlabel('taille en bp')
    plt.title('Distribution de la taille des reads-MINION')
    plt.savefig('TailleReadsMinion.png')
    plt.show()
    x = np.linspace(0, 1)
    plt.plot(x, b(x, sum(tab), len(tab)))
    plt.title('Distribution beta 1/g*(n)*(1-x)^(n-1) g : taille génome, n=nombre de fragments')
    plt.savefig('beta.png')
    plt.show()
    N = np.random.uniform(0, sum(tab)/1000, 1000)
    t = sum(tab) / N
    plt.hist(t, range=(0, 10000))
    plt.xlabel('taille de N fragment où N suit loi unif')
    plt.title('DIstribution de la taille de N fragment')
    plt.savefig('tailleNfragment.png')
    plt.show()
    print(sum(tab))
    print(len(tab))
    id = []
    tab = []


#plt.hist(nbr)
#plt.show()

err1=open("../Documents/MINION_FASTQ/ERR1111170_1.fastq", "r")
err2=open("../Documents/MINION_FASTQ/ERR1111170_2.fastq", "r")
err3=open("../Documents/MINION_FASTQ/ERR1111171_1.fastq", "r")
err4=open("../Documents/MINION_FASTQ/ERR1111171_2.fastq", "r")

err=[err1, err2, err3, err4]


fich1=open("../Documents/MINION_FASTQ/DRR129261.fastq", "r")
fich2=open("../Documents/MINION_FASTQ/DRR129651.fastq", "r")
fich3=open("../Documents/MINION_FASTQ/DRR129656.fastq", "r")
fich4=open("../Documents/MINION_FASTQ/DRR129658.fastq", "r")
fich5=open("../Documents/MINION_FASTQ/DRR129660.fastq", "r")
fich6=open("../Documents/MINION_FASTQ/DRR129662.fastq", "r")
fich7=open("../Documents/MINION_FASTQ/DRR154113.fastq", "r")
fich8=open("../Documents/MINION_FASTQ/DRR154114.fastq", "r")
fich9=open("../Documents/MINION_FASTQ/DRR154115.fastq", "r")
fich10=open("../Documents/MINION_FASTQ/DRR154116.fastq", "r")
fich11=open("../Documents/MINION_FASTQ/DRR155408.fastq", "r")
fich12=open("../Documents/MINION_FASTQ/DRR155409.fastq", "r")
fich13=open("../Documents/MINION_FASTQ/DRR161061.fastq", "r")
fich14=open("../Documents/MINION_FASTQ/DRR164908.fastq", "r")
fich15=open("../Documents/MINION_FASTQ/DRR164909.fastq", "r")
fich16=open("../Documents/MINION_FASTQ/DRR164910.fastq", "r")
fich17=open("../Documents/MINION_FASTQ/DRR164911.fastq", "r")
fich18=open("../Documents/MINION_FASTQ/DRR164912.fastq", "r")
fich19=open("../Documents/MINION_FASTQ/DRR164913.fastq", "r")
fich20=open("../Documents/MINION_FASTQ/DRR164914.fastq", "r")
fich21=open("../Documents/MINION_FASTQ/DRR164915.fastq", "r")

drr = [fich1, fich2, fich3, fich4, fich5, fich6, fich7, fich8, fich9, fich10, fich11, fich12, fich13, fich14, fich15,
       fich16, fich17, fich18, fich19, fich20, fich21]

"""
def afficher_hist(drr):
    #x = np.linspace(0, 1)
    taille = []
    f = []
    for fich in drr:
        for line in fich:
            fields = line.strip().split()
            for idx, word in enumerate(fields):
                f.append(word)

        for i, element in enumerate(f):
            if element.endswith("/1") and i + 1 < len(f):
                taille.append(len(f[i+1]))

        fig, ax = plt.subplots(1, 1)
        ax.hist(taille, range=(0, 10000))
        #plt.savefig('TailleReadsMinionDRR.png')
        plt.show()
        x=np.linspace(0, 1)
        #x = np.linspace(beta.ppf(0.01, 1, len(taille)), beta.ppf(0.99, 1, len(taille)), 100000)
        #plt.plot(x, b(x, sum(taille), len(taille)))
        #plt.show()
        #print(stats.kstest(taille, 'beta', args=(1, len(taille))))
        print(taille)
        print(len(taille))
        #N = np.random.uniform(0, sum(taille) / 1000, 1000)
        #t = sum(taille) / N
        #plt.hist(t, range=(0, 10000))
        #plt.xlabel('taille de N fragment où N suit loi unif')
        #plt.title('DIstribution de la taille de N fragment')
        #plt.savefig('tailleNfragment.png')
        #beta.fit(taille, 1, len(taille), scale=10000)

        size = max(taille)
        x = np.linspace(0, max(taille), len(taille))
        # creating the dummy sample (using beta distribution)

        # creating the histogram
        h = plt.hist(taille, range=(0, 10000))


        dist = getattr(scipy.stats, 'beta')
        param = dist.fit(taille)
        pdf_fitted = dist.pdf(x, 1, len(taille), loc=0, scale=1) * max(taille)
        plt.plot(pdf_fitted, label='beta')
        plt.xlim(0, 10000)
        plt.show()





        plt.show()



        #nbr.append(len(taille))
        f = []
        taille = []



pacbio1=open("../Documents/PACBIO/DRR008620_subreads.fastq", "r")
pacbio2=open("../Documents/PACBIO/DRR013332_subreads.fastq", "r")
pacbio3=open("../Documents/PACBIO/DRR013333_subreads.fastq", "r")
pacbio4=open("../Documents/PACBIO/DRR013334_subreads.fastq", "r")
pacbio5=open("../Documents/PACBIO/DRR013335_subreads.fastq", "r")
pacbio6=open("../Documents/PACBIO/DRR013336_subreads.fastq", "r")
pacbio7=open("../Documents/PACBIO/DRR013337_subreads.fastq", "r")
pacbio8=open("../Documents/PACBIO/DRR013338_subreads.fastq", "r")
pacbio9=open("../Documents/PACBIO/DRR013339_subreads.fastq", "r")
pacbio10=open("../Documents/PACBIO/DRR013340_subreads.fastq", "r")
pacbio11=open("../Documents/PACBIO/DRR013341_subreads.fastq", "r")
pacbio12=open("../Documents/PACBIO/DRR013342_subreads.fastq", "r")
pacbio13=open("../Documents/PACBIO/DRR013343_subreads.fastq", "r")
pacbio14=open("../Documents/PACBIO/DRR013344_subreads.fastq", "r")
pacbio15=open("../Documents/PACBIO/DRR013345_subreads.fastq", "r")
pacbio16=open("../Documents/PACBIO/DRR015080_subreads.fastq", "r")
pacbio17=open("../Documents/PACBIO/DRR015123_subreads.fastq", "r")
pacbio18=open("../Documents/PACBIO/DRR016446_subreads.fastq", "r")
pacbio19=open("../Documents/PACBIO/DRR016447_subreads.fastq", "r")
pacbio20=open("../Documents/PACBIO/DRR017980_subreads.fastq", "r")
pacbio21=open("../Documents/PACBIO/DRR018508_subreads.fastq", "r")
pacbio22=open("../Documents/PACBIO/DRR018509_subreads.fastq", "r")
pacbio23=open("../Documents/PACBIO/DRR028956_subreads.fastq", "r")
pacbio24=open("../Documents/PACBIO/DRR074681_subreads.fastq", "r")
pacbio25=open("../Documents/PACBIO/ERR070563_1.fastq", "r")
pacbio26=open("../Documents/PACBIO/ERR070564_1.fastq", "r")
pacbio27=open("../Documents/PACBIO/ERR070565_1.fastq", "r")
pacbio28=open("../Documents/PACBIO/ERR070566_1.fastq", "r")
pacbio29=open("../Documents/PACBIO/ERR070567_1.fastq", "r")


pacbio = [pacbio1, pacbio2, pacbio3, pacbio4, pacbio5, pacbio6, pacbio7, pacbio8, pacbio9, pacbio10, pacbio11, pacbio12,
          pacbio13, pacbio14, pacbio15, pacbio16, pacbio17, pacbio18, pacbio19, pacbio20, pacbio21, pacbio22, pacbio23,
          pacbio24, pacbio25, pacbio26, pacbio27, pacbio28, pacbio29]
#x = np.linspace(0, 1)
#plt.plot(x, b(x, 11110, 114))
#plt.show()
#afficher_hist(drr)


afficher_hist(pacbio)

"""
print(nbr)
plt.hist(nbr)
plt.xlabel('nombres de fragment')
plt.title('Distribution du nombre de fragment pour 40 génomes')
plt.savefig('DistributionN40genomes.png')
plt.show()

afficher_hist(err)
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import random
import time
import re

import sys


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

Fichiers = [monFichier ,monFichier1 ,monFichier2 ,monFichier3 ,monFichier4 ,monFichier5,
monFichier6 ,monFichier7 ,monFichier8,
monFichier9,monFichier10,monFichier11 ,monFichier12,monFichier13,monFichier14 ,monFichier15 ,monFichier16 ,monFichier17 ,
monFichier18 ]

id = []
nbr=[]
tab=[]


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
    print((tab))
    plt.hist(tab, range=(0, 100000), bins=20)
    plt.savefig('image.png')
    plt.show()
    id=[]
    tab=[]


plt.hist(nbr)
plt.show()


taille=[]
fich1=open("../Documents/MINION_FASTQ/DRR129261.fastq", "r")

for line in fich1:
    fields = line.strip().split()
    for idx, word in enumerate(fields):
        if word.endswith("/1"):
            taille.append(fields[idx+1])


print(taille)





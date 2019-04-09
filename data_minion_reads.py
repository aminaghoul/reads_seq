import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import random
import time
import re





monFichier = open("../Documents/longueur_readSRR7548026.txt", "r")
monFichier1 = open("../Documents/longueur_readSRR6472704.txt", "r")
monFichier2 = open("../Documents/longueur_readSRR7989207.txt", "r")
monFichier3 = open("../Documents/longueur_read3.txt", "r")


Fichiers=[monFichier, monFichier1, monFichier2, monFichier3]

z=[]

for fichier in Fichiers:
    for i in fichier:
        z.append(i.strip())

    plt.hist(z)
    print(len(z))
    plt.show()
    z=[]



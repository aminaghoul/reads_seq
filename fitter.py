import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import random
import time
import re




def b(x, g, n):
    """
    fonction loi beta qui correspond à la distribution de la taille des reads
    :param x: taille des reads
    :param g: taille du génome
    :param n: nombre de fragments
    :return: fonction loi beta
    """
    return (1 / g) * n * (1 - x) ** (n - 1)


monFichier = open("../Documents/Donnees/SRR6472704.fastq", "r")
id = []
taille = []
norm = []
x=np.linspace(0, 1)
tmp = 0


for line in monFichier:
    fields = line.strip().split()
    for word in fields:
        if word.startswith("length="):
            id.append(word.strip("length="))
id = list(map(int, id))

for i in id:
    if i != tmp:
        taille.append(i)
    tmp = i
plt.hist(taille, normed=True)
plt.xlabel('taille en bp')
plt.title('Distribution de la taille des reads-MINION')
plt.show()

for i in taille:
    norm.append(i / max(taille))

weights = np.ones_like(norm) / float(len(norm))
plt.hist(norm, weights=weights)
plt.plot(x, b(x, 1000000, len(taille)))
plt.show()
print(len(taille))
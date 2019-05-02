import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fmin
from scipy.stats import beta
from scipy.special import gamma as gammaf
from scipy import stats
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
n = len(taille)

xx = np.linspace(0, max(taille))

plt.hist(taille)
plt.plot(xx, beta.pdf(xx, 1, len(taille)))
plt.xlabel('taille en bp')
plt.title('Distribution de la taille des reads-MINION')
plt.show()
print(stats.kstest(taille, 'beta', args=(1, len(taille))))

plt.plot(xx, (1/1000000)*n * (1 - xx) ** (n - 1))

plt.show()


for i in taille:
    norm.append(i / max(taille))

def betaNLL(param,*args):
    """
    Negative log likelihood function for beta
    <param>: list for parameters to be fitted.
    <args>: 1-element array containing the sample data.

    Return <nll>: negative log-likelihood to be minimized.
    """

    a, b = param
    data = args[0]
    pdf = beta.pdf(data,a,b,loc=0,scale=1)
    lg = np.log(pdf)
    mask = np.isfinite(lg)
    nll = -lg[mask].sum()
    return nll

mean=np.mean(norm)
var=np.var(norm,ddof=1)
alpha1=mean**2*(1-mean)/var-mean
beta1=alpha1*(1-mean)/mean

#------------------Fit using mle------------------
result=fmin(betaNLL,[1,1],args=(norm,))
alpha2,beta2=result

#----------------Fit using beta.fit----------------
alpha3,beta3,xx,yy=beta.fit(norm)

print('\n# alpha,beta from moments:',alpha1,beta1)
print('# alpha,beta from mle:',alpha2,beta2)
print('# alpha,beta from beta.fit:',alpha3,beta3)

#-----------------------Plot-----------------------
weights = np.ones_like(norm) / float(len(norm))
plt.hist(norm, weights=weights)
fitted=lambda x,a,b:gammaf(a+b)/gammaf(a)/gammaf(b)*x**(a-1)*(1-x)**(b-1) #pdf of beta

xx=np.linspace(0,max(norm))

plt.plot(xx,fitted(xx,alpha1,beta1),'g')
plt.plot(xx,fitted(xx,alpha2,beta2),'b')
plt.plot(xx,fitted(xx,alpha3,beta3),'r')

plt.show()


weights = np.ones_like(norm) / float(len(norm))
plt.hist(norm, weights=weights)
plt.plot(x, b(x, 1000000, len(taille)))
plt.show()
print(len(taille))
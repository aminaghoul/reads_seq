import matplotlib
import matplotlib.pyplot as plt
import numpy
import numpy as np
from scipy.optimize import fmin
from scipy.stats import beta
from scipy.special import gamma as gammaf
from scipy import stats
import random
import time
import re

from data_minion_reads import lire_tailles, lire_taille


def b(x, g, n):
	"""
	fonction loi beta qui correspond à la distribution de la taille des reads
	:param x: taille des reads
	:param g: taille du génome
	:param n: nombre de fragments
	:return: fonction loi beta
	"""
	print(locals())
	return (1 / g) * n * (1 - x) ** (n - 1)


'''monFichier = open("../Documents/Donnees/SRR6472704.fastq", "r")
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


for i in tailles:
	if i != tmp:
		taille.append(i)
	tmp = i

			'''
# id = list(map(int, id))

taille = lire_taille()

nombre_fragments = len(taille)

axe_abscisse = np.linspace(0, max(taille))

# plt.hist(taille)
# plt.show()
"""
fig, plots = plt.subplots(1, 2, )

axe_abscisse = np.linspace(0, 1, 10) # len(set(taille))
plots[0].plot(axe_abscisse, beta.pdf(axe_abscisse, 1, len(axe_abscisse)))
print(beta.pdf(axe_abscisse, 1, len(axe_abscisse)))
plt.show()


plots[1].set_label('taille en bp')
plots[1].set_title('Distribution de la taille des reads-MINION')
print(stats.kstest(taille, 'beta', args=(1, len(taille))))
plots[1].plot(np.linspace(0, 1, len(set(taille))), b(np.linspace(0, 1, len(set(taille))), 3000000, nombre_fragments, ))
print(b(np.linspace(0, 1, len(set(taille))), 3000000, nombre_fragments, ))
# plt.plot(axe_abscisse, (1 / 1000000) * nombre_fragments * (1 - axe_abscisse) ** (nombre_fragments - 1))"""

"""plt.xlabel('taille en bp')
plt.title('Distribution de la taille des reads-MINION')
print(stats.kstest(taille, 'beta', args=(1, len(taille))))
plt.plot(np.linspace(0, 1, len(set(taille))), b(np.linspace(0, 1, len(set(taille))), 3000000, nombre_fragments, ))
print(b(np.linspace(0, 1, len(set(taille))), 3000000, nombre_fragments, ))
# plt.plot(axe_abscisse, (1 / 1000000) * nombre_fragments * (1 - axe_abscisse) ** (nombre_fragments - 1))
 plt.show()"""

t = max(taille)
norm = [i / t for i in taille]


def betaNLL(param, *args):
	"""
	Negative log likelihood function for beta
	<param>: list for parameters to be fitted.
	<args>: 1-element array containing the sample data.

	Return <nll>: negative log-likelihood to be minimized.
	"""

	a, b = param
	data = args[0]
	pdf = beta.pdf(data, a, b, loc=0, scale=1)
	lg = numpy.log(pdf)
	# -----Replace -inf with 0s------
	lg = numpy.where(lg == -numpy.inf, 0, lg)
	print(lg)
	nll = -1 * numpy.sum(lg)
	return nll


"""mean = np.mean(norm)
var = np.var(norm, ddof=1)
alpha1 = mean ** 2 * (1 - mean) / var - mean
beta1 = alpha1 * (1 - mean) / mean"""

fitted = lambda x, a, b: gammaf(a + b) / gammaf(a) / gammaf(b) * x ** (a - 1) * (1 - x) ** (b - 1)  # pdf of beta

# ------------------Fit using mle------------------
# result = fmin(betaNLL, [1, 1], args=(norm,))
# alpha2, beta2 = result

# ----------------Fit using beta.fit----------------
s = sum(taille)
norm = [i / s for i in taille]
alpha3, beta3, axe_abscisse, yy = beta.fit(norm)
plt.plot(axe_abscisse, fitted(axe_abscisse, alpha3, beta3), 'r')
plt.show()

# print('\n# alpha,beta from moments:', alpha1, beta1)
# print('# alpha,beta from mle:', alpha2, beta2)
print('# alpha,beta from beta.fit:', alpha3, beta3)

# -----------------------Plot-----------------------
weights = np.ones_like(norm) / float(len(norm))
plt.hist(norm, weights=weights)

axe_abscisse = np.linspace(0, max(norm))

plt.plot(axe_abscisse, fitted(axe_abscisse, alpha1, beta1), 'g')
# plt.plot(axe_abscisse, fitted(axe_abscisse, alpha2, beta2), 'b')


weights = np.ones_like(norm) / float(len(norm))
plt.hist(norm, weights=weights)
plt.plot(x, b(x, 1000000, len(taille)))
plt.show()
print(len(taille))
if __name__ == '__main__':
	pass

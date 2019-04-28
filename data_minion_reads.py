import json
import os
import typing

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
importer les fichiers et les mettre dans un tableau
"""

'''monFichier = open("../Documents/MINION_FASTQ/SRR6472704.fastq", "r")
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

err1 = open("../Documents/MINION_FASTQ/ERR1111170_1.fastq", "r")
err2 = open("../Documents/MINION_FASTQ/ERR1111170_2.fastq", "r")
err3 = open("../Documents/MINION_FASTQ/ERR1111171_1.fastq", "r")
err4 = open("../Documents/MINION_FASTQ/ERR1111171_2.fastq", "r")

fich1 = open("../Documents/MINION_FASTQ/DRR129261.fastq", "r")
fich2 = open("../Documents/MINION_FASTQ/DRR129651.fastq", "r")
fich3 = open("../Documents/MINION_FASTQ/DRR129656.fastq", "r")
fich4 = open("../Documents/MINION_FASTQ/DRR129658.fastq", "r")
fich5 = open("../Documents/MINION_FASTQ/DRR129660.fastq", "r")
fich6 = open("../Documents/MINION_FASTQ/DRR129662.fastq", "r")
fich7 = open("../Documents/MINION_FASTQ/DRR154113.fastq", "r")
fich8 = open("../Documents/MINION_FASTQ/DRR154114.fastq", "r")
fich9 = open("../Documents/MINION_FASTQ/DRR154115.fastq", "r")
fich10 = open("../Documents/MINION_FASTQ/DRR154116.fastq", "r")
fich11 = open("../Documents/MINION_FASTQ/DRR155408.fastq", "r")
fich12 = open("../Documents/MINION_FASTQ/DRR155409.fastq", "r")
fich13 = open("../Documents/MINION_FASTQ/DRR161061.fastq", "r")
fich14 = open("../Documents/MINION_FASTQ/DRR164908.fastq", "r")
fich15 = open("../Documents/MINION_FASTQ/DRR164909.fastq", "r")
fich16 = open("../Documents/MINION_FASTQ/DRR164910.fastq", "r")
fich17 = open("../Documents/MINION_FASTQ/DRR164911.fastq", "r")
fich18 = open("../Documents/MINION_FASTQ/DRR164912.fastq", "r")
fich19 = open("../Documents/MINION_FASTQ/DRR164913.fastq", "r")
fich20 = open("../Documents/MINION_FASTQ/DRR164914.fastq", "r")
fich21 = open("../Documents/MINION_FASTQ/DRR164915.fastq", "r")

pacbio1 = open("../Documents/PACBIO/DRR008620_subreads.fastq", "r")
pacbio2 = open("../Documents/PACBIO/DRR013332_subreads.fastq", "r")
pacbio3 = open("../Documents/PACBIO/DRR013333_subreads.fastq", "r")
pacbio4 = open("../Documents/PACBIO/DRR013334_subreads.fastq", "r")
pacbio5 = open("../Documents/PACBIO/DRR013335_subreads.fastq", "r")
pacbio6 = open("../Documents/PACBIO/DRR013336_subreads.fastq", "r")
pacbio7 = open("../Documents/PACBIO/DRR013337_subreads.fastq", "r")
pacbio8 = open("../Documents/PACBIO/DRR013338_subreads.fastq", "r")
pacbio9 = open("../Documents/PACBIO/DRR013339_subreads.fastq", "r")
pacbio10 = open("../Documents/PACBIO/DRR013340_subreads.fastq", "r")
pacbio11 = open("../Documents/PACBIO/DRR013341_subreads.fastq", "r")
pacbio12 = open("../Documents/PACBIO/DRR013342_subreads.fastq", "r")
pacbio13 = open("../Documents/PACBIO/DRR013343_subreads.fastq", "r")
pacbio14 = open("../Documents/PACBIO/DRR013344_subreads.fastq", "r")
pacbio15 = open("../Documents/PACBIO/DRR013345_subreads.fastq", "r")
pacbio16 = open("../Documents/PACBIO/DRR015080_subreads.fastq", "r")
pacbio17 = open("../Documents/PACBIO/DRR015123_subreads.fastq", "r")
pacbio18 = open("../Documents/PACBIO/DRR016446_subreads.fastq", "r")
pacbio19 = open("../Documents/PACBIO/DRR016447_subreads.fastq", "r")
pacbio20 = open("../Documents/PACBIO/DRR017980_subreads.fastq", "r")
pacbio21 = open("../Documents/PACBIO/DRR018508_subreads.fastq", "r")
pacbio22 = open("../Documents/PACBIO/DRR018509_subreads.fastq", "r")
pacbio23 = open("../Documents/PACBIO/DRR028956_subreads.fastq", "r")
pacbio24 = open("../Documents/PACBIO/DRR074681_subreads.fastq", "r")
pacbio25 = open("../Documents/PACBIO/ERR070563_1.fastq", "r")
pacbio26 = open("../Documents/PACBIO/ERR070564_1.fastq", "r")
pacbio27 = open("../Documents/PACBIO/ERR070565_1.fastq", "r")
pacbio28 = open("../Documents/PACBIO/ERR070566_1.fastq", "r")
pacbio29 = open("../Documents/PACBIO/ERR070567_1.fastq", "r")

pacbio = [
	pacbio1, pacbio2, pacbio3, pacbio4, pacbio5, pacbio6, pacbio7, pacbio8, pacbio9, pacbio10, pacbio11, pacbio12,
	pacbio13, pacbio14, pacbio15, pacbio16, pacbio17, pacbio18, pacbio19, pacbio20, pacbio21, pacbio22, pacbio23,
	pacbio24, pacbio25, pacbio26, pacbio27, pacbio28, pacbio29]

drr = [
	fich1, fich2, fich3, fich4, fich5, fich6, fich7, fich8, fich9, fich10, fich11, fich12, fich13, fich14, fich15,
	fich16, fich17, fich18, fich19, fich20, fich21]

err = [err1, err2, err3, err4]

srr = [
	monFichier, monFichier1, monFichier2, monFichier3, monFichier4, monFichier5, monFichier6, monFichier7,
	monFichier8, monFichier9, monFichier10, monFichier11, monFichier12, monFichier13, monFichier14,
	monFichier15, monFichier16, monFichier17, monFichier18]
'''


def b(x, g, n):
    """
    fonction loi beta qui correspond à la distribution de la taille des reads
    :param x: taille des reads
    :param g: taille du génome
    :param n: nombre de fragments
    :return: fonction loi beta
    """
    return (1 / g) * n * (1 - x) ** (n - 1)


def afficherhistdrr(srr):
    """
    affiche l'histogramme de la distribution des tailles de reads pour chaque fichier srr
    affiche l'histogramme normalisé avec la courbe b correspondante
    affiche l'histogramme des tailles de reads avec N : le nombre de fragments suivant une loi connue (uniforme)
    :param srr: liste de plusieurs fichiers SRR
    :return: nbr : liste du nombre de fragments pour chaque génome
    """

    id = []
    nbr = []
    taille = []
    norm = []

    tmp = 0

    for fichier in srr:
        for line in fichier:
            fields = line.strip().split()
            for word in fields:
                if word.startswith("length="):
                    id.append(word.strip("length="))
        id = list(map(int, id))

        for i in id:
            if i != tmp:
                taille.append(i)
            tmp = i

        m = max(taille)
        s = sum(taille)
        n = len(taille)

        nbr.append(n)

        plt.hist(taille, normed=True)
        plt.xlabel('taille en bp')
        plt.title('Distribution de la taille des reads-MINION')
        plt.show()

        for i in taille:
            norm.append(i / m)

        weights = np.ones_like(norm) / float(len(norm))
        plt.hist(norm, weights=weights)
        x = np.linspace(0, 1)
        plt.plot(x, b(x, sum(norm), n))
        plt.title('Distribution beta 1/g*(n)*(1-x)^(n-1) g : taille génome, n=nombre de fragments')
        plt.show()

        weights = np.ones_like(norm) / float(len(norm))
        plt.hist(norm, weights=weights)
        dist = getattr(scipy.stats, 'beta')
        param = dist.fit(norm)
        x = scipy.arange(100)
        pdf_fitted = dist.pdf(x, param[0], param[1], loc=param[-2], scale=param[-1])
        plt.plot(pdf_fitted, label='beta')
        plt.xlim(0, 1)
        plt.show()

        print(stats.kstest(norm, 'beta', args=(1, len(norm))))

        N = np.random.uniform(0, sum(taille) / 1000, 1000)
        t = sum(taille) / N
        plt.hist(t, normed=True)
        plt.xlabel('taille de N fragment où N suit loi unif')
        plt.title('DIstribution de la taille de N fragment')
        plt.show()

        id = []
        taille = []


def square(length):
    x, y = length, length
    while x * y > length:
        x -= 1
        if x * y < length:
            return x + 1, y

        y -= 1
        if x * y < length:
            return x, y + 1
    return x, y


def afficher_hist(drr, normalise=False, inclure_theorique=True):
    """
    affiche l'histogramme de la distribution des tailles de reads pour chaque fichier drr
    affiche l'histogramme normalisé avec la courbe b correspondante
    affiche l'histogramme des tailles de reads avec N : le nombre de fragments suivant une loi connue (uniforme)
    :param drr: liste de plusieurs fichiers MINION ou pacbio
    :param nbr: liste du nombre de fragments pour chaque génome
    :return: nbr: liste du nombre de fragments pour chaque génome complétée
    """
    sq = square(len(drr))
    fig, _plots = plt.subplots(*sq)
    plots = []
    for x in _plots:
        plots.extend(x)
    tailles = [(extraire_fichier(path), path) for path in drr]
    for plot, (i, path) in zip(plots, sorted(tailles, key=lambda x: x[0][-1])):
        t, m, s, n = i

        norm = [i / s for i in t]
        if inclure_theorique:
            plot.plot(norm if normalise else t, [b(x, s, n) * s for x in norm])

        if not normalise:
            plot.hist(t)  # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)
            plot.set_xlabel('taille en bp')
            plot.set_title('Distribution de %r' % path)
        else:
            plot.hist(norm)  # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)
            plot.set_xlabel('taille normalisée en bp')
            plot.set_title('Distribution normalisée de %r' % path)
    plt.show()


def human(size: int) -> str:
    suffixes, index = ['o', 'Ko', 'Mo', 'Go', 'To'], 0
    while size > 1024 and index < len(suffixes):
        index += 1
        size /= 1024
    return '%.2f%s' % (size, suffixes[index])


def extraire_fichier(chemin: str) -> typing.Tuple[typing.List[int], int, int, int]:
    """

    :param chemin: fichier MINION ou pacbio
    :return tailles, m, s, n: Une liste de tailles de fragments, le max, la somme et le nombre de fragments
    """
    try:
        with open(chemin + '_length.json', 'rb') as f:
            taille = json.load(f)
    except FileNotFoundError:
        print('Chargement du fichier %r...' % os.path.relpath(chemin), end='', flush=True)

        if chemin.endswith('.fastq'):

            taille, f, index = [], open(chemin), 0
            taille_fichier = os.path.getsize(chemin)
            _t, _taille = time.time(), 0

            ligne, index = f.readline(), 0
            while ligne:
                index += len(ligne)
                match = re.match(pattern=r'\s*(A|T|C|G)+\s*', string=ligne)
                if match is not None:
                    taille.append(len(match.group().strip()))
                ligne, index = f.readline(), index + 1
                if _t + 0.3 < time.time():
                    txt = '%s/%s (%.2f%%)' % (human(index), human(taille_fichier), (100 * index / taille_fichier))
                    padding = ' ' * max(0, _taille - len(txt)) + '\b' * max(0, _taille - len(txt))
                    print('\b' * _taille + txt + padding, end='', flush=True)
                    _taille = len(txt)
                    _t = time.time()

            dumps = json.dumps(taille)
            with open(chemin + '_length.json', 'w') as f:
                f.write(dumps)
        else:
            print('erreur')

            raise NotImplementedError('Fichier %r non supporté' % chemin)
        txt = '%s/%s (%.2f%%)' % (human(taille_fichier), human(taille_fichier), 100.00)
        padding = ' ' * max(0, _taille - len(txt)) + '\b' * max(0, _taille - len(txt))
        print('\b' * _taille + txt + padding, end='', flush=True)

    # IMPORTANT: il faut pouvoir avoir des x croissants (pour comparer avec 𝜷)
    return list(sorted(taille)), max(taille), sum(taille), len(taille)


def extraire_fichiers(chemins: typing.List[str]) -> typing.List[typing.Tuple[typing.List[int], int, int, int]]:
    return [extraire_fichier(chemin) for chemin in chemins]


def extraire_tailles(chemins: typing.List[str]) -> typing.List[int]:
    return [len(i[0]) for i in extraire_fichiers(chemins)]


def main():
    # Chemin du dossier contenant les fichiers à traiter
    origine = os.path.abspath(os.path.join(__file__, '..', '..', 'Données'))

    # Chemins des fichiers à traiter (pour l'instant uniquement les fichiers fastq sont supportés
    fichiers = [os.path.abspath(i.path) for i in os.scandir(origine) if i.path.endswith('.fastq')]

    # On ajoute les fichiers déjà extrait
    suf = '_length.json'
    fichiers.extend([os.path.abspath(i.path.replace(suf, '')) for i in os.scandir(origine) if i.path.endswith(suf)])
    fichiers = sorted(set(fichiers))

    # Extraire la liste des nombres de fragments
    nombres_fragments = extraire_tailles(fichiers)
    print(nombres_fragments)

    # Afficher les tailles, sous forme d'un histogramme
    fig, plot = plt.subplots(1, 1)
    plot.hist(nombres_fragments)  # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)

    # Affiche les tableaux
    afficher_hist(fichiers)

    return 0


if __name__ == '__main__':
    exit(main())

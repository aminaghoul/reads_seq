import json
import os
import random
import re
import time
import typing

import matplotlib.pyplot as plt
from numpy.matlib import randn
from tensorflow.python.ops.distributions import normal
import numpy as np

"""
importer les fichiers et les mettre dans un tableau
"""


def b(x, g, n):
    """
    fonction loi beta qui correspond à la distribution de la taille des reads
    :param x: taille des reads
    :param g: taille du génome
    :param n: nombre de fragments
    :return: fonction loi beta
    """
    return (1 / g) * n * (1 - x) ** (n - 1)


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


def afficher_hist(normalise=False, inclure_theorique=True):
    """
    affiche l'histogramme de la distribution des tailles de reads pour chaque fichier drr
    affiche l'histogramme normalisé avec la courbe b correspondante
    affiche l'histogramme des tailles de reads avec N : le nombre de fragments suivant une loi connue (uniforme)
    :return: nbr: liste du nombre de fragments pour chaque génome complétée
    """
    tailles = sorted(lire_tailles().items(), key=lambda y: y[1][0])[-9:]
    sq = square(len(tailles))
    fig, _plots = plt.subplots(*sq)
    plots = []
    for x in _plots:
        plots.extend(x)
    for plot, (path, taille) in zip(plots, tailles):
        s = sum(taille)
        n = len(taille)

        norm = [i / s for i in taille]
        if inclure_theorique:
            plot.plot(norm if normalise else taille, [b(x, s, n) * s for x in norm])

        if not normalise:
            plot.hist(taille)  # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)
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


cache = os.path.abspath(os.path.join(__file__, '..', 'cache'))
os.makedirs(cache, exist_ok=True)
os.makedirs(os.path.join(cache, 'tailles'), exist_ok=True)
cached = {}


def scanner_fichier(chemin):
    """Scanne un fichier et retourne les tailles"""
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

    else:
        print('erreur')

        raise NotImplementedError('Fichier %r non supporté' % chemin)
    txt = '%s/%s (%.2f%%)' % (human(taille_fichier), human(taille_fichier), 100.00)
    padding = ' ' * max(0, _taille - len(txt)) + '\b' * max(0, _taille - len(txt))
    print('\b' * _taille + txt + padding, flush=True)
    return list(sorted(taille))


def extraire_fichier(chemin):
    parent, nom = os.path.split(chemin)
    # On essaie de charger la valeur depuis le cache
    try:
        cached[nom]
    except KeyError:
        try:
            cached[nom] = json.load(open(os.path.join(cache, 'tailles', nom + '.json')))
        except (FileNotFoundError, json.decoder.JSONDecodeError):
            tailles = scanner_fichier(chemin)
            # On convertis la liste en une chaine de caractères, pour
            #   pouvoir les écrires dans un fichier
            dumps = json.dumps(tailles)
            with open(os.path.join(cache, 'tailles', nom + '.json'), 'w') as f:
                f.write(dumps)
            cached[nom] = tailles
    return cached[nom]


def recherche_fichiers(*origines):
    """Scanne les dossier à la recherche de données, rajoute les tailles dans le cache"""
    fichiers = []
    for origine in origines:
        for element in os.scandir(origine):
            if element.path.endswith('.fastq'):
                fichiers.append(os.path.abspath(element.path))
    fichiers = list(set(sorted(fichiers)))
    for chemin in fichiers:
        extraire_fichier(chemin)


def lire_tailles() -> typing.Dict[str, typing.List[int]]:
    recherche_fichiers(
        os.path.abspath(os.path.join(__file__, '..', '..', 'Données'))
    )
    return dict((item.name, json.load(open(item.path))) for item in os.scandir(os.path.join(cache, 'tailles')))


def main():
    # Extraire la liste des nombres de fragments
    nombres_fragments = list(lire_tailles().values())

    # Afficher les tailles, sous forme d'un histogramme
    fig, plot = plt.subplots(1, 1)
    plot.hist(list(sorted(map(len, nombres_fragments))))  # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)
    # Affiche les tableaux
    afficher_hist()

    return 0


def generer_genes(n=1000, g=1000000) -> typing.List[int]:
    """
    Génère une liste de `n` échantillons d'ADN de taille `g`
    :param n: Le nombre d'échantillons
    :param g: Taille des échantillons
    :return: Liste de tailles de génomes (pour des raisons de performances)
    """
    return [g] * n


# ''.join(random.sample('ATCG' * g, g))

def fracture(echantillons, n=65952, sigma=100):
    """
    Fracture une liste d'échantillons en utilisant un nombre de fractures variables de loi n = Norm(N, sigma)
    :param sigma:
    :param n: Ecart-type
    :param echantillons:
    :return:
    """
    #
    échantillons_fracturés = []

    # randn nous donne la loi normale centrée, on la convertis donc avec cette ligne
    n = n + np.random.randn(len(echantillons)) * sigma

    # On prends chaque échantillon et lui associe un nombre de fractures, converti en int à la volée
    for i, taille_échantillon in zip(map(int, n), echantillons):
        if i > taille_échantillon:
            raise ValueError(
                '\n\n\n'
                'ERREUR: Le nombre de fragments demandé est trop grand, il est impossible de continuer\n'
                'Vous pouvez: \n'
                '\t- Réduire le nombre de fragments demandé, en jouant sur N et sigma\n'
                '\t- Augmenter la taille des échantillons à fracturer\n'
            )
        # Il y aura n - 1 fractures pour obtenir n elements
        nombre_fractures = i - 1

        # Une fracture sera placée avant l'index désigné, donc une fracture
        #   localisée à l'index 1 donnera un morceau de 1 et un de taille - 1
        fractures = set()
        while len(fractures) < nombre_fractures:
            fractures.update(np.random.randint(1, taille_échantillon - 1, nombre_fractures))
        print(nombre_fractures)
        if len(fractures) != nombre_fractures:
            if len(fractures) < nombre_fractures:
                raise RuntimeError('ERREUR: Une erreur s\'est produite: pas assez de fragments')
            fractures = set(list(fractures)[:nombre_fractures])
        fractures = [0] + list(sorted(fractures)) + [taille_échantillon]
        échantillon_fracturé = [
            fractures[i + 1] - fractures[i] for i in range(nombre_fractures + 1)
        ]
        if sum(échantillon_fracturé) != taille_échantillon:
            raise RuntimeError(
                "ERREUR: Une erreur s'est produite lors de la fragmentation, \n"
                f"Nous avons {abs(sum(échantillon_fracturé) - taille_échantillon)}"
                f" fragments de différence entre ce qu'on attendais et ce qu'on a"
            )
        échantillon_fracturé = list(i for i in sorted(échantillon_fracturé) if i < 100)
        échantillons_fracturés.append(échantillon_fracturé)

    return échantillons_fracturés


if __name__ == '__main__':
    recherche_fichiers(
        os.path.abspath(os.path.join(__file__, '..', '..', 'Données')),
        os.path.abspath(os.path.join(__file__, '..', '..')),
        os.path.abspath(os.path.join(__file__, '..'))
    )
    for folder in os.scandir('../Documents/'):
        recherche_fichiers(os.path.abspath(folder.path))


    # ech = fracture(generer_genes(10, 1000000), 65952, 1)[0]
    # fig, plot = plt.subplots(1, 1);
    # for ech in fracture(generer_genes(1, 1000000), 65952, 10000):
    #     plot.hist(ech)  # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)
    # plot.hist(generer_genes(1000, 1000)) # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)
    # plot.hist(generer_genes(1000, 50)) # , bins=np.arange(min(taille), m + 0.2, 0.2), rwidth=0.5)
    # Affiche les tableaux
    exit(plt.show())
    exit(main())

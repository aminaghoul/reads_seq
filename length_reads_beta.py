import matplotlib.pyplot as plt
import numpy as np

import random
import time
from scipy.stats import beta


def get_parts(length, number_of_chunks):
    cuts = set()
    while len(cuts) < number_of_chunks - 1:
        cuts.add(random.randint(1, length - 1))
    cuts = [0] + list(sorted(cuts)) + [length]
    return [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]


def print_parts(parts):
    print('|' + '|'.join(''.join(random.sample('ATCG' * i, i)) for i in parts) + '|')


def show_plot(parts):
    # https://matplotlib.org/gallery/lines_bars_and_markers/simple_plot.html
    # Data for plotting
    t = sorted(set(parts))
    s = [parts.count(i) for i in t]

    fig, ax = plt.subplots()
    ax.plot(t, s)
    ax.set(xlabel='Chunk size', ylabel='Number of occurences',
           title='%s nucleotide long chain divided into %s chunks' % (sum(parts), len(parts)))
    ax.grid()

    fig.savefig("test.png")
    plt.show()


def show_hist(parts):
    # Data for plotting
    t = sorted(set(parts))
    s = [parts.count(i) for i in t]

    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(s, density=10)

    ax.set(xlabel='Chunk size', ylabel='Number of occurences',
           title='%s nucleotide long chain divided into %s chunks' % (sum(parts), len(parts)))
    ax.grid()

    plt.show()



# print_parts(get_parts(100, 10))
show_plot(get_parts(10000000, 5000000))
show_hist(get_parts(1000000, 500000))


# print_parts(get_parts(100, 10))

def main():
    # Data for plotting

    fig, ax = plt.subplots()
    for i in range(6):
        print('Loading %s...' % i)
        parts = get_parts(10 ** i, (10 ** i) / 2)
        T = time.time()
        t = sorted(set(parts))
        s = [parts.count(i) for i in t]
        s = [i / max(s) for i in s]

        ax.plot(t, s)
        ax.set(xlabel='Chunk size', ylabel='Number of occurences',
               title='%s nucleotide long chain divided into %s chunks' % (sum(parts), len(parts)))
        ax.grid()

        plt.savefig("plot-%s.png" % i)
        print('Loading %s...done in %s' % (i, time.time() - T))

    plt.show()





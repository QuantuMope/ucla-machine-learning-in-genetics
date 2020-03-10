import numpy as np
import sys
import itertools
import time


class haplotypes:
    def __init__(self):
        self.haplos = []

    def append(self, x):
        self.haplos.append(x)

    def len(self):
        return len(self.haplos)

    def cutout(self):
        self.haplos = np.unique(self.haplos).tolist()


unique_haplotypes = haplotypes()
sys.setrecursionlimit(20000)


def get_unique_hap_recurse(individual, index, haplotype):
    if index == len(individual):
        unique_haplotypes.append(haplotype)
        return

    snp = individual[index]
    while snp != 1 and index < len(individual):
        if snp == 0:
            haplotype += "0"
        elif snp == 2:
            haplotype += "1"
        index += 1
        if index == len(individual): break
        snp = individual[index]
    if index != len(individual) and snp == 1:
        new_haplo = haplotype
        haplotype += "1"
        get_unique_hap_recurse(individual, index + 1, haplotype)
        new_haplo = new_haplo + "0"
        get_unique_hap_recurse(individual, index + 1, new_haplo)

    if index == len(individual):
        unique_haplotypes.append(haplotype)


def get_unique_loop(data, num_individuals):
    x = [[] for _ in range(num_individuals)]
    for g in range(num_individuals):
        sample_genotype = data[:, g]
        for i in range(4000):
            if i == 3999:
                block_genotype = sample_genotype[i * 10:]
            else:
                block_genotype = sample_genotype[i * 10: (i + 1) * 10]
            num_ones = len(np.where(block_genotype == 1)[0])
            # print("num_ones: {}".format(num_ones))
            for y in range(2 ** num_ones):
                x[g].append(bin(y)[2:])
                # keep 0 and 2 in original same, and replace 1 by whatever is in x
            # print("iterations {}".format(i))


def clarks_algo(num_blocks, num_individuals, block_size):
    for i in range(num_blocks):
        # Clark's algo here
        # Define starting haplos
        unique_haplos = []
        for g in range(num_individuals):
            # Try to match with current one or create a new error
            for b in range(i*block_size, (i+1)*block_size):
                unique_haplos.append()


def main():
    data = np.loadtxt(fname=sys.argv[1])
    num_snps, num_individuals = data.shape
    sample_genotype = data[:, 0]
    num_ones = len(np.where(sample_genotype == 1)[0])

    block_size = 2
    num_blocks = num_snps/block_size


if __name__ == '__main__':
    main()

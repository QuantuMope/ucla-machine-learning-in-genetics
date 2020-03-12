import numpy as np
import sys
from math import ceil
from time import time


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


def clarks_algo(data, num_individuals, num_snps, num_blocks, block_size):
    phased_haplotypes = np.zeros((num_individuals, num_snps*2))
    phased_haplotypes[:] = 9  # no entry should be a 9 at the end
    unique_haplos = []
    remainder = num_snps % block_size
    for b in range(num_blocks):
        unique_haplos.clear()
        h1s = b*block_size
        h1e = (b+1)*block_size
        h2s = b*block_size + num_snps
        h2e = (b+1)*block_size + num_snps
        if remainder != 0 and b == num_blocks-1:
            h1e = h1s + remainder
            h2e = h2s + remainder
            block_size = remainder
        for i in range(num_individuals):
            geno_snp = data[i, h1s:h1e]
            assert geno_snp.shape == (block_size,)  # will need to reshape if this is false
            found_match = False
            if i == 0:
                haplo_1, haplo_2 = phase_from_geno(geno_snp, unique_haplos, block_size)
                phased_haplotypes[i, h1s:h1e] = haplo_1
                phased_haplotypes[i, h2s:h2e] = haplo_2
                continue
            # Try to find a match with existing haplotypes
            for x, haplo_1 in enumerate(unique_haplos):
                for y, haplo_2 in enumerate(unique_haplos, start=x):
                    if np.array_equal(geno_snp, haplo_1 + haplo_2):
                        phased_haplotypes[i, h1s:h1e] = haplo_1
                        phased_haplotypes[i, h2s:h2e] = haplo_2
                        found_match = True
                        break
                if found_match:
                    break
            # If no matches with current list.
            for haplo in unique_haplos:
                potential_haplo = geno_snp - haplo
                if not np.any(potential_haplo == 2) and not np.any(potential_haplo == -1):
                    print('half match')
                    unique_haplos.append(potential_haplo)
                    phased_haplotypes[i, h1s:h1e] = haplo
                    phased_haplotypes[i, h2s:h2e] = potential_haplo
                    found_match = True
                    break
            if not found_match:
                haplo_1, haplo_2 = phase_from_geno(geno_snp, unique_haplos, block_size)
                phased_haplotypes[i, h1s:h1e] = haplo_1
                phased_haplotypes[i, h2s:h2e] = haplo_2

        print(len(unique_haplos))

    assert not np.any(phased_haplotypes == 9)
    assert not np.any(phased_haplotypes == 2)

    return phased_haplotypes


def phase_from_geno(geno_snp, unique_haplos, block_size):
    haplo_1 = np.zeros((block_size,))
    haplo_2 = np.zeros((block_size,))
    for i in range(len(geno_snp)):
        allele = geno_snp[i]
        if allele == 0:
            haplo_1[i] = 0
            haplo_2[i] = 0
        elif allele == 2:
            haplo_1[i] = 1
            haplo_2[i] = 1
        elif allele == 1:
            haplo_1[i] = 0
            haplo_2[i] = 1
    unique_haplos.append(haplo_1)
    unique_haplos.append(haplo_2)
    return haplo_1, haplo_2



def main():
    data = np.loadtxt(fname=sys.argv[1])
    data = data.T
    num_individuals, num_snps = data.shape

    size = 238
    sample_genotype = data[:, :size]
    block_size = 5
    num_blocks = ceil(size / block_size)
    print(num_blocks)
    num_snps = size
    phased_haplotypes = clarks_algo(sample_genotype, num_individuals, num_snps, num_blocks, block_size)
    print(phased_haplotypes)

    np.savetxt('phased', data.T, fmt='%d', delimiter=' ')


if __name__ == '__main__':
    main()

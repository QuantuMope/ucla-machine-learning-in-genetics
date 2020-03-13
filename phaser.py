import numpy as np
import sys
from math import ceil
from time import time


"""
    Phases imputed genotypes into haplotypes.
    
    Uses a tiling variant of Clark's algorithm.
    A minimum parsimony approach is used to get
    the haplotype phases for blocks of 12 SNPs
    at a time.
"""


def update_phased_haplos(phased_haplotypes, i, h1s, h1e, haplo_1, haplo_2):
    phased_haplotypes[2 * i, h1s:h1e] = haplo_1
    phased_haplotypes[2 * i + 1, h1s:h1e] = haplo_2


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


def clarks_algo(data, num_individuals, num_snps, num_blocks, block_size):
    phased_haplotypes = np.zeros((num_individuals*2, num_snps))
    phased_haplotypes[:] = 9  # no entry should be a 9 at the end
    unique_haplos = []
    remainder = num_snps % block_size
    for b in range(num_blocks):
        unique_haplos.clear()
        h1s = b*block_size
        h1e = (b+1)*block_size
        if remainder != 0 and b == num_blocks-1:
            h1e = h1s + remainder
            block_size = remainder
        for i in range(num_individuals):
            geno_snp = data[i, h1s:h1e]
            found_match = False
            if i == 0:
                haplo_1, haplo_2 = phase_from_geno(geno_snp, unique_haplos, block_size)
                update_phased_haplos(phased_haplotypes, i, h1s, h1e, haplo_1, haplo_2)
                continue
            # Try to find a match with existing haplotypes
            for x, haplo_1 in enumerate(unique_haplos):
                for y, haplo_2 in enumerate(unique_haplos, start=x):
                    if np.array_equal(geno_snp, haplo_1 + haplo_2):
                        update_phased_haplos(phased_haplotypes, i, h1s, h1e, haplo_1, haplo_2)
                        found_match = True
                        break
                if found_match: break
            if found_match: continue
            # If no matches with current list, find new complement.
            for haplo in unique_haplos:
                potential_haplo = geno_snp - haplo
                if not np.any(potential_haplo == 2) and not np.any(potential_haplo == -1):
                    unique_haplos.append(potential_haplo)
                    update_phased_haplos(phased_haplotypes, i, h1s, h1e, haplo, potential_haplo)
                    found_match = True
                    break
            # If no possible initial pair, add two new haplos.
            if not found_match:
                haplo_1, haplo_2 = phase_from_geno(geno_snp, unique_haplos, block_size)
                update_phased_haplos(phased_haplotypes, i, h1s, h1e, haplo_1, haplo_2)

    assert not np.any(phased_haplotypes == 9)
    assert not np.any(phased_haplotypes == 2)
    assert not np.any(phased_haplotypes == -1)

    return phased_haplotypes


def main():
    data = np.loadtxt(fname=sys.argv[1])
    data = data.T
    num_individuals, num_snps = data.shape

    block_size = 12

    num_blocks = ceil(num_snps / block_size)

    start = time()
    phased_haplotypes = clarks_algo(data, num_individuals, num_snps, num_blocks, block_size)
    print('Phasing took {} seconds'.format(round(time() - start, 2)))

    np.savetxt('phased', phased_haplotypes.T, fmt='%d', delimiter=' ')


if __name__ == '__main__':
    main()

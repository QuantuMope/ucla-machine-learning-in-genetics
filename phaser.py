import numpy as np
import sys
import csv


def main():
    data = np.loadtxt(fname=sys.argv[1])
    num_individuals, num_snps = data.shape

    allele_counts = np.zeros((num_snps, 4))

    for i in range(num_snps):
        _, counts = np.unique(data[:, i], return_counts=True)
        allele_counts[i, :] = counts

    print(allele_counts)


if __name__ == '__main__':
    main()
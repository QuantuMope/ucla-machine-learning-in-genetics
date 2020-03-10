import numpy as np
import sys
import csv
import random


def main():
    data = np.loadtxt(fname=sys.argv[1])
    num_individuals, num_snps = data.shape

    allele_counts = np.zeros((num_snps, 4))

    # Transition probabilities of next SNP
    # Probability that the next SNP is either a 0, 1, 2
    # Nine values in total. Every 3 must sum up to 1.
    transition_probabilities = np.zeros((num_snps, 9))

    for i in range(num_snps):
        _, counts = np.unique(data[:, i], return_counts=True)
        allele_counts[i, :] = counts

    for i in range(num_snps):
        if i == 0:
            masked = np.where(data[:, 0] == 9)[0]
            for mask in masked:
                data[mask, 0] = random.choice([0, 2])
            prev_zeroes = np.where(data[:, 0] == 0)[0]
            prev_ones = np.where(data[:, 0] == 1)[0]
            prev_twos = np.where(data[:, 0] == 2)[0]
            continue

        next_from_zeros = data[prev_zeroes, i]
        next_from_ones = data[prev_ones, i]
        next_from_twos = data[prev_twos, i]

        zeroes_after_zeroes = len(np.where(next_from_zeros == 0)[0])
        ones_after_zeroes = len(np.where(next_from_zeros == 1)[0])
        twos_after_zeroes = len(np.where(next_from_zeros == 2)[0])
        total_zeroes = zeroes_after_zeroes + ones_after_zeroes + twos_after_zeroes

        zeroes_after_ones = len(np.where(next_from_ones == 0)[0])
        ones_after_ones = len(np.where(next_from_ones == 1)[0])
        twos_after_ones = len(np.where(next_from_ones == 2)[0])
        total_ones = zeroes_after_ones + ones_after_ones + twos_after_ones

        zeroes_after_twos = len(np.where(next_from_twos == 0)[0])
        ones_after_twos = len(np.where(next_from_twos == 1)[0])
        twos_after_twos = len(np.where(next_from_twos == 2)[0])
        total_twos = zeroes_after_twos + ones_after_twos + twos_after_twos

        transition_probabilities[i, :3] = [zeroes_after_zeroes/total_zeroes,
                                           ones_after_zeroes/total_zeroes,
                                           twos_after_zeroes/total_zeroes]
        transition_probabilities[i, 3:6] = [zeroes_after_ones/total_ones,
                                            ones_after_ones/total_ones,
                                            twos_after_ones/total_ones]
        transition_probabilities[i, 6:] = [zeroes_after_twos/total_twos,
                                           ones_after_twos/total_twos,
                                           twos_after_twos/total_twos]

        print(transition_probabilities[i, :])

        trp = transition_probabilities

        z_assign = 0 if trp[i, 0] > trp[i, 2] else 2
        o_assign = 0 if trp[i, 3] > trp[i, 5] else 2
        t_assign = 0 if trp[i, 6] > trp[i, 8] else 2

        masked = np.where(data[:, i] == 9)[0]

        for mask in masked:
            if data[mask, i-1] == 0:
                data[mask, i] = z_assign
            elif data[mask, i-1] == 1:
                data[mask, i] = o_assign
            elif data[mask, i-1] == 2:
                data[mask, i] = t_assign

    np.savetxt('output.txt', data, fmt='%d', delimiter=' ')





if __name__ == '__main__':
    main()
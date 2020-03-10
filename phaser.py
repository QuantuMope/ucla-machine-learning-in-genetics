import numpy as np
import sys
import csv
import random


def main():
    data = np.loadtxt(fname=sys.argv[1])
    num_individuals, num_snps = data.shape

    # Transition probabilities of next SNP
    # Probability that the next SNP is either a 0 or 2 when prev SNP was 0, 1, 2
    trans_probs_1 = np.zeros((num_snps, 6))

    # trans probs of next SNP when prev two SNPs were either 00, 01, 02, 10, 11, 12, 20, 21, 22
    trans_probs_2 = np.zeros((num_snps, 18))

    # for i in range(num_snps):
    #     _, counts = np.unique(data[:, i], return_counts=True)
    #     allele_counts[i, :] = counts

    for i in range(num_snps):
        if i == 0:
            masked = np.where(data[:, 0] == 9)[0]
            for mask in masked:
                data[mask, 0] = random.choice([0, 2])
            prev_zeroes = np.where(data[:, 0] == 0)[0]
            prev_ones = np.where(data[:, 0] == 1)[0]
            prev_twos = np.where(data[:, 0] == 2)[0]

        elif i == 1:
            next_from_zeros = data[prev_zeroes, i]
            next_from_ones = data[prev_ones, i]
            next_from_twos = data[prev_twos, i]

            zeroes_after_zeroes = len(np.where(next_from_zeros == 0)[0])
            twos_after_zeroes = len(np.where(next_from_zeros == 2)[0])
            total_zeroes = zeroes_after_zeroes + twos_after_zeroes

            zeroes_after_ones = len(np.where(next_from_ones == 0)[0])
            twos_after_ones = len(np.where(next_from_ones == 2)[0])
            total_ones = zeroes_after_ones + twos_after_ones

            zeroes_after_twos = len(np.where(next_from_twos == 0)[0])
            twos_after_twos = len(np.where(next_from_twos == 2)[0])
            total_twos = zeroes_after_twos + twos_after_twos

            trans_probs_1[i, :2] = [zeroes_after_zeroes/total_zeroes,
                                    twos_after_zeroes/total_zeroes]
            trans_probs_1[i, 2:4] = [zeroes_after_ones/total_ones,
                                     twos_after_ones/total_ones]
            trans_probs_1[i, 4:] = [zeroes_after_twos/total_twos,
                                    twos_after_twos/total_twos]

            trp = trans_probs_1

            z_assign = 0 if trp[i, 0] > trp[i, 1] else 2
            o_assign = 0 if trp[i, 2] > trp[i, 3] else 2
            t_assign = 0 if trp[i, 4] > trp[i, 5] else 2

            masked = np.where(data[:, i] == 9)[0]

            for mask in masked:
                if data[mask, i-1] == 0:
                    data[mask, i] = z_assign
                elif data[mask, i-1] == 1:
                    data[mask, i] = o_assign
                elif data[mask, i-1] == 2:
                    data[mask, i] = t_assign

        elif i == 2:
            prev_indexes = []
            # zz, zo, zt, oz, oo, ot, tz, to, tt
            for p in range(3):
                for j in range(3):
                    prev_indexes.append(np.where(np.logical_and(data[:, i-2] == p, data[:, i-1] == j))[0])

            snps = []
            for prev in prev_indexes:
                snps.append(data[prev, i])

            counts = []
            for snp in snps:
                counts.append(len(np.where(snp == 0)[0]))
                counts.append(len(np.where(snp == 2)[0]))

            for j in range(0, len(counts), 2):
                total = counts[j] + counts[j+1]
                trans_probs_2[i, j] = counts[j] / total
                trans_probs_2[i, j+1] = counts[j+1] / total

            trp = trans_probs_2

            assigns = []
            for j in range(0, len(trp[i]), 2):
                assigns.append(0 if trp[i, j] > trp[i, j+1] else 2)

            masked = np.where(data[:, i] == 9)[0]

            for mask in masked:
                counter = 0
                pp, p = data[mask, i-2:i]
                for z in range(3):
                    for j in range(3):
                        if pp == z and p == j:
                            data[mask, i] = assigns[counter]
                        counter += 1
        else:
            # 27 total possible permutations
            prev_indexes = []
            for p in range(3):
                for j in range(3):
                    for t in range(3):
                        prev_indexes.append(np.where(np.logical_and.reduce((data[:, i-3] == p,
                                                                            data[:, i-2] == j,
                                                                            data[:, i-1] == t)))[0])

            snps = []
            for prev in prev_indexes:
                snps.append(data[prev, i])

            counts = []
            for snp in snps:
                counts.append(len(np.where(snp == 0)[0]))
                counts.append(len(np.where(snp == 2)[0]))

            trp = []
            for j in range(0, len(counts), 2):
                total = counts[j] + counts[j + 1]
                trp.append(counts[j] / total)
                trp.append(counts[j+1] / total)

            assigns = []
            for j in range(0, len(trp), 2):
                assigns.append(0 if trp[j] > trp[j+1] else 2)
                assert abs((trp[j] + trp[j+1]) - 1.) < 0.01

            masked = np.where(data[:, i] == 9)[0]

            for mask in masked:
                counter = 0
                ppp, pp, p = data[mask, i-3:i]
                for z in range(3):
                    for j in range(3):
                        for t in range(3):
                            if ppp == z and pp == j and p == t:
                                data[mask, i] = assigns[counter]
                            counter += 1

        print("Processing SNP {}".format(i), flush=True)

    np.savetxt('output', data, fmt='%d', delimiter=' ')


if __name__ == '__main__':
    main()
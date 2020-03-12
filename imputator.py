import numpy as np
import sys
import random
import sys
from time import time


def main():
    data = np.loadtxt(fname=sys.argv[1])
    data = data.T
    num_individuals, num_snps = data.shape
    start = time()

    for i in range(num_snps):
        if i == 0:
            masked = np.where(data[:, 0] == 9)[0]
            for mask in masked:
                data[mask, 0] = random.choice([0, 2])

        elif i == 1:
            prev_indexes = []
            for j in range(3):
                prev_indexes.append(np.where(data[:, i-1] == j)[0])

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
                if total == 0:
                    trp.append(0)
                    trp.append(0)
                    continue
                trp.append(counts[j] / total)
                trp.append(counts[j + 1] / total)

            assigns = []
            for j in range(0, len(trp), 2):
                assigns.append(0 if trp[j] > trp[j+1] else 2)

            masked = np.where(data[:, i] == 9)[0]

            for mask in masked:
                counter = 0
                p = data[mask, i-1]
                for j in range(3):
                    if p == j:
                        data[mask, i] = assigns[counter]
                    counter += 1

        # elif i == 2:
        else:
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

            trp = []
            for j in range(0, len(counts), 2):
                total = counts[j] + counts[j+1]
                if total == 0:
                    trp.append(0)
                    trp.append(0)
                    continue
                trp.append(counts[j] / total)
                trp.append(counts[j+1] / total)

            assigns = []
            for j in range(0, len(trp), 2):
                assigns.append(0 if trp[j] > trp[j+1] else 2)

            masked = np.where(data[:, i] == 9)[0]

            for mask in masked:
                counter = 0
                pp, p = data[mask, i-2:i]
                for z in range(3):
                    for j in range(3):
                        if pp == z and p == j:
                            data[mask, i] = assigns[counter]
                        counter += 1
        # else:
        #     # 27 total possible permutations
        #     prev_indexes = []
        #     for p in range(3):
        #         for j in range(3):
        #             for t in range(3):
        #                 prev_indexes.append(np.where(np.logical_and.reduce((data[:, i-3] == p,
        #                                                                     data[:, i-2] == j,
        #                                                                     data[:, i-1] == t)))[0])
        #
        #     snps = []
        #     for prev in prev_indexes:
        #         snps.append(data[prev, i])
        #
        #     counts = []
        #     for snp in snps:
        #         counts.append(len(np.where(snp == 0)[0]))
        #         counts.append(len(np.where(snp == 2)[0]))
        #
        #     trp = []
        #     for j in range(0, len(counts), 2):
        #         total = counts[j] + counts[j + 1]
        #         if total == 0:
        #             trp.append(0)
        #             trp.append(0)
        #             continue
        #         trp.append(counts[j] / total)
        #         trp.append(counts[j+1] / total)
        #
        #     assigns = []
        #     for j in range(0, len(trp), 2):
        #         assigns.append(0 if trp[j] > trp[j+1] else 2)
        #
        #     masked = np.where(data[:, i] == 9)[0]
        #
        #     for mask in masked:
        #         counter = 0
        #         ppp, pp, p = data[mask, i-3:i]
        #         for z in range(3):
        #             for j in range(3):
        #                 for t in range(3):
        #                     if ppp == z and pp == j and p == t:
        #                         data[mask, i] = assigns[counter]
        #                     counter += 1
        #
        # elif i == 4:
        #     # 64 total possible permutations
        #     prev_indexes = []
        #     for p in range(3):
        #         for j in range(3):
        #             for t in range(3):
        #                 for z in range(3):
        #                     prev_indexes.append(np.where(np.logical_and.reduce((data[:, i-4] == p,
        #                                                                         data[:, i-3] == j,
        #                                                                         data[:, i-2] == t,
        #                                                                         data[:, i-1] == z)))[0])
        #
        #     snps = []
        #     for prev in prev_indexes:
        #         snps.append(data[prev, i])
        #
        #     counts = []
        #     for snp in snps:
        #         counts.append(len(np.where(snp == 0)[0]))
        #         counts.append(len(np.where(snp == 2)[0]))
        #
        #     trp = []
        #     for j in range(0, len(counts), 2):
        #         total = counts[j] + counts[j + 1]
        #         if total == 0:
        #             trp.append(0)
        #             trp.append(0)
        #             continue
        #         trp.append(counts[j] / total)
        #         trp.append(counts[j+1] / total)
        #
        #     assigns = []
        #     for j in range(0, len(trp), 2):
        #         assigns.append(0 if trp[j] > trp[j+1] else 2)
        #
        #     masked = np.where(data[:, i] == 9)[0]
        #
        #     for mask in masked:
        #         counter = 0
        #         pppp, ppp, pp, p = data[mask, i-4:i]
        #         for z in range(3):
        #             for j in range(3):
        #                 for t in range(3):
        #                     for n in range(3):
        #                         if pppp == z and ppp == j and pp == t and p == n:
        #                             data[mask, i] = assigns[counter]
        #                         counter += 1
        # else:
        #     # 125 total possible permutations
        #     prev_indexes = []
        #     for p in range(3):
        #         for j in range(3):
        #             for t in range(3):
        #                 for z in range(3):
        #                     for y in range(3):
        #                         prev_indexes.append(np.where(np.logical_and.reduce((data[:, i-5] == p,
        #                                                                             data[:, i-4] == j,
        #                                                                             data[:, i-3] == t,
        #                                                                             data[:, i-2] == z,
        #                                                                             data[:, i-1] == y)))[0])
        #
        #     snps = []
        #     for prev in prev_indexes:
        #         snps.append(data[prev, i])
        #
        #     counts = []
        #     for snp in snps:
        #         counts.append(len(np.where(snp == 0)[0]))
        #         counts.append(len(np.where(snp == 2)[0]))
        #
        #     trp = []
        #     for j in range(0, len(counts), 2):
        #         total = counts[j] + counts[j + 1]
        #         if total == 0:
        #             trp.append(0)
        #             trp.append(0)
        #             continue
        #         trp.append(counts[j] / total)
        #         trp.append(counts[j+1] / total)
        #
        #     assigns = []
        #     for j in range(0, len(trp), 2):
        #         assigns.append(0 if trp[j] > trp[j+1] else 2)
        #
        #     masked = np.where(data[:, i] == 9)[0]
        #
        #     for mask in masked:
        #         counter = 0
        #         ppppp, pppp, ppp, pp, p = data[mask, i-5:i]
        #         for z in range(3):
        #             for j in range(3):
        #                 for t in range(3):
        #                     for n in range(3):
        #                         for y in range(3):
        #                             if ppppp == z and pppp == j and ppp == t and pp == n and p == y:
        #                                 data[mask, i] = assigns[counter]
        #                             counter += 1

        # print("Processing SNP {}".format(i), flush=True)

    np.savetxt('output', data.T, fmt='%d', delimiter=' ')
    print("Imputation took {} seconds".format(round(time() - start, 2)))


if __name__ == '__main__':
    main()
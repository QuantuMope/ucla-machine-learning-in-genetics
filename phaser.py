import numpy as np
import sys
import csv


def main():
    if len(sys.argv) != 2:
        print('Program expects exactly one argument but {} were given.'.format(len(sys.argv)-1))
        sys.exit(1)

    file_name = sys.argv[1]

    try:
        file = open(file_name)
    except FileNotFoundError:
        print('File {} could not be found'.format(file_name))
    except IOError:
        print('File {} could not be opened/read'.format(file_name))

    file_reader = csv.reader(file, delimiter=' ')
    for entry in file:
        print(entry)
    file.close()


if __name__ == '__main__':
    main()
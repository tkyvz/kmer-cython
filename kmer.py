import argparse
import bfcounter
import dsk
import heapq
import math
import time
import sys


def check_positive(value):
    """
    Checks value is a positive integer, returns True if so, else raise error.
    :param  value: value to be checked
    """
    try:
        ivalue = int(value)
        if ivalue <= 0:
            # is int but non-positive
            raise argparse.ArgumentTypeError(
                '{} is an invalid positive integer value'.format(value))
        return ivalue
    except ValueError:
        # not int
        raise argparse.ArgumentTypeError('{} is not an integer'.format(value))


def check_between_zero_one(value):
    """
    Checks value is between 0 and 1, returns True if so, else raise error.
    :param  value: value to be checked
    """
    try:
        fvalue = float(value)
        if 0 <= fvalue < 1:
            return fvalue
        # is float but not between 0 and 1
        raise argparse.ArgumentTypeError(
            'should be between 0 and 1'.format(value))
    except ValueError:
        # not float
        raise argparse.ArgumentTypeError(
            '{} is not a floating point number'.format(value))


def count_kmers(file_name, k, verbose=False):
    """
    Counts how many k-mers exists in a given file
    :param  file_name: Fastq file to be counted
    :param  k: K-mer size
    """
    if verbose:
        start = time.time()
        print('Counting kmers in {}'.format(file_name))
    total_kmers = 0
    with open(file_name, 'r') as f:
        line_num = 0
        for line in f:
            if line_num % 4 == 1:  # dna sequence
                total_kmers += len(line) - k  # eliminate new-line
            line_num += 1
    if verbose:
        end = time.time()
        print('{} kmers are counted in {:.2f} seconds'.format(
            total_kmers, end - start))
    return total_kmers


def parameters(total_kmers, target_disk, target_memory, k, verbose=False):
    """
    Calculates paramters
    :param  total_kmers: Number of total kmers in the file
    :param  target_disk: Target disk space
    :param  target_memory: Target memory
    :param  k: K-mer size
    :return number_of_iterations: Number of iterations for DSK
            number_of_partitions: Number of partitions for DSK
            bloom_filter_capacity: Bloom Filter Capacity
            dsk: True if DSK should be used, False if BFCounter should be used
    """
    if verbose:
        print('Calculating paramters...')
        print('Total K-mers: {}'.format(total_kmers))
        print('Target Disk: {:.2f}GB'.format(target_disk / 1024**3 / 8))
        print('Target Memory: {:.2f}GB'.format(target_memory / 1024**3 / 8))

    b = (k + sys.getsizeof('')) * 8  # size of kmer in bits
    b_disk = (k + 1) * 8  # size of kmer in disk in bits
    # max unique kmers
    if (4 ** k) < total_kmers:
        v = 4 ** k
    else:
        v = total_kmers
    # DSK Paramters
    number_of_iterations = int(
        math.ceil(float(total_kmers) * b_disk / target_disk))
    number_of_partitions = int(
        math.ceil((v * b) / (0.7 * target_memory * number_of_iterations)))
    use_dsk = 0.7 * target_memory < b * v
    # Bloom Filter Capacity
    if dsk:
        capacity = total_kmers / (number_of_iterations *
                                  number_of_partitions)
    else:
        capacity = total_kmers

    if verbose:
        print('Parameters are calculated')
        if use_dsk:
            print('Algorithm: DSK')
            print('# of iterations: {}'.format(number_of_iterations))
            print('# of partitions: {}'.format(number_of_partitions))
        else:
            print('Algorithm: BFCounter')
        print('Bloom Filter Capacity: {}'.format(capacity))

    return number_of_iterations, number_of_partitions, capacity, use_dsk


def run():
    desc = 'A program for counting most frequent k-mers in .fastq files.'
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.RawTextHelpFormatter)
    # file name
    parser.add_argument('-f',
                        '--file-name',
                        required=True,
                        help='Name of the .fastq file to be processed',
                        type=str)
    # kmer size
    parser.add_argument('-k',
                        '--kmer-size',
                        required=True,
                        help='Length of k-mers',
                        type=check_positive)
    # most frequent count
    parser.add_argument('-n',
                        '--most-frequent',
                        required=True,
                        help='Number of most frequent k-mers to be outputted',
                        type=check_positive)
    # Bloom Filter error rate for BFCounter
    parser.add_argument('-e',
                        '--error-rate',
                        help='Bloom Filter error rate for BFCounter algorithm',
                        default=1e-2,
                        type=check_between_zero_one)
    # DSK target disk space
    parser.add_argument('-d',
                        '--target-disk',
                        help='Target Disk Space in GB for DSK algorithm',
                        default=50,
                        type=check_positive)
    # DSK target disk memory
    parser.add_argument('-m',
                        '--target-memory',
                        help='Target Memory in GB for DSK algorithm',
                        default=4,
                        type=check_positive)
    # Verbose
    parser.add_argument('-v',
                        '--verbose',
                        help='Verbose',
                        action='store_true')

    # Parameters
    cli_args = parser.parse_args()
    verbose = cli_args.verbose  # verbose
    error_rate = cli_args.error_rate  # bloom filter error rate
    target_disk = cli_args.target_disk * (1024**3) * 8  # target disk
    target_memory = cli_args.target_memory * (1024**3) * 8  # target memory
    n = cli_args.most_frequent  # n, the number of kmers to be returned
    k = cli_args.kmer_size  # k, size of kmer
    file_name = cli_args.file_name  # fastq file name

    start = time.time()

    # Count total k-mers
    total_kmers = count_kmers(file_name, k, verbose=verbose)
    # Calculate paramters
    iters, parts, capacity, is_dsk = parameters(total_kmers,
                                                target_disk,
                                                target_memory,
                                                k,
                                                verbose=verbose)
    if is_dsk:  # DSK Algorithm implementation
        heap = dsk.dsk(
            file_name, k, n, capacity, error_rate, iters, parts, verbose
        )
    else:  # BFCounter Algorithm implementation
        heap = bfcounter.bfcounter(
            file_name, k, n, capacity, error_rate, verbose
        )

    for count, kmer in heapq.nlargest(n, heap):
        print('{}: {}'.format(kmer, count))

    end = time.time()
    print('Duration: {:.2f} seconds'.format(end - start))


if __name__ == '__main__':
    run()

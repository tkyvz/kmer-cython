from collections import defaultdict
import heapq
import math
import os
import time
import sys

try:
    from pybloomfilter import BloomFilter
except ImportError:
    raise ImportError('pybloomfilter module is required. ' +
                      'Use pip install pybloomfiltermmap3')


def bfcounter(file_name, int k, int n, long capacity, double error_rate, verbose=False):
    """
    Implementation of Bloom Filter k-mer Counting algorithm

    Creates a Bloom Filter and check the k-mer is previously encountered or
    not. Only previously encountered k-mers are added to the Hash Table, which
    drastically reduced the size of the Hash Table.

    Referenced Paper:
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-333
    """
    if verbose:
        start = time.time()
        print('BFCounter started.')

    heap = []
    for i in range(n):
        heap.append((0, ''))

    bf = BloomFilter(capacity, error_rate, 'kmer_bf')

    kmer_counter = defaultdict(lambda: 1)

    # Assign functions to local variables for performance improvement
    add_to_bf = bf.add
    heap_pushpop = heapq.heappushpop

    with open(file_name, 'r') as f:
        line_num = 0
        for line in f:
            if line_num % 4 == 1:  # dna sequence
                kmer_count = len(line) - k
                for i in range(kmer_count):
                    kmer = line[i:i + k]
                    if kmer not in bf:  # not in Bloom Filter
                        add_to_bf(kmer)
                    else:  # in Bloom Filter
                        kmer_counter[kmer] += 1
            line_num += 1
    if verbose:
        end_hash = time.time()
        hash_table_size = sys.getsizeof(kmer_counter) / (1024 ** 2)
        print('Hash table is created in {:.2f} seconds.'.format(
            end_hash - start))
        print('Hash table size: {:.2f} MB.'.format(hash_table_size))
        start_populate = time.time()
        print('Populating the heap...')

    for kmer, count in kmer_counter.items():
        # insert to the heap if count is bigger than minimum
        if count > heap[0][0]:
            heap_pushpop(heap, (count, kmer))

    if verbose:
        end_populate = time.time()
        print('Heap is populated in {:.2f} seconds.'.format(
            end_populate - start_populate
        ))

    os.remove('kmer_bf')
    if verbose:
        end = time.time()
        print('BFCounter is completed in {:.2f} seconds.'.format(end - start))

    return heap

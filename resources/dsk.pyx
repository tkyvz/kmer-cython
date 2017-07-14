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

try:
    import mmh3
except ImportError:
    raise ImportError('mmh3 module is required. Use pip install mmh3')


def dsk(file_name, int k, int n, long capacity, double error_rate,
        int iters, int parts, verbose=False):
  """
  Implementation of DSK, k-mer counting with very low memory algorithm

  Hashes each kmer and puts them into different files according to their
  hash values. By using target disk and memory spaces, determines how many
  distinct files should be used and in how many iterations the program needs
  to perform.

  Referenced Paper:
  http://minia.genouest.org/dsk/
  """
  if verbose:
      start = time.time()

  # Assign functions to local variables for performance improvement
  hash_function = mmh3.hash
  heap_pushpop = heapq.heappushpop

  cdef int CHUNK_LIMIT = math.floor(capacity / 10)  # write approximately in 10 calls

  heap = []
  for i in range(n):
      heap.append((0, ''))
  for it in range(iters):  # iteration
      if verbose:
          start_iter = time.time()
          print('Iteration#{} started.'.format(it + 1))
      files = [open('{}'.format(j), 'w') for j in range(parts)]  # open files

      # Write to files in chunks to have less file.write calls
      chunks = [[] for j in range(parts)]

      # Assign functions to local variables for performance improvement
      writers = [files[j].write for j in range(parts)]
      chunk_appender = [chunks[j].append for j in range(parts)]
      chunk_cleaner = [chunks[j].clear for j in range(parts)]
      chunk_joiner = ''.join

      with open(file_name, 'r') as f:
          line_num = 0
          for line in f:
              if line_num % 4 == 1:  # dna sequence
                  kmer_count = len(line) - k
                  for i in range(kmer_count):
                      kmer = line[i:i + k]
                      h = hash_function(kmer)
                      if h % iters == it:  # belongs to this iteration
                          j = (h / iters) % parts
                          _j = int(j)
                          chunk_appender[_j](kmer + '\n')
                          if len(chunks[_j]) == CHUNK_LIMIT:
                              # write to file
                              writers[_j](chunk_joiner(chunks[_j]))
                              chunk_cleaner[_j]()
              line_num += 1

      # Write remaining kmers
      for j in range(parts):
          writers[j](chunk_joiner(chunks[j]))

      for f in files:
          f.close()  # close files

      del chunks

      if verbose:
          end_disk_write = time.time()
          print('Disk write is completed in {:.2f} seconds.'.format(
              end_disk_write - start_iter
          ))

      for j in range(parts):
          bf = BloomFilter(capacity, error_rate, 'kmer_bf')

          kmer_counter = defaultdict(lambda: 1)

          # Assign functions to local variables for performance improvement
          add_to_bf = bf.add

          if verbose:
              start_partition = time.time()
              print('Reading partition#{} started.'.format(j + 1))

          with open(str(j), 'r') as f:
              for kmer in f:
                  if kmer not in bf:  # not in Bloom Filter
                      add_to_bf(kmer)
                  else:  # in Bloom Filter
                      kmer_counter[kmer] += 1

          if verbose:
              end_partition = time.time()
              print('Reading partition#{} is completed '.format(j + 1) +
                    'in {:.2f} seconds.'.format(
                    end_partition - start_partition))
              print('Hash table size: {:.2f}MB.'.format(sys.getsizeof(kmer_counter) / 1024**2))
              start_populate = time.time()
              print('Populating the heap...')

          for kmer, count in kmer_counter.items():
              # insert to the heap if count is bigger than minimum
              if count > heap[0][0]:
                  heap_pushpop(heap, (count, kmer.rstrip()))

          if verbose:
              end_populate = time.time()
              print('Heap is populated in {:.2f} seconds.'.format(
                  end_populate - start_populate
              ))

          os.remove(str(j))
          os.remove('kmer_bf')

  if verbose:
      end = time.time()
      print('DSK Duration: {:.2f} seconds.'.format(end - start))
  return heap

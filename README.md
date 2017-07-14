# K-mer Counter
K-mer counter counts the most frequent `n` `k-mers` in a given FASTQ file. This project has been developed and tested with Python 3.5.2.

## Getting Started

### Algorithms
K-mer counter uses two different algorithms:
1. [BFCounter](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-333)
   BFCounter uses Bloom Filter for eliminating the unique k-mers in the FASTQ file. This operation reduces the number of items to be put in the hash table drastically. Thus, minimizes the memory usage.
2. [DSK](http://minia.genouest.org/dsk/)
   DSK uses hash values of the k-mers for grouping and partitioning the k-mers according to the specified target disk and memory spaces. K-mers are written in separate files, according to their hash value, for making the memory usage within the target value. For making the disk usage within the target value, this operation is done in multiple iterations.

### Prerequisites
K-mer counter uses Python 3.5.2, and usage of virtual environment is encouraged.
If you don't have virtual environment, you can use
```
pip install virtualenv
```
for installing the virtual environment module. Then create and activate the virtual environment.
```
virtualenv -p ${PYTHON_PATH} venv
source venv\bin\activate
```
For installing required modules, use following commands while virtual environment is activated.
```
pip install -r requirements.txt
```

#### Used Packages
1. **pybloomfiltermmap3:** Bloom filter implementation for Python 3
2. **mmh3:** MurmurHash implementation for Python
3. **Cython:** Required by *pybloomfiltermmap3* module

There are also other packages used for code styling an linter purposes, such as
*flake8*, *pep8* etc.

## Running
Before running the project, Cython dependencies should be built. This can be achieved by `make` command.

Project can be run as follows
```
python kmer.py --file-name ${FASTQ_FILE} --kmer-size ${KMER_SIZE} --most-frequent ${MOST_FREQUENT} --error-rate ${ERROR_RATE} --target-disk ${TARGET_DISK} --target-memory ${TARGET_MEMORY} --algorithm ${ALGORITHM} --verbose
```

After running the program, if you want to delete built files, use `make clean` command.

### Parameters
* **FASTQ_FILE:** Name of the .fastq file in which k-mers will be counted. (required)
* **KMER_SIZE:** K (required)
* **MOST_FREQUENT:** N (required)
* **ERROR_RATE:** Bloom Filter error rate, used only for *BFCounter* algorithm. (DEFAULT=0.01)
* **TARGET_DISK:** Target disk space that will be used in Gigabytes, used only for *DSK* algorithm. (DEFAULT=50)
* **TARGET_MEMORY:** Target memory space that will be used in Gigabytes, used only for *DSK* algorithm. (DEFAULT=4)
* **VERBOSE:** For printing elapsed time and the hash table memory usage.

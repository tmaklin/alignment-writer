# alignment-writer
Pack/unpack [Themisto](https://github.com/algbio/themisto)
or [Fulgor](https://github.com/jermp/fulgor)
pseudoalignment files into/from a compact representation using the
[BitMagic](https://github.com/tlk00/BitMagic) library.

# Installation
## Compiling from source
### Requirements
- C++17 compliant compiler.
- cmake
- git
- OpenMP (optional)

### How-to
- Clone the repository, enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
- This will compile the alignment-writer executable in build/bin/.

# Usage
Default options assume that the alignment is written in the Themisto
format. Add the `--format fulgor` toggle to read in alignments from
Fulgor.

## Read from a file
All calls print the results to cout.

Pack a themisto pseudoalignment file `alignment.txt.gz` containing `2000000` reads pseudoaligned against `1000` reference
sequences, and write the results to `alignment.aln`
```
alignment-writer -f alignment.txt -n 1000 -r 2000000 > alignment.aln
```

Unpack the packed alignment `alignment.aln` with the `-d` toggle
```
alignment-writer -d -f alignment.aln > alignment.tsv
```

The produced `alignment.tsv` file will contain the pseudoalignments
from the original `alignment.txt.gz` file sorted according to the
first column (read_id). This is equivalent to using the
`--sort-output` toggle when running themisto.

## Read from cin
Omitting the `-f` option sets alignment-writer to read input from
cin. This can be used to pack the output from a pseudoaligner without first writing it to disk
```
themisto pseudoalign -q query_reads.fastq -i index --temp-dir tmp | alignment-writer -n <number of reference sequences> -r <number of reads> > alignment.aln
```

## More options
alignment-writer accepts the following flags
```
Usage: alignment-writer -f <input-file>
-f	Pseudoalignment file, packed or unpacked, read from cin if not supplied.
-r	Input reads in .fastq(.gz) format (required for packing).
-l	List containing reference names (1 per line, required for packing).
-d	Unpack pseudoalignment.
--buffer-size	Buffer size for buffered packing (default: 100000
--format	Input file format (one of `themisto` (default), `fulgor`)
--help	Print the help message.
```

# File format
Alignment-writer writes the packed pseudoalignments in chunks
containing ~100000 pseudoalignments by default (can be controlled with
the `--buffer-size` option). The first line in the file contains the
number of reads and the number of reference sequences as
comma-separated integers.

An example packed file with `120000` reads aligned against `8` reference sequences would like the following
```
120000,8
22398
<first chunk: binary representation of the alignments spanning multiple lines>
22297
<second chunk>
8007
<last chunk>
```
The values before each chunk correspond to the size required to read the chunk in as an `unsigned char` array. 

## Reading the file format
Alignment-writer header `unpack.hpp` provides the `Unpack` and
`ParallelUnpack` functions to read the file format into memory.

### Single-threaded
Use the `alignment-writer::Unpack` function to read a file using a single thread.

### Multi-threaded
The `alignment-writer::ParallelUnpack` function can be used to read a
file using multiple threads. The parallelization is implemented using
[OpenMP](https://www.openmp.org/). The caller is responsible for
setting the number of threads using `omp_set_num_threads(NR_THREADS)`
before calling `alignment-writer::ParallelUnpack`.

Using multiple threads incurs a memory overhead that depends on the
`--buffer-size` parameter given when packing the alignment. Larger
buffer sizes result in better multi-thread performance at the cost of
increased memory consumption.

# License
alignment-writer is licensed under the [BSD-3-Clause license](https://opensource.org/licenses/BSD-3-Clause). A copy of the license is supplied with the project, or can alternatively be obtained from [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

## Dependencies
[BitMagic](https://github.com/tlk00/BitMagic) is licensed under the [Apache-2.0 license](https://opensource.org/licenses/Apache-2.0).

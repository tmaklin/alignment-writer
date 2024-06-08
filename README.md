# alignment-writer
Pack/unpack pseudoalignments into/from a compact representation using the
[BitMagic](https://github.com/tlk00/BitMagic) library.

Supported pseudoaligners:
- [Themisto](https://github.com/algbio/themisto)
- [Fulgor](https://github.com/jermp/fulgor)
- [Bifrost](https://github.com/pmelsted/bifrost)
- [Metagraph](https://github.com/ratschlab/metagraph)

Converting .sam files produced by other tools is also supported.

## Installation
### Prebuilt binaries
Download prebuilt binaries for
- Linux x86_64
- macOS arm64
- macOS x86_64

from the [Releases page](https://github.com/tmaklin/alignment-writer/releases).

### Compiling from source
#### Requirements
- C++17 compliant compiler.
- liblzma
- cmake
- git
- OpenMP (optional)

#### How-to
Clone the repository, enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
This will compile the alignment-writer executable in build/bin/.

## Usage
Default options assume that the alignment is written in the Themisto
format. Add the `--format fulgor`, `--format bifrost`, `--format metagraph`
toggle to read in alignments from other pseudoaligners, or `--format sam` to
read alignments from a .sam file.

### Read from a file
#### Pack an alignment
Pack a themisto pseudoalignment file `alignment.txt` created from aligning the reads `reads.fastq.gz` against a Themisto index constructed from the assemblies listed in `target_assemblies_list.txt`
```
alignment-writer -l target_assemblies_list.txt -r reads.fastq.gz alignment.txt
```
This writes the packed results to `alignment.txt.aln`

#### Unpack an alignment
Unpack the packed alignment `alignment.txt.aln` with the `-d` toggle
```
alignment-writer -d alignment.txt.aln
```
The produced `alignment.txt` file will contain the pseudoalignments
from the original `alignment.txt` file sorted according to the
order they occur in the `reads.fastq.gz` file. This is equivalent to using the
`--sort-output` toggle when running themisto.

### Read from terminal
Omitting the positional argument `alignment.txt` option sets alignment-writer to read input from
terminal. This can be used to pack the output from a pseudoaligner without first writing it to disk
```
themisto pseudoalign -q query_reads.fastq -i index --temp-dir tmp | alignment-writer -l target_assemblies_list.txt -r query_reads.fastq > alignment.txt.aln
```

### More options
alignment-writer accepts the following flags
```
  alignment-writer [options] [files]

  -z, --compress         Compress file(s).
  -d, --decompress       Decompress file(s).
  -r, --reads arg        Reads used in the input alignment.
  -l, --target-list arg  File listing the input alignment targets.
      --format arg       Input/output format (default: themisto)
  -k, --keep             Keep input file(s) instead of deleting.
  -f, --force            Force overwrite output file(s).
  -c, --stdout           Write to standard out, keep files.
  -T, --threads arg      Use `arg` threads, 0 = all available. (default: 1)
  -b, --buffer-size arg  Buffer writes every `arg` hits. (default: 256000)
  -h, --help             Print this message and quit.
  -V, --version          Print the version and quit.
```

## File format
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

### Reading the file format
Alignment-writer header `unpack.hpp` provides the `Unpack` and
`ParallelUnpack` functions to read the file format into memory.

#### Single-threaded
Use the `alignment-writer::Unpack` function to read a file using a single thread.

#### Multi-threaded
The `alignment-writer::ParallelUnpack` function can be used to read a
file using multiple threads. The parallelization is implemented using
[OpenMP](https://www.openmp.org/). The caller is responsible for
setting the number of threads using `omp_set_num_threads(NR_THREADS)`
before calling `alignment-writer::ParallelUnpack`.

Using multiple threads incurs a memory overhead that depends on the
`--buffer-size` parameter given when packing the alignment. Larger
buffer sizes result in better multi-thread performance at the cost of
increased memory consumption.

## License
alignment-writer is licensed under the [BSD-3-Clause license](https://opensource.org/licenses/BSD-3-Clause). A copy of the license is supplied with the project, or can alternatively be obtained from [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

### Dependencies
- [BitMagic](https://github.com/tlk00/BitMagic) is licensed under the [Apache-2.0 license](https://opensource.org/licenses/Apache-2.0).
- [bxzstr](https://github.com/tmaklin/bxzstr) is licensed under the [MPL-2.0 license](https://opensource.org/licenses/MPL-2.0).
- [cxxargs](https://github.com/tmaklin/cxxargs) is licensed under the [MPL-2.0 license](https://opensource.org/licenses/MPL-2.0).
- [kseq++](https://github.com/cartoonist/kseqpp) is licensed under the [MIT license](https://opensource.org/licenses/MIT).
- [JSON for Modern C++](https://github.com/nlohmann/json) is licensed under the [MIT license](https://opensource.org/licenses/MIT).

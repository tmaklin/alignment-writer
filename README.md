# alignment-writer
Pack/unpack [Themisto](https://github.com/algbio/themisto)
pseudoalignment files into/from a compact representation using the
[BitMagic](https://github.com/tlk00/BitMagic) library.

# Installation
## Compiling from source
### Requirements
- C++17 compliant compiler.
- cmake
- git

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
All calls print the results to cout.

Pack a themisto pseudoalignment file `alignment.txt` (can be
compressed) produced from a themisto index containing `1000` reference
sequences, and write the results to `alignment.aln`
```
alignment-writer -f alignment.txt -n 1000 > alignment.aln
```

Unpack the packed alignment `alignment.aln` by adding the `-d` toggle and the number of reads `-r` in the alignment
```
alignment-writer -f alignment.aln -d -n 1000 -r 100000 > alignment.tsv
```

The produced `alignment.tsv` file will contain the pseudoalignment
from the original `alignment.txt` file in sorted output (equivalent to
using the `--sort-output` toggle when running themisto).

Intersect packed alignments `alignment_1.aln` and `alignment_2.aln` and write the results to `merged.aln`
```
alignment-writer -f alignment_1.aln --merge alignment_2.aln -n 1000 > merged.aln
```

# License
alignment-writer is licensed under the [BSD-3-Clause license](https://opensource.org/licenses/BSD-3-Clause). A copy of the license is supplied with the project, or can alternatively be obtained from [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

## Dependencies
[BitMagic](https://github.com/tlk00/BitMagic) is licensed under the [Apache-2.0 license](https://opensource.org/licenses/Apache-2.0).

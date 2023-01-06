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
cin. This can be used to pack the output from
[Themisto](https://github.com/algbio/themisto) without first writing it to disk
```
themisto pseudoalign -q query_reads.fastq -i index --temp-dir tmp | alignment-writer -n <number of reference sequences> -r <number of reads> > alignment.aln
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
Alignment-writer `unpack.cpp` provides functions to read the file
format. The following C++ code can be incorporated in other projects
to read the file into memory
```
#include <fstream>
#include <string>
#include <cstddef>

#include "unpack.hpp" // Include the alignment-writer unpacking functions
#include "bm64.h" // BitMagic headers

bm::bvector<> ReadPackedAlignment(std::istream *infile) {
	std::string next_line;

	// Read the number of reads and reference sequences from the first line
	std::getline(*infile, next_line);
	size_t n_reads;
	size_t n_refs;
	alignment_writer::ReadHeader(next_line, &n_reads, &n_refs);

	// Read the chunks into `pseudoalignment`
	bm::bvector<> pseudoalignment(n_reads*n_refs);
	while (std::getline(*infile, next_line)) {
		size_t next_buffer_size = std::stoul(next_line); // Read the size of the next chunk
		alignment_writer::DeserializeBuffer(next_buffer_size, infile, &pseudoalignment); // Read the next chunk
	}

	// Return the `n_reads x n_refs` contiguously stored matrix containing the pseudoalignment.
	// The pseudoalignment for the `n`th read against the `k`th reference sequence is contained
	// at position `n*n_refs + k` assuming indexing starts at 0.
	return pseudoalignment;
}

```

# License
alignment-writer is licensed under the [BSD-3-Clause license](https://opensource.org/licenses/BSD-3-Clause). A copy of the license is supplied with the project, or can alternatively be obtained from [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).

## Dependencies
[BitMagic](https://github.com/tlk00/BitMagic) is licensed under the [Apache-2.0 license](https://opensource.org/licenses/Apache-2.0).

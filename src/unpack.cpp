// alignment-writer: pack/unpack Themisto pseudoalignment files
// https://github.com/tmaklin/alignment-writer
// Copyright (c) 2022 Tommi Mäklin (tommi@maklin.fi)
//
// BSD-3-Clause license
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     (1) Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//
//     (2) Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in
//     the documentation and/or other materials provided with the
//     distribution.
//
//     (3)The name of the author may not be used to
//     endorse or promote products derived from this software without
//     specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "unpack.hpp"

#include <string>
#include <cmath>
#include <sstream>
#include <exception>

#include "bmserial.h"

#include "alignment-writer_openmp_config.hpp"

namespace alignment_writer {
void ReadHeader(const std::string &header_line, size_t *n_reads, size_t *n_refs) {
    std::stringstream header(header_line);
    std::string line;
    std::getline(header, line, ',');
    (*n_reads) = std::stoul(line); // First value is number of reads
    std::getline(header, line, ',');
    (*n_refs) = std::stoul(line); // Second value is number of references
}

void DeserializeBuffer(const size_t buffer_size, std::istream *in, bm::bvector<> *out) {
  // Allocate space for the block
  char* cbuf = new char[buffer_size];

  // Read the next block into buf
  in->read(cbuf, buffer_size);
  unsigned char* buf = reinterpret_cast<unsigned char*>(const_cast<char*>(cbuf));

  // Deserialize block (OR with old data in bits)
  bm::deserialize((*out), buf);

  delete[] cbuf;
}

void Print(std::istream *in, std::ostream *out) {
    // Read size of alignment from the file
    size_t n_reads;
    size_t n_refs;
    std::string line;
    std::getline(*in, line);
    ReadHeader(line, &n_reads, &n_refs);

    // Deserialize the buffer
    bm::bvector<> bits(n_reads*n_refs, bm::BM_GAP);

    while (std::getline(*in, line)) { // Read size of next block
	size_t next_buffer_size = std::stoul(line);
	DeserializeBuffer(next_buffer_size, in, &bits);
    }

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    for (size_t i = 0; i < n_reads; ++i) {
	// Write read id (data compressed with Pack() is sorted so read id is just the iterator id)
	*out << i << ' ';
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    while (*en < i*n_refs + n_refs && en < en_end) {
		*out<< (*en) - i*n_refs << ' ';
		++en;
	    }
	}
	*out << '\n';
    }
    out->flush(); // Flush
}

void StreamingUnpack(std::istream *in, std::ostream *out) {
    // Read size of alignment from the file
    size_t n_reads;
    size_t n_refs;
    std::string line;
    std::getline(*in, line);
    ReadHeader(line, &n_reads, &n_refs);

    while (std::getline(*in, line)) { // Read size of next block
	bm::bvector<> bits;
	size_t next_buffer_size = std::stoul(line);
	DeserializeBuffer(next_buffer_size, in, &bits);

	// Use an enumerator to traverse the pseudoaligned bits
	bm::bvector<>::enumerator en = bits.first();
	bm::bvector<>::enumerator en_end = bits.end();

	bool first = true;
	while (en < en_end) {
	    size_t read_id = std::floor((*en)/n_refs);
	    if (!first) {
		*out << read_id << ' ';
	    }
	    while (std::floor((*en)/n_refs) == read_id && en < en_end) {
		if (first) {
		    *out << read_id << ' ';
		    first = false;
		}
		*out << (*en) - read_id*n_refs << ' ';
		++en;
	    }
	    *out << '\n';
	}
    }

    out->flush(); // Flush
}

bm::bvector<> Unpack(std::istream *infile, size_t *n_reads, size_t *n_refs) {
    std::string next_line;

    // Read the number of reads and reference sequences from the first line
    std::getline(*infile, next_line);
    alignment_writer::ReadHeader(next_line, n_reads, n_refs);

    // Read the chunks into `pseudoalignment`
    bm::bvector<> pseudoalignment((*n_reads)*(*n_refs));
    while (std::getline(*infile, next_line)) {
        size_t next_buffer_size = std::stoul(next_line); // Read the size of the next chunk
	alignment_writer::DeserializeBuffer(next_buffer_size, infile, &pseudoalignment); // Read the next chunk
    }

    // Return the `n_reads x n_refs` contiguously stored matrix containing the pseudoalignment.
    // The pseudoalignment for the `n`th read against the `k`th reference sequence is contained
    // at position `n*n_refs + k` assuming indexing starts at 0.
    return pseudoalignment;
}

bm::bvector<> ParallelUnpack(std::istream *infile, size_t *n_reads, size_t *n_refs) {
#if defined(ALIGNMENTWRITER_OPENMP_SUPPORT) && (ALIGNMENTWRITER_OPENMP_SUPPORT) == 1
    std::string next_line;
    size_t n_threads;

    // Read the number of reads and reference sequences from the first line
    std::getline(*infile, next_line);
    alignment_writer::ReadHeader(next_line, n_reads, n_refs);

    // Get number of threads
    #pragma omp parallel
    {
      n_threads = omp_get_num_threads();
    }

    bm::bvector<> pseudoalignment((*n_reads)*(*n_refs));

    // Read the chunks into `pseudoalignment` in parallel.
    // This loop reads in four chunks at a time using a single thread and then deserializes them in parallel.
    while (std::getline(*infile, next_line)) { // Read size of next block
        std::vector<std::basic_string<unsigned char>> vals;
	for (size_t i = 0; i < n_threads && *infile; ++i) {
	    size_t next_buffer_size = std::stoul(next_line);

	    // Allocate space for the block
	    char* cbuf = new char[next_buffer_size];

	    // Read the next block into buf
	    infile->read(cbuf, next_buffer_size);

	    // Store the block in `vals`
	    unsigned char* buf = reinterpret_cast<unsigned char*>(const_cast<char*>(cbuf));
	    vals.emplace_back(std::basic_string<unsigned char>(buf, next_buffer_size));

	    delete[] cbuf;

	    // Check if there is more to read
	    infile->peek();
	    if (i < n_threads - 1 && *infile) {
	        // If there is more, get the size of the next buffer.
	      std::getline(*infile, next_line);
	    }
	}
	// Deserialize the blocks in parallel (reduce by ORring into output variable)
#pragma omp parallel shared(vals) reduction(bm_bvector_or : pseudoalignment)
	{
	    size_t thread_id = omp_get_thread_num();
	    // If we are at end of file there might be fewer blocks read than there are worker threads
	    if (thread_id < vals.size()) {
	        bm::deserialize(pseudoalignment, vals[thread_id].c_str());
	    }
	}
    }
    // Return the `n_reads x n_refs` contiguously stored matrix containing the pseudoalignment.
    // The pseudoalignment for the `n`th read against the `k`th reference sequence is contained
    // at position `n*n_refs + k` assuming indexing starts at 0.
    return pseudoalignment;
#else
    throw std::runtime_error("Error in alignment-writer::ParallelUnpack: Alignment-writer was not compiled with OpenMP support.");
#endif
}
}

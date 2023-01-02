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

#include "bmserial.h"

namespace alignment_writer {

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

void UnpackBuffered(const size_t &n_refs, const size_t &n_reads, std::istream *in, std::ostream *out) {
    // Deserialize the buffer
    bm::bvector<> bits(n_reads*n_refs, bm::BM_GAP);

    std::string line;
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

void StreamingUnpackBuffered(const size_t &n_refs, const size_t &n_reads, std::istream *in, std::ostream *out) {
    std::string line;
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


void UnpackPlain(const size_t &n_refs, const size_t &n_reads, std::istream *in, std::ostream *out) {
    // Read *in into a buffer
    std::string contents;
    contents.assign(std::istreambuf_iterator<char>(*in),
		    std::istreambuf_iterator<char>());
    const char* cbuf = contents.c_str();
    unsigned char* buf = reinterpret_cast<unsigned char*>(const_cast<char*>(cbuf));

    // Deserialize the buffer
    bm::bvector<> bits(n_reads*n_refs, bm::BM_GAP);
    bm::deserialize(bits, buf);

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
}

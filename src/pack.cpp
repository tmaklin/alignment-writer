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
#include "pack.hpp"

#include <sstream>
#include <exception>

#include "bm64.h"
#include "bmserial.h"

namespace alignment_writer {
void CheckInput(const size_t n_refs, const size_t n_reads) {
    size_t aln_size = (size_t)(n_reads * n_refs);
    if (aln_size > (size_t)2^47) {
	throw std::length_error("Input size exceeds maximum capacity (number of reads x number of references > 2^(48 - 1)).");
    }
}

void WriteHeader(const size_t n_refs, const size_t n_reads, std::ostream *out) {
    // Write the header line of the packed format
    *out << n_reads << ',' << n_refs << '\n';
}

void WriteBuffer(const bm::bvector<> &bits, bm::serializer<bm::bvector<>> &bvs, std::ostream *out) {
    // Use serialization buffer class (automatic RAI, freed on destruction)
    bm::serializer<bm::bvector<>>::buffer sbuf;
    bvs.serialize(bits, sbuf);

    //  Write to *out
    auto sz = sbuf.size();
    *out << sz << '\n';
    unsigned char* buf = sbuf.data();
    for (size_t i = 0; i < sz; ++i) {
	*out << buf[i];
    }
}

void BufferedPack(const size_t n_refs, const size_t n_reads, const size_t &buffer_size, std::istream *in, std::ostream *out) {
    // Buffered read + packing from a stream
    // Write info about the pseudoalignment
    CheckInput(n_refs, n_reads);
    WriteHeader(n_refs, n_reads, out);

    // Next settings provide the lowest size (see BitMagic documentation/examples)
    bm::serializer<bm::bvector<>> bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);

    bm::bvector<> bits;
    bits.set_new_blocks_strat(bm::BM_GAP);
    bm::bvector<>::bulk_insert_iterator it(bits);

    size_t n_in_buffer = 0;
    std::string line;
    while (std::getline(*in, line)) {
	std::stringstream stream(line);
	std::string part;
	std::getline(stream, part, ' ');
	size_t read_id = std::stoul(part); // First column is a numerical ID for the read
	while(std::getline(stream, part, ' ')) {
	    // Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	    it = read_id*n_refs + std::stoul(part);
	    ++n_in_buffer;
	}
	if (n_in_buffer > buffer_size) {
  	    // Force flush on the inserter to ensure everything is saved
	    it.flush();

	    WriteBuffer(bits, bvs, out);
	    bits.clear(true);
	    bits.set_new_blocks_strat(bm::BM_GAP);
	    n_in_buffer = 0;
	}
    }

    // Write the remaining bits
    it.flush();
    WriteBuffer(bits, bvs, out);
    out->flush(); // Flush
}

void Pack(const bm::bvector<> &bits, const size_t n_refs, const size_t n_reads, std::ostream *out) {
    // Pack a pseudoalignment that has been stored in memory
    // Write info about the pseudoalignment
    CheckInput(n_refs, n_reads);
    WriteHeader(n_refs, n_reads, out);

    // Next settings provide the lowest size (see BitMagic documentation/examples)
    bm::serializer<bm::bvector<>> bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);

    WriteBuffer(bits, bvs, out);
    out->flush(); // Flush
}
}

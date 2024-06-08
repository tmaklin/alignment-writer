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
#ifndef ALIGNMENT_WRITER_PACK_HPP
#define ALIGNMENT_WRITER_PACK_HPP

#include <cstddef>
#include <istream>
#include <ostream>
#include <unordered_map>

#include "bm64.h"
#include "bmserial.h"

#include "parser.hpp"

namespace alignment_writer {
// Pack a pseudoalignment already in memory as a BitMagic vector
void Pack(const bm::bvector<> &bits, const size_t n_refs, const size_t n_reads, std::ostream *out);

// Pack a pseudoalignment stored in some iterable container
template <typename T>
inline void PackContainer(T &bits, const size_t n_refs, const size_t n_reads, std::ostream *out) {
    bm::bvector<> bm_bits;
    bm_bits.set_new_blocks_strat(bm::BM_GAP);
    bm::bvector<>::bulk_insert_iterator bm_it(bm_bits);

    bm_bits.resize((size_t)(n_reads * n_refs));
    typename T::iterator iter = bits.begin();
    size_t pos = 0;
    while (iter < bits.end()) {
	if (*iter) {
	    bm_it = pos;
	}
	++iter;
	++pos;
    }
    Pack(bits, n_refs, n_reads, out);
}

// Write a single chunk
void WriteBlock(const bm::bvector<> &bits, const std::string &query_info,
		std::ostream *out, bm::serializer<bm::bvector<>> *bvs);
// Write a single chunk stored in some iterable container
template <typename T>
void WriteContainerBlock(T &bits, const std::string &query_info, std::ostream *out) {
    bm::bvector<> bm_bits;
    bm_bits.set_new_blocks_strat(bm::BM_GAP);
    bm::bvector<>::bulk_insert_iterator bm_it(bm_bits);

    bm_bits.resize((size_t)(bits.end() - bits.begin()));
    typename T::iterator iter = bits.begin();
    size_t pos = 0;
    while (iter < bits.end()) {
	if (*iter) {
	    bm_it = pos;
	}
	++iter;
	++pos;
    }
    bm::serializer<bm::bvector<>> bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);
    WriteBlock(bm_bits, query_info, out, &bvs);
}

// Buffered read of a pseudoalignment from a stream and packing
void BufferedPack(const Format &format,
		  const std::unordered_map<std::string, size_t> &query_to_position,
		  const std::unordered_map<std::string, size_t> &ref_to_position,
		  const size_t &buffer_size, std::istream *in, std::ostream *out);
}

#endif

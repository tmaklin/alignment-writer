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

#include <lzma.h>

#include <sstream>
#include <exception>
#include <functional>

#include "bm64.h"
#include "bmserial.h"
#include "bxzstr.hpp"

namespace alignment_writer {
void CheckInput(const size_t n_refs, const size_t n_reads) {
    size_t aln_size = (size_t)(n_reads * n_refs);
    if (aln_size > (size_t)140737488355328) {
	throw std::length_error("Input size exceeds maximum capacity (number of reads x number of references > 2^(48 - 1)).");
    }
}

void WriteHeader(const std::unordered_map<std::string, size_t> &query_to_position,
		 const std::unordered_map<std::string, size_t> &ref_to_position,
		 std::ostream *out) {
    // Write the header line of the packed format
    size_t n_reads = query_to_position.size();
    size_t n_refs = ref_to_position.size();

    std::stringbuf buf;
    bxz::ostream lzma(&buf, bxz::lzma, 1);
    lzma << "{";
    lzma << "\"n_queries\":" << n_reads << ',' << "\"n_targets\":" << n_refs;
    lzma << ',';
    lzma << "\"targets\":[";
    size_t n_written = 0;
    for (auto kv : ref_to_position) {
	lzma << "{\"target\":\"" << kv.first << '"' << ',' << "\"pos\":" << kv.second << '}'<< ((n_written == n_refs - 1) ? ']' : ',');
	++n_written;
    }
    lzma << '}';
    lzma.flush();
    *out << buf.str().size() << '\n';
    *out << buf.str();
    out->flush();
}

void WriteBufferHeader(const std::unordered_map<size_t, std::string> &position_to_query,
		       const std::vector<size_t> &queries_in_buffer,
		       const bm::serializer<bm::bvector<>>::buffer &sbuf, std::ostream *out) {
    // Write the header line of the packed format
    size_t n_reads = queries_in_buffer.size();

    auto sz = sbuf.size();

    std::stringbuf buf;
    bxz::ostream lzma(&buf, bxz::lzma, 1);
    lzma << "{";
    lzma << "\"queries\":[";
    size_t n_written = 0;
    for (auto query_pos : queries_in_buffer) {
	lzma << "{\"query\":\"" << position_to_query.at(query_pos) << '"' << ',' << "\"pos\":" << query_pos << '}'<< ((n_written == n_reads - 1) ? ']' : ',');
	++n_written;
    }
    lzma << '}';
    lzma.flush();
    *out << '{' << "\"header_size\":" << buf.str().size() << ',' << "\"block_size\":" << sz << '}' << '\n';
    *out << buf.str();
    out->flush();
}

void WriteBuffer(bm::serializer<bm::bvector<>>::buffer &sbuf, std::ostream *out) {
    //  Write to *out
    auto sz = sbuf.size();
    unsigned char* buf = sbuf.data();
    for (size_t i = 0; i < sz; ++i) {
	*out << buf[i];
    }
    out->flush();
}

void WriteBlock(const bm::bvector<> &bits, const std::unordered_map<size_t, std::string> &pos_to_query,
		const std::vector<size_t> &reads_in_buffer, std::ostream *out,
		bm::serializer<bm::bvector<>> *bvs) {

    // Use serialization buffer class (automatic RAI, freed on destruction)
    bm::serializer<bm::bvector<>>::buffer sbuf;
    bvs->serialize(bits, sbuf);

    WriteBufferHeader(pos_to_query, reads_in_buffer, sbuf, out);
    WriteBuffer(sbuf, out);
}

void BufferedPack(const Format &format, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, const size_t &buffer_size, std::istream *in, std::ostream *out) {
    // Buffered read + packing from a stream
    // Write info about the pseudoalignment
    size_t n_reads = query_to_position.size();
    size_t n_refs = ref_to_position.size();
    CheckInput(n_refs, n_reads);
    WriteHeader(query_to_position, ref_to_position, out);

    // Next settings provide the lowest size (see BitMagic documentation/examples)
    bm::serializer<bm::bvector<>> bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);

    bm::bvector<> bits;
    bits.set_new_blocks_strat(bm::BM_GAP);
    bm::bvector<>::bulk_insert_iterator it(bits);

    std::unordered_map<size_t, std::string> pos_to_query;
    for (auto kv : query_to_position) {
	pos_to_query.insert(std::make_pair(kv.second, kv.first));
    }

    std::function<size_t(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::vector<size_t> *reads_in_buffer)> parser;
    if (format == themisto) {
	parser = ThemistoParser;
    } else if (format == fulgor) {
	parser = FulgorParser;
    } else if (format == bifrost) {
	// Consume first line in Bifrost format (contains ref names)
	std::string header;
	std::getline(*in, header);
	parser = BifrostParser;
    } else if (format == metagraph) {
	parser = MetagraphParser;
    } else if (format == sam) {
	// Consume the header lines from the SAM format
	std::string header;
	while (in->peek() == '@') {
	    std::getline(*in, header);
	}
	parser = SAMParser;
    } else {
	throw std::runtime_error("Unrecognized input format.");
    }

    size_t n_in_buffer = 0;
    std::vector<size_t> reads_in_buffer;
    std::string line;
    while (std::getline(*in, line)) {
	// Parse the line
	n_in_buffer += parser(line, query_to_position, ref_to_position, &it, &reads_in_buffer);

	if (n_in_buffer > buffer_size) {
  	    // Force flush on the inserter to ensure everything is saved
	    it.flush();
	    WriteBlock(bits, pos_to_query, reads_in_buffer, out, &bvs);

	    bits.clear(true);
	    bits.set_new_blocks_strat(bm::BM_GAP);
	    reads_in_buffer.clear();
	    n_in_buffer = 0;
	}
    }

    if (reads_in_buffer.size() > 0) {
	// Write the remaining bits
	it.flush();
	WriteBlock(bits, pos_to_query, reads_in_buffer, out, &bvs);
    }
}

void Pack(const bm::bvector<> &bits, const std::unordered_map<std::string, size_t> &query_to_position,
	  const std::unordered_map<std::string, size_t> &ref_to_position,
	  const size_t n_refs, const size_t n_reads, std::ostream *out) {
    // Pack a pseudoalignment that has been stored in memory
    // Write info about the pseudoalignment
    CheckInput(n_refs, n_reads);
    WriteHeader(query_to_position, ref_to_position, out);

    // Next settings provide the lowest size (see BitMagic documentation/examples)
    bm::serializer<bm::bvector<>> bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);

    std::vector<size_t> queries_in_buffer(query_to_position.size());
    for (size_t i = 0; i < query_to_position.size(); ++i) {
	queries_in_buffer[i] = i;
    }

    std::unordered_map<size_t, std::string> pos_to_query;
    for (auto kv : query_to_position) {
	pos_to_query.insert(std::make_pair(kv.second, kv.first));
    }

    WriteBlock(bits, pos_to_query, queries_in_buffer, out, &bvs);
}
}

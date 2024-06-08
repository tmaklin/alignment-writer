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

#include "nlohmann/json.hpp"
#include "bmserial.h"
#include "bxzstr.hpp"

#include "alignment-writer_openmp_config.hpp"

namespace alignment_writer {
void ReadHeader(std::istream *in, size_t *n_reads, size_t *n_refs) {
    std::string first_line;
    std::getline(*in, first_line);
    size_t buffer_size = std::stoul(first_line);

    std::stringbuf buffer;
    for (size_t i = 0; i < buffer_size; ++i) {
	buffer.sputc(in->get());
    }

    bxz::istream instr(&buffer);
    nlohmann::json_abi_v3_11_3::json header = nlohmann::json_abi_v3_11_3::json::parse(instr);
    *n_reads = header["n_queries"];
    *n_refs = header["n_targets"];
}

nlohmann::json_abi_v3_11_3::json DeserializeBlockHeader(std::stringbuf &buffer) {
    bxz::istream instr(&buffer);
    return nlohmann::json_abi_v3_11_3::json::parse(instr);
}

std::stringbuf ReadBlockHeader(const std::string &line, std::istream *in, size_t *block_size) {
    auto header_data = nlohmann::json_abi_v3_11_3::json::parse(line);
    size_t header_buffer_size = (size_t)header_data["header_size"];

    std::stringbuf buffer;
    for (size_t i = 0; i < header_buffer_size; ++i) {
	buffer.sputc(in->get());
    }

    *block_size = (size_t)header_data["block_size"];
    return buffer; // Use DeserializeBlockHeader to read contents is needed
}


std::basic_string<unsigned char> ReadBytes(const size_t bytes, std::istream *in) {
    // Allocate space for the block
    std::basic_string<char> cbuf(bytes, '=');
    // char* cbuf = new char[bytes];

    // Read the next block into buf
    in->read(cbuf.data(), bytes);

    // Store the block in `vals`
    unsigned char* buf = reinterpret_cast<unsigned char*>(cbuf.data());

    return std::basic_string<unsigned char>(buf, bytes);
}

void ReadBlock(const std::string &line, std::istream *in, bm::bvector<> *bits_out) {
    size_t block_size = 0;
    ReadBlockHeader(line, in, &block_size);
    std::basic_string<unsigned char> buf = ReadBytes(block_size, in);

    // Deserialize block (OR with old data in bits)
    bm::deserialize((*bits_out), buf.data());
}

void Print(std::istream *in, std::ostream *out) {
    // Read size of alignment from the file
    size_t n_reads;
    size_t n_refs;
    ReadHeader(in, &n_reads, &n_refs);

    // Deserialize the buffer
    bm::bvector<> bits(n_reads*n_refs, bm::BM_GAP);

    std::string line;
    while (std::getline(*in, line)) {
	ReadBlock(line, in, &bits);
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
    ReadHeader(in, &n_reads, &n_refs);

    std::string line;
    while (std::getline(*in, line)) { // Read size of next block
	bm::bvector<> bits;
	ReadBlock(line, in, &bits);

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

bm::bvector<> Unpack(std::istream *in, size_t *n_reads, size_t *n_refs) {
    // Read the number of reads and reference sequences from the first line
    alignment_writer::ReadHeader(in, n_reads, n_refs);

    // Read the chunks into `pseudoalignment`
    bm::bvector<> bits((*n_reads)*(*n_refs));
    std::string line;
    while (std::getline(*in, line)) {
	ReadBlock(line, in, &bits);
    }

    // Return the `n_reads x n_refs` contiguously stored matrix containing the pseudoalignment.
    // The pseudoalignment for the `n`th read against the `k`th reference sequence is contained
    // at position `n*n_refs + k` assuming indexing starts at 0.
    return bits;
}

void ParallelUnpackData(std::istream *infile, bm::bvector<> &pseudoalignment) {
    // Read the chunks into `pseudoalignment` in parallel.
#if defined(ALIGNMENTWRITER_OPENMP_SUPPORT) && (ALIGNMENTWRITER_OPENMP_SUPPORT) == 1
    // Get number of threads
    size_t n_threads;
    #pragma omp parallel
    {
      n_threads = omp_get_num_threads();
    }

    std::string next_line;
    // This loop reads in four chunks at a time using a single thread and then deserializes them in parallel.
    while (std::getline(*infile, next_line)) { // Read size of next block
        std::vector<std::basic_string<unsigned char>> vals;
	for (size_t i = 0; i < n_threads && *infile; ++i) {

	    // Read the block header
	    size_t block_data_size = 0;
	    ReadBlockHeader(next_line, infile, &block_data_size);

	    // Read the block data
	    std::basic_string<unsigned char> block_data = std::move(ReadBytes(block_data_size, infile));

	    // Store data in `vals`
	    vals.emplace_back(block_data);

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
#else
    throw std::runtime_error("Error in alignment-writer::ParallelUnpack: Alignment-writer was not compiled with OpenMP support.");
#endif
}

bm::bvector<> ParallelUnpack(std::istream *infile, size_t *n_reads, size_t *n_refs) {
    // Read the number of reads and reference sequences from the first line
    alignment_writer::ReadHeader(infile, n_reads, n_refs);

    bm::bvector<> pseudoalignment((*n_reads)*(*n_refs));
    ParallelUnpackData(infile, pseudoalignment);

    // Return the `n_reads x n_refs` contiguously stored matrix containing the pseudoalignment.
    // The pseudoalignment for the `n`th read against the `k`th reference sequence is contained
    // at position `n*n_refs + k` assuming indexing starts at 0.
    return pseudoalignment;
}
}

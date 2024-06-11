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

#include "Alignment.hpp"
#include "printer.hpp"

namespace alignment_writer {
bool check_xz_header(std::istream *in, std::stringbuf *out) {
    (*out) = std::stringbuf();
    for (size_t i = 0; i < 6; ++i) {
	out->sputc(in->get());
    }
    unsigned char b0 = out->str()[0];
    unsigned char b1 = out->str()[1];
    unsigned char b2 = out->str()[2];
    unsigned char b3 = out->str()[3];
    unsigned char b4 = out->str()[4];
    unsigned char b5 = out->str()[5];

    bool is_xz_header = (b0 == 0xFD && b1 == 0x37 && b2 == 0x7A
			&& b3 == 0x58 && b4 == 0x5A && b5 == 0x00);
    return is_xz_header;
}

bool read_until_xz_end(std::istream *in, std::stringbuf *out) {
    bool stream_end = false;
    while (!stream_end && in->good()) {
	unsigned char next = in->get();
	out->sputc(next);
	if (next == 0x59) {
	    // Possibly found stream end
	    unsigned char nextnext = in->get();
	    out->sputc(nextnext);
	    stream_end = (nextnext == 0x5A);
	}
    }
    return stream_end && in->good();
}

nlohmann::json_abi_v3_11_3::json ReadHeader(std::istream *in) {
    std::stringbuf buffer;
    if (!check_xz_header(in, &buffer)) {
	throw std::runtime_error("Input file does not start with a .xz stream header.");
    }

    if (!read_until_xz_end(in, &buffer)) {
	throw std::runtime_error("Unexpected end of input.");
    }

    bxz::istream instr(&buffer);
    const nlohmann::json_abi_v3_11_3::json &header = nlohmann::json_abi_v3_11_3::json::parse(instr);

    return header;
}

nlohmann::json_abi_v3_11_3::json DeserializeBlockHeader(std::stringbuf &buffer) {
    bxz::istream instr(&buffer);
    return nlohmann::json_abi_v3_11_3::json::parse(instr);
}

std::stringbuf ReadBlockHeader(std::istream *in, size_t *block_size) {
    const auto &header_data = ReadHeader(in);
    size_t header_buffer_size = (size_t)header_data["header_size"];
    *block_size = (size_t)header_data["block_size"];

    std::stringbuf ret_buffer;
    for (size_t i = 0; i < header_buffer_size; ++i) {
	ret_buffer.sputc(in->get());
    }

    return ret_buffer; // Use DeserializeBlockHeader to read contents is needed
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

std::stringbuf ReadBlock(std::istream *in, bm::bvector<> *bits_out) {
    size_t block_size = 0;
    std::stringbuf block_header = std::move(ReadBlockHeader(in, &block_size));
    std::basic_string<unsigned char> buf = ReadBytes(block_size, in);

    // Deserialize block (OR with old data in bits)
    bm::deserialize((*bits_out), buf.data());

    return block_header;
}

Alignment DecompressStream(std::istream *in) {
    // Read size of alignment from the file
    const nlohmann::json_abi_v3_11_3::json &file_header = ReadHeader(in);
    size_t n_reads = file_header["n_queries"];
    size_t n_refs = file_header["n_targets"];

    // Consume stream
    std::string line;
    bool first = true;
    Alignment bits(file_header);
    while (in->good() && in->peek() != EOF) {
	std::stringbuf header = std::move(ReadBlock(in, &bits));
	auto block_headers = DeserializeBlockHeader(header);
	bits.annotate(block_headers);
    }
    return bits;
}

void Print(const Format &format, std::istream *in, std::ostream *out) {
    // Deserialize the file
    const Alignment &alignment = DecompressStream(in);

    // Initialize formatter
    Printer printer(format);

    // Print results
    const std::stringbuf &ret = printer.format(alignment);
    *out << ret.str();
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
	    ReadBlockHeader(infile, &block_data_size);

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
    const auto &file_header = alignment_writer::ReadHeader(infile);

    Alignment pseudoalignment(file_header);
    ParallelUnpackData(infile, pseudoalignment);

    // Return the `n_reads x n_refs` contiguously stored matrix containing the pseudoalignment.
    // The pseudoalignment for the `n`th read against the `k`th reference sequence is contained
    // at position `n*n_refs + k` assuming indexing starts at 0.
    return pseudoalignment.raw();
}
}

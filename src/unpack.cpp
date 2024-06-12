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
#include "BS_thread_pool.hpp"

#include "Alignment.hpp"
#include "printer.hpp"

namespace alignment_writer {
bool read_xz_header(std::istream *in, std::basic_string<unsigned char> *out) {
    for (size_t i = 0; i < 6 && in->good() && in->peek() != EOF; ++i) {
	out->push_back(in->get());
    }
    if (out->size() < 6) {
	return 0;
    }


    bool is_xz_header = ((*out)[0] == 0xFD);
    is_xz_header &= ((*out)[1] == 0x37);
    is_xz_header &= ((*out)[2] == 0x7A);
    is_xz_header &= ((*out)[3] == 0x58);
    is_xz_header &= ((*out)[4] == 0x5A);
    is_xz_header &= ((*out)[5] == 0x00);

    return is_xz_header;
}

bool read_until_xz_end(std::istream *in, std::basic_string<unsigned char> *out) {
    bool stream_end = false;
    while (!stream_end && in->good()) {
	unsigned char next = in->get();
	out->push_back(next);
	if (next == 0x59) {
	    // Possibly found stream end
	    unsigned char nextnext = in->get();
	    out->push_back(nextnext);
	    stream_end = (nextnext == 0x5A);
	}
    }
    return stream_end && in->good();
}

std::basic_string<char> decompress_xz(const std::basic_string<unsigned char> &buffer) {
    std::basic_stringstream<unsigned char> str(buffer);
    bxz::istream instr(reinterpret_cast<std::stringbuf*>(str.rdbuf()));
    std::string out;
    std::string read;
    while (instr.good() && instr.peek() != EOF) {
	out += instr.get();
    }

    return out;
}

std::basic_string<char> ReadHeader(std::istream *in) {
    std::basic_string<unsigned char> buffer;
    if (!read_xz_header(in, &buffer)) {
	throw std::runtime_error("Input file does not start with a .xz stream header.");
    }

    if (!read_until_xz_end(in, &buffer)) {
	throw std::runtime_error("Unexpected end of input.");
    }

    std::basic_string<char> header = decompress_xz(buffer);

    return header;
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

std::basic_string<unsigned char> ReadBlockHeader(std::istream *in, size_t *block_size) {
    const json &header_data = json::parse(ReadHeader(in));
    size_t header_buffer_size = (size_t)header_data["header_size"];
    *block_size = (size_t)header_data["block_size"];

    std::basic_string<unsigned char> ret = ReadBytes(header_buffer_size, in);

    return ret; // Use DecompressXZ to read contents is needed
}

void ReadBlock(std::istream *in, std::basic_string<unsigned char> *block_header, std::basic_string<unsigned char> *block) {
    size_t block_size = 0;
    (*block_header) = std::move(ReadBlockHeader(in, &block_size));
    (*block) = std::move(ReadBytes(block_size, in));
}

void DecompressBlock(const std::basic_string<unsigned char> &block_header_bytes, const std::basic_string<unsigned char> &block_bytes, Alignment *bits) {
    const json &block_headers = json::parse(decompress_xz(block_header_bytes));
    bits->annotate(block_headers["queries"]);
    bm::deserialize(*bits, block_bytes.data());
}

Alignment DecompressBlock2(const json &file_header, const std::basic_string<unsigned char> &block_header_bytes, const std::basic_string<unsigned char> &block_bytes) {
    Alignment bits(file_header);
    const json &block_headers = json::parse(decompress_xz(block_header_bytes));
    bits.annotate(block_headers["queries"]);
    bm::deserialize(bits, block_bytes.data());
    return bits;
}

void Print(const Format &format, std::istream *in, std::ostream *out, size_t n_threads) {
    // Initialize formatter
    Printer printer(format);

    // Read size of alignment from the file
    const json &file_header = json::parse(ReadHeader(in));
    size_t n_reads = file_header["n_queries"];
    size_t n_refs = file_header["n_targets"];

    BS::thread_pool pool(n_threads);

    std::vector<std::basic_string<unsigned char>> block_header_bytes;
    std::vector<std::basic_string<unsigned char>> block_bytes;
    std::vector<std::future<Alignment>> futures;
    // Stream the output
    while (in->good() && in->peek() != EOF) {
	if (block_bytes.size() < n_threads) {
	    // Read blocks until each thread has something to process
	    block_header_bytes.emplace_back(std::basic_string<unsigned char>());
	    block_bytes.emplace_back(std::basic_string<unsigned char>());
	    ReadBlock(in, &block_header_bytes.back(), &block_bytes.back());
	    futures.emplace_back(pool.submit(DecompressBlock2, file_header, block_header_bytes.back(), block_bytes.back()));
	} else {
	    // Process all blocks
	    for (size_t i = 0; i < futures.size(); ++i) {
		// Deserialize and print contents
		const Alignment &alignment = futures[i].get();
		const std::stringbuf &ret = printer.format(alignment);
		*out << ret.str();
	    }
	    // Clear thread buffers
	    block_header_bytes.clear();
	    block_bytes.clear();
	    futures.clear();
	}
    }
    // Process remaining blocks
    if (futures.size() > 0) {
	for (size_t i = 0; i < futures.size(); ++i) {
	    const Alignment &alignment = futures[i].get();
	    const std::stringbuf &ret = printer.format(alignment);
	    *out << ret.str();
	}
    }
}

Alignment Read(std::istream *in, size_t n_threads) {
    // Read size of alignment from the file
    const json &file_header = json::parse(ReadHeader(in));
    size_t n_reads = file_header["n_queries"];
    size_t n_refs = file_header["n_targets"];

    Alignment alignment(file_header);
    BS::thread_pool pool(n_threads);

    std::vector<std::basic_string<unsigned char>> block_header_bytes;
    std::vector<std::basic_string<unsigned char>> block_bytes;
    std::vector<std::future<Alignment>> futures;
    // Stream the output
    while (in->good() && in->peek() != EOF) {
	if (block_bytes.size() < n_threads) {
	    // Read blocks until each thread has something to process
	    block_header_bytes.emplace_back(std::basic_string<unsigned char>());
	    block_bytes.emplace_back(std::basic_string<unsigned char>());
	    ReadBlock(in, &block_header_bytes.back(), &block_bytes.back());
	    futures.emplace_back(pool.submit(DecompressBlock2, file_header, block_header_bytes.back(), block_bytes.back()));
	} else {
	    // Process all blocks
	    for (size_t i = 0; i < futures.size(); ++i) {
		// Deserialize and print contents
		const Alignment &got = futures[i].get();
		alignment.bit_or(got);
		alignment.annotate(got.annotation());
	    }
	    // Clear thread buffers
	    block_header_bytes.clear();
	    block_bytes.clear();
	    futures.clear();
	}
    }
    // Process remaining blocks
    if (futures.size() > 0) {
	for (size_t i = 0; i < futures.size(); ++i) {
	    const Alignment &got = futures[i].get();
	    alignment.bit_or(got);
	    alignment.annotate(got.annotation());
	}
    }

    return alignment;
}
}

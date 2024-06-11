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

std::basic_string<unsigned char> ReadBlock(std::istream *in, std::stringbuf *block_header) {
    size_t block_size = 0;
    (*block_header) = std::move(ReadBlockHeader(in, &block_size));
    const std::basic_string<unsigned char> &buf = ReadBytes(block_size, in);
    return buf;
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
    std::vector<std::basic_string<unsigned char>> buffers;
    while (in->good() && in->peek() != EOF) {
	std::stringbuf block_header_ser;
	buffers.emplace_back(ReadBlock(in, &block_header_ser));
	auto block_headers = DeserializeBlockHeader(block_header_ser);
	bits.annotate(block_headers);
    }

    // Deserialize
    for (auto buffer : buffers) {
	bm::deserialize(bits, buffer.data());
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
}

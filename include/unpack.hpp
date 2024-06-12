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
#ifndef ALIGNMENT_WRITER_UNPACK_HPP
#define ALIGNMENT_WRITER_UNPACK_HPP

#include <istream>
#include <ostream>
#include <cstddef>
#include <string>

#include "bm64.h"

#include "Alignment.hpp"
#include "parser.hpp"

namespace alignment_writer {

// Read compressed files
//
// Print compressed data to std::ostream
void Print(const Format &format, std::istream *in, std::ostream *out, size_t n_threads=1);
//
// Read compressed data into memory
Alignment Read(std::istream *in, size_t n_threads=1);


// Read parts (blocks) from compressed files
//
// Read and decompress file or block header
std::basic_string<char> ReadHeader(std::istream *in);
//
// Read the bytes comprising a block
void ReadBlock(std::istream *in, std::basic_string<unsigned char> *block_header, std::basic_string<unsigned char> *block);
//
// Decompress the bytes comprising a block by OR:ing the contents into `*bits`
void DecompressBlock(const std::basic_string<unsigned char> &block_header_bytes, const std::basic_string<unsigned char> &block_bytes, Alignment *bits);
//
// Decompress the bytes comprising a block into a returned `Alignment`
Alignment DecompressBlock2(const json &file_header, const std::basic_string<unsigned char> &block_header_bytes, const std::basic_string<unsigned char> &block_bytes);

}

#endif

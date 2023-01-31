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

#include <cstddef>
#include <istream>
#include <ostream>
#include <vector>

#include "bm64.h"

namespace alignment_writer {
// Print data that has been written using BufferedPack
void Print(std::istream *in, std::ostream *out);
void StreamingPrint(std::istream *in, std::ostream *out);

// Read in pseudoalignment data written using BufferedPack
bm::bvector<> Unpack(std::istream *infile, size_t *n_reads, size_t *n_refs);
// Parallel read
bm::bvector<> ParallelUnpack(std::istream *infile, size_t *n_reads, size_t *n_refs);

// Deserialize one section of data written with BufferedPack
void DeserializeBuffer(const size_t buffer_size, std::istream *in, bm::bvector<> *out);

// Function for reading the header line of the alignment file
void ReadHeader(const std::string &header_line, size_t *n_reads, size_t *n_refs);
}

#endif

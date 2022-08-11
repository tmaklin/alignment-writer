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
#include <string>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <exception>
#include <memory>

#include "bm64.h"
#include "bmserial.h"
#include "cxxargs.hpp"
#include "bxzstr.hpp"

#include "version.h"

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_short_argument<std::string>('f', "Pseudoalignment file, packed or unpacked, read from cin if not supplied.", "");
  args.add_short_argument<bool>('d', "Unpack pseudoalignment.", false);
  args.add_short_argument<size_t>('n', "Number of reference sequences in the pseudoalignment (required).");
  args.add_short_argument<size_t>('r', "Number of reads in the pseudoalignment (required for unpacking).");
  args.add_long_argument<size_t>("buffer-size", "Buffer size for buffered packing (default: 100000", (size_t)100000);
  if (!CmdOptionPresent(argv, argv+argc, "-d")) {
      args.set_not_required('r');
  }
  args.add_long_argument<bool>("help", "Print the help message.", false);
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
      std::cerr << "\n" + args.help() << '\n' << '\n';
  }
  args.parse(argc, argv);
}

void Dump(bm::bvector<>::bulk_insert_iterator &it, bm::bvector<> &bits, bm::serializer<bm::bvector<>> &bvs, std::ostream *out) {
    // Force flush on the inserter to ensure everything is saved
    it.flush();

    // Use serialization buffer class (automatic RAI, freed on destruction)
    bm::serializer<bm::bvector<>>::buffer sbuf;
    bvs.optimize_serialize_destroy(bits, sbuf);

    //  Write to *out
    auto sz = sbuf.size();
    *out << sz << std::endl;
    unsigned char* buf = sbuf.data();
    for (size_t i = 0; i < sz; ++i) {
	*out << buf[i];
    }
}

void BufferedPack(const size_t &n_refs, const size_t &buffer_size, std::istream *in, std::ostream *out) {
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
	    Dump(it, bits, bvs, out);
	    bits.clear(true);
	    bits.set_new_blocks_strat(bm::BM_GAP);
	    n_in_buffer = 0;
	}
    }

    Dump(it, bits, bvs, out);
    out->flush(); // Flush
}

void UnpackBuffered(const size_t &n_refs, const size_t &n_reads, std::istream *in, std::ostream *out) {
    // Deserialize the buffer
    bm::bvector<> bits(n_reads*n_refs, bm::BM_GAP);

    std::string line;
    while (std::getline(*in, line)) { // Read size of next block
	size_t next_buffer_size = std::stoul(line);

	// Allocate space for the block
	char* cbuf = new char[next_buffer_size];

	// Read the next block into buf
	in->read(cbuf, next_buffer_size);
	unsigned char* buf = reinterpret_cast<unsigned char*>(const_cast<char*>(cbuf));

	// Deserialize block (OR with old data in bits)
	bm::deserialize(bits, buf);

	delete[] cbuf;
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

void Pack(const size_t &n_refs, std::istream *in, std::ostream *out) {
    bm::bvector<> bits;
    bits.set_new_blocks_strat(bm::BM_GAP);
    bm::bvector<>::bulk_insert_iterator it(bits);

    size_t n_reads = 0; // Count the reads
    std::string line;
    while (std::getline(*in, line)) {
	std::stringstream stream(line);
	std::string part;
	std::getline(stream, part, ' ');
	size_t read_id = std::stoul(part); // First column is a numerical ID for the read
	while(std::getline(stream, part, ' ')) {
	    // Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	    it = read_id*n_refs + std::stoul(part);
	}
	++n_reads;
    }

    // Force flush on the inserter to ensure everything is saved
    it.flush();

    // Compress the pseudoalignment matrix
    bits.resize(n_reads*n_refs); // Final size is now known
    bits.optimize();
    bits.freeze();

    // Next settings provide the lowest size (see BitMagic documentation/examples)
    bm::serializer<bm::bvector<>> bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);

    // Use serialization buffer class (automatic RAI, freed on destruction)
    bm::serializer<bm::bvector<>>::buffer sbuf;
    bvs.serialize(bits, sbuf);

    //  Write to *out
    auto sz = sbuf.size();
    unsigned char* buf = sbuf.data();
    for (size_t i = 0; i < sz; ++i) {
	*out << buf[i];
    }
    out->flush(); // Flush
}

void Unpack(const size_t &n_refs, const size_t &n_reads, std::istream *in, std::ostream *out) {
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

int main(int argc, char* argv[]) {
    cxxargs::Arguments args("alignment-writer-" + std::string(ALIGNMENT_WRITER_BUILD_VERSION), "Usage: alignment-writer -f <input-file>");
    try {
	parse_args(argc, argv, args);
    } catch (std::exception &e) {
	std::cerr << "Parsing arguments failed:\n"
		  << std::string("\t") + std::string(e.what()) + "\n"
		  << "\trun alignment-writer with the --help option for usage instructions.\n";
	std::cerr << std::endl;
	return 1;
    }

    std::unique_ptr<std::istream> in;
    if (args.value<std::string>('f').empty()) {
	in = std::unique_ptr<std::istream>(&std::cin);
    } else {
	const std::string &infile = args.value<std::string>('f');
	in = std::unique_ptr<std::istream>(new bxz::ifstream(infile));
    }

    if (args.value<bool>('d')) {
	UnpackBuffered(args.value<size_t>('n'), args.value<size_t>('r'), in.get(), &std::cout);
    } else {
	BufferedPack(args.value<size_t>('n'), args.value<size_t>("buffer-size"), in.get(), &std::cout);
    }
    if (args.value<std::string>('f').empty()) {
	in.release(); // Release ownership of std::cout so we don't try to free it
    }

    return 0;
}

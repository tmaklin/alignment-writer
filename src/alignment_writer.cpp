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
#include <cstddef>
#include <string>
#include <algorithm>
#include <iostream>
#include <exception>
#include <memory>
#include <unordered_map>

#include "cxxargs.hpp"
#include "bxzstr.hpp"
#include "kseq++/seqio.hpp"

#include "version.h"
#include "unpack.hpp"
#include "pack.hpp"

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_short_argument<std::string>('f', "Pseudoalignment file, packed or unpacked, read from cin if not supplied.", "");
  args.add_short_argument<std::string>('r', "Input reads in .fastq(.gz) format (required for packing)");
  args.add_short_argument<bool>('d', "Unpack pseudoalignment.", false);
  args.add_short_argument<size_t>('n', "Number of reference sequences in the pseudoalignment (required for packing).");
  args.add_long_argument<size_t>("buffer-size", "Buffer size for buffered packing (default: 100000", (size_t)100000);
  args.add_long_argument<std::string>("format", "Input file format (one of `themisto` (default), `fulgor`, `bifrost`)", "themisto");
  if (CmdOptionPresent(argv, argv+argc, "-d")) {
      args.set_not_required('n');
  }
  args.add_long_argument<bool>("help", "Print the help message.", false);
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
      std::cerr << "\n" + args.help() << '\n' << '\n';
  }
  args.parse(argc, argv);
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

    alignment_writer::Format format;
    if (args.value<std::string>("format") == "themisto") {
	format = alignment_writer::themisto;
    } else if (args.value<std::string>("format") == "fulgor") {
	format = alignment_writer::fulgor;
    } else if (args.value<std::string>("format") == "bifrost") {
	format = alignment_writer::bifrost;
    } else {
	throw std::runtime_error("Unrecognized input format.");
    }

    std::unique_ptr<std::istream> in;
    if (args.value<std::string>('f').empty()) {
	in = std::unique_ptr<std::istream>(&std::cin);
    } else {
	const std::string &infile = args.value<std::string>('f');
	in = std::unique_ptr<std::istream>(new bxz::ifstream(infile));
    }

    if (args.value<bool>('d')) {
	alignment_writer::Print(in.get(), &std::cout);
    } else {

	std::unordered_map<std::string, size_t> query_to_position;
	try {
	    klibpp::KSeq record;
	    klibpp::SeqStreamIn iss(args.value<std::string>('r').c_str());
	    size_t read_pos = 0;
	    while (iss >> record) {
		query_to_position.insert(std::make_pair(std::string(record.name), read_pos));
		++read_pos;
	    }
	} catch (const std::exception &e) {
	    std::cerr << "Reading from input `-r " << args.value<std::string>('r') << "` failed:" << e.what() << std::endl;
	}

	if (query_to_position.size() == 0) {
	    std::cerr << "Input: `-r " << args.value<std::string>('r') << "` has no reads!" << std::endl;
	    return 1;
	}

	try {
	    alignment_writer::BufferedPack(format, query_to_position, args.value<size_t>('n'), args.value<size_t>("buffer-size"), in.get(), &std::cout);
	} catch (const std::invalid_argument &e) {
	    std::cerr << "Reading the alignment failed: " << e.what() << " (is `--format " << args.value<std::string>("format") << "` correct?)" << std::endl;
	    return 1;
	} catch (const std::exception &e) {
	    std::cerr << "Reading the alignment failed: " << e.what() << '.' << std::endl;
	}
    }
    if (args.value<std::string>('f').empty()) {
	in.release(); // Release ownership of std::cout so we don't try to free it
    }

    return 0;
}

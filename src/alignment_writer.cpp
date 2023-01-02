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

#include "cxxargs.hpp"
#include "bxzstr.hpp"

#include "version.h"
#include "unpack.hpp"
#include "pack.hpp"
#include "merge.hpp"

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_short_argument<std::string>('f', "Pseudoalignment file, packed or unpacked, read from cin if not supplied.", "");
  args.add_short_argument<bool>('d', "Unpack pseudoalignment.", false);
  args.add_short_argument<size_t>('n', "Number of reference sequences in the pseudoalignment (required).");
  args.add_short_argument<size_t>('r', "Number of reads in the pseudoalignment (required).");
  args.add_long_argument<size_t>("buffer-size", "Buffer size for buffered packing (default: 100000", (size_t)100000);
  args.add_long_argument<std::string>("merge", "Intersect the pseudoalignment from this file with the file in -f.", "");
  if (CmdOptionPresent(argv, argv+argc, "--merge")) {
      args.set_not_required('r');
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

    std::unique_ptr<std::istream> in;
    if (args.value<std::string>('f').empty()) {
	in = std::unique_ptr<std::istream>(&std::cin);
    } else {
	const std::string &infile = args.value<std::string>('f');
	in = std::unique_ptr<std::istream>(new bxz::ifstream(infile));
    }

    if (args.value<bool>('d')) {
	alignment_writer::UnpackBuffered(in.get(), &std::cout);
    } else if (!args.value<std::string>("merge").empty()) {
      const std::string &in_2 = args.value<std::string>("merge");
      std::unique_ptr<std::istream> pair(new bxz::ifstream(in_2));
      std::vector<std::istream*> infiles({ in.get(), pair.get() });
      alignment_writer::Merge(telescope::get_mode("intersection"), args.value<size_t>('n'), infiles, &std::cout);
    } else {
      alignment_writer::BufferedPack(args.value<size_t>('n'), args.value<size_t>('r'), args.value<size_t>("buffer-size"), in.get(), &std::cout);
    }
    if (args.value<std::string>('f').empty()) {
	in.release(); // Release ownership of std::cout so we don't try to free it
    }

    return 0;
}

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
#include <filesystem>

#include "cxxopts.hpp"
#include "bxzstr.hpp"
#include "kseq++/seqio.hpp"

#include "version.h"
#include "unpack.hpp"
#include "pack.hpp"
#include "parser.hpp"

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

bool parse_args(int argc, char* argv[], cxxopts::Options &options) {
    options.add_options()
	("z,compress", "Compress file(s).", cxxopts::value<bool>()->default_value("false"))
	("d,decompress", "Decompress file(s).", cxxopts::value<bool>()->default_value("false"))
	("r,reads", "Reads used in the input alignment.", cxxopts::value<std::string>())
	("l,target-list", "File listing the input alignment targets.", cxxopts::value<std::string>())
	("format", "Input/output format", cxxopts::value<std::string>()->default_value("themisto"))
	("k,keep", "Keep input file(s) instead of deleting.", cxxopts::value<bool>()->default_value("false"))
	("f,force", "Force overwrite output file(s).", cxxopts::value<bool>()->default_value("false"))
	("c,stdout", "Write to standard out, keep files.", cxxopts::value<bool>()->default_value("false"))
	("T,threads", "Use `arg` threads, 0 = all available.", cxxopts::value<size_t>()->default_value("1"))
	("b,buffer-size", "Buffer writes every `arg` hits.", cxxopts::value<size_t>()->default_value("256000"))
	("h,help", "Print this message and quit.", cxxopts::value<bool>()->default_value("false"))
	("V,version", "Print the version and quit.", cxxopts::value<bool>()->default_value("false"))
	("filenames", "Input files as positional arguments", cxxopts::value<std::vector<std::string>>()->default_value(""));

    // options.allow_unrecognised_options(); // Handle the levels
    options.positional_help("[files]");
    options.custom_help("[options]");
    options.parse_positional({ "filenames" });

    if (CmdOptionPresent(argv, argv+argc, "--help") || CmdOptionPresent(argv, argv+argc, "-h")) {
	std::cerr << options.help() << std::endl;
	return true; // quit
    } else if (CmdOptionPresent(argv, argv+argc, "--version") || CmdOptionPresent(argv, argv+argc, "-V")) {
	std::cerr << "alignment-writer " << ALIGNMENT_WRITER_BUILD_VERSION << std::endl;
	return true; // quit
    }

    return false; // parsing successfull
}

bool file_exists(const std::string &file_path) {
    std::filesystem::path check_file{ file_path };
    return std::filesystem::exists(check_file);
}

void read_compression_inputs(const std::string &reads_file, const std::string &refs_file,
			     std::unordered_map<std::string, size_t> *query_to_position,
			     std::unordered_map<std::string, size_t> *ref_to_position) {
    klibpp::KSeq record;
    klibpp::SeqStreamIn iss(reads_file.c_str());
    size_t read_pos = 0;
    while (iss >> record) {
	query_to_position->insert(std::make_pair(std::string(record.name), read_pos));
	++read_pos;
    }

    if (query_to_position->size() == 0) {
	throw std::runtime_error("Input `--reads " + reads_file + "` has no reads!");
    }

    bxz::ifstream in(refs_file);
    std::string line;
    size_t ref_pos = 0;
    while (std::getline(in, line)) {
	ref_to_position->insert(std::make_pair(line, ref_pos));
	++ref_pos;
    }
    in.close();

    if (ref_to_position->size() == 0) {
	throw std::runtime_error("Input `--target-list " + refs_file + "` is empty!");
    }
}

int main(int argc, char* argv[]) {
    std::string program_name = "alignment-writer";
    std::string fileformat_suffix = ".aln";
    cxxopts::Options opts(program_name, program_name + ": compress or decompress pseudoalignment files.");
    try {
	bool quit = parse_args(argc, argv, opts);
	if (quit) {
	    return 0;
	}
    } catch (std::exception &e) {
	std::cerr << program_name + ": parsing arguments failed: " << std::string(e.what()) << std::endl;
	return 1;
    }

    cxxopts::ParseResult args;
    try {
	args = opts.parse(argc, argv);
    } catch (const cxxopts::exceptions::no_such_option &e) {
	std::cerr << program_name + ": " << std::string(e.what()) << std::endl;
	return 1;
    }

    alignment_writer::Format format;
    if (args["format"].as<std::string>() == "themisto") {
	format = alignment_writer::themisto;
    } else if (args["format"].as<std::string>() == "fulgor") {
	format = alignment_writer::fulgor;
    } else if (args["format"].as<std::string>() == "bifrost") {
	format = alignment_writer::bifrost;
    } else if (args["format"].as<std::string>() == "metagraph") {
	format = alignment_writer::metagraph;
    } else if (args["format"].as<std::string>() == "sam") {
	format = alignment_writer::sam;
    } else {
	throw std::runtime_error("Unrecognized input format.");
    }

    // 0 == read from cin
    const std::vector<std::string> &input_files = args["filenames"].as<std::vector<std::string>>();
    size_t n_input_files = input_files.size();

    if (n_input_files > 1 && !args["force"].as<bool>() && !args["decompress"].as<bool>()) {
	std::cerr << program_name + ": refusing to compress more than 1 input files. Use -f to force compression\nNote: multiple inputs must have the same --reads and --target-list." << std::endl;
	return 1;
    }

    // Don't write compressed data to terminal unless asked
    if (!args["force"].as<bool>() && !args["decompress"].as<bool>() && input_files[0].empty()) {
	// Refuse to write to terminal without -f or -c
	if (isatty(fileno(stdout))) {
	    std::cerr << program_name + ": refusing to write compressed data to terminal. Use -f to force write.\n" + program_name + ": try `" + program_name + "--help` for help." << std::endl;
	    return 1;
	}
    }

    // Decompress from and to terminal
    if (!isatty(fileno(stdin))) {
	// Compress from cin to cout
	if (args["decompress"].as<bool>()) {
	    try {
		alignment_writer::Print(format, &std::cin, &std::cout);
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": error in reading compressed data from terminal." << std::endl;
		return 1;
	    }
	} else {
	    std::string reads_file;
	    std::string refs_file;
	    try {
		reads_file = args["reads"].as<std::string>().c_str();
	    } catch (const cxxopts::exceptions::option_has_no_value &e) {
		std::cerr << program_name + ": ERROR: option --reads has no value!" << std::endl;
		return 1;
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": ERROR: " << e.what() << std::endl;
		return 1;
	    }
	    try {
		refs_file = args["target-list"].as<std::string>().c_str();
	    } catch (const cxxopts::exceptions::option_has_no_value &e) {
		std::cerr << program_name + ": ERROR: option --target-list has no value!" << std::endl;
		return 1;
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": ERROR: " << e.what() << std::endl;
		return 1;
	    }
	    std::unordered_map<std::string, size_t> query_to_position;
	    std::unordered_map<std::string, size_t> ref_to_position;
	    try {
		read_compression_inputs(reads_file, refs_file, &query_to_position, &ref_to_position);
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": ERROR: " << e.what() << std::endl;
		return 1;
	    }
	    try {
		alignment_writer::BufferedPack(format, query_to_position, ref_to_position, args["buffer-size"].as<size_t>(), &std::cin, &std::cout);
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": error in reading data from terminal." << std::endl;
		return 1;
	    }
	}
    }

    // Decompress input files
    if (!input_files[0].empty()) {
	// Compress/decompress all positional arguments
	if (!args["decompress"].as<bool>() || args["compress"].as<bool>()) {
	    // Run compression loop
	    std::string reads_file;
	    std::string refs_file;
	    try {
		reads_file = args["reads"].as<std::string>().c_str();
	    } catch (const cxxopts::exceptions::option_has_no_value &e) {
		std::cerr << program_name + ": ERROR: option --reads has no value!" << std::endl;
		return 1;
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": ERROR: \"" << e.what() << '"' << std::endl;
		return 1;
	    }
	    try {
		refs_file = args["target-list"].as<std::string>().c_str();
	    } catch (const cxxopts::exceptions::option_has_no_value &e) {
		std::cerr << program_name + ": ERROR: option --target-list has no value!" << std::endl;
		return 1;
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": ERROR: " << e.what() << std::endl;
		return 1;
	    }
	    std::unordered_map<std::string, size_t> query_to_position;
	    std::unordered_map<std::string, size_t> ref_to_position;
	    try {
		read_compression_inputs(reads_file, refs_file, &query_to_position, &ref_to_position);
	    } catch (const std::exception &e) {
		std::cerr << program_name + ": ERROR: " << e.what() << std::endl;
		return 1;
	    }

	    for (size_t i = 0; i < n_input_files; ++i) {
		const std::string &infile = input_files[i];
		if (!file_exists(infile)) {
		    std::cerr << program_name + ": " << infile << ": no such file or directory." << std::endl;
		    return 1;
		}
		bxz::ifstream in_stream(infile); // Handle compressed inputs
		if (!args["stdout"].as<bool>()) {
		    // Add .aln suffix to infile name
		    const std::string &outfile = infile + fileformat_suffix;
		    if (file_exists(outfile) && !args["force"].as<bool>()) {
			std::cerr << program_name + ": " << outfile << ": file exists; use `--force` to overwrite." << std::endl;
			return 1;
		    }
		    std::ofstream out_stream(outfile);

		    try {
			alignment_writer::BufferedPack(format, query_to_position, ref_to_position, args["buffer-size"].as<size_t>(), &in_stream, &out_stream);
		    } catch (const std::exception &e) {
			std::cerr << program_name + ": error in compressing file " << infile << " to file " << outfile << '.' << std::endl;
			return 1;
		    }
		} else {
		    // Compress to cout
		    try {
			alignment_writer::BufferedPack(format, query_to_position, ref_to_position, args["buffer-size"].as<size_t>(), &in_stream, &std::cout);
		    } catch (const std::exception &e) {
			std::cerr << program_name + ": error in compressing file " << infile << '.' << std::endl;
			return 1;
		    }
		}

		if (!args["keep"].as<bool>() && !args["stdout"].as<bool>()) {
		    std::filesystem::path remove_file{ infile };
		    std::filesystem::remove(remove_file);
		}
	    }
	} else if (args["decompress"].as<bool>() && !args["compress"].as<bool>()) {
	    // Run decompressor loop
	    for (size_t i = 0; i < n_input_files; ++i) {
		const std::string &infile = input_files[i];
		if (!file_exists(infile)) {
		    std::cerr << program_name + ": " << infile << ": no such file or directory." << std::endl;
		    return 1;
		}

		// Remove the .aln suffix from the infile name
		size_t lastindex = infile.find_last_of(".");

		// Decompresses to cout if `outfile` equals empty
		std::string outfile = (args["stdout"].as<bool>() ? "" : infile.substr(0, lastindex));

		if (file_exists(outfile) && !args["force"].as<bool>()) {
		    std::cerr << program_name + ": " << outfile << ": file exists; use `--force` to overwrite." << std::endl;
		    return 1;
		}

		bxz::ifstream in_stream(infile);
		if (args["stdout"].as<bool>()) {
		    try {
			alignment_writer::Print(format, &in_stream, &std::cout);
		    } catch (const std::exception &e) {
			std::cerr << program_name + ": error in reading compressed file " + infile << ':' << e.what() << std::endl;
			return 1;
		    }
		} else {
		    std::ofstream out_stream(outfile);
		    try {
			alignment_writer::Print(format, &in_stream, &out_stream);
		    } catch (const std::exception &e) {
			std::cerr << program_name + ": error in decompressing file " + infile << " to file " << outfile << '.' << std::endl;
			return 1;
		    }
		}

		if (!args["keep"].as<bool>() && !args["stdout"].as<bool>()) {
		    std::filesystem::path remove_file{ infile };
		    std::filesystem::remove(remove_file);
		}
	    }
	}
    }

    return 0;
}

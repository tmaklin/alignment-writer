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
#ifndef ALIGNMENT_WRITER_PARSER_HPP
#define ALIGNMENT_WRITER_PARSER_HPP

#include <cstddef>
#include <string>
#include <unordered_map>
#include <sstream>
#include <vector>

#include "bm64.h"

namespace alignment_writer {
enum Format { themisto, fulgor, bifrost, metagraph, sam };

inline size_t ThemistoParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::vector<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Themisto* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = ' ';
    std::stringstream stream(line);
    std::string part;
    std::getline(stream, part, separator);
    size_t read_id = std::stoul(part); // First column is a numerical ID for the read
    reads_in_buffer->push_back(read_id);
    size_t n_alignments = 0;
    while(std::getline(stream, part, separator)) {
	// Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	(*it) = read_id*n_refs + std::stoul(part);
	++n_alignments;
    }
    return n_alignments;
}

inline size_t FulgorParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::vector<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Fulgor* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = '\t';
    std::stringstream stream(line);
    std::string part;
    std::string query_name;
    std::getline(stream, query_name, separator); // First column is the fragment name
    std::getline(stream, part, separator); // Second column is the number of alignments
    size_t n_alignments = std::stoul(part);
    size_t read_id = query_to_position.at(query_name);
    reads_in_buffer->push_back(read_id);
    while(std::getline(stream, part, separator)) {
	// Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	(*it) = read_id*n_refs + std::stoul(part);
    }
    return n_alignments;
}

inline size_t BifrostParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::vector<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Bifrost* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = '\t';
    std::stringstream stream(line);
    std::string part;
    std::string query_name;
    std::getline(stream, query_name, separator); // First column is the fragment name
    size_t n_alignments = 0;
    size_t ref_id = 0;
    size_t read_id = query_to_position.at(query_name);
    reads_in_buffer->push_back(read_id);
    while(std::getline(stream, part, separator)) {
	// Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	bool aligned = std::stoul(part) == 1;
	if (aligned) {
	    (*it) = read_id*n_refs + ref_id;
	    ++n_alignments;
	}
	++ref_id;
    }
    return n_alignments;
}

inline size_t MetagraphParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::vector<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Bifrost* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = '\t';
    std::stringstream stream(line);
    std::string part;
    std::getline(stream, part, separator); // first column is the read position in input
    std::string query_name;
    std::getline(stream, query_name, separator); // second column is the fragment name
    size_t read_id = query_to_position.at(query_name);
    reads_in_buffer->push_back(read_id);
    size_t n_alignments = 0;
    std::getline(stream, part, separator); // Third column contains `:` separated alignments
    separator = ':';
    std::stringstream alns(part);
    std::string ref_name;
    while(std::getline(alns, ref_name, separator)) {
	// Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	size_t ref_id = ref_to_position.at(ref_name);
	(*it) = read_id*n_refs + ref_id;
	++n_alignments;
    }
    return n_alignments;
}

inline size_t SAMParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::vector<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Bifrost* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = '\t';
    std::stringstream stream(line);
    std::string part;

    // First column is query name
    std::string query_name;
    std::getline(stream, query_name, separator);
    size_t read_id = query_to_position.at(query_name);
    reads_in_buffer->push_back(read_id);

    // Second column contains the sam bitwise flags
    std::getline(stream, part, separator);

    // Third column contains the reference name
    std::string ref_name;
    std::getline(stream, ref_name, separator);
    if (ref_name == "*") {
	// Unmapped
	return 0;
    }
    size_t ref_id = ref_to_position.at(ref_name);

    (*it) = read_id*n_refs + ref_id;

    return 1;
}
}

#endif

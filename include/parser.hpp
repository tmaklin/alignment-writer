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
#include <unordered_set>

#include "bm64.h"
#include "nlohmann/json.hpp"

#include "version.h"

namespace alignment_writer {
enum Format { themisto, fulgor, bifrost, metagraph, sam };
inline std::string format_to_string(const Format &format) {
    switch (format) {
	    case themisto:   return "themisto";
	    case fulgor:     return "fulgor";
	    case bifrost:    return "bifrost";
	    case metagraph:  return "metagraph";
	    case sam:        return "SAM";
	    default:         return "themisto";
    }
}

inline size_t ThemistoParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::unordered_set<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Themisto* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = ' ';
    std::stringstream stream(line);
    std::string part;
    std::getline(stream, part, separator);
    size_t read_id = std::stoul(part); // First column is a numerical ID for the read
    reads_in_buffer->insert(read_id);
    size_t n_alignments = 0;
    while(std::getline(stream, part, separator)) {
	// Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	(*it) = read_id*n_refs + std::stoul(part);
	++n_alignments;
    }
    return n_alignments;
}

inline void ThemistoPrinter(const bm::bvector<> &bits, const nlohmann::json_abi_v3_11_3::json &header, const nlohmann::json_abi_v3_11_3::json &block_headers, std::ostream *out) {
    size_t n_reads = header["n_queries"];
    size_t n_refs = header["n_targets"];

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

inline size_t FulgorParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::unordered_set<size_t> *reads_in_buffer) {
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
    reads_in_buffer->insert(read_id);
    while(std::getline(stream, part, separator)) {
	// Buffered insertion to contiguously stored n_reads x n_refs pseudoalignment matrix
	(*it) = read_id*n_refs + std::stoul(part);
    }
    return n_alignments;
}

inline void FulgorPrinter(const bm::bvector<> &bits, const nlohmann::json_abi_v3_11_3::json &header, const nlohmann::json_abi_v3_11_3::json &block_headers, std::ostream *out) {
    size_t n_reads = header["n_queries"];
    size_t n_refs = header["n_targets"];

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    nlohmann::json_abi_v3_11_3::json query_info = block_headers["queries"];
    std::unordered_map<size_t, std::string> query_map;
    for (auto kv : query_info) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }
    for (size_t i = 0; i < n_reads; ++i) {
	*out << query_map.at(i) << '\t';
	std::vector<size_t> mapped_targets;
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    while (*en < i*n_refs + n_refs && en < en_end) {
		mapped_targets.emplace_back((*en) - i*n_refs);
		++en;
	    }
	}
	size_t n_mapped = mapped_targets.size();
	*out << n_mapped;
	if (n_mapped > 0) {
	    *out << '\t';
	    for (size_t i = 0; i < n_mapped; ++i) {
		*out << mapped_targets[i] << ((i == n_mapped - 1) ? '\n' : '\t');
	    }
	} else {
	    *out << '\n';
	}
    }
    out->flush(); // Flush
}

inline size_t BifrostParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::unordered_set<size_t> *reads_in_buffer) {
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
    reads_in_buffer->insert(read_id);
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

inline void BifrostPrinter(const bm::bvector<> &bits, const nlohmann::json_abi_v3_11_3::json &header, const nlohmann::json_abi_v3_11_3::json &block_headers, std::ostream *out) {
    size_t n_reads = header["n_queries"];
    size_t n_refs = header["n_targets"];

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    nlohmann::json_abi_v3_11_3::json query_info = block_headers["queries"];
    std::unordered_map<size_t, std::string> query_map;
    for (auto kv : query_info) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }

    std::vector<std::string> targets(n_refs);
    for (size_t i = 0; i < n_refs; ++i) {
	size_t pos = header["targets"].at(i)["pos"];
	targets[pos] = header["targets"].at(i)["target"];
    }

    *out << "query_name" << '\t';
    for (size_t i = 0; i < n_refs; ++i) {
	*out << targets[i] << ((i == n_refs - 1) ? '\n' : '\t');
    }

    for (size_t i = 0; i < n_reads; ++i) {
	*out << query_map.at(i) << '\t';
	std::vector<bool> alignment(n_refs, false);
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    while (*en < i*n_refs + n_refs && en < en_end) {
		alignment[(*en) - i*n_refs] = true;
		++en;
	    }
	}
	for (size_t j = 0; j < n_refs; ++j) {
	    *out << (alignment[j] ? "1" : "0") << ((j == n_refs - 1) ? '\n' : '\t');
	}
    }
    out->flush(); // Flush
}

inline size_t MetagraphParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::unordered_set<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Bifrost* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = '\t';
    std::stringstream stream(line);
    std::string part;
    std::getline(stream, part, separator); // first column is the read position in input
    std::string query_name;
    std::getline(stream, query_name, separator); // second column is the fragment name
    size_t read_id = query_to_position.at(query_name);
    reads_in_buffer->insert(read_id);
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

inline void MetagraphPrinter(const bm::bvector<> &bits, const nlohmann::json_abi_v3_11_3::json &header, const nlohmann::json_abi_v3_11_3::json &block_headers, std::ostream *out) {
    size_t n_reads = header["n_queries"];
    size_t n_refs = header["n_targets"];

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    nlohmann::json_abi_v3_11_3::json query_info = block_headers["queries"];
    std::unordered_map<size_t, std::string> query_map;
    for (auto kv : query_info) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }

    std::unordered_map<size_t, std::string> target_map;
    for (size_t i = 0; i < n_refs; ++i) {
	size_t pos = header["targets"].at(i)["pos"];
	std::string target = header["targets"].at(i)["target"];
	target_map.insert(std::make_pair(pos, target));
    }

    for (size_t i = 0; i < n_reads; ++i) {
	*out << i << '\t' << query_map.at(i) << '\t';
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    bool first = true;
	    while (*en < i*n_refs + n_refs && en < en_end) {
		if (!first) {
		    *out << ':';
		}
		first = false;
		*out << target_map.at((*en) - i*n_refs);
		++en;
	    }
	}
	*out << '\n';
    }
    out->flush(); // Flush
}

inline size_t SAMParser(const std::string &line, const std::unordered_map<std::string, size_t> &query_to_position, const std::unordered_map<std::string, size_t> &ref_to_position, bm::bvector<>::bulk_insert_iterator *it, std::unordered_set<size_t> *reads_in_buffer) {
    // Reads a pseudoalignment line stored in the *Bifrost* format and returns the number of pseudoalignments on the line
    size_t n_refs = ref_to_position.size();
    char separator = '\t';
    std::stringstream stream(line);
    std::string part;

    // First column is query name
    std::string query_name;
    std::getline(stream, query_name, separator);
    size_t read_id = query_to_position.at(query_name);
    reads_in_buffer->insert(read_id);

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

inline void SAMPrinter(const bm::bvector<> &bits, const nlohmann::json_abi_v3_11_3::json &header, const nlohmann::json_abi_v3_11_3::json &block_headers, std::ostream *out) {
    size_t n_reads = header["n_queries"];
    size_t n_refs = header["n_targets"];

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    nlohmann::json_abi_v3_11_3::json query_info = block_headers["queries"];
    std::unordered_map<size_t, std::string> query_map;
    for (auto kv : query_info) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }

    std::unordered_map<size_t, std::string> target_map;
    for (size_t i = 0; i < n_refs; ++i) {
	size_t pos = header["targets"].at(i)["pos"];
	std::string target = header["targets"].at(i)["target"];
	target_map.insert(std::make_pair(pos, target));
    }

    for (size_t i = 0; i < n_refs; ++i) {
	*out << "@SQ" << '\t' << "SN:" << target_map.at(i) << '\n';
    }
    *out << "@PG" << '\t' << "ID:" << std::string(header["input_format"]) << '\t' << "PN:alignment-writer" << '\t' << "VN:" << ALIGNMENT_WRITER_BUILD_VERSION << '\n';

    for (size_t i = 0; i < n_reads; ++i) {
	std::vector<size_t> targets;
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    while (*en < i*n_refs + n_refs && en < en_end) {
		targets.emplace_back((*en) - i*n_refs);
		++en;
	    }
	}
	if (targets.size() > 0) {
	    for (size_t j = 0; j < targets.size(); ++j) {
		*out << query_map.at(i) << '\t' << "0" << '\t' << target_map.at(targets[j]) << '\t';
		*out << "1" << '\t' << "255" << '\t' << '*' << '\t' << '*' << '\t' << '0' << '\t' << '0' << '\t';
		*out << '*' << '\t' << '*' << '\n';
	    }
	} else {
	    *out << query_map.at(i) << '\t' << "0" << '\t' << '*' << '\t';
	    *out << "0" << '\t' << "255" << '\t' << '*' << '\t' << '*' << '\t' << '0' << '\t' << '0' << '\t';
	    *out << '*' << '\t' << '*' << '\n';
	}
    }
    out->flush(); // Flush
}
}

#endif

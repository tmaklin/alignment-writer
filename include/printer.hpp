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
#ifndef ALIGNMENT_WRITER_PRINTER_HPP
#define ALIGNMENT_WRITER_PRINTER_HPP

#include <cstddef>
#include <string>
#include <unordered_map>
#include <map>
#include <sstream>
#include <unordered_set>
#include <functional>

#include "bm64.h"
#include "nlohmann/json.hpp"

#include "Alignment.hpp"
#include "version.h"

namespace alignment_writer {
inline std::stringbuf ThemistoPrinter(const Alignment &bits) {
    size_t n_reads = bits.queries();
    size_t n_refs = bits.targets();

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    std::string out;
    for (size_t i = 0; i < n_reads; ++i) {
	// Write read id (data compressed with Pack() is sorted so read id is just the iterator id)
	out += std::to_string(i);
	out += ' ';
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    while (*en < i*n_refs + n_refs && en < en_end) {
		out += std::to_string((*en) - i*n_refs);
		out += ' ';
		++en;
	    }
	}
	out += '\n';
    }
    return std::stringbuf(out);
}

inline std::stringbuf FulgorPrinter(const Alignment &bits) {
    size_t n_reads = bits.queries();
    size_t n_refs = bits.targets();

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    std::unordered_map<size_t, std::string> query_map;
    for (auto kv : bits.annotation()) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }

    std::string out;
    for (size_t i = 0; i < n_reads; ++i) {
	out += query_map.at(i);
	out += '\t';
	std::vector<size_t> mapped_targets;
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    while (*en < i*n_refs + n_refs && en < en_end) {
		mapped_targets.emplace_back((*en) - i*n_refs);
		++en;
	    }
	}
	size_t n_mapped = mapped_targets.size();
	out += std::to_string(n_mapped);
	if (n_mapped > 0) {
	    out += '\t';
	    for (size_t i = 0; i < n_mapped; ++i) {
		out += std::to_string(mapped_targets[i]);
		out += ((i == n_mapped - 1) ? '\n' : '\t');
	    }
	} else {
	    out += '\n';
	}
    }
    return std::stringbuf(out);
}

inline std::stringbuf BifrostPrinter(const Alignment &bits) {
    size_t n_reads = bits.queries();
    size_t n_refs = bits.targets();

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    std::unordered_map<size_t, std::string> query_map;
    for (auto kv : bits.annotation()) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }

    const std::vector<std::string> &targets = bits.target_names();

    std::string out;
    out += "query_name";
    out += '\t';
    for (size_t i = 0; i < n_refs; ++i) {
	out += targets[i];
	out += ((i == n_refs - 1) ? '\n' : '\t');
    }

    for (size_t i = 0; i < n_reads; ++i) {
	out += query_map.at(i);
	out += '\t';
	std::vector<bool> alignment(n_refs, false);
	if (*en < i*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    while (*en < i*n_refs + n_refs && en < en_end) {
		alignment[(*en) - i*n_refs] = true;
		++en;
	    }
	}
	for (size_t j = 0; j < n_refs; ++j) {
	    out += (alignment[j] ? "1" : "0");
	    out += ((j == n_refs - 1) ? '\n' : '\t');
	}
    }
    return std::stringbuf(out);
}

inline std::stringbuf MetagraphPrinter(const Alignment &bits) {
    size_t n_reads = bits.queries();
    size_t n_refs = bits.targets();

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    std::map<size_t, std::string> query_map; // Need to iterate in order
    for (auto kv : bits.annotation()) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }

    std::unordered_map<size_t, std::string> target_map;
    for (size_t i = 0; i < n_refs; ++i) {
	target_map.insert(std::make_pair(i, bits.target_names()[i]));
    }

    std::string out;
    for (auto query : query_map) {
	out += std::to_string(query.first);
	out += '\t';
	out += query.second;
	out += '\t';
	if (en != en_end && *en < query.first*n_refs + n_refs) { // Next pseudoalignment is for this read
	    // Write found pseudoalignments using the enumerator
	    bool first = true;
	    while (*en < query.first*n_refs + n_refs && en < en_end) {
		if (!first) {
		    out += ':';
		}
		first = false;
		out += target_map.at((*en) - query.first*n_refs);
		++en;
	    }
	}
	out += '\n';
    }
    return std::stringbuf(out);
}

inline std::stringbuf SAMPrinter(const Alignment &bits) {
    size_t n_reads = bits.queries();
    size_t n_refs = bits.targets();

    // Use an enumerator to traverse the pseudoaligned bits
    bm::bvector<>::enumerator en = bits.first();
    bm::bvector<>::enumerator en_end = bits.end();

    std::unordered_map<size_t, std::string> query_map;
    for (auto kv : bits.annotation()) {
	query_map.insert(std::make_pair(kv["pos"], kv["query"]));
    }

    std::unordered_map<size_t, std::string> target_map;
    for (size_t i = 0; i < n_refs; ++i) {
	target_map.insert(std::make_pair(i, bits.target_names()[i]));
    }

    std::string out;
    for (size_t i = 0; i < n_refs; ++i) {
	out += "@SQ";
	out += '\t';
	out += "SN:";
	out += target_map.at(i);
	out += '\n';
    }
    out += "@PG\tID:";
    out += std::string(bits.input_format());
    out += "\tPN:alignment-writer\tVN:";
    out += ALIGNMENT_WRITER_BUILD_VERSION;
    out += '\n';

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
		out += query_map.at(i);
		out +="\t0\t";
		out += target_map.at(targets[j]);
		out += "\t'1\t255\t*\t*\t0\t0\t";
		out += "*\t*\n";
	    }
	} else {
	    out += query_map.at(i);
	    out += "\t0\t*\t";
	    out += "0\t255\t*\t*\t0\t0\t";
	    out += "*\t*\n";
	}
    }
    return std::stringbuf(out);
}

class Printer {
public:
    Printer(const Format &_format) {
	if (_format == themisto) {
	    this->format = ThemistoPrinter;
	} else if (_format == fulgor) {
	    this->format = FulgorPrinter;
	} else if (_format == bifrost) {
	    this->format = BifrostPrinter;
	} else if (_format == metagraph) {
	    this->format = MetagraphPrinter;
	} else if (_format == sam) {
	    this->format = SAMPrinter;
	}
    }

    std::function<std::stringbuf(Alignment &bits)> format;
};
}

#endif

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
#ifndef ALIGNMENT_WRITER_ALIGNMENT_HPP
#define ALIGNMENT_WRITER_ALIGNMENT_HPP

#include <cstddef>
#include <vector>
#include <string>

#include "bm64.h"
#include "nlohmann/json.hpp"

#include "parser.hpp"

namespace alignment_writer {
using json = nlohmann::json_abi_v3_11_3::json;

class Alignment : public bm::bvector<> {
private:
    size_t n_queries;
    size_t n_targets;
    std::vector<std::string> targets_names;
    bm::bvector<> bitvector;
    std::string format;

public:
    Alignment(size_t _n_queries, const std::vector<std::string> &_targets_names) {
	this->n_queries = _n_queries;
	this->n_targets = _targets_names.size();
	this->targets_names = _targets_names;
	this->bitvector = bm::bvector<>((size_t)(n_queries * n_targets));
    }

    Alignment(size_t _n_queries, const std::unordered_map<std::string, size_t>& _targets_names_to_pos) {
	this->n_queries = _n_queries;
	this->n_targets = _targets_names_to_pos.size();
	this->targets_names = std::vector<std::string>(_targets_names_to_pos.size());
	for (auto kv : _targets_names_to_pos) {
	    this->targets_names[kv.second] = kv.first;
	}
	this->bitvector = bm::bvector<>((size_t)(n_queries * n_targets));
    }

    Alignment(const json& file_header) {
	this->n_queries = file_header["n_queries"];
	this->n_targets = file_header["n_targets"];
	this->targets_names = std::vector<std::string>(this->n_targets);
	for (auto kv : file_header["targets"]) {
	    this->targets_names[kv["pos"]] = kv["target"];
	}
	this->bitvector = bm::bvector<>((size_t)(this->n_queries * this->n_targets));
	this->format = file_header["input_format"];
    }

    const size_t queries() const { return this->n_queries; }
    const size_t targets() const { return this->n_targets; }

    const std::vector<std::string> target_names() const { return this->targets_names; }

    const std::string input_format() const { return this->format; }

    const bm::bvector<> raw() const { return this->bitvector; }
};
}

#endif

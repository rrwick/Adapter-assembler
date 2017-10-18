// Copyright 2017 Ryan Wick

// This file is part of Adapter-assembler

// Adapter-assembler is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
// version.

// Adapter-assembler is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.

// You should have received a copy of the GNU General Public License along with Adapter-assembler.  If not, see
// <http://www.gnu.org/licenses/>.


#include "kmers.h"

#include <iostream>
#include <algorithm>
#include <zlib.h>
#include "kseq.h"
#include "misc.h"

KSEQ_INIT(gzFile, gzread)


Kmers::Kmers(int kmer_size) {
    m_kmer_size = size_t(kmer_size);
}


void Kmers::add_fastq(std::string filename, bool start, int margin) {

    std::cerr << "Hashing " << m_kmer_size << "-mers from " << filename;
    if (start)
        std::cerr << " starts\n";
    else  // end
        std::cerr << " ends\n";

    int l;
    int sequence_count = 0;

    long long base_count = 0;
    long long last_progress = 0;

    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t * seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        if (l == -3)
            std::cerr << "Error reading " << filename << "\n";
        else {
            ++sequence_count;

            base_count += seq->seq.l;
            char * sequence = seq->seq.s;

            int range_start, range_end;
            if (start) {
                range_start = 0;
                range_end = std::min(margin, int(seq->seq.l)) + 1 - int(m_kmer_size);
            }
            else {  // end
                range_start = std::max(int(seq->seq.l) - margin, 0);
                range_end = int(seq->seq.l) + 1 - int(m_kmer_size);
            }

            for (int i = range_start; i < range_end; ++i)
                add_kmer(kmer_to_bits(sequence + i));

            if (base_count - last_progress >= 483611) {  // a big prime number so progress updates don't round off
                last_progress = base_count;
                print_hash_progress(filename, base_count);
            }
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    print_hash_progress(filename, base_count);

    std::cerr << "\n  " << int_to_string(sequence_count) << " reads, "
              << int_to_string(int(m_kmers.size())) << " " << m_kmer_size << "-mers\n\n";
}


void Kmers::add_kmer(uint32_t kmer) {
    if (!is_kmer_present(kmer))
        m_kmers[kmer] = 1;
    else
        ++m_kmers[kmer];
}


bool Kmers::is_kmer_present(uint32_t kmer) {
    return m_kmers.find(kmer) != m_kmers.end();
}


uint32_t Kmers::base_to_bits(char base) {
    switch (base) {
        case 'A':
            return 0;  // 00000000000000000000000000000000
        case 'C':
            return 1;  // 00000000000000000000000000000001
        case 'G':
            return 2;  // 00000000000000000000000000000010
        case 'T':
            return 3;  // 00000000000000000000000000000011
        case 'a':
            return 0;
        case 'c':
            return 1;
        case 'g':
            return 2;
        case 't':
            return 3;
        default:
            return 0;
    }
}

char Kmers::bits_to_base(uint32_t bits) {
    switch (bits) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            return 'A';
    }
}


uint32_t Kmers::kmer_to_bits(char * sequence) {
    uint32_t kmer = 0;
    for (size_t i = 0; i < m_kmer_size; ++i) {
        kmer <<= 2;
        kmer |= base_to_bits(sequence[i]);
    }
    return kmer;
}


uint32_t Kmers::kmer_to_bits(std::string sequence) {
    uint32_t kmer = 0;
    for (size_t i = 0; i < m_kmer_size; ++i) {
        kmer <<= 2;
        kmer |= base_to_bits(sequence[i]);
    }
    return kmer;
}


std::string Kmers::bits_to_kmer(uint32_t bits) {
    std::string kmer;
    for (size_t i = 0; i < m_kmer_size; ++i) {
        kmer.insert(0, 1, bits_to_base(bits % 4));
        bits >>= 2;
    }
    return kmer;
}


void Kmers::remove_low_depth_kmers(int min_depth) {
    std::vector<uint32_t> kmers_to_remove;
    for (auto kv : m_kmers) {
        uint32_t kmer = kv.first;
        int count = kv.second;
        if (count < min_depth)
            kmers_to_remove.push_back(kmer);
    }
    for (auto kmer : kmers_to_remove)
        m_kmers.erase(kmer);
}


void Kmers::output_gfa() {
    std::vector<uint32_t> all_kmers;
    for (auto kv : m_kmers)
        all_kmers.push_back(kv.first);
    std::sort(all_kmers.begin(), all_kmers.end());

    for (auto kmer : all_kmers)
        print_segment_line(kmer);
    for (auto kmer : all_kmers) {
        for (auto next : get_downstream_kmers(kmer))
            print_link_line(kmer, next);
    }
}


void Kmers::print_segment_line(uint32_t kmer) {
    std::cout << "S\t" << kmer << "\t" << bits_to_kmer(kmer) << "\tdp:f:" << m_kmers[kmer] << "\n";
}


void Kmers::print_link_line(uint32_t kmer_1, uint32_t kmer_2) {
    std::cout << "L\t";
    std::cout << kmer_1 << "\t+\t";
    std::cout << kmer_2 << "\t+\t";
    uint32_t overlap = m_kmer_size - 1;
    std::cout << overlap << "M\t\n";
}


std::vector<uint32_t> Kmers::get_upstream_kmers(uint32_t kmer) {
    std::string kmer_string = bits_to_kmer(kmer);
    kmer_string.pop_back();

    std::string next_1 = 'A' + kmer_string;
    std::string next_2 = 'C' + kmer_string;
    std::string next_3 = 'G' + kmer_string;
    std::string next_4 = 'T' + kmer_string;

    uint32_t next_1_bits = kmer_to_bits(next_1);
    uint32_t next_2_bits = kmer_to_bits(next_2);
    uint32_t next_3_bits = kmer_to_bits(next_3);
    uint32_t next_4_bits = kmer_to_bits(next_4);

    std::vector<uint32_t> upstream_kmers;

    if (is_kmer_present(next_1_bits))
        upstream_kmers.push_back(next_1_bits);
    if (is_kmer_present(next_2_bits))
        upstream_kmers.push_back(next_2_bits);
    if (is_kmer_present(next_3_bits))
        upstream_kmers.push_back(next_3_bits);
    if (is_kmer_present(next_4_bits))
        upstream_kmers.push_back(next_4_bits);

    return upstream_kmers;
}


std::vector<uint32_t> Kmers::get_downstream_kmers(uint32_t kmer) {
    std::string kmer_string = bits_to_kmer(kmer);
    kmer_string.erase(0, 1);

    std::string next_1 = kmer_string + 'A';
    std::string next_2 = kmer_string + 'C';
    std::string next_3 = kmer_string + 'G';
    std::string next_4 = kmer_string + 'T';

    uint32_t next_1_bits = kmer_to_bits(next_1);
    uint32_t next_2_bits = kmer_to_bits(next_2);
    uint32_t next_3_bits = kmer_to_bits(next_3);
    uint32_t next_4_bits = kmer_to_bits(next_4);

    std::vector<uint32_t> downstream_kmers;

    if (is_kmer_present(next_1_bits))
        downstream_kmers.push_back(next_1_bits);
    if (is_kmer_present(next_2_bits))
        downstream_kmers.push_back(next_2_bits);
    if (is_kmer_present(next_3_bits))
        downstream_kmers.push_back(next_3_bits);
    if (is_kmer_present(next_4_bits))
        downstream_kmers.push_back(next_4_bits);

    return downstream_kmers;
}


int Kmers::get_max_depth() {
    int max_depth = 0;
    for (auto kv : m_kmers)
        max_depth = std::max(max_depth, kv.second);
    return max_depth;
}


void Kmers::remove_tips() {
    std::vector<uint32_t> kmers_to_remove;
    for (auto kv : m_kmers) {
        uint32_t kmer = kv.first;
        int count = kv.second;

        std::vector<uint32_t> upstream = get_upstream_kmers(kmer);
        std::vector<uint32_t> downstream = get_downstream_kmers(kmer);

        if (downstream.empty()) {
            int max_upstream_count = 0;
            for (auto upstream_kmer : upstream)
                max_upstream_count = std::max(max_upstream_count, m_kmers[upstream_kmer]);
            if (max_upstream_count > count * 2)
                kmers_to_remove.push_back(kmer);
        }

        else if (upstream.empty()) {
            int max_downstream_count = 0;
            for (auto downstream_kmer : downstream)
                max_downstream_count = std::max(max_downstream_count, m_kmers[downstream_kmer]);
            if (max_downstream_count > count * 2)
                kmers_to_remove.push_back(kmer);
        }
    }

    for (auto kmer : kmers_to_remove)
        m_kmers.erase(kmer);
}


void Kmers::remove_large_diff() {
    std::vector<uint32_t> kmers_to_remove;
    for (auto kv : m_kmers) {
        uint32_t kmer = kv.first;
        int count = kv.second;

        std::vector<uint32_t> neighbours = get_upstream_kmers(kmer);
        std::vector<uint32_t> downstream = get_downstream_kmers(kmer);
        neighbours.insert(neighbours.end(), downstream.begin(), downstream.end());

        int max_neighbour_count = 0;
        for (auto neighbour : neighbours)
            max_neighbour_count = std::max(max_neighbour_count, m_kmers[neighbour]);
        if (max_neighbour_count > count * 5)
            kmers_to_remove.push_back(kmer);
    }

    for (auto kmer : kmers_to_remove)
        m_kmers.erase(kmer);
}


void Kmers::remove_singletons() {
    std::vector<uint32_t> kmers_to_remove;
    for (auto kv : m_kmers) {
        uint32_t kmer = kv.first;

        std::vector<uint32_t> neighbours = get_upstream_kmers(kmer);
        std::vector<uint32_t> downstream = get_downstream_kmers(kmer);
        neighbours.insert(neighbours.end(), downstream.begin(), downstream.end());

        int num_neighbours = 0;
        for (auto neighbour : neighbours) {
            if (neighbour != kmer)
                num_neighbours += 1;
        }
        if (num_neighbours == 0)
            kmers_to_remove.push_back(kmer);
    }

    for (auto kmer : kmers_to_remove)
        m_kmers.erase(kmer);
}

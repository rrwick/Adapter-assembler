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

#ifndef KMERS_H
#define KMERS_H


#include <string>
#include <vector>
#include <unordered_map>


class Kmers
{
public:
    Kmers(int kmer_size);

    void add_fastq(std::string filename, bool start, int margin);
    void remove_low_depth_kmers(int min_depth);
    void output_gfa();
    bool is_kmer_present(uint32_t kmer);

    uint32_t kmer_to_bits(char * sequence);
    uint32_t kmer_to_bits(std::string sequence);
    uint32_t base_to_bits(char base);

    std::string bits_to_kmer(uint32_t kmer);
    char bits_to_base(uint32_t kmer);

private:
    size_t m_kmer_size;
    std::unordered_map<uint32_t, int> m_kmers;

    std::vector<uint32_t> get_upstream_kmers(uint32_t kmer);
    std::vector<uint32_t> get_downstream_kmers(uint32_t kmer);

    void print_segment_line(uint32_t kmer);
    void print_link_line(uint32_t kmer_1, uint32_t kmer_2);

    void add_kmer(uint32_t kmer);
};


#endif // KMERS_H

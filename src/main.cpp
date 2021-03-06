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


#include <iostream>
#include <vector>
#include <unordered_map>

#include "arguments.h"
#include "kmers.h"
#include "misc.h"

#define PROGRAM_VERSION "0.1.0"


int main(int argc, char **argv)
{
    Arguments args(argc, argv);
    if (args.parsing_result == BAD)
        return 1;
    else if (args.parsing_result == HELP)
        return 0;
    else if (args.parsing_result == VERSION) {
        std::cout << "Adapter-assembler v" << PROGRAM_VERSION << "\n";
        return 0;
    }

    std::cerr << "\n";

    Kmers kmers(args.kmer);
    for (auto read_file : args.input_reads)
        kmers.add_fastq(read_file, args.start, args.margin);

    int max_depth = kmers.get_max_depth();
    std::cerr << "Maximum depth: " << max_depth << "\n";
    auto filter_depth = int(max_depth * args.filter_depth);
    std::cerr << "Filter depth:  " << filter_depth << "\n\n";


    std::cerr << "Cleaning step                      Remaining " << args.kmer << "-mers\n";
    std::cerr << "-------------------------------------------------------\n";

    std::cerr << "remove low-depth nodes             ";
    kmers.remove_low_depth_kmers(filter_depth);
    std::cerr << int_to_string(kmers.get_kmer_count()) << "\n";

    std::cerr << "prune tips                         ";
    kmers.remove_tips();
    std::cerr << int_to_string(kmers.get_kmer_count()) << "\n";

    std::cerr << "remove large differences           ";
    kmers.remove_large_diff();
    std::cerr << int_to_string(kmers.get_kmer_count()) << "\n";

    std::cerr << "remove singletons                  ";
    kmers.remove_singletons();
    std::cerr << int_to_string(kmers.get_kmer_count()) << "\n";

    kmers.output_gfa();

    std::cerr << "\n";
    return 0;
}

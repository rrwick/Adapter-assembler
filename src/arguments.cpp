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


#include "arguments.h"

#include <iostream>
#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>

#include "args.h"


struct DoublesReader
{
    void operator()(const std::string &name, const std::string &value, double &destination) {
        try {
            if (value.find_first_not_of("0123456789.") != std::string::npos)
                throw std::invalid_argument("");
            destination = std::stod(value);
        }
        catch ( ... ) {
            std::ostringstream problem;
            problem << "Error: argument '" << name << "' received invalid value type '" << value << "'";
            throw args::ParseError(problem.str());
        }
    }
};


typedef args::ValueFlag<double, DoublesReader> d_arg;
typedef args::ValueFlag<long long> i_arg;
typedef args::ValueFlag<std::string> s_arg;
typedef args::Flag f_arg;


Arguments::Arguments(int argc, char **argv) {

    args::ArgumentParser parser("Adapter-assembler: a tool for extracting adapter sequences from long reads",
                                "For more information, go to: https://github.com/rrwick/Adapter-assembler");
    parser.LongSeparator(" ");

    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    int terminal_width = w.ws_col;

    int indent_size;
    if (terminal_width > 120)
        indent_size = 4;
    else if (terminal_width > 80)
        indent_size = 3;
    else if (terminal_width > 60)
        indent_size = 2;
    else
        indent_size = 1;

    parser.helpParams.showTerminator = false;
    parser.helpParams.progindent = 0;
    parser.helpParams.descriptionindent = 0;
    parser.helpParams.width = terminal_width;
    parser.helpParams.flagindent = indent_size;
    parser.helpParams.eachgroupindent = indent_size;


    i_arg kmer_arg(parser, "int",
                   "k-mer size used for assembly (default: 10)",
                   {'k', "kmer"}, 10);
    i_arg min_depth_arg(parser, "int",
                   "k-mers with depth lower than this will be filtered out (default: 10)",
                   {'m', "min_depth"}, 10);

    args::Positional<std::string> input_reads_arg(parser, "input_reads",
                                      "input long reads for adapter assembly");

    f_arg version_arg(parser, "version",
                      "display the program version and quit",
                      {"version"});

    args::HelpFlag help(parser, "help",
                        "display this help menu",
                        {'h', "help"});


    parsing_result = GOOD;
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        std::cerr << parser;
        parsing_result = HELP;
        return;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << "\n";
        parsing_result = BAD;
        return;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << "\n";
        parsing_result = BAD;
        return;
    }
    if (argc == 1) {
        std::cerr << parser;
        parsing_result = HELP;
        return;
    }
    if (args::get(version_arg)) {
        parsing_result = VERSION;
        return;
    }

    input_reads = args::get(input_reads_arg);
    if (input_reads.empty()) {
        std::cerr << "Error: input reads are required" << "\n";
        parsing_result = BAD;
        return;
    }
    if (!does_file_exist(input_reads)) {
        std::cerr << "Error: cannot find file: " << input_reads << "\n";
        parsing_result = BAD;
        return;
    }
}


bool Arguments::does_file_exist(std::string filename){
    std::ifstream infile(filename);
    return infile.good();
}
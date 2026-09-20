#include <zlib.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct Fragment {
    std::string seq;
    std::string qual;
};

static std::vector<std::string> split_tab(const std::string &line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        const size_t end = line.find('\t', start);
        fields.emplace_back(line.substr(start, end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return fields;
}

static std::string molecule_name(const std::vector<std::string> &fields) {
    for (size_t i = 11; i < fields.size(); ++i) {
        if (fields[i].rfind("MI:Z:", 0) == 0) return fields[i].substr(5);
    }
    std::string name = fields[0];
    for (int i = 0; i < 2; ++i) {
        const size_t colon = name.rfind(':');
        if (colon == std::string::npos) return fields[0];
        name.resize(colon);
    }
    return name;
}

static char complement(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G';
        case 'G': return 'C'; case 'T': return 'A';
        case 'a': return 't'; case 'c': return 'g';
        case 'g': return 'c'; case 't': return 'a';
        default: return c;
    }
}

static void restore_original_orientation(Fragment &fragment, int flag) {
    if ((flag & 0x10) == 0) return;
    std::reverse(fragment.seq.begin(), fragment.seq.end());
    std::transform(fragment.seq.begin(), fragment.seq.end(),
                   fragment.seq.begin(), complement);
    std::reverse(fragment.qual.begin(), fragment.qual.end());
}

static void write_text(gzFile out, const std::string &text) {
    if (gzwrite(out, text.data(), static_cast<unsigned int>(text.size())) !=
        static_cast<int>(text.size())) {
        throw std::runtime_error("failed to write gzip output");
    }
}

static uint64_t emit_pairs(const std::string &molecule,
                           const std::vector<Fragment> &fragments,
                           gzFile r1, gzFile r2) {
    uint64_t count = 0;
    for (size_t i = 0; i < fragments.size(); ++i) {
        for (size_t j = i + 1; j < fragments.size(); ++j) {
            ++count;
            const std::string name = molecule + "_porec_pair_" + std::to_string(count);
            write_text(r1, "@" + name + "/1\n" + fragments[i].seq +
                           "\n+\n" + fragments[i].qual + "\n");
            write_text(r2, "@" + name + "/2\n" + fragments[j].seq +
                           "\n+\n" + fragments[j].qual + "\n");
        }
    }
    return count;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <name-sorted-monomers.bam> <R1.fastq.gz> <R2.fastq.gz>\n";
        return 2;
    }

    const std::string command = "samtools view \"" + std::string(argv[1]) + "\"";
    FILE *input = popen(command.c_str(), "r");
    gzFile r1 = gzopen(argv[2], "wb1");
    gzFile r2 = gzopen(argv[3], "wb1");
    if (!input || !r1 || !r2) {
        std::cerr << "ERROR: cannot open input or output\n";
        return 1;
    }
    gzbuffer(r1, 1 << 20);
    gzbuffer(r2, 1 << 20);

    char *buffer = nullptr;
    size_t capacity = 0;
    std::string current;
    std::vector<Fragment> fragments;
    uint64_t molecules = 0, pairs = 0;

    try {
        while (getline(&buffer, &capacity, input) >= 0) {
            std::string line(buffer);
            if (!line.empty() && line.back() == '\n') line.pop_back();
            const auto fields = split_tab(line);
            if (fields.size() < 11 || fields[9] == "*") continue;
            const std::string molecule = molecule_name(fields);
            if (!current.empty() && molecule != current) {
                pairs += emit_pairs(current, fragments, r1, r2);
                ++molecules;
                fragments.clear();
            }
            current = molecule;
            Fragment fragment{fields[9], fields[10]};
            if (fragment.qual == "*") fragment.qual.assign(fragment.seq.size(), 'I');
            if (fragment.seq.size() != fragment.qual.size())
                throw std::runtime_error("sequence/quality mismatch for " + fields[0]);
            restore_original_orientation(fragment, std::stoi(fields[1]));
            fragments.emplace_back(std::move(fragment));
        }
        if (!current.empty()) {
            pairs += emit_pairs(current, fragments, r1, r2);
            ++molecules;
        }
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        free(buffer); gzclose(r1); gzclose(r2); pclose(input);
        return 1;
    }

    free(buffer);
    const int input_status = pclose(input);
    const int r1_status = gzclose(r1), r2_status = gzclose(r2);
    if (input_status != 0 || r1_status != Z_OK || r2_status != Z_OK) return 1;
    std::cerr << "Completed: " << molecules << " molecules, " << pairs << " pairs\n";
    return 0;
}

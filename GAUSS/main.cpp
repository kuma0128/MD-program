#include "GAUSS/gaussian.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::size_t samples{1'000'000};
    int bins{200};
    double width{0.02};
    std::uint64_t seed{80};
    std::string output{"output.d"};
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            std::cout << "Usage: gaussian [--samples N] [--bins N] [--width X] "
                         "[--seed N] [--output FILE]\n";
            std::exit(0);
        }
        if (index + 1 >= argc) {
            throw std::invalid_argument("missing value for " + argument);
        }
        const std::string value = argv[++index];
        if (argument == "--samples") {
            options.samples = std::stoull(value);
        } else if (argument == "--bins") {
            options.bins = std::stoi(value);
        } else if (argument == "--width") {
            options.width = std::stod(value);
        } else if (argument == "--seed") {
            options.seed = std::stoull(value);
        } else if (argument == "--output") {
            options.output = value;
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.samples == 0 || options.bins <= 0 || !(options.width > 0.0)) {
        throw std::invalid_argument("samples, bins, and width must be positive");
    }
    return options;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        md::GaussianGenerator gaussian(options.seed);
        std::vector<std::size_t> histogram(static_cast<std::size_t>(2 * options.bins + 1));
        double sum = 0.0;
        double squared_sum = 0.0;

        for (std::size_t sample = 0; sample < options.samples; ++sample) {
            const double value = gaussian();
            int bin = static_cast<int>(std::llround(value / options.width));
            bin = std::clamp(bin, -options.bins, options.bins);
            ++histogram[static_cast<std::size_t>(bin + options.bins)];
            sum += value;
            squared_sum += value * value;
        }

        const double mean = sum / static_cast<double>(options.samples);
        const double variance = squared_sum / static_cast<double>(options.samples) - mean * mean;
        std::ofstream output(options.output);
        if (!output) {
            throw std::runtime_error("cannot open output file: " + options.output);
        }
        output << std::setprecision(12);
        for (int bin = -options.bins; bin <= options.bins; ++bin) {
            output << static_cast<double>(bin) * options.width << ' '
                   << histogram[static_cast<std::size_t>(bin + options.bins)] << '\n';
        }
        output << "# mean variance " << mean << ' ' << variance << '\n';
        std::cout << "wrote " << options.samples << " samples to " << options.output << '\n';
    } catch (const std::exception& error) {
        std::cerr << "gaussian: " << error.what() << '\n';
        return 1;
    }
    return 0;
}

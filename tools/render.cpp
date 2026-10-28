// Portable CPU entry point for the preserved FRAKTALFAKTUR v1.8.6 engine.
// All dynamics, sampling, colours, presets and file formats come from src/.
#include "../src/catalogue.hpp"
#include "../src/focus_presets.hpp"
#include "../src/render_workers.hpp"
#include <iostream>
#include <optional>

namespace {
int integer(const std::string& text, int minimum, int maximum, const std::string& name) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error(name + " requires an integer");
    unsigned long long value = 0;
    try { value = std::stoull(text); }
    catch (...) { throw std::runtime_error(name + " is out of range"); }
    if (value < static_cast<unsigned>(minimum) || value > static_cast<unsigned>(maximum))
        throw std::runtime_error(name + " must be in [" + std::to_string(minimum) + ", " +
                                 std::to_string(maximum) + "]");
    return static_cast<int>(value);
}
void usage() {
    std::cout << "Fractalengine portable CPU renderer (FRAKTALFAKTUR v1.8.6)\n"
        "  fractal-render --list\n"
        "  fractal-render --renderer ID --output image.bmp [options]\n"
        "  fractal-render --view saved.view --output image.bmp [options]\n"
        "Options:\n"
        "  --width N       1..8192 (default 800)\n"
        "  --height N      1..8192 (default 600); at most 16777216 pixels\n"
        "  --preset N      1..4; native Currents bank for RH drivers\n"
        "  --iterations N  32..1024; default from the engine/view\n"
        "  --palette N     0..2; default from the engine/view\n"
        "  --samples N     1..4 per axis (default 1)\n"
        "  --threads N     1..64 (default up to 16 hardware threads)\n"
        "  --save-view P   save resulting parameters in native .view format\n"
        "  --help\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 1) { usage(); return 0; }
        std::string renderer, inputView, output, savedView;
        int width = 800, height = 600, samples = 1;
        unsigned threads = std::min(16u, std::max(1u, std::thread::hardware_concurrency()));
        std::optional<int> selectedPreset, iterations, palette;
        bool list = false, help = false;
        std::set<std::string> seen;
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (!seen.insert(option).second) throw std::runtime_error("Duplicate option: " + option);
            if (option == "--list") { list = true; continue; }
            if (option == "--help") { help = true; continue; }
            if (option != "--renderer" && option != "--view" && option != "--output" &&
                option != "--save-view" && option != "--width" && option != "--height" &&
                option != "--samples" && option != "--threads" && option != "--preset" &&
                option != "--iterations" && option != "--palette")
                throw std::runtime_error("Unknown option: " + option);
            if (++i == argc || std::string(argv[i]).empty() || std::string(argv[i]).rfind("--", 0) == 0)
                throw std::runtime_error("Missing value for " + option);
            const std::string value = argv[i];
            if (option == "--renderer") renderer = value;
            else if (option == "--view") inputView = value;
            else if (option == "--output") output = value;
            else if (option == "--save-view") savedView = value;
            else if (option == "--width") width = integer(value, 1, 8192, option);
            else if (option == "--height") height = integer(value, 1, 8192, option);
            else if (option == "--samples") samples = integer(value, 1, 4, option);
            else if (option == "--threads") threads = unsigned(integer(value, 1, 64, option));
            else if (option == "--preset") selectedPreset = integer(value, 1, 4, option) - 1;
            else if (option == "--iterations") iterations = integer(value, 32, 1024, option);
            else if (option == "--palette") palette = integer(value, 0, 2, option);
        }
        if (help || list) {
            if (argc != 2) throw std::runtime_error("Use --help or --list on its own");
            if (help) { usage(); return 0; }
            for (const auto& entry : catalogue::builtins())
                std::cout << entry.id << '\t' << entry.title << '\n';
            return 0;
        }
        if (renderer.empty() == inputView.empty())
            throw std::runtime_error("Choose exactly one of --renderer and --view");
        if (output.empty()) throw std::runtime_error("--output is required");
        if (std::filesystem::path(output).extension() != ".bmp")
            throw std::runtime_error("--output must have a .bmp extension");
        if (!savedView.empty() && std::filesystem::path(output).lexically_normal() ==
                                      std::filesystem::path(savedView).lexically_normal())
            throw std::runtime_error("Image and saved view paths must differ");
        const std::size_t pixelCount = std::size_t(width) * std::size_t(height);
        if (pixelCount > 16777216) throw std::runtime_error("Image exceeds 16777216 pixels");
        const auto entries = catalogue::builtins();
        amp::State state;
        if (!inputView.empty()) {
            state = amp::readView(inputView);
            const auto found = std::find_if(entries.begin(), entries.end(), [&](const auto& entry) {
                return entry.family == state.family && entry.generation == state.generation;
            });
            if (found == entries.end()) throw std::runtime_error("View is not in the current catalogue");
            renderer = found->id;
        } else {
            const auto found = std::find_if(entries.begin(), entries.end(), [&](const auto& entry) {
                return entry.id == renderer;
            });
            if (found == entries.end()) throw std::runtime_error("Unknown renderer ID: " + renderer);
            state = amp::defaults(found->family, found->generation);
            if (state.family == 20) focus_presets::apply(state, 1, 0);
        }
        if (selectedPreset) {
            if (state.family == 20) focus_presets::apply(state, 1, *selectedPreset);
            else amp::preset(state, *selectedPreset);
        }
        if (iterations) state.iterations = *iterations;
        if (palette) state.palette = *palette;
        const auto scene = amp::scene(state);
        std::vector<uint32_t> pixels(pixelCount);
        std::atomic<int> nextRow{0};
        renderworkers::Pool workers;
        workers.run(std::min(threads, unsigned(height)), [&] {
            for (int y; (y = nextRow.fetch_add(1)) < height;)
                for (int x = 0; x < width; ++x)
                    pixels[std::size_t(y) * width + x] = amp::pixel(x, y, width, height, scene, samples);
        });
        amp::writeBmp(output, pixels, width, height);
        if (!savedView.empty()) amp::writeView(savedView, state);
        std::cout << "RENDER=PASS BACKEND=CPU_DOUBLE RENDERER=" << renderer
                  << " WIDTH=" << width << " HEIGHT=" << height
                  << " THREADS=" << workers.lastWorkerCount() << " OUTPUT=" << output << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RENDER=FAIL " << error.what() << '\n';
        return 1;
    }
}

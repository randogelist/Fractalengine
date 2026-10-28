#pragma once
#include "amp_core.hpp"
#include "frk_extra_catalogue.hpp"

namespace catalogue {
struct Entry {
    std::string id, title, description;
    int family = 0, generation = 1;
};
inline std::vector<Entry> builtins() {
    std::vector<Entry> entries;
    entries.reserve(rhfresh::CandidateCount + frk_extra::driver_count);
    const auto& candidates=rhfresh::registry();
    for(int i=0;i<rhfresh::CandidateCount;++i) {
        const auto& e=candidates[i];
        entries.push_back({e.id,e.name,e.description,20,i+1});
    }
    for (int i=0;i<frk_extra::driver_count;++i) {
        const auto& a=frk_extra::extra_catalogue[i];
        const std::string id="extra_"+std::to_string(a.driver_id);
        entries.push_back({id,std::to_string(19+i)+" / "+a.label,
                           i>=8 ? "QQRH finite Hardy geometry / N=12 / GPU + CPU double" :
                           "Native complex-dynamics generation / GPU + CPU double",25,i+1});
    }
    return entries;
}
inline Entry read(const std::filesystem::path &p) {
    std::ifstream f(p);
    std::string header;
    std::getline(f, header);
    if (!header.empty() && header.back() == '\r')
        header.pop_back();
    if (header != "FRAKTALFAKTUR_RENDERER 1")
        throw std::runtime_error("Bad renderer descriptor: " + p.filename().string());
    Entry e;
    std::string line;
    std::set<std::string> seen;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty() || line[0] == '#')
            continue;
        auto space = line.find(' ');
        if (space == std::string::npos)
            throw std::runtime_error("Malformed renderer descriptor");
        auto key = line.substr(0, space), value = line.substr(space + 1);
        if (!seen.insert(key).second)
            throw std::runtime_error("Duplicate descriptor field");
        if (key == "id")
            e.id = value;
        else if (key == "title")
            e.title = value;
        else if (key == "description")
            e.description = value;
        else if (key == "family") {
            if (value == "amp")
                e.family = 0;
            else if (value == "gp")
                e.family = 1;
            else if (value == "percolation")
                e.family = 2;
            else if (value == "restart")
                e.family = 4;
            else if (value == "cayley")
                e.family = 5;
            else if (value == "dirichlet")
                e.family = 6;
            else if (value == "grh_cayley")
                e.family = 7;
            else if (value == "primitive_pair")
                e.family = 8;
            else if (value == "grh_observable")
                e.family = 9;
            else if (value == "grh_phase_feedback")
                e.family = 10;
            else if (value == "spectral_jet")
                e.family = 11;
            else if (value == "quantum_measurement")
                e.family = 12;
            else if (value == "spherical_fractals")
                e.family = 13;
            else if (value == "spin_chaos")
                e.family = 14;
            else if (value == "project_fractals")
                e.family = 15;
            else if (value == "qqrh_direct")
                e.family = 24;
            else if (value == "qqrh")
                e.family = 23;
            else if (value == "bkqr_rh")
                e.family = 22;
            else if (value == "bkqr_lambda_legacy")
                e.family = 21;
            else if (value == "rh_fresh" || value == "rh_master")
                e.family = 20;
            else if (value == "frk_extra")
                e.family = 25;
            else if (value == "moment_triple")
                e.family = 19;
            else if (value == "moment_quotient")
                e.family = 18;
            else if (value == "moment_cubic")
                e.family = 17;
            else if (value == "riemann_witness")
                e.family = 16;
            else if (value == "crown")
                e.family = 3;
            else
                throw std::runtime_error("Unknown renderer family");
        } else if (key == "generation") {
            size_t end = 0;
            e.generation = std::stoi(value, &end);
            if (end != value.size())
                throw std::runtime_error("Invalid generation");
        } else
            throw std::runtime_error("Unknown descriptor field");
    }
    if (seen.size() != 5 || e.id.empty() || e.title.empty() || e.id.size() > 64 || e.title.size() > 100 ||
        e.description.size() > 400 || e.generation < 1 ||
        e.generation > (e.family == 25 ? frk_extra::driver_count : e.family == 24 ? 6 : e.family == 23 ? 5 : e.family == 22 ? 3 : e.family == 21 ? 1 : e.family == 20 ? rhfresh::CandidateCount : e.family == 19 ? 3 : e.family == 18 ? 1 : e.family == 17 ? 1 : e.family == 16 ? 7 : e.family == 15 ? 5 : e.family == 14 ? 3 : e.family == 13 ? 3 : e.family == 12 ? 1 : e.family == 11 ? 2 : (e.family == 1 || e.family >= 3) ? 4
                        : e.family == 2                  ? 6
                                                         : 1))
        throw std::runtime_error("Invalid renderer descriptor");
    if (e.id.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_-") != std::string::npos)
        throw std::runtime_error("Invalid renderer ID");
    return e;
}
inline std::vector<Entry> load(const std::filesystem::path &root,
                               std::vector<std::string> *warnings = nullptr) {
    std::vector<Entry> entries;
    std::set<std::string> ids;
    std::vector<std::filesystem::path> paths;
    if (std::filesystem::exists(root))
        for (const auto &p : std::filesystem::directory_iterator(root))
            if (p.is_regular_file() && p.path().extension() == ".renderer")
                paths.push_back(p.path());
    // Catalogue prefixes are numbers: 100 must follow 99, not precede 10.
    // Compare digit strings without integer conversion or overflow limits.
    std::sort(paths.begin(), paths.end(), [](const auto &a, const auto &b) {
        const auto an=a.filename().wstring(), bn=b.filename().wstring();
        auto prefix=[](const std::wstring &name) {
            const auto end=name.find_first_not_of(L"0123456789");
            if (end==0 || end==std::wstring::npos || name[end]!=L'_') return std::wstring{};
            const auto first=name.find_first_not_of(L'0');
            return first>=end ? std::wstring(L"0") : name.substr(first,end-first);
        };
        const auto ap=prefix(an), bp=prefix(bn);
        if (ap.empty()!=bp.empty()) return !ap.empty();
        if (ap.size()!=bp.size()) return ap.size()<bp.size();
        if (ap!=bp) return ap<bp;
        return an<bn;
    });
    for (const auto &p : paths) {
        try {
            auto e = read(p);
            // This fresh catalogue accepts only the registered new candidates.
            // Old descriptors can remain beside a copied runner without bringing
            // any previous catalogue entry back into the menu.
            const bool oldCandidate=e.family==20 && e.generation>=1 &&
                e.generation<=rhfresh::CandidateCount &&
                e.id==rhfresh::registry()[e.generation-1].id;
            const bool newCandidate=e.family==25 && e.generation>=1 &&
                e.generation<=frk_extra::driver_count &&
                e.id=="extra_"+std::to_string(frk_extra::extra_catalogue[e.generation-1].driver_id);
            if(!oldCandidate && !newCandidate)
                throw std::runtime_error("Descriptor is not a registered fractal driver: "+p.filename().string());
            if (!ids.insert(e.id).second)
                throw std::runtime_error("Duplicate renderer ID");
            entries.push_back(e);
        } catch (const std::exception &e) {
            if (warnings)
                warnings->push_back(e.what());
        }
    }
    if(entries.empty())return builtins();
    const auto all=builtins();
    for(const auto& entry:all)
        if(ids.insert(entry.id).second)entries.push_back(entry);
    return entries;
}
} // namespace catalogue

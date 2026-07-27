// nga_range_stats.hpp  --  header-only dynamic-range / product-magnitude stats
// ---------------------------------------------------------------------------
// mp-spice / NGA side-task: instrument GROMACS force kernels to record the
// magnitude distribution of the numeric operands and products they form, so the
// same value-distribution / dynamic-range histogram we build for SPICE circuit
// matrices can be built for an MD force computation.
//
// This is the GROMACS-side twin of the SPICE-side DynamicRangeStats /
// product_magnitude_stats in include/sw/mp_spice/klu_study.hpp. Keep the field
// names here aligned with that header so both apps land on the same axes when
// you overlay them (charge/mass/LJ/Coulomb operands and per-interaction
// products vs. the exponent reach of float32 / float16 / posit<16,1> /
// posit<32,2>).
//
// Design goals:
//   * header-only, C++17, no dependencies beyond the standard library
//   * ZERO overhead when NGA_INSTRUMENT is not defined (macros compile to void)
//   * cheap + thread-safe under GROMACS's OpenMP nonbonded loops: each thread
//     accumulates into thread_local storage; a global registry merges on flush
//   * records exactly what a posit regime selection cares about: a log2 (binade)
//     histogram of |value|, plus abs_min / abs_max / signed extremes / zeros
//
// Usage in a kernel (see instrumentation/gromacs_nbnxm_instrument.md):
//     NGA_RECORD("coulomb_qq",   qq);          // operand product q_i*q_j
//     NGA_RECORD("coulomb_term", qq * rinv);   // full Coulomb term magnitude
//     NGA_RECORD("lj_c6",        c6);
//     NGA_RECORD("lj_term",      FrLJ12 - FrLJ6);
// and once, after the mdrun integration loop finishes:
//     NGA_FLUSH("nga_gromacs_stats");   // writes nga_gromacs_stats.{json,csv}
// ---------------------------------------------------------------------------
#ifndef NGA_RANGE_STATS_HPP
#define NGA_RANGE_STATS_HPP

#ifdef NGA_INSTRUMENT

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace nga {

// log2 binade histogram spanning exponents [-MINEXP, +MAXEXP]. IEEE double
// normal exponents run about [-1022, 1023]; we bucket a generous, MD-relevant
// window and lump anything outside into underflow/overflow guard bins.
constexpr int kMinExp = -64;   // 2^-64  ~ 5.4e-20
constexpr int kMaxExp = 64;    // 2^+64  ~ 1.8e+19
constexpr int kNBins  = kMaxExp - kMinExp + 3;  // +3: underflow, overflow, zero

struct DynamicRangeStats {
    std::string name;
    uint64_t    count      = 0;
    uint64_t    count_zero = 0;
    uint64_t    count_nan  = 0;
    double      abs_min    =  INFINITY;   // smallest nonzero |x|
    double      abs_max    =  0.0;        // largest |x|
    double      signed_min =  INFINITY;
    double      signed_max = -INFINITY;
    std::array<uint64_t, kNBins> log2_hist{};  // index 0 = zero bin

    inline void record(double x) {
        ++count;
        if (x == 0.0) { ++count_zero; ++log2_hist[0]; return; }
        if (std::isnan(x)) { ++count_nan; return; }
        if (x < signed_min) signed_min = x;
        if (x > signed_max) signed_max = x;
        double a = std::fabs(x);
        if (a < abs_min) abs_min = a;
        if (a > abs_max) abs_max = a;
        int e;
        std::frexp(a, &e);              // a = m * 2^e, m in [0.5,1) -> e is binade
        int idx;
        if (e < kMinExp)       idx = kNBins - 2;   // underflow guard
        else if (e > kMaxExp)  idx = kNBins - 1;   // overflow guard
        else                   idx = 1 + (e - kMinExp);
        ++log2_hist[idx];
    }

    void merge(const DynamicRangeStats& o) {
        count += o.count; count_zero += o.count_zero; count_nan += o.count_nan;
        if (o.abs_min < abs_min) abs_min = o.abs_min;
        if (o.abs_max > abs_max) abs_max = o.abs_max;
        if (o.signed_min < signed_min) signed_min = o.signed_min;
        if (o.signed_max > signed_max) signed_max = o.signed_max;
        for (int i = 0; i < kNBins; ++i) log2_hist[i] += o.log2_hist[i];
    }

    double decades() const {
        return (abs_min > 0 && std::isfinite(abs_min) && abs_max > 0)
                   ? std::log10(abs_max / abs_min) : 0.0;
    }
};

// ---- registry -------------------------------------------------------------
// Per-thread tables avoid any locking on the hot path. A process-wide list of
// all thread tables is merged at flush time.
class Registry {
  public:
    static Registry& global() { static Registry r; return r; }

    // thread-local table; registered with the global registry on first use
    std::map<std::string, DynamicRangeStats>& local_table() {
        thread_local std::map<std::string, DynamicRangeStats>* t = [this] {
            auto* m = new std::map<std::string, DynamicRangeStats>();
            std::lock_guard<std::mutex> lk(mu_);
            tables_.push_back(m);
            return m;
        }();
        return *t;
    }

    void flush(const std::string& stem) {
        std::lock_guard<std::mutex> lk(mu_);
        std::map<std::string, DynamicRangeStats> merged;
        for (auto* t : tables_)
            for (auto& kv : *t) {
                auto& d = merged[kv.first];
                d.name = kv.first;
                d.merge(kv.second);
            }
        write_json(stem + ".json", merged);
        write_csv(stem + ".csv", merged);
        std::fprintf(stderr, "[nga] wrote %s.json / .csv (%zu operands)\n",
                     stem.c_str(), merged.size());
    }

  private:
    std::mutex mu_;
    std::vector<std::map<std::string, DynamicRangeStats>*> tables_;

    static void write_json(const std::string& path,
                           std::map<std::string, DynamicRangeStats>& m) {
        FILE* f = std::fopen(path.c_str(), "w");
        if (!f) return;
        std::fprintf(f, "{\n  \"minexp\": %d, \"maxexp\": %d,\n  \"operands\": {\n",
                     kMinExp, kMaxExp);
        size_t k = 0;
        for (auto& kv : m) {
            auto& d = kv.second;
            std::fprintf(f, "    \"%s\": {\n", kv.first.c_str());
            std::fprintf(f, "      \"count\": %llu, \"count_zero\": %llu, \"count_nan\": %llu,\n",
                         (unsigned long long)d.count, (unsigned long long)d.count_zero,
                         (unsigned long long)d.count_nan);
            std::fprintf(f, "      \"abs_min\": %.9g, \"abs_max\": %.9g, \"decades\": %.4f,\n",
                         d.abs_min, d.abs_max, d.decades());
            std::fprintf(f, "      \"signed_min\": %.9g, \"signed_max\": %.9g,\n",
                         d.signed_min, d.signed_max);
            std::fprintf(f, "      \"log2_hist\": [");
            for (int i = 0; i < kNBins; ++i)
                std::fprintf(f, "%llu%s", (unsigned long long)d.log2_hist[i],
                             i + 1 < kNBins ? "," : "");
            std::fprintf(f, "]\n    }%s\n", (++k < m.size()) ? "," : "");
        }
        std::fprintf(f, "  }\n}\n");
        std::fclose(f);
    }

    static void write_csv(const std::string& path,
                          std::map<std::string, DynamicRangeStats>& m) {
        FILE* f = std::fopen(path.c_str(), "w");
        if (!f) return;
        // one row per (operand, binade) so plot_instrumented.py can pivot it
        std::fprintf(f, "operand,binade_exp,count\n");
        for (auto& kv : m) {
            auto& d = kv.second;
            for (int i = 1; i < kNBins - 2; ++i) {
                if (d.log2_hist[i] == 0) continue;
                std::fprintf(f, "%s,%d,%llu\n", kv.first.c_str(),
                             kMinExp + (i - 1), (unsigned long long)d.log2_hist[i]);
            }
        }
        std::fclose(f);
    }
};

inline void record(const char* name, double x) {
    auto& tbl = Registry::global().local_table();
    auto it = tbl.find(name);
    if (it == tbl.end()) it = tbl.emplace(name, DynamicRangeStats{}).first;
    it->second.record(x);
}
inline void flush(const std::string& stem) { Registry::global().flush(stem); }

}  // namespace nga

#define NGA_RECORD(name, value) ::nga::record((name), (double)(value))
#define NGA_FLUSH(stem)         ::nga::flush((stem))

#else  // NGA_INSTRUMENT not defined -> compile to nothing

#define NGA_RECORD(name, value) ((void)0)
#define NGA_FLUSH(stem)         ((void)0)

#endif  // NGA_INSTRUMENT
#endif  // NGA_RANGE_STATS_HPP

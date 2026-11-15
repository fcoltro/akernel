// Minimal test support: counted checks and a deterministic random source.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#if defined(_MSC_VER) && defined(_DEBUG)
#include <crtdbg.h>
#endif

namespace akt {

#if defined(_MSC_VER) && defined(_DEBUG)
/// In a Debug build a failed CRT assertion (a checked iterator, say) opens a
/// modal dialog and the program waits for a click, so an unattended test run
/// hangs. Report on stderr and stop with an error code instead.
inline const bool kNoCrtDialogs = [] {
    for (int type : {_CRT_WARN, _CRT_ERROR, _CRT_ASSERT}) {
        _CrtSetReportMode(type, _CRTDBG_MODE_FILE | _CRTDBG_MODE_DEBUG);
        _CrtSetReportFile(type, _CRTDBG_FILE_STDERR);
    }
    _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
    return true;
}();
#endif

inline int& failures() { static int n = 0; return n; }
inline int& checks()   { static int n = 0; return n; }

inline int finish(const char* name) {
    std::printf("%s: %d checks, %d failures\n", name, checks(), failures());
    return failures() == 0 ? 0 : 1;
}

/// SplitMix64. Same sequence on every platform, unlike <random> distributions.
class Rng {
    std::uint64_t s_;
public:
    explicit Rng(std::uint64_t seed) : s_(seed) {}
    std::uint64_t next() {
        std::uint64_t z = (s_ += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    /// Uniform double in [0, 1) with 53 random bits.
    double unit() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
    double range(double a, double b) { return a + (b - a) * unit(); }
};

} // namespace akt

#define AK_CHECK(cond)                                                        \
    do {                                                                      \
        ++akt::checks();                                                      \
        if (!(cond)) {                                                        \
            ++akt::failures();                                                \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);       \
        }                                                                     \
    } while (false)

#define AK_CHECK_NEAR(a, b, tol)                                              \
    do {                                                                      \
        ++akt::checks();                                                      \
        const double ak_a_ = (a), ak_b_ = (b), ak_t_ = (tol);                 \
        if (!(std::fabs(ak_a_ - ak_b_) <= ak_t_)) {                           \
            ++akt::failures();                                                \
            std::printf("FAIL %s:%d: |%s - %s| = %.3e > %.3e\n", __FILE__,    \
                        __LINE__, #a, #b, std::fabs(ak_a_ - ak_b_), ak_t_);   \
        }                                                                     \
    } while (false)

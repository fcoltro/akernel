// Tests for ak::num::on_undecided (include/akernel/num/escalation_policy.hpp).
//
// The same cases are checked twice: at compile time with static_assert, and
// at run time through volatile inputs so the compiler cannot fold them. The
// two can differ (MSVC's constant evaluator orders NaN below everything), and
// the policy must give the same answer in both.
//
// Build and run from the project root, in a Visual Studio developer prompt:
//   cl /nologo /std:c++20 /W4 /WX /EHsc /I include tests\num\escalation_policy_test.cpp
//   escalation_policy_test.exe
// Exit code 0 means every case passed.

#include <akernel/num/escalation_policy.hpp>

#include <cstdio>
#include <limits>

namespace {

using namespace ak::num;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr std::uint32_t kCap = 4096;

struct Case {
    const char*    name;
    PredicateState state;
    Action         expected;
};

// Fields of PredicateState: tier, precision_bits, max_precision_bits,
// enclosure_width, local_feature_size, shared_provenance, symbolic_tried.
constexpr Case kCases[] = {
    // 1. Shared ancestor, symbolic proof not yet tried: prove first, at any tier.
    {"shared ancestor at Double",
     {Tier::Double, 0, kCap, 1e-12, 1e-3, true, false}, Action::ProveBySymbol},
    {"shared ancestor at the precision cap",
     {Tier::MultiPrecision, kCap, kCap, 1e-12, 1e-3, true, false}, Action::ProveBySymbol},

    // 2. Precision still available: escalate.
    {"Double, no shared ancestor",
     {Tier::Double, 0, kCap, 1e-12, 1e-3, false, false}, Action::Escalate},
    {"Expansion, symbolic proof already failed",
     {Tier::Expansion, 0, kCap, 1e-12, 1e-3, true, true}, Action::Escalate},
    {"MultiPrecision below the cap",
     {Tier::MultiPrecision, 1024, kCap, 1e-12, 1e-3, false, false}, Action::Escalate},

    // 3. Precision exhausted, symbolic proof never tried: try it once.
    {"cap reached, symbolic untried",
     {Tier::MultiPrecision, kCap, kCap, 1e-12, 1e-3, false, false}, Action::ProveBySymbol},

    // 4. Nothing proved: merge only strictly inside half the feature size.
    {"enclosure well inside the bound",
     {Tier::MultiPrecision, kCap, kCap, 1e-12, 1e-3, false, true}, Action::MergeRecorded},
    {"enclosure outside the bound",
     {Tier::MultiPrecision, kCap, kCap, 6e-4, 1e-3, false, true}, Action::Fail},
    {"enclosure exactly on the bound (bound is strict)",
     {Tier::MultiPrecision, kCap, kCap, 5e-4, 1e-3, false, true}, Action::Fail},
    {"feature size zero",
     {Tier::MultiPrecision, kCap, kCap, 1e-12, 0.0, false, true}, Action::Fail},
    {"feature size NaN",
     {Tier::MultiPrecision, kCap, kCap, 1e-12, kNaN, false, true}, Action::Fail},
    {"enclosure width NaN",
     {Tier::MultiPrecision, kCap, kCap, kNaN, 1e-3, false, true}, Action::Fail},

    // Termination: past the cap the policy must not ask for more precision.
    {"precision above the cap, everything tried",
     {Tier::MultiPrecision, 2 * kCap, kCap, 1e-12, 1e-3, true, true}, Action::MergeRecorded},
};

constexpr bool all_cases_pass_at_compile_time() {
    for (const Case& c : kCases) {
        if (on_undecided(c.state) != c.expected) return false;
    }
    return true;
}
static_assert(all_cases_pass_at_compile_time(),
              "on_undecided disagrees with the case table at compile time");

const char* to_string(Action a) {
    switch (a) {
        case Action::Escalate:      return "Escalate";
        case Action::ProveBySymbol: return "ProveBySymbol";
        case Action::MergeRecorded: return "MergeRecorded";
        case Action::Fail:          return "Fail";
    }
    return "?";
}

} // namespace

int main() {
    int failures = 0;
    for (const Case& c : kCases) {
        // Copy through a volatile object so the call is evaluated at run time.
        volatile PredicateState runtime_state = c.state;
        const PredicateState s = {
            runtime_state.tier,               runtime_state.precision_bits,
            runtime_state.max_precision_bits, runtime_state.enclosure_width,
            runtime_state.local_feature_size, runtime_state.shared_provenance,
            runtime_state.symbolic_tried};
        const Action got = on_undecided(s);
        if (got != c.expected) {
            ++failures;
            std::printf("FAIL  %s: expected %s, got %s\n", c.name,
                        to_string(c.expected), to_string(got));
        }
    }
    const int total = static_cast<int>(sizeof(kCases) / sizeof(kCases[0]));
    std::printf("%d of %d run-time cases passed\n", total - failures, total);
    return failures == 0 ? 0 : 1;
}

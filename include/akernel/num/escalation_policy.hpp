// akernel — escalation policy for undecided predicates.
//
// See docs/BLUEPRINT.md §3.2 (filter cascade) and §3.3 step 5.
//
// Every geometric predicate returns a proven Sign or Sign::Undecided.
// When a tier cannot prove the sign, the cascade asks this policy what to
// do next. The policy is a pure function of PredicateState so that it can be
// unit-tested and so that the kernel's behaviour is reproducible.
#pragma once

#include <cstdint>

namespace ak::num {

/// Arithmetic tiers of the filter cascade, cheapest first.
enum class Tier : std::uint8_t {
    Double,         ///< double with a static error bound
    Interval,       ///< interval arithmetic, outward rounding
    Expansion,      ///< double-double / quad-double
    MultiPrecision, ///< ball arithmetic, precision doubles on each retry
};

/// Result of a predicate. Undecided means "not proven", never "probably zero".
enum class Sign : std::int8_t { Negative = -1, Zero = 0, Positive = 1, Undecided = 2 };

/// What the cascade knows at the moment a tier returns Sign::Undecided.
struct PredicateState {
    Tier          tier;               ///< tier that just failed to decide
    std::uint32_t precision_bits;     ///< working precision (MultiPrecision only)
    std::uint32_t max_precision_bits; ///< configured cap for MultiPrecision
    double        enclosure_width;    ///< width of the enclosure that straddles zero
    double        local_feature_size; ///< distance to the nearest non-incident entity
    bool          shared_provenance;  ///< operands have a common constructive ancestor
    bool          symbolic_tried;     ///< provenance / exact algebra already attempted
};

/// What the cascade should do next.
enum class Action : std::uint8_t {
    Escalate,      ///< retry at the next tier, or double precision_bits
    ProveBySymbol, ///< ask the provenance engine / exact algebra for identity
    MergeRecorded, ///< declare the operands coincident and write enclosure_width
                   ///< to the tolerance ledger of the affected entity
    Fail,          ///< return Indeterminate; the enclosing transaction rolls back
};

/// Decide the next step for a predicate that is still undecided.
///
/// Invariants the kernel relies on (the audit will assert them):
///   - MergeRecorded is only legal when
///       enclosure_width < 0.5 * local_feature_size      (blueprint T3)
///   - The policy must terminate: it may not return Escalate once
///       tier == MultiPrecision && precision_bits >= max_precision_bits
///
/// Order of preference:
///   1. Operands share an ancestor: ask for a symbolic proof first. It is
///      cheap, and for a true coincidence more precision can never decide.
///   2. Precision is still available: escalate.
///   3. Precision is exhausted: try symbolic proof once, even without a
///      known shared ancestor (exact algebra may still settle analytic pairs).
///   4. Nothing proved: merge only inside the half-feature-size bound, and
///      record it. Otherwise fail. A NaN or zero feature size fails.
[[nodiscard]] constexpr Action on_undecided(const PredicateState& s) noexcept {
    if (s.shared_provenance && !s.symbolic_tried) return Action::ProveBySymbol;

    const bool precision_left =
        s.tier != Tier::MultiPrecision || s.precision_bits < s.max_precision_bits;
    if (precision_left) return Action::Escalate;

    if (!s.symbolic_tried) return Action::ProveBySymbol;

    // NaN is rejected explicitly. MSVC's constant evaluator orders NaN as
    // less than everything, so relying on `nan < x == false` would give a
    // different answer at compile time than at run time.
    const bool bound_known = s.enclosure_width == s.enclosure_width &&
                             s.local_feature_size == s.local_feature_size;
    const bool within_bound =
        bound_known && s.enclosure_width < 0.5 * s.local_feature_size;
    return within_bound ? Action::MergeRecorded : Action::Fail;
}

} // namespace ak::num

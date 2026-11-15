// Expansion arithmetic must be exact. Tested with identities whose exact
// value is known and which plain doubles get wrong.
#include <akernel/num/expansion.hpp>

#include "check.hpp"

using ak::num::Expansion;

int main() {
    // Catastrophic cancellation in doubles, exact here.
    {
        const double big = 0x1.0p60;
        AK_CHECK(((big + 1.0) - big) != 1.0); // the double result is wrong
        const Expansion r = (Expansion(big) + Expansion(1.0)) - Expansion(big);
        AK_CHECK(r.sign() == 1);
        AK_CHECK(r.estimate() == 1.0);
        AK_CHECK(r.size() == 1);
    }
    {
        const Expansion r = Expansion(1e16) + Expansion(1.0) - Expansion(1e16) - Expansion(1.0);
        AK_CHECK(r.sign() == 0);
    }

    // Zero handling.
    AK_CHECK(Expansion(0.0).sign() == 0);
    AK_CHECK((Expansion(3.0) * Expansion(0.0)).sign() == 0);
    AK_CHECK((Expansion(3.0) - Expansion(3.0)).sign() == 0);
    AK_CHECK((-Expansion(2.0)).sign() == -1);

    // Tiny differences keep their sign.
    {
        const double tiny = 0x1.0p-200;
        AK_CHECK(((Expansion(1.0) + Expansion(tiny)) - Expansion(1.0)).sign() == 1);
        AK_CHECK(((Expansion(1.0) - Expansion(tiny)) - Expansion(1.0)).sign() == -1);
    }

    akt::Rng rng(7);
    for (int i = 0; i < 5000; ++i) {
        const double a = rng.range(-1e6, 1e6), b = rng.range(-1e-3, 1e-3);
        const double c = rng.range(-1e3, 1e3);
        const Expansion A(a), B(b), C(c);

        // (a + b)^2 - a^2 - 2ab - b^2 == 0
        const Expansion s = A + B;
        const Expansion binom = s * s - A * A - Expansion(2.0) * A * B - B * B;
        AK_CHECK(binom.sign() == 0);

        // (a - b)(a + b) == a^2 - b^2
        AK_CHECK(((A - B) * (A + B) - (A * A - B * B)).sign() == 0);

        // Distributivity and associativity hold exactly.
        AK_CHECK((A * (B + C) - (A * B + A * C)).sign() == 0);
        AK_CHECK(((A + B) + C - (A + (B + C))).sign() == 0);
        AK_CHECK(((A * B) * C - A * (B * C)).sign() == 0);

        // Sign agrees with the obvious answer for a well-separated case.
        AK_CHECK((A - Expansion(a + 1.0)).sign() == -1);
    }

    // A determinant that is exactly zero: rows are exact multiples.
    {
        const double x = 0.1, y = 0.7; // not exactly representable in decimal, exact as doubles
        const Expansion det = Expansion(x) * Expansion(2.0 * y) - Expansion(y) * Expansion(2.0 * x);
        AK_CHECK(det.sign() == 0);
    }

    return akt::finish("expansion_test");
}

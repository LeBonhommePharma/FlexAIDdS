// Regression: the map_par RESERVATION INVARIANT that prevents dangling
// atoms[].par pointers.
//
// WHAT THIS TEST COVERS, AND WHAT IT DOES NOT. add2_optimiz_vec() stores RAW
// POINTERS into FA->map_par (atoms[].par at add2_optimiz_vec.cpp:71/136/164/176/
// 187/215, plus map_par_sidechain_first/last), and realloc_par() RELOCATES that
// array. The function is called up to four times per setup ("" x N, then "SC",
// "NM"), so a growth on any later call orphans every pointer the earlier calls
// captured. reserve_par(FA, MAX_PAR) before the first call makes the array
// immovable for the whole setup; that is a LIFETIME INVARIANT held by call-site
// ordering, not by the type system, since atoms[].par is a raw pointer
// (flexaid.h:303). This test pins that invariant.
//
// IT IS NOT A CRASH REGRESSION TEST, and the distinction was established by
// measurement rather than by reading the existing comments. atoms[].par is
// assigned 8 times and DEREFERENCED ZERO times anywhere in LIB, so a stale
// atoms[].par cannot fault. The exit-139 crash on 1R55/1R58 was ASan-localised
// to atoms[].optres -- a different array, dereferenced at vcfunction.cpp:796-798
// and GISTEvaluator.cpp:206, relocated by build_rotamers.cpp:308. The
// "dangling atoms[].par dereferenced in populate_chromosomes()" narrative in
// top.cpp:2655-2668 is NOT the crash mechanism, and this file deliberately does
// not repeat it.
//
// So the optres lifetime -- the one that actually crashed -- remains UNCOVERED.
// A test for it needs reserve_optres plus the build_rotamers realloc path and is
// the more important missing test of the two.
//
// WHY IT EXISTS *NOW* AND NOT THEN. An external audit found a THIRD call site,
// direct_input.cpp, with the same ""/"SC"/"NM" sequence and NO reservation at
// all. It was found by reading source, not by a test: reserve_par,
// reserve_optres and atoms_with_optres appear in ZERO of this tree's test files,
// so nothing could have caught an unguarded site. This test closes that gap for
// reserve_par.
//
// MUTATION-TESTED. The second case is a negative control that drives the
// assertion to FAIL, using the engine's own FLEXAIDDS_PAR_RESERVE=0 escape
// (add2_optimiz_vec.cpp:322) plus direct stock growth. A test that cannot fail
// is not evidence, and this project has shipped several of those.

#include <gtest/gtest.h>

// flexaid.h defines E (Euler) as a macro — must come after gtest headers.
#include "flexaid.h"

#include <cstdlib>
#include <cstring>

// Defined in LIB/add2_optimiz_vec.cpp
void reserve_par(FA_Global* FA, int need);
void realloc_par(FA_Global* FA, int* MIN_PAR);

namespace {

// Minimal FA_Global carrying only the par-array fields. realloc() on a null
// pointer is malloc(), so a zeroed struct is a valid starting state.
struct ParFixture {
    FA_Global FA{};
    ParFixture() {
        std::memset(&FA, 0, sizeof(FA));
        FA.MIN_PAR = 0;
        FA.map_par = nullptr;
        FA.npar = 0;
    }
    ~ParFixture() { std::free(FA.map_par); }
};

// The invariant the two guarded call sites rely on: once reserve_par has been
// asked for MAX_PAR, no later demand at or below MAX_PAR may relocate the array.
TEST(ParReservation, ReserveMaxParMakesTheArrayImmovable) {
    ParFixture f;
    ASSERT_EQ(f.FA.map_par, nullptr) << "fixture must start unallocated";

    reserve_par(&f.FA, MAX_PAR);
    ASSERT_NE(f.FA.map_par, nullptr) << "reserve_par must allocate";
    ASSERT_GT(f.FA.MIN_PAR, MAX_PAR)
        << "reserve_par's loop runs while MIN_PAR <= need, so MIN_PAR must end above need";

    const optmap* base_after_reserve = f.FA.map_par;
    const int min_par_after_reserve = f.FA.MIN_PAR;

    // Every later demand a setup can make, up to the engine's declared gene
    // ceiling. This is the sequence the ""/"SC"/"NM" calls generate.
    for (int need = 1; need <= MAX_PAR; ++need) {
        reserve_par(&f.FA, need);
        ASSERT_EQ(f.FA.map_par, base_after_reserve)
            << "map_par RELOCATED at need=" << need
            << " — any atoms[].par captured before this point is now dangling";
    }
    EXPECT_EQ(f.FA.MIN_PAR, min_par_after_reserve)
        << "no growth should have been necessary after reserving MAX_PAR";
}

// NEGATIVE CONTROL. Without the reservation, stock growth relocates the array —
// so the assertion above is capable of failing, and this documents the exposure
// the fix removes. If this test ever passes trivially (base never moves), the
// allocator is masking the defect and the positive test proves nothing.
TEST(ParReservation, WithoutReservationStockGrowthRelocatesTheArray) {
    ParFixture f;

    // Stock behaviour: start small, grow in fixed steps, exactly as unpatched
    // runs did before the reservation existed.
    realloc_par(&f.FA, &f.FA.MIN_PAR);
    ASSERT_NE(f.FA.map_par, nullptr);
    const optmap* first_base = f.FA.map_par;

    bool relocated = false;
    for (int i = 0; i < 32 && !relocated; ++i) {
        realloc_par(&f.FA, &f.FA.MIN_PAR);
        if (f.FA.map_par != first_base) relocated = true;
    }
    EXPECT_TRUE(relocated)
        << "stock growth never moved map_par in 32 steps — the positive test above "
           "would then be vacuous, because there is nothing for reserve_par to prevent";
}

// The escape hatch must actually disable the reservation, otherwise the
// baseline-exposure measurements taken with FLEXAIDDS_PAR_RESERVE=0 were
// measuring the patched path and are void.
TEST(ParReservation, ParReserveZeroDisablesTheReservation) {
    ParFixture f;
    setenv("FLEXAIDDS_PAR_RESERVE", "0", 1);
    reserve_par(&f.FA, MAX_PAR);
    const bool disabled = (f.FA.map_par == nullptr && f.FA.MIN_PAR == 0);
    unsetenv("FLEXAIDDS_PAR_RESERVE");
    EXPECT_TRUE(disabled)
        << "FLEXAIDDS_PAR_RESERVE=0 must make reserve_par a no-op; it allocated instead, "
           "so any 'baseline exposure' measured with it is confounded by the reservation";
}

}  // namespace

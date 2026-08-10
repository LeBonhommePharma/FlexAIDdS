// tests/test_contacts_epoch.cpp
// Regression coverage for the FLEXAIDDS_CONTACTS_EPOCH stamp-buffer invariant.
//
// THE BUG THIS PINS
// -----------------
// vcfunction() clears a per-atom "already visited this eval" flag array
// (FA_Global::contacts) on every evaluation. Under FLEXAIDDS_CONTACTS_EPOCH the
// O(MAX_ATOM_NUMBER) memset is replaced by an O(1) epoch bump: a slot stamped
// with the current epoch means "visited".
//
// The epoch counter originally lived in FA_Global. But the threaded GA
// evaluation path re-snapshots its per-thread FA_Global from the master
// *every generation* (LIB/gaboom.cpp, `tl_fa[t] = *FA;`) while the stamp buffer
// it re-points at (ParEvalWS::tl_contacts) is allocated ONCE and stays resident
// across generations. The counter was therefore rewound at every generation
// boundary against a buffer that still held the previous generation's
// high-water stamps. Stale stamps then aliased live epoch values, atoms were
// wrongly treated as already-visited, contacts were skipped, and CF was
// silently wrong — no crash, no warning, no diagnostic.
//
// The fix moves the counter INSIDE the buffer it stamps (slot
// CONTACTS_EPOCH_SLOT), so a struct copy cannot separate them.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

// PROOF THAT THIS TEST HAS TEETH
// ------------------------------
// Compiling this file with -DCONTACTS_EPOCH_PREFIX_LAYOUT=1 re-creates the
// pre-fix layout (epoch carried by the per-thread struct instead of by the
// buffer). SurvivesPerGenerationSnapshotRefresh MUST fail in that mode.
// CMake builds that second binary as `test_contacts_epoch_prefix_layout` and
// registers it with WILL_FAIL, so ctest continuously proves both directions:
//   - shipping layout  -> the invariant holds;
//   - pre-fix layout   -> this test detects the corruption.
// A "fix" that silently stopped detecting the bug would turn the WILL_FAIL
// target green-when-it-should-be-red and break the suite.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/flexaid.h"

#include <vector>

#ifndef CONTACTS_EPOCH_PREFIX_LAYOUT
#define CONTACTS_EPOCH_PREFIX_LAYOUT 0
#endif

namespace {

// ---------------------------------------------------------------------------
// The per-thread snapshot that the GA refreshes at every generation boundary.
// ---------------------------------------------------------------------------
// In the shipping layout this is exactly FA_Global: the epoch is not part of it,
// it lives in the buffer `fa.contacts` points at. In the pre-fix layout the
// snapshot also carries the counter — and copying the snapshot then rewinds it.
struct Snapshot {
	FA_Global fa{};
#if CONTACTS_EPOCH_PREFIX_LAYOUT
	int contacts_epoch = 0;
#endif
};

// Begin one evaluation; returns the epoch to stamp with.
inline int begin_eval(Snapshot& s)
{
#if CONTACTS_EPOCH_PREFIX_LAYOUT
	return ++s.contacts_epoch;                    // pre-fix: counter in the struct
#else
	return contacts_epoch_begin(s.fa.contacts);   // shipping: counter in the buffer
#endif
}

}  // namespace

namespace {

// ---------------------------------------------------------------------------
// Shared workload model
// ---------------------------------------------------------------------------
// Deterministic stand-in for a GA run: kGenerations generations, each of
// kEvalsPerGen evaluations, each evaluation touching a *different* subset of
// atoms. The varying subset matters — it is what lets a stamp written late in
// generation N survive untouched into generation N+1 and collide with the
// re-issued epoch there. A model that stamped every atom on every eval would
// overwrite the evidence and pass even on the broken implementation.
//
// Partition rule: atom `a` is touched by exactly one eval per generation,
// namely eval (a % kEvalsPerGen).

constexpr int kGenerations  = 4;
constexpr int kEvalsPerGen  = 32;
constexpr int kAtoms        = 256;

inline bool touched_by(int atom, int eval) { return (atom % kEvalsPerGen) == eval; }

}  // namespace

// ---------------------------------------------------------------------------
// 1. THE REGRESSION TEST — production code path.
// ---------------------------------------------------------------------------
// Mirrors LIB/gaboom.cpp: a resident stamp buffer, and an FA_Global snapshot
// refreshed from the master at every generation boundary. No stale stamp may
// ever alias the live epoch.
TEST(ContactsEpoch, SurvivesPerGenerationSnapshotRefresh)
{
	// Allocated ONCE and resident across every generation, exactly like
	// ParEvalWS::tl_contacts.
	std::vector<int> resident(CONTACTS_BUFFER_SIZE, 0);

	Snapshot master{};   // the master FA; its scalars do not advance during GA
	Snapshot tl{};       // the per-thread snapshot

	int evals = 0;

	for (int gen = 0; gen < kGenerations; ++gen) {
		// ── the generation-boundary refresh (LIB/gaboom.cpp) ──
		tl = master;
		tl.fa.contacts = resident.data();

		for (int eval = 0; eval < kEvalsPerGen; ++eval) {
			const int epoch = begin_eval(tl);
			int* contacts = tl.fa.contacts;
			++evals;

			// Start of a fresh evaluation: NOTHING may read as visited.
			for (int a = 0; a < kAtoms; ++a) {
				ASSERT_NE(contacts[a], epoch)
					<< "stale stamp aliased the live epoch — a contact would be "
					   "silently skipped and CF would be wrong. gen=" << gen
					<< " eval=" << eval << " atom=" << a << " epoch=" << epoch;
			}

			// Stamp this evaluation's subset and confirm it reads back visited.
			for (int a = 0; a < kAtoms; ++a) {
				if (!touched_by(a, eval)) continue;
				contacts[a] = epoch;
				ASSERT_EQ(contacts[a], epoch) << "atom=" << a;
			}
		}
	}

	// The counter advanced once per evaluation and was never rewound by the
	// per-generation struct copy.
	EXPECT_EQ(contacts_epoch_current(resident.data()), evals);
	EXPECT_EQ(evals, kGenerations * kEvalsPerGen);
}

// ---------------------------------------------------------------------------
// 2. MECHANISM PIN — the pre-fix layout, reproduced.
// ---------------------------------------------------------------------------
// Same workload, same resident buffer, but with the epoch stored in the
// FA-like struct (the pre-fix layout). This MUST still corrupt. If it ever
// stops corrupting, the workload model no longer exercises the bug and test 1
// above has become vacuous — which is precisely the failure mode a regression
// test is supposed to be protected against.
namespace {
struct LegacyFA {
	int* contacts       = nullptr;
	int  contacts_epoch = 0;   // <-- the defect: rewound by `tl = master`
};
}  // namespace

TEST(ContactsEpoch, LegacyEpochInStructRewindsAndAliases)
{
	std::vector<int> resident(CONTACTS_BUFFER_SIZE, 0);

	LegacyFA master;   // master epoch stays at 0 for the whole GA
	LegacyFA tl;

	int aliases = 0;
	int first_alias_gen = -1;

	for (int gen = 0; gen < kGenerations; ++gen) {
		tl = master;                    // <-- rewinds contacts_epoch to 0
		tl.contacts = resident.data();

		for (int eval = 0; eval < kEvalsPerGen; ++eval) {
			const int epoch = ++tl.contacts_epoch;

			for (int a = 0; a < kAtoms; ++a) {
				if (tl.contacts[a] == epoch) {   // false "already visited"
					++aliases;
					if (first_alias_gen < 0) first_alias_gen = gen;
				}
			}
			for (int a = 0; a < kAtoms; ++a) {
				if (touched_by(a, eval)) tl.contacts[a] = epoch;
			}
		}
	}

	EXPECT_GT(aliases, 0)
		<< "the pre-fix layout no longer reproduces the stale-stamp alias; the "
		   "workload model above has stopped exercising the bug and "
		   "SurvivesPerGenerationSnapshotRefresh is now vacuous";
	EXPECT_EQ(first_alias_gen, 1)
		<< "corruption must appear at the FIRST generation boundary";
}

// ---------------------------------------------------------------------------
// 3. A struct copy cannot rewind the epoch.
// ---------------------------------------------------------------------------
TEST(ContactsEpoch, EpochTravelsWithBufferNotWithStruct)
{
	std::vector<int> buf(CONTACTS_BUFFER_SIZE, 0);

	FA_Global fa{};
	fa.contacts = buf.data();

	for (int i = 0; i < 100; ++i) contacts_epoch_begin(fa.contacts);
	ASSERT_EQ(contacts_epoch_current(fa.contacts), 100);

	// A shallow copy of a *stale* FA_Global — the pre-fix rewind vector.
	FA_Global stale{};
	FA_Global snapshot = stale;
	snapshot.contacts  = buf.data();

	// The copy sees the buffer's epoch, not the stale struct's.
	EXPECT_EQ(contacts_epoch_current(snapshot.contacts), 100);
	EXPECT_EQ(contacts_epoch_begin(snapshot.contacts), 101);
	EXPECT_EQ(contacts_epoch_current(fa.contacts), 101);
}

// ---------------------------------------------------------------------------
// 4. Wraparound reset zeroes the stamps it invalidates (invariant 3).
// ---------------------------------------------------------------------------
TEST(ContactsEpoch, WraparoundResetZeroesTheStamps)
{
	std::vector<int> buf(CONTACTS_BUFFER_SIZE, 0);

	// Drive the counter to the edge and leave stamps at the high-water mark.
	buf[CONTACTS_EPOCH_SLOT] = INT_MAX - 2;
	const int last = contacts_epoch_begin(buf.data());
	ASSERT_EQ(last, INT_MAX - 1);
	for (int a = 0; a < 64; ++a) buf[a] = last;

	// Next eval must wrap. The reset has to take the stamps with it.
	const int epoch = contacts_epoch_begin(buf.data());
	EXPECT_EQ(epoch, 1);
	for (int a = 0; a < 64; ++a) {
		ASSERT_EQ(buf[a], 0) << "stale stamp survived the wraparound reset, atom=" << a;
		ASSERT_NE(buf[a], epoch);
	}
}

// ---------------------------------------------------------------------------
// 5. Independent buffers keep independent counters.
// ---------------------------------------------------------------------------
// Each worker thread owns its own contacts buffer; one worker's progress must
// not be visible to another.
TEST(ContactsEpoch, IndependentBuffersDoNotInterfere)
{
	std::vector<int> a(CONTACTS_BUFFER_SIZE, 0);
	std::vector<int> b(CONTACTS_BUFFER_SIZE, 0);

	for (int i = 0; i < 10; ++i) contacts_epoch_begin(a.data());
	for (int i = 0; i < 3;  ++i) contacts_epoch_begin(b.data());

	EXPECT_EQ(contacts_epoch_current(a.data()), 10);
	EXPECT_EQ(contacts_epoch_current(b.data()), 3);

	const int ea = contacts_epoch_current(a.data());
	a[7] = ea;
	EXPECT_EQ(b[7], 0);   // not visited in b
	EXPECT_NE(b[7], contacts_epoch_current(b.data()));
}

// ---------------------------------------------------------------------------
// 6. Allocation contract.
// ---------------------------------------------------------------------------
// The epoch slot must sit immediately past the addressable atom range, so it
// can never be clobbered by a legitimate atom stamp. Every allocator of a
// contacts buffer must size it CONTACTS_BUFFER_SIZE.
TEST(ContactsEpoch, EpochSlotSitsPastTheAtomRange)
{
	EXPECT_EQ(CONTACTS_EPOCH_SLOT, MAX_ATOM_NUMBER);
	EXPECT_EQ(CONTACTS_BUFFER_SIZE, MAX_ATOM_NUMBER + 1);

	// A freshly allocated buffer has never been stamped.
	std::vector<int> buf(CONTACTS_BUFFER_SIZE, 0);
	EXPECT_EQ(contacts_epoch_current(buf.data()), 0);

	// The first epoch is 1, so a zeroed slot can never read as visited.
	EXPECT_EQ(contacts_epoch_begin(buf.data()), 1);
	for (int a = 0; a < 32; ++a) EXPECT_NE(buf[a], contacts_epoch_current(buf.data()));
}

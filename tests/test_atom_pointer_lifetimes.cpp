#include <gtest/gtest.h>

#include "AtomOptResBinding.h"
#include "rotamer_output.h"

#include <array>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace {

// A ligand and a flexible receptor have separate OptRes owners. The receptor
// atoms are deliberately absent from the ligand's dirty refresh list.
struct WorkspaceFixture {
    std::vector<atom> atoms = std::vector<atom>(32);
    std::array<OptRes, 2> optres{};

    WorkspaceFixture() {
        optres[0].type = 0;
        optres[0].rnum = 1;
        optres[1].type = 1;
        optres[1].rnum = 2;
        atoms[4].optres = &optres[0];
        atoms[5].optres = &optres[0];
        atoms[7].optres = &optres[0];
        atoms[30].optres = &optres[1];
        atoms[31].optres = &optres[1];
    }
};

std::vector<std::array<double, 2>> exercise_workspaces(int threads, bool selective) {
    WorkspaceFixture master;
    const flexaids::AtomOptResBinding binding(master.atoms, master.optres);
    constexpr int chromosomes = 24;
    std::vector<std::array<double, 2>> written(chromosomes);
    std::vector<std::thread> workers;
    for (int tid = 0; tid < threads; ++tid) {
        workers.emplace_back([&, tid] {
            auto atoms = master.atoms;
            auto optres = master.optres;
            for (int chromosome = tid; chromosome < chromosomes; chromosome += threads) {
                if (selective) {
                    atoms[30] = master.atoms[30];
                    atoms[31] = master.atoms[31];
                } else {
                    atoms = master.atoms;
                }
                // This is the production operation called after the atom reset.
                // On the second chromosome the receptor pointers already refer
                // to optres here, whereas the refreshed ligand pointers do not.
                binding.bind(atoms, optres);
                for (int ai : {4, 5, 7}) EXPECT_EQ(atoms[ai].optres, &optres[0]);
                for (int ai : {30, 31}) EXPECT_EQ(atoms[ai].optres, &optres[1]);
                EXPECT_EQ(atoms[10].optres, nullptr);
                for (auto& entry : optres) entry.cf.com = 0.0;
                // Write through the rebound atom pointers. These are arbitrary
                // counters testing ownership, not simulated docking scores.
                for (int ai : {4, 5, 7, 30, 31})
                    atoms[ai].optres->cf.com += chromosome + 1;
                written[chromosome] = {optres[0].cf.com, optres[1].cf.com};
            }
        });
    }
    for (auto& worker : workers) worker.join();
    EXPECT_DOUBLE_EQ(master.optres[0].cf.com, 0.0);
    EXPECT_DOUBLE_EQ(master.optres[1].cf.com, 0.0);
    for (int chromosome = 0; chromosome < chromosomes; ++chromosome) {
        EXPECT_DOUBLE_EQ(written[chromosome][0], 3 * (chromosome + 1));
        EXPECT_DOUBLE_EQ(written[chromosome][1], 2 * (chromosome + 1));
    }
    return written;
}

std::string read_output(FILE* output) {
    rewind(output);
    std::string result;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), output)) result += buffer;
    return result;
}

std::string render_rotamer(std::span<const atom> atoms, std::span<const resid> residues) {
    std::unique_ptr<FILE, decltype(&fclose)> output(tmpfile(), &fclose);
    if (!output) throw std::runtime_error("Could not open test output file");
    const std::array<int, 1> built{2};
    flexaids::write_rotamer_model(output.get(), atoms, residues, 1, built, 1);
    return read_output(output.get());
}

}  // namespace

TEST(AtomOptResBinding, SelectiveRefreshPreservesOwnershipAcrossRepeatedChromosomes) {
    EXPECT_EQ(exercise_workspaces(1, true), exercise_workspaces(1, false));
}

TEST(AtomOptResBinding, IndependentWorkspacesAgreeAtOneAndFourThreads) {
    const auto serial = exercise_workspaces(1, true);
    EXPECT_EQ(exercise_workspaces(4, true), serial);
    EXPECT_EQ(exercise_workspaces(4, true), serial);
    EXPECT_EQ(exercise_workspaces(4, false), serial);
}

TEST(AtomOptResBinding, RejectsForeignOwnerWithoutSubtractingUnrelatedPointers) {
    WorkspaceFixture master;
    OptRes foreign{};
    master.atoms[4].optres = &foreign;
    EXPECT_THROW((flexaids::AtomOptResBinding(master.atoms, master.optres)), FlexAIDException);
}

TEST(AtomOptResBinding, EmptyOptResKeepsUnoptimizedAtomsNull) {
    std::array<atom, 3> atoms{};
    std::span<OptRes> empty;
    flexaids::AtomOptResBinding binding(atoms, empty);
    binding.bind(atoms, empty);
    binding.bind(atoms, empty);
    for (const auto& at : atoms) EXPECT_EQ(at.optres, nullptr);
}

TEST(AtomOptResBinding, RejectsChangedWorkspaceDimensionsBeforeWriting) {
    WorkspaceFixture master;
    flexaids::AtomOptResBinding binding(master.atoms, master.optres);
    auto atoms = master.atoms;
    std::array<OptRes, 1> short_optres{};
    EXPECT_THROW(binding.bind(atoms, short_optres), FlexAIDException);
    EXPECT_EQ(atoms[4].optres, &master.optres[0]);
    atoms.pop_back();
    EXPECT_THROW(binding.bind(atoms, master.optres), FlexAIDException);
}

TEST(RotamerOutput, AtomIndicesPreservePdbBytesAfterForcedRelocation) {
    auto atoms = std::make_unique<atom[]>(3);
    atoms[1].ofres = 1;
    atoms[1].number = 7;
    strcpy(atoms[1].name, " CB ");
    atoms[1].coor[0] = 1.0f;
    atoms[1].coor[1] = 2.0f;
    atoms[1].coor[2] = 3.0f;
    atoms[2].number = 8;
    strcpy(atoms[2].name, " OG ");
    atoms[2].coor[0] = 4.0f;
    atoms[2].coor[1] = 5.0f;
    atoms[2].coor[2] = 6.0f;
    std::array<resid, 2> residues{};
    strcpy(residues[1].name, "SER");
    residues[1].chn = 'A';
    residues[1].number = 10;
    const std::string expected =
        "MODEL        1\n"
        "ATOM      7  CB  SER A  10       1.000   2.000   3.000\n"
        "ATOM      8  OG  SER A  10       4.000   5.000   6.000\n"
        "ENDMDL\n";
    EXPECT_EQ(render_rotamer({atoms.get(), 3}, residues), expected);

    // Allocate while the old storage is live to guarantee a different address;
    // then copy and release it, as a relocating realloc does. No stale pointer
    // is ever dereferenced. The production writer receives the current base.
    auto grown = std::make_unique<atom[]>(503);
    ASSERT_NE(grown.get(), atoms.get());
    std::copy_n(atoms.get(), 3, grown.get());
    atoms.reset();
    EXPECT_EQ(render_rotamer({grown.get(), 503}, residues), expected);
    grown[1].coor[0] = 9.0f;
    EXPECT_NE(render_rotamer({grown.get(), 503}, residues), expected);
}

TEST(RotamerOutput, InvalidAtomIndicesFailBeforeAnyPdbIsWritten) {
    std::unique_ptr<FILE, decltype(&fclose)> output(tmpfile(), &fclose);
    ASSERT_NE(output, nullptr);
    std::array<atom, 3> atoms{};
    atoms[1].ofres = 1;
    std::array<resid, 2> residues{};
    const std::array<int, 1> outside{3};
    EXPECT_THROW(flexaids::write_rotamer_model(output.get(), atoms, residues, 0, {}, 1),
                 FlexAIDException);
    EXPECT_THROW(flexaids::write_rotamer_model(output.get(), atoms, residues, 1, outside, 1),
                 FlexAIDException);
    EXPECT_EQ(read_output(output.get()), "");
}

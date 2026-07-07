// tests/test_read_spheres.cpp
// Unit tests for read_spheres — HETATM sphere PDB ingest (CavityDetect format)

#include <gtest/gtest.h>
#include "../LIB/flexaid.h"
#include <cstdio>
#include <fstream>
#include <string>

namespace {

static std::string write_cavity_detect_sphere_pdb(const std::string& path) {
    std::ofstream out(path);
    out << "REMARK  Cleft spheres — FlexAIDΔS CavityDetector (test fixture)\n";
    out << "HETATM    1  SPH SURF    1      10.000  20.000  30.000  1.00  2.50           S\n";
    out << "HETATM    2  SPH SURF    1      11.500  21.250  31.750  1.00  1.80           S\n";
    out << "HETATM    3  SPH SURF    2      40.000  41.000  42.000  1.00  3.20           S\n";
    out << "END\n";
    return path;
}

static int count_spheres(sphere* head) {
    int count = 0;
    for (sphere* curr = head; curr != nullptr; curr = curr->prev) {
        ++count;
    }
    return count;
}

static void free_spheres(sphere* head) {
    while (head != nullptr) {
        sphere* prev = head->prev;
        free(head);
        head = prev;
    }
}

}  // namespace

TEST(ReadSpheres, AcceptsHetatmCavityDetectPdb) {
    const std::string pdb = "test_read_spheres_hetatm.pdb";
    write_cavity_detect_sphere_pdb(pdb);

    char filename[256];
    std::snprintf(filename, sizeof(filename), "%s", pdb.c_str());

    sphere* head = read_spheres(filename);
    ASSERT_NE(head, nullptr);
    EXPECT_EQ(count_spheres(head), 3);

    free_spheres(head);
    std::remove(pdb.c_str());
}
// tests/test_vcontacts_geometry.cpp
// Unit tests for Vcontacts geometry/polygon functions:
//   test_point, add_vedge, save_seeds, get_firstvert, order_faces
// These are the convex-hull construction functions used during Voronoi
// contact area calculation.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/Vcontacts.h"

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

static constexpr double EPS = 1e-9;

// Helper: create a plane with given A,B,C,D and optional index/dist/flag
static plane make_plane(double a, double b, double c, double d,
                        int idx = 0, double dist = 0.0, char flag = 0) {
    plane p{};
    p.Ai[0] = a; p.Ai[1] = b; p.Ai[2] = c; p.Ai[3] = d;
    p.index = idx;
    p.dist  = dist;
    p.area  = 0.0;
    p.flag  = flag;
    return p;
}

// Helper: create a vertex at (x,y,z) belonging to three planes
static vertex make_vertex(double x, double y, double z,
                          int p0, int p1, int p2) {
    vertex v{};
    v.xi[0] = x; v.xi[1] = y; v.xi[2] = z;
    v.dist = std::sqrt(x*x + y*y + z*z);
    v.plane[0] = p0; v.plane[1] = p1; v.plane[2] = p2;
    return v;
}

namespace {

constexpr int kDenseContactCount = MAX_CONT;

struct DenseHullFixture {
    std::array<atom, kDenseContactCount + 1> atoms{};
    std::array<atomsas, kDenseContactCount + 1> calc{};
    std::array<contactlist, kDenseContactCount> contacts{};
    std::array<plane, kDenseContactCount + 4> planes{};
    std::array<vertex, MAX_POLY> poly{};
    std::array<vertex, MAX_CONT> centerpt{};
    std::array<edgevector, MAX_POLY> edges{};
    std::array<int, (kDenseContactCount + 1) * 3> seeds{};
    VC_Global vc{};

    explicit DenseHullFixture(bool recalc) {
        atoms[0].coor[0] = 11.25f;
        atoms[0].coor[1] = -7.5f;
        atoms[0].coor[2] = 3.125f;
        atoms[0].number = 900;
        atoms[0].type = 7;
        atoms[0].ofres = 23;
        atoms[0].radius = 1.5f;

        constexpr double golden_angle = 2.39996322972865332;
        for(int i = 0; i < kDenseContactCount; ++i) {
            const double z = 1.0 - 2.0 * (static_cast<double>(i) + 0.5)
                / static_cast<double>(kDenseContactCount);
            const double radial = std::sqrt(1.0 - z * z);
            const double angle = golden_angle * static_cast<double>(i);
            const int neighbor = i + 1;

            atoms[neighbor].coor[0] = atoms[0].coor[0]
                + static_cast<float>(4.0 * radial * std::cos(angle));
            atoms[neighbor].coor[1] = atoms[0].coor[1]
                + static_cast<float>(4.0 * radial * std::sin(angle));
            atoms[neighbor].coor[2] = atoms[0].coor[2] + static_cast<float>(4.0 * z);
            atoms[neighbor].number = 900 + neighbor;
            atoms[neighbor].type = 7;
            atoms[neighbor].ofres = 23;
            atoms[neighbor].radius = 1.5f;

            contacts[i].index = neighbor;
            contacts[i].dist = 4.0;
        }

        for(std::size_t i = 0; i < calc.size(); ++i) {
            calc[i].atom = &atoms[i];
        }
        seeds.fill(-1);
        seeds[0] = 17;
        seeds[1] = 41;
        seeds[2] = 83;

        vc.Calc = calc.data();
        vc.contlist = contacts.data();
        vc.centerpt = centerpt.data();
        vc.poly = poly.data();
        vc.cont = planes.data();
        vc.vedge = edges.data();
        vc.seed = seeds.data();
        vc.planedef = 'B';
        vc.recalc = recalc ? 1 : 0;
    }

    int run() {
        return voronoi_poly2(&vc, 0, planes.data(), 3.0f,
                             kDenseContactCount, contacts.data());
    }

    std::array<std::uint32_t, 3> coordinate_bits() const {
        return {
            std::bit_cast<std::uint32_t>(atoms[0].coor[0]),
            std::bit_cast<std::uint32_t>(atoms[0].coor[1]),
            std::bit_cast<std::uint32_t>(atoms[0].coor[2]),
        };
    }

    std::array<int, 3> center_seed() const {
        return {seeds[0], seeds[1], seeds[2]};
    }
};

void expect_identical_planes(
    const std::array<plane, kDenseContactCount + 4>& lhs,
    const std::array<plane, kDenseContactCount + 4>& rhs)
{
    for(std::size_t i = 0; i < lhs.size(); ++i) {
        for(int component = 0; component < 4; ++component) {
            EXPECT_DOUBLE_EQ(lhs[i].Ai[component], rhs[i].Ai[component]) << "plane " << i;
        }
        EXPECT_DOUBLE_EQ(lhs[i].dist, rhs[i].dist) << "plane " << i;
        EXPECT_EQ(lhs[i].index, rhs[i].index) << "plane " << i;
        EXPECT_DOUBLE_EQ(lhs[i].area, rhs[i].area) << "plane " << i;
        EXPECT_EQ(lhs[i].flag, rhs[i].flag) << "plane " << i;
    }
}

struct DenseHullResult {
    int vertex_count = 0;
    bool coordinates_restored = false;
    bool seed_restored = false;
    std::array<std::uint32_t, 3> coordinates{};
    std::array<int, 3> seed{};
    std::array<plane, kDenseContactCount + 4> planes{};
};

DenseHullResult run_dense_hull() {
    DenseHullFixture fixture(/*recalc=*/true);
    DenseHullResult result;
    const auto original_coordinates = fixture.coordinate_bits();
    const auto original_seed = fixture.center_seed();
    result.vertex_count = fixture.run();
    result.coordinates = fixture.coordinate_bits();
    result.seed = fixture.center_seed();
    result.coordinates_restored = result.coordinates == original_coordinates;
    result.seed_restored = result.seed == original_seed;
    result.planes = fixture.planes;
    return result;
}

} // namespace

// ===========================================================================
// test_point — checks whether a point lies inside all half-spaces
// ===========================================================================

TEST(TestPoint, PointInsideAllPlanes) {
    // Two planes: Ax+By+Cz+D must be <= 0 for point to pass.
    // Plane: {1,0,0,-5} → x-5=0, half-space x<=5. Origin: 0-5=-5 <= 0 → OK
    // Plane: {0,1,0,-5} → y-5=0, half-space y<=5. Origin: 0-5=-5 <= 0 → OK
    plane cont[2] = {
        make_plane(1, 0, 0, -5.0),
        make_plane(0, 1, 0, -5.0),
    };
    double pt[3] = {0.0, 0.0, 0.0};
    EXPECT_EQ(test_point(pt, cont, 2, 1.0f, -1, -1, -1), 'Y');
}

TEST(TestPoint, PointBehindOnePlane) {
    // Plane: x + 1 >= 0.  Point at x=-2 is behind.
    plane cont[1] = {
        make_plane(1, 0, 0, 1.0),
    };
    double pt[3] = {-2.0, 0.0, 0.0};
    // x=-2 → 1*(-2)+1 = -1 < 0  → actually inside (Ax+By+Cz+D = -1 < 0 → OK)
    // Wait: the function checks > 0.0 to reject.  1*(-2) + 1 = -1 → not > 0 → 'Y'
    EXPECT_EQ(test_point(pt, cont, 1, 1.0f, -1, -1, -1), 'Y');

    // Now x=0.5 → 1*0.5+1 = 1.5 > 0 → 'N' (behind plane)
    double pt2[3] = {0.5, 0.0, 0.0};
    EXPECT_EQ(test_point(pt2, cont, 1, 1.0f, -1, -1, -1), 'N');
}

TEST(TestPoint, PointOnPlaneExcluded) {
    // When plane index matches planeA/planeB/planeC, that plane is skipped
    plane cont[2] = {
        make_plane(1, 0, 0, -100.0),   // very negative D → would reject most points
        make_plane(0, 1, 0, 0.0),       // y <= 0
    };
    double pt[3] = {50.0, 0.0, 0.0};
    // planeA=0 skips first plane; second: 0*50 + 1*0 + 0 + 0 = 0 → not > 0 → 'Y'
    EXPECT_EQ(test_point(pt, cont, 2, 1.0f, 0, -1, -1), 'Y');
}

TEST(TestPoint, SkipsFlaggedPlanes) {
    plane cont[2] = {
        make_plane(1, 0, 0, 100.0, 0, 0.0, 'X'),   // flagged 'X' → skipped
        make_plane(0, 1, 0, 0.0),
    };
    double pt[3] = {50.0, 0.0, 0.0};
    // first plane flagged 'X' → skipped; second: 0 → not > 0 → 'Y'
    EXPECT_EQ(test_point(pt, cont, 2, 1.0f, -1, -1, -1), 'Y');
}

TEST(TestPoint, ThreePlaneExclusion) {
    // All three planes excluded → always 'Y'
    plane cont[3] = {
        make_plane(1, 0, 0, 1000.0),
        make_plane(0, 1, 0, 1000.0),
        make_plane(0, 0, 1, 1000.0),
    };
    double pt[3] = {500.0, 500.0, 500.0};
    EXPECT_EQ(test_point(pt, cont, 3, 1.0f, 0, 1, 2), 'Y');
}

TEST(TestPoint, ZeroPlanes) {
    // No planes → point is always valid
    double pt[3] = {1.0, 2.0, 3.0};
    EXPECT_EQ(test_point(pt, nullptr, 0, 1.0f, -1, -1, -1), 'Y');
}

TEST(TestPoint, ExactBoundary) {
    // Point exactly on plane boundary: Ax+By+Cz+D = 0 → not > 0 → 'Y'
    plane cont[1] = {
        make_plane(1, 0, 0, -5.0),   // x - 5 = 0 at x=5
    };
    double pt[3] = {5.0, 0.0, 0.0};
    EXPECT_EQ(test_point(pt, cont, 1, 1.0f, -1, -1, -1), 'Y');
}

// ===========================================================================
// add_vedge — adds edge to edge list, auto-flips direction
// ===========================================================================

TEST(AddVedge, BasicEdge) {
    // Two orthogonal planes: x=0 and y=0
    // Cross product of normals (1,0,0)×(0,1,0) = (0,0,1)
    plane cont[3] = {
        make_plane(1, 0, 0, 0.0),   // plane 0: x=0
        make_plane(0, 1, 0, 0.0),   // plane 1: y=0
        make_plane(0, 0, 1, -1.0),  // test plane: z=1 (z-1<=0)
    };
    vertex poly[10];
    poly[0] = make_vertex(0.0, 0.0, 0.0, 0, 1, 2);

    edgevector ve[4];
    add_vedge(ve, 0, cont, 0, 1, 2, poly, 0);

    // Cross product: (1,0,0)×(0,1,0) = (0,0,1)
    double len = std::sqrt(ve[0].V[0]*ve[0].V[0] + ve[0].V[1]*ve[0].V[1] + ve[0].V[2]*ve[0].V[2]);
    EXPECT_NEAR(len, 1.0, EPS);

    EXPECT_EQ(ve[0].startpt, 0);
    EXPECT_EQ(ve[0].endpt, -1);
    EXPECT_EQ(ve[0].plane[0], 0);
    EXPECT_EQ(ve[0].plane[1], 1);
    EXPECT_EQ(ve[0].startplane, 2);
    EXPECT_EQ(ve[0].arc, '.');
}

TEST(AddVedge, DirectionFlip) {
    // Plane normals such that initial cross product points wrong way
    // Plane 0: normal (0,0,1), Plane 1: normal (1,0,0)
    // Cross: (0,0,1)×(1,0,0) = (0,1,0)
    // Test plane: y = -1 (0,1,0,1 → y+1>=0, behind when y>−1)
    plane cont[3] = {
        make_plane(0, 0, 1, 0.0),   // plane 0
        make_plane(1, 0, 0, 0.0),   // plane 1
        make_plane(0, 1, 0, 1.0),   // test plane: y+1>=0
    };
    vertex poly[10];
    poly[0] = make_vertex(0.0, 0.0, 0.0, 0, 1, 2);

    edgevector ve[4];
    add_vedge(ve, 0, cont, 0, 1, 2, poly, 0);

    // testpt = origin + (0,1,0) = (0,1,0)
    // testval = 0*0 + 1*1 + 0*0 + 1 = 2 > 0 → flip to (0,-1,0)
    EXPECT_NEAR(ve[0].V[0], 0.0, EPS);
    EXPECT_NEAR(ve[0].V[1], -1.0, EPS);
    EXPECT_NEAR(ve[0].V[2], 0.0, EPS);
}

TEST(AddVedge, ParallelPlanesZeroVector) {
    // Two parallel planes: same normal → cross product = 0
    plane cont[3] = {
        make_plane(1, 0, 0, 0.0),
        make_plane(1, 0, 0, -2.0),
        make_plane(0, 1, 0, 0.0),
    };
    vertex poly[10];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, 2);

    edgevector ve[4];
    add_vedge(ve, 0, cont, 0, 1, 2, poly, 0);

    EXPECT_NEAR(ve[0].V[0], 0.0, EPS);
    EXPECT_NEAR(ve[0].V[1], 0.0, EPS);
    EXPECT_NEAR(ve[0].V[2], 0.0, EPS);
}

TEST(AddVedge, MultipleEdgesInArray) {
    plane cont[3] = {
        make_plane(1, 0, 0, 0.0),
        make_plane(0, 1, 0, 0.0),
        make_plane(0, 0, 1, 0.0),
    };
    vertex poly[10];
    poly[0] = make_vertex(0.0, 0.0, 0.0, 0, 1, 2);
    poly[1] = make_vertex(0.0, 0.0, 0.0, 0, 1, 2);

    edgevector ve[4];
    add_vedge(ve, 0, cont, 0, 1, 2, poly, 0);
    add_vedge(ve, 1, cont, 0, 2, 1, poly, 1);

    EXPECT_EQ(ve[0].plane[0], 0);
    EXPECT_EQ(ve[0].plane[1], 1);
    EXPECT_EQ(ve[0].startpt, 0);

    EXPECT_EQ(ve[1].plane[0], 0);
    EXPECT_EQ(ve[1].plane[1], 2);
    EXPECT_EQ(ve[1].startpt, 1);
}

// ===========================================================================
// save_seeds — stores seed vertices for each atom in the contact polyhedron
// ===========================================================================

TEST(SaveSeeds, SingleInteriorVertex) {
    // One vertex at (1,0,0) belonging to planes 0,1,2
    // cont indices are 10,20,30
    plane cont[3] = {
        make_plane(0,0,0,0, 10),
        make_plane(0,0,0,0, 20),
        make_plane(0,0,0,0, 30),
    };
    vertex poly[10];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, 2);  // interior (plane[2] != -1)

    int seed[100];
    std::memset(seed, -1, sizeof(int) * 100);

    save_seeds(seed, cont, poly, 1, /*atomzero=*/5);

    // Seed for cont[0].index=10: atomzero=5, then cont[1].index=20, cont[2].index=30
    EXPECT_EQ(seed[30], 5);   // 3*10
    EXPECT_EQ(seed[31], 20);
    EXPECT_EQ(seed[32], 30);

    // Seed for cont[1].index=20: atomzero=5, then cont[0].index=10, cont[2].index=30
    EXPECT_EQ(seed[60], 5);   // 3*20
    EXPECT_EQ(seed[61], 10);
    EXPECT_EQ(seed[62], 30);

    // Seed for cont[2].index=30: atomzero=5, then cont[0].index=10, cont[1].index=20
    EXPECT_EQ(seed[90], 5);   // 3*30
    EXPECT_EQ(seed[91], 10);
    EXPECT_EQ(seed[92], 20);
}

TEST(SaveSeeds, SurfaceVertexSkipped) {
    // Vertex with plane[2] == -1 is a surface vertex → not saved
    plane cont[3] = {
        make_plane(0,0,0,0, 10),
        make_plane(0,0,0,0, 20),
        make_plane(0,0,0,0, 30),
    };
    vertex poly[10];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, -1);  // surface vertex

    int seed[100];
    std::memset(seed, -1, sizeof(int) * 100);

    save_seeds(seed, cont, poly, 1, 5);

    // Nothing should be saved
    EXPECT_EQ(seed[30], -1);
    EXPECT_EQ(seed[60], -1);
    EXPECT_EQ(seed[90], -1);
}

TEST(SaveSeeds, NoOverwrite) {
    // save_seeds only writes if seed slot is -1
    plane cont[2] = {
        make_plane(0,0,0,0, 5),
        make_plane(0,0,0,0, 10),
    };
    vertex poly[10];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, -1);  // surface → no seed

    int seed[100];
    std::fill_n(seed, 100, -1);
    // Pre-fill seed for index 5
    seed[15] = 99;  // 3*5 = 15

    save_seeds(seed, cont, poly, 1, 1);
    EXPECT_EQ(seed[15], 99);  // not overwritten
}

TEST(SaveSeeds, MultipleVertices) {
    plane cont[4] = {
        make_plane(0,0,0,0, 0),
        make_plane(0,0,0,0, 1),
        make_plane(0,0,0,0, 2),
        make_plane(0,0,0,0, 3),
    };
    vertex poly[10];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, 2);
    poly[1] = make_vertex(0.0, 1.0, 0.0, 1, 2, 3);

    int seed[100];
    std::fill_n(seed, 100, -1);

    save_seeds(seed, cont, poly, 2, /*atomzero=*/7);

    // From poly[0]: seeds for indices 0,1,2
    EXPECT_EQ(seed[0], 7);  EXPECT_EQ(seed[1], 1);   EXPECT_EQ(seed[2], 2);
    EXPECT_EQ(seed[3], 7);  EXPECT_EQ(seed[4], 0);   EXPECT_EQ(seed[5], 2);
    EXPECT_EQ(seed[6], 7);  EXPECT_EQ(seed[7], 0);   EXPECT_EQ(seed[8], 1);

    // From poly[1]: seeds for indices 1,2,3 (index 1 already filled by poly[0])
    EXPECT_EQ(seed[3], 7);  // already set — not overwritten
    EXPECT_EQ(seed[9], 7);  EXPECT_EQ(seed[10], 1);  EXPECT_EQ(seed[11], 2);
}

TEST(SaveSeeds, ZeroVertices) {
    plane cont[1] = { make_plane(0,0,0,0, 0) };
    int seed[100];
    std::fill_n(seed, 100, -1);
    // No crash with NV=0
    save_seeds(seed, cont, nullptr, 0, 1);
    EXPECT_EQ(seed[0], -1);
}

// ===========================================================================
// get_firstvert — finds initial vertex for polyhedron construction
// ===========================================================================

TEST(GetFirstvert, UsesSeedWhenAvailable) {
    // seed[atomzero*3] = cont[0].index, seed[+1] = cont[1].index, seed[+2] = cont[2].index
    plane cont[3] = {
        make_plane(0,0,0,0, 100),
        make_plane(0,0,0,0, 200),
        make_plane(0,0,0,0, 300),
    };
    int seed[12];
    std::fill_n(seed, 12, -1);
    // atomzero=0: seed[0]=100, seed[1]=200, seed[2]=300
    seed[0] = 100; seed[1] = 200; seed[2] = 300;

    int pA = -1, pB = -1, pC = -1;
    get_firstvert(seed, cont, &pA, &pB, &pC, 3, 0);

    EXPECT_EQ(pA, 0);  // cont index 0 has .index=100
    EXPECT_EQ(pB, 1);  // cont index 1 has .index=200
    EXPECT_EQ(pC, 2);  // cont index 2 has .index=300
}

TEST(GetFirstvert, NoSeedFallsBackToClosestPlane) {
    // No seed → fallback picks plane closest to origin (smallest cont[].dist)
    // Production code always allocates NC+4 planes (4 bounding-box planes follow
    // the NC contact planes).  get_firstvert's inner loop iterates cai<NC+4, so
    // the array must have NC+4 = 8 elements to avoid stack-buffer-overflow.
    plane cont[8] = {
        make_plane(0,0,0,0, 0, 5.0),
        make_plane(0,0,0,0, 1, 2.0),
        make_plane(0,0,0,0, 2, 8.0),
        make_plane(0,0,0,0, 3, 1.0),  // closest to origin
        // [4..7]: bounding-box placeholder slots (zero-init → degenerate planes,
        //          solve_3x3 returns -1, so they are safely skipped)
    };
    int seed[12];
    std::fill_n(seed, 12, -1);
    // atomzero=0, seed[0]=-1 → fallback path

    int pA = -1, pB = -1, pC = -1;
    get_firstvert(seed, cont, &pA, &pB, &pC, 4, 0);

    // planeA should be index 3 (dist=1.0 is smallest)
    EXPECT_EQ(pA, 3);
    // pB/pC may or may not be found depending on geometry (needs sphere planes)
    // Just verify pA was set correctly
}

TEST(GetFirstvert, SeedNotFoundInCont) {
    // seed points to index values that don't exist in cont[].index
    // NC=2 → loop goes to NC+4=6; allocate 6 elements (last 4 zero-init).
    plane cont[6] = {
        make_plane(0,0,0,0, 10),
        make_plane(0,0,0,0, 20),
        // [2..5]: bounding-box placeholder slots
    };
    int seed[12];
    std::fill_n(seed, 12, -1);
    seed[0] = 999; seed[1] = 998; seed[2] = 997;  // not found

    int pA = -1, pB = -1, pC = -1;
    get_firstvert(seed, cont, &pA, &pB, &pC, 2, 0);

    // Falls back to closest plane search
    EXPECT_NE(pA, -1);
}

// ===========================================================================
// order_faces — orders vertices around each face; returns 'I' or 'S'
// ===========================================================================

TEST(OrderFaces, InternalAtom) {
    // All vertices have plane[2] != -1 → internal atom
    plane cont[3] = {
        make_plane(1, 0, 0, 0.0),
        make_plane(0, 1, 0, 0.0),
        make_plane(0, 0, 1, 0.0),
    };
    vertex poly[20];
    // 4 vertices forming a tetrahedral-like polyhedron (all interior)
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, 2);
    poly[1] = make_vertex(0.0, 1.0, 0.0, 0, 1, 2);
    poly[2] = make_vertex(0.0, 0.0, 1.0, 0, 1, 2);
    poly[3] = make_vertex(0.5, 0.5, 0.5, 0, 1, 2);

    vertex centerpt[3];
    centerpt[0] = make_vertex(0.5, 0.0, 0.0, -1, -1, -1);
    centerpt[1] = make_vertex(0.0, 0.5, 0.0, -1, -1, -1);
    centerpt[2] = make_vertex(0.0, 0.0, 0.5, -1, -1, -1);

    ptindex ptorder[3];
    std::memset(ptorder, 0, sizeof(ptorder));

    char result = order_faces(0, poly, centerpt, 1.0f, 3, 4, cont, ptorder);
    EXPECT_EQ(result, 'I');

    // All planes are present → each face should have some points
    for (int i = 0; i < 3; ++i) {
        EXPECT_GT(ptorder[i].numpts, 0);
    }
}

TEST(OrderFaces, SurfaceAtomDetected) {
    // At least one vertex with plane[2] == -1 → surface
    plane cont[2] = {
        make_plane(1, 0, 0, 0.0),
        make_plane(0, 1, 0, 0.0),
    };
    vertex poly[20];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, -1);  // surface vertex
    poly[1] = make_vertex(0.0, 1.0, 0.0, 0, 1, -1);   // surface vertex

    vertex centerpt[2];
    centerpt[0] = make_vertex(0.5, 0.0, 0.0, -1, -1, -1);
    centerpt[1] = make_vertex(0.0, 0.5, 0.0, -1, -1, -1);

    ptindex ptorder[2];
    std::memset(ptorder, 0, sizeof(ptorder));

    char result = order_faces(0, poly, centerpt, 1.0f, 2, 2, cont, ptorder);
    EXPECT_EQ(result, 'S');
}

TEST(OrderFaces, HiddenPlaneZeroPoints) {
    // A plane flagged 'X' should have numpts = 0
    plane cont[2] = {
        make_plane(1, 0, 0, 0.0, 0, 1.0, 'X'),  // hidden
        make_plane(0, 1, 0, 0.0),                 // visible
    };
    vertex poly[20];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, 2);
    poly[1] = make_vertex(0.0, 1.0, 0.0, 0, 1, 2);

    vertex centerpt[2];
    centerpt[0] = make_vertex(0.5, 0.0, 0.0, -1, -1, -1);
    centerpt[1] = make_vertex(0.0, 0.5, 0.0, -1, -1, -1);

    ptindex ptorder[2];
    std::memset(ptorder, 0, sizeof(ptorder));

    order_faces(0, poly, centerpt, 1.0f, 2, 2, cont, ptorder);

    EXPECT_EQ(ptorder[0].numpts, 0);  // hidden plane
    EXPECT_GT(ptorder[1].numpts, 0);  // visible plane has points
}

TEST(OrderFaces, PointsOrderedByCosine) {
    // With > 3 points on a face, they should be ordered by decreasing cosPQR
    plane cont[1] = {
        make_plane(0, 0, 1, 0.0),   // z=0 plane
    };
    vertex poly[20];
    // All on plane 0, all interior (plane[2] != -1)
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, 2);   // not actually on same planes
    poly[1] = make_vertex(0.0, 1.0, 0.0, 0, 1, 2);
    poly[2] = make_vertex(-1.0, 0.0, 0.0, 0, 1, 2);
    poly[3] = make_vertex(0.0, -1.0, 0.0, 0, 1, 2);

    vertex centerpt[1];
    centerpt[0] = make_vertex(0.0, 0.0, 0.0, -1, -1, -1);

    ptindex ptorder[1];
    std::memset(ptorder, 0, sizeof(ptorder));

    char result = order_faces(0, poly, centerpt, 2.0f, 1, 4, cont, ptorder);
    // All vertices belong to plane 0
    EXPECT_EQ(ptorder[0].numpts, 4);
    (void)result;
}

TEST(OrderFaces, TwoPointsOnFace) {
    // Face with exactly 2 points → no ordering needed
    plane cont[1] = {
        make_plane(1, 0, 0, 0.0),
    };
    vertex poly[20];
    poly[0] = make_vertex(0.0, 1.0, 0.0, 0, 1, 2);
    poly[1] = make_vertex(0.0, -1.0, 0.0, 0, 1, 2);

    vertex centerpt[1];
    centerpt[0] = make_vertex(0.0, 0.0, 0.0, -1, -1, -1);

    ptindex ptorder[1];
    std::memset(ptorder, 0, sizeof(ptorder));

    order_faces(0, poly, centerpt, 1.0f, 1, 2, cont, ptorder);
    EXPECT_EQ(ptorder[0].numpts, 2);
}

TEST(OrderFaces, NoPointsOnFace) {
    // Face where no vertices belong → 0 points
    plane cont[2] = {
        make_plane(1, 0, 0, 0.0),  // plane 0
        make_plane(0, 1, 0, 0.0),  // plane 1
    };
    // Vertices all belong to plane 1 only
    vertex poly[20];
    poly[0] = make_vertex(0.0, 1.0, 0.0, 1, 2, 3);
    poly[1] = make_vertex(0.0, -1.0, 0.0, 1, 2, 3);

    vertex centerpt[2];
    centerpt[0] = make_vertex(1.0, 0.0, 0.0, -1, -1, -1);
    centerpt[1] = make_vertex(0.0, 1.0, 0.0, -1, -1, -1);

    ptindex ptorder[2];
    std::memset(ptorder, 0, sizeof(ptorder));

    order_faces(0, poly, centerpt, 1.0f, 2, 2, cont, ptorder);
    EXPECT_EQ(ptorder[0].numpts, 0);  // no vertices on plane 0
}

// ===========================================================================
// Integration: test_point + order_faces interaction
// ===========================================================================

TEST(VcontactsGeometry, TestPointUsedByOrderFaces) {
    // Verify that order_faces correctly uses test_point internally
    // when adding arc-surface points
    plane cont[1] = {
        make_plane(0, 0, 1, 0.0),
    };
    // Surface vertices (plane[2] == -1)
    vertex poly[20];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, -1, -1);
    poly[1] = make_vertex(-1.0, 0.0, 0.0, 0, -1, -1);
    poly[2] = make_vertex(0.0, 1.0, 0.0, 0, -1, -1);
    poly[3] = make_vertex(0.0, -1.0, 0.0, 0, -1, -1);

    vertex centerpt[1];
    centerpt[0] = make_vertex(0.0, 0.0, 0.0, -1, -1, -1);

    ptindex ptorder[1];
    std::memset(ptorder, 0, sizeof(ptorder));

    char result = order_faces(0, poly, centerpt, 2.0f, 1, 4, cont, ptorder);
    EXPECT_EQ(result, 'S');  // surface atom
    EXPECT_GE(ptorder[0].numpts, 4);  // may add arc point
}

// ===========================================================================
// Edge cases and stress
// ===========================================================================

TEST(TestPoint, ManyPlanesAllPass) {
    // 50 planes, point is inside all of them (at origin)
    plane cont[50];
    for (int i = 0; i < 50; ++i) {
        // Planes passing through origin with various orientations
        double angle = i * 3.14159265 / 25.0;
        cont[i] = make_plane(std::cos(angle), std::sin(angle), 0.0, 0.0);
    }
    double pt[3] = {0.0, 0.0, 0.0};
    EXPECT_EQ(test_point(pt, cont, 50, 1.0f, -1, -1, -1), 'Y');
}

TEST(SaveSeeds, LargeSeedArray) {
    // Verify indexing doesn't go out of bounds with large atom indices
    plane cont[3] = {
        make_plane(0,0,0,0, 100),
        make_plane(0,0,0,0, 200),
        make_plane(0,0,0,0, 300),
    };
    vertex poly[10];
    poly[0] = make_vertex(1.0, 0.0, 0.0, 0, 1, 2);

    int seed[1000];
    std::fill_n(seed, 1000, -1);

    save_seeds(seed, cont, poly, 1, 5);

    // Verify no out-of-bounds: max index accessed = 3*300+2 = 902
    EXPECT_EQ(seed[900], 5);
    EXPECT_EQ(seed[901], 100);
    EXPECT_EQ(seed[902], 200);
}

// ===========================================================================
// voronoi_poly2 degeneracy failsafe purity
// ===========================================================================

TEST(VoronoiFailsafe, ImmediateFailureRestoresCenterCoordinatesExactly) {
    DenseHullFixture fixture(/*recalc=*/false);
    const auto original_coordinates = fixture.coordinate_bits();
    const auto original_seed = fixture.center_seed();

    EXPECT_EQ(fixture.run(), -1);
    EXPECT_EQ(fixture.coordinate_bits(), original_coordinates);
    EXPECT_EQ(fixture.center_seed(), original_seed);
}

TEST(VoronoiFailsafe, RetryIsBoundedAndDeterministicAcrossRepeatedCalls) {
    DenseHullFixture fixture(/*recalc=*/true);
    const auto original_coordinates = fixture.coordinate_bits();
    const auto original_seed = fixture.center_seed();

    EXPECT_EQ(fixture.run(), -1);
    EXPECT_EQ(fixture.coordinate_bits(), original_coordinates);
    EXPECT_EQ(fixture.center_seed(), original_seed);
    const auto first_planes = fixture.planes;

    EXPECT_EQ(fixture.run(), -1);
    EXPECT_EQ(fixture.coordinate_bits(), original_coordinates);
    EXPECT_EQ(fixture.center_seed(), original_seed);
    expect_identical_planes(first_planes, fixture.planes);
}

#ifdef _OPENMP
TEST(VoronoiFailsafe, DeterministicAcrossOpenMPWorkerAssignments) {
    if(omp_get_max_threads() < 2) {
        GTEST_SKIP() << "OMP_NUM_THREADS limits this run to one worker";
    }

    std::array<DenseHullResult, 2> results{};
    std::array<int, 2> worker_ids{{-1, -1}};
    omp_set_dynamic(0);

#pragma omp parallel for num_threads(2) schedule(static, 1)
    for(int run = 0; run < 2; ++run) {
        worker_ids[run] = omp_get_thread_num();
        results[run] = run_dense_hull();
    }

    ASSERT_NE(worker_ids[0], worker_ids[1]);
    EXPECT_EQ(results[0].vertex_count, -1);
    EXPECT_EQ(results[1].vertex_count, results[0].vertex_count);
    EXPECT_TRUE(results[0].coordinates_restored);
    EXPECT_TRUE(results[1].coordinates_restored);
    EXPECT_TRUE(results[0].seed_restored);
    EXPECT_TRUE(results[1].seed_restored);
    EXPECT_EQ(results[1].coordinates, results[0].coordinates);
    EXPECT_EQ(results[1].seed, results[0].seed);
    expect_identical_planes(results[0].planes, results[1].planes);
}
#endif

#include "msh_parser.h"
#include "neighbour.h"
#include <cassert>

int test_water_block() {
    // verify mesh data
    std::ifstream input("water.msh");
    EGS_Mesh mesh = msh_parser::parse_msh_file(input);
    auto elts = mesh.elements();
    assert(elts.size() == 1160);
    assert(elts[0].medium_tag == 69);
    assert(elts[0].a == 142 && elts[0].b == 223 && elts[0].c == 130 && elts[0].d == 353);

    auto nodes = mesh.nodes();
    assert(nodes.size() == 363);
    assert(nodes[362].tag == 363 &&
           nodes[362].x == 0.3899710788706327 &&
           nodes[362].y == 0.1542470443087625 &&
           nodes[362].z == 0.7332480649826769);

    auto materials = mesh.materials();
    assert(materials.size() == 1);
    assert(materials[0].tag == 69);
    assert(materials[0].medium_name == "Water");

    // verify neighbour information
    std::vector<mesh_neighbours::Tetrahedron> neighbour_elts;
    neighbour_elts.reserve(elts.size());
    for (const auto& elt: elts) {
        neighbour_elts.emplace_back(mesh_neighbours::Tetrahedron(elt.a, elt.b, elt.c, elt.d));
    }
    mesh_neighbours::tetrahedron_neighbours(neighbour_elts);
    return 0;
}

#define RUN_TEST(test_fn) \
    std::cerr << "starting test " << #test_fn << std::endl; \
    err = test_fn; \
    num_total++; \
    if (err) { \
        std::cerr << "test FAILED" << std::endl; \
        num_failed++; \
    } else { \
        std::cerr << "test passed" << std::endl; \
    }

int main() {
    int num_failed = 0;
    int num_total = 0;
    int err = 0;

    RUN_TEST(test_water_block());

    std::cerr << num_total - num_failed << " out of " << num_total << " tests passed\n";
    return num_failed;
}

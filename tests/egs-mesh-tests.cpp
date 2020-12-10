#include "msh_parser.h"
#include "neighbour.h"
#include <cassert>

// O(n2) neighbour finding function to verify our implementation
std::vector<std::array<std::size_t, 4>> naive_neighbours(const std::vector<mesh_neighbours::Tetrahedron>& elements) {
    using mesh_neighbours::NONE;
    std::vector<std::array<std::size_t, 4>> nbrs(elements.size(), {NONE, NONE, NONE, NONE});
    for (std::size_t i = 0; i < elements.size(); i++) {
        auto elt_faces = elements[i].faces();
        for (std::size_t f = 0; f < 4; f++) {
            if (nbrs[i][f] != NONE) {
                continue;
            }
            for (std::size_t j = 0; j < elements.size(); j++) {
                if (i == j) {
                    continue;
                }
                auto other_faces = elements[j].faces();
                for (std::size_t fj = 0; fj < 4; fj++) {
                    if (elt_faces[f] == other_faces[fj]) {
                        nbrs[i][f] = j;
                        nbrs[j][fj] = i;
                        break;
                    }
                }
            }
        }
    }
    return nbrs;
}

int test_water_block() {
    // verify mesh data
    std::ifstream input("water.msh");
    EGS_Mesh mesh = msh_parser::parse_msh_file(input);
    auto elts = mesh.elements();
    assert(elts.size() == 1160);
    assert(elts[0].medium_tag == 1);
    assert(elts[0].a == 142 && elts[0].b == 223 && elts[0].c == 130 && elts[0].d == 353);

    auto nodes = mesh.nodes();
    assert(nodes.size() == 363);
    assert(nodes[362].tag == 363 &&
           nodes[362].x == 0.3899710788706327 &&
           nodes[362].y == 0.1542470443087625 &&
           nodes[362].z == 0.7332480649826769);

    auto materials = mesh.materials();
    assert(materials.size() == 1);
    assert(materials[0].tag == 1);
    assert(materials[0].medium_name == "Water");

    // verify neighbour information
    std::vector<mesh_neighbours::Tetrahedron> neighbour_elts;
    neighbour_elts.reserve(elts.size());
    for (const auto& elt: elts) {
        neighbour_elts.emplace_back(mesh_neighbours::Tetrahedron(elt.a, elt.b, elt.c, elt.d));
    }
    auto nbrs = mesh_neighbours::tetrahedron_neighbours(neighbour_elts);
    std::cout << "element 1 has neighbours "
        << nbrs[0][0] + 1 << " " << nbrs[0][1] + 1 << " " << nbrs[0][2] + 1 << " " << nbrs[0][3] + 1 << "\n";
    assert(nbrs == naive_neighbours(neighbour_elts));

    // check that we have no isolated tetrahedrons
    for (std::size_t i = 0; i < elts.size(); i++) {
        assert(!(nbrs[i][0] == mesh_neighbours::NONE &&
                 nbrs[i][1] == mesh_neighbours::NONE &&
                 nbrs[i][2] == mesh_neighbours::NONE &&
                 nbrs[i][3] == mesh_neighbours::NONE));
    }
    return 0;
}

int test_water10000_block() {
    // verify mesh data
    std::ifstream input("water10000.msh");
    EGS_Mesh mesh = msh_parser::parse_msh_file(input);
    auto elts = mesh.elements();
    assert(elts.size() == 9280);
    assert(elts[0].medium_tag == 1);
    assert(elts[0].a == 142 && elts[0].b == 364 && elts[0].c == 366 && elts[0].d == 367);

    auto nodes = mesh.nodes();
    assert(nodes.size() == 2197);
    assert(nodes[2196].tag == 2197 &&
           nodes[2196].x == 0.8045166131834418 &&
           nodes[2196].y == 0.175446902578746 &&
           nodes[2196].z == 0.1343100687184781);

    auto materials = mesh.materials();
    assert(materials.size() == 1);
    assert(materials[0].tag == 1);
    assert(materials[0].medium_name == "Water");

    // verify neighbour information
    std::vector<mesh_neighbours::Tetrahedron> neighbour_elts;
    neighbour_elts.reserve(elts.size());
    for (const auto& elt: elts) {
        neighbour_elts.emplace_back(mesh_neighbours::Tetrahedron(elt.a, elt.b, elt.c, elt.d));
    }
    auto nbrs = mesh_neighbours::tetrahedron_neighbours(neighbour_elts);
    std::cout << "element 1 has neighbours "
        << nbrs[0][0] + 1 << " " << nbrs[0][1] + 1 << " " << nbrs[0][2] + 1 << " " << nbrs[0][3] + 1 << "\n";
    assert(nbrs == naive_neighbours(neighbour_elts));

    // check that we have no isolated tetrahedrons
    for (std::size_t i = 0; i < elts.size(); i++) {
        assert(!(nbrs[i][0] == mesh_neighbours::NONE &&
                 nbrs[i][1] == mesh_neighbours::NONE &&
                 nbrs[i][2] == mesh_neighbours::NONE &&
                 nbrs[i][3] == mesh_neighbours::NONE));
    }

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
    RUN_TEST(test_water10000_block());

    std::cerr << num_total - num_failed << " out of " << num_total << " tests passed\n";
    return num_failed;
}

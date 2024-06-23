// Computational Fabrication Assignment #1
//  By David Levin 2014
#include <cmath>
#include <iostream>
#include <vector>

#include "../include/Voxelizer.h"

int main(int argc, char **argv) {
    unsigned int dim = 32;      // Default dimension of voxel grid
    unsigned int num_rays = 2;  // Default number of rays for surface area estimation
    unsigned int distance = 3;  // Default distance to mesh surface when carving
    double jitter = 0.1;        // Default jittering factor for ray directions
    int num_iterations = 10;    // Default number of iterations for base expansion
    int base_distance = 1;      // Default balance distance from base to barycenter

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--dim" && i + 1 < argc) {
            dim = std::atoi(argv[++i]);
            if (dim <= 0) {
                std::cout << "Invalid dimension\n";
                return 0;
            }
        } else if (arg == "--num-rays" && i + 1 < argc) {
            num_rays = std::atoi(argv[++i]);
            if (num_rays <= 0) {
                std::cout << "Invalid number of rays\n";
                return 0;
            }
        } else if (arg == "--distance" && i + 1 < argc) {
            distance = std::atoi(argv[++i]);
        } else if (arg == "--jitter" && i + 1 < argc) {
            jitter = std::atof(argv[++i]);
        } else if (arg == "--base-distance" && i + 1 < argc) {
            base_distance = std::atoi(argv[++i]);
        } else if (arg == "--max-base-dilation" && i + 1 < argc) {
            num_iterations = std::atoi(argv[++i]);
        } else if (arg == "--help") {
            std::cout << "make_it_stand: A tool for modifying 3D models to ensure they can stand stably.\n"
                      << "This program adjusts models by carving the inside to reduce weight and dilating\n"
                      << "the base in a specific direction to increase the footprint, thereby enhancing stability.\n"
                      << "It utilizes a voxel grid representation for precise modifications.\n\n"
                      << "Usage: make_it_stand [OPTIONS]\n"
                      << "Options:\n"
                      << "  --dim <dimension>            Dimension of voxel grid (default: 32)\n"
                      << "  --num-rays <number>          Number of rays per voxel (default: 2)\n"
                      << "  --jitter <jitter>            Jittering factor for ray directions (default: 0.1)\n"
                      << "  --distance <distance>        Distance to mesh surface when carving (default: 3)\n"
                      << "  --base-distance <distance>   Distance from base to barycenter (default: 1)\n"
                      << "  --max-base-dilation <number> Number of iterations for base expansion (default: 10)\n"
                      << "  --help                       Display this help and exit\n";
            return 0;
        }
    }

    // print all settings
    std::cout << "Running with the following settings: \n"
              << "\tDimension of voxel grid: " << dim << "x" << dim << "x" << dim << "\n"
              << "\tNumber of rays per voxel: " << num_rays << "\n"
              << "\tJittering factor for ray directions: " << jitter << "\n"
              << "\tDistance to mesh surface when carving: " << distance << "\n"
              << "\tDistance from base to barycenter: " << base_distance << "\n"
              << "\tNumber of iterations for base expansion: " << num_iterations << std::endl;

    std::cout << "Load Mesh : " << argv[1] << "\n";

    // loadVoxels("../data/voxels.txt");
    loadMesh(argv[1], dim);
    voxelize(num_rays, jitter, distance);
    fill_holes();
    mark_inner_voxels(distance);
    saveInnerVoxelsToObj("../data/inner_before.obj");

    CompFab::VoxelGrid *base = compute_base();
    // fill base
    CompFab::VoxelGrid *base_filled = fill_base(base);
    CompFab::VoxelGrid *base_eroded = erode2d(base_filled, base_distance);
    saveVoxelsToObj("../data/base_eroded.obj", base_eroded);

    CompFab::Vec3 barycenter = std::get<0>(compute_barycenter());
    std::cout << "barycenter: " << barycenter.m_x << " " << barycenter.m_y << " " << barycenter.m_z << std::endl;

    std::cout << is_balanced(barycenter, base_eroded) << std::endl;
    CompFab::Vec3 dir = compute_direction_to_base(barycenter, base_eroded);
    std::cout << "direction: " << dir.m_x << " " << dir.m_y << " " << dir.m_z << std::endl;
    std::vector<CompFab::Vec3> sorted_voxels = sort_voxels(barycenter, dir);
    bool isbalanced = carve(base_eroded);
    if (!isbalanced) {
        expand_base(base_eroded, num_iterations);
    }

    save_barycenter(barycenter, "../data/barycenter.obj", base_eroded);

    barycenter = std::get<0>(compute_barycenter());
    save_barycenter(barycenter, "../data/barycenter_after.obj", base_eroded);

    saveVoxelsToObj(argv[2]);
    saveVoxelsToObj("../data/base.obj", base);
    saveVoxelsToObj("../data/base_filled.obj", base_filled);
    // saveVoxels("../data/teapot_voxels.txt");
    saveInnerVoxelsToObj("../data/inner.obj");
}

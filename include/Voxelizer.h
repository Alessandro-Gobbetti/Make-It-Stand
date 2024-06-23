#include <cmath>
#include <vector>

#include "../include/CompFab.h"
#include "../include/Mesh.h"

// Ray-Triangle Intersection
// Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle);

// Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

// Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir);

// random sampling of directions (for better results)
std::vector<CompFab::Vec3> generateDirections(unsigned int num_rays, double jitter);

void mark_inner_voxels(unsigned int voxel_distance);
void mark_inner_voxels(CompFab::VoxelGrid *voxels, unsigned int voxel_distance);

CompFab::VoxelGrid *voxelize(unsigned int num_rays, double jitter, unsigned int voxel_distance);

CompFab::VoxelGrid *compute_base();
CompFab::VoxelGrid *fill_base(CompFab::VoxelGrid *base);

std::tuple<CompFab::Vec3, CompFab::Vec3, unsigned int> compute_barycenter();
bool is_balanced(CompFab::Vec3 barycenter, CompFab::VoxelGrid *base);
void save_barycenter(CompFab::Vec3 barycenter, const char *outfile = "../data/barycenter.obj", CompFab::VoxelGrid *base = NULL);

CompFab::Vec3 compute_direction_to_base(CompFab::Vec3 barycenter, CompFab::VoxelGrid *base_filled);
std::vector<CompFab::Vec3> sort_voxels(CompFab::Vec3 barycenter, CompFab::Vec3 direction);
bool carve(CompFab::VoxelGrid *base);
CompFab::VoxelGrid *erode2d(CompFab::VoxelGrid *base, int iterations);
void expand_base(CompFab::VoxelGrid *base, int max_iterations);

void fill_holes();

bool loadMesh(char *filename, unsigned int dim);

void saveVoxelsToObj(const char *outfile);
void saveVoxelsToObj(const char *outfile, CompFab::VoxelGrid *voxels);
void saveInnerVoxelsToObj(const char *outfile);
void saveInnerVoxelsToObj(const char *outfile, CompFab::VoxelGrid *voxels);

void saveVoxels(const char *outfile, CompFab::VoxelGrid *voxels);
void saveVoxels(const char *outfile);

void loadVoxels(const char *filename);

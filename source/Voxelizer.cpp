#include "../include/Voxelizer.h"

#include <omp.h>

#include <iostream>

#include "../include/ConvexHull.h"

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;
// CompFab::VoxelGrid *g_voxelGridBase;
int base_layer = -1;

// #include "../include/CompFab.h"

// Ray-Triangle Intersection
// Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle) {
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle,
     * 0 otherwise */

    CompFab::Vec3 edge1 = triangle.m_v2 - triangle.m_v1;
    CompFab::Vec3 edge2 = triangle.m_v3 - triangle.m_v1;

    // compute determinant
    CompFab::Vec3 pvec = ray.m_direction % edge2;  // cross product
    float det = edge1 * pvec;                      // dot product

    // if determinant is near zero, ray lies in plane of triangle
    if (det < EPSILON && det > -EPSILON) {
        return 0;
    }

    // compute barycentric coordinates u and v and test bounds
    float invDet = 1.0f / det;
    CompFab::Vec3 svec = ray.m_origin - triangle.m_v1;

    float u = invDet * (svec * pvec);
    if (u < 0 || u > 1) {  // if u is out of bounds, no intersection
        return 0;
    }

    CompFab::Vec3 qvec = svec % edge1;
    float v = invDet * (ray.m_direction * qvec);
    if (v < 0 || u + v > 1) {  // if v is out of bounds, no intersection
        return 0;
    }

    // compute t, distance along ray to triangle
    float t = invDet * (edge2 * qvec);

    if (t > 0) {
        // intersection point is ray_origin + t * ray_direction
        return 1;
    }

    // this means that there is a line intersection but not a ray intersection (the ray is pointing away from the triangle)
    return 0;
}

// Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir) {
    unsigned int numHits = 0;

    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir,
     * from voxel center voxelPos intersects the surface */

    CompFab::Ray ray(voxelPos, dir);
    for (unsigned int i = 0; i < g_triangleList.size(); i++) {
        numHits += rayTriangleIntersection(ray, g_triangleList[i]);
    }

    return numHits;
}

// random sampling of directions (for better results)
std::vector<CompFab::Vec3> generateDirections(unsigned int num_rays, double jitter) {
    // fibonacci sphere algorithm
    std::vector<CompFab::Vec3> directions;

    const double goldenRatio = (1 + sqrt(5)) / 2.0;
    const double angleIncrement = M_PI * 2 * goldenRatio;

    for (unsigned int i = 0; i < num_rays; i++) {
        double t = (double)i / (double)num_rays;
        double inclination = acos(1 - 2 * t);
        double azimuth = angleIncrement * i;

        // add jitter to inclination and azimuth
        // note that we add twice more randomness to inclination than azimuth since the range of inclination is smaller
        inclination += 0.5 * jitter * ((double)rand() / RAND_MAX - (double)rand() / RAND_MAX);
        azimuth += jitter * ((double)rand() / RAND_MAX - (double)rand() / RAND_MAX);

        // convert from spherical to cartesian coordinates
        double x = sin(inclination) * cos(azimuth);
        double y = sin(inclination) * sin(azimuth);
        double z = cos(inclination);

        CompFab::Vec3 direction(x, y, z);
        direction.normalize();

        directions.push_back(direction);
    }

    return directions;
}

CompFab::VoxelGrid *voxelize(unsigned int num_rays, double jitter, unsigned int voxel_distance = 0) {
    // Cast ray, check if voxel is inside or outside
    // even number of surface intersections = outside (OUT then IN then OUT)
    //  odd number = inside (IN then OUT)
    CompFab::Vec3 voxelPos;
    // CompFab::Vec3 direction(1.0,0.0,0.0);

    // unsigned int num_rays = 10;

    double spacing = g_voxelGrid->m_spacing;

/********* ASSIGNMENT *********/
/* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
 * surface defined by the triangles in g_triangleList */

// print pragma number of cores
#pragma omp parallel
    {
#pragma omp master
        {
            std::cout << "Executing with " << omp_get_num_threads() << " threads\n";
        }
    }

    unsigned int iterations_count = 0;

// parallelize with openmp
#pragma omp parallel for private(voxelPos) schedule(dynamic)
    for (unsigned int i = 0; i < g_voxelGrid->m_dimX; i++) {
// percentage of completion
// atomically update the iterations count
#pragma omp atomic
        iterations_count++;

// print percentage of completion
#pragma omp critical
        std::cout << "\rVoxelizing: " << (iterations_count * 100) / g_voxelGrid->m_dimX << "%" << std::flush;

        for (unsigned int j = 0; j < g_voxelGrid->m_dimY; j++) {
            for (unsigned int k = 0; k < g_voxelGrid->m_dimZ; k++) {
                voxelPos = g_voxelGrid->m_lowerLeft +
                           CompFab::Vec3(i * spacing, j * spacing, k * spacing);

                std::vector<CompFab::Vec3> directions = generateDirections(num_rays, jitter);
                // std::vector<CompFab::Vec3> directions;
                // directions.push_back(CompFab::Vec3(0.0, 0.0, 1.0));

                unsigned int count = 0;
                for (unsigned int l = 0; l < directions.size(); l++) {
                    int numIntersections = numSurfaceIntersections(voxelPos, directions[l]);
                    count += numIntersections % 2;
                }

                g_voxelGrid->isInside(i, j, k) = count > directions.size() / 2.0;
            }
        }
    }
    std::cout << std::endl;

    // if (voxel_distance > 0) {
    //     mark_inner_voxels(g_voxelGrid, voxel_distance);
    // }

    return g_voxelGrid;
}

void mark_inner_voxels(unsigned int voxel_distance) {
    mark_inner_voxels(g_voxelGrid, voxel_distance);
}

void mark_inner_voxels(CompFab::VoxelGrid *voxels, unsigned int voxel_distance) {
    // mark inner voxels
    for (unsigned int i = voxel_distance; i < voxels->m_dimX - voxel_distance; i++) {
        std::cout << "\rMarking Inner Voxels with Distance " << voxel_distance << ": " << ((i + 1) * 100) / voxels->m_dimX << "%" << std::flush;
        for (unsigned int j = voxel_distance; j < voxels->m_dimY - voxel_distance; j++) {
            for (unsigned int k = voxel_distance; k < voxels->m_dimZ - voxel_distance; k++) {
                // if (i+voxel_distance >= voxels->m_dimX || j+voxel_distance >= voxels->m_dimY || k+voxel_distance >= voxels->m_dimZ ||
                //     i < voxel_distance || j < voxel_distance || k < voxel_distance) {
                //     voxels->isInner(i,j,k) = false;
                //     continue;
                // }
                bool inner = true;
                for (unsigned int ii = i - voxel_distance; ii <= i + voxel_distance; ii++) {
                    for (unsigned int jj = j - voxel_distance; jj <= j + voxel_distance; jj++) {
                        for (unsigned int kk = k - voxel_distance; kk <= k + voxel_distance; kk++) {
                            inner = inner && voxels->isInside(ii, jj, kk);
                        }
                    }
                }
                voxels->isInner(i, j, k) = inner;
            }
        }
    }
    std::cout << std::endl;
}

// keeps only the base of the object
CompFab::VoxelGrid *compute_base() {
    CompFab::VoxelGrid *g_voxelGridBase = new CompFab::VoxelGrid(g_voxelGrid->m_lowerLeft, g_voxelGrid->m_dimX, g_voxelGrid->m_dimY, g_voxelGrid->m_dimZ, g_voxelGrid->m_spacing);

    base_layer = -1;
    for (unsigned int j = 0; j < g_voxelGrid->m_dimY; j++) {
        for (unsigned int i = 0; i < g_voxelGrid->m_dimX; i++) {
            for (unsigned int k = 0; k < g_voxelGrid->m_dimZ; k++) {
                if (!g_voxelGrid->isInside(i, j, k)) {
                    continue;
                }

                if (base_layer != -1 && j > base_layer) {
                    return g_voxelGridBase;
                }

                base_layer = j;
                g_voxelGridBase->isInside(i, j, k) = true;
            }
        }
    }

    return g_voxelGridBase;
}

void dilate(CompFab::VoxelGrid *voxels, int layer, CompFab::Vec3 direction) {
    direction.normalize();
    // copy the voxels
    CompFab::VoxelGrid *old_voxels = new CompFab::VoxelGrid(voxels->m_lowerLeft, voxels->m_dimX, voxels->m_dimY, voxels->m_dimZ, voxels->m_spacing);
    for (unsigned int i = 0; i < voxels->m_dimX; i++) {
        for (unsigned int j = 0; j < voxels->m_dimY; j++) {
            for (unsigned int k = 0; k < voxels->m_dimZ; k++) {
                old_voxels->isInside(i, j, k) = voxels->isInside(i, j, k);
            }
        }
    }

    for (unsigned int i = 0; i < voxels->m_dimX; i++) {
        for (unsigned int k = 0; k < voxels->m_dimZ; k++) {
            if (!old_voxels->isInside(i, layer, k)) {
                continue;
            }

            // Calculate the next voxel position based on the direction vector
            int nextI = static_cast<int>(std::round(i + direction[0]));
            int nextK = static_cast<int>(std::round(k + direction[2]));

            // Check bounds and set the next voxel as inside if within bounds
            if (nextI >= 0 && nextI < voxels->m_dimX && nextK >= 0 && nextK < voxels->m_dimZ) {
                if (layer == base_layer || voxels->isInside(nextI, layer - 1, nextK)) {
                    voxels->isInside(nextI, layer, nextK) = true;
                }
            }
        }
    }
}

CompFab::VoxelGrid *erode2d(CompFab::VoxelGrid *base, int iterations) {
    int layer = base_layer;
    CompFab::VoxelGrid *eroded = new CompFab::VoxelGrid(base->m_lowerLeft, base->m_dimX, base->m_dimY, base->m_dimZ, base->m_spacing);

    for (unsigned int i = 1; i < base->m_dimX - 1; i++) {
        for (unsigned int k = 1; k < base->m_dimZ - 1; k++) {
            if (!base->isInside(i, layer, k)) {
                continue;
            }

            bool inside = true;
            for (unsigned int ii = i - 1; ii <= i + 1; ii++) {
                for (unsigned int kk = k - 1; kk <= k + 1; kk++) {
                    inside = inside && base->isInside(ii, layer, kk);
                }
            }
            eroded->isInside(i, layer, k) = inside;
        }
    }

    if (iterations > 1) {
        CompFab::VoxelGrid *old_eroded = eroded;
        eroded = erode2d(old_eroded, iterations - 1);
        delete old_eroded;
    }

    return eroded;
}

CompFab::VoxelGrid *fill_base(CompFab::VoxelGrid *base) {
    CompFab::VoxelGrid *convex_hull = new CompFab::VoxelGrid(base->m_lowerLeft, base->m_dimX, base->m_dimY, base->m_dimZ, base->m_spacing);

    // array of points
    std::vector<Point> base_points;

    // get base points
    int base_layer = -1;
    for (int j = 0; j < base->m_dimY; j++) {
        for (int i = 0; i < base->m_dimX; i++) {
            for (int k = 0; k < base->m_dimZ; k++) {
                if (!base->isInside(i, j, k)) {
                    continue;
                }

                if (base_layer != -1 && j > base_layer) {
                    break;
                }

                base_layer = j;
                base_points.push_back({i, k});
            }
            if (base_layer != -1 && j > base_layer) break;
        }
        if (base_layer != -1 && j > base_layer) break;
    }

    int n = base_points.size();
    std::vector<Point> vertices = convexHull(base_points, n);

    std::cout << "Base Points: " << n << std::endl;
    std::cout << "Convex Hull Points: " << vertices.size() << std::endl;
    for (int i = 0; i < vertices.size(); i++) {
        std::cout << "(" << vertices[i].x << " " << vertices[i].y << ") ";
    }
    std::cout << std::endl;

    // fill base

    // int i, j, c = 0;
    // for (i = 0, j = nvert-1; i < nvert; j = i++) {
    //     if ( ((verty[i]>testy) != (verty[j]>testy)) &&
    //     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
    //     c = !c;
    // }

    for (int i = 0; i < base->m_dimX; i++) {
        for (int k = 0; k < base->m_dimZ; k++) {
            if (base->isInside(i, base_layer, k)) {
                convex_hull->isInside(i, base_layer, k) = true;
                continue;
            }

            int l, m = 0;
            bool inside = false;
            int nvert = vertices.size();
            for (l = 0, m = nvert - 1; l < nvert; m = l++) {
                if (((vertices[l].y > k) != (vertices[m].y > k)) && (i < (vertices[m].x - vertices[l].x) * (k - vertices[l].y) / (vertices[m].y - vertices[l].y) + vertices[l].x)) {
                    inside = !inside;
                }
            }
            convex_hull->isInside(i, base_layer, k) = inside;
        }
    }

    return convex_hull;
}

void fill_holes() {
    for (unsigned int j = 0; j < g_voxelGrid->m_dimY; j++) {
        for (unsigned int i = 0; i < g_voxelGrid->m_dimX; i++) {
            for (unsigned int k = 0; k < g_voxelGrid->m_dimZ; k++) {
                if (g_voxelGrid->isInside(i, j, k)) {
                    continue;
                }
                // check if corners
                if (i == 0 || i == g_voxelGrid->m_dimX - 1 || j == 0 || j == g_voxelGrid->m_dimY - 1 || k == 0 || k == g_voxelGrid->m_dimZ - 1) {
                    continue;
                }

                // check for neighbors

                bool top_bottom = g_voxelGrid->isInside(i, j - 1, k) && g_voxelGrid->isInside(i, j + 1, k);
                bool left_right = g_voxelGrid->isInside(i - 1, j, k) && g_voxelGrid->isInside(i + 1, j, k);
                bool front_back = g_voxelGrid->isInside(i, j, k - 1) && g_voxelGrid->isInside(i, j, k + 1);
                if (top_bottom && left_right && front_back) {
                    std::cout << "Filling hole at " << i << " " << j << " " << k << std::endl;
                    g_voxelGrid->isInside(i, j, k) = true;
                }
            }
        }
    }
}

// compute barycenter of the voxels
std::tuple<CompFab::Vec3, CompFab::Vec3, unsigned int> compute_barycenter() {
    double x, y, z;
    x = y = z = 0;
    unsigned int count = 0;

    // iterate over all voxels
    for (unsigned int i = 0; i < g_voxelGrid->m_dimX; i++) {
        for (unsigned int j = 0; j < g_voxelGrid->m_dimY; j++) {
            for (unsigned int k = 0; k < g_voxelGrid->m_dimZ; k++) {
                if (!g_voxelGrid->isInside(i, j, k)) {
                    continue;
                }
                x += i + 0.5;
                y += j + 0.5;
                z += k + 0.5;
                count++;
            }
        }
    }
    CompFab::Vec3 barycenter = CompFab::Vec3(x / count, y / count, z / count);
    return std::tuple<CompFab::Vec3, CompFab::Vec3, unsigned int>(barycenter, CompFab::Vec3(x, y, z), count);
}

bool is_balanced(CompFab::Vec3 barycenter, CompFab::VoxelGrid *base) {
    double x = barycenter[0];
    double y = barycenter[1];
    double z = barycenter[2];

    // chck if barycenter is inside the base
    return base->isInside(x, base_layer, z);
}

double compute_distance_2d(CompFab::Vec3 a, CompFab::Vec3 b) {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[2] - b[2]) * (a[2] - b[2]));
}

// find direction from barycenter to closest point on the base
CompFab::Vec3 compute_direction_to_base(CompFab::Vec3 barycenter, CompFab::VoxelGrid *base_filled) {
    double x = barycenter[0];
    double z = barycenter[2];
    barycenter = CompFab::Vec3(x, base_layer, z);
    // find closest point on the base
    double min_distance = 1e9;
    // CompFab::Vec3 closest_point;
    std::vector<std::pair<CompFab::Vec3, double>> closest_points;
    for (unsigned int i = 0; i < base_filled->m_dimX; i++) {
        for (unsigned int k = 0; k < base_filled->m_dimZ; k++) {
            if (!base_filled->isInside(i, base_layer, k)) {
                continue;
            }
            CompFab::Vec3 point(i + 0.5, base_layer, k + 0.5);

            double distance = compute_distance_2d(point, barycenter);
            closest_points.push_back(std::make_pair(point, distance));
        }
    }

    // sort the closest points
    std::sort(closest_points.begin(), closest_points.end(), [](const std::pair<CompFab::Vec3, double> &a, const std::pair<CompFab::Vec3, double> &b) {
        return a.second < b.second;
    });

    // get the n closest points
    int n = std::min(1, (int)closest_points.size());
    std::vector<CompFab::Vec3> closest_points_n;
    for (unsigned int i = 0; i < n; i++) {
        closest_points_n.push_back(closest_points[i].first);
    }
    // average the n closest points
    CompFab::Vec3 closest_point(0, 0, 0);
    for (unsigned int i = 0; i < n; i++) {
        closest_point += closest_points_n[i];
    }
    closest_point[0] /= n;
    closest_point[1] /= n;
    closest_point[2] /= n;

    CompFab::Vec3 direction = closest_point - barycenter;
    direction.normalize();
    return direction;
}

// sort the voxels in increasing order of their signed distance to the cutting plane
std::vector<CompFab::Vec3> sort_voxels(CompFab::Vec3 barycenter, CompFab::Vec3 direction) {
    std::vector<std::pair<CompFab::Vec3, double>> voxels;

    for (unsigned int i = 0; i < g_voxelGrid->m_dimX; i++) {
        for (unsigned int j = 0; j < g_voxelGrid->m_dimY; j++) {
            for (unsigned int k = 0; k < g_voxelGrid->m_dimZ; k++) {
                if (!g_voxelGrid->isInner(i, j, k)) {
                    continue;
                }

                CompFab::Vec3 voxel(i, j, k);
                CompFab::Vec3 barycenter_to_voxel = voxel - barycenter;
                double distance = barycenter_to_voxel * direction;
                if (distance < 0) {
                    // voxels that are on the base side of the cutting plane will be ignored
                    // since the direction is from the barycenter to the base, we need to invert the sign of the distance
                    voxels.push_back(std::make_pair(voxel, -distance));
                }
            }
        }
    }

    // sort the voxels vector and distances vector in decreasing order of distances
    std::sort(voxels.begin(), voxels.end(), [](const std::pair<CompFab::Vec3, double> &a, const std::pair<CompFab::Vec3, double> &b) {
        if (a.second == b.second) {
            return a.first[1] < b.first[1];
        }

        return a.second < b.second;
    });

    std::vector<CompFab::Vec3> sorted_voxels;
    for (unsigned int i = 0; i < voxels.size(); i++) {
        sorted_voxels.push_back(voxels[i].first);
    }

    // print closest pixels
    // for (unsigned int i = 0; i < 20; i++) {
    //     std::cout << "Sorted: " << sorted_voxels[i][0] << " " << sorted_voxels[i][1] << " " << sorted_voxels[i][2] << " | " << voxels[i].second << std::endl;
    // }

    return sorted_voxels;
}

void expand_base(CompFab::VoxelGrid *base, int max_iterations) {
    CompFab::Vec3 barycenter = std::get<0>(compute_barycenter());
    CompFab::Vec3 direction = compute_direction_to_base(barycenter, base);

    std::cout << "Expanding Base" << std::endl;
    for (unsigned int i = 0; i < max_iterations; i++) {
        dilate(base, base_layer, -direction);
        dilate(g_voxelGrid, base_layer, -direction);
        for (unsigned int j = 1; j < i; j++) {
            dilate(base, base_layer + j, -direction);
            dilate(g_voxelGrid, base_layer + j, -direction);
        }
        barycenter = std::get<0>(compute_barycenter());
        if (is_balanced(barycenter, base)) {
            std::cout << "Balanced with base dilations at iteration " << i << std::endl;
            break;
        }
    }
    if (!is_balanced(barycenter, base)) {
        std::cout << "Not Balanced after " << max_iterations << " iterations" << std::endl;
    }
}

bool carve(CompFab::VoxelGrid *base) {
    int i = 0;

    auto [barycenter, sum, count] = compute_barycenter();
    CompFab::Vec3 direction = compute_direction_to_base(barycenter, base);
    std::vector<CompFab::Vec3> sorted_voxels = sort_voxels(barycenter, direction);
    std::cout << "Carving " << sorted_voxels.size() << " voxels" << std::endl;

    while (sorted_voxels.size() > 0) {
        i++;
        CompFab::Vec3 voxel = sorted_voxels.back();
        sorted_voxels.pop_back();
        if (!g_voxelGrid->isInside(voxel[0], voxel[1], voxel[2]) || !g_voxelGrid->isInner(voxel[0], voxel[1], voxel[2])) {
            std::cout << "Error: Voxel not inside" << std::endl;
            continue;
        }
        g_voxelGrid->isInside(voxel[0], voxel[1], voxel[2]) = false;
        g_voxelGrid->isInner(voxel[0], voxel[1], voxel[2]) = false;

        // update barycenter
        sum[0] -= voxel[0] + 0.5;
        sum[1] -= voxel[1] + 0.5;
        sum[2] -= voxel[2] + 0.5;
        count--;

        barycenter = CompFab::Vec3(sum[0] / count, sum[1] / count, sum[2] / count);
        if (is_balanced(barycenter, base)) {
            std::cout << "Balanced at iteration " << i << std::endl;
            break;
        } else if (i % 10 == 0) {
            direction = compute_direction_to_base(barycenter, base);
            sorted_voxels = sort_voxels(barycenter, compute_direction_to_base(barycenter, base));
        }
    }

    bool isbalanced = is_balanced(barycenter, base);

    if (!isbalanced) {
        std::cout << "Not Balanced" << std::endl;
    }

    return isbalanced;
}

void save_barycenter(CompFab::Vec3 barycenter, const char *outfile, CompFab::VoxelGrid *base) {
    double x = barycenter[0];
    double y = barycenter[1];
    double z = barycenter[2];

    // save point to file
    CompFab::VoxelGrid *g_voxelBarycenter = new CompFab::VoxelGrid(base->m_lowerLeft, base->m_dimX, base->m_dimY, base->m_dimZ, base->m_spacing);
    g_voxelBarycenter->isInside(x, y, z) = true;

    // copy base to barycenter
    for (unsigned int i = 0; i < base->m_dimX; i++) {
        for (unsigned int k = 0; k < base->m_dimZ; k++) {
            if (base->isInside(i, base_layer, k)) {
                g_voxelBarycenter->isInside(i, base_layer, k) = true;
                // std::cout << "Base: " << i << " " << j << " " << k << std::endl;
            }
        }
    }

    // project barycenter to base
    g_voxelBarycenter->isInside(x, base_layer, z) = false;

    saveVoxelsToObj(outfile, g_voxelBarycenter);

    // barycenter /= g_voxelGrid->m_size;
    std::cout << "barycenter: " << x << " " << y << " " << z << std::endl;
    // convert to voxel coordinates
    x = 0.5 + x * g_voxelGrid->m_spacing;
    y = 0.5 + y * g_voxelGrid->m_spacing;
    z = 0.5 + z * g_voxelGrid->m_spacing;
    std::cout << "barycenter: " << x << " " << y << " " << z << std::endl;
}

bool loadMesh(char *filename, unsigned int dim) {
    g_triangleList.clear();

    Mesh *tempMesh = new Mesh(filename, true);

    CompFab::Vec3 v1, v2, v3;

    // copy triangles to global list
    for (unsigned int tri = 0; tri < tempMesh->t.size(); ++tri) {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1, v2, v3));
    }

    // Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);

    // Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;

    if (bbX > bbY && bbX > bbZ) {
        spacing = bbX / (double)(dim - 2);
    } else if (bbY > bbX && bbY > bbZ) {
        spacing = bbY / (double)(dim - 2);
    } else {
        spacing = bbZ / (double)(dim - 2);
    }

    CompFab::Vec3 hspacing(0.5 * spacing, 0.5 * spacing, 0.5 * spacing);

    g_voxelGrid = new CompFab::VoxelGrid(bbMin - hspacing, dim, dim, dim, spacing);

    delete tempMesh;

    return true;
}

void saveVoxelsToObj(const char *outfile, CompFab::VoxelGrid *g_voxelGrid) {
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;

    CompFab::Vec3 hspacing(0.5 * spacing, 0.5 * spacing, 0.5 * spacing);

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if (!g_voxelGrid->isInside(ii, jj, kk)) {
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii) * spacing, 0.5f + ((double)jj) * spacing, 0.5f + ((double)kk) * spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}

void saveInnerVoxelsToObj(const char *outfile, CompFab::VoxelGrid *g_voxelGrid) {
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;

    CompFab::Vec3 hspacing(0.5 * spacing, 0.5 * spacing, 0.5 * spacing);

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if (!g_voxelGrid->isInner(ii, jj, kk)) {
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii) * spacing, 0.5f + ((double)jj) * spacing, 0.5f + ((double)kk) * spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}

void saveInnerVoxelsToObj(const char *outfile) {
    saveInnerVoxelsToObj(outfile, g_voxelGrid);
}

void saveVoxelsToObj(const char *outfile) {
    saveVoxelsToObj(outfile, g_voxelGrid);
}

void saveVoxels(const char *outfile) {
    saveVoxels(outfile, g_voxelGrid);
}

void saveVoxels(const char *outfile, CompFab::VoxelGrid *voxels) {
    // save cpp object to file
    std::ofstream file(outfile);
    if (file.is_open()) {
        file << voxels->m_dimX << " " << voxels->m_dimY << " " << voxels->m_dimZ << " " << voxels->m_spacing << std::endl;
        for (unsigned int i = 0; i < voxels->m_dimX; i++) {
            for (unsigned int j = 0; j < voxels->m_dimY; j++) {
                for (unsigned int k = 0; k < voxels->m_dimZ; k++) {
                    file << voxels->isInside(i, j, k) << " ";
                }
            }
        }
        for (unsigned int i = 0; i < voxels->m_dimX; i++) {
            for (unsigned int j = 0; j < voxels->m_dimY; j++) {
                for (unsigned int k = 0; k < voxels->m_dimZ; k++) {
                    file << voxels->isInner(i, j, k) << " ";
                }
            }
        }
        file.close();
    }
}

void loadVoxels(const char *filename) {
    // load cpp object from file
    std::ifstream file(filename);
    if (file.is_open()) {
        unsigned int dimX, dimY, dimZ;
        double spacing;
        file >> dimX >> dimY >> dimZ >> spacing;
        g_voxelGrid = new CompFab::VoxelGrid(CompFab::Vec3(0, 0, 0), dimX, dimY, dimZ, spacing);
        for (unsigned int i = 0; i < dimX; i++) {
            for (unsigned int j = 0; j < dimY; j++) {
                for (unsigned int k = 0; k < dimZ; k++) {
                    int value;
                    file >> value;
                    g_voxelGrid->isInside(i, j, k) = value;
                }
            }
        }
        for (unsigned int i = 0; i < dimX; i++) {
            for (unsigned int j = 0; j < dimY; j++) {
                for (unsigned int k = 0; k < dimZ; k++) {
                    int value;
                    file >> value;
                    g_voxelGrid->isInner(i, j, k) = value;
                }
            }
        }
        file.close();
    }
}

// set deformation
// deformation origin
// CompFab::Vec3 deformation_origin;
// double deformation_radius = 0.0;
// void set
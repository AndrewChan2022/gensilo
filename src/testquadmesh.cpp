#include <silo.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>


static int
test_quadmesh(DBfile *dbfile)
{
    int			nerrors=0;
    int dims[] = {5, 5, 5};
    float coords0[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    float coords1[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    float coords2[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    // static float	*coords[] = {coords0, coords1, coords2};

    const int nx = 5;  // number of nodes in x
    const int ny = 5;  // number of nodes in y
    const int nz = 5;  // number of nodes in z

    // Create coordinate arrays (uniform spacing from 0 to 1)
    std::vector<float> x_coords(nx);
    std::vector<float> y_coords(ny);
    std::vector<float> z_coords(nz);

    // for (int i = 0; i < nx; i++) {
    //     coords0[i] = static_cast<float>(i) / (nx - 1) * 0.1;
    //     coords1[i] = static_cast<float>(i) / (ny - 1) * 0.1;
    //     coords2[i] = static_cast<float>(i) / (nz - 1) * 0.1;
    // }

    // Prepare coordinate arrays for Silo
    float* coords[3] = {coords0, coords1, coords2};
    // int dims[3] = {nx, ny, nz};

    // static double	varA[] = {10, 11, 12, 13, 14,
	// 			  15, 16, 17, 18, 19,
	// 			  20, 21, 22, 23, 24,
	// 			  25, 26, 27, 28, 29,
	// 			  30, 31, 32, 33, 34};
    // static double	varB[] = {35, 36, 37, 38, 39,
	// 			  40, 41, 42, 43, 44,
	// 			  45, 46, 47, 48, 49,
	// 			  50, 51, 52, 53, 54,
	// 			  55, 56, 57, 58, 59};
    // static double	varC[] = {60, 61, 62, 63, 64,
	// 			  65, 66, 67, 68, 69,
	// 			  70, 71, 72, 73, 74,
	// 			  75, 76, 77, 78, 79,
	// 			  80, 81, 82, 83, 84};
    // static double	*vars[] = {varA, varB, varC};
    // static char		*varnames[] = {"varA", "varB", "varC"};
	
    puts("=== Quadmesh ===");

    // DBMkDir(dbfile, "/quad");
    // DBSetDir(dbfile, "/quad");

    // Optional: Set mesh extents explicitly
    double min_extent = 0.1;
    double max_extent = 0.5;
    double extents[6] = {
        min_extent, min_extent, min_extent,
        max_extent, max_extent, max_extent};
    DBoptlist* optlist = DBMakeOptlist(1);
    DBAddOption(optlist, DBOPT_EXTENTS_SIZE, extents);

    if (DBPutQuadmesh(dbfile, "qm1", nullptr, coords, dims, 3, DB_FLOAT, DB_COLLINEAR, optlist)<0) {
	    puts("DBPutQuadmesh(qm1) failed");
	    nerrors++;
    }
    DBFreeOptlist(optlist);

    // if (DBPutQuadvar(dbfile, "qv1", "qm1", 3, (DBCAS_t) varnames,
    //     vars, dims, 3, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL)<0) {
	// puts("DBPutQuadmesh(qv1) failed");
	// nerrors++;
    // }
    

    return nerrors;
}

int test_quadmesh_main() {
    // Initialize random seed
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    // Create a new Silo file (using HDF5 driver)
    DBfile* file = DBCreate("rectilinear_grid6.silo", DB_CLOBBER, DB_LOCAL, 
                           "10x10x10 rectilinear grid with random scalars", DB_HDF5);

    if (!file) {
        std::cerr << "Failed to create Silo file!" << std::endl;
        return 1;
    }

    test_quadmesh(file);

    // // Define grid dimensions (10x10x10 zones means 11x11x11 nodes)
    // const int nx = 11;  // number of nodes in x
    // const int ny = 11;  // number of nodes in y
    // const int nz = 11;  // number of nodes in z

    // // Create coordinate arrays (uniform spacing from 0 to 1)
    // std::vector<float> x_coords(nx);
    // std::vector<float> y_coords(ny);
    // std::vector<float> z_coords(nz);

    // for (int i = 0; i < nx; i++) {
    //     x_coords[i] = static_cast<float>(i) / (nx - 1);
    //     y_coords[i] = static_cast<float>(i) / (ny - 1);
    //     z_coords[i] = static_cast<float>(i) / (nz - 1);
    // }

    // // Prepare coordinate arrays for Silo
    // float* coords[3] = {x_coords.data(), y_coords.data(), z_coords.data()};

    // // Define mesh dimensions (number of zones in each direction)
    // int dims[3] = {nx, ny, nz};

    // // Write the rectilinear mesh
    // DBPutQuadmesh(file, "rectilinear_mesh", nullptr, coords, dims, 3,
    //               DB_FLOAT, DB_COLLINEAR, nullptr);

    // // Create random scalar values for each zone (10x10x10 zones)
    // const int n_zones = 10 * 10 * 10;
    // std::vector<float> scalar_values(n_zones);

    // for (int i = 0; i < n_zones; i++) {
    //     scalar_values[i] = static_cast<float>(std::rand()) / RAND_MAX; // 0.0 to 1.0
    // }

    // // Define variable dimensions (number of zones)
    // int var_dims[3] = {10, 10, 10};  // One less than node dimensions

    // // Write the scalar variable
    // DBPutQuadvar1(file, "scalar_field", "rectilinear_mesh", scalar_values.data(),
    //               var_dims, 3, nullptr, 0, DB_FLOAT, DB_ZONECENT, nullptr);

    // Close the file
    DBClose(file);

    std::cout << "Successfully created rectilinear_grid.silo with 10x10x10 grid" << std::endl;

    return 0;
}
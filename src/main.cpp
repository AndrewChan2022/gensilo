
int testmat(int argc, char **argv);
int test_quadmesh_main();
int testvoxelmat(int argc, char **argv);


#include <silo.h>
#include <stdio.h>
#include <memory>
#include <vector>

void WriteRectLinear3D(DBfile *dbfile) {
    const int NX = 4;
    const int NY = 3;
    const int NZ = 2;
    // Write a rectilinear mesh
    float x[] = {0.f, 1.f, 3.f, 6.f};
    float y[] = {0.f, 1.f, 3.25f};
    float z[] = {0.f, 1.f};
    int dims[] = {NX, NY, NZ};
    int ndims = 3;
    float *coords[] = {x, y, z};
    DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims,
        DB_FLOAT, DB_COLLINEAR, NULL);
}

void WriteCurveLinear2D(DBfile *dbfile) {
    // Write a curvilinear mesh.
    #define NX 4
    #define NY 3
    float x[NY][NX] = {{0., 1., 3., 3.5}, {0., 1., 2.5, 3.5},
        {0.7, 1.3, 2.3, 3.5}};
    float y[NY][NX] = {{0., 0., 0., 0.}, {1.5, 1.5, 1.25, 1.5},
        {3., 2.75, 2.75, 3.}};
    int dims[] = {NX, NY};
    int ndims = 2;
    float *coords[] = {(float*)x, (float*)y};
    DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims,
        DB_FLOAT, DB_NONCOLLINEAR, NULL);
}

void WriteCurveLinear3D(DBfile *dbfile) {
    // Write a curvilinear mesh.
    #define NX 4
    #define NY 3
    #define NZ 2
    float x[NZ][NY][NX] = {
        {{0.,1.,2.,3.},{0.,1.,2.,3.}, {0.,1.,2.,3.}},
        {{0.,1.,2.,3.},{0.,1.,2.,3.}, {0.,1.,2.,3.}}
    };
    float y[NZ][NY][NX] = {
        {{0.5,0.,0.,0.5},{1.,1.,1.,1.}, {1.5,2.,2.,1.5}},
        {{0.5,0.,0.,0.5},{1.,1.,1.,1.}, {1.5,2.,2.,1.5}}
    };
    float z[NZ][NY][NX] = {
        {{0.,0.,0.,0.},{0.,0.,0.,0.},{0.,0.,0.,0.}},
        {{1.,1.,1.,1.},{1.,1.,1.,1.},{1.,1.,1.,1.}}
    };
    int dims[] = {NX, NY, NZ};
    int ndims = 3;
    float *coords[] = {(float*)x, (float*)y, (float*)z};
    DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims,
        DB_FLOAT, DB_NONCOLLINEAR, NULL);
}

void WriteUgrid3D(DBfile *dbfile) {
    // Node coordinates
    float x[] = {0.,2.,2.,0.,0.,2.,2.,0.,0.,2.,2.,0.,1.,2.,4.,4.};
    float y[] = {0.,0.,0.,0.,2.,2.,2.,2.,4.,4.,4.,4.,6.,0.,0.,0.};
    float z[] = {2.,2.,0.,0.,2.,2.,0.,0.,2.,2.,0.,0.,1.,4.,2.,0.};
    float *coords[] = {x, y, z};

    // Connectivity 
    int nodelist[] = {
        1,2,3,4,5,6,7,8,    // hex, zone 1
        5,6,7,8,9,10,11,12, // hex, zone 2
        9,10,11,12,13,      // pyramid, zone 3
        2,3,16,15,6,7,      // prism, zone 4
        2,15,14,6           // tet, zone 5
    };
    int lnodelist = sizeof(nodelist) / sizeof(int);
    // shape type 1 has 8 nodes (hex)
    // shape type 2 has 5 nodes (pyramid)
    // shape type 3 has 6 nodes (prism)
    // shape type 4 has 4 nodes (tet)
    int shapesize[] = {8,5,6,4};

    
    // We have 2 hex, 1 pyramid, 1 prism, 1 tet
    int shapecounts[] = {2,1,1,1};
    int nshapetypes = 4;
    int nnodes = 16;
    int nzones = 5;
    int ndims = 3;
    // Write out connectivity information.
    DBPutZonelist(dbfile, "zonelist", nzones, ndims, nodelist, lnodelist,
                1, shapesize, shapecounts, nshapetypes);
    // Write an unstructured mesh.
    DBPutUcdmesh(dbfile, "mesh", ndims, NULL, coords, nnodes, nzones,
                "zonelist", NULL, DB_FLOAT, NULL);
}

#define NX 4
#define NY 3
#define NZ 2

// Write a zone-centered variable.
void write_zonecent_quadvar(DBfile *dbfile) {
    int i, dims[3], ndims = 3;
    int ncells = (NX-1)*(NY-1)*(NZ-1);
    float *data = (float *)malloc(sizeof(float)*ncells);
    for(i = 0; i < ncells; ++i)
        data[i] = (float)i;
    dims[0] = NX-1; dims[1] = NY-1; dims[2] = NZ-1;
    DBPutQuadvar1(dbfile, "zonal", "quadmesh", data, dims,
                  ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);
    free(data);
}

// Write a node-centered variable.
void write_nodecent_quadvar(DBfile *dbfile) {
    int i, dims[3], ndims = 3;
    int nnodes = NX*NY*NZ;
    float *data = (float *)malloc(sizeof(float)*nnodes);
    for(i = 0; i < nnodes; ++i)
        data[i] = (float)i;
    dims[0] = NX; dims[1] = NY; dims[2] = NZ;
    DBPutQuadvar1(dbfile, "nodal", "quadmesh", data, dims,
                  ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
    free(data);
}

void write_u_nodecent_quadvar(DBfile *dbfile, int nnodes, int nzones) {

    float *nodal = (float *)malloc(sizeof(float)*nnodes);
    for(int i = 0; i < nnodes; ++i) {
        nodal[i] = (float)i;
    }
    // Write a node-centered variable.
    DBPutUcdvar1(dbfile, "nodal", "mesh", nodal, nnodes, NULL, 0,
        DB_FLOAT, DB_NODECENT, NULL);

    
    // Write a zone-centered variable.
    float *zonal = (float *)malloc(sizeof(float)*nzones);
    for(int i = 0; i < nzones; ++i) {
        zonal[i] = (float)i;
    }
    DBPutUcdvar1(dbfile, "zonal", "mesh", zonal, nzones, NULL, 0,
        DB_FLOAT, DB_ZONECENT, NULL);
    
    free(nodal);
    free(zonal);
}

void write_cgrid3d_with_node_scalar(DBfile *dbfile) {
    WriteCurveLinear3D(dbfile);
    write_nodecent_quadvar(dbfile);
}

void write_rgrid3d_with_node_scalar(DBfile *dbfile) {
    WriteRectLinear3D(dbfile);
    write_nodecent_quadvar(dbfile);
}

void write_ugrid3d_with_node_scalar(DBfile *dbfile) {
    WriteUgrid3D(dbfile);
    write_u_nodecent_quadvar(dbfile, 16, 5);
}

void write_rgrid_material(DBfile *dbfile) {
    const int nx = 101;
    const int ny = 101;
    const int nz = 101;
    char *meshname = "mesh";
    char *matname = "mat";

    // mesh
    std::vector<float> x(nx);
    std::vector<float> y(ny);
    std::vector<float> z(nz);
    for (size_t i = 0; i < nx; i++) {
        x[i] = static_cast<float>(nx-i)/nx;
        printf("x[i]:%f\n", x[i]);
    }
    for (size_t i = 0; i < ny; i++) {
        y[i] = static_cast<float>(ny-i)/ny;
    }
    for (size_t i = 0; i < nz; i++) {
        z[i] = static_cast<float>(nz-i)/nz;
    }
    int dims[] = {nx, ny, nz};
    int ndims = 3;
    float *coords[] = {x.data(), y.data(), z.data()};
    DBoptlist *optList = DBMakeOptlist(6);
    DBAddOption(optList, DBOPT_XLABEL, (void*)"Width");
    DBAddOption(optList, DBOPT_YLABEL, (void*)"Height");
    DBAddOption(optList, DBOPT_ZLABEL, (void*)"Depth");
    DBAddOption(optList, DBOPT_XUNITS, (void*)"parsec");
    DBAddOption(optList, DBOPT_YUNITS, (void*)"parsec");
    DBAddOption(optList, DBOPT_ZUNITS, (void*)"parsec");
    DBPutQuadmesh(dbfile, meshname, NULL, coords, dims, ndims,
        DB_FLOAT, DB_COLLINEAR, NULL);
    DBFreeOptlist(optList);
    
    // material at zones
    int nmats = 2;
    int *matnos = new int[nmats];
    for (int i = 0 ; i < nmats ; i++) {
        matnos[i] = i+1;
    }

    int nzones = (nx-1)*(ny-1)*(nz-1);
    int *matlist = new int[nzones];
    dims[0] = nx-1;
    dims[1] = ny-1;
    dims[2] = nz-1;
    // matlist is matno/mixindex of each zone
    for (int i = 0 ; i < dims[0] ; i++) {
        for (int j = 0 ; j < dims[1] ; j++) {
            for (int k = 0 ; k < dims[2] ; k++) {
                int index = k*ny*nx + j*nx + i;
                int index2 = k*(ny-1)*(nx-1) + j*(nx-1) + i;
                int mat = 0;
                if (i < nx/2) {
                    mat = matnos[0];
                } else if (i > nx/2) {
                    mat = matnos[1];
                } else {
                    // 
                    mat = matnos[1];
                }
                matlist[index2] = mat;
            }
        }
    }

    // DBPutMaterial(dbfile, matname, meshname, nmats, matnos, matlist, dims, 3,
    //               NULL, NULL, NULL, NULL, 0, DB_FLOAT, NULL);

    // Clean up
    delete[] matnos;
    delete[] matlist;
}

int test_silo_(int argc, char *argv[]) {
    DBfile *dbfile = NULL;


    // Open the Silo file
    dbfile = DBCreate("rmat.silo", DB_CLOBBER, DB_LOCAL,
                      "Comment about the data", DB_HDF5);
    if(dbfile == NULL) {
        fprintf(stderr, "Could not create Silo file!\n");
        return -1;
    }

    // Add other Silo calls here.
    write_rgrid_material(dbfile);

    // Close the Silo file.
    DBClose(dbfile);

    return 0;
}


/// rectlinear save error: min_extents/max_extents of yz always 0
int main(int argc, char *argv[]) {

    printf("hello\n");

    testvoxelmat(argc, argv);

    return 0;
}

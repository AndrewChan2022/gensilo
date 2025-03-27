#include <stdlib.h>
#include <string.h>
#include "silo.h"
#include "std.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

const std::string input_file = "D:/data/github/mc/trinity/data/voxel_test.txt";
const std::string output_file = "voxel_mat5.pdb.silo";
#define zx 144
#define zy 75
#define zz 273
#define Z_CHANGING_FAST 1

constexpr int nx = (zx + 1);
constexpr int ny = (zy + 1);
constexpr int nz = (zz + 1);

class ReadMaterial {
public:
    ReadMaterial();
    void read(const std::string& file);

public:

    // dims
    std::vector<int> dims = {nx, ny, nz};

    // node
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;

    // material
    int nmat=3;
    std::vector<int> matnos;    // Dynamically populated

    std::vector<int> matlist;   // matno or mix index

    int mixlen;                 // size of mix_next
    std::vector<int> mix_mat;   // matno
    std::vector<int> mix_next;  // list of each zone
    std::vector<float> mix_vf;  // vf 
    std::vector<int> mix_zone;  // zone
};

ReadMaterial::ReadMaterial() {
    // Size for non-collinear mesh: full 3D coordinate arrays
    const size_t total_nodes = nx * ny * nz;
    x.resize(total_nodes);
    y.resize(total_nodes);
    z.resize(total_nodes);

    // Initialize coordinates for a 3D grid
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                size_t idx = i + j * nx + k * nx * ny;
                x[idx] = static_cast<float>(i);
                y[idx] = static_cast<float>(j);
                z[idx] = static_cast<float>(k);
            }
        }
    }

    // Initialize matlist
    matlist.resize(zx * zy * zz);
}

void ReadMaterial::read(const std::string& file) {

    // file is text file
    // each line is fraction of mat1, mat1, mat2.  fraction is 0~129, 0 is 0%, 129 is 100%, mat1 mat2 together 100%

    std::ifstream inFile(file);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file " << file << std::endl;
        return;
    }

    int line_idx = 0;  // Line number from the file (0 to 26)
    std::string line;
    int mix_idx = 0;   // Index into mix arrays
    const int total_zones = zx * zy * zz;
    int last_percent = -1;  // Track the last printed percentage

    // Clear mix arrays and matnos in case of re-read
    mix_mat.clear();
    mix_next.clear();
    mix_vf.clear();
    mix_zone.clear();
    matnos.clear();  // We'll build this dynamically

    // Read each line: fraction of mat1 (0-129), mat1_no, mat2_no
    while (std::getline(inFile, line) && line_idx < total_zones) {
        int frac1, mat1_no, mat2_no;
        if (sscanf(line.c_str(), "%d %d %d", &frac1, &mat1_no, &mat2_no) != 3) {
            std::cerr << "Error: Invalid line format at line " << line_idx + 1 << std::endl;
            continue;
        }

        // Validate fraction
        if (frac1 < 0 || frac1 > 129) {
            std::cerr << "Error: Fraction out of range (0-129) at line " << line_idx << std::endl;
            continue;
        }

        // Convert fractions to volume fractions (0.0 to 1.0)
        float vf1 = static_cast<float>(frac1) / 129.0f;
        float vf2 = static_cast<float>(129 - frac1) / 129.0f;

        // Find or add material numbers to matnos
        int mat1_idx = -1, mat2_idx = -1;
        for (size_t i = 0; i < matnos.size(); i++) {
            if (matnos[i] == mat1_no) mat1_idx = i;
            if (matnos[i] == mat2_no) mat2_idx = i;
        }
        if (mat1_idx == -1) {
            mat1_idx = matnos.size();
            matnos.push_back(mat1_no);
        }
        if (mat2_idx == -1 && mat2_no != mat1_no) {  // Avoid duplicates
            mat2_idx = matnos.size();
            matnos.push_back(mat2_no);
        } else if (mat2_no == mat1_no) {
            mat2_idx = mat1_idx;
        }

        // Update nmat
        nmat = matnos.size();

        // Map line_idx to zone_idx with z changing fastest
#if Z_CHANGING_FAST
        int k = line_idx % zz;              // z-index
        int j = (line_idx / zz) % zy;       // y-index
        int i = (line_idx / (zz * zy));     // x-index
        int zone_idx = i + j * zx + k * zx * zy;
#else
        int i = line_idx % zx;              // x-index
        int j = (line_idx / zx) % zy;       // y-index
        int k = (line_idx / (zx * zy));     // x-index
        int zone_idx = i + j * zx + k * zx * zy;
#endif

        // test
        // frac1 = 129;

        // Determine if this is a pure or mixed zone
        if(mat1_idx != 0)  {
            // printf("index:%d\n", mat1_idx);
        }
        if (frac1 == 129) {
            matlist[zone_idx] = mat1_no;  // Pure mat1
        } else if (frac1 == 0) {
            matlist[zone_idx] = mat2_no;  // Pure mat2
        } else {
            // Mixed zone
            matlist[zone_idx] = -(mix_idx + 1);  // Negative index (1-based) into mix arrays

            // Add mat1
            mix_mat.push_back(mat1_no);
            mix_vf.push_back(vf1);
            mix_zone.push_back(zone_idx);
            mix_next.push_back(mix_idx + 2);  // Point to mat2
            mix_idx++;

            // Add mat2
            if (mat2_no == 255) {
                // printf("mat2_no = 255\n");
            }
            mix_mat.push_back(mat2_no);
            mix_vf.push_back(vf2);
            mix_zone.push_back(zone_idx);
            mix_next.push_back(0);  // End of mix for this zone
            mix_idx++;
        }
        line_idx++;

        // printf("ijk:%d %d %d matlist:%d, mat1:%d mat2:%d fraction:%d\n", i, j, k, matlist[zone_idx], mat1_no, mat2_no, frac1);

        int current_percent = static_cast<int>(static_cast<float>(line_idx) / total_zones * 100.0f);
        if (current_percent != last_percent) {
            std::cout << "\rReading file: " << current_percent << "% (" 
                      << line_idx << "/" << total_zones << ")" << std::flush;
            last_percent = current_percent;
        }
    }

    mixlen = mix_idx;
    inFile.close();

    if (line_idx != zx * zy * zz) {
        std::cerr << "Warning: Expected " << zx * zy * zz << " zones, but read " << line_idx << std::endl;
    }

    // for (size_t k = 0; k < zz; k++) {
    //     printf("\nk:%d\n", k);
    //     for (size_t j = 0; j < zy; j++) {
    //         printf("\nj:%d\n", j);
    //         for (size_t i = 0; i < zx; i++) {
    //             size_t index = k * zy * zx + j * zx + i;
    //             printf(" %08x", matlist[index]);
    //             if (index % 8 == 0) {
    //                 printf("\n");
    //             }
    //         }
    //     }
    // }
}


#define ALLOC_N(t,n) (t*)(calloc(n, sizeof(n)));

//--------------------//
//    Main Program    //
//--------------------//
int testvoxelmat(int argc, char **argv) {
    DBfile *db;
    int            i, driver = DB_PDB;
    const char     *filename = output_file.c_str();
    int            show_all_errors = FALSE;
    int            custom_mat = FALSE;
    char const * const coordnames[3] = {"x", "y", "z"};
    float *coord[3];

    for (i=1; i<argc; i++) {
        if (!strncmp(argv[i], "DB_PDB",6)) {
            driver = StringToDriver(argv[i]);
            filename = "mat3d_3across.pdb";
        } else if (!strncmp(argv[i], "DB_HDF5", 7)) {
            driver = StringToDriver(argv[i]);
            filename = "mat3d_3across.h5";
        } else if (!strcmp(argv[i], "show-all-errors")) {
            show_all_errors = 1;
        } else if (!strcmp(argv[i], "custom-mat")) {
            custom_mat = 1;
	} else if (argv[i][0] != '\0') {
            fprintf(stderr, "%s: ignored argument `%s'\n", argv[0], argv[i]);
        }
    }

    if (show_all_errors) DBShowErrors(DB_ALL_AND_DRVR, 0);

    ReadMaterial reader;
    reader.read(input_file);
    coord[0] = reader.x.data();
    coord[1] = reader.y.data();
    coord[2] = reader.z.data();

    

    db=DBCreate(filename, DB_CLOBBER, DB_LOCAL,
                "Mixed zone 3d test", driver);

    DBPutQuadmesh(db, "mesh", coordnames, coord, reader.dims.data(), 3, 
                    DB_FLOAT, DB_NONCOLLINEAR, NULL);

    int dims[3];
    dims[0]=zx;
    dims[1]=zy;
    dims[2]=zz;

    int nmat = reader.nmat;
    int* matnos = reader.matnos.data();
    int* matlist = reader.matlist.data();
    int* mix_mat = reader.mix_mat.data();
    int* mix_next = reader.mix_next.data();
    float* mix_vf = reader.mix_vf.data();
    int*  mix_zone = reader.mix_zone.data();
    int mixlen = static_cast<int>(reader.mix_next.size());
    
    DBPutMaterial(db, "material", "mesh", nmat, matnos, matlist, dims, 3, 
        mix_next, mix_mat, mix_zone, mix_vf, mixlen, DB_FLOAT, NULL);

    if (custom_mat)
    {
        long count[3]; 
        DBobject *udef_matobj = DBMakeObject("userdef_material", DB_MATERIAL, 20);

        /* Standard material stuf (in order of args to DBPutMaterial but that doesn't matter) */
        DBAddStrComponent(udef_matobj, "meshid", "mesh");
        DBAddIntComponent(udef_matobj, "nmat", nmat);
        count[0] = nmat;
        DBWriteComponent(db, udef_matobj, "matnos", "userdef_material", "integer", matnos, 1, count);
        count[0] = dims[0]; count[1] = dims[1]; count[2] = dims[2];
        DBWriteComponent(db, udef_matobj, "matlist", "userdef_material", "integer", matlist, 3, count);
        count[0] = 3;
        DBWriteComponent(db, udef_matobj, "dims", "userdef_material", "integer", dims, 1, count);
        DBAddIntComponent(udef_matobj, "ndims", 3);
        count[0] = mixlen;
        DBWriteComponent(db, udef_matobj, "mix_next", "userdef_material", "integer", mix_next, 1, count);
        DBWriteComponent(db, udef_matobj, "mix_mat", "userdef_material", "integer", mix_mat, 1, count);
        DBWriteComponent(db, udef_matobj, "mix_zone", "userdef_material", "integer", mix_zone, 1, count);
        DBWriteComponent(db, udef_matobj, "mix_vf", "userdef_material", "float", mix_vf, 1, count);
        DBAddIntComponent(udef_matobj, "mixlen", mixlen);
        DBAddIntComponent(udef_matobj, "datatype", DB_FLOAT);

        /* Ok, lets write some extra arrays with interesting stuff */
        {
            char *strArray[] = {"mark","sandy","fred","steve","sue","JayLo"};
            char *tmpList = 0;
            int len;

            /* Add a simple integer valued scalar member */
            DBAddIntComponent(udef_matobj, "foo", 42);

            /* Add a simple double valued scalar member */
            DBAddDblComponent(udef_matobj, "M_PI", 3.1415926);

            /* Add a string valued component */
            DBAddStrComponent(udef_matobj, "make", "Toyota");

            /* Add an array of strings (Katie's case) */
            DBStringArrayToStringList((DBCAS_t) strArray, 6, &tmpList, &len);
            count[0] = len;
            DBWriteComponent(db, udef_matobj, "Katies_Names", "userdef_material", "char", tmpList, 1, count);
            free(tmpList);
        }

        /* Finally, write the generic object to the file */
        DBWriteObject(db, udef_matobj, 0);

        DBFreeObject(udef_matobj);
    }


    DBClose(db);

    printf("write to file:%s\n", output_file.c_str());
    
    CleanupDriverStuff();
    return 0;
}

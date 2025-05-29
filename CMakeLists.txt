cmake_minimum_required(VERSION 3.10)

project(silogen)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (MSVC)
    # msvc parallel build
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

    # vcpkg toolchain
    # set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")

    # cgal gmp/mpfr lib
    # set(GMP_INCLUDE_DIR "C:/dev/gmp/include" CACHE PATH "Path to GMP include directory")
    # set(GMP_LIBRARIES "C:/dev/gmp/lib/libgmp-10.lib" CACHE FILEPATH "Path to GMP library")
    # set(MPFR_INCLUDE_DIR "C:/dev/gmp/include" CACHE PATH "Path to GMP include directory")
    # set(MPFR_LIBRARIES "C:/dev/gmp/lib/libmpfr-4.lib" CACHE FILEPATH "Path to GMP library")
endif()


########################### libs for windows embedded inside


########################### libs
 
function(print_all_var)
    message(STATUS "********** begin print all var")
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
    message(STATUS "********** end  print all var")
endfunction()


if (MSVC)
    set(CMAKE_PREFIX_PATH  
        "E:/local/silo4.10.3" 
        "E:/local/hdf5/1.8.19" 
        ${CMAKE_PREFIX_PATH} 
    )
endif()


find_package(Silo REQUIRED)
message(STATUS "Silo_FOUND: ${Silo_FOUND}")
message(STATUS "Silo_VERSION: ${Silo_VERSION}")
message(STATUS "SILO_INCLUDE_DIR: ${SILO_INCLUDE_DIR}")
message(STATUS "SILO_LIBRARY: ${SILO_LIBRARY}")
find_library(SILO_LIBRARY NAMES silo siloh5 PATHS "${SILO_INCLUDE_DIR}/../lib")
# find_library(SILO_LIBRARY NAMES silo siloh5 PATHS "C:/Users/amcax/LLNL/VisIt3.4.1/lib")
message(STATUS "SILO_LIBRARY: ${SILO_LIBRARY}")

# message(STATUS "###############################################################################")
# message(STATUS "########################### after #############################################")
# message(STATUS "###############################################################################")
# print_all_var()

########################### src

file(GLOB GEN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)

source_group("gensilo" FILES ${GEN_SRC})

add_executable(${PROJECT_NAME}
    ${GEN_SRC} 
)

########################### include
target_include_directories(${PROJECT_NAME} PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${SILO_INCLUDE_DIR}
)

########################### link
target_link_libraries(${PROJECT_NAME} PRIVATE 
    ${SILO_LIBRARY} 
)

if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_ENVIRONMENT "PATH=E:/local/silo4.10.3/bin;%PATH%")
endif()
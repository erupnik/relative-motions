cmake_minimum_required(VERSION 2.8)
project(cov_in_motions)

set(BUILD_PATH_BIN  ${PROJECT_SOURCE_DIR}/bin CACHE PATH "binaires install path")
set(INSTALL_DIR ${BUILD_PATH_BIN})


find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

#find_package(CXSparse REQUIRED)
#include_directories(${CXSparse_INCLUDE_DIRS})

#find_package(SuiteSparse 4.0 COMPONENTS CHOLMOD SPQR)
#include_directories(${SuiteSparse_INCLUDE_DIRS})

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(gflags)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories("include")
include_directories("src")


file(GLOB_RECURSE SOURCE_FILES
    "src/main.cpp"
    "src/triplet_utils.cpp"
    "src/cov_in_motions.cpp"
    "src/bundle.cpp"
    "src/check_propagation.cpp"
)


add_executable(${PROJECT_NAME} "src/main.cpp" "src/triplet_utils.cpp" "src/cov_in_motions.cpp" "src/bundle.cpp" )
add_executable("triplet_utils" "src/triplet_utils_main.cpp" "src/triplet_utils.cpp" "src/cov_in_motions.cpp" "src/bundle.cpp")
#add_executable("bench" "src/check_propagation.cpp" "src/cov_in_motions.cpp" "src/bundle.cpp" )


target_link_libraries(${PROJECT_NAME} ${CERES_LIBRARIES} ${Boost_LIBRARIES} ) #
target_link_libraries("triplet_utils" ${CERES_LIBRARIES} ) #
#target_link_libraries("bench" ${CERES_LIBRARIES} ) #



install (TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_DIR})
install (TARGETS "triplet_utils" DESTINATION ${INSTALL_DIR})
#install (TARGETS "bench" DESTINATION ${INSTALL_DIR})


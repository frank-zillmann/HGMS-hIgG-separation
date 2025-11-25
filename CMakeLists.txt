cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME "HGMS_hIgG_separation")
project(${PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add FS3 installation to CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/external/FS3/install")

# Find FS3 package
find_package(FS3 REQUIRED)

# Add your own sources as a static library
add_library(${PROJECT_NAME} STATIC
    src/HGMS_hIgG_separation.cpp
)

# Include directories for your code
# (FS3 exports its own include dirs automatically)
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME} PUBLIC FS3::FS3)

# Executable for normal run
add_executable(HGMS_hIgG_separation_main_run
    src/main_run.cpp
)
target_link_libraries(HGMS_hIgG_separation_main_run PRIVATE ${PROJECT_NAME})

# Executable for hyperparameter search
add_executable(HGMS_hIgG_separation_hyperparameter_grid_search
    src/hyperparameter_grid_search.cpp
)
target_link_libraries(HGMS_hIgG_separation_hyperparameter_grid_search PRIVATE ${PROJECT_NAME})

# Installation
include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(TARGETS HGMS_hIgG_separation_main_run HGMS_hIgG_separation_hyperparameter_grid_search
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY src/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)


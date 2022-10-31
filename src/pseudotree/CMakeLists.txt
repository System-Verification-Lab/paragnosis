# CMake file for DAOOPT solver
# by Lars Otten, 2012

cmake_minimum_required(VERSION 2.6)
project(DAOOPT)

# To select the two master modes
option(DAOOPT_WORKER "Force worker binary" OFF)
option(DAOOPT_MASTER_STATIC "Force static master binary" OFF)
option(DAOOPT_MASTER_DYNAMIC "Force dynamic master binary" OFF)

# To enable static linking of the final daoopt binary
option(DAOOPT_LINK_STATIC "Link binary statically" OFF)

# General Compiler flags
add_definitions(-Wall)

if(WIN32)
  add_definitions(-DWINDOWS)
else()
  add_definitions(-DLINUX)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -DDEBUG -O0 -g3")
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -DDEBUG -O0 -g3")

set(CMAKE_CXX_FLAGS_RELEASE "-w -O3")
set(CMAKE_C_FLAGS_RELEASE "-w -O3")

if(DAOOPT_LINK_STATIC)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

# Optional forced version: static over dynamic master, worker over both
if(DAOOPT_WORKER)
  add_definitions(-DNOTHREADS)
else()
  if(DAOOPT_MASTER_STATIC)
    add_definitions(-DPARALLEL_STATIC)
  endif()
  if(DAOOPT_MASTER_DYNAMIC AND NOT DAOOPT_MASTER_STATIC)
    add_definitions(-DPARALLEL_DYNAMIC)
    find_package(Threads)  # includes pthread
    find_package( Boost REQUIRED COMPONENTS system )
  endif()
endif()

find_package( Boost REQUIRED COMPONENTS program_options thread )
include_directories( ${Boost_INCLUDE_DIRS} )

# Included folders for header files
include_directories(
  ${DAOOPT_SOURCE_DIR}/include
  ${DAOOPT_SOURCE_DIR}/lib
)

# Libraries in the lib folder via recursive CMakeLists.txt
add_subdirectory(lib)
set(LIBS Gzstream Sls4mpe Mex Cvo Zlib)

# Main daoopt source files
file(GLOB FILES "${CMAKE_CURRENT_LIST_DIR}/source/*.cpp")

# Main executable and library dependencies
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    get_filename_component(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
endif()
set(INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${INCLUDE_DIR} ${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

set(PROJECT "daoopt")

add_library(${PROJECT}-base STATIC ${FILES})

include_directories(${CMAKE_CURRENT_DIR})
set(DAOOPT_LIBRARY "libdaoopt.a")
add_custom_command(
    OUTPUT ${DAOOPT_LIBRARY}
    COMMAND echo create ${DAOOPT_LIBRARY}           >  build-lib.mri
    COMMAND echo addlib $<TARGET_FILE:Sls4mpe>      >> build-lib.mri
    COMMAND echo addlib $<TARGET_FILE:Mex>          >> build-lib.mri
    COMMAND echo addlib $<TARGET_FILE:Cvo>          >> build-lib.mri
    COMMAND echo addlib $<TARGET_FILE:daoopt-base>  >> build-lib.mri
    COMMAND echo addlib $<TARGET_FILE:Zlib>         >> build-lib.mri
    COMMAND echo addlib $<TARGET_FILE:Gzstream>     >> build-lib.mri
    COMMAND echo save                               >> build-lib.mri
    COMMAND echo end                                >> build-lib.mri
    COMMAND ar -M < build-lib.mri
    DEPENDS Gzstream Sls4mpe Mex Cvo ${PROJECT}-base
    COMMENT "Create ${DAOOPT_LIBRARY}"
)
add_custom_target(${PROJECT}-gen DEPENDS ${DAOOPT_LIBRARY})
add_library(${PROJECT}-lib STATIC IMPORTED)
set_property(TARGET ${PROJECT}-lib PROPERTY IMPORTED_LOCATION ${DAOOPT_LIBRARY})
add_dependencies(${PROJECT}-lib ${PROJECT}-gen)

add_executable(${PROJECT}-bin daoopt.cpp)
set_target_properties(${PROJECT}-bin PROPERTIES OUTPUT_NAME ${PROJECT})
target_link_libraries(${PROJECT}-bin ${PROJECT}-lib ${CMAKE_THREAD_LIBS_INIT} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY})

file(GLOB HEADERS "${INCLUDE_DIR}/*.h")
file(GLOB CVO_HEADERS "${CMAKE_CURRENT_LIST_DIR}/lib/cvo/*.hxx" "${CMAKE_CURRENT_LIST_DIR}/lib/cvo/*.h")
file(GLOB MEX_HEADERS "${CMAKE_CURRENT_LIST_DIR}/lib/mex/*.h")
file(GLOB SLS4MPE_HEADERS "${CMAKE_CURRENT_LIST_DIR}/lib/sls4mpe/*.h")
file(GLOB ZLIB_HEADERS "${CMAKE_CURRENT_LIST_DIR}/lib/zlib/*.h")
file(GLOB GZSTREAM_HEADERS "${CMAKE_CURRENT_LIST_DIR}/lib/*.h")

install(FILES $<TARGET_FILE:daoopt-lib> DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(FILES ${HEADERS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT}")
install(FILES ${CVO_HEADERS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT}/cvo")
install(FILES ${MEX_HEADERS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT}/mex")
install(FILES ${SLS4MPE_HEADERS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT}/sls4mpe")
install(FILES ${ZLIB_HEADERS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT}/zlib")
install(FILES ${GZSTREAM_HEADERS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT}")


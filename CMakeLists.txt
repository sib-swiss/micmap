#########################################################################################
# Header Message
#
MESSAGE(STATUS "+--------------------------------------------------------------------+")
MESSAGE(STATUS "|                              MICMAP                                |")
MESSAGE(STATUS "+--------------------------------------------------------------------+")
MESSAGE(STATUS "| Copyright (C) 2015-2019 SIB Swiss Institute of Bioinformatics      |")
MESSAGE(STATUS "| Copyright (C) 2019-2022 University of Lausanne, Switzerland        |")
MESSAGE(STATUS "| Copyright (C) 2020-2022 EPFL, Lausanne, Switzerland                |")
MESSAGE(STATUS "|         Nicolas Guex       (nicolas.guex@unil.ch)                  |")
MESSAGE(STATUS "|         Christian Iseli    (christian.iseli@epfl.ch)               |")
MESSAGE(STATUS "|         Thierry Schuepbach                                         |")
MESSAGE(STATUS "|                                                                    |")
MESSAGE(STATUS "+--------------------------------------------------------------------+")

#########################################################################################
# CMAKE
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.8 FATAL_ERROR)
PROJECT(GTL C CXX)
IF(COMMAND cmake_policy)
  IF(POLICY CMP0003)
		CMAKE_POLICY(SET CMP0003 NEW)
  ENDIF(POLICY CMP0003)
  IF(POLICY CMP0053)
		CMAKE_POLICY(SET CMP0053 NEW)
  ENDIF(POLICY CMP0053)
ENDIF(COMMAND cmake_policy)

LIST(APPEND CMAKE_MODULE_PATH ${GTL_SOURCE_DIR}/cmake)
enable_testing()

#########################################################################################
# Version with GIT
#
SET(MAJOR_VERSION "2")
SET(MINOR_VERSION "4")
SET(VERSION "${MAJOR_VERSION}.${MINOR_VERSION}")

FIND_PACKAGE(Git REQUIRED)
cmake_host_system_information(RESULT buildhost QUERY FQDN)
FILE(WRITE ${CMAKE_BINARY_DIR}/gtlVersion.h.in
"\#define GTL_MAJOR_VERSION ${MAJOR_VERSION}\n"
"\#define GTL_MINOR_VERSION ${MINOR_VERSION}\n"
"\#define GTL_PATCH_VERSION @PATCH_VERSION@\n"
"\#define GTL_VERSION \"${MAJOR_VERSION}.${MINOR_VERSION}-@PATCH_VERSION@\"\n"
"\#define GTL_VERSION_FULL_STRING \"GTL version ${MAJOR_VERSION}.${MINOR_VERSION}-@PATCH_VERSION@ built \" __DATE__ \", \" __TIME__ \" on ${buildhost}\"\n"
)
FILE(WRITE ${CMAKE_BINARY_DIR}/version.cmake
"EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE} rev-parse --verify --short HEAD
     OUTPUT_VARIABLE PATCH_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
 )
 CONFIGURE_FILE(\${SRC} \${DST} @ONLY)
")
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
ADD_CUSTOM_TARGET(
    git_version
    ${CMAKE_COMMAND} -D SRC=${CMAKE_BINARY_DIR}/gtlVersion.h.in
                     -D DST=${CMAKE_BINARY_DIR}/gtlVersion.h
                     -P ${CMAKE_BINARY_DIR}/version.cmake
)

EXECUTE_PROCESS(
     COMMAND ${GIT_EXECUTABLE} rev-parse --verify --short HEAD
     OUTPUT_VARIABLE PATCH_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
)

#########################################################################################
# Architectures
#
INCLUDE(CMakeDetermineSystem)
#INCLUDE(CMakePrintSystemInformation)
IF(CMAKE_SYSTEM_PROCESSOR)
    MESSAGE(STATUS "Compilation on architecture ${CMAKE_SYSTEM_PROCESSOR}.")
ENDIF(CMAKE_SYSTEM_PROCESSOR)

#########################################################################################
# What type of build
#
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#########################################################################################
# Determine C Compiler
#
IF(NOT CMAKE_COMPILER_IS_GNUCC)
    EXEC_PROGRAM(${CMAKE_C_COMPILER} ARGS ${CMAKE_BOOT_C_FLAGS} -E "\"${PFTOOLS_SOURCE_DIR}/cmake/CMakeTestIntel.c\"" OUTPUT_VARIABLE CMAKE_COMPILER_OUTPUT RETURN_VALUE CMAKE_COMPILER_RETURN)
    SET(CMAKE_COMPILER_IS_GNUCC_RUN 1)
    IF(NOT CMAKE_COMPILER_RETURN)
        IF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*THIS_IS_INTEL.*" )
             INCLUDE(${PFTOOLS_SOURCE_DIR}/cmake/Linux-icc.cmake)
        ENDIF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*THIS_IS_INTEL.*" )
    ENDIF(NOT CMAKE_COMPILER_RETURN)
ENDIF(NOT CMAKE_COMPILER_IS_GNUCC)

INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
CHECK_C_COMPILER_FLAG(-msse HANDLE_SSE)
IF(HANDLE_SSE)
  SET(CMAKE_C_SSE_FLAGS "-msse")
ELSE(HANDLE_SSE)
  MESSAGE(FATAL_ERROR "C compiler does not support SSE")
ENDIF(HANDLE_SSE)
CHECK_C_COMPILER_FLAG(-msse2 HANDLE_SSE2)
IF(HANDLE_SSE2)
  SET(CMAKE_C_SSE2_FLAGS "-msse2")
ELSE(HANDLE_SSE2)
  MESSAGE(FATAL_ERROR "C compiler does not support SSE 2")
ENDIF(HANDLE_SSE2)
SET(SSE41_FLAGS "-msse4.1" "-msse41")
UNSET(CMAKE_C_SSE41_FLAGS)
FOREACH(flag IN LISTS SSE41_FLAGS)
  MESSAGE(STATUS "testing flag ${flag}...")
  CHECK_C_COMPILER_FLAG("${flag}" HANDLE_SSE41)
  IF(HANDLE_SSE41)
    SET(CMAKE_C_SSE41_FLAGS "${flag}")
    BREAK()
  ENDIF(HANDLE_SSE41)
ENDFOREACH(flag)
IF(NOT DEFINED CMAKE_C_SSE41_FLAGS)
  MESSAGE(FATAL_ERROR "C compiler does not support SSE 4.1")
ENDIF(NOT DEFINED CMAKE_C_SSE41_FLAGS)
CHECK_C_COMPILER_FLAG(-mavx HANDLE_AVX)
IF(HANDLE_AVX)
  SET(CMAKE_C_AVX_FLAGS "-mavx")
ELSE(HANDLE_AVX)
  MESSAGE(FATAL_ERROR "C compiler does not support AVX")
ENDIF(HANDLE_AVX)

CHECK_C_COMPILER_FLAG(-mavx2 HANDLE_AVX2)
IF(HANDLE_AVX2)
  SET(CMAKE_C_AVX2_FLAGS "-mavx2")
ELSE(HANDLE_AVX2)
  MESSAGE(WARNING "C compiler does not support AVX2")
ENDIF(HANDLE_AVX2)

IF (${CMAKE_C_COMPILER} MATCHES "icc.*$")
	CHECK_C_COMPILER_FLAG(-xMIC-AVX512 HANDLE_AVX512F)
	IF(HANDLE_AVX512F)
		SET(CMAKE_C_AVX512F_FLAGS "-xMIC-AVX512")
		SET(CMAKE_REQUIRED_FLAGS ${CMAKE_C_AVX512F_FLAGS})
		SET(CMAKE_REQUIRED_INCLUDES "immintrin.h")
		CHECK_FUNCTION_EXISTS(_mm512_int2mask HAS_AVX512_INT2MASK)
	ELSE(HANDLE_AVX512F)
		MESSAGE(WARNING "C compiler does not support AVX512-F")
	ENDIF(HANDLE_AVX512F)

	CHECK_C_COMPILER_FLAG(-xCORE-AVX512 HANDLE_AVX512CD)
	IF(HANDLE_AVX512CD)
		SET(CMAKE_C_AVX512CD_FLAGS "-xCORE-AVX512")
		SET(CMAKE_C_AVX512BW_FLAGS "-xCORE-AVX512")
		SET(HANDLE_AVX512BW ON)
	ELSE(HANDLE_AVX512CD)
		MESSAGE(WARNING "C compiler does not support AVX512-CD")
		MESSAGE(WARNING "C compiler does not support AVX512-BW")
	ENDIF(HANDLE_AVX512CD)

ELSE (${CMAKE_C_COMPILER} MATCHES "icc.*$")
	CHECK_C_COMPILER_FLAG(-mavx512f HANDLE_AVX512F)
	IF(HANDLE_AVX512F)
		SET(CMAKE_C_AVX512F_FLAGS "-mavx512f")
		SET(CMAKE_REQUIRED_FLAGS ${CMAKE_C_AVX512F_FLAGS})
		#SET(CMAKE_REQUIRED_INCLUDES "immintrin.h")
		CHECK_SYMBOL_EXISTS(_mm512_int2mask "immintrin.h" HAS_AVX512_INT2MASK)
	ELSE(HANDLE_AVX512F)
		MESSAGE(WARNING "C compiler does not support AVX512-F")
	ENDIF(HANDLE_AVX512F)

	CHECK_C_COMPILER_FLAG(-mavx512cd HANDLE_AVX512CD)
	IF(HANDLE_AVX512CD)
		SET(CMAKE_C_AVX512CD_FLAGS "-mavx512cd")
	ELSE(HANDLE_AVX512CD)
		MESSAGE(WARNING "C compiler does not support AVX512-CD")
	ENDIF(HANDLE_AVX512CD)

	CHECK_C_COMPILER_FLAG(-mavx512bw HANDLE_AVX512BW)
	IF(HANDLE_AVX512BW)
		SET(CMAKE_C_AVX512BW_FLAGS "-mavx512bw")
	ELSE(HANDLE_AVX512BW)
		MESSAGE(WARNING "C compiler does not support AVX512-BW")
	ENDIF(HANDLE_AVX512BW)
ENDIF (${CMAKE_C_COMPILER} MATCHES "icc.*$")

CHECK_C_COMPILER_FLAG("-std=c99" HANDLE_C99)
IF(NOT HANDLE_C99)
  MESSAGE(FATAL_ERROR "GTL require a C compiler with C99 support")
ELSE(NOT HANDLE_C99)
  IF(CMAKE_C_FLAGS)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
  ELSE(CMAKE_C_FLAGS)
		SET(CMAKE_C_FLAGS "-std=c99")
  ENDIF(CMAKE_C_FLAGS)
  MESSAGE(STATUS "Add -std=c99 to C compiler options")
ENDIF(NOT HANDLE_C99)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_SSE2_FLAGS}")
MESSAGE(STATUS "Add SSE2 to C compiler options")

IF(CMAKE_COMPILER_IS_GNUCC)
  CHECK_C_COMPILER_FLAG("-fms-extensions" MS_EXTENSION)
  IF(MS_EXTENSION)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fms-extensions")
  ENDIF(MS_EXTENSION)
ENDIF(CMAKE_COMPILER_IS_GNUCC)

if (${CMAKE_C_COMPILER} MATCHES "icc.*$")
  GET_FILENAME_COMPONENT(ICC_PATH ${CMAKE_C_COMPILER} DIRECTORY)
  GET_FILENAME_COMPONENT(INTEL_LIB_DIR ${ICC_PATH}/../lib/intel64 ABSOLUTE)
  SET(CMAKE_EXE_LINKER_FLAGS "-L${INTEL_LIB_DIR} ${CMAKE_EXE_LINKER_FLAGS}")
  MESSAGE(STATUS "Intel Compiler at use, adding library path ${INTEL_LIB_DIR}")
endif ()

#########################################################################################
# External Dependencies
#
SET(FIND_LIBRARY_USE_LIB64_PATHS 1)

# Check for pthread
FIND_PACKAGE(Threads REQUIRED)
IF( ${CMAKE_USE_WIN32_THREADS_INIT})
  MESSAGE(STATUS "Using Windows API for threads")
  SET(USE_WINAPI 1)
ENDIF(${CMAKE_USE_WIN32_THREADS_INIT})

#########################################################################################
# Third party directory
#
SET(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third-party)
SET(EXTERNAL_INSTALL_LOCATION ${PROJECT_BINARY_DIR}/external)

#########################################################################################
# Potentially Internal Dependencies
#
OPTION(STANDALONE "This will compile executable in static mode." OFF)

#########################################################################################
# EXTRA HEADER FILE
#
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (emmintrin.h HAVE_SSE2_H)
CHECK_INCLUDE_FILES (smmintrin.h HAVE_SSE_41_H)
IF (NOT HAVE_SSE2_H OR NOT HAVE_SSE_41_H)
#   MESSAGE(FATAL_ERROR "PFTOOLS require a compiler with both SSE 2 and 4.1 support")
ENDIF(NOT HAVE_SSE2_H OR NOT HAVE_SSE_41_H)

CHECK_INCLUDE_FILES (mm_malloc.h HAVE_MM_MALLOC_H)
CHECK_INCLUDE_FILES (alloca.h HAVE_ALLOCA_H)

SET(MMAP ON)

#########################################################################################
# Compiling options
#
# Use RPATH?
OPTION(USE_RPATH "Use -rpath when linking libraries, executables" ON)
IF(USE_RPATH)
	SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags")
	# use, i.e. don't skip the full RPATH for the build tree
	SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

	# when building, don't use the install RPATH already
	# (but later on when installing)
	SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

	# add the automatically determined parts of the RPATH
	# which point to directories outside the build tree to the install RPATH
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

	# the RPATH to be used when installing, but only if it's not a system directory
	LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
	IF("${isSystemDir}" STREQUAL "-1")
		SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:/opt/intel/compilers_and_libraries/linux/lib/intel64")
	ENDIF("${isSystemDir}" STREQUAL "-1")
ENDIF(USE_RPATH)

OPTION(USE_AFFINITY "CPU core affinity will be used for threads." ON)

OPTION(BUILD_SHARED_LIBS "Builds shared development library." OFF)
OPTION(BUILD_STATIC_LIBS "Builds static development library." OFF)

#########################################################################################
# Enabling TESTING
include(CTest)
IF(BUILD_TESTING)
#	SET(GTL_TEST_FILES_PATH "/home/tschuepb/Data/Caller/HG001/" CACHE PATH "Set path to GTL testing files")
#	SET(GTL_TEST_GENOME_PATH "/home/tschuepb/Data/Genomes/" CACHE PATH "Set path to Genome testing files")
	SET(GTL_TESTS_PATH ${CMAKE_BINARY_DIR}/Tests CACHE PATH "Set output path for tests")
	IF (NOT EXISTS "${GTL_TESTS_PATH}" )
		FILE(MAKE_DIRECTORY "${GTL_TESTS_PATH}" )
	ENDIF (NOT EXISTS "${GTL_TESTS_PATH}" )

	IF (EXISTS "${GTL_TEST_FILES_PATH}")
		MESSAGE(STATUS "GTL test files in ${GTL_TEST_FILES_PATH}")
		IF (EXISTS "${GTL_TEST_GENOME_PATH}")
			MESSAGE(STATUS "Genome test files in ${GTL_TEST_GENOME_PATH}")
		ELSE(EXISTS "${GTL_TEST_GENOME_PATH}")
			MESSAGE(STATUS "Genome test files not found")
		ENDIF(EXISTS "${GTL_TEST_GENOME_PATH}")
	ELSE(EXISTS "${GTL_TEST_FILES_PATH}")
		MESSAGE(STATUS "GTL test files not found")
	ENDIF(EXISTS "${GTL_TEST_FILES_PATH}")
ENDIF(BUILD_TESTING)

#########################################################################################
# CONFIGURATION FILE
#
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
ADD_DEFINITIONS("-I ${CMAKE_CURRENT_BINARY_DIR}/")

#########################################################################################
# STATIC BUILD
#
IF(STANDALONE)

  #######################################################################################
  ## CHECK TO SEE IF PCRE HAS STATIC LIBRARY
  FIND_LIBRARY(ZLIB_STATIC libz.a)
  IF(NOT ZLIB_STATIC)
    MESSAGE(FATAL_ERROR "Standalone static building is not possible without static library zlib\n")
  ENDIF(NOT ZLIB_STATIC)

  ADD_DEFINITIONS("-DNO_USERNAME")
  SET(USE_RPATH OFF)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
#  SET (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")

	SET(CMAKE_SKIP_INSTALL_RPATH ON)
	SET(CMAKE_SKIP_RPATH ON)
	SET(BUILD_STATIC_LIBS ON)

ENDIF(STANDALONE)

#########################################################################################
## PATHs
INCLUDE_DIRECTORIES(${GTL_SOURCE_DIR}/include;${CMAKE_BINARY_DIR})

#########################################################################################
# Now treats subdirectories
#
ADD_SUBDIRECTORY(${GTL_SOURCE_DIR}/utils)
ADD_SUBDIRECTORY(${GTL_SOURCE_DIR}/src)

#########################################################################################
# Enabling TESTING
include(CTest)
IF(BUILD_TESTING)
#       SET(GTL_TEST_FILES_PATH "/home/tschuepb/Data/Caller/HG001/" CACHE PATH "Set path to GTL testing files")
#       SET(GTL_TEST_GENOME_PATH "/home/tschuepb/Data/Genomes/" CACHE PATH "Set path to Genome testing files")
        SET(GTL_TESTS_PATH ${CMAKE_BINARY_DIR}/Tests CACHE PATH "Set output path for tests")
        IF (NOT EXISTS "${GTL_TESTS_PATH}" )
                FILE(MAKE_DIRECTORY "${GTL_TESTS_PATH}" )
        ENDIF (NOT EXISTS "${GTL_TESTS_PATH}" )
	INCLUDE(cmake/tests.cmake)
ENDIF(BUILD_TESTING)

#########################################################################################
## Data files
# ADD_SUBDIRECTORY(data)

#########################################################################################
## Documentation
#ADD_SUBDIRECTORY(doc)

#########################################################################################
## Install script
ADD_SUBDIRECTORY(scripts/)

#########################################################################################
## Packaging
include(cmake/Packaging.cmake)

#########################################################################################
## Output global build information
include(summary)
summary()

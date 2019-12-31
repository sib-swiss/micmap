
# htslib
#FIND_PACKAGE(htslib QUIET)
IF (NOT htslib_FOUND)
	find_package(ZLIB REQUIRED)
	ExternalProject_Add(
			htslib
			PREFIX  ${EXTERNAL_INSTALL_LOCATION}/htslib
			GIT_REPOSITORY https://gitlab.isb-sib.ch/tschuepb/htslib.git
			TIMEOUT 10
			UPDATE_COMMAND ${GIT_EXECUTABLE} pull
			CONFIGURE_COMMAND autoreconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl --disable-s3 --disable-gcs
			BUILD_COMMAND make -j 4
			BUILD_IN_SOURCE 1
			INSTALL_COMMAND ""
			LOG_DOWNLOAD OFF
	)
	ExternalProject_Get_Property(htslib source_dir)
	MESSAGE(STATUS "htslib downloaded into ${source_dir}")
	SET(HTSLIB_DIR "${source_dir}" CACHE PATH "Path to HTS library directory")
	SET(HTSLIB_INCLUDE_DIRS "${source_dir}" CACHE PATH "Path to HTS include directory")
	
	IF(NOT STANDALONE)
		SET(HTSLIB_LIBRARIES "${source_dir}/libhts.so;-lz;-lpthread" CACHE STRING "HTS libraries")
		
		INSTALL(FILES "${source_dir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}" "${source_dir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}.2"
						DESTINATION lib
						PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
												GROUP_EXECUTE GROUP_READ 
												WORLD_EXECUTE WORLD_READ
		)
 	ELSE(NOT STANDALONE)
		SET(HTSLIB_LIBRARIES "${source_dir}/libhts.a;-lz;-lpthread" CACHE STRING "HTS libraries")
 	ENDIF(NOT STANDALONE)
 	mark_as_advanced(HTSLIB_DIR HTSLIB_INCLUDE_DIRS HTSLIB_LIBRARIES)
ENDIF(NOT htslib_FOUND)

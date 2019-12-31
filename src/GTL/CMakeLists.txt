
INCLUDE_DIRECTORIES(Compresser)
SET(GTL_COMPRESSER_SRCS
		Compresser/lz4.c
		Compresser/lz4frame.c
		Compresser/lz4hc.c
		Compresser/rangecod.c
		Compresser/xxhash.c
)

SET(GTL_GENOME_SRCS Genome/virt_chr.c Genome/genome.c)

ADD_LIBRARY(gtl_shared SHARED
						GTL.c
						GTLdispatch.c
						threadpool.c
						${GTL_COMPRESSER_SRCS}
						${GTL_GENOME_SRCS}
)
set_target_properties( gtl_shared PROPERTIES LINKER_LANGUAGE C LIBRARY_OUTPUT_NAME "gtl")
target_link_libraries(gtl_shared z;bz2;pthread)
target_include_directories(gtl_shared PUBLIC ${PROJECT_SOURCE_DIR}/include)
ADD_DEPENDENCIES(gtl_shared git_version)
INSTALL(TARGETS gtl_shared LIBRARY DESTINATION lib)

IF (BUILD_STATIC_LIBS)
	ADD_LIBRARY(gtl_static STATIC
							GTL.c
							GTLdispatch.c
							threadpool.c
							${GTL_COMPRESSER_SRCS}
							${GTL_GENOME_SRCS}
	)
	set_target_properties( gtl_static PROPERTIES LINKER_LANGUAGE C ARCHIVE_OUTPUT_NAME "gtl")
	target_link_libraries(gtl_static z;bz2;pthread)
	target_include_directories(gtl_static PUBLIC ${PROJECT_SOURCE_DIR}/include)
	ADD_DEPENDENCIES(gtl_static git_version)
	INSTALL(TARGETS gtl_static ARCHIVE DESTINATION lib)
ENDIF (BUILD_STATIC_LIBS)

ADD_SUBDIRECTORY(tools)

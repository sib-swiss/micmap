OPTION(BUILD_MPI "Build MPI softwares" OFF)
FIND_PACKAGE(MPI)
IF(MPI_FOUND AND BUILD_MPI)
	SET(ALIGNERS $<TARGET_OBJECTS:sse41_Aligner>)
	IF(HANDLE_AVX512F)
		LIST(APPEND ALIGNERS $<TARGET_OBJECTS:avx512_Aligner>)
	ENDIF(HANDLE_AVX512F)

	ADD_LIBRARY(MPILink OBJECT mpi_transfer.c)
	target_include_directories(MPILink PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../;${MPI_C_HEADER_DIR})
	
	ADD_EXECUTABLE(MpiMap ../main.c
                       $<TARGET_OBJECTS:Reader>
                       $<TARGET_OBJECTS:Decoder>
                       $<TARGET_OBJECTS:Writer>
                       $<TARGET_OBJECTS:Aligner>
	               ${ALIGNERS}
	               ${TOPOLOGY}
                       ../system.c
                       $<TARGET_OBJECTS:Statistics>
                       $<TARGET_OBJECTS:MPILink>)
                       
	TARGET_LINK_LIBRARIES(MpiMap -lpthread;gtl_shared;dl;${MPI_C_LIBRARIES})
	target_compile_definitions(MpiMap PRIVATE USE_MPI)
	target_include_directories(MpiMap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../;${MPI_C_HEADER_DIR})
	
ENDIF(MPI_FOUND AND BUILD_MPI)

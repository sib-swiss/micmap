ADD_DEFINITIONS("-DUSE_AFFINITY")
ADD_LIBRARY(Reader OBJECT readFASTQ.c)
IF(HANDLE_AVX2)
	ADD_LIBRARY(Decoder OBJECT decode_common.c decode.c decode_avx2.c)
	SET_PROPERTY(SOURCE decode_avx2.c APPEND PROPERTY COMPILE_FLAGS ${CMAKE_C_AVX2_FLAGS})
ELSEIF(HANDLE_AVX2)
	ADD_LIBRARY(Decoder OBJECT decode_common.c decode.c)
ENDIF(HANDLE_AVX2)

IF (HUGETLB_FS)
	ADD_EXECUTABLE(LoadDecodingTable LoadDecodingTable.c decode_common.c)
	TARGET_LINK_LIBRARIES(LoadDecodingTable rt;pthread)
	ADD_EXECUTABLE(ReadDecodingTable ReadDecode.c)
	TARGET_LINK_LIBRARIES(ReadDecodingTable rt;hugetlbfs)
	INSTALL(TARGETS LoadDecodingTable ReadDecodingTable DESTINATION bin) 
ENDIF(HUGETLB_FS)

ADD_LIBRARY(Writer OBJECT write.c)
ADD_LIBRARY(System OBJECT system.c system.h)
ADD_EXECUTABLE(SystemInfo system.c)
target_compile_definitions(SystemInfo PRIVATE _SYSTEM_TEST USE_AFFINITY )
target_link_libraries(SystemInfo dl)

ADD_LIBRARY(Statistics OBJECT statistics.c statistics.h)
ADD_LIBRARY(Aligner 
					  OBJECT 
							align_common.c
							align_cpu.c
)

ADD_SUBDIRECTORY(SSE4.1/)
SET(SPECIFIC_ALIGNERS $<TARGET_OBJECTS:sse41_Aligner>)

IF (HANDLE_AVX2)
 	ADD_SUBDIRECTORY(AVX2)
ENDIF(HANDLE_AVX2)

IF (HANDLE_AVX512F)
	ADD_SUBDIRECTORY(AVX512F)
	LIST(APPEND SPECIFIC_ALIGNERS "$<TARGET_OBJECTS:avx512_Aligner>")
ENDIF(HANDLE_AVX512F)

SET(TOPOLOGY "")
IF(USE_AFFINITY)
	ADD_LIBRARY(Topology OBJECT 
							topology.c
							${THIRD_PARTY_DIR}/json-c/arraylist.c
							${THIRD_PARTY_DIR}/json-c/debug.c
							${THIRD_PARTY_DIR}/json-c/json_object.c
							${THIRD_PARTY_DIR}/json-c/json_pointer.c
							${THIRD_PARTY_DIR}/json-c/json_tokener.c
							${THIRD_PARTY_DIR}/json-c/json_util.c
							${THIRD_PARTY_DIR}/json-c/linkhash.c
							${THIRD_PARTY_DIR}/json-c/printbuf.c
							${THIRD_PARTY_DIR}/json-c/random_seed.c
							${THIRD_PARTY_DIR}/json-c/strerror_override.c)
	SET_PROPERTY(SOURCE ${THIRD_PARTY_DIR}/json-c/arraylist.c ${THIRD_PARTY_DIR}/json-c/debug.c
							${THIRD_PARTY_DIR}/json-c/json_object.c ${THIRD_PARTY_DIR}/json-c/json_pointer.c
							${THIRD_PARTY_DIR}/json-c/json_tokener.c ${THIRD_PARTY_DIR}/json-c/json_util.c
							${THIRD_PARTY_DIR}/json-c/linkhash.c ${THIRD_PARTY_DIR}/json-c/printbuf.c
							${THIRD_PARTY_DIR}/json-c/random_seed.c ${THIRD_PARTY_DIR}/json-c/strerror_override.c
							APPEND PROPERTY COMPILE_DEFINITIONS _GNU_SOURCE )
	TARGET_INCLUDE_DIRECTORIES(Topology PRIVATE ${THIRD_PARTY_DIR}/json-c)
	SET(TOPOLOGY $<TARGET_OBJECTS:Topology>)
ENDIF(USE_AFFINITY)

IF(CMAKE_C_FLAGS_<CONFIG> MATCHES "-DNDEBUG")
	SET(STREAM_FOR_DEBUG "")
ELSE(CMAKE_C_FLAGS_<CONFIG> MATCHES "-DNDEBUG")
	SET(STREAM_FOR_DEBUG "stream.c")
ENDIF(CMAKE_C_FLAGS_<CONFIG> MATCHES "-DNDEBUG")
ADD_EXECUTABLE(MicMap main.c 
              $<TARGET_OBJECTS:Reader>
              $<TARGET_OBJECTS:Decoder>
              $<TARGET_OBJECTS:Writer>
              $<TARGET_OBJECTS:Aligner>
              ${SPECIFIC_ALIGNERS}
              $<TARGET_OBJECTS:System>
              $<TARGET_OBJECTS:Statistics>
              ${TOPOLOGY}
              ${STREAM_FOR_DEBUG}
)


TARGET_LINK_LIBRARIES(MicMap -lpthread;gtl_shared;dl;rt)

INSTALL(TARGETS MicMap SystemInfo DESTINATION bin)

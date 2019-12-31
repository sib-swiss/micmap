add_library(FastAlignLib OBJECT
            fastalign_common.c
            fastalign_cpu.c
            )
TARGET_INCLUDE_DIRECTORIES(FastAlignLib PUBLIC "${THIRD_PARTY_DIR}/json-c;${CMAKE_CURRENT_SOURCE_DIR}")
# SET_PROPERTY(SOURCE fastalign_common.c 
# 						 APPEND PROPERTY COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/../")

SET_PROPERTY(SOURCE fastalign_cpu.c APPEND PROPERTY COMPILE_FLAGS ${CMAKE_C_SSE41_FLAGS})
add_executable(FastAlign
							 FastAlignTest.c
							 $<TARGET_OBJECTS:FastAlignLib>
 							 $<TARGET_OBJECTS:json-c>
)

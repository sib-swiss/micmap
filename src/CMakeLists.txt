###################################################################################
## THIRD PARTY PACKAGES 
ADD_LIBRARY(json-c OBJECT
    ${THIRD_PARTY_DIR}/json-c/arraylist.c
    ${THIRD_PARTY_DIR}/json-c/debug.c
    ${THIRD_PARTY_DIR}/json-c/json_object.c
    ${THIRD_PARTY_DIR}/json-c/json_pointer.c
    ${THIRD_PARTY_DIR}/json-c/json_tokener.c
    ${THIRD_PARTY_DIR}/json-c/json_util.c
    ${THIRD_PARTY_DIR}/json-c/linkhash.c
    ${THIRD_PARTY_DIR}/json-c/printbuf.c
    ${THIRD_PARTY_DIR}/json-c/random_seed.c
    ${THIRD_PARTY_DIR}/json-c/strerror_override.c
)
SET_PROPERTY(TARGET json-c APPEND PROPERTY COMPILE_DEFINITIONS _GNU_SOURCE)
TARGET_INCLUDE_DIRECTORIES(json-c PUBLIC ${THIRD_PARTY_DIR}/json-c)

###################################################################################
## SOURCES 
ADD_SUBDIRECTORY(FastAlignLib)
ADD_SUBDIRECTORY(GTL)
ADD_SUBDIRECTORY(MICMAP)
ADD_SUBDIRECTORY(Callers)

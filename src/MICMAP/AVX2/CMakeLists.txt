ADD_LIBRARY(avx2align OBJECT fastalign.c)
SET_PROPERTY(SOURCE fastalign.c PROPERTY COMPILE_FLAGS "${CMAKE_C_AVX2_FLAGS}")

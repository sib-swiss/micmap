ADD_EXECUTABLE(mkblob lzBlob.c)
TARGET_LINK_LIBRARIES(mkblob z)

ADD_EXECUTABLE(BedSplitter bedSplitter.c)

ADD_EXECUTABLE(VCFCorrector VCFCorrector.c)

ADD_EXECUTABLE(TestAlignerQueue testAlignerQueue.c)
SET_PROPERTY(SOURCE TestAlignerQueue.c APPEND PROPERTY COMPILE_FLAGS "-std=c11")
TARGET_LINK_LIBRARIES(TestAlignerQueue pthread)

INSTALL(TARGETS BedSplitter VCFCorrector
        DESTINATION bin )

IF(STANDALONE)
	SET(GTL_LIB gtl_static)
ELSE(STANDALONE)
	SET(GTL_LIB gtl_shared)
ENDIF(STANDALONE)

ADD_EXECUTABLE(GTLlist GTLlist.c)
TARGET_LINK_LIBRARIES(GTLlist ${GTL_LIB})

ADD_EXECUTABLE(GTLinfo GTLinfo.c)
TARGET_LINK_LIBRARIES(GTLinfo ${GTL_LIB})

ADD_EXECUTABLE(GTLfetch GTLfetch.c)
TARGET_LINK_LIBRARIES(GTLfetch ${GTL_LIB})

ADD_EXECUTABLE(GTLcnt GTLcnt.c)
TARGET_LINK_LIBRARIES(GTLcnt ${GTL_LIB})

ADD_EXECUTABLE(GTLcoverage GTLcoverage.c)
TARGET_LINK_LIBRARIES(GTLcoverage ${GTL_LIB})

ADD_EXECUTABLE(GTLsubCoverage GTLsubCoverage.c)
TARGET_LINK_LIBRARIES(GTLsubCoverage ${GTL_LIB})

ADD_EXECUTABLE(GTLindex GTLindex.c)
TARGET_LINK_LIBRARIES(GTLindex ${GTL_LIB})

ADD_EXECUTABLE(GTLmerge GTLmerge.c)
TARGET_LINK_LIBRARIES(GTLmerge ${GTL_LIB})

ADD_EXECUTABLE(GTLrecompress GTLrecompress.c)
TARGET_LINK_LIBRARIES(GTLrecompress ${GTL_LIB})

ADD_EXECUTABLE(GTLsimpleSort GTLsimpleSort.c)
TARGET_LINK_LIBRARIES(GTLsimpleSort ${GTL_LIB})

ADD_EXECUTABLE(GTLamplicon GTLamplicon.c)
TARGET_LINK_LIBRARIES(GTLamplicon ${GTL_LIB})

ADD_EXECUTABLE(GTLencodegenome encode_genome.c)
TARGET_LINK_LIBRARIES(GTLencodegenome ${GTL_LIB})

ADD_EXECUTABLE(GTLencodeUEM encode_UEM.c)
TARGET_LINK_LIBRARIES(GTLencodeUEM ${GTL_LIB})

ADD_EXECUTABLE(GTLdecompress GTLdecompress.c $<TARGET_OBJECTS:FastAlignLib> $<TARGET_OBJECTS:json-c>)
TARGET_INCLUDE_DIRECTORIES(GTLdecompress PRIVATE $<TARGET_PROPERTY:FastAlignLib,INTERFACE_INCLUDE_DIRECTORIES>)  
TARGET_LINK_LIBRARIES(GTLdecompress ${GTL_LIB})

ADD_EXECUTABLE(GTLregionAnalysis  GTLregionAnalysis.c)
TARGET_LINK_LIBRARIES(GTLregionAnalysis ${GTL_LIB})
#TARGET_INCLUDE_DIRECTORIES(GTLregionAnalysis PUBLIC "${HTSLIB_INCLUDE_DIRS}")

INSTALL(TARGETS GTLamplicon
                GTLcnt
                GTLcoverage
                GTLsubCoverage
                GTLdecompress
                GTLencodeUEM
                GTLencodegenome
                GTLindex
                GTLinfo
                GTLlist
                GTLfetch
                GTLmerge
                GTLrecompress
                GTLsimpleSort
                GTLregionAnalysis
        DESTINATION bin)

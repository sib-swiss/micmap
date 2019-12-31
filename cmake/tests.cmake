IF(BUILD_TESTING)
	OPTION(USE_HUGEPAGE_TLB_FS "Do use Huge Page TLB fs" ON)
	FIND_PATH(1GB_TLB_FS_PATH 
	          NAMES pagesize-1GB
	          HINTS "/var/lib/hugetlbfs" "/var/lib/hugetlbfs/global"
	          NO_DEFAULT_PATH)
	IF (1GB_TLB_FS_PATH)
		IF(NOT 1GB_TLB_FS_PATH MATCHES ".*pagesize-1GB")
			SET(1GB_TLB_FS_PATH "${1GB_TLB_FS_PATH}/pagesize-1GB" CACHE PATH "Huge TLB 1Gb fs path" FORCE)
		ENDIF(NOT 1GB_TLB_FS_PATH MATCHES ".*pagesize-1GB")
		MESSAGE(STATUS "Huge page TLB File System in ${1GB_TLB_FS_PATH}")
	ELSE(1GB_TLB_FS_PATH)
		MESSAGE(STATUS "Huge page TLB File System not found")
		SET(USE_HUGEPAGE_TLB_FS OFF CACHE PATH "Do use Huge Page TLB fs" FORCE)
	ENDIF(1GB_TLB_FS_PATH)
	
	IF (UNIX)
		EXECUTE_PROCESS(COMMAND grep -l '^flags[[:space:]]*:.*avx512f' /proc/cpuinfo OUTPUT_VARIABLE CPU_HAS_AVX512F)
	ELSE(UNIX)
		SET(CPU_HAS_AVX512F FALSE)
	ENDIF(UNIX)
	
	IF(NOT CPU_HAS_AVX512F)
		MESSAGE(STATUS "AVX512F tests will not be performed")
	ENDIF(NOT CPU_HAS_AVX512F)
	
	INCLUDE (FindUnixCommands)
	# ---------------------------------------------- MACROS -------------------------------------------------------------	
	MACRO(TestGenomeTLBLoading Genome_file Table_prefix Ref)
		SET(SUFFIXES A.2nt.bin;A.chr.bin;A.strand.bin;A.valid.bin)
		GET_FILENAME_COMPONENT(Genome ${Genome_file} NAME)
		GET_FILENAME_COMPONENT(table ${Table_prefix} NAME)
		SET(GenomeFiles "${Genome}")
		FOREACH(extension IN ITEMS ${SUFFIXES})
			SET(GenomeFiles "${GenomeFiles},${table}${extension}")
		ENDFOREACH(extension IN ITEMS ${SUFFIXES})
		ADD_TEST(NAME "LoadingTLBFS_${Ref}"
		         COMMAND ${BASH} -c "rm -f ${1GB_TLB_FS_PATH}/*; $<TARGET_FILE:MicMap> -g ${Genome_file} -p ${Table_prefix} -s ${1GB_TLB_FS_PATH}")
		SET_TESTS_PROPERTIES("LoadingTLBFS_${Ref}" PROPERTIES TIMEOUT 480)
	ENDMACRO(TestGenomeTLBLoading)
	
	MACRO(DecompressExistingGTL gtlFileName FASTQFileName ChrNumber Arguments)
		STRING(REPLACE ".gtl" "" CHRName "${gtlFileName}")
		ADD_TEST(NAME "Decompress${CHRName}"
		         COMMAND ${BASH} -c "($<TARGET_FILE:GTLdecompress> -g ${GENOME_B37_FILE} -p ${GENOME_B37_TABLE_PREFIX} -C ${ChrNumber} -r ${GTL_TEST_FILES_PATH}/${gtlFileName} ${Arguments} -d | sort -n | sed 's/^[0-9][0-9]* //'>${GTL_TESTS_PATH}/${FASTQFileName}_1.fq) 2>&1 | sort -n | sed 's/^[0-9][0-9]* //'>${GTL_TESTS_PATH}/${FASTQFileName}_2.fq" )
	ENDMACRO(DecompressExistingGTL)

	MACRO(TestMapping ExecutableName Type FASTQBaseFileName Arguments)
		SET(RESULTING_FILES "${GTL_TESTS_PATH}/${Type}_${ExecutableName}_${FASTQBaseFileName}")
		IF(USE_HUGEPAGE_TLB_FS)
			SET(GENOME_FILE "${1GB_TLB_FS_PATH}/b37.bin")
			SET(TABLE_FILE "${1GB_TLB_FS_PATH}/start0.")
		ELSE(USE_HUGEPAGE_TLB_FS)
			SET(GENOME_FILE ${GENOME_B37_FILE})
			SET(TABLE_FILE ${GENOME_B37_TABLE_PREFIX})
		ENDIF(USE_HUGEPAGE_TLB_FS)
		ADD_TEST(NAME ${Type}_${ExecutableName}_${FASTQBaseFileName}
		         COMMAND ${BASH} -c "rm -f ${RESULTING_FILES}*.gtl && $<TARGET_FILE:${ExecutableName}> -g ${GENOME_FILE} -p ${TABLE_FILE} -r ${RESULTING_FILES} ${Arguments} -1 ${GTL_TESTS_PATH}/${FASTQBaseFileName}_1.fq -2 ${GTL_TESTS_PATH}/${FASTQBaseFileName}_2.fq" )
		IF(USE_HUGEPAGE_TLB_FS)
			SET_TESTS_PROPERTIES(${Type}_${ExecutableName}_${FASTQBaseFileName} PROPERTIES DEPENDS LoadingTLBFS_b37)
		ENDIF(USE_HUGEPAGE_TLB_FS)
				
		ADD_TEST(NAME ${Type}_${ExecutableName}_${FASTQBaseFileName}_Decompress
		         COMMAND ${BASH} -c "($<TARGET_FILE:GTLdecompress> -g ${GENOME_B37_FILE} -r <(cat ${RESULTING_FILES}*.gtl) -p -n -m -a -u -h -c -d | sort -n | sed 's/^[0-9][0-9]* //'>${RESULTING_FILES}_1.fq) 2>&1 | sort -n | sed 's/^[0-9][0-9]* //'>${RESULTING_FILES}_2.fq" )
		IF(USE_HUGEPAGE_TLB_FS)
			SET_TESTS_PROPERTIES(${Type}_${ExecutableName}_${FASTQBaseFileName} PROPERTIES DEPENDS LoadingTLBFS_b37)
		ENDIF(USE_HUGEPAGE_TLB_FS)
		ADD_TEST(NAME ${Type}_${ExecutableName}_${FASTQBaseFileName}_RAW_Decompress
		         COMMAND ${BASH} -c "$<TARGET_FILE:GTLdecompress> -g ${GENOME_B37_FILE} -r <(cat ${RESULTING_FILES}*.gtl) -a -p -c -h -u -n -m -o raw -d | sort -n >${RESULTING_FILES}.raw")
		IF(USE_HUGEPAGE_TLB_FS)
			SET_TESTS_PROPERTIES(${Type}_${ExecutableName}_${FASTQBaseFileName} PROPERTIES DEPENDS LoadingTLBFS_b37)
		ENDIF(USE_HUGEPAGE_TLB_FS)
		ADD_TEST(NAME ${Type}_${ExecutableName}_${FASTQBaseFileName}_CmpPair1
		         COMMAND ${CMAKE_COMMAND} -E compare_files ${RESULTING_FILES}_1.fq ${GTL_TESTS_PATH}/${FASTQBaseFileName}_1.fq)
		ADD_TEST(NAME ${Type}_${ExecutableName}_${FASTQBaseFileName}_CmpPair2
		         COMMAND ${CMAKE_COMMAND} -E compare_files ${RESULTING_FILES}_2.fq ${GTL_TESTS_PATH}/${FASTQBaseFileName}_2.fq)
		IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
			SET_TESTS_PROPERTIES(${Type}_${ExecutableName}_${FASTQBaseFileName} PROPERTIES TIMEOUT 480)
		ELSE(CMAKE_BUILD_TYPE STREQUAL "Debug")
			SET_TESTS_PROPERTIES(${Type}_${ExecutableName}_${FASTQBaseFileName} PROPERTIES TIMEOUT 240)
		ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
	ENDMACRO(TestMapping)
	
	MACRO(SimpleTest FASTQBaseFileName Genome Table Ref)
		IF(USE_HUGEPAGE_TLB_FS)
			GET_FILENAME_COMPONENT(Genome_Name ${Genome} NAME)
			GET_FILENAME_COMPONENT(Table_Name ${Table} NAME)
			SET(GenomeFile "${1GB_TLB_FS_PATH}/${Genome_Name}")
			SET(TablePrefix "${1GB_TLB_FS_PATH}/${Table_Name} ")
		ELSE(USE_HUGEPAGE_TLB_FS)
			SET(GenomeFile "${Genome}")
                        SET(TablePrefix "${Table} ")
		ENDIF(USE_HUGEPAGE_TLB_FS)
		ADD_TEST(NAME SSE_${Ref}_${FASTQBaseFileName} COMMAND ${BASH} -c "rm -f ${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}*.gtl \
			&& $<TARGET_FILE:MicMap> -g ${GenomeFile} -p ${TablePrefix} -G 448 -a 1 -1 ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_1.fq -2 ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_2.fq -r ${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName} \
			&& ( $<TARGET_FILE:GTLdecompress> -g ${Genome} -r <(cat ${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_*.gtl) -p -n -m -a -u -h -c -d | sort -n | sed 's/^[0-9][0-9]* //' >${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_r1.fq ) 2>&1 | sort -n | sed 's/^[0-9][0-9]* //'>${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_r2.fq \
			&& cmp ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_1.fq ${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_r1.fq \
			&& cmp ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_2.fq ${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_r2.fq \
			&& echo \"SSE_${Ref}_${FASTQBaseFileName} matches original fastq files\" \
			&& $<TARGET_FILE:GTLdecompress> -g ${Genome} -r <(cat ${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_*.gtl) -p -n -m -a -u -h -c -o raw | sort -n >${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_raw.txt \
			&& cmp ${CMAKE_SOURCE_DIR}/data/${Ref}_${FASTQBaseFileName}_raw.txt ${GTL_TESTS_PATH}/SSE_${Ref}_${FASTQBaseFileName}_raw.txt \
			&& echo \"SSE_${Ref}_${FASTQBaseFileName} matches expected raw file\"")
		IF(USE_HUGEPAGE_TLB_FS)
			SET_TESTS_PROPERTIES(SSE_${Ref}_${FASTQBaseFileName} PROPERTIES DEPENDS LoadingTLBFS_${Ref})
		ENDIF(USE_HUGEPAGE_TLB_FS)
		IF(CPU_HAS_AVX512F)
			ADD_TEST(NAME AVX512F_${Ref}_${FASTQBaseFileName} COMMAND ${BASH} -c "rm -f ${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName}*.gtl \
				&& $<TARGET_FILE:MicMap> -g ${GenomeFile} -p ${TablePrefix} -G 448 -M 1 -x AVX512F -1 ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_1.fq -2 ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_2.fq -r ${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName} \
				&& ( $<TARGET_FILE:GTLdecompress> -g ${Genome} -r <(cat ${GTL_TESTS_PATH}/MIC_${Ref}_${FASTQBaseFileName}_*.gtl) -p -n -m -a -u -h -c -d | sort -n | sed 's/^[0-9][0-9]* //' >${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName}_r1.fq ) 2>&1 | sort -n | sed 's/^[0-9][0-9]* //'>${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName}_r2.fq \
				&& cmp ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_1.fq ${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName}_r1.fq \
				&& cmp ${CMAKE_SOURCE_DIR}/data/${FASTQBaseFileName}_2.fq ${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName}_r2.fq \
				&& echo \"AVX512F_${Ref}_${FASTQBaseFileName} matches original fastq files\" \
				&& $<TARGET_FILE:GTLdecompress> -g ${Genome} -r <(cat ${GTL_TESTS_PATH}/MIC_${Ref}_${FASTQBaseFileName}_*.gtl) -p -n -m -a -u -h -c -o raw | sort -n >${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName}_raw.txt \
				&& cmp ${CMAKE_SOURCE_DIR}/data/${Ref}_${FASTQBaseFileName}_raw.txt ${GTL_TESTS_PATH}/AVX512F_${Ref}_${FASTQBaseFileName}_raw.txt \
				&& echo \"AVX512F_${Ref}_${FASTQBaseFileName} matches expected raw file\"")
			SET_TESTS_PROPERTIES(AVX512F_${Ref}_${FASTQBaseFileName} PROPERTIES TIMEOUT 180)
			IF(USE_HUGEPAGE_TLB_FS)
				SET_TESTS_PROPERTIES(AVX512F_${Ref}_${FASTQBaseFileName} PROPERTIES DEPENDS LoadingTLBFS_${Ref})
			ENDIF(USE_HUGEPAGE_TLB_FS)
		ENDIF(CPU_HAS_AVX512F)
		# FIXME - need to add verification against a known-good reference
	ENDMACRO(SimpleTest)
		
	# ------------------------------------------ GENOMES --------------------------------------------------------------
	SET(GENOME_B37 0)
	FIND_PATH(GENOME_B37_PATH NAMES b37/b37.bin
	          PATHS /home/tschuepb/Data/Genomes /data6 "${GTL_TEST_GENOME_PATH}"
	          NO_DEFAULT_PATH )
	IF (NOT GENOME_B37_PATH-NOTFOUND)
		SET(SUFFIXES A.2nt.bin;A.chr.bin;A.strand.bin;A.valid.bin)
		SET(ALL_FOUND 1)
		FOREACH(extension IN ITEMS ${SUFFIXES})
			MESSAGE(STATUS "Testing for ${GENOME_B37_PATH}/b37/tbl_18nt/start0.${extension}")
			IF(NOT EXISTS "${GENOME_B37_PATH}/b37/tbl_18nt/start0.${extension}")
				SET(ALL_FOUND 0)
			ENDIF(NOT EXISTS "${GENOME_B37_PATH}/b37/tbl_18nt/start0.${extension}")
		ENDFOREACH(extension IN ITEMS ${SUFFIXES})
		IF(ALL_FOUND)
			MESSAGE(STATUS "Genome b37 found in ${GENOME_B37_PATH}")
			SET(GENOME_B37 1)
			SET(GENOME_B37_FILE "${GENOME_B37_PATH}/b37/b37.bin")
			SET(GENOME_B37_TABLE_PREFIX "${GENOME_B37_PATH}/b37/tbl_18nt/start0.")
		ENDIF(ALL_FOUND)		
	ENDIF(NOT GENOME_B37_PATH-NOTFOUND)

	SET(GENOME_B38 0)
        FIND_PATH(GENOME_B38_PATH NAMES b38/b38.bin
                  PATHS /home/tschuepb/Data/Genomes /data6 "${GTL_TEST_GENOME_PATH}"
                  NO_DEFAULT_PATH )
        IF (NOT GENOME_B38_PATH-NOTFOUND)
                SET(SUFFIXES A.2nt.bin;A.chr.bin;A.strand.bin;A.valid.bin)
                SET(ALL_FOUND 1)
                FOREACH(extension IN ITEMS ${SUFFIXES})
                        MESSAGE(STATUS "Testing for ${GENOME_B38_PATH}/b38/tbl_18nt/start0.${extension}")
                        IF(NOT EXISTS "${GENOME_B38_PATH}/b38/tbl_18nt/start0.${extension}")
                                SET(ALL_FOUND 0)
                        ENDIF(NOT EXISTS "${GENOME_B38_PATH}/b38/tbl_18nt/start0.${extension}")
                ENDFOREACH(extension IN ITEMS ${SUFFIXES})
                IF(ALL_FOUND)
                        MESSAGE(STATUS "Genome b38 found in ${GENOME_B38_PATH}")
                        SET(GENOME_B38 1)
                        SET(GENOME_B38_FILE "${GENOME_B38_PATH}/b38/b38.bin")
                        SET(GENOME_B38_TABLE_PREFIX "${GENOME_B38_PATH}/b38/tbl_18nt/start0.")
                ENDIF(ALL_FOUND)
        ENDIF(NOT GENOME_B38_PATH-NOTFOUND)

	# ------------------------------------------ GTL FILES ------------------------------------------------------------
	FIND_PATH(GTL_TEST_FILES_PATH NAMES chr9c_p.gtl
	          PATHS /home/tschuepb/Data/Caller/HG001/ /data6/FDA/Schupi/HG001
	          NO_DEFAULT_PATH )


	# ------------------------------------------ TESTS ----------------------------------------------------------------
	IF (GENOME_B38)	
		IF(USE_HUGEPAGE_TLB_FS)
			TestGenomeTLBLoading(${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		ENDIF(USE_HUGEPAGE_TLB_FS)
		SimpleTest(al_sch ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(badCode ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(badPos ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(bad_S ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(bug_I_D ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(bug_sa1 ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(bug_sa2 ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(bug_sa3 ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(chim ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(cigarbug ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(lowq ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(mis_a ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(mis_a2 ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(mis_m ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(mis_n ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(mis_p ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(mis_p2 ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(mis_rev_gap ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(miss_a ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(miss_b ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(refIsN ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(strange ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(test ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
		SimpleTest(wrong ${GENOME_B38_FILE} ${GENOME_B38_TABLE_PREFIX} b38)
	ENDIF(GENOME_B38)
	
	IF (GENOME_B37)
		IF(USE_HUGEPAGE_TLB_FS)
			TestGenomeTLBLoading(${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		ENDIF(USE_HUGEPAGE_TLB_FS)
		SimpleTest(al_sch ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(badCode ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(badPos ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(bad_S ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(bug_I_D ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(bug_sa1 ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(bug_sa2 ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(bug_sa3 ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(chim ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(cigarbug ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(lowq ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(mis_a ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(mis_a2 ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(mis_m ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(mis_n ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(mis_p ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(mis_p2 ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(mis_rev_gap ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(miss_a ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(miss_b ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(refIsN ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(strange ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(test ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(wrong ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)
		SimpleTest(trim_err ${GENOME_B37_FILE} ${GENOME_B37_TABLE_PREFIX} b37)	
	
		IF(EXISTS "${GTL_TEST_FILES_PATH}/chr9c_p.gtl")
			DecompressExistingGTL(chr9c_p.gtl chr9c_p 9 "-p")
			TestMapping(MicMap    "SSE" chr9c_p "-G 800 -A 8")
		ENDIF(EXISTS "${GTL_TEST_FILES_PATH}/chr9c_p.gtl")

		IF(EXISTS "${GTL_TEST_FILES_PATH}/chr8c_p.gtl")
			DecompressExistingGTL(chr8c_p.gtl chr8c_p 8 "-p")
			TestMapping(MicMap    "SSE" chr8c_p "-G 800 -A 8")
		ENDIF(EXISTS "${GTL_TEST_FILES_PATH}/chr8c_p.gtl")

		IF(EXISTS "${GTL_TEST_FILES_PATH}/chr21a_g.gtl")
			DecompressExistingGTL(chr21a_g.gtl chr21a_g 21 "-a")
			TestMapping(MicMap    "SSE" chr21a_g "-G 800 -A 8")

			IF(EXISTS "${GTL_TEST_FILES_PATH}/chr21a_g.gtl.idr")
				IF(EXISTS "${GTL_TESTS_PATH}/chr21a_g.gtl.idr")
					FILE(REMOVE "${GTL_TESTS_PATH}/chr21a_g.gtl.idr")
				ENDIF(EXISTS "${GTL_TESTS_PATH}/chr21a_g.gtl.idr")
				ADD_TEST(NAME IndexChr21
								COMMAND ${BASH} -c "$<TARGET_FILE:GTLindex> -g ${GTL_TEST_GENOME_PATH}/b37/b37.bin ${GTL_TEST_FILES_PATH}/chr21a_g.gtl -f -o ${GTL_TESTS_PATH}" )
				ADD_TEST(NAME IndexChr21Cmp
								COMMAND ${CMAKE_COMMAND} -E compare_files "${GTL_TEST_FILES_PATH}/chr21a_g.gtl.idr" "${GTL_TESTS_PATH}/chr21a_g.gtl.idr" )
			ENDIF(EXISTS "${GTL_TEST_FILES_PATH}/chr21a_g.gtl.idr")

		ENDIF(EXISTS "${GTL_TEST_FILES_PATH}/chr21a_g.gtl")	
	ENDIF(GENOME_B37)
ENDIF(BUILD_TESTING)

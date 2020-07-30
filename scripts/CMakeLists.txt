INSTALL (FILES
	   analyzeVarRegions
           annotateVCF
	   ExomeDiffCall
	   genRGBsnp
	   grabAllContacts
	   grabAllNonContacts
	   mergeVCFandSelect
	   mergeVCFdetails
	   MicMapChIP2bam
	   MicMapExome2bam
	   MicMapExomeUMI2bam
	   MicMapHiC2table
	   MicMapRNA2bam
	   MicMapRNAse2bam
	   sortCompressCall
	   sortVCF.sh
	   trimAndSplit
	   tromer2sam
         DESTINATION bin
         PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

INSTALL (FILES
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
	   sortCompressCall
	   sortVCF.sh
	   tromer2sam
         DESTINATION bin
         PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

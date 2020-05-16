INSTALL (FILES
           annotateVCF
	   genRGBsnp
	   grabAllContacts
	   grabAllNonContacts
	   MicMapChIP2bam
	   MicMapExomeUMI2bam
	   MicMapHiC2table
	   MicMapRNA2bam
	   sortCompressCall
	   sortVCF.sh
	   tromer2sam
         DESTINATION bin
         PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

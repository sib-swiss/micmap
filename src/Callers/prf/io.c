/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland      2019 Nicolas Guex and Christian Iseli
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *
 *      Code:       Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *      Contacts:   Nicolas.Guex@unil.ch and Christian.Iseli@unil.ch
 *      Repository: https://github.com/sib-swiss/micmap
 *
 * ------------------------------------------------------------------------------------------------------------------------
 */
/*******************************************************
                        PFTOOLS
 *******************************************************
  Jan 18, 2011 io.c
 *******************************************************
 (C) 2011 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <errno.h>
#define _GNU_SOURCE
#include <string.h>
#include "pfProfile.h"
#include "pfIOInline.h"
#ifdef PRF_CORE_PCRE
#include "pfRegexp.h"
#endif
const char NormalizationModeName[3][16] = { 
  "linear\0\0\0\0\0\0\0\0\0\0",
  "gle_zscore\0\0\0\0\0\0",
  "gle_zscave\0\0\0\0\0\0"
};

#ifdef _TEST
#define _VERBOSE_
unsigned int UseColor = 0;
unsigned int out_profile = 0;
#endif


/* TODO:
 *       - Profile length is zero based. We therefore added an extra line to account for it.
 *         Nevertheless, no memory space checking is so far implemented.
 *         So potentially we may completely screw up upon incorrect profiles and segfault!
 */

#define SUB_COMMAND_MAX_SIZE 16
#ifndef _BEST_IS_NEGATIVE_
#define MAXMLOW(a) (StoredIntegerFormat) ( ((a)>MLOW) ? (a) : MLOW )
#define CHECK_AND_SET(a,b,c) {\
  register const int Lscore = ((int) (b) ) + ((int) (c));\
  if (Lscore < STORED_INT_MIN || Lscore > STORED_INT_MAX) {\
    fprintf(stderr, "Profile contains out of bound values for StoredIntegerFormat : %i arising from the addition of %i with %i\n", Lscore, (int)(b), (int)(c));\
    return 1;\
  } else {\
    a = (StoredIntegerFormat) ( (Lscore > MLOW) ? Lscore : MLOW );\
  }\
}
#else
#define MAXMLOW(a) (StoredIntegerFormat) ( (a<MLOW) ? a : MLOW )
#define CHECK_AND_SET(a,b,c) {\
  register const int Lscore = ((int) b ) + ((int) c);\
  if (Lscore < STORED_INT_MIN || Lscore > STORED_INT_MAX) {\
    fprintf(stderr, "Profile contains out of bound values for StoredIntegerFormat : %i arising from the addition of %i with %i\n", Lscore, (int)b, (int)c);\
    return 1;\
  } else {\
    a = (StoredIntegerFormat) ( (Lscore < MLOW) ? Lscore : MLOW );\
  }\
}
#endif

/*
 * This structure is fed by the Line Analyzer which parses the input, transforms delimiting symbols
 * into '\0' end character and sets the starting pointer within the ProfileLine structure.
 */
struct ProfileLine { 
  char * command;
  char * subcommand;
  char * keywords[32];
  char * values[32];
  size_t counter;
  _Bool isMatrix;
};

// static char previousSubCommand[SUB_COMMAND_MAX_SIZE];

/*
 * Reads a line of profile and corrects the ending with a '\0' character.
 */
static inline size_t GetLine(FILE * const stream, char * const restrict Destination, const size_t MaxLength)
{
  char * const read = fgets(Destination, MaxLength, stream);
  if (read != NULL) {
    size_t length = strlen(Destination);
    if (Destination[length-1] == '\n') Destination[--length] = '\0';
    return length;
  }
  else 
    return 0;
}
//---------------------------------------------------------------

/*
 * Move to the first non space or tab character while within memory up to MaxMempoint.
 */
static inline char * FirstChar(const char * restrict Line, const uintptr_t MaxMemPoint)
{
  char * pos = (char*) Line; 
  while ( (*pos == ' ') || (*pos == '\t') ) {
    if ( (uintptr_t) ++pos == MaxMemPoint) break;
  }
  return pos;
}
//---------------------------------------------------------------

/*
 * Copy the string Line into Destination until a given symbol is found. Max length applies.
 * WARNING: Destination should be cleared before since no '\0' will be added.
 */
static inline char * CopyUptoSymbol(char * const restrict Destination, const char * restrict Line, const char Symbol, const size_t MaxLength)
{
   for (size_t i=0; i<MaxLength; i++) {
    if ( Line[i] == Symbol){
      return (char*) &Line[i+1];
    }
    else 
      Destination[i] = Line[i];
  }
  return NULL;
}
//---------------------------------------------------------------

/*
 * Copy the string Line into Destination until a given symbol is found. Max length applies.
 */
static inline char * CopyUptoSymbolTerminated(char * const restrict Destination, const char * restrict Line, const char Symbol, const size_t MaxLength)
{
   for (size_t i=0; i<MaxLength; i++) {
    if ( Line[i] == Symbol){
      Destination[i] = '\0';
      return (char*) &Line[i+1];
    }
    else 
      Destination[i] = Line[i];
  }
  Destination[MaxLength] = '\0';
  return NULL;
}
//---------------------------------------------------------------

static inline char * CopySubCommand(char * const restrict Destination, char * restrict Line, const char Symbol, const uintptr_t MaxLength)
{
#ifdef _DEBUG_VERBOSE_
  fputs(" COPY_SUB_COMMAND : '",stdout);
#endif
  for (size_t i=0; i<SUB_COMMAND_MAX_SIZE-1; ++i) {
    if ((uintptr_t) &Line[i] > MaxLength) break;

    if ( Line[i] == Symbol){
      Destination[i] = '\0';
      Line[i] = '\0';
#ifdef _DEBUG_VERBOSE_
      puts("'");
#endif
      return (char*) &Line[i+1];
    }
    Destination[i] = Line[i];
#ifdef _DEBUG_VERBOSE_
    fputc((int) Line[i], stdout);
#endif
  }
  return NULL;
}
//---------------------------------------------------------------

/*
 * Trim the string from the right
 */
static inline void CleanEndString(char * const restrict String, const size_t Length)
{
  char * pos = &String[Length-1];
  while ( pos >= String && *pos == ' ') *pos-- = '\0';
}
//---------------------------------------------------------------

/*
 * Move up to symbol while staying within memory limit.
 */
static inline char * GoUptoSymbol(char * const restrict String, const char Symbol, const uintptr_t MaxMemPoint)
{
  char * pos = String;
  while ( (uintptr_t) pos <= MaxMemPoint && *pos != Symbol) ++pos;
  if (*pos == Symbol) {
    *pos = '\0';
    ++pos;
  }
  return pos;
}
//---------------------------------------------------------------

/*
 * Count and replace symbol with '\0' character within given string. Overall is bounded by a maximum memory position
 * Pay attention that not all count will give correct amount of tags. Indeed commas are n-1 tags.
 */
static inline size_t CountAndReplaceSymbol(char * const restrict String, const char Symbol, const uintptr_t MaxMemPoint)
{
  char * pos = String;
#ifdef _DEBUG_VERBOSE_
  fputs("CountAndReplaceSymbol : ",stdout); fflush(stdout); int c=0;
  while ( (uintptr_t) pos <= MaxMemPoint ) {
    printf("%i @0x%16.16x of 0x%16.16x\t", c++, (uintptr_t) pos, MaxMemPoint);
    fputc((int) *pos,stdout);
    fputc((int) '\n',stdout);
    pos++;
    fflush(stdout);
  }
  fputs("\n", stdout);
  fputs("CountAndReplaceSymbol : ",stdout);fflush(stdout);
#endif
  pos = String;
  size_t count = 0;
  while ( (uintptr_t) pos <= MaxMemPoint ) {
    if (*pos == Symbol) {
      ++count;
#ifdef _DEBUG_VERBOSE_
      if (UseColor) fputs("\e[31;40m",stdout);
      fputc((int) *pos,stdout);
      if (UseColor) fputs("\e[0m",stdout);
#endif
      *pos = '\0';
    }
#ifdef _DEBUG_VERBOSE_
    else {
      fputc((int) *pos,stdout);
    } 
#endif
    ++pos;
  }
#ifdef _DEBUG_VERBOSE_
  printf("\t count = %lu\n", count);
#endif
  return count;
}
//---------------------------------------------------------------

/*
 * Count given symbol within given string. Overall is bounded by a maximum memory position
 * Pay attention that not all count will give correct amount of tags. Indeed commas are n-1 tags.
 */
static inline size_t CountSymbol(const char * const restrict String, const char Symbol, const uintptr_t MaxMemPoint)
{
  const char * pos = String;
#ifdef _DEBUG_VERBOSE_
  fputs("CountSymbol : ",stdout); fflush(stdout); int c=0;
  while ( (uintptr_t) pos <= MaxMemPoint ) {
    printf("%i @0x%16.16x of 0x%16.16x\t", c++, (uintptr_t) pos, MaxMemPoint);
    fputc((int) *pos,stdout);
    fputc((int) '\n',stdout);
    pos++;
    fflush(stdout);
  }
  fputs("\n", stdout);
  fputs("CountSymbol : ",stdout);fflush(stdout);
#endif
  pos = String;
  size_t count = 0;
  while ( (uintptr_t) pos <= MaxMemPoint ) {
    if (*pos == Symbol) {
      ++count;
#ifdef _DEBUG_VERBOSE_
      if (UseColor) fputs("\e[31;40m",stdout);
      fputc((int) *pos,stdout);
      if (UseColor) fputs("\e[0m",stdout);
#endif
    }
#ifdef _DEBUG_VERBOSE_
    else {
      fputc((int) *pos,stdout);
    } 
#endif
    ++pos;
  }
#ifdef _DEBUG_VERBOSE_
  printf("\t count = %lu\n", count);
#endif
  return count;
}
//---------------------------------------------------------------

/*
 * Profile Line parser.
 */
static int AnalyzeLine(char * restrict currentLine, struct ProfileLine * const restrict prfLine,
                       _Bool * restrict MultipleLine, char * const restrict previousSubCommand)
{
  const size_t LineSize       = strlen(currentLine);
  const uintptr_t MaxMemPoint = (uintptr_t) &currentLine[LineSize-1];
#ifdef _VERBOSE_
  if (out_profile) {
    if (UseColor) 
      printf("\e[31;40mLine : %s\e[0m\n", currentLine);
    else
      printf("Line : %s\n", currentLine);
  }
#endif
  /* Set Multiple line to false */
  *MultipleLine = false;
  
  /* Clear result first */
  memset(prfLine, 0, sizeof(struct ProfileLine)-SUB_COMMAND_MAX_SIZE);

  /* Get Line command */
  char * Position  = FirstChar(currentLine, MaxMemPoint);
  prfLine->command = Position;
  const _Bool isMatrix = (Position[0] == 'M' && Position[1] == 'A') ? true : false;
  prfLine->isMatrix = isMatrix;
  
  /* Get to end of command and replace space by '\0' */
  while ( *Position != ' ' && (uintptr_t)Position <= MaxMemPoint) ++Position;
  while ( *Position == ' ' && (uintptr_t)Position <= MaxMemPoint) *Position++ = '\0';
  if ((uintptr_t) Position == MaxMemPoint) return 1;
#ifdef _VERBOSE_
  if (out_profile) printf("\t'%s'", prfLine->command);
#endif
  if (!isMatrix) {
    prfLine->subcommand = Position;
    for (size_t i=0; i<SUB_COMMAND_MAX_SIZE;++i) previousSubCommand[i] = '\0';
    CleanEndString(currentLine, LineSize);
#ifdef _VERBOSE_
    if (out_profile) printf("\t\t'%s'\n", prfLine->subcommand);
#endif
  } else {
    if (*Position == '/') {
      prfLine->subcommand = ++Position;
      Position = CopySubCommand(previousSubCommand, Position, ':', MaxMemPoint);
      //Position = GoUptoSymbol(Position, ':', MaxMemPoint);
    } else {
#ifdef _VERBOSE_
  if (out_profile) {
    if (UseColor) 
      printf(" \e[31;40mMULTIPLE MA LINES DETECTED\e[0m");
    else
      printf(" MULTIPLE MA LINES DETECTED");
  }
#endif 
      prfLine->subcommand = previousSubCommand;
      *MultipleLine = true;
    }
#ifdef _VERBOSE_
    if (out_profile) printf("\t\t'%s'\n", prfLine->subcommand);
#endif
    //if ((uintptr_t)Position >= MaxMemPoint) return 1; causes bug when line is /MA M:
    // in any case not needed
    const size_t counter = CountAndReplaceSymbol(Position, ';', MaxMemPoint);
#ifdef _VERBOSE_
    if (out_profile) printf("\t\t%lu tags found\n", counter);
#endif
    if (counter > 32) return 2;
    for (size_t i=0; i<counter; ++i) {
      /* TODO: Avoid checking from the beginning */
        char * tpos = GoUptoSymbol(Position, '=', MaxMemPoint);
        if ((uintptr_t) tpos >= MaxMemPoint) {
          fprintf(stderr, "\nError no '=' symbol detected in input line %s\n", currentLine);
          return 3;
        } else {
          prfLine->values[i] = tpos;
	  /* get back to = symbol */
// 	  while ( *tpos != '\0' && tpos > Position) --tpos;
	  const char * EqualPosition = tpos - 1;
	  /* get back to ; symbol */
          tpos -= 2;
          while ( *tpos != '\0' && tpos >= Position) --tpos;
	  ++tpos;
	  /* trim early space */
	  while ( *tpos == ' ' && tpos < EqualPosition) ++tpos;
          prfLine->keywords[i] = tpos;
        }
#ifdef _VERBOSE_
#ifdef _DEBUG_VERBOSE_
	printf("\t\t\t tag:'%s' @ 0x%8.8x \t\t value:'%s' @ 0x%8.8x\n", prfLine->keywords[i],prfLine->keywords[i], prfLine->values[i], prfLine->values[i]);
#else
        if (out_profile) printf("\t\t\t tag:'%s' \t\t value:'%s'\n", prfLine->keywords[i], prfLine->values[i] );
#endif
#endif
    }
    prfLine->counter = counter;
  }
  return 0;
}
//---------------------------------------------------------------

/*
 * The following functions read the score(s) from a given string and set the corresponding
 * destination(s) with the StoredIntegerFormat (16 or 32 bits) values read.
 * Upon error, errno should hold the operating POSIX error.
 */
static inline size_t ReadScore(const char * const restrict String, StoredIntegerFormat * const restrict Score)
{
  errno = 0;
  if (String[0] == '*') {
    *Score = NLOW;
    return 0;
  } else {
    const long int Lscore =  strtol(String, NULL, 10);
    if (Lscore < STORED_INT_MIN || Lscore > STORED_INT_MAX) {
	fprintf(stderr, "Profile contains out of bound values for StoredIntegerFormat : %li\n", Lscore);
	exit(1);
    }
    *Score = (StoredIntegerFormat) Lscore;
    return errno;
  }
}
//---------------------------------------------------------------

static inline size_t ReadScores(const char * const restrict String, StoredIntegerFormat * const restrict Scores, const size_t count)
{
  const char * restrict pos = String;

  errno = 0;
  if (*pos == '\'') ++pos;
  for (size_t i=0; i<count; ++i) {
    if (*pos == '*') {
      Scores[i] = NLOW;
      pos += 2;
    } else {
      char * newpos;
      const long int Lscore =  strtol(pos, &newpos, 10);
      if (Lscore < STORED_INT_MIN || Lscore > STORED_INT_MAX) {
	  fprintf(stderr, "Profile contains out of bound values for StoredIntegerFormat : %li\n", Lscore);
	  exit(1);
      }
      Scores[i] = (StoredIntegerFormat) Lscore;
      if (errno == 0) {
        while(*newpos != '\0') ++newpos;
        pos = (const char*) newpos + 1;
      } else {
        return 1;
      }
    }
//     if ( count > 1) fprintf(stderr, "%4i ", Scores[i]);
  }
//   if ( count > 1) fputs("\n", stderr);
  return 0;
}
//---------------------------------------------------------------

static int internalReadProfile(FILE* prfStream, struct Profile * const prf, const _Bool SetExtraTable,
															 const char * const FileName, size_t * const restrict LineOffset,
															 const _Bool CompleteCycleOnly)
{
  StoredIntegerFormat IIPD_ALPHABET       [      ALPHABET_MEMORY_SIZE] __attribute__((aligned(16)));
  StoredIntegerFormat IIPD_BOUNDARIES     [ INSERTION_BOUNDARIES_SIZE] __attribute__((aligned(16)));
  StoredIntegerFormat IMPD_ALPHABET       [      ALPHABET_MEMORY_SIZE] __attribute__((aligned(16)));
  TransitionScores IIPD_TRANSITIONS __attribute__((aligned(16)));
  
  char currentLine[PROFILE_MAX_LINE_SIZE] __attribute__((aligned(16)));
  char previousSubCommand[SUB_COMMAND_MAX_SIZE];

  struct ProfileLine AnalyzedLine;

  union Scores DefaultScores = { 
    .Insertion = { IIPD_ALPHABET,
		  IIPD_BOUNDARIES,
		  &IIPD_TRANSITIONS,
		  0,
		  0,
		  0,
		  { IMPD_ALPHABET } }
  };
  union Scores WorkingScores;
  char * ProfileSequence = 0;
  size_t Line = *LineOffset, Alphabet_Length=0, Length=0;
  
  int InsertionCounter=-1, MatchCounter=0;
  int res;
  _Bool LZCO = false;
  _Bool MultipleLine;
  char DefaultMatchSymbol;
  char DefaultInsertionSymbol;

   
  /*
   * Initialize position-independent profile parameters
   *
   *   - general specification
   */
  
  prf->isCircular = false;
	prf->CompleteCycleOnly = CompleteCycleOnly;
  prf->Length     = 0;

  /*   - disjoint mode */
  prf->DisjointData.MDIS = 1;
  strcpy(prf->DisjointData.CDIS[0], "UNIQUE\0");
  strcpy(prf->DisjointData.CDIS[1], "PROTECT\0");
  prf->DisjointData.JDIP[0] = 0;
  prf->DisjointData.JDIP[1] = 2;
  

  /*   - normalization modes */
  prf->NormalizationData.JNOR = 0;
  strcpy(prf->NormalizationData.CNOR[0], "LINEAR\0");
  strcpy(prf->NormalizationData.CNOR[1], "GLE_ZSCORE\0");
  strcpy(prf->NormalizationData.CNOR[2], "GLE_ZSCAVE\0");
  prf->NormalizationData.JNOP[0] = 2;
  prf->NormalizationData.JNOP[1] = 5;
  prf->NormalizationData.JNOP[2] = 5;

  for (int i=0; i<MAXN; ++i) {
      prf->NormalizationData.Values[i].NNOR = i;
      prf->NormalizationData.Values[i].NNPR = i;
  }

  /*   - cut-off */
  prf->CutOffData.JCUT = 0;
  
  /* Pattern */
  prf->Pattern = NULL;
  
  InitializeDefault(&DefaultScores, &DefaultMatchSymbol, &DefaultInsertionSymbol);
  
  memset(previousSubCommand, 0, sizeof(char)*SUB_COMMAND_MAX_SIZE);
  
  /*
   * Read profile
   */
#ifdef _VERBOSE_
  _Bool FirstCall = true;
#endif
  
  while (!feof(prfStream)) {   
    const size_t length = GetLine(prfStream, currentLine, PROFILE_MAX_LINE_SIZE);
    if (length == 0) continue;
    ++Line;
#ifdef _VERBOSE_
    if (out_profile) {
      if (FirstCall) {
	printf("============================== NEW PROFILE ================================================== Input line %lu\n", Line);
	FirstCall = false;
      }
      else
	printf("--------------------------------------------------------------------------------------------- Input line %lu\n", Line);
      if (UseColor)
	printf("\e[33;40m\t Match: %i\t\tInsertion: %i\e[0m\n", MatchCounter, InsertionCounter);
      else
	printf("\t Match: %i\t\tInsertion: %i\n", MatchCounter, InsertionCounter);
    }
#endif

    if ((res=AnalyzeLine(currentLine, &AnalyzedLine, &MultipleLine, previousSubCommand)) != 0) {
      fprintf(stderr,"Error %i in analysis at line %lu of %s\n", res, Line, FileName);
      return 1;
    }
    
    /* Header or matrix */
    if (!AnalyzedLine.isMatrix)
    {
       if ( (AnalyzedLine.command[0] == 'I' && AnalyzedLine.command[1] == 'D' ) ) {
         char * tpos = CopyUptoSymbol(prf->Identification, AnalyzedLine.subcommand, ';', 64);
         if (strstr(tpos, "MATRIX") != NULL) {
				   prf->Type = PF_MATRIX;
				 }
				 else if (strstr(tpos, "PATTERN") != NULL) {
				   prf->Type = PF_PATTERN;
				 }
				 else {
				   goto MissingType;  
				 }
       } 
       else if ( (AnalyzedLine.command[0] == 'A' && AnalyzedLine.command[1] == 'C' ) ) {
          CopyUptoSymbol(prf->AC_Number, AnalyzedLine.subcommand, ';', 64);
       }
       else if ( (AnalyzedLine.command[0] == 'D' && AnalyzedLine.command[1] == 'T' ) ) {
          strncpy(prf->Date, AnalyzedLine.subcommand, 127);
       }
       else if ( (AnalyzedLine.command[0] == 'D' && AnalyzedLine.command[1] == 'E' ) ) {
          strncpy(prf->Description, AnalyzedLine.subcommand, 255);
       }
       else if ( (AnalyzedLine.command[0] == 'P' && AnalyzedLine.command[1] == 'A' ) ) {
#ifdef PRF_CORE_PCRE
			  const size_t len = strlen(AnalyzedLine.subcommand);
			  char * tmpPattern;
			  if (prf->Pattern) {

			    const size_t len2 = strlen(prf->Pattern);
			    tmpPattern = (char*) malloc((1+len+len2)*sizeof(char));
			    if (prf->Pattern == NULL) goto AllocationError;
			    for (size_t i=0; i<len2; ++i) tmpPattern[i] = prf->Pattern[i];
			    free(prf->Pattern);
			    prf->Pattern = tmpPattern;
			    tmpPattern += len2;
			  }
			  else {
			    prf->Pattern = (char*) malloc((1+len)*sizeof(char));
			    if (prf->Pattern == NULL) goto AllocationError;
			    tmpPattern = prf->Pattern;
			  }
			  CopyUptoSymbolTerminated(tmpPattern, AnalyzedLine.subcommand, '.', len);
#ifdef _VERBOSE_
			  char *const TranslatedRegex = PatternToRegex(prf->Pattern);
			  printf("\t\t\'%s\' -> \'%s\'\n", prf->Pattern, TranslatedRegex);
			  free(TranslatedRegex);
#endif
#else
				fputs("Pattern profile have been deactivated at compile time!\n", stderr);
				goto FIN;
#endif
	     }
       else if ( (AnalyzedLine.command[0] == '/' && AnalyzedLine.command[1] == '/' ) ) {
#ifdef _VERBOSE_
		 if (out_profile) {
		  if (UseColor)
		    puts("\e[33;41m\tEND OF PROFILE REACHED\e[0m");
		  else
		    puts("\tEND OF PROFILE REACHED");
		 }
#endif
	  goto END_OF_PROFILE;
       }
#ifdef _VERBOSE_
       else {
	  if (out_profile) {
	    if (UseColor)
	      puts("\e[33;41m\t!!! UNTREATED COMMAND !!!\e[0m");
	    else
	      puts("\t!!! UNTREATED COMMAND !!!");
	  }
       }
#endif
    }
    else
    {
      if (prf->Type == PF_PATTERN) goto WrongType;
      /* ----------------------- GENERAL SPECIFICATIONS -------------------------------*/
      if (strcmp(AnalyzedLine.subcommand, "GENERAL_SPEC")==0) {
        for (size_t keys=0; keys<AnalyzedLine.counter;++keys) {
          if        (strcmp(AnalyzedLine.keywords[keys], "ALPHABET")==0) {
            Alphabet_Length = strlen(AnalyzedLine.values[keys]) - 2;
            if ( Alphabet_Length > ALPHABET_SIZE ) goto AlphabetSizeTooLarge;

            /* Map all character to unknown set as 0 */
            memset (prf->Alphabet_Mapping, 0, (ALPHABET_SIZE+1)*sizeof(char));
            memset (prf->CABC, 'X', (ALPHABET_SIZE+1)*sizeof(char));
				    prf->CABC[ALPHABET_SIZE+1] = '\0';

            // Update the mapping
            for (size_t i=1; i<=Alphabet_Length; ++i) {
              prf->CABC[i] = AnalyzedLine.values[keys][i];
              const size_t index = (size_t) ((unsigned char) AnalyzedLine.values[keys][i] - (unsigned char) 'A');
              prf->Alphabet_Mapping[index] = (unsigned char) i;
            }

            Alphabet_Length = ALPHABET_SIZE;
            prf->Alphabet_Length = Alphabet_Length;
            
          } else if (strcmp(AnalyzedLine.keywords[keys], "LENGTH")==0) {
            Length = (size_t) atoi(AnalyzedLine.values[keys]);
            prf->Length = Length;
	               
            // Allocates memory
            prf->Sequence = (char*) calloc((Length+1), sizeof(char));
				    if (prf->Sequence == NULL) goto AllocationError;
				    ProfileSequence = prf->Sequence;
				    
				    //if ( AllocateScores(&DefaultScores, Alphabet_Length, Length) != 0 ) goto AllocationError; THIS IS IN THE STACK
            if ( AllocateScores(&WorkingScores, Alphabet_Length, Length) != 0 ) goto AllocationError;

            /* Copy score pointers to profile */
            memcpy(&(prf->Scores), &WorkingScores, sizeof(union Scores));

          } else if (strcmp(AnalyzedLine.keywords[keys], "TOPOLOGY")==0) {
						if (strcmp(AnalyzedLine.values[keys], "CIRCULAR") == 0) {
							prf->isCircular = true;
						}

          } else if (strcmp(AnalyzedLine.keywords[keys], "LOG_BASE")==0) {

          } else if (strcmp(AnalyzedLine.keywords[keys], "P0")==0) {

          } else if (strcmp(AnalyzedLine.keywords[keys], "P")==0) {

          }
        }
      }
      /* -----------------------      DISJOINT          -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "DISJOINT")==0) {
        SDisjoint * const djt = &(prf->DisjointData);
        /* Set some default in case data is missing */
        djt->NDIP[0] = 1;
        djt->NDIP[1] = Length;
        for (size_t keys=0; keys<AnalyzedLine.counter;++keys) {
          if        (strcmp(AnalyzedLine.keywords[keys], "DEFINITION")==0) {
            const char * const tval = AnalyzedLine.values[keys];
            for (size_t i=0; i<KDIS; ++i) {
              if ( strcmp(tval, djt->CDIS[i])==0) djt->MDIS = (int) i;
            }
          } else if (AnalyzedLine.keywords[keys][0] == 'N' && AnalyzedLine.keywords[keys][1] == '1' ) {
            djt->NDIP[0] = atoi(AnalyzedLine.values[keys]);
          } else if (AnalyzedLine.keywords[keys][0] == 'N' && AnalyzedLine.keywords[keys][1] == '2' ) {
            djt->NDIP[1] = atoi(AnalyzedLine.values[keys]);
          }
        }
      }
      /* -----------------------   NORMALIZATION        -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "NORMALIZATION")==0) {
        SNormalization * const nrm = &(prf->NormalizationData);
        register const size_t JNOR = nrm->JNOR;
        if (JNOR >= MAXN) goto TooManyNormalization;
				register SNormalizationItem * const nrmItem = &(nrm->Values[JNOR]);
        nrmItem->CNTX[0] = ' ';
        
        for (size_t keys=0; keys<AnalyzedLine.counter;++keys) {
          if        (strcmp(AnalyzedLine.keywords[keys], "FUNCTION")==0) {
	    char ctmp2[] = "GLE_ZSCORE";
            if (strcmp(AnalyzedLine.values[keys], "GRIBSKOV") == 0)
              AnalyzedLine.values[keys] = ctmp2;

            int index = -1;
            for (int i=0; i<KNOR; ++i) {
              if (strcmp(AnalyzedLine.values[keys], nrm->CNOR[i])==0) index = i;
            }
            if ( index < 0 ) goto NormalizationError;
            nrmItem->MNOR = index;

          } else if (strcmp(AnalyzedLine.keywords[keys], "MODE")==0) {
            nrmItem->NNOR = atoi(AnalyzedLine.values[keys]);
          } else if (strcmp(AnalyzedLine.keywords[keys], "PRIORITY")==0) {
            nrmItem->NNPR = atoi(AnalyzedLine.values[keys]);
          } else if (strcmp(AnalyzedLine.keywords[keys], "TEXT")==0) {
	    // Move to the first quote
	    const char * pos = AnalyzedLine.values[keys];
            const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
	    pos = GoUptoSymbol(AnalyzedLine.values[keys], '\'', MaxMemory);
	    // Copy up to the second quote
	    CopyUptoSymbol(nrmItem->CNTX, pos, '\'', 32);
          } else {
            if ( AnalyzedLine.keywords[keys][0] == 'R') {
              // get the number
              const size_t index = (size_t) ( (unsigned char) AnalyzedLine.keywords[keys][1] - (unsigned char) '1');
	      if (index >= KNPM) {
		fprintf(stderr, "Normalization R index out of bound (%lu).\n", index);
		return 1;
	      }
              nrmItem->RNOP[index] = (float) atof(AnalyzedLine.values[keys]);
	      {
		char * endptr;
		nrmItem->RNOP[index] = (float) strtof(AnalyzedLine.values[keys], &endptr);
		if ( endptr == NULL ) {
		    fprintf(stderr, "Line %lu : unable to convert %s to float\n", Line, AnalyzedLine.values[keys]);
		    return 1;
		} 
		if ( endptr == AnalyzedLine.values[keys] ) {
		    fprintf(stderr, "Line %lu : conversion error of %s to float\n"
				    "         : error code %i = %s",
			    Line, AnalyzedLine.values[keys], errno, strerror(errno));
		    
		    return 1;
		} 
	      }
	    }
          }
        }
        nrm->JNOR = JNOR + 1;
      }
      /* -----------------------        CUTOFF          -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "CUT_OFF")==0) {
        register SCutOff * const ct = &(prf->CutOffData);
        register const size_t JCUT = (size_t) ct->JCUT;
        if (JCUT >= MAXC) goto TooManyCutOff;
	register SCutOffItem * const restrict ctItem = &(ct->Values[JCUT]);

        for (size_t keys=0; keys<AnalyzedLine.counter;++keys) {
          if        (strcmp(AnalyzedLine.keywords[keys], "LEVEL")==0) {
            register const int itmp = atoi(AnalyzedLine.values[keys]);
            ctItem->MCLE = itmp;
            if (itmp == 0 ) LZCO = true;
          } else if (strcmp(AnalyzedLine.keywords[keys], "SCORE")==0) {
            ctItem->ICUT = atoi(AnalyzedLine.values[keys]);
          } else if (strcmp(AnalyzedLine.keywords[keys], "H_SCORE")==0) {
            ctItem->HCUT = (unsigned int) atoi(AnalyzedLine.values[keys]);
          } else if (strcmp(AnalyzedLine.keywords[keys], "TEXT")==0) {
	    // Move to the first quote
	    const char * pos = AnalyzedLine.values[keys];
            const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
	    pos = GoUptoSymbol(AnalyzedLine.values[keys], '\'', MaxMemory);
	    // Copy up to the second quote
	    CopyUptoSymbol(ctItem->CCUT, pos, '\'', 32);
          } else if (strcmp(AnalyzedLine.keywords[keys], "N_SCORE")==0) {
            const char * pos = AnalyzedLine.values[keys];
            const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
            const size_t count = 1 + CountAndReplaceSymbol(AnalyzedLine.values[keys], ',', MaxMemory);
            if (count > MAXN) goto TooManyScores;
            
            for (size_t i=0; i<count; ++i) {
              char * newpos;
              ctItem->RCUT[i] = strtof(pos, &newpos);
	      if (newpos == pos || newpos == 0) {
		fputs("Unable to read N_SCORE values\n" , stderr);
		return 1;
	      }
              while(*newpos != '\0') ++newpos;
              pos = newpos + 1; 
            }
            ctItem->JCNM = (int) count;
          } else if (strcmp(AnalyzedLine.keywords[keys], "MODE")==0) {
            const char * pos = AnalyzedLine.values[keys];
            const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
            const size_t count = 1 + CountAndReplaceSymbol(AnalyzedLine.values[keys], ',', MaxMemory);
            if (count > MAXN) goto TooManyModes;

            for (size_t i=0; i<count; ++i) {
              char * newpos;
              ctItem->MCUT[i] = (int) strtol(pos, &newpos, 10);
              while(*newpos != '\0') ++newpos;
              pos = newpos + 1;
            }
          }
        }
        ct->JCUT = JCUT + 1;
      }
      /* -----------------------       DEFAULT          -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "DEFAULT")==0) {
        if (AnalyzedLine.counter == 0) {
          InitializeDefault(&DefaultScores, &DefaultMatchSymbol, &DefaultInsertionSymbol);
        }
         else {
          for (size_t keys=0; keys<AnalyzedLine.counter;++keys) {
            const char * key = AnalyzedLine.keywords[keys];
				    if ( key[0] == 'S' && key[1] == 'Y') {
				      if (key[3] == 'M') {
								DefaultMatchSymbol = AnalyzedLine.values[keys][1];
				      } 
				      else if (key[3] == 'I') {
								DefaultInsertionSymbol = AnalyzedLine.values[keys][1];
				      }
				    }
            else if ( key[0] == 'M' && key[1] == '0' ) {
              if ( ReadScore(AnalyzedLine.values[keys], &(DefaultScores.Match.Alphabet[0])) != 0 ) goto ReadError;
            } 
            else if ( key[0] == 'M' && key[1] == '\0' ) {
              const char * pos = AnalyzedLine.values[keys];
              const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
              const size_t count = CountAndReplaceSymbol(AnalyzedLine.values[keys], ',', MaxMemory) + 1;
              if (count == 1) {
                StoredIntegerFormat data;
                if ( ReadScore(AnalyzedLine.values[keys], &data) != 0 ) goto ReadError;
                for (size_t i=0; i<Alphabet_Length; ++i) DefaultScores.Match.Alphabet[i+1] = data;
              } 
              else {
                if ( ReadScores(AnalyzedLine.values[keys], &(DefaultScores.Match.Alphabet[1]), count) != 0 ) goto ReadError;
              }
            }
            else if ( key[0] == 'D' && key[1] == '\0' ) {
              if ( ReadScore(AnalyzedLine.values[keys], &(DefaultScores.Match.Alphabet[Alphabet_Length+1])) != 0 ) goto ReadError;
            }
            else {
              StoredIntegerFormat *ptr;
              const size_t type = GetInsertionMemory(key, &DefaultScores.Insertion, &ptr);
              switch (type) {
                case (0):
                  if ( ReadScore(AnalyzedLine.values[keys], ptr) != 0 ) goto ReadError;
                  break;
                case (1):
                  {
                    const char * pos = AnalyzedLine.values[keys];
                    const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
                    const size_t count = CountAndReplaceSymbol(AnalyzedLine.values[keys], ',', MaxMemory) + 1;
                    if (count == 1) {
                      StoredIntegerFormat data;
                      if ( ReadScore(AnalyzedLine.values[keys], &data) != 0 ) goto ReadError;
                      for (size_t i=0; i<Alphabet_Length; ++i) ptr[i] = data;
                    }
                    else {
                      if ( ReadScores(AnalyzedLine.values[keys], ptr, count) != 0 ) goto ReadError;
                    }
                  }
                  break;
                default:
                  goto UnknownKey;
              }
            }
          }
        }
      }
      /* -----------------------        INSERTIONS      -------------------------------*/
      else if (AnalyzedLine.subcommand[0] == 'I') {
				if (MultipleLine) {
					--InsertionCounter;
					PreviousInsertionProfile(&WorkingScores.Insertion);
				}
				else {
					/* Copy default values */
					memcpy(WorkingScores.Insertion.Alphabet,    DefaultScores.Insertion.Alphabet,
								 (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
					memcpy(WorkingScores.Insertion.Boundaries,  DefaultScores.Insertion.Boundaries,
								 (INSERTION_BOUNDARIES_SIZE)*sizeof(StoredIntegerFormat));
					memcpy(WorkingScores.Insertion.Transitions, DefaultScores.Insertion.Transitions,
								 sizeof(TransitionScores));
				}
				/* Check whether there is an implicit /M:, if so insert it */
				if (MatchCounter < InsertionCounter + 1) {
//           fprintf(stderr, "Line %lu : implicit Match : %lu < %lu\n", Line, MatchCounter, InsertionCounter);
#ifdef _VERBOSE_
					if (out_profile) {
						if (UseColor)
							puts(" \e[0;32m IMPLICIT INSERTION OF M\e[0m\n");
						else
							puts("  IMPLICIT INSERTION OF M\n");
					}
#endif
					*ProfileSequence = DefaultMatchSymbol; //ProfileSequence[-1];
					++ProfileSequence;
					memcpy(WorkingScores.Match.Alphabet, DefaultScores.Match.Alphabet, (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
          NextMatchProfile(&WorkingScores.Match);
          ++MatchCounter;
        }

        /* Read the scores */
        for (size_t keys=0; keys<AnalyzedLine.counter;++keys) {
          const char * key = AnalyzedLine.keywords[keys];
					if ( key[0] == 'S' && key[1] == 'Y' ) continue; 
					StoredIntegerFormat *ptr;
					const size_t type = GetInsertionMemory(key, &WorkingScores.Insertion, &ptr);
					switch (type) {
						case (0):
							if ( ReadScore(AnalyzedLine.values[keys], ptr) != 0 ) goto ReadError;
							break;
						case (1):
							{
								const char * pos = AnalyzedLine.values[keys];
								const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
								const size_t count = CountAndReplaceSymbol(AnalyzedLine.values[keys], ',', MaxMemory) + 1;
								if (count == 1) {
									StoredIntegerFormat data;
									if ( ReadScore(AnalyzedLine.values[keys], &data) != 0 ) goto ReadError;
									for (size_t i=0; i<Alphabet_Length; ++i) ptr[i] = data;
								} else {
									if ( ReadScores(AnalyzedLine.values[keys], ptr, count) != 0 ) goto ReadError;
								}
							}
							break;
						default:
							goto UnknownKey;
					}
				}

        /* Increment I counter */
        NextInsertionProfile(&WorkingScores.Insertion);
        ++InsertionCounter;
      }
      /* -----------------------         MATCHES       --------------------------------*/
      else if (AnalyzedLine.subcommand[0] == 'M') {
				if (MultipleLine) {
					--MatchCounter;
					PreviousMatchProfile(&WorkingScores.Match);
				} else {
					/* Copy default values */
					memcpy(WorkingScores.Match.Alphabet, DefaultScores.Match.Alphabet, (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
				}
        /* Check whether there is an implicit /I:, if so insert it */
        if (InsertionCounter < MatchCounter ) {
//           fprintf(stderr, "Line %lu : implicit Insertion : %lu <= %lu\n", Line, InsertionCounter, MatchCounter);
#ifdef _VERBOSE_
					if (out_profile) {
						if (UseColor) 
							puts(" \e[0;32m IMPLICIT INSERTION OF I\e[0m\n");
						else
							puts("  IMPLICIT INSERTION OF I\n");
					}
#endif
					memcpy(WorkingScores.Insertion.Alphabet,    DefaultScores.Insertion.Alphabet,
					       (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
					memcpy(WorkingScores.Insertion.Boundaries,  DefaultScores.Insertion.Boundaries,
					       (INSERTION_BOUNDARIES_SIZE)*sizeof(StoredIntegerFormat));
					memcpy(WorkingScores.Insertion.Transitions, DefaultScores.Insertion.Transitions,
					       sizeof(TransitionScores));

//        CopyPreviousInsertionProfile(&WorkingScores.Insertion);
	  
					NextInsertionProfile(&WorkingScores.Insertion);
          ++InsertionCounter;
        }

        /* Read the scores */
				_Bool HasGivenSequence = false;
        for (size_t keys=0; keys<AnalyzedLine.counter;++keys) {
          const char * key = AnalyzedLine.keywords[keys];
					if ( key[0] == 'S' && key[1] == 'Y' ) {
						const char * AminoAcid = AnalyzedLine.values[keys];
						while (*AminoAcid == '\'') AminoAcid++;
						*ProfileSequence = *AminoAcid;
						++ProfileSequence;
						HasGivenSequence = true;
					}
					else if ( key[0] == 'M' && key[1] == '0' ) {
						if ( ReadScore(AnalyzedLine.values[keys], &(WorkingScores.Match.Alphabet[0])) != 0 ) goto ReadError;
					}
					else if ( key[0] == 'M' && key[1] == '\0' ) {
						const char * pos = AnalyzedLine.values[keys];
						const uintptr_t MaxMemory = (uintptr_t) pos + strlen(pos);
						const size_t count = CountAndReplaceSymbol(AnalyzedLine.values[keys], ',', MaxMemory) + 1;
						if (count == 1) {
							StoredIntegerFormat data;
							if ( ReadScore(AnalyzedLine.values[keys], &data) != 0 ) goto ReadError;
							for (size_t i=0; i<Alphabet_Length; ++i) WorkingScores.Match.Alphabet[i+1] = data;
						} else {
							if ( ReadScores(AnalyzedLine.values[keys], &(WorkingScores.Match.Alphabet[1]), count) != 0 ) goto ReadError;
						}
					}
					else if ( key[0] == 'D' && key[1] == '\0' ) {
						if ( ReadScore(AnalyzedLine.values[keys], &WorkingScores.Match.Alphabet[Alphabet_Length+1]) != 0 ) goto ReadError;
					}
        }
        if (!HasGivenSequence) {
					*ProfileSequence = DefaultMatchSymbol;
					++ProfileSequence;
				}

        /* Increment M counter */
        NextMatchProfile(&WorkingScores.Match);
        ++MatchCounter;
      }
      /* -----------------------       UNKNOWN          -------------------------------*/
      else {
#ifdef _VERBOSE_
	if (out_profile) {
	  if (UseColor)
	    puts("\e[33;41m\t!!! UNTREATED COMMAND !!!\e[0m");
	  else
	    puts("\t!!! UNTREATED COMMAND !!!");
	}
#endif
      }
    }
  }
  
  /* 
   * End of file reached without end of profile keyword. This typically happens when profile
   * has empty lines after the last end of profile keyword.
   * Let us check this and return here with special value -1.
   */ 
  if (Length < 1) return -1;
  
END_OF_PROFILE:
  if (prf->Type == PF_MATRIX) {
    /* Insert possible missing I when profile starts with M */
    if (InsertionCounter == Length-1 ) {
      /* Copy default values */
      memcpy(WorkingScores.Insertion.Alphabet,    DefaultScores.Insertion.Alphabet, (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
      memcpy(WorkingScores.Insertion.Boundaries,  DefaultScores.Insertion.Boundaries,  (INSERTION_BOUNDARIES_SIZE)*sizeof(StoredIntegerFormat));
      memcpy(WorkingScores.Insertion.Transitions, DefaultScores.Insertion.Transitions, sizeof(TransitionScores));
      ++InsertionCounter;
    }

    /* Insert possible missing M line at last */
    if (MatchCounter == Length) {
      /* Copy default values */
      memcpy(WorkingScores.Match.Alphabet, DefaultScores.Match.Alphabet, (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
    }
    
    /* CHECK CONSISTENCY */
    if (MatchCounter != Length) {
      if ( MatchCounter < Length )
				fprintf(stderr, "There is not enough match lines %i <> %lu\n", MatchCounter, Length );
      else
				fprintf(stderr, "There is too many match lines within profile, %i <> %lu\n", MatchCounter, Length );
      return -1;
    }
    if (InsertionCounter != Length) {
      if ( InsertionCounter < Length )
				fprintf(stderr, "There is not enough insertion lines %i <> %lu\n", InsertionCounter, Length );
      else
				fprintf(stderr, "There is too many insertion lines within profile, %i <> %lu\n", InsertionCounter, Length );
      return -1;
    }
#ifndef _BEST_IS_NEGATIVE_
    if (LZCO != true || prf->CutOffData.JCUT == 0) {
	fputs("No level 0 CUT-OFF data block in profile.\n", stderr);
	return 1;
    }
#endif
    if (Length < 1) {
      fputs("Unexpected end of profile. Profile has zero length.\n", stderr);
      return -1;
    }

    if (prf->DisjointData.NDIP[0] <= 0 || prf->DisjointData.NDIP[0] > Length) {
      fprintf(stderr,
	      "Warning: Disjointness parameter 1 (%i) out of bound. Parameter set to acceptable value.\n",
	      prf->DisjointData.NDIP[0]);
      prf->DisjointData.NDIP[0] = 1;
    }
    if (prf->DisjointData.NDIP[1] <= 0 || prf->DisjointData.NDIP[1] > Length) {
      fprintf(stderr,
	      "Warning: Disjointness parameter 2 (%i) out of bound. Parameter set to acceptable value (%lu).\n",
	      prf->DisjointData.NDIP[1], Length);
      prf->DisjointData.NDIP[1] = (int) Length;    
    }

    if (prf->DisjointData.NDIP[1] < prf->DisjointData.NDIP[0]) {
      const int tmp = prf->DisjointData.NDIP[1];
      prf->DisjointData.NDIP[1] = prf->DisjointData.NDIP[0];
      prf->DisjointData.NDIP[0] = tmp;
    }
    
    /* 
		 * Have to do this no earlier than here, otherwise the consitency check and modifications are
		 * not brought correctly, nasty bug 
		 */
    if (SetExtraTable) 
    {
      const int ret = PrepareExtraTable(prf);
      if (ret) return ret;
    }
    
  }
  
  *LineOffset = Line;
  return 0;
  
  /*
   * ERRORS
   */
  
// MissingSymbol:
//    fprintf(stderr, "Missing symbol %c at line %lu of %s\n\tLine: %s\n",Symbol,Line, FileName, prfLine);
//    return 1;
   
MissingType:
   if (prf->Identification[0] != '\0') {
    fprintf(stderr, "Missing MATRIX or PATTERN keyword at line %lu of profile %s\n\tLine: %s\n", Line, prf->Identification, currentLine);
   } else {
    fprintf(stderr, "Missing MATRIX or PATTERN keyword at line %lu of %s\n\tLine: %s\n", Line, FileName, currentLine);
   }
  goto FIN;
   
AlphabetSizeTooLarge:
   fprintf(stderr, "Alphabet size exceeds hard defined size: %u > %lu\n", ALPHABET_SIZE, prf->Alphabet_Length);
  goto FIN;

NormalizationError:
  fprintf(stderr, "Error within normalization section at line %lu of %s\nFUNCTION value matches none of the following.\n",
          Line, FileName);
  for (int i=0; i<KNOR; ++i) fprintf(stderr, "%s ", prf->NormalizationData.CNOR[i]);
  fputs("\n", stderr);
  
 goto FIN;
    
TooManyNormalization:
  fprintf(stderr, "Too many normalization parameters at line %lu of %s\n\tMaximum is %i.\n", Line, FileName, 0);
  goto FIN;

TooManyCutOff:
  fprintf(stderr, "Too many cutoffs parameters at line %lu of %s\n\tMaximum is %i.\n", Line, FileName, 0);
  goto FIN;

TooManyModes:
   fprintf(stderr, "Too many modes parameters at line %lu of %s\n\tMaximum is %i.\n", Line, FileName, 0);
  goto FIN;

TooManyScores:
   fprintf(stderr, "Too many scores at line %lu of %s\n\tMaximum is %i.\n", Line, FileName, 0);
  goto FIN;

ReadError:
   fprintf(stderr, "Error reading one parameter at line %lu of %s\n\tLine: %s\n", Line, FileName, currentLine);
  goto FIN;

UnknownKey:
   fprintf(stderr, "Unknown keyword found at line %lu of %s\n\tLine: %s\n", Line, FileName, currentLine);
  goto FIN;

AllocationError:
  fputs("Unable to allocate sufficient memory\n", stderr);
  if (prf->Sequence) free(prf->Sequence);
 goto FIN;
   
WrongType:
  fprintf(stderr, "MA keyword should not appear in PATTERN profile at line %lu of %s\n", Line, FileName);
  
FIN:
    return 1;
}
//---------------------------------------------------------------

/*
 * Exported functions to work on profile structure.
 */
#ifndef _BEST_IS_NEGATIVE_
int PrepareExtraTable(struct Profile * const prf)
{
  /*
   * Operate on the Insertion Score matrix to place boundaries scores within the structure at
   * location reserved for them.
   * NOTE: This was initially the IIPX matrix filled up within pfsearch.
   */
  /*
     IIPX( XM,I1) = MAX(MLOW,IIPP( B1,I1) + IIPP( BM,I1))
     IIPX( XI,I1) = MAX(MLOW,IIPP( B1,I1) + IIPP( BI,I1))
     IIPX( XD,I1) = MAX(MLOW,IIPP( B1,I1) + IIPP( BD,I1))

     IIPX( YM,I1) = MAX(MLOW,IIPP( B0,I1) + IIPP( BM,I1))
     IIPX( YI,I1) = MAX(MLOW,IIPP( B0,I1) + IIPP( BI,I1))
     IIPX( YD,I1) = MAX(MLOW,IIPP( B0,I1) + IIPP( BD,I1))

     IIPX( MX,I1) = MAX(MLOW,IIPP( E1,I1) + IIPP( ME,I1))
     IIPX( IX,I1) = MAX(MLOW,IIPP( E1,I1) + IIPP( IE,I1))
     IIPX( DX,I1) = MAX(MLOW,IIPP( E1,I1) + IIPP( DE,I1))

     IIPX( MY,I1) = MAX(MLOW,IIPP( E0,I1) + IIPP( ME,I1))
     IIPX( IY,I1) = MAX(MLOW,IIPP( E0,I1) + IIPP( IE,I1))
     IIPX( DY,I1) = MAX(MLOW,IIPP( E0,I1) + IIPP( DE,I1))
  */
  struct Profile * restrict lprf = prf;
  do {
    const StoredIntegerFormat * restrict const InsertionBoundaries = lprf->Scores.Insertion.Boundaries;
    StoredIntegerFormat * restrict const InsertionScores           = lprf->Scores.Insertion.Transitions->Element;
    ScoreTuple * restrict const FirstScores                        = lprf->Scores.Insertion.FirstSequenceProtein;
    ScoreTuple * restrict const LastScores                         = lprf->Scores.Insertion.LastSequenceProtein;
    const int MLOW = NLOW/4*3;
    const size_t Length = lprf->Length;
  
    if ( lprf->isCircular && lprf->CompleteCycleOnly) {
      for (size_t i=0; i<=Length; ++i) {
				register const size_t offset  = INSERTION_TRANSITIONS_SIZE*i;
				register const size_t Boffset = INSERTION_BOUNDARIES_SIZE*i;
				/*
				  * WARNING: POTENTIAL ISSUE WITH _MX or _MY
				  */
				CHECK_AND_SET(InsertionScores[offset + _XM], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(InsertionScores[offset + _XI], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(InsertionScores[offset + _XD], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BD]);
				
				// Minimize dummy element 
				InsertionScores[offset + _DUMMY] = NLOW;
				CHECK_AND_SET(InsertionScores[offset + _MX], InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(InsertionScores[offset + _IX], InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(InsertionScores[offset + _DX], InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _DE]);
				
				CHECK_AND_SET(FirstScores[i].To[MATCH]    , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(FirstScores[i].To[INSERTION], InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(FirstScores[i].To[DELETION] , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BD]);

				CHECK_AND_SET(LastScores[i].From[MATCH]    , InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(LastScores[i].From[INSERTION], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(LastScores[i].From[DELETION] , InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _DE]);
      }
    }
    else {
      const size_t NDIP1 = lprf->DisjointData.NDIP[0];
      const size_t NDIP2 = lprf->DisjointData.NDIP[1];
      
      for (size_t i=0; i<NDIP1; ++i) {
				register const size_t offset  = INSERTION_TRANSITIONS_SIZE*i;
				register const size_t Boffset = INSERTION_BOUNDARIES_SIZE*i;

				CHECK_AND_SET(InsertionScores[offset + _XM], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(InsertionScores[offset + _XI], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(InsertionScores[offset + _XD], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BD]);
				
				// Minimize dummy element 
				InsertionScores[offset + _DUMMY] = NLOW;
				CHECK_AND_SET(InsertionScores[offset + _MX], NLOW, InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(InsertionScores[offset + _IX], NLOW, InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(InsertionScores[offset + _DX], NLOW, InsertionBoundaries[Boffset + _DE]);

				CHECK_AND_SET(FirstScores[i].To[MATCH]    , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(FirstScores[i].To[INSERTION], InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(FirstScores[i].To[DELETION] , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BD]);

				CHECK_AND_SET(LastScores[i].From[MATCH]    , NLOW, InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(LastScores[i].From[INSERTION], NLOW, InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(LastScores[i].From[DELETION] , NLOW, InsertionBoundaries[Boffset + _DE]);
      }
      for (size_t i=NDIP1; i<NDIP2; ++i) {
				register const size_t offset  = INSERTION_TRANSITIONS_SIZE*i;
				register const size_t Boffset = INSERTION_BOUNDARIES_SIZE*i;

				CHECK_AND_SET(InsertionScores[offset + _XM], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(InsertionScores[offset + _XI], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(InsertionScores[offset + _XD], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BD]);
				
				// Minimize dummy element 
				InsertionScores[offset + _DUMMY] = NLOW;
				CHECK_AND_SET(InsertionScores[offset + _MX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(InsertionScores[offset + _IX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(InsertionScores[offset + _DX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _DE]);

				CHECK_AND_SET(FirstScores[i].To[MATCH]    , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(FirstScores[i].To[INSERTION], InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(FirstScores[i].To[DELETION] , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BD]);

				CHECK_AND_SET(LastScores[i].From[MATCH]    , InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(LastScores[i].From[INSERTION], InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(LastScores[i].From[DELETION] , InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _DE]);
      }
      for (size_t i=NDIP2; i<=Length; ++i) {
				register const size_t offset  = INSERTION_TRANSITIONS_SIZE*i;
				register const size_t Boffset = INSERTION_BOUNDARIES_SIZE*i;

				CHECK_AND_SET(InsertionScores[offset + _XM], NLOW, InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(InsertionScores[offset + _XI], NLOW, InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(InsertionScores[offset + _XD], NLOW, InsertionBoundaries[Boffset + _BD]);
				
				// Minimize dummy element 
				InsertionScores[offset + _DUMMY] = NLOW;
				CHECK_AND_SET(InsertionScores[offset + _MX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(InsertionScores[offset + _IX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(InsertionScores[offset + _DX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _DE]);

				CHECK_AND_SET(FirstScores[i].To[MATCH]    , NLOW, InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(FirstScores[i].To[INSERTION], NLOW, InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(FirstScores[i].To[DELETION] , NLOW, InsertionBoundaries[Boffset + _BD]);

				CHECK_AND_SET(LastScores[i].From[MATCH]    , InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(LastScores[i].From[INSERTION], InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(LastScores[i].From[DELETION] , InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _DE]);
      }
    }
    lprf = lprf->next;
  } while (lprf != NULL);
  return 0;
}
#else
int PrepareExtraTable(struct Profile * const prf)
{
  /*
   * Operate on the Insertion Score matrix to place boundaries scores within the structure at
   * location reserved for them.
   * NOTE: This was initially the IIPX matrix filled up within pfsearch.
   */
  struct Profile * restrict lprf = prf;
  do {
    const StoredIntegerFormat * restrict const InsertionBoundaries = lprf->Scores.Insertion.Boundaries;
    StoredIntegerFormat * restrict const InsertionScores           = lprf->Scores.Insertion.Transitions->Element;
    ScoreTuple * restrict const FirstScores                        = lprf->Scores.Insertion.FirstSequenceProtein;
    ScoreTuple * restrict const LastScores                         = lprf->Scores.Insertion.LastSequenceProtein;
    const int MLOW = NLOW/4*3;
    const size_t Length = lprf->Length;
  
    if ( lprf->isCircular) {
      fputs("No coverage on circular profile", stderr);
      return 1;
    } 
    else {
      const size_t NDIP1 = lprf->DisjointData.NDIP[0];
      const size_t NDIP2 = lprf->DisjointData.NDIP[1];
      
      for (size_t i=0; i<NDIP1; ++i) {
				register const size_t offset  = INSERTION_TRANSITIONS_SIZE*i;
				register const size_t Boffset = INSERTION_BOUNDARIES_SIZE*i;

				CHECK_AND_SET(InsertionScores[offset + _XM], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(InsertionScores[offset + _XI], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(InsertionScores[offset + _XD], InsertionBoundaries[Boffset + _B1], InsertionBoundaries[Boffset + _BD]);
				
				// Minimize dummy element 
				InsertionScores[offset + _DUMMY] = NLOW;
				CHECK_AND_SET(InsertionScores[offset + _MX], NLOW, InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(InsertionScores[offset + _IX], NLOW, InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(InsertionScores[offset + _DX], NLOW, InsertionBoundaries[Boffset + _DE]);

				CHECK_AND_SET(FirstScores[i].To[MATCH]    , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(FirstScores[i].To[INSERTION], InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(FirstScores[i].To[DELETION] , InsertionBoundaries[Boffset + _B0], InsertionBoundaries[Boffset + _BD]);

				CHECK_AND_SET(LastScores[i].From[MATCH]    , NLOW, InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(LastScores[i].From[INSERTION], NLOW, InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(LastScores[i].From[DELETION] , NLOW, InsertionBoundaries[Boffset + _DE]);
      }
      for (size_t i=NDIP1; i<NDIP2; ++i) {
				register const size_t offset  = INSERTION_TRANSITIONS_SIZE*i;
				register const size_t Boffset = INSERTION_BOUNDARIES_SIZE*i;

				CHECK_AND_SET(InsertionScores[offset + _XM], NLOW, InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(InsertionScores[offset + _XI], NLOW, InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(InsertionScores[offset + _XD], NLOW, InsertionBoundaries[Boffset + _BD]);
				
				// Minimize dummy element 
				InsertionScores[offset + _DUMMY] = NLOW;
				CHECK_AND_SET(InsertionScores[offset + _MX], NLOW, InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(InsertionScores[offset + _IX], NLOW, InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(InsertionScores[offset + _DX], NLOW, InsertionBoundaries[Boffset + _DE]);

				CHECK_AND_SET(FirstScores[i].To[MATCH]    , NLOW, InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(FirstScores[i].To[INSERTION], NLOW, InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(FirstScores[i].To[DELETION] , NLOW, InsertionBoundaries[Boffset + _BD]);

				CHECK_AND_SET(LastScores[i].From[MATCH]    , NLOW, InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(LastScores[i].From[INSERTION], NLOW, InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(LastScores[i].From[DELETION] , NLOW, InsertionBoundaries[Boffset + _DE]);
      }
      for (size_t i=NDIP2; i<=Length; ++i) {
				register const size_t offset  = INSERTION_TRANSITIONS_SIZE*i;
				register const size_t Boffset = INSERTION_BOUNDARIES_SIZE*i;

				CHECK_AND_SET(InsertionScores[offset + _XM], NLOW, InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(InsertionScores[offset + _XI], NLOW, InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(InsertionScores[offset + _XD], NLOW, InsertionBoundaries[Boffset + _BD]);
				
				// Minimize dummy element 
				InsertionScores[offset + _DUMMY] = NLOW;
				CHECK_AND_SET(InsertionScores[offset + _MX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(InsertionScores[offset + _IX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(InsertionScores[offset + _DX], InsertionBoundaries[Boffset + _E1], InsertionBoundaries[Boffset + _DE]);

				CHECK_AND_SET(FirstScores[i].To[MATCH]    , NLOW, InsertionBoundaries[Boffset + _BM]);
				CHECK_AND_SET(FirstScores[i].To[INSERTION], NLOW, InsertionBoundaries[Boffset + _BI]);
				CHECK_AND_SET(FirstScores[i].To[DELETION] , NLOW, InsertionBoundaries[Boffset + _BD]);

				CHECK_AND_SET(LastScores[i].From[MATCH]    , InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _ME]);
				CHECK_AND_SET(LastScores[i].From[INSERTION], InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _IE]);
				CHECK_AND_SET(LastScores[i].From[DELETION] , InsertionBoundaries[Boffset + _E0], InsertionBoundaries[Boffset + _DE]);
      }
    }
    lprf = lprf->next;
  } while (lprf != NULL);
  return 0;
}
//---------------------------------------------------------------
#endif

struct Profile * ReverseProfile(const struct Profile * const restrict inprf)
{
  /* Allocates memory */
  struct Profile * const restrict outprf = _mm_malloc(sizeof(struct Profile), 16);
  if (outprf == NULL) return NULL;
  register const size_t prfLength = inprf->Length;
  
  /* Copy parameters */
#if 1
  memcpy(outprf, inprf, sizeof(struct Profile));
#else
  memcpy(outprf->Identification, inprf->Identification, 64*sizeof(char));
  memcpy(outprf->AC_Number, inprf->AC_Number, 64*sizeof(char));
  memcpy(outprf->Date, inprf->Date, 128*sizeof(char));
  memcpy(outprf->Description, inprf->Description, 256*sizeof(char));
  outprf->Alphabet_Length = inprf->Alphabet_Length;
  outprf->isCircular      = inprf->isCircular;
  outprf->Length          = prfLength;
  memcpy(outprf->Alphabet_Mapping, inprf->Alphabet_Mapping, (ALPHABET_SIZE+1)*sizeof(char));
  memcpy(outprf->CABC, inprf->CABC, (ALPHABET_SIZE+2)*sizeof(char));
  outprf->Level = inprf->Level;
  outprf->Mode  = inprf->Mode;
#endif
//   fputs("Testing with NDIP offset by 1 to circumvent null value\n",stderr);
  outprf->DisjointData.NDIP[0] = (int) prfLength - inprf->DisjointData.NDIP[1] + 1;
  outprf->DisjointData.NDIP[1] = (int) prfLength - inprf->DisjointData.NDIP[0] + 1;
//   fprintf(stderr, "NDIP %i %i -> %i %i\n", inprf->DisjointData.NDIP[0], inprf->DisjointData.NDIP[1],
// 					   outprf->DisjointData.NDIP[0], outprf->DisjointData.NDIP[1]);
  
  /* Allocates memory for profile sequence */
  char * Sequence = malloc(prfLength*sizeof(char));
  if (Sequence == NULL) {
    _mm_free(outprf);
    return NULL;
  }
  outprf->Sequence = Sequence;
  for (size_t iprf=0; iprf<prfLength; ++iprf) {
    Sequence[iprf] = inprf->Sequence[prfLength-1-iprf];
  }
  
  /* Allocates memory for score tables*/
  if (AllocateScores(&(outprf->Scores), inprf->Alphabet_Length, prfLength) != 0) {
    free(Sequence);
    _mm_free(outprf);
    return NULL;
  }
  
  /* Transfer match and insertion alphabets */
  {
    register const size_t AlignedStep = inprf->Scores.Insertion.AlignStep;
    register const StoredIntegerFormat * restrict inMatch = &inprf->Scores.Match.Alphabet[AlignedStep*(inprf->Length-1)];
    register StoredIntegerFormat * restrict outMatch      = outprf->Scores.Match.Alphabet;
    for (size_t iprf=0; iprf<prfLength; ++iprf) {
      memcpy(outMatch,inMatch, AlignedStep*sizeof(StoredIntegerFormat));
      outMatch     += AlignedStep;
      inMatch      -= AlignedStep;
    }
    
    register const StoredIntegerFormat * restrict inInsertion = &inprf->Scores.Insertion.Alphabet[AlignedStep*(inprf->Length)];
    register StoredIntegerFormat * restrict outInsertion      = outprf->Scores.Insertion.Alphabet;
    for (size_t iprf=0; iprf<=prfLength; ++iprf) {
      memcpy(outInsertion, inInsertion, AlignedStep*sizeof(StoredIntegerFormat));
      outInsertion += AlignedStep;
      inInsertion  -= AlignedStep;
    }
  }
  
  /* Transfer boundaries */
  register const StoredIntegerFormat (*inBoundaries)[INSERTION_BOUNDARIES_SIZE] = (const StoredIntegerFormat (*)[INSERTION_BOUNDARIES_SIZE]) inprf->Scores.Insertion.Boundaries;
  register StoredIntegerFormat (*outBoundaries)[INSERTION_BOUNDARIES_SIZE]      = (StoredIntegerFormat (*)[INSERTION_BOUNDARIES_SIZE]) outprf->Scores.Insertion.Boundaries;
  for (size_t iprf=0; iprf<=prfLength; ++iprf) {
    outBoundaries[iprf][_B0] = inBoundaries[prfLength-iprf][_E0];
    outBoundaries[iprf][_B1] = inBoundaries[prfLength-iprf][_E1];
    outBoundaries[iprf][_E0] = inBoundaries[prfLength-iprf][_B0];
    outBoundaries[iprf][_E1] = inBoundaries[prfLength-iprf][_B1];
    outBoundaries[iprf][_BI] = inBoundaries[prfLength-iprf][_IE];
    outBoundaries[iprf][_IE] = inBoundaries[prfLength-iprf][_BI];
    outBoundaries[iprf][_BM] = inBoundaries[prfLength-iprf][_ME];
    outBoundaries[iprf][_ME] = inBoundaries[prfLength-iprf][_BM];
    outBoundaries[iprf][_BD] = inBoundaries[prfLength-iprf][_DE];
    outBoundaries[iprf][_DE] = inBoundaries[prfLength-iprf][_BD];
    outBoundaries[iprf][_BE] = inBoundaries[prfLength-iprf][_BE];
	}
  
  /* Transfer transitions */
  register const TransitionScores *inTransitions = (TransitionScores*) inprf->Scores.Insertion.Transitions;
  register TransitionScores *outTransitions      = (TransitionScores*) outprf->Scores.Insertion.Transitions;
  for (size_t iprf=0; iprf<=prfLength; ++iprf) {
     outTransitions[iprf].Element[_DD] = inTransitions[prfLength-iprf].Element[_DD];
     outTransitions[iprf].Element[_DI] = inTransitions[prfLength-iprf].Element[_ID];
     outTransitions[iprf].Element[_DM] = inTransitions[prfLength-iprf].Element[_MD];
// 		 outTransitions[iprf].Element[_DX] = inTransitions[prfLength-iprf].Element[_XD];
     
     outTransitions[iprf].Element[_MD] = inTransitions[prfLength-iprf].Element[_DM];
     outTransitions[iprf].Element[_MM] = inTransitions[prfLength-iprf].Element[_MM];
     outTransitions[iprf].Element[_MI] = inTransitions[prfLength-iprf].Element[_IM];
// 		 outTransitions[iprf].Element[_MX] = inTransitions[prfLength-iprf].Element[_XM];
     
     outTransitions[iprf].Element[_ID] = inTransitions[prfLength-iprf].Element[_DI];
     outTransitions[iprf].Element[_IM] = inTransitions[prfLength-iprf].Element[_MI];
     outTransitions[iprf].Element[_II] = inTransitions[prfLength-iprf].Element[_II];
// 		 outTransitions[iprf].Element[_IX] = inTransitions[prfLength-iprf].Element[_XI];
		 
// 		 outTransitions[iprf].Element[_XD] = inTransitions[prfLength-iprf].Element[_DX];
//      outTransitions[iprf].Element[_XM] = inTransitions[prfLength-iprf].Element[_MX];
//      outTransitions[iprf].Element[_XI] = inTransitions[prfLength-iprf].Element[_IX];
  }
  outprf->isReversed = true;
  return outprf;
}
//---------------------------------------------------------------

/*
 * ReadProfile returns the number of profile read, or -1 on error
 */
int ReadProfile(const char * const restrict FileName, struct Profile * const prf,
                const _Bool SetExtraTable, const _Bool CompleteCycleOnly)
{
	int nprf = 0;
	size_t Line = 0;
	struct Profile * newPrf = prf;
	
	/*
	* Try to open the file, upon failure emmit error
	*/
	FILE* prfStream = fopen(FileName,"r");
	if (prfStream == NULL) {
		return -1;
	}

	/*
	* Read all internal data 
	*/
	
	/* Clean Profile structure */
	memset(newPrf, 0, sizeof(struct Profile));
	while (!feof(prfStream)) {   
		/* Read one profile structure at a time*/
		const int res = internalReadProfile(prfStream, newPrf, SetExtraTable, FileName, &Line, CompleteCycleOnly);
		++nprf;
		if (res == 0 ) {
			struct Profile * const tmpPrf = (struct Profile*) _mm_malloc(sizeof(struct Profile), 16);
			if (tmpPrf) {
				/* Clean Profile structure */
				memset(tmpPrf, 0, sizeof(struct Profile));

				tmpPrf->previous = newPrf;
				newPrf->next = tmpPrf;
				newPrf = tmpPrf;
			}
			else {
				nprf = -1;
				break;
			}
		}
		else {
			if (res == -1) {
				/* Clean the last one as its length is null meaning extra empty lines within previous profile */
				struct Profile * const tmpPrf = newPrf->previous;
				tmpPrf->next = NULL;
				_mm_free(newPrf);
				--nprf;
			}
			else {
				nprf = -1;
			}
			break; 
		}
	}

	/*
	 * Close the file 
	 */
	fclose(prfStream);

	return nprf;
}
//---------------------------------------------------------------

int WriteProfile(const char * const restrict FileNameIn, struct Profile * const prf,
		 const char * const restrict AdditionalComments, FILE * const prfStreamOut)
{
  char currentLine[PROFILE_MAX_LINE_SIZE] __attribute__((aligned(16)));
  char currentLineCopy[PROFILE_MAX_LINE_SIZE] __attribute__((aligned(16)));
	char previousSubCommand[SUB_COMMAND_MAX_SIZE];
	
  _Bool NormalizationDone = false, CutOffDone = false;
  struct ProfileLine AnalyzedLine;
  /*
   * Prepare basic structures
   */
  size_t CommentStorageSize = 4096;
  char * restrict CommentStorage = (char*) malloc(CommentStorageSize*sizeof(char));
  if (CommentStorage == NULL) {
      fputs("Unable to allocate memory\n", stderr);
      return(1);
  }
  size_t CommentStorageUsed = 0;
  /*
   * Try to open the file for reading, upon failure emmit error
   */
  FILE* prfStreamIn = fopen(FileNameIn,"r");
  if (prfStreamIn == NULL) {
    return 1;
  }
  
  /*
   * Write down the header data 
   */
  fprintf(prfStreamOut,
	  "ID   %s; MATRIX.\n"
	  "AC   %s;\n"
	  "DT   %s\n"
	  "DE   %s\n",
	  prf->Identification, prf->AC_Number, prf->Date, prf->Description
	 );
  
  /*
   * Follow the initial profile and parse/reuse
   */
  size_t Line = 0;
  int res;
  _Bool MultipleLine;
  memset(previousSubCommand, 0, sizeof(char)*SUB_COMMAND_MAX_SIZE);
  
  while (!feof(prfStreamIn)) {   
    const size_t length = GetLine(prfStreamIn, currentLine, PROFILE_MAX_LINE_SIZE);
    ++Line;
    if (length == 0) continue;

    /* Make a copy */
    memcpy(currentLineCopy, currentLine, length+1);
    currentLineCopy[length] = '\0';
    
    if ((res=AnalyzeLine(currentLine, &AnalyzedLine, &MultipleLine, previousSubCommand)) != 0) {
      fprintf(stderr,"Error %i in analysis at line %lu\n", res, Line);
      return 1;
    }
    
    /* Header or matrix */
    if (!AnalyzedLine.isMatrix)
    {
       if ( (AnalyzedLine.command[0] == 'C' && AnalyzedLine.command[1] == 'C' ) ) {
         if (CommentStorageUsed+length+2 >= CommentStorageSize) {
	    CommentStorageSize += 4096;
	    CommentStorage = realloc(CommentStorage, CommentStorageSize);
	    if (CommentStorage == NULL) {
		fputs("Unable to grow memory size\n",stderr);
		fclose(prfStreamIn);
		return 1;
	    }
	 }
	 memcpy(&CommentStorage[CommentStorageUsed], currentLineCopy, length);
	 CommentStorageUsed += length;
	 CommentStorage[CommentStorageUsed++] = '\n';
	 CommentStorage[CommentStorageUsed] = '\0';
       }
       else if ( (AnalyzedLine.command[0] == 'D' && AnalyzedLine.command[1] == 'R' ) || 
	         (AnalyzedLine.command[0] == 'N' && AnalyzedLine.command[1] == 'R' ) ||
	         (AnalyzedLine.command[0] == '3' && AnalyzedLine.command[1] == 'D' ) ||
	         (AnalyzedLine.command[0] == 'P' && AnalyzedLine.command[1] == 'R' ) ||
	         (AnalyzedLine.command[0] == 'D' && AnalyzedLine.command[1] == 'O' ) ) {
	  fprintf(prfStreamOut,"%s\n",currentLineCopy);
       }
       else if ( (AnalyzedLine.command[0] == '/' && AnalyzedLine.command[1] == '/' ) ) {
	  goto END_OF_PROFILE;
       }
//        else {
// 	fprintf(stderr,"The following unknown line has been avoided\n %s\n",currentLineCopy);
//       }
    }
    else {
      /* ----------------------- GENERAL SPECIFICATIONS -------------------------------*/
      if (strcmp(AnalyzedLine.subcommand, "GENERAL_SPEC")==0) {
        fprintf(prfStreamOut,"%s\n",currentLineCopy);
      }
      /* -----------------------      DISJOINT          -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "DISJOINT")==0) {
	fprintf(prfStreamOut,"%s\n",currentLineCopy);
      }
      /* -----------------------   NORMALIZATION        -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "NORMALIZATION")==0) {
	if (NormalizationDone) continue;
	const int Count = prf->NormalizationData.JNOR;
        if (Count) NormalizationDone = true;
	for (int Norm=0; Norm<Count; ++Norm) {
	  register const SNormalizationItem * const nrmItem = &(prf->NormalizationData.Values[Norm]);
	  const char * Function = prf->NormalizationData.CNOR[nrmItem->MNOR];
	  const int nCoef       = prf->NormalizationData.JNOP[nrmItem->MNOR];
	  const int Mode        = nrmItem->NNOR;
	  const int Priority    = nrmItem->NNPR;
	  const char * Text     = nrmItem->CNTX;
	  
	  fprintf(prfStreamOut, "MA   /NORMALIZATION: MODE=%i; FUNCTION=%s;", Mode, Function);
	  for (int iCoef=1; iCoef<=nCoef; ++iCoef) {
	      fprintf(prfStreamOut, " R%1.1i=%f;", iCoef, nrmItem->RNOP[iCoef-1]);
	  }
	  if (Priority) fprintf(prfStreamOut, " PRIORITY=%i;", Priority);
	  if (strlen(Text)) fprintf(prfStreamOut, " TEXT=\'%s\';", Text);
	  fputc('\n', prfStreamOut);
	}
      }
      /* -----------------------        CUTOFF          -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "CUT_OFF")==0) {
	if (CutOffDone) continue;
	const int Count = prf->CutOffData.JCUT;
	if (Count) CutOffDone = true;
	for (int Cut=0; Cut<Count; ++Cut) {
	  register const SCutOffItem * const cutItem = &(prf->CutOffData.Values[Cut]);
	  const int Level   = cutItem->MCLE;
	  const int RScore  = cutItem->ICUT;
	  const unsigned int HScore  = cutItem->HCUT;
	  const char * Text = cutItem->CCUT;
	  const int nModes  = cutItem->JCNM;
	  
	  fprintf(prfStreamOut, "MA   /CUT_OFF: LEVEL=%i; SCORE=%i;", Level, RScore);
	  if (HScore>0) {
	      fprintf(prfStreamOut, " H_SCORE=%u;", HScore);
	  }
	  if (nModes) {
	    int iMode = 0;
	    fputs(" N_SCORE=", prfStreamOut);
	    do {
	      if (iMode>0) fputc(',', prfStreamOut);
	      fprintf(prfStreamOut, "%.1f", cutItem->RCUT[iMode]);
	    } while (++iMode < nModes);
	    fputs("; MODE=", prfStreamOut);
	    iMode = 0;
	    do {
	      if (iMode>0) fputc(',', prfStreamOut);
	      fprintf(prfStreamOut, "%i", cutItem->MCUT[iMode]);
	    } while (++iMode < nModes);
	    fputc(';', prfStreamOut);
	  }
	  if (strlen(Text)) fprintf(prfStreamOut, " TEXT=\'%s\';", Text);
	  fputc('\n', prfStreamOut);
	}
      }
      /* -----------------------       DEFAULT          -------------------------------*/
      else if (strcmp(AnalyzedLine.subcommand, "DEFAULT")==0) {
	fprintf(prfStreamOut,"%s\n",currentLineCopy);
      }
      /* -----------------------        INSERTIONS      -------------------------------*/
      else if (AnalyzedLine.subcommand[0] == 'I') {
	fprintf(prfStreamOut,"%s\n",currentLineCopy);
      }
      /* -----------------------         MATCHES       --------------------------------*/
      else if (AnalyzedLine.subcommand[0] == 'M') {
	fprintf(prfStreamOut,"%s\n",currentLineCopy);
      }
      /* -----------------------       UNKNOWN          -------------------------------*/
      else {
	fprintf(stderr,"The following unknown line has been avoided\n->%s\n",currentLineCopy);
      }
      
    }
  }
  
  END_OF_PROFILE:
  
  /*
   * Write down the comments now 
   */
  if (CommentStorageUsed>0) fputs(CommentStorage, prfStreamOut);
  
  /*
   * Write down the additional comments now 
   */
  if (AdditionalComments) fputs(AdditionalComments, prfStreamOut);
  
  /*
   * Terminate profile
   */
  fputs("//\n", prfStreamOut);
  
  fclose(prfStreamIn);
  free(CommentStorage);
  return 0;
}
//---------------------------------------------------------------

void FreeProfile(struct Profile * const prf, const _Bool IsPointer)
{
  struct Profile * nextPrf = prf->next;
  while ( nextPrf != NULL ) {
    FreeScores(&(nextPrf->Scores));
    free(nextPrf->Sequence);
    FreeAverage(&nextPrf->Average);
    struct Profile * tmpPrf = nextPrf->next;
    _mm_free(nextPrf);
    nextPrf = tmpPrf;
  }
  
  FreeScores(&(prf->Scores));
  free(prf->Sequence);
  FreeAverage(&prf->Average);
  if (prf->Pattern) free(prf->Pattern);
  
  if (IsPointer) {
    _mm_free(prf);
  } else {
    memset(prf, 0, sizeof(struct Profile));
  }
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */

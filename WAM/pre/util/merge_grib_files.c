#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <float.h>
#include "bds.h"
#include "gdstool.h"
#include "gds.h"
#include "gribw.h"
#include "pds4.h"

/*
 * Compile: make -f merge_grib_files.c
 *  
 * Requires: libgribw
 *
 * Oyvind.Breivik@met.no
 */

/* Constants */
#define USAGE "Usage: %s file1.grb file2.grb out.grb par\n"
#define VERSION "Version: %s v0.1, 2011-04-07, Oyvind.Breivik@met.no\n"
#define HELP "\n\
A GRIB hack to merge two files\n\
 file2.grb contains one field which is weaved in after fields with parameter par\n\
 in file1.grb. The time is adjusted to match the forecast and analysis time of\n\
 the previous field in file1.grb\n\
"

int main(int argc, char **argv) {

   /* Vars */
   unsigned char *pds, *gds, *bms, *bds;
   unsigned char *pds2, *gds2, *bms2, *bds2;
   unsigned char pdstmp[100];

   int nrec=0;     /* record counter */

   long int len_grib, len_grib2, pos=0, pos2=0; /* GRIB file position counters */

   int i;   // index
   int par; // parameter number

   FILE *input1, *input2, *output;

   /* Initialize */

   /* Help line */
   if (argc == 1) {
      fprintf(stderr, USAGE, argv[0]);
      fprintf(stderr, VERSION, argv[0]);
      fprintf(stderr, HELP);
      exit(8);
   }

   /* Open files */
   if (argc < 5) {
      fprintf(stderr, USAGE, argv[0]);
      exit(8);
   }
   if ((input1 = fopen(argv[1],"rb")) == NULL) {
      fprintf(stderr,"Could not open file: %s\n", argv[1]); 
      exit(7);                                              
   }                                                         
   if ((input2 = fopen(argv[2],"rb")) == NULL) {             
      fprintf(stderr,"Could not open file: %s\n", argv[2]); 
      exit(7);                                              
   }                                                         
   if ((output = fopen(argv[3],"wb")) == NULL) {             
      fprintf(stderr,"Could not open file: %s\n", argv[3]); 
      exit(7);                                              
   }                                                         
   par = atoi(argv[4]);   /* parameter number */

  // Read file2 grid

   len_grib2 = rd_grib_msg(input2, pos2, &pds2, &gds2, &bms2, &bds2);

  // Loop over file1 input records

   while ((len_grib = rd_grib_msg(input1, pos, &pds, &gds, &bms, &bds)) > 0) {

      nrec++; /* count records */

      /* Write GRIB record from file 1 */
      wrt_grib_msg(output, pds, gds, bms, bds);

      // Write modified GRIB record from file 2 if after par field in file1 
      if (PDS_PARAM(pds) == par) {
        // Save date code from last pds field of file1 */
         for (i=12; i<=20; i++) {
            pdstmp[i] = pds[i];
         }
        // Must re-read every time ... crap
         len_grib2 = rd_grib_msg(input2, pos2, &pds2, &gds2, &bms2, &bds2);
        // Copy date code from last field of file1 to pds of file2 */
         for (i=12; i<=20; i++) {
            pds2[i] = pdstmp[i];
         }
         wrt_grib_msg(output, pds2, gds2, bms2, bds2);
         nrec++;
      }

      pos += len_grib; /* file position */
   } /* end for nrec */

   /* Close files */
   fclose(input1);
   fclose(input2);
   fclose(output);

   /* Output */
   if (nrec > 1) {
      printf("Wrote %d records\n", nrec);
   }
   else if (nrec == 1) {
      printf("Wrote one record\n");
   }
   else {
      fprintf(stderr,"ERROR: no records written\n");
      exit(4);
   } /* end if */

   return 0;
} /* end main */

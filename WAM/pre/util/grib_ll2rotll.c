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

#define new_RotLL_GDS(pds,nlon,nlat,lon1,lat1,lon2,lat2,dlon,dlat,) \
	GDStool(pds,NULL,g_init,32,0,g_2bytes,6,nlon,g_2bytes,8,nlat, \
	g_s3bytes,10,(int) (1000.0*(lat1)),g_s3bytes,13,(int) (1000.0*(lon1)), \
	g_s3bytes,17,(int) (1000.0*(lat2)),g_s3bytes,20,(int) (1000.0*(lon2)), \
	g_s2bytes,23,(int) (1000.0*(dlon)+0.5),g_s2bytes,25,(int) (1000.0*(dlat)+0.5), \
	g_s3bytes,32,(int) (1000.0*(lasp), g_s2bytes,35,(int) (1000.0*(losp)),38,(int) flt2ibm(rotang), \
	g_or, 27, ((lat1) < (lat2)) * 64, \
	g_byte,16,128,g_byte,4,255,g_end)
/*

/*
 *
 * Handles rotated lat-lon grids.
 * 
 * Compile: make -f grib_switch_scan.make
 *  
 * Usage: type bilingrb without arguments
 * 
 * Requires: libgribw
 *
 * Oyvind.Breivik@met.no
 */

/* Constants */
#define USAGE "Usage: %s in.grb out.grb scan\n"
#define VERSION "Version: %s v0.1, 2010-12-03, Oyvind.Breivik@met.no\n"
#define HELP "\n\
grib_switch_scan -- A GRIB hack to flip the internal arrangement of arrays\n\
 Handles grids in rotated lat-lon and grids periodic in longitude.\n\
\n\
 scan: the scan direction - 64 is south to north, 0 is north to south\n"
#define MDEG 1000.0
#define N 100000

int main(int argc, char **argv) {

   /* Function prototype */
   int gribindex(int, int, int, int, int);

   /* Vars */
   unsigned char *pds, *gds, *bms, *bds, *gdsnew;
   float *arr, *brr; /* tmp work arrays */

   int nrec=0;     /* record counter */
   int nxny;       /* grid dimension */

   long int len_grib, pos=0; /* GRIB file position counters */

   int nx, ny;           /* grid dimensions, old and new */
   int i, k, l, m, n;    /* indices, 1D and 2D */
   int idir, jdir, iadj; /* logical vars */
   int idir2, jdir2;     /* logical vars */

   float dlon, dlat;                           /* from old GDS */
   float wlon, elon, slat, nlat;               /* old grid boundaries */
   float firstlon, lastlon, firstlat, lastlat; /* for new GDS */

   unsigned char mode; /* resolution and mode flags, GDS 17 */
   int scan, newscan;           /* GRIB record scan direction, GDS 28 */

   FILE *input, *output;

   /* Initialize */

   /* Help line */
   if (argc == 1) {
      fprintf(stderr, USAGE, argv[0]);
      fprintf(stderr, VERSION, argv[0]);
      fprintf(stderr, HELP);
      exit(8);
   }

   /* Open files */
   if (argc < 3) {
      fprintf(stderr, USAGE, argv[0]);
      exit(8);
   }
   if ((input = fopen(argv[1],"rb")) == NULL) {
      fprintf(stderr,"Could not open file: %s\n", argv[1]); 
      exit(7);                                              
   }                                                         
   if ((output = fopen(argv[2],"wb")) == NULL) {             
      fprintf(stderr,"Could not open file: %s\n", argv[2]); 
      exit(7);                                              
   }                                                         
   /* New grid scan direction */
   newscan = atoi(argv[3]);   /* new scan mode */
   idir2 = ((newscan&128)==128);
   idir2 = 1-idir2;          /* West to east is "true" */
   jdir2 = ((newscan&64)==64); /* South to north is "true" */

   /* Preallocate BIG array for speed */
   arr = (float *) malloc(N * sizeof(float));
   
   /* Loop over input records */
   while ((len_grib = rd_grib_msg(input, pos, &pds, &gds, &bms, &bds)) > 0) {

      nrec++; /* count records */

      /* Reallocate arrays if necessary */
      nxny = get_nxny(pds, gds, bms, bds);
      if (nxny > N) {
         free(arr);
         arr = (float *) malloc((nxny) * sizeof(float));
         if (arr == NULL) {
            fprintf(stderr,"Memory allocation failed\n");
         }
         brr = (float *) malloc((nxny) * sizeof(float));
         if (brr == NULL) {
            fprintf(stderr,"Memory allocation failed\n");
         }
      } /* end if nxny */

      /* Unpack data */
      unpk_bds(arr, pds, gds, bms, bds, nxny);

      /* Is grid geographic or rotated geographic? If not, bitch and quit */
      if (!GDS_LatLon(gds) && !GDS_RotLL(gds)) {
         fprintf(stderr,\
         "Error: GRIB rec no %4d not in geographic or rotated geographic co-ords\n",\
          nrec);
         exit(6);
      } 

      /* Old grid parameters */
      scan = GDS_LatLon_scan(gds);
      nx = GDS_LatLon_nx(gds);
      ny = GDS_LatLon_ny(gds);
      dlon = (float) GDS_LatLon_dx(gds)/MDEG;
      dlat = (float) GDS_LatLon_dy(gds)/MDEG;

      /* Grid scan direction */
      idir = ((scan&128)==128);
      idir = 1-idir;          /* West to east is "true" */
      jdir = ((scan&64)==64); /* South to north is "true" */
      iadj = ((scan&32)==32);
      iadj = 1-iadj;          /* Zonal points adjacent is "true" */

      if (idir) {
         wlon = (float) GDS_LatLon_Lo1(gds)/MDEG; /* grid running W to E (default)? */
         elon = (float) GDS_LatLon_Lo2(gds)/MDEG;
      }
      else {
         elon = (float) GDS_LatLon_Lo1(gds)/MDEG; /* or E to W? */
         wlon = (float) GDS_LatLon_Lo2(gds)/MDEG;
      }
         
      if (jdir) {
         slat = (float) GDS_LatLon_La1(gds)/MDEG; /* grid running S to N? */
         nlat = (float) GDS_LatLon_La2(gds)/MDEG;
      }
      else {
         nlat = (float) GDS_LatLon_La1(gds)/MDEG; /* or N to S (default)? */
         slat = (float) GDS_LatLon_La2(gds)/MDEG;
      }

      if (nlat < slat) {
         fprintf(stderr,"Error: scan direction mismatches with south and "
         "north grid boundaries, scan=%d\n",scan);
         exit(5);
      }

      /* Loop over new grid */
      for (m=1; m<=nx; m++) {    /* lons */
         for (n=1; n<=ny; n++) { /* lats */

            /* 1D indices in old and new grid arrays */
            k = gribindex(m,n,nx,ny,scan);
            l = gribindex(m,n,nx,ny,newscan);

            /* Rearrange new array according to new scan direction */
            brr[l] = arr[k];

         } /* end for n */
      } /* end for m */

      /* Compute new grid boundaries in (rotated) lat-lon */
      if (idir2) {
         firstlon = wlon;  /* grid starts on western rim (default) */
         lastlon = elon;
      }
      else {
         firstlon = elon;
         lastlon = wlon;
      }
         
      if (jdir2) {
         firstlat = slat;
         lastlat = nlat;
      }
      else {
         firstlat = nlat;  /* grid starts on northern rim (default) */
         lastlat = slat;
      }

      /* Remake grid definition section, GDS */
      mode = gds[16];
      if ((gdsnew = malloc(GDS_LEN(gds))) == NULL) {
         fprintf(stderr,"malloc failure gdsnew\n");
         exit(8);
      }

      /* Copy old GDS ... */
      for (i=0; i<GDS_LEN(gds); i++) gdsnew[i] = gds[i];

      /* ... change some (works for rotated lat-lon too) ... */
      gds = new_LatLon_GDS(pds,nx,ny,firstlon,firstlat,lastlon,
                           lastlat,dlon,dlat);
      /* ... and copy again */
      for (i=6; i<GDS_LEN(gds); i++) gdsnew[i] = gds[i];
      gdsnew[16] = mode; /* Restore the mode flags */

      /* Remake bitmap section, BMS */
      bms = mk_BMS(pds, brr, &nxny, UNDEFINED_LOW, UNDEFINED_HIGH);

      /* Remake binary data section, BDS */
      if (BDS_BinScale(bds)) {               /* ECMWF style? */
         set_BDSMinBits(BDS_NumBits(bds));
         set_BDSMaxBits(BDS_NumBits(bds));
      }
      bds = mk_BDS(pds, brr, nxny);

      /* Write new GRIB record */
      wrt_grib_msg(output, pds, gdsnew, bms, bds);

      /* Clean up */
      if (gds != NULL) {
        free(gds);
        gds = NULL;
      }
      if (gdsnew != NULL) {
        free(gdsnew);
        gdsnew = NULL;
      }
      if (bms != NULL) {
        free(bms);
        bms = NULL;
      }
      if (bds != NULL) {
        free(bds);
        bds = NULL;
      }

      pos += len_grib; /* file position */
   } /* end for nrec */

   /* Close files */
   fclose(input);
   fclose(output);

   /* Clean up */
   if (arr != NULL) {
     free(arr);
     arr = NULL;
   }
   if (brr != NULL) {
     free(brr);
     brr = NULL;
   }

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


int gribindex(int i, int j, int nx, int ny, int scan)
{
   /* Compute 1D GRIB record index from two independent grid indices (east and
    * north indices) */

   int k;
   int ii, jj, idir, jdir, iadj;


   /* Grid orientation, GDS 28 */
   idir = ((scan&128)==128);
   idir = 1-idir;          /* west to east */
   jdir = ((scan&64)==64); /* south to north */
   iadj = ((scan&32)==32);
   iadj = 1-iadj;          /* zonal points adjacent */

   /* Positive direction index i? */
   ii = i;
   if (!idir) {
      ii = nx-i+1;     /* not default */
   }

   /* Positive direction index j? */
   jj = j;
   if (!jdir) {
      jj = ny-j+1;     /* default */
   }

   /* Index i adjacent (grid running east-west)? */
   if (iadj) {
      jj = (jj-1)*nx;  /* default */
   }
   else {
      ii = (ii-1)*ny;  /* index j adjacent, running north-south */
   }

   k = ii+jj-1;

   return k;
} /* end function gribindex */

#include "test.h"
#include <fitsio.h>
#define STR_LEN_0	127

_Bool doTest_wf = TRUE;

/*....................................................................*/
void
_writeWCS(fitsfile *fptr, const int i, int axesOrder[4], float cdelt[4]\
  , double crpix[4], double crval[4], char ctype[4][9], char cunit[4][9]){

  char myStr[9];
  int status = 0;

  sprintf(myStr, "CTYPE%d  ", i+1);
  fits_write_key(fptr, TSTRING, myStr, &ctype[axesOrder[i]], "", &status);
  sprintf(myStr, "CDELT%d  ", i+1);
  fits_write_key(fptr, TFLOAT, myStr, &cdelt[axesOrder[i]], "", &status);
  sprintf(myStr, "CRPIX%d  ", i+1);
  fits_write_key(fptr, TDOUBLE, myStr, &crpix[axesOrder[i]], "", &status);
  sprintf(myStr, "CRVAL%d  ", i+1);
  fits_write_key(fptr, TDOUBLE, myStr, &crval[axesOrder[i]], "", &status);
  sprintf(myStr, "CUNIT%d  ", i+1);
  fits_write_key(fptr, TSTRING, myStr, &cunit[axesOrder[i]], "", &status);	
}

/*....................................................................*/
unsigned long flattenImageIndices(const int sizes[3], const int xi, const int yi, const int zi){
  unsigned long ppi;
  ppi = ((unsigned long)zi*sizes[1] + yi)*sizes[0] + xi;
  return ppi;
}

/*....................................................................*/
void 
write3Dfits(char *fileName, double *cube, const int sizes[3], double xLos[3], double xHis[3]){
  /*
****.
  */
  const int numAxes=3;
  double bscale,bzero;
  int axesOrder[] = {0,1,2};
  char ctype[numAxes][9],cunit[numAxes][9];
  double crpix[numAxes],crval[numAxes];
  float cdelt[numAxes];
  float *row;
  int xi,yi,zi,i;
  fitsfile *fptr;
  int status = 0;
  int naxis=numAxes, bitpix=-32;
  long naxes[numAxes];
  long int fpixels[numAxes],lpixels[numAxes];
  char negfile[100]="! ";
  unsigned long ppi;

const int testXi=0,testYi=1,testZi=0;

  row = malloc(sizeof(*row)*sizes[0]);

  for(i=0;i<numAxes;i++)
    naxes[axesOrder[i]] = sizes[i];

  fits_create_file(&fptr, fileName, &status);

  if(status!=0){
    printf("Overwriting existing fits file.\n");
    status=0;
    strcat(negfile, fileName);
    fits_create_file(&fptr, negfile, &status);
  }

  /* Write FITS header */ 
  fits_create_img(fptr, bitpix, naxis, naxes, &status);

  i = 0;
  sprintf(ctype[axesOrder[i]], "X      ");
  cdelt[axesOrder[i]] = (float)(xHis[i] - xLos[i])/(float)(sizes[i] - 1);
  crpix[axesOrder[i]] = 0.5;
  crval[axesOrder[i]] = xLos[i];
  sprintf(cunit[axesOrder[i]], "M      ");

  i = 1;
  sprintf(ctype[axesOrder[i]], "Y      ");
  cdelt[axesOrder[i]] = (float)(xHis[i] - xLos[i])/(float)(sizes[i] - 1);
  crpix[axesOrder[i]] = 0.5;
  crval[axesOrder[i]] = xLos[i];
  sprintf(cunit[axesOrder[i]], "M      ");

  i = 2;
  sprintf(ctype[axesOrder[i]], "Z      ");
  cdelt[axesOrder[i]] = (float)(xHis[i] - xLos[i])/(float)(sizes[i] - 1);
  crpix[axesOrder[i]] = 0.5;
  crval[axesOrder[i]] = xLos[i];
  sprintf(cunit[axesOrder[i]], "M      ");

  bscale  =1.0e0;
  bzero   =0.0e0;

  for(i=0;i<numAxes;i++)
    _writeWCS(fptr, i, axesOrder, cdelt, crpix, crval, ctype, cunit);

  fits_write_key(fptr, TDOUBLE, "BSCALE  ", &bscale,        "", &status);
  fits_write_key(fptr, TDOUBLE, "BZERO   ", &bzero,         "", &status);

  /* Write FITS data
   */
  for(zi=0;zi<sizes[2];zi++){
    for(yi=0;yi<sizes[1];yi++){
      for(xi=0;xi<sizes[0];xi++){
        ppi = flattenImageIndices(sizes, xi, yi, zi);
        row[xi] = (float)cube[ppi];
if(doTest_wf && xi==testXi && yi==testYi && zi==testZi){
  printf("In write_fits.write3Dfits(). Hypercube flat index at [%d,%d,%d] is %lu\n", testXi, testYi, testZi, ppi);
}
      }
      fpixels[axesOrder[0]] = 1;
      fpixels[axesOrder[1]] = yi+1;
      fpixels[axesOrder[2]] = zi+1;
      lpixels[axesOrder[0]] = sizes[0];
      lpixels[axesOrder[1]] = yi+1;
      lpixels[axesOrder[2]] = zi+1;

if(doTest_wf && zi==testZi && yi==testYi){
  printf("  ...Value at [%d,%d,%d] is %e\n", testXi, testYi, testZi, row[testXi]);
}
      fits_write_subset(fptr, TFLOAT, fpixels, lpixels, row, &status);
    }
  }

  fits_close_file(fptr, &status);

  free(row);
}


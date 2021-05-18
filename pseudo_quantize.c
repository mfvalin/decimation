#include <stdio.h>
#include <unistd.h>
#include <math.h>

void QuantizeRestore(float *z, int ni, int li, int nj, int nbits){
  float min, max, range, scale ;
  float *t ;
  int i, j, quant, mask, limit ;
  int *irange = (int *) &range ;

  max = z[0] ; min = max ;               // find minimum and maximum value
  for(j=0, t=z ; j<nj ; j++,t+=li){
    for(i=0 ; i<ni ; i++){
      min = (min > t[i]) ? t[i] : min ;
      max = (max < t[i]) ? t[i] : max ;
    }
  }
  mask = 1 << nbits ;
  limit = mask - 1 ;
  range = max - min ;
  printf("0 - min = %f, max = %f, range = %f\n",min,max,range) ;
  *irange = *irange - 1 ;           // if we have an exact power of 2, bump down ;
  *irange = *irange & 0x7F800000 ;  // flush mantissa, power of 2 <= range
  range *= 2.0f ;                   // power of 2 >= mantissa
  scale = 1.0f / range ;
  scale *= mask ;
  printf("1 - min = %f, max = %f, range = %f, quantum = %f\n",min,max,range,1.0f/scale) ;
  for(j=0, t=z ; j<nj ; j++,t+=li){
    for(i=0 ; i<ni ; i++){
      quant = (t[i] - min) * scale + .5 ;    // quantize
      quant = (quant > limit) ? limit : quant ;
      t[i] = min + quant / scale ;           // restore
    }
  }
  max = z[0] ; min = max ;
  for(j=0, t=z ; j<nj ; j++,t+=li){
    for(i=0 ; i<ni ; i++){
      min = (min > t[i]) ? t[i] : min ;
      max = (max < t[i]) ? t[i] : max ;
    }
  }
  printf("2 - min = %f, max = %f, range = %f\n",min,max,max-min) ;
  
}

#define NI 31
#define NJ 21
#define LI 37
#define NBITS 15

int main(){
  float z[LI*NJ] ;
  int i, j, ij ;
  ij = 0 ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      z[ij+i] = (i - NI*.5f)*(i - NI*.5f) + (j - NJ*.5f)*(j - NJ*.5f) ;
      z[ij+i] = sqrtf(z[ij+i]) ;
    }
    ij += LI ;
  }
//   z[0] = 32.707107f ;
  QuantizeRestore(z, NI, LI, NJ, NBITS) ;
  printf("\n") ;
  QuantizeRestore(z, NI, LI, NJ, NBITS) ;
  printf("\n") ;
  QuantizeRestore(z, NI, LI, NJ, NBITS) ;
  return 0;
}

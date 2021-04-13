/*
 * Copyright (C) 2021  Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * Author:
 *     M. Valin,   Recherche en Prevision Numerique, April 2021
 */
//
// linear decimation and inverse decimation (restore to initial dimensions) routines
// support is for (preferably large) float (Fortran real *4) data elements
//
// decimation by a factor of 2,3, or 4 / restore 
//
// decimation is performed with simple N x N averaging (N = 2/3/4)
// restoration is performed using linear interpolation
//
// a bi-linear array will be restored without loss (other than float rounding errors)
//
// number of elements after decimation of npts elements by nd
int N_decimated(int npts, int nd){
  int ntup = (npts-2) / nd ;  // number of averaged points
  if(ntup == 1) return npts ; // there MUST be more than 1 group of averaged elements
  return (npts - (ntup * nd) + ntup) ;
}
//
// ============================= decimation by 2 =============================
//
//  +---+     +---+---+---+---+     +---+---+---+
//  | 1 | ... | i |i+1|i+2|i+3| ... |   |   | n |  before decimation : n elements
//  +---+     +---+---+---+---+     +---+---+---+
//
//  +---+     +-------+-------+     +-------+---+
//  | 1 | ... |  i+1  |  i+2  | ... |   j   | d |  after decimation : d = (n-2)/2 + 2 + mod(n-2, 2) elements
//  +---+     +-------+-------+     +-------+---+                     i = 1 + mod(n-2, 2)
//
//  the first and last elements are always kept verbatim, if n-2 is not a multiple of 2, 
//  mod(n-2, 2) elements are kept verbatim after the first element
//
// decimation by 2 of a 1D array
static inline int Decimate_by2_1d(float *src1, float *dst, int ni, int li){
  int i, j, pre ;
  int ntup = (ni-2) / 2 ;    // number of averaged groups of rows
  float *src2 = src1 + li ;  // second row to be averaged (same as src1 if li is 0)

  pre = ni - (2 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = (src1[i] + src2[i]) * .5f;     // first pre point(s) verbatim

  for(i = pre, j = pre ; i < pre + ntup ; i++ , j+=2) {              // averaged tuples (pairs)
    dst[i] = ( src1[j] + src1[j+1] + src2[j] + src2[j+1] ) * .25f ;
  }

  dst[i] = ( src1[j] + src2[j] ) * .5f ;                             // last point verbatim

  return i + 1 ;
}

// decimation by 2 of a 2D array
int Decimate_by2_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/2 ;                   // number of averaged pairs along i
  int jtup = (nj-2)/2 ;                   // number of averaged pairs along j
  int ipre = ni - (2 * itup) - 1 ;
  int jpre = nj - (2 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;             // total number of decimated points along i
  int j;

  if(itup <= 1 || jtup <= 1) return -1 ;  // one or both dimensions too small

  for(j = 0 ; j < jpre ; j ++){
    Decimate_by2_1d(src, dst, ni, 0) ;    // first jpre row(s) (only decimated along i)
    src += li ;
    dst += nid ;
  }
  for(j=1 ; j <= jtup ; j++) {            // jtup  blocks of 2 rows
    Decimate_by2_1d(src, dst, ni, li) ;
    src += (li+li) ;
    dst += nid ;
  }
  Decimate_by2_1d(src, dst, ni, 0) ;      // last row (only decimated along i)
  return 0 ;
}

// inverse decimation by 2 of a 1D array (restore to original dimension)
static inline void UnDecimate_by2_1d(float *src, float *dst, int ni){
  int i, j, pre ;
  int ntup = (ni-2) / 2 ;    // number of averaged groups of rows
  float third  = 1.0f / 3.0f ;
  float third2 = 2.0f / 3.0f ;

  pre = ni - (2 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = src[i] ;                 // first non averaged points

  dst[pre  ] = src[pre  ] * third2 +  src[pre-1] * third ;     // first group of averaged points
  dst[pre+1] = src[pre  ] * .75f  +  src[pre+1] * .25f ;

  for(i = pre+1, j = pre+2 ; i < pre + ntup - 1 ; i++ , j+=2) {
    dst[j  ] = src[i-1] * .25f + src[i] * .75f ;               // (ntup -2) middle groups
    dst[j+1] = src[i+1] * .25f + src[i] * .75f ;
  }

  dst[j++] = src[i-1] * .25f + src[i] * .75f ;                 // last group of averaged points
  dst[j++] = src[i] * third2 + src[i+1] * third ;

  dst[j  ] = src[i+1] ;                                        // last point
}

// inverse decimation by 2 of a 2D array (restore to original dimensions)
int UnDecimate_by2_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/2 ;                  // number of averaged pairs along i
  int jtup = (nj-2)/2 ;                  // number of averaged pairs along j
  int ipre = ni - (2 * itup) - 1 ;
  int jpre = nj - (2 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;            // total number of decimated points along i
  int i, j;
  float third  = 1.0f / 3.0f ;
  float third2 = 2.0f / 3.0f ;
  float *topd = dst + li * (nj -1 ) ;    // top row of undecimated array used as temporary 

  for(j = 0 ; j < jpre ; j ++){
    UnDecimate_by2_1d(src, dst, ni) ;    // undecimate first pre rows (no processing along j)
    src += nid ;
    dst += li ;
  }
  // undecimate first group of rows
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * third2 + src[i-nid] * third ;  // along j
  UnDecimate_by2_1d(topd, dst, ni) ;                                         // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * .75f   + src[i+nid] * .25f ;   // along j
  UnDecimate_by2_1d(topd, dst + li, ni) ;                                    // along i
  dst += (li+li) ;
  src += nid ;
  for(j=2 ; j < jtup ; j++){ // undecimate second group -> next to last  group (jtup -2 groups)
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * .75f + src[i-nid] * .25f ;   // along j
    UnDecimate_by2_1d(topd, dst, ni) ;                                       // along i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * .75f + src[i+nid] * .25f ;   // along j
    UnDecimate_by2_1d(topd, dst + li, ni) ;                                  // along i
    dst += (li+li) ;
    src += nid ;
  }
  // undecimate last group
  dst[0] = 0.1f ; dst[li] = 0.2f ;
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * .75f   + src[i-nid] * .25f ;   // along j
  UnDecimate_by2_1d(topd, dst, ni) ;                                         // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * third2 + src[i+nid] * third ;  // along j
  UnDecimate_by2_1d(topd, dst + li, ni) ;                                    // along i
  dst += (li+li) ;
  src += nid ;
  // undecimate last row (no processing along j)
  UnDecimate_by2_1d(src, dst, ni) ;
  return 0 ;
}
//
// ============================= decimation by 3 =============================
//
//  +---+     +---+---+---+     +---+---+---+---+
//  | 1 | ... | i |i+1|i+2| ... |   |   |   | n |  before decimation : n elements
//  +---+     +---+---+---+     +---+---+---+---+
//
//  +---+     +-----------+     +-----------+---+
//  | 1 | ... |     i     | ... |     j     | d |  after decimation : d = (n-2)/3 + 2 + mod(n-2, 3) elements
//  +---+     +-----------+     +-----------+---+                     i = 1 + mod(n-2, 3)
//
//  the first and last elements are always kept verbatim, if n-2 is not a multiple of 3, 
//  mod(n-2, 3) elements are kept verbatim after the first element
//
// decimation by 3 along i of one row
static inline int Decimate_by3_1d(float *src1, float *dst, int ni, int li){
  int i, j, pre ;
  int ntup = (ni-2) / 3 ;             // number of averaged groups of rows
  float *src2 = src1 + li ;
  float *src3 = src2 + li ;
  float third = 1.0f / 3.0f ;
  float ninth = 1.0f / 9.0f ;

  pre = ni - (3 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = (src1[i] + src2[i] + src3[i]) * third;     // first pre point(s) verbatim

  for(i = pre, j = pre ; i < pre + ntup ; i++ , j+=3) {                          // averaged tuples (trios)
    dst[i] = ( src1[j] + src1[j+1] + src1[j+2] + 
               src2[j] + src2[j+1] + src2[j+2] + 
               src3[j] + src3[j+1] + src3[j+2]    ) * ninth ;
  }

  dst[i] = ( src1[j] + src2[j] + src3[j] ) * third ;                             // last point verbatim

  return i + 1 ;
}

// decimation by 3 of a 2D array
int Decimate_by3_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/3 ;                   // number of averaged trios along i
  int jtup = (nj-2)/3 ;                   // number of averaged trios along j
  int ipre = ni - (3 * itup) - 1 ;
  int jpre = nj - (3 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;             // total number of decimated points along i
  int j;

  if(itup <= 1 || jtup <= 1) return -1 ;  // one or both dimensions too small

  for(j = 0 ; j < jpre ; j ++){
    Decimate_by3_1d(src, dst, ni, 0) ;    // first row
    src += li ;
    dst += nid ;
  }
  for(j=1 ; j <= jtup ; j++) {            // jtup  blocks of 3 rows
    Decimate_by3_1d(src, dst, ni, li) ;
    src += (li+li+li) ;
    dst += nid ;
  }
  Decimate_by3_1d(src, dst, ni, 0) ;      // last row
  return 0 ;
}

// inverse decimation by 3 along i of a 1D array (restore to original dimension)
static inline void UnDecimate_by3_1d(float *src, float *dst, int ni){
  int i, j, pre ;
  int ntup = (ni-2) / 3 ;             // number of averaged groups of rows
  float third  = 1.0f / 3.0f ;
  float third2 = 2.0f / 3.0f ;

  pre = ni - (3 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = src[i] ;                   // first non averaged points

  dst[pre  ] = src[pre  ] * .5f +  src[pre-1] * .5f ;            // first group of averaged points
  dst[pre+1] = src[pre  ] ;
  dst[pre+2] = src[pre  ] * third2  +  src[pre+1] * third ;
  for(i = pre+1, j = pre+3 ; i < pre + ntup - 1 ; i++ , j+=3) {
    dst[j  ] = src[i-1] * third + src[i] * third2;               // ntup middle groups
    dst[j+1] = src[i] ;
    dst[j+2] = src[i+1] * third + src[i] * third2 ;
  }
  dst[j++] = src[i-1] * third + src[i] * third2 ;                // last group of averaged points
  dst[j++] = src[i] ;
  dst[j++] = src[i] * .5f + src[i+1] * .5f ;

  dst[j  ] = src[i+1] ;                                          // last point
}

// inverse decimation by 3 of a 2D array (restore to original dimensions)
int UnDecimate_by3_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/3 ;                  // number of averaged trios along i
  int jtup = (nj-2)/3 ;                  // number of averaged trios along j
  int ipre = ni - (3 * itup) - 1 ;
  int jpre = nj - (3 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;            // total number of decimated points along i
  int i, j;
  float third  = 1.0f / 3.0f ;
  float third2 = 2.0f / 3.0f ;
  float *topd = dst + li * (nj -1 ) ;    // top row of undecimated array used as temporary

  for(j = 0 ; j < jpre ; j ++){
    UnDecimate_by3_1d(src, dst, ni) ;   // undecimate first pre rows (no processing along j)
    src += nid ;
    dst += li ;
  }
  // undecimate first group of rows
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * .5f + src[i-nid] * .5f ;          // along j
  UnDecimate_by3_1d(topd, dst, ni) ;                                            // along i
  UnDecimate_by3_1d(src, dst + li, ni) ;                                        // along j and i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * third2   + src[i+nid] * third ;   // along j
  UnDecimate_by3_1d(topd, dst + li + li, ni) ;                                  // along i
  dst += (li+li+li) ;
  src += nid ;
  // undecimate second group -> next to last  group (jtup -2 groups)
  for(j=2 ; j < jtup ; j++){
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * third2 + src[i-nid] * third ;   // along j
    UnDecimate_by3_1d(topd, dst, ni) ;                                          // along i
    UnDecimate_by3_1d(src, dst + li, ni) ;                                      // along j and i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * third2 + src[i+nid] * third ;   // along j
    UnDecimate_by3_1d(topd, dst + li + li, ni) ;                                // along i
    dst += (li+li+li) ;
    src += nid ;
  }
  // undecimate last group
  dst[0] = 0.1f ; dst[li] = 0.2f ;
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * third2   + src[i-nid] * third ;   // along j
  UnDecimate_by3_1d(topd, dst, ni) ;                                            // along i
  UnDecimate_by3_1d(src, dst + li, ni) ;                                        // along j and i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * .5f + src[i+nid] * .5f ;          // along j
  UnDecimate_by3_1d(topd, dst + li + li, ni) ;                                  // along i
  dst += (li+li+li) ;
  src += nid ;
  // undecimate last row (no processing along j)
  UnDecimate_by3_1d(src, dst, ni) ;
  return 0 ;
}
//
// ============================= decimation by 4 =============================
//
//  +---+     +---+---+---+---+     +---+---+---+---+---+
//  | 1 | ... | i |i+1|i+2|i+3| ... |   |   |   |   | n |  before decimation : n elements
//  +---+     +---+---+---+---+     +---+---+---+---+---+
//
//  +---+     +---------------+     +---------------+---+
//  | 1 | ... |       i       | ... |       j       | d |  after decimation : d = (n-2)/4 + 2 + mod(n-2, 4) elements
//  +---+     +---------------+     +---------------+---+                     i = 1 + mod(n-2, 4)
//
//  the first and last elements are always kept verbatim, if n-2 is not a multiple of 4, 
//  mod(n-2, 4) elements are kept verbatim after the first element
//
// decimation by 4 of a 1D array
static inline int Decimate_by4_1d(float *src1, float *dst, int ni, int li){
  int i, j, pre ;
  int ntup = (ni-2) / 4 ;    // number of averaged groups of rows (quads)
  float *src2 = src1 + li ;  // second row to be averaged (same as src1 if li is 0)
  float *src3 = src2 + li ;  // third row to be averaged (same as src1 if li is 0)
  float *src4 = src3 + li ;  // fourth row to be averaged (same as src1 if li is 0)
  float sixteenth = 1.0f/16.0f ;

  pre = ni - (4 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = (src1[i] + src2[i] + src3[i] + src4[i]) * .25f;     // first pre point(s) verbatim

  for(i = pre, j = pre ; i < pre + ntup ; i++ , j+=4) {                                   // averaged tuples (quads)
    dst[i] = ( src1[j] + src1[j+1] + src1[j+2] + src1[j+3] + 
               src2[j] + src2[j+1] + src2[j+2] + src2[j+3] +
               src3[j] + src3[j+1] + src3[j+2] + src3[j+3] +
               src4[j] + src4[j+1] + src4[j+2] + src4[j+3]   ) * sixteenth ;
  }

  dst[i] = ( src1[j] + src2[j] + src3[j] + src4[j] ) * .25f ;                             // last point verbatim

  return i + 1 ;
}

// decimation by 4 of a 2D array
int Decimate_by4_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/4 ;                   // number of averaged pairs along i
  int jtup = (nj-2)/4 ;                   // number of averaged pairs along j
  int ipre = ni - (4 * itup) - 1 ;
  int jpre = nj - (4 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;             // total number of decimated points along i
  int j;

  if(itup <= 1 || jtup <= 1) return -1 ;  // one or both dimensions too small

  for(j = 0 ; j < jpre ; j ++){
    Decimate_by4_1d(src, dst, ni, 0) ;    // first row
    src += li ;
    dst += nid ;
  }
  for(j=1 ; j <= jtup ; j++) {            // jtup  blocks of 4 rows
    Decimate_by4_1d(src, dst, ni, li) ;
    src += (li+li+li+li) ;
    dst += nid ;
  }
  Decimate_by4_1d(src, dst, ni, 0) ;      // last row
  return 0 ;
}

// inverse decimation by 4 of a 1D array (restore to original dimension)
static inline void UnDecimate_by4_1d(float *src, float *dst, int ni){
  int i, j, pre ;
  int ntup = (ni-2) / 4 ;         // number of averaged groups of rows (quads)
  float eighth  = 1.0f / 8.0f ;
  float eighth3 = 3.0f / 8.0f ;
  float eighth5 = 5.0f / 8.0f ;
  float eighth7 = 7.0f / 8.0f ;
  float fifth   = 1.0f / 5.0f ;
  float fifth2  = 2.0f / 5.0f ;
  float fifth3  = 3.0f / 5.0f ;
  float fifth4  = 4.0f / 5.0f ;

  pre = ni - (4 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = src[i] ;                   // first non averaged points

  dst[pre  ] = src[pre] * fifth2  +  src[pre-1] * fifth3 ;       // first group of averaged points
  dst[pre+1] = src[pre] * fifth4  +  src[pre-1] * fifth ;
  dst[pre+2] = src[pre] * eighth7 +  src[pre+1] * eighth ;
  dst[pre+3] = src[pre] * eighth5 +  src[pre+1] * eighth3 ;

  for(i = pre+1, j = pre+4 ; i < pre + ntup - 1 ; i++ , j+=4) {
    dst[j  ] = src[i-1] * eighth3 + src[i] * eighth5;            // (ntup -2) middle groups
    dst[j+1] = src[i-1] * eighth  + src[i] * eighth7 ;
    dst[j+2] = src[i+1] * eighth  + src[i] * eighth7 ;
    dst[j+3] = src[i+1] * eighth3 + src[i] * eighth5 ;
  }

  dst[j++] = src[i-1] * eighth3 + src[i]   * eighth5 ;           // last group of averaged points
  dst[j++] = src[i-1] * eighth  + src[i]   * eighth7 ;
  dst[j++] = src[i]   * fifth4  + src[i+1] * fifth ;
  dst[j++] = src[i]   * fifth2  + src[i+1] * fifth3 ;

  dst[j  ] = src[i+1] ;                                          // last point
}

// inverse decimation by 4 of a 2D array (restore to original dimensions)
int UnDecimate_by4_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/4 ;                  // number of averaged pairs along i
  int jtup = (nj-2)/4 ;                  // number of averaged pairs along j
  int ipre = ni - (4 * itup) - 1 ;
  int jpre = nj - (4 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;            // total number of decimated points along i
  int i, j;
  float eighth  = 1.0f / 8.0f ;
  float eighth3 = 3.0f / 8.0f ;
  float eighth5 = 5.0f / 8.0f ;
  float eighth7 = 7.0f / 8.0f ;
  float fifth   = 1.0f / 5.0f ;
  float fifth2  = 2.0f / 5.0f ;
  float fifth3  = 3.0f / 5.0f ;
  float fifth4  = 4.0f / 5.0f ;
  float *topd = dst + li * (nj -1 ) ;    // top row of undecimated array used as temporary

  for(j = 0 ; j < jpre ; j ++){
    UnDecimate_by4_1d(src, dst, ni) ;   // undecimate first pre rows (no processing along j)
    src += nid ;
    dst += li ;
  }

  // undecimate first group of rows
  dst[0] = 0.1f ; dst[li] = 0.2f ;
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth2 + src[i-nid] * fifth3 ;      // along j
  UnDecimate_by4_1d(topd, dst, ni) ;                                              // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth4 + src[i-nid] * fifth  ;      // along j
  UnDecimate_by4_1d(topd, dst + li, ni) ;                                         // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth7   + src[i+nid] * eighth  ;  // along j
  UnDecimate_by4_1d(topd, dst + li + li, ni) ;                                    // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth5   + src[i+nid] * eighth3 ;  // along j
  UnDecimate_by4_1d(topd, dst + li + li + li, ni) ;                               // along i
  dst += (li+li+li+li) ;
  src += nid ;

  // undecimate second group -> next to last  group (jtup -2 groups)
  for(j=2 ; j < jtup ; j++){
  dst[0] = 0.3f ; dst[li] = 0.4f ;
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth5 + src[i-nid] * eighth3 ;  // along j
    UnDecimate_by4_1d(topd, dst, ni) ;                                            // along i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth7 + src[i-nid] * eighth ;   // along j
    UnDecimate_by4_1d(topd, dst + li, ni) ;                                       // along i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth7 + src[i+nid] * eighth ;   // along j
    UnDecimate_by4_1d(topd, dst + li + li, ni) ;                                  // along i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth5 + src[i+nid] * eighth3 ;  // along j
    UnDecimate_by4_1d(topd, dst + li + li + li, ni) ;                             // along i
    dst += (li+li+li+li) ;
    src += nid ;
  }

  // undecimate last group
  dst[0] = 0.1f ; dst[li] = 0.2f ;
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth5   + src[i-nid] * eighth3 ;  // along j
  UnDecimate_by4_1d(topd, dst, ni) ;                                              // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * eighth7   + src[i-nid] * eighth  ;  // along j
  UnDecimate_by4_1d(topd, dst + li, ni) ;                                         // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth4 + src[i+nid] * fifth ;       // along j
  UnDecimate_by4_1d(topd, dst + li + li, ni) ;                                    // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth2 + src[i+nid] * fifth3 ;      // along j
  UnDecimate_by4_1d(topd, dst + li + li + li, ni) ;                               // along i
  dst += (li+li+li+li) ;
  src += nid ;

  // undecimate last row (no processing along j)
  UnDecimate_by4_1d(src, dst, ni) ;
  return 0 ;
}
//
// ============================= decimation by 4 =============================
//
//  +---+     +---+---+---+---+---+     +---+---+---+---+---+---+
//  | 1 | ... | i |i+1|i+2|i+3|i+4| ... |   |   |   |   |   | n |  before decimation : n elements
//  +---+     +---+---+---+---+---+     +---+---+---+---+---+---+
//
//  +---+     +-------------------+     +-------------------+---+
//  | 1 | ... |         i         | ... |         j         | d |  after decimation : d = (n-2)/5 + 2 + mod(n-2, 5) elements
//  +---+     +-------------------+     +-------------------+---+                     i = 1 + mod(n-2, 5)
//
//  the first and last elements are always kept verbatim, if n-2 is not a multiple of 5, 
//  mod(n-2, 5) elements are kept verbatim after the first element
//
// decimation by 5 of a 1D array
static inline int Decimate_by5_1d(float *src1, float *dst, int ni, int li){
  int i, j, pre ;
  int ntup = (ni-2) / 5 ;    // number of averaged groups of rows (quads)
  float *src2 = src1 + li ;  // second row to be averaged (same as src1 if li is 0)
  float *src3 = src2 + li ;  // third row to be averaged (same as src1 if li is 0)
  float *src4 = src3 + li ;  // fourth row to be averaged (same as src1 if li is 0)
  float *src5 = src4 + li ;  // fifth row to be averaged (same as src1 if li is 0)
  float scale = 1.0f / 25.0f ;
  float fifth = 1.0f / 5.0f ;

  pre = ni - (5 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = (src1[i] + src2[i] + src3[i] + src4[i] + src5[i]) * fifth;     // first pre point(s) verbatim

  for(i = pre, j = pre ; i < pre + ntup ; i++ , j+=5) {                                   // averaged tuples (quads)
    dst[i] = ( src1[j] + src1[j+1] + src1[j+2] + src1[j+3] + src1[j+4] + 
               src2[j] + src2[j+1] + src2[j+2] + src2[j+3] + src2[j+4] +
               src3[j] + src3[j+1] + src3[j+2] + src3[j+3] + src3[j+4] +
               src4[j] + src4[j+1] + src4[j+2] + src4[j+3] + src4[j+4] +
               src5[j] + src5[j+1] + src5[j+2] + src5[j+3] + src5[j+4]   ) * scale ;
  }

  dst[i] = ( src1[j] + src2[j] + src3[j] + src4[j] + src5[j] ) * fifth ;                   // last point verbatim

  return i + 1 ;
}

// decimation by 5 of a 2D array
int Decimate_by5_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/5 ;                   // number of averaged pairs along i
  int jtup = (nj-2)/5 ;                   // number of averaged pairs along j
  int ipre = ni - (5 * itup) - 1 ;
  int jpre = nj - (5 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;             // total number of decimated points along i
  int j;

  if(itup <= 1 || jtup <= 1) return -1 ;  // one or both dimensions too small

  for(j = 0 ; j < jpre ; j ++){
    Decimate_by5_1d(src, dst, ni, 0) ;    // first row
    src += li ;
    dst += nid ;
  }
  for(j=1 ; j <= jtup ; j++) {            // jtup  blocks of 5 rows
    Decimate_by5_1d(src, dst, ni, li) ;
    src += (li+li+li+li+li) ;
    dst += nid ;
  }
  Decimate_by5_1d(src, dst, ni, 0) ;      // last row
  return 0 ;
}

// inverse decimation by 5 of a 1D array (restore to original dimension)
static inline void UnDecimate_by5_1d(float *src, float *dst, int ni){
  int i, j, pre ;
  int ntup = (ni-2) / 5 ;         // number of averaged groups of rows (quads)
  float sixth2  = 2.0f / 6.0f ;
  float sixth4  = 4.0f / 6.0f ;
  float fifth   = 1.0f / 5.0f ;
  float fifth2  = 2.0f / 5.0f ;
  float fifth3  = 3.0f / 5.0f ;
  float fifth4  = 4.0f / 5.0f ;

  pre = ni - (5 * ntup) - 1 ;
  for(i = 0 ; i < pre ; i++) dst[i] = src[i] ;                   // first non averaged points

  dst[pre  ] = src[pre] * sixth2  +  src[pre-1] * sixth4 ;       // first group of averaged points
  dst[pre+1] = src[pre] * sixth4  +  src[pre-1] * sixth2 ;
  dst[pre+2] = src[pre] ;
  dst[pre+3] = src[pre] * fifth4  +  src[pre+1] * fifth ;
  dst[pre+4] = src[pre] * fifth3  +  src[pre+1] * fifth2 ;

  for(i = pre+1, j = pre+5 ; i < pre + ntup - 1 ; i++ , j+=5) {
    dst[j  ] = src[i-1] * fifth2  + src[i] * fifth3;            // (ntup -2) middle groups
    dst[j+1] = src[i-1] * fifth   + src[i] * fifth4 ;
    dst[j+2] = src[i] ;
    dst[j+3] = src[i+1] * fifth   + src[i] * fifth4 ;
    dst[j+4] = src[i+1] * fifth2  + src[i] * fifth3 ;
  }

  dst[j++] = src[i-1] * fifth2  + src[i]   * fifth3 ;           // last group of averaged points
  dst[j++] = src[i-1] * fifth   + src[i]   * fifth4 ;
  dst[j++] = src[i] ;
  dst[j++] = src[i]   * sixth4  + src[i+1] * sixth2 ;
  dst[j++] = src[i]   * sixth2  + src[i+1] * sixth4 ;

  dst[j  ] = src[i+1] ;                                          // last point
}

// inverse decimation by 5 of a 2D array (restore to original dimensions)
int UnDecimate_by5_2d(float *src, float *dst, int ni, int li, int nj){
  int itup = (ni-2)/5 ;                  // number of averaged pairs along i
  int jtup = (nj-2)/5 ;                  // number of averaged pairs along j
  int ipre = ni - (5 * itup) - 1 ;
  int jpre = nj - (5 * jtup) - 1 ;
  int nid = ipre + 1 + itup ;            // total number of decimated points along i
  int i, j;
  float sixth2  = 2.0f / 6.0f ;
  float sixth4  = 4.0f / 6.0f ;
  float fifth   = 1.0f / 5.0f ;
  float fifth2  = 2.0f / 5.0f ;
  float fifth3  = 3.0f / 5.0f ;
  float fifth4  = 4.0f / 5.0f ;
  float *topd = dst + li * (nj -1 ) ;    // top row of undecimated array used as temporary

  for(j = 0 ; j < jpre ; j ++){
    UnDecimate_by5_1d(src, dst, ni) ;   // undecimate first pre rows (no processing along j)
    src += nid ;
    dst += li ;
  }

  // undecimate first group of rows
  dst[0] = 0.1f ; dst[li] = 0.2f ;
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * sixth2 + src[i-nid] * sixth4 ;      // along j
  UnDecimate_by5_1d(topd, dst, ni) ;                                              // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * sixth4 + src[i-nid] * sixth2  ;     // along j
  UnDecimate_by5_1d(topd, dst + li, ni) ;                                         // along i
  UnDecimate_by5_1d(src, dst + li + li, ni) ;                                     // along j and i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth4   + src[i+nid] * fifth  ;    // along j
  UnDecimate_by5_1d(topd, dst + li + li + li, ni) ;                               // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth3   + src[i+nid] * fifth2 ;    // along j
  UnDecimate_by5_1d(topd, dst + li + li + li+ li, ni) ;                           // along i
  dst += (li+li+li+li+li) ;
  src += nid ;

  // undecimate second group -> next to last  group (jtup -2 groups)
  for(j=2 ; j < jtup ; j++){
  dst[0] = 0.3f ; dst[li] = 0.4f ;
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth3 + src[i-nid] * fifth2 ;    // along j
    UnDecimate_by5_1d(topd, dst, ni) ;                                            // along i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth4 + src[i-nid] * fifth ;     // along j
    UnDecimate_by5_1d(topd, dst + li, ni) ;                                       // along i
    UnDecimate_by5_1d(src, dst + li + li, ni) ;                                   // along j and i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth4 + src[i+nid] * fifth ;     // along j
    UnDecimate_by5_1d(topd, dst + li + li + li, ni) ;                             // along i
    for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth3 + src[i+nid] * fifth2 ;    // along j
    UnDecimate_by5_1d(topd, dst + li + li + li + li, ni) ;                             // along i
    dst += (li+li+li+li+li) ;
    src += nid ;
  }

  // undecimate last group
  dst[0] = 0.1f ; dst[li] = 0.2f ;
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth3   + src[i-nid] * fifth2 ;    // along j
  UnDecimate_by5_1d(topd, dst, ni) ;                                              // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * fifth4   + src[i-nid] * fifth  ;    // along j
  UnDecimate_by5_1d(topd, dst + li, ni) ;                                         // along i
  UnDecimate_by5_1d(src, dst + li + li, ni) ;                                     // along j and i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * sixth4 + src[i+nid] * sixth2 ;      // along j
  UnDecimate_by5_1d(topd, dst + li + li + li, ni) ;                               // along i
  for(i=0 ; i < nid ; i++) topd[i] = src[i] * sixth2 + src[i+nid] * sixth4 ;      // along j
  UnDecimate_by5_1d(topd, dst + li + li + li + li, ni) ;                          // along i
  dst += (li+li+li+li+li) ;
  src += nid ;

  // undecimate last row (no processing along j)
  UnDecimate_by5_1d(src, dst, ni) ;
  return 0 ;
}

// 2D general decimation function
int Decimate_2d(float *src, int factor, float *dst, int ni, int li, int nj){
  switch(factor)
  {
    case 2:
      return Decimate_by2_2d(src, dst, ni, li, nj) ;
    case 3:
      return Decimate_by3_2d(src, dst, ni, li, nj) ;
    case 4:
      return Decimate_by4_2d(src, dst, ni, li, nj) ;
    case 5:
      return Decimate_by5_2d(src, dst, ni, li, nj) ;
    default:
      return -1 ;
  }
}

// 2D general inverses decimation (restore) function
int Undecimate_2d(float *src, int factor, float *dst, int ni, int li, int nj){
  switch(factor)
  {
    case 2:
      return UnDecimate_by2_2d(src, dst, ni, li, nj) ;
    case 3:
      return UnDecimate_by3_2d(src, dst, ni, li, nj) ;
    case 4:
      return UnDecimate_by4_2d(src, dst, ni, li, nj) ;
    case 5:
      return UnDecimate_by5_2d(src, dst, ni, li, nj) ;
    default:
      return -1 ;
  }
}

#if defined(SELF_TEST)
#if ! defined(NI)
#define NI 12
#endif
#if ! defined(NJ)
#define NJ 12
#endif

#include <stdio.h>
int main(int argc, char **argv){
  float src1[NI] ;
  float src2[NI] ;
  float dst1[NI * NJ] ;
  float srca[NJ * NI] ;
  float srcb[NJ * NI] ;
  int np, nd ;
  int i, j ;

  for(nd = 2 ; nd < 6 ; nd++){
    np = N_decimated(NI, nd) ;
    printf("decimating NI = %d points by %d => %d elements\n",NI,nd,np) ;
    np = N_decimated(NJ, nd) ;
    printf("decimating NJ = %d points by %d => %d elements\n",NJ,nd,np) ;
  }

  for(i=0 ; i < NI ; i++) src1[i] = i + 1.0f ;
  for(i=0 ; i < NI ; i++) src2[i] = 0.0f ;
  printf("========= original array =========\n") ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      srca[j*NI +i] = i + j + 1.0f ;
      printf("%4.1f ",srca[j*NI +i]) ;
    }
    printf("\n");
  }
//
// test decimation by 2
//
  printf("\n");
//   for(i=0 ; i < NI ; i++) printf("%4.1f ",src1[i]) ;
//   printf("\n");
  printf(" 1 D decimation by 2 (first row)\n");
  np = Decimate_by2_1d(srca, dst1, NI, 0) ;
  for(i=0 ; i < np ; i++) printf("%4.1f ",dst1[i]) ;
  printf("\n");
  UnDecimate_by2_1d(dst1, src2, NI) ;
  for(i=0 ; i < NI ; i++) printf("%4.1f ",src2[i]) ;
  printf("\n\n");

  printf(" 2 D decimation by 2\n");
  for(i=0 ; i<NI*NJ ; i++) dst1[i] = 0.0 ;
//   Decimate_by2_2d(srca, dst1, NI, NI, NJ) ;
  Decimate_2d(srca, 2, dst1, NI, NI, NJ) ;
  for(j=0 ; j<N_decimated(NJ,2) ; j++){
    for(i=0 ; i<N_decimated(NI,2) ; i++){
      printf("%4.1f ",dst1[j*N_decimated(NI,2) +i]) ;
    }
    printf("\n");
  }
  printf("\n");
  for(i=0 ; i<NI*NJ ; i++) srcb[i] = 0.0 ;
  UnDecimate_by2_2d(dst1, srcb, NI, NI, NJ) ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      printf("%4.1f ",srcb[j*NI +i]) ;
    }
    printf("\n");
  }
  printf("\n");
//
// test decimation by 3
//
  printf(" 1 D decimation by 3 (first row)\n");
  for(i=0 ; i < NI ; i++) dst1[i] = 0.0f ;
  np = Decimate_by3_1d(src1, dst1, NI, 0) ;
  for(i=0 ; i < np ; i++) printf("%4.1f ",dst1[i]) ;
  printf("\n");

  for(i=0 ; i < NI ; i++) src2[i] = 0.0f ;
  UnDecimate_by3_1d(dst1, src2, NI) ;
  for(i=0 ; i < NI ; i++) printf("%4.1f ",src2[i]) ;
  printf("\n\n");

  printf(" 2 D decimation by 3\n");
  for(i=0 ; i<NI*NJ ; i++) dst1[i] = 0.0 ;
//   Decimate_by3_2d(srca, dst1, NI, NI, NJ) ;
  Decimate_2d(srca, 3, dst1, NI, NI, NJ) ;
  for(j=0 ; j<N_decimated(NJ,3) ; j++){
    for(i=0 ; i<N_decimated(NI,3) ; i++){
      printf("%4.1f ",dst1[j*N_decimated(NI,3) +i]) ;
    }
    printf("\n");
  }
  printf("\n");

  for(i=0 ; i<NI*NJ ; i++) srcb[i] = 0.0 ;
  UnDecimate_by3_2d(dst1, srcb, NI, NI, NJ) ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      printf("%4.1f ",srcb[j*NI +i]) ;
    }
    printf("\n");
  }
  printf("\n");
//
// test decimation by 4
//
  printf(" 1 D decimation by 4 (first row)\n");
  for(i=0 ; i < NI ; i++) dst1[i] = 0.0f ;
  np = Decimate_by4_1d(src1, dst1, NI, 0) ;
  for(i=0 ; i < np ; i++) printf("%4.1f ",dst1[i]) ;
  printf("\n");

  for(i=0 ; i < NI ; i++) src2[i] = 0.0f ;
  UnDecimate_by4_1d(dst1, src2, NI) ;
  for(i=0 ; i < NI ; i++) printf("%4.1f ",src2[i]) ;
  printf("\n\n");

  printf(" 2 D decimation by 4\n");
  for(i=0 ; i<NI*NJ ; i++) dst1[i] = 0.0 ;
//   Decimate_by4_2d(srca, dst1, NI, NI, NJ) ;
  Decimate_2d(srca, 4, dst1, NI, NI, NJ) ;
  for(j=0 ; j<N_decimated(NJ,4) ; j++){
    for(i=0 ; i<N_decimated(NI,4) ; i++){
      printf("%4.1f ",dst1[j*N_decimated(NI,4) +i]) ;
    }
    printf("\n");
  }
  printf("\n");

  for(i=0 ; i<NI*NJ ; i++) srcb[i] = 0.0 ;
  UnDecimate_by4_2d(dst1, srcb, NI, NI, NJ) ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      printf("%4.1f ",srcb[j*NI +i]) ;
    }
    printf("\n");
  }
  printf("\n");
//
// test decimation by 5
//
  printf(" 1 D decimation by 5 (first row)\n");
  for(i=0 ; i < NI ; i++) dst1[i] = 0.0f ;
  np = Decimate_by5_1d(src1, dst1, NI, 0) ;
  for(i=0 ; i < np ; i++) printf("%4.1f ",dst1[i]) ;
  printf("\n");

  for(i=0 ; i < NI ; i++) src2[i] = 0.0f ;
  UnDecimate_by5_1d(dst1, src2, NI) ;
  for(i=0 ; i < NI ; i++) printf("%4.1f ",src2[i]) ;
  printf("\n\n");

  printf(" 2 D decimation by 5\n");
  for(i=0 ; i<NI*NJ ; i++) dst1[i] = 0.0 ;
//   Decimate_by5_2d(srca, dst1, NI, NI, NJ) ;
  Decimate_2d(srca, 5, dst1, NI, NI, NJ) ;
  for(j=0 ; j<N_decimated(NJ,5) ; j++){
    for(i=0 ; i<N_decimated(NI,5) ; i++){
      printf("%4.1f ",dst1[j*N_decimated(NI,5) +i]) ;
    }
    printf("\n");
  }
  printf("\n");
  for(i=0 ; i<NI*NJ ; i++) srcb[i] = 0.0 ;
  UnDecimate_by5_2d(dst1, srcb, NI, NI, NJ) ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      printf("%4.1f ",srcb[j*NI +i]) ;
    }
    printf("\n");
  }
  printf("\n");

  return 0 ;
}
#endif

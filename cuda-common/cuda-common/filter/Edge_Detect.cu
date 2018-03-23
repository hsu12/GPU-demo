# include <time .h>
# include <stdlib .h>
# include <stdio .h>
# include <string .h>
# include <math .h>
# include <cuda .h>
# include <cutil .h>
# include <ctime >


unsigned int width , height ;

 int Gx [3][3] = { -1 , 0 , 1 ,-2 , 0 , 2 ,-1 , 0 , 1};

 int Gy [3][3] = {1 ,2 ,1 , 0 ,0 ,0 , -1 , -2 , -1};

 int getPixel ( unsigned char * org , int col , int row) {

    int sumX , sumY ;
      sumX = sumY = 0;

        for (int i= -1; i <= 1; i++) {
            for (int j= -1; j <=1; j++) {
              int curPixel = org [( row + j) * width + (col + i) ];
              sumX += curPixel * Gx[i +1][ j +1];
              sumY += curPixel * Gy[i +1][ j +1];
          }
      }
      int sum = abs( sumY ) + abs( sumX ) ;
      if (sum > 255) sum = 255;
      if (sum < 0) sum = 0;
      return sum ;
 }

 void h_EdgeDetect ( unsigned char * org , unsigned char * result ) {
     int offset = 1 * width ;
     for (int row =1; row < height -2; row ++) {
     for (int col =1; col <width -2; col ++) {
        result [ offset + col ] = getPixel (org , col , row ) ;
      }
      offset += width ;
     }
  }


 __global__ void d_EdgeDetect ( unsigned char *org , unsigned char *result , int width , int height ) {
    int col = blockIdx .x * blockDim .x + threadIdx .x;
    int row = blockIdx .y * blockDim .y + threadIdx .y;

    if (row < 2 || col < 2 || row >= height -3 || col >= width -3 )
      return ;

    int Gx [3][3] = { -1 , 0 , 1 , -2 , 0 , 2 ,-1 , 0 , 1};

    int Gy [3][3] = {1 ,2 ,1 , 0 ,0 ,0 ,-1 , -2 , -1};

    int sumX , sumY ;
    sumX = sumY = 0;

    for (int i= -1; i <= 1; i++) {
      for (int j= -1; j <=1; j++) {
        int curPixel = org [( row + j) * width + (col + i) ];
        sumX += curPixel * Gx[i +1][ j +1];
        sumY += curPixel * Gy[i +1][ j +1];
      }
    }

    int sum = abs( sumY ) + abs( sumX ) ;
    if (sum > 255) sum = 255;
    if (sum < 0) sum = 0;

    result [row * width + col ] = sum ;

 }

 int main ( int argc , char ** argv )
 {
   printf (" Starting program \n") ;



   unsigned char * d_resultPixels ;
   unsigned char * h_resultPixels ;
   unsigned char * h_pixels = NULL ;
   unsigned char * d_pixels = NULL ;

   char * srcPath = "/ Developer /GPU Computing /C/src / EdgeDetection /
   image / cartoon .pgm";
   char * h_ResultPath = "/ Developer /GPU Computing /C/src /
   EdgeDetection / output / h_cartoon .pgm";
   char * d_ResultPath = "/ Developer /GPU Computing /C/src /
   EdgeDetection / output / d_cartoon .pgm";

   cutLoadPGMub ( srcPath , & h_pixels , &width , & height ) ;

   int ImageSize = sizeof ( unsigned char ) * width * height ;

   h_resultPixels = ( unsigned char *) malloc ( ImageSize ) ;
   cudaMalloc (( void **) & d_pixels , ImageSize ) ;
   cudaMalloc (( void **) & d_resultPixels , ImageSize ) ;
   cudaMemcpy ( d_pixels , h_pixels , ImageSize , cudaMemcpyHostToDevice) ;


   clock_t starttime , endtime , difference ;

   printf (" Starting host processing \n") ;
   starttime = clock () ;
   h_EdgeDetect ( h_pixels , h_resultPixels ) ;
   endtime = clock () ;
   printf (" Completed host processing \n") ;

   difference = ( endtime - starttime ) ;

   double interval = difference / ( double ) CLOCKS_PER_SEC ;
   printf ("CPU execution time = %f ms\n", interval * 1000) ;
   cutSavePGMub ( h_ResultPath , h_resultPixels , width , height ) ;

   dim3 block (16 ,16) ;
   dim3 grid ( width /16 , height /16) ;
   unsigned int timer = 0;





    printf (" Press enter to exit ...\ n") ;
      getchar () ;
 }

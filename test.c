#include "mat.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#pragma GCC optimize(3,"Ofast","inline")

int main(int argc, char *argv[]) {

  int matsize = 8000;

  char *ch1 = random_mat_file(

      matsize, matsize,

      get_mat_filename(matsize));
  char *ch2 = random_mat_file(

      matsize, matsize,

      get_mat_filename(matsize));

  Mat *mat1 = (Mat *)malloc(sizeof(Mat));
  Mat *mat2 = (Mat *)malloc(sizeof(Mat));
  mat1 = initMat_file(mat1, matsize, matsize, ch1);
  mat2 = initMat_file(mat2, matsize, matsize, ch2);

  // 矩阵乘法速度对比
  Mat *res = (Mat *)malloc(sizeof(Mat));
  res = initMat_random(res, matsize, matsize);
  struct timeval start, end;

  // //生成float数组进行测试

  // float *p1 = (float *)(aligned_alloc(256, sizeof(float) * 32 * 32));
  // float *p2 = (float *)(aligned_alloc(256, sizeof(float) * 32 * 32));
  // float *pres = (float *)(aligned_alloc(256, sizeof(float) * 32 * 32));
  // //随机初始化float数组
  // for (int i = 0; i < 32 * 32; i++) {
  //   p1[i] = (float)rand() / (float)RAND_MAX;
  //   p2[i] = (float)rand() / (float)RAND_MAX;
  // }

  // float_vectorAdd_avx2(pres, p1, p2, 32 * 32);
  // float_vectorSub_avx2(pres, p1, p2, 32 * 32);
  // float_vectorMul_avx2(pres, p1, p2, 32 * 32);

  // loat_matAdd_avx2(res, mat1, mat2);
  //  matAdd(res, mat1, mat2);
  //  float_matSub_avx2(res, mat1, mat2);
  matmul_ikj_avx2(res, mat1, mat2);
  // matmul_plain(res, mat1, mat2);

  // gettimeofday(&start, NULL);
  // matmul_ikj_avx2(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("ikj_avx2 Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));

  printMat(res);
  // for (int i = 0; i < 32 * 32; i++) {
  //   printf("%f ", pres[i]);
  // }

  // free(p1);
  // free(p2);
  // free(pres);

  return 0;
}
// gcc -o test test.c mat.c matmul.c -lm -fgnu89-inline -mavx2 -DWITH_AVX2 && ./test





//   void version5(int mat1[N][N], int mat2[N][N], int result[N][N])
// {
// __m256i vec_multi_res = _mm256_setzero_si256(); //Initialize vector to zero
// __m256i vec_mat1 = _mm256_setzero_si256(); //Initialize vector to zero
// __m256i vec_mat2 = _mm256_setzero_si256(); //Initialize vector to zero

// int i, j, k;
// for (i = 0; i < N; i++)
// {
//     for (j = 0; j < N; ++j)
//     {
//         //Stores one element in mat1 and use it in all computations needed before proceeding
//         //Stores as vector to increase computations per cycle
//         vec_mat1 = _mm256_set1_epi32(mat1[i][j]);

//         for (k = 0; k < N; k += 8)
//         {
//             vec_mat2 = _mm256_loadu_si256((__m256i*)&mat2[j][k]); //Stores row of second matrix (eight in each iteration)
//             vec_multi_res = _mm256_loadu_si256((__m256i*)&result[i][k]); //Loads the result matrix row as a vector
//             vec_multi_res = _mm256_add_epi32(vec_multi_res ,_mm256_mullo_epi32(vec_mat1, vec_mat2));//Multiplies the vectors and adds to th the result vector

//             _mm256_storeu_si256((__m256i*)&result[i][k], vec_multi_res); //Stores the result vector into the result array
//         }
//     }
//     }
// }

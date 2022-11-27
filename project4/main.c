#include "mat.h"
#include </opt/OpenBLAS/include/cblas.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#pragma GCC optimize(3, "Ofast", "inline")

int main(int argc, char *argv[]) {

  int matsize = 1024;

  printf("Matrix size: %d\n", matsize);

  //     char *ch1 = random_mat_file(

  //         matsize, matsize,

  //         get_mat_filename(matsize));
  //     char *ch2 = random_mat_file(

  //         matsize, matsize,

  //         get_mat_filename(matsize));

  char *ch1 = get_mat_filename(matsize);
  char *ch2 = get_mat_filename(matsize);

  Mat *mat1 = (Mat *)malloc(sizeof(Mat));
  Mat *mat2 = (Mat *)malloc(sizeof(Mat));
  mat1 = initMat_file(mat1, matsize, matsize, ch1);
  mat2 = initMat_file(mat2, matsize, matsize, ch2);

  // 矩阵乘法速度对比
  Mat *res = (Mat *)malloc(sizeof(Mat));
  Mat *res2 = (Mat *)malloc(sizeof(Mat));
  Mat *res3 = (Mat *)malloc(sizeof(Mat));
  res = initMat(res, matsize, matsize);
  res2 = initMat(res2, matsize, matsize);
  res3 = initMat(res3, matsize, matsize);
  struct timeval start, end;

  // // 朴素算法
  // gettimeofday(&start, NULL);
  // matmul_plain(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("plain Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //        start.tv_usec));
  // printMat(res);
  //   访存优化ikj
  //   gettimeofday(&start, NULL);
  //   matmul_ikj(res, mat1, mat2);
  //   gettimeofday(&end, NULL);
  //   printf("ikj Time cost: %ld us\n",
  //          (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //          start.tv_usec));
  //   printMat(res);
  // // 访存优化kij
  // gettimeofday(&start, NULL);
  // matmul_kij(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("kij Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //        start.tv_usec));
  // // 分块算法
  // gettimeofday(&start, NULL);
  // matmul_block(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("block Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //        start.tv_usec));

  //   // 分块算法,访存优化
  //   gettimeofday(&start, NULL);
  //   matmul_block_opt(res, mat1, mat2);
  //   gettimeofday(&end, NULL);
  //   printf("block_opt Time cost: %ld us\n",
  //          (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //          start.tv_usec));

  // Strassen算法
  //   gettimeofday(&start, NULL);
  //   matmul_strassen(res, mat1, mat2);
  //   gettimeofday(&end, NULL);
  //   printf("strassen Time cost: %ld us\n",
  //          (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //          start.tv_usec));
  // printMat(res);

  // // 访存优化ikj，使用avx2指令集
  // gettimeofday(&start, NULL);
  // matmul_ikj_avx2(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("ikj_avx2 Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //        start.tv_usec));
  // //    printMat(res);

  // // 访存优化ikj，使用avx2指令集,OPENMP并行
  // gettimeofday(&start, NULL);
  // matmul_ikj_avx2_omp(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("ikj_avx2_omp Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //        start.tv_usec));
  // //  printMat(res);

  // // Strassen算法，使用avx2指令集
  // gettimeofday(&start, NULL);
  // matmul_strassen_avx2(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("strassen_avx2 Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //        start.tv_usec));
  // //  printMat(res);

  // // Strassen算法，使用avx2指令集,OPENMP并行
  // gettimeofday(&start, NULL);
  // matmul_strassen_avx2_omp(res, mat1, mat2);
  // gettimeofday(&end, NULL);
  // printf("strassen_avx2_omp Time cost: %ld us\n",
  //        (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec -
  //        start.tv_usec));
  //  printMat(res);

  // 矩阵乘法，访存优化，AVX2指令集，OpenMP并行，分块,循环展开
  gettimeofday(&start, NULL);
  matmul_block_jik_avx2_omp(res, mat1, mat2);
  gettimeofday(&end, NULL);
  printf("block_ijk_avx2_omp Time cost: %ld us\n",
         (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
  // printMat(res);

  // OPENBLAS
  gettimeofday(&start, NULL);
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mat1->row, mat2->col,
              mat1->col, 1, mat1->data, mat1->col, mat2->data, mat2->col, 0,
              res2->data, res2->col);
  gettimeofday(&end, NULL);
  printf("openblas Time cost: %ld us\n",
         (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));
  // printMat(res2);

//   matSub(res3, res2, res);
//   float_vectorDiv_avx2_omp(res3->data, res3->data, res2->data, res2->row * res2->col);
//   printf("The Max difference rate is %f\n", matMax(res3));
// //   printMat(res3);

  freeMat(res);
  freeMat(res2);
  freeMat(res3);
  freeMat(mat1);
  freeMat(mat2);

  return 0;
}
// gcc -o main main.c mat.c matmul.c -lm -fgnu89-inline -mavx2 -DWITH_AVX2 -fopenmp -L/opt/OpenBLAS/lib -lopenblas  && ./main

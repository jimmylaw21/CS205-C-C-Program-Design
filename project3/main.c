#include "mat.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#pragma GCC optimize(2)
int main(int argc, char *argv[]) {
  char *ch1 = random_mat_file(

      1024, 1024,

      "/mnt/c/Users/jimmylaw21/OneDrive - "

      "南方科技大学/桌面/主要文件/CPP-main/project4/mat-A-1024.txt");
  char *ch2 = random_mat_file(

      1024, 1024,

      "/mnt/c/Users/jimmylaw21/OneDrive - "

      "南方科技大学/桌面/主要文件/CPP-main/project4/mat-B-1024.txt");
      
  float *mata = random_mat_ptr(1024, 1024);
  float *matb = random_mat_ptr(1024, 1024);
  Mat *mat1 = (Mat *)malloc(sizeof(Mat));
  Mat *mat2 = (Mat *)malloc(sizeof(Mat));
  mat1 = initMat_file(mat1, 1024, 1024, ch1);
  mat2 = initMat_file(mat2, 1024, 1024, ch2);
  Mat *mat3 = (Mat *)malloc(sizeof(Mat));
  Mat *mat4 = (Mat *)malloc(sizeof(Mat));
  mat3 = initMat_array(mat3, 4, 4, mata);
  mat4 = initMat_array(mat4, 4, 4, matb);
  Mat *mat5 = (Mat *)malloc(sizeof(Mat));
  Mat *mat6 = (Mat *)malloc(sizeof(Mat));
  mat5 = initMat_random(mat5, 4, 4);
  mat6 = initMat_random(mat6, 4, 4);
  // Testcase 1 三种方法生成浮点数矩阵
  // printMat(mat1);
  // printf("\n");
  // printMat(mat2);
  // printf("\n");
  // printMat(mat3);
  // printf("\n");
  // printMat(mat4);
  // printf("\n");
  // printMat(mat5);
  // printf("\n");
  // printMat(mat6);
  // Testcase 2 复制，释放和删除矩阵
  // if (cpyMat(mat5, mat6)) {
  // printf("Copy successfully\n");
  // printMat(mat5);
  // printf("\n");
  // printMat(mat6);
  // } else {
  // printf("Copy failed\n");
  // }
  // if (delMat(&mat5)) {
  // printf("Delete successfully\n");
  // } else {
  // printf("Delete failed\n");
  // }
  // printMat(mat5);
  // if (freeMat(mat6)) {
  // printf("Free successfully\n");
  // } else {
  // printf("Free failed\n");
  // }
  // printMat(mat6);
  // Testcase 3 矩阵求极值
  // printf("Max of mat1 is %f\n", matMax(mat1));
  // printf("Min of mat1 is %f\n", matMin(mat1));
  // Testcase 4 矩阵和标量加减乘除
  // Mat *res = (Mat *)malloc(sizeof(Mat));
  // res = initMat_random(res, 4, 4);
  // printMat(res);
  // matAddScalar(res,mat1, 1.0f);
  // printMat(res);
  // matSubScalar(res, mat1, 1.0f);
  // printMat(res);
  // matMulScalar(res, mat1, 2.0f);
  // printMat(res);
  // matDivScalar(res, mat1, 2.0f);
  // printMat(res);
  // Testcase 5 访存硬件级优化对矩阵乘法速度对比
  Mat *res = (Mat *)malloc(sizeof(Mat));
  res = initMat_random(res, 1024, 1024);
  // printMat(mat1);
  // printMat(mat2);
  struct timeval start1, end1;
  struct timeval start2, end2;

  gettimeofday(&start1, NULL);
  matMul(res, mat1, mat2);
  gettimeofday(&end1, NULL);
  printf("OPT Time cost: %ld us\n", (end1.tv_sec - start1.tv_sec) * 1000000 +

                                        (end1.tv_usec - start1.tv_usec));
  // printMat(res);
  gettimeofday(&start2, NULL);
  matMulNaive(res, mat1, mat2);
  gettimeofday(&end2, NULL);
  printf("Naive Time cost: %ld us\n", (end2.tv_sec - start2.tv_sec) * 1000000 +

                                          (end2.tv_usec - start2.tv_usec));
  // printMat(res);
  // Testcase6 矩阵转置
  // Mat *res = (Mat *)malloc(sizeof(Mat));
  // res = initMat_random(res, 4, 4);
  // printMat(mat1);
  // matTranspose(res, mat1);
  // printMat(res);
  // Testcase7 矩阵求逆
  // Mat *res = (Mat *)malloc(sizeof(Mat));
  // res = initMat_random(res, 4, 4);
  // printMat(mat1);
  // matInverse(res, mat1);
  // printMat(res);
  return 0;
}
// gcc -o main main.c mat.c -lm -fgnu89-inline && ./main

#pragma once
#include <stddef.h>
#ifndef _MATMUL_F
#include <math.h>
#include <stdio.h>

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#define _MATMUL_F
#pragma GCC optimize(3,"Ofast","inline")

typedef struct Matrix {
  size_t row;
  size_t col;
  float *data;
} Mat;

Mat *initMat(Mat *src, int r, int c); // 初始化矩阵

Mat *initMat_file(Mat *src, int r, int c,
                  char *filename); // 用文件初始化矩阵

Mat *initMat_array(Mat *mat, int row, int col,
                   float *data); // 用float指针初始化矩阵

Mat *initMat_random(Mat *mat, int row, int col); // 用随机数初始化矩阵

Mat *readMat_file(Mat *mat, int r, int c, char *filename); // 从文件读取矩阵

char *get_mat_filename(int matsize);

bool freeMat(Mat *mat); // 释放矩阵

bool delMat(Mat **mat); // 删除矩阵

bool cpyMat(Mat *mat1, Mat *mat2); // 复制矩阵

float matMax(const Mat *mat); // 获取矩阵最大值

float matMin(const Mat *mat); // 获取矩阵最小值

bool matAdd(Mat *result, const Mat *mat1,
            const Mat *mat2); // 矩阵加法

bool matSub(Mat *result, const Mat *mat1,
            const Mat *mat2); // 矩阵减法

bool matMul(Mat *result, const Mat *mat1,
            const Mat *mat2); // 矩阵乘法，访存优化

bool matAddScalar(Mat *result, const Mat *mat,
                  float scalar); // 矩阵加标量

bool matSubScalar(Mat *result, const Mat *mat,
                  float scalar); // 矩阵减标量

bool matMulScalar(Mat *result, const Mat *mat,
                  float scalar); // 矩阵乘标量

bool matDivScalar(Mat *result, const Mat *mat,
                  float scalar); // 矩阵除标量

// int matSize_file(char str[]); // 通过文件名获取矩阵规模

void random_mat(Mat *mat); // 给矩阵随机赋值

char *random_mat_file(int row, int col, char *filename); // 随机生成矩阵数据文件

float *random_mat_ptr(int row, int col); // 随机生成矩阵数据指针

void printMat(const Mat *mat); // 打印矩阵

bool isMatNull(const Mat *mat); // 判断矩阵是否为空

bool matTranspose(Mat *result, const Mat *mat); // 矩阵转置

float matDet(const Mat *mat); // 矩阵行列式

bool isMatSingular(const Mat *mat); // 判断矩阵是否奇异

bool matInverse(Mat *result, const Mat *mat); // 矩阵求逆

bool getMatUnit(Mat *result, const int n); // 获取单位矩阵

inline float Max(const float a, const float b); // 获取最大值

inline float Min(const float a, const float b); // 获取最小值

// void splitMat(Mat *result, Mat *mat, int row, int col, int rowLen, int colLen); // 切分矩阵

// void mergeMat(Mat *mat, Mat *subMat, int row, int col); // 合并矩阵

bool matmul_plain(Mat *result, const Mat *mat1,
                  const Mat *mat2); // 矩阵乘法，朴素算法

bool matmul_ikj(Mat *result, const Mat *mat1,
                const Mat *mat2); // 矩阵乘法，访存优化ikj

bool matmul_kij(Mat *result, const Mat *mat1,
                const Mat *mat2); // 矩阵乘法，访存优化kij

bool matmul_block(Mat *result, const Mat *mat1,
                  const Mat *mat2); // 矩阵乘法，分块算法

bool matmul_block_opt(Mat *result, const Mat *mat1,
                      const Mat *mat2); // 矩阵乘法，分块算法，访存优化

bool matmul_strassen(Mat *result, const Mat *mat1,
                     const Mat *mat2); // 矩阵乘法，Strassen算法

bool matmul_ikj_avx2(Mat *result, const Mat *mat1,
                       const Mat *mat2); // 矩阵乘法，访存优化，AVX2

bool matmul_ikj_avx2_omp(Mat *result, const Mat *mat1, const Mat *mat2);// 矩阵乘法，访存优化，AVX2，OpenMP


bool matmul_strassen_avx2(Mat *result, const Mat *mat1,
                          const Mat *mat2); // 矩阵乘法，Strassen算法，AVX2指令集

bool matmul_strassen_avx2_omp(Mat *result, const Mat *mat1,
                            const Mat *mat2); // 矩阵乘法，Strassen算法，AVX2指令集，OpenMP并行

bool matmul_block_jik_avx2_omp( Mat *result, const Mat *mat1, const Mat *mat2); // 矩阵乘法，访存优化，AVX2指令集，OpenMP并行，分块

float *float_vectorAdd_avx2(float *res, const float *p1, const float *p2,size_t n);

float *float_vectorSub_avx2(float *res, const float *p1, const float *p2,size_t n);

float *float_vectorMul_avx2(float *res, const float *p1, const float *p2,size_t n);

float *float_vectorAdd_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n);

float *float_vectorSub_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n);

float *float_vectorMul_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n); 

float *float_vectorDiv_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n);

bool float_matAdd_avx2(Mat *result, const Mat *mat1, const Mat *mat2);

bool float_matSub_avx2(Mat *result, const Mat *mat1, const Mat *mat2);

bool float_matMul_avx2(Mat *result, const Mat *mat1, const Mat *mat2);


#endif
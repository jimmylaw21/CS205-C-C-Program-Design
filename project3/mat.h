#ifndef _MATMUL_F
#include <math.h>
#include <stdio.h>

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#define _MATMUL_F

typedef struct Matrix {
  int row;
  int col;
  float *data;
} Mat;

Mat *initMat(Mat *src, int r, int c); // 初始化矩阵

Mat *initMat_file(Mat *src, int r, int c,
                  char *filename); // 用文件初始化矩阵

Mat *initMat_array(Mat *mat, int row, int col,
                   float *data); // 用float指针初始化矩阵

Mat *initMat_random(Mat *mat, int row, int col); // 用随机数初始化矩阵

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

bool matMulNaive(Mat *result, const Mat *mat1,
                 const Mat *mat2); // 矩阵乘法，朴素算法

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

#endif
#include "mat.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

// 初始化矩阵
Mat *initMat(Mat *mat, int row, int col) {
  if (row <= 0 || col <= 0) {
    printf("The line or colomn input is improper!\n");
    return NULL;
  }
  // 为矩阵结构各结构体变量赋值
  mat->row = row;
  mat->col = col;
  // 动态分配二维数组的内存
  mat->data = (float *)malloc(sizeof(float) * row * col);
  return mat;
}

// 初始化矩阵，从文件中读取矩阵
Mat *initMat_file(Mat *mat, int row, int col, char *filename) {
  if (row <= 0 || col <= 0) {
    printf("The line or colomn input is improper!\n");
    return NULL;
  }
  if (filename == NULL) {
    printf("The file name input is improper!\n");
    return NULL;
  }
  // 初始化矩阵
  mat = initMat(mat, row, col);
  // 给矩阵动态分配空间之后，读取相应的数据存储到矩阵中
  if (mat != NULL && mat->row == row && mat->col == col) {
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
      printf("File open failed!\n");
      exit(1);
      return NULL;
    }
    for (int i = 0; i < row * col; i++) {
      fscanf(fp, "%f", &mat->data[i]);
    }
    fclose(fp);
  } else {
    printf("The matrix is not initialized properly!\n");
  }
  return mat;
}

// 初始化矩阵，并通过float指针读取矩阵
Mat *initMat_array(Mat *mat, int row, int col, float *data) {
  if (row <= 0 || col <= 0) {
    printf("The line or colomn input is improper!\n");
    return NULL;
  }
  if (data == NULL) {
    printf("The array data input is improper!\n");
    return NULL;
  }
  // 初始化矩阵
  mat = initMat(mat, row, col);
  // 给矩阵动态分配空间之后，读取相应的数据存储到矩阵中
  if (mat != NULL && mat->row == row && mat->col == col) {
    for (int i = 0; i < row * col; i++) {
      mat->data[i] = data[i];
    }
  } else {
    printf("The matrix is not initialized properly!\n");
  }
  return mat;
}

// 初始化矩阵，并生成随机矩阵
Mat *initMat_random(Mat *mat, int row, int col) {
  if (row <= 0 || col <= 0) {
    printf("The line or colomn input is improper!\n");
    return NULL;
  }
  // 初始化矩阵
  mat = initMat(mat, row, col);
  mat->data = random_mat_ptr(row, col);
  return mat;
}

// 释放矩阵
bool freeMat(Mat *mat) {
  if (mat == NULL || mat->data == NULL) {
    return false;
  }
  free(mat->data);
  mat->data = NULL;
  free(mat);
  mat = NULL;
  return true;
}

// 删除矩阵
bool delMat(Mat **mat) {
  if (*mat == NULL || (*mat)->data == NULL) {
    return false;
  }
  free((*mat)->data);
  (*mat)->data = NULL;
  free(*mat);
  *mat = NULL;
  return true;
}

// 复制矩阵,将source矩阵复制到target矩阵中
bool cpyMat(Mat *target, Mat *source) {
  if (source == NULL || source->row <= 0 || source->col <= 0) {
    return false;
  }
  // 如果target矩阵为空，或者target矩阵的行列数与source矩阵不同，则重新初始化target矩阵
  if (target == NULL || source->row != target->row ||
      source->col != target->col) {
    if (!freeMat(target)) {
      printf("Fail freeing matrix");
    }
    Mat *tar = (Mat *)malloc(sizeof(Mat));
    target = initMat(tar, source->row, source->col);
  }
  // 将source矩阵的数据复制到target矩阵中
  for (int i = 0; i < source->row * source->col; i++) {
    target->data[i] = source->data[i];
  }
  return true;
}

float matMax(const Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 获取矩阵的最大值
  float max = mat->data[0];
  for (int i = 0; i < mat->row; i++) {
    for (int j = 0; j < mat->col; j++) {
      max = Max(mat->data[i * mat->col + j], max);
    }
  }
  return max;
}

float matMin(const Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 获取矩阵的最小值
  float min = mat->data[0];
  for (int i = 0; i < mat->row; i++) {
    for (int j = 0; j < mat->col; j++) {
      min = Min(mat->data[i * mat->col + j], min);
    }
  }
  return min;
}

// 矩阵的加法
bool matAdd(Mat *result, const Mat *mat1, const Mat *mat2) {

  // 检查矩阵为空
  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }
  // 检查矩阵的行列数是否相同
  if (mat1->row != mat2->row || mat1->col != mat2->col) {
    printf("The matrix is not the same size!\n");
    return false;
  }
  // 初始化结果矩阵
  Mat *res = (Mat *)malloc(sizeof(Mat));
  result = initMat(res, mat1->row, mat1->col);
  // 矩阵加法
  for (int i = 0; i < mat1->row * mat1->col; i++) {
    result->data[i] = mat1->data[i] + mat2->data[i];
  }
  return true;
}

// 矩阵的减法
bool matSub(Mat *result, const Mat *mat1, const Mat *mat2) {
  // 检查矩阵为空
  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }
  // 检查矩阵的行列数是否相同
  if (mat1->row != mat2->row || mat1->col != mat2->col) {
    printf("The matrix is not the same size!\n");
    return false;
  }
  // 初始化结果矩阵
  Mat *res = (Mat *)malloc(sizeof(Mat));
  result = initMat(res, mat1->row, mat1->col);
  // 矩阵减法
  for (int i = 0; i < mat1->row * mat1->col; i++) {
    result->data[i] = mat1->data[i] - mat2->data[i];
  }
  return true;
}

// 矩阵的乘法,访存优化
bool matMul(Mat *result, const Mat *mat1, const Mat *mat2) {
  // 检查矩阵为空
  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }

  // 检查矩阵的行列数是否相同
  if (mat1->col != mat2->row) {
    printf("The matrix is not the same size!\n");
    return false;
  }

  // 初始化结果矩阵
  if (result == NULL || result->row != mat1->row || result->col != mat2->col) {
    if (!freeMat(result)) {
      printf("Fail freeing matrix");
    }
    Mat *res = (Mat *)malloc(sizeof(Mat));
    result = initMat(res, mat1->row, mat2->col);
  }

  // 矩阵乘法,访存优化
  int row = mat1->row;
  int col = mat1->col;

  for (int i = 0; i < row; i++) {
    for (int k = 0; k < row; k++) {
      float tmp = mat1->data[i * col + k];
      for (int j = 0; j < col; j++) {
        result->data[i * col + j] += tmp * mat2->data[k * col + j];
      }
    }
  }

  return true;
}

// 矩阵的乘法，朴素算法
bool matMulNaive(Mat *result, const Mat *mat1, const Mat *mat2) {
  // 检查矩阵为空
  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }
  // 检查矩阵的行列数是否相同
  if (mat1->col != mat2->row) {
    printf("The matrix is not the same size!\n");
    return false;
  }
  // 初始化结果矩阵
  if (result == NULL || result->row != mat1->row || result->col != mat2->col) {
    if (!freeMat(result)) {
      printf("Fail freeing matrix");
    }
    Mat *res = (Mat *)malloc(sizeof(Mat));
    result = initMat(res, mat1->row, mat2->col);
  }

  // 矩阵乘法
  int row = mat1->row;
  int col = mat1->col;
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      for (int k = 0; k < col; k++) {
        result->data[i * col + j] += mat1->data[i * col + k] * mat2->data[k * col + j];
      }
    }
  }
  return true;
}

// 矩阵与标量相加
bool matAddScalar(Mat *result, const Mat *mat, float scalar) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 矩阵与标量相加
  for (int i = 0; i < mat->row * mat->col; i++) {
    result->data[i] = mat->data[i] + scalar;
  }
  return true;
}

// 矩阵与标量相减
bool matSubScalar(Mat *result, const Mat *mat, float scalar) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 矩阵与标量相减
  for (int i = 0; i < mat->row * mat->col; i++) {
    result->data[i] = mat->data[i] - scalar;
  }
  return true;
}

// 矩阵与标量相乘
bool matMulScalar(Mat *result, const Mat *mat, float scalar) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 矩阵与标量相乘
  for (int i = 0; i < mat->row * mat->col; i++) {
    result->data[i] = mat->data[i] * scalar;
  }
  return true;
}

// 矩阵与标量相除
bool matDivScalar(Mat *result, const Mat *mat, float scalar) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 矩阵与标量相除
  for (int i = 0; i < mat->row * mat->col; i++) {
    result->data[i] = mat->data[i] / scalar;
  }
  return true;
}

// 给矩阵随机赋值
void random_mat(Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return;
  }
  for (int i = 0; i < mat->row * mat->col; i++) {
    mat->data[i] = (rand() % (100 * 10 - 1)) / 10.0;
  }
}

// 随机生成矩阵数据文件
char *random_mat_file(int row, int col, char *filename) {
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("Open file failed!\n");
    return NULL;
  }
  // fprintf(fp, "%d %d  ", row, col);
  for (int i = 0; i < row * col; i++) {
    fprintf(fp, "%f ", (rand() % (100 * 10 - 1)) / 10.0);
  }
  fclose(fp);
  return filename;
}

// 随机生成矩阵数据float指针
float *random_mat_ptr(int row, int col) {
  float *data = (float *)malloc(sizeof(float) * row * col);
  for (int i = 0; i < row * col; i++) {
    data[i] = (rand() % (100 * 10 - 1)) / 10.0;
  }
  return data;
}

// 打印矩阵
void printMat(const Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return;
  }
  // 打印矩阵
  for (int i = 0; i < mat->row; i++) {
    for (int j = 0; j < mat->col; j++) {
      printf("%f ", mat->data[i * mat->col + j]);
    }
    printf("\n");
  }
  printf("\n");
}

// 获取数组的大小
// int matSize(char str[]){
//     int size = 0;
//     for(int i = 0; str[i] != '\0'; i++){
//         if(str[i] == ' '){
//             size++;
//         }
//     }
//     return size;
// }

// 矩阵转置
bool matTranspose(Mat *result, const Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 矩阵转置
  for (int i = 0; i < mat->row; i++) {
    for (int j = 0; j < mat->col; j++) {
      result->data[j * mat->row + i] = mat->data[i * mat->col + j];
    }
  }
  return true;
}

// 检查矩阵是否为空
bool isMatNull(const Mat *mat) {
  // 检查指针
  if (mat == NULL) {
    printf("The pointer to matrix is NULL!\n");
    return true;
  }
  // 检查矩阵是否为空
  if (mat->data == NULL) {
    printf("The pointer to matrix data is NULL!\n");
    return true;
  }
  // 检查矩阵是否为空
  if (mat->row <= 0 || mat->col <= 0) {
    printf("The matrix is NULL!\n");
    return true;
  }
  return false;
}

// 求矩阵行列式
float matDet(const Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return 0;
  }
  // 检查矩阵是否为方阵
  if (mat->row != mat->col) {
    printf("The matrix is not square matrix!\n");
    return 0;
  }
  // 求矩阵行列式
  float det = 0;       // 行列式
  if (mat->row == 1) { // 1阶行列式
    det = mat->data[0];
  } else if (mat->row == 2) { // 2阶行列式
    det = mat->data[0] * mat->data[3] - mat->data[1] * mat->data[2];
  } else {
    for (int i = 0; i < mat->col; i++) {                // 3阶以上行列式
      Mat *temp = (Mat *)malloc(sizeof(Mat));           // 临时矩阵
      temp = initMat(temp, mat->row - 1, mat->col - 1); // 初始化临时矩阵
      int k = 0;                                        // 临时矩阵的行
      for (int m = 1; m < mat->row; m++) {
        for (int n = 0; n < mat->col; n++) {
          if (n != i) {
            temp->data[k] = mat->data[m * mat->col + n]; // 赋值
            k++;                                         // 行加1
          }
        }
      }
      det += pow(-1, i) * mat->data[i] * matDet(temp);
      delMat(&temp);
    }
  }
  return det;
}

// 检查矩阵是否为奇异矩阵
bool isMatSingular(const Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 检查矩阵是否为方阵
  if (mat->row != mat->col) {
    printf("The matrix is not square matrix!\n");
    return false;
  }
  // 检查矩阵是否为奇异矩阵
  if (matDet(mat) == 0) {
    return true;
  }
  return false;
}

// 矩阵求逆
bool matInverse(Mat *result, const Mat *mat) {
  // 检查矩阵为空
  if (isMatNull(mat)) {
    return false;
  }
  // 检查矩阵是否为方阵
  if (mat->row != mat->col) {
    printf("The matrix is not square matrix!\n");
    return false;
  }
  // 检查矩阵是否为奇异矩阵
  if (isMatSingular(mat)) {
    printf("The matrix is singular matrix!\n");
    return false;
  }
  // 初始化结果矩阵
  result = initMat(result, mat->row, mat->col);
  // 矩阵求逆
  float det = matDet(mat);
  for (int i = 0; i < mat->row; i++) {
    for (int j = 0; j < mat->col; j++) {
      Mat *temp = (Mat *)malloc(sizeof(Mat));
      temp = initMat(temp, mat->row - 1, mat->col - 1);
      int k = 0;
      for (int m = 0; m < mat->row; m++) {
        for (int n = 0; n < mat->col; n++) {
          if (m != i && n != j) {
            temp->data[k] = mat->data[m * mat->col + n];
            k++;
          }
        }
      }
      result->data[j * mat->row + i] = pow(-1, i + j) * matDet(temp) / det;
    }
  }
  return true;
}

// 获得单位矩阵
bool getMatUnit(Mat *result, const int n) {
  // 检查矩阵为空
  if (isMatNull(result)) {
    return false;
  }
  // 检查矩阵是否为方阵
  if (result->row != result->col) {
    printf("The matrix is not square matrix!\n");
    return false;
  }
  // 检查矩阵大小
  if (result->row != n) {
    printf("The matrix is not unit matrix!\n");
    return false;
  }
  // 获得单位矩阵
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        result->data[i * n + j] = 1;
      } else {
        result->data[i * n + j] = 0;
      }
    }
  }
  return true;
}

// 求最大值
inline float Max(const float a, const float b) { return a > b ? a : b; }

// 求最小值
inline float Min(const float a, const float b) { return a < b ? a : b; }
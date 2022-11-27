#include "mat.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#pragma GCC optimize(3, "Ofast", "inline")

#define UNROLL 4
#define BLOCKSIZE 32

#ifdef WITH_AVX2
#include <immintrin.h>
#endif

#ifdef WITH_NEON
#include <arm_neon.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// 矩阵乘法，朴素算法
bool matmul_plain(Mat *result, const Mat *mat1, const Mat *mat2) {
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
  // if (result == NULL || result->row != mat1->row || result->col != mat2->col)
  // {
  //   if (!freeMat(result)) {
  //     printf("Fail freeing matrix");
  //   }
  //   Mat *res = (Mat *)malloc(sizeof(Mat));
  //   result = initMat(res, mat1->row, mat2->col);
  // }

  // 矩阵乘法
  int row = mat1->row;
  int col = mat1->col;
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      for (int k = 0; k < col; k++) {
        result->data[i * col + j] +=
            mat1->data[i * col + k] * mat2->data[k * col + j];
      }
    }
  }
  return true;
}

// 矩阵乘法,访存优化ikj
bool matmul_ikj(Mat *result, const Mat *mat1, const Mat *mat2) {
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
  // if (result == NULL || result->row != mat1->row || result->col != mat2->col)
  // {
  //   if (!freeMat(result)) {
  //     printf("Fail freeing matrix");
  //   }
  //   Mat *res = (Mat *)malloc(sizeof(Mat));
  //   result = initMat(res, mat1->row, mat2->col);
  // }

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

// 矩阵乘法,访存优化kij
bool matmul_kij(Mat *result, const Mat *mat1, const Mat *mat2) {
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
  // if (result == NULL || result->row != mat1->row || result->col != mat2->col)
  // {
  //   if (!freeMat(result)) {
  //     printf("Fail freeing matrix");
  //   }
  //   Mat *res = (Mat *)malloc(sizeof(Mat));
  //   result = initMat(res, mat1->row, mat2->col);
  // }

  // 矩阵乘法,访存优化
  int row = mat1->row;
  int col = mat1->col;

  for (int k = 0; k < row; k++) {
    for (int i = 0; i < row; i++) {
      float tmp = mat1->data[i * col + k];
      for (int j = 0; j < col; j++) {
        result->data[i * col + j] += tmp * mat2->data[k * col + j];
      }
    }
  }

  return true;
}

// 矩阵乘法,分块算法
bool matmul_block(Mat *result, const Mat *mat1, const Mat *mat2) {
  // 检查矩阵为空
  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }

  // 检查矩阵的行列数是否相同
  if (mat1->col != mat2->row) {
    printf("The matrix is not the same size!\n");
    return false;
  }

  // 矩阵乘法,分块算法
  int row = mat1->row;
  int col = mat1->col;
  for (int i = 0; i < row; i += BLOCKSIZE) {
    for (int j = 0; j < col; j += BLOCKSIZE) {
      for (int k = 0; k < col; k += BLOCKSIZE) {
        for (int ii = i; ii < i + BLOCKSIZE; ii++) {
          for (int jj = j; jj < j + BLOCKSIZE; jj++) {
            for (int kk = k; kk < k + BLOCKSIZE; kk++) {
              result->data[ii * col + jj] +=
                  mat1->data[ii * col + kk] * mat2->data[kk * col + jj];
            }
          }
        }
      }
    }
  }

  return true;
}

// 矩阵乘法,分块算法,访存优化
bool matmul_block_opt(Mat *result, const Mat *mat1, const Mat *mat2) {
  // 检查矩阵为空
  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }

  // 检查矩阵的行列数是否相同
  if (mat1->col != mat2->row) {
    printf("The matrix is not the same size!\n");
    return false;
  }

  // 矩阵乘法,分块算法
  int row = mat1->row;
  int col = mat1->col;

  for (int i = 0; i < row; i += BLOCKSIZE) {
    for (int j = 0; j < col; j += BLOCKSIZE) {
      for (int k = 0; k < col; k += BLOCKSIZE) {
        for (int ii = i; ii < i + BLOCKSIZE; ii++) {
          for (int kk = k; kk < k + BLOCKSIZE; kk++) {
            float tmp = mat1->data[ii * col + kk];
            for (int jj = j; jj < j + BLOCKSIZE; jj++) {
              result->data[ii * col + jj] += tmp * mat2->data[kk * col + jj];
            }
          }
        }
      }
    }
  }

  return true;
}

// 矩阵乘法,Strassen算法
bool matmul_strassen(Mat *result, const Mat *mat1, const Mat *mat2) {
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
  // if (result == NULL || result->row != mat1->row || result->col != mat2->col)
  // {
  //   if (!freeMat(result)) {
  //     printf("Fail freeing matrix");
  //   }
  //   Mat *res = (Mat *)malloc(sizeof(Mat));
  //   result = initMat(res, mat1->row, mat2->col);
  // }

  // 矩阵乘法,Strassen算法
  int row = mat1->row;
  int col = mat1->col;
  int size = mat1->row / 2;

  if (size <= 32) {
    matmul_ikj(result, mat1, mat2);
    return true;
  }

  Mat *A11 = (Mat *)malloc(sizeof(Mat));
  Mat *A12 = (Mat *)malloc(sizeof(Mat));
  Mat *A21 = (Mat *)malloc(sizeof(Mat));
  Mat *A22 = (Mat *)malloc(sizeof(Mat));

  initMat(A11, size, size);
  initMat(A12, size, size);
  initMat(A21, size, size);
  initMat(A22, size, size);

  Mat *B11 = (Mat *)malloc(sizeof(Mat));
  Mat *B12 = (Mat *)malloc(sizeof(Mat));
  Mat *B21 = (Mat *)malloc(sizeof(Mat));
  Mat *B22 = (Mat *)malloc(sizeof(Mat));

  initMat(B11, size, size);
  initMat(B12, size, size);
  initMat(B21, size, size);
  initMat(B22, size, size);

  Mat *C11 = (Mat *)malloc(sizeof(Mat));
  Mat *C12 = (Mat *)malloc(sizeof(Mat));
  Mat *C21 = (Mat *)malloc(sizeof(Mat));
  Mat *C22 = (Mat *)malloc(sizeof(Mat));

  initMat(C11, size, size);
  initMat(C12, size, size);
  initMat(C21, size, size);
  initMat(C22, size, size);

  Mat *P1 = (Mat *)malloc(sizeof(Mat));
  Mat *P2 = (Mat *)malloc(sizeof(Mat));
  Mat *P3 = (Mat *)malloc(sizeof(Mat));
  Mat *P4 = (Mat *)malloc(sizeof(Mat));
  Mat *P5 = (Mat *)malloc(sizeof(Mat));
  Mat *P6 = (Mat *)malloc(sizeof(Mat));
  Mat *P7 = (Mat *)malloc(sizeof(Mat));

  initMat(P1, size, size);
  initMat(P2, size, size);
  initMat(P3, size, size);
  initMat(P4, size, size);
  initMat(P5, size, size);
  initMat(P6, size, size);
  initMat(P7, size, size);

  Mat *S1 = (Mat *)malloc(sizeof(Mat));
  Mat *S2 = (Mat *)malloc(sizeof(Mat));
  Mat *S3 = (Mat *)malloc(sizeof(Mat));
  Mat *S4 = (Mat *)malloc(sizeof(Mat));
  Mat *S5 = (Mat *)malloc(sizeof(Mat));
  Mat *S6 = (Mat *)malloc(sizeof(Mat));
  Mat *S7 = (Mat *)malloc(sizeof(Mat));
  Mat *S8 = (Mat *)malloc(sizeof(Mat));
  Mat *S9 = (Mat *)malloc(sizeof(Mat));
  Mat *S10 = (Mat *)malloc(sizeof(Mat));

  initMat(S1, size, size);
  initMat(S2, size, size);
  initMat(S3, size, size);
  initMat(S4, size, size);
  initMat(S5, size, size);
  initMat(S6, size, size);
  initMat(S7, size, size);
  initMat(S8, size, size);
  initMat(S9, size, size);
  initMat(S10, size, size);

  Mat *resA = (Mat *)malloc(sizeof(Mat));
  Mat *resB = (Mat *)malloc(sizeof(Mat));

  initMat(resA, size, size);
  initMat(resB, size, size);

  // 对矩阵进行分块
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      A11->data[i * size + j] = mat1->data[i * row + j];
      A12->data[i * size + j] = mat1->data[i * row + j + size];
      A21->data[i * size + j] = mat1->data[(i + size) * row + j];
      A22->data[i * size + j] = mat1->data[(i + size) * row + j + size];

      B11->data[i * size + j] = mat2->data[i * col + j];
      B12->data[i * size + j] = mat2->data[i * col + j + size];
      B21->data[i * size + j] = mat2->data[(i + size) * col + j];
      B22->data[i * size + j] = mat2->data[(i + size) * col + j + size];
    }
  }

  // 计算S1~S10，P1~P7，C11~C22
  matSub(S1, B12, B22);
  matAdd(S2, A11, A12);
  matAdd(S3, A21, A22);
  matSub(S4, B21, B11);
  matAdd(S5, A11, A22);
  matAdd(S6, B11, B22);
  matSub(S7, A12, A22);
  matAdd(S8, B21, B22);
  matSub(S9, A11, A21);
  matAdd(S10, B11, B12);

  matmul_strassen(P1, A11, S1);
  matmul_strassen(P2, S2, B22);
  matmul_strassen(P3, S3, B11);
  matmul_strassen(P4, A22, S4);
  matmul_strassen(P5, S5, S6);
  matmul_strassen(P6, S7, S8);
  matmul_strassen(P7, S9, S10);

  matAdd(resA, P4, P5);
  matSub(resB, resA, P2);
  matAdd(C11, resB, P6);
  matAdd(C12, P1, P2);
  matAdd(C21, P3, P4);
  matAdd(resA, P5, P1);
  matSub(resB, resA, P3);
  matSub(C22, resB, P7);

  // 将C11~C22合并到C中
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      result->data[i * col + j] = C11->data[i * size + j];
      result->data[i * col + j + size] = C12->data[i * size + j];
      result->data[(i + size) * col + j] = C21->data[i * size + j];
      result->data[(i + size) * col + j + size] = C22->data[i * size + j];
    }
  }

  // 释放内存
  freeMat(A11);
  freeMat(A12);
  freeMat(A21);
  freeMat(A22);
  freeMat(B11);
  freeMat(B12);
  freeMat(B21);
  freeMat(B22);
  freeMat(C11);
  freeMat(C12);
  freeMat(C21);
  freeMat(C22);
  freeMat(P1);
  freeMat(P2);
  freeMat(P3);
  freeMat(P4);
  freeMat(P5);
  freeMat(P6);
  freeMat(P7);
  freeMat(S1);
  freeMat(S2);
  freeMat(S3);
  freeMat(S4);
  freeMat(S5);
  freeMat(S6);
  freeMat(S7);
  freeMat(S8);
  freeMat(S9);
  freeMat(S10);
  freeMat(resA);
  freeMat(resB);

  return true;
}

float *float_vectorAdd_avx2(float *res, const float *p1, const float *p2,
                            size_t n) {
#ifdef WITH_AVX2
  if (n % 8 != 0) {
    printf("The size n must be a multiple of 8.");
    return res;
  }

  __m256 a, b;
  for (size_t i = 0; i < n; i += 8) {
    a = _mm256_load_ps(p1 + i);
    b = _mm256_load_ps(p2 + i);
    _mm256_store_ps(res + i, _mm256_add_ps(a, b));
  }
  return res;
#else
  printf("AVX2 is not supported");
  return res;
#endif
}

float *float_vectorSub_avx2(float *res, const float *p1, const float *p2,
                            size_t n) {
#ifdef WITH_AVX2
  if (n % 8 != 0) {
    printf("The size n must be a multiple of 8.");
    return res;
  }

  __m256 a, b;
  for (size_t i = 0; i < n; i += 8) {
    a = _mm256_load_ps(p1 + i);
    b = _mm256_load_ps(p2 + i);
    _mm256_store_ps(res + i, _mm256_sub_ps(a, b));
  }
  return res;
#else
  printf("AVX2 is not supported");
  return res;
#endif
}

float *float_vectorMul_avx2(float *res, const float *p1, const float *p2,
                            size_t n) {
#ifdef WITH_AVX2
  if (n % 8 != 0) {
    printf("The size n must be a multiple of 8.");
    return res;
  }

  __m256 a, b;
  for (size_t i = 0; i < n; i += 8) {
    a = _mm256_load_ps(p1 + i);
    b = _mm256_load_ps(p2 + i);
    _mm256_store_ps(res + i, _mm256_mul_ps(a, b));
  }
  return res;
#else
  printf("AVX2 is not supported");
  return res;
#endif
}

bool float_matAdd_avx2(Mat *result, const Mat *mat1, const Mat *mat2) {
#ifdef WITH_AVX2

  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }
  if (mat1->row != mat2->row || mat1->col != mat2->col) {
    printf("The size of two matrix must be the same.");
    return false;
  }

  float *p1 = mat1->data;
  float *p2 = mat2->data;
  float *res = result->data;
  size_t n = mat1->row * mat1->col;

  float_vectorAdd_avx2(res, p1, p2, n);
  return true;
#else
  printf("AVX2 is not supported");
  return false;
#endif
}

bool float_matSub_avx2(Mat *result, const Mat *mat1, const Mat *mat2) {
#ifdef WITH_AVX2
  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }
  if (mat1->row != mat2->row || mat1->col != mat2->col) {
    printf("The size of two matrix must be the same.");
    return false;
  }

  float *p1 = mat1->data;
  float *p2 = mat2->data;
  float *res = result->data;
  size_t n = mat1->row * mat1->col;

  float_vectorSub_avx2(res, p1, p2, n);
  return true;
#else
  printf("AVX2 is not supported");
  return false;
#endif
}

// 矩阵乘法，访存优化，使用AVX2指令集
bool matmul_ikj_avx2(Mat *result, const Mat *mat1, const Mat *mat2) {
#ifdef WITH_AVX2
  if (isMatNull(mat1) || isMatNull(mat2) || isMatNull(result)) {
    return false;
  }
  if (mat1->col != mat2->row) {
    printf("The size of two matrix is not suitable.");
    return false;
  }
  if (mat1->row != result->row || mat2->col != result->col) {
    printf("The size of result matrix is not suitable.");
    return false;
  }

  float *p1 = mat1->data;
  float *p2 = mat2->data;
  float *res = result->data;
  size_t n = mat1->row;

  __m256 a, b, c;

  for (size_t i = 0; i < n; i++) {
    for (size_t k = 0; k < n; k++) {
      a = _mm256_set1_ps(p1[i * n + k]);
      for (size_t j = 0; j < n; j += 8) {
        b = _mm256_loadu_ps(p2 + k * n + j);
        c = _mm256_loadu_ps(res + i * n + j);
        _mm256_storeu_ps(res + i * n + j,
                         _mm256_add_ps(c, _mm256_mul_ps(a, b)));
      }
    }
  }

  return true;
#else
  printf("AVX2 is not supported");
  return false;
#endif
}

// 矩阵乘法，访存优化，使用AVX2指令集，使用OPENMP多线程
bool matmul_ikj_avx2_omp(Mat *result, const Mat *mat1, const Mat *mat2) {
#ifdef WITH_AVX2
  if (isMatNull(mat1) || isMatNull(mat2) || isMatNull(result)) {
    return false;
  }
  if (mat1->col != mat2->row) {
    printf("The size of two matrix is not suitable.");
    return false;
  }
  if (mat1->row != result->row || mat2->col != result->col) {
    printf("The size of result matrix is not suitable.");
    return false;
  }

  size_t THREAD_NUM = 8;

  float *p1 = mat1->data;
  float *p2 = mat2->data;
  float *res = result->data;
  size_t n = mat1->row;

  __m256 a, b, c;

#pragma omp parallel for num_threads(THREAD_NUM) private(a, b, c)
  for (size_t i = 0; i < n; i++) {
    for (size_t k = 0; k < n; k++) {
      a = _mm256_set1_ps(p1[i * n + k]);
      for (size_t j = 0; j < n; j += 8) {
        b = _mm256_loadu_ps(p2 + k * n + j);
        c = _mm256_loadu_ps(res + i * n + j);
        _mm256_storeu_ps(res + i * n + j,
                         _mm256_add_ps(c, _mm256_mul_ps(a, b)));
      }
    }
  }

  return true;
#else
  printf("AVX2 is not supported");
  return false;
#endif
}

bool matmul_block_jik_avx2_omp(Mat *result, const Mat *mat1, const Mat *mat2) {
#ifdef WITH_AVX2
  if (isMatNull(mat1) || isMatNull(mat2) || isMatNull(result)) {
    return false;
  }
  if (mat1->col != mat2->row) {
    printf("The size of two matrix is not suitable.");
    return false;
  }
  if (mat1->row != result->row || mat2->col != result->col) {
    printf("The size of result matrix is not suitable.");
    return false;
  }
  float *A = mat1->data;
  float *B = mat2->data;
  float *C = result->data;
  size_t n = mat1->row;

#pragma omp parallel for
  for (int sj = 0; sj < n; sj += BLOCKSIZE) {
    for (int si = 0; si < n; si += BLOCKSIZE) {
      for (int sk = 0; sk < n; sk += BLOCKSIZE) {
        for (int i = si; i < si + BLOCKSIZE; i += UNROLL * 8) {
          for (int j = sj; j < sj + BLOCKSIZE; j++) {
            __m256 c[UNROLL];
            for (int x = 0; x < UNROLL; x++) {
              c[x] = _mm256_load_ps(C + i + x * 8 + j * n);
            }
            for (int k = sk; k < sk + BLOCKSIZE; k++) {
              __m256 b = _mm256_broadcast_ss(B + k + j * n);
              for (int x = 0; x < UNROLL; x++) {
                c[x] = _mm256_add_ps(
                    c[x],
                    _mm256_mul_ps(_mm256_load_ps(A + n * k + x * 8 + i), b));
              }
            }

            for (int x = 0; x < UNROLL; x++) {
              _mm256_store_ps(C + i + x * 8 + j * n, c[x]);
            }
          }
        }
      }
    }
  }
  return true;
#else
  printf("AVX2 is not supported");
  return false;
#endif
}

// 矩阵乘法,Strassen算法,使用AVX2指令集
bool matmul_strassen_avx2(Mat *result, const Mat *mat1, const Mat *mat2) {

  if (isMatNull(mat1) || isMatNull(mat2)) {
    return false;
  }
  if (mat1->col != mat2->row) {
    printf("The size of two matrix is not suitable.");
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

  int row = mat1->row;
  int col = mat1->col;
  int size = mat1->row / 2;

  if (size <= 32) {
    matmul_ikj_avx2(result, mat1, mat2);
    return true;
  }

  Mat *A11 = (Mat *)malloc(sizeof(Mat));
  Mat *A12 = (Mat *)malloc(sizeof(Mat));
  Mat *A21 = (Mat *)malloc(sizeof(Mat));
  Mat *A22 = (Mat *)malloc(sizeof(Mat));

  initMat(A11, size, size);
  initMat(A12, size, size);
  initMat(A21, size, size);
  initMat(A22, size, size);

  Mat *B11 = (Mat *)malloc(sizeof(Mat));
  Mat *B12 = (Mat *)malloc(sizeof(Mat));
  Mat *B21 = (Mat *)malloc(sizeof(Mat));
  Mat *B22 = (Mat *)malloc(sizeof(Mat));

  initMat(B11, size, size);
  initMat(B12, size, size);
  initMat(B21, size, size);
  initMat(B22, size, size);

  Mat *C11 = (Mat *)malloc(sizeof(Mat));
  Mat *C12 = (Mat *)malloc(sizeof(Mat));
  Mat *C21 = (Mat *)malloc(sizeof(Mat));
  Mat *C22 = (Mat *)malloc(sizeof(Mat));

  initMat(C11, size, size);
  initMat(C12, size, size);
  initMat(C21, size, size);
  initMat(C22, size, size);

  Mat *P1 = (Mat *)malloc(sizeof(Mat));
  Mat *P2 = (Mat *)malloc(sizeof(Mat));
  Mat *P3 = (Mat *)malloc(sizeof(Mat));
  Mat *P4 = (Mat *)malloc(sizeof(Mat));
  Mat *P5 = (Mat *)malloc(sizeof(Mat));
  Mat *P6 = (Mat *)malloc(sizeof(Mat));
  Mat *P7 = (Mat *)malloc(sizeof(Mat));

  initMat(P1, size, size);
  initMat(P2, size, size);
  initMat(P3, size, size);
  initMat(P4, size, size);
  initMat(P5, size, size);
  initMat(P6, size, size);
  initMat(P7, size, size);

  Mat *S1 = (Mat *)malloc(sizeof(Mat));
  Mat *S2 = (Mat *)malloc(sizeof(Mat));
  Mat *S3 = (Mat *)malloc(sizeof(Mat));
  Mat *S4 = (Mat *)malloc(sizeof(Mat));
  Mat *S5 = (Mat *)malloc(sizeof(Mat));
  Mat *S6 = (Mat *)malloc(sizeof(Mat));
  Mat *S7 = (Mat *)malloc(sizeof(Mat));
  Mat *S8 = (Mat *)malloc(sizeof(Mat));
  Mat *S9 = (Mat *)malloc(sizeof(Mat));
  Mat *S10 = (Mat *)malloc(sizeof(Mat));

  initMat(S1, size, size);
  initMat(S2, size, size);
  initMat(S3, size, size);
  initMat(S4, size, size);
  initMat(S5, size, size);
  initMat(S6, size, size);
  initMat(S7, size, size);
  initMat(S8, size, size);
  initMat(S9, size, size);
  initMat(S10, size, size);

  Mat *resA = (Mat *)malloc(sizeof(Mat));
  Mat *resB = (Mat *)malloc(sizeof(Mat));

  initMat(resA, size, size);
  initMat(resB, size, size);

  // 对矩阵进行分块
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      A11->data[i * size + j] = mat1->data[i * row + j];
      A12->data[i * size + j] = mat1->data[i * row + j + size];
      A21->data[i * size + j] = mat1->data[(i + size) * row + j];
      A22->data[i * size + j] = mat1->data[(i + size) * row + j + size];

      B11->data[i * size + j] = mat2->data[i * col + j];
      B12->data[i * size + j] = mat2->data[i * col + j + size];
      B21->data[i * size + j] = mat2->data[(i + size) * col + j];
      B22->data[i * size + j] = mat2->data[(i + size) * col + j + size];
    }
  }

  // 计算S1~S10，P1~P7，C11~C22
  float_matSub_avx2(S1, B12, B22);
  float_matAdd_avx2(S2, A11, A12);
  float_matAdd_avx2(S3, A21, A22);
  float_matSub_avx2(S4, B21, B11);
  float_matAdd_avx2(S5, A11, A22);
  float_matAdd_avx2(S6, B11, B22);
  float_matSub_avx2(S7, A12, A22);
  float_matAdd_avx2(S8, B21, B22);
  float_matSub_avx2(S9, A11, A21);
  float_matAdd_avx2(S10, B11, B12);

  matmul_strassen_avx2(P1, A11, S1);
  matmul_strassen_avx2(P2, S2, B22);
  matmul_strassen_avx2(P3, S3, B11);
  matmul_strassen_avx2(P4, A22, S4);
  matmul_strassen_avx2(P5, S5, S6);
  matmul_strassen_avx2(P6, S7, S8);
  matmul_strassen_avx2(P7, S9, S10);

  float_matAdd_avx2(resA, P4, P5);
  float_matSub_avx2(resB, resA, P2);
  float_matAdd_avx2(C11, resB, P6);
  float_matAdd_avx2(C12, P1, P2);
  float_matAdd_avx2(C21, P3, P4);
  float_matAdd_avx2(resA, P5, P1);
  float_matSub_avx2(resB, resA, P3);
  float_matSub_avx2(C22, resB, P7);

  // 将C11~C22合并到C中
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      result->data[i * col + j] = C11->data[i * size + j];
      result->data[i * col + j + size] = C12->data[i * size + j];
      result->data[(i + size) * col + j] = C21->data[i * size + j];
      result->data[(i + size) * col + j + size] = C22->data[i * size + j];
    }
  }

  // 释放内存
  freeMat(A11);
  freeMat(A12);
  freeMat(A21);
  freeMat(A22);
  freeMat(B11);
  freeMat(B12);
  freeMat(B21);
  freeMat(B22);
  freeMat(C11);
  freeMat(C12);
  freeMat(C21);
  freeMat(C22);
  freeMat(P1);
  freeMat(P2);
  freeMat(P3);
  freeMat(P4);
  freeMat(P5);
  freeMat(P6);
  freeMat(P7);
  freeMat(S1);
  freeMat(S2);
  freeMat(S3);
  freeMat(S4);
  freeMat(S5);
  freeMat(S6);
  freeMat(S7);
  freeMat(S8);
  freeMat(S9);
  freeMat(S10);
  freeMat(resA);
  freeMat(resB);

  return true;
}

float *float_vectorAdd_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n) {
#ifdef WITH_AVX2
  if (n % 8 != 0) {
    printf("The size n must be a multiple of 8.");
    return res;
  }

  __m256 a, b;
  size_t THREAD_NUM = 8;
#pragma omp parallel for num_threads(THREAD_NUM) private(a, b)
  for (size_t i = 0; i < n; i += 8) {
    a = _mm256_load_ps(p1 + i);
    b = _mm256_load_ps(p2 + i);
    _mm256_store_ps(res + i, _mm256_add_ps(a, b));
  }
  return res;

#else
  printf("AVX2 is not supported");
  return res;
#endif
}

float *float_vectorSub_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n) {
#ifdef WITH_AVX2
  if (n % 8 != 0) {
    printf("The size n must be a multiple of 8.");
    return res;
  }
  __m256 a, b;
  size_t THREAD_NUM = 8;
#pragma omp parallel for num_threads(THREAD_NUM) private(a, b)
  for (size_t i = 0; i < n; i += 8) {
    a = _mm256_load_ps(p1 + i);
    b = _mm256_load_ps(p2 + i);
    _mm256_store_ps(res + i, _mm256_sub_ps(a, b));
  }
  return res;
#else
  printf("AVX2 is not supported");
  return res;
#endif
}

float *float_vectorMul_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n) {
#ifdef WITH_AVX2
  if (n % 8 != 0) {
    printf("The size n must be a multiple of 8.");
    return res;
  }
  __m256 a, b;
  size_t THREAD_NUM = 8;
#pragma omp parallel for num_threads(THREAD_NUM) private(a, b)
  for (size_t i = 0; i < n; i += 8) {
    a = _mm256_load_ps(p1 + i);
    b = _mm256_load_ps(p2 + i);
    _mm256_store_ps(res + i, _mm256_mul_ps(a, b));
  }
  return res;
#else
  printf("AVX2 is not supported");
  return res;
#endif
}

float *float_vectorDiv_avx2_omp(float *res, const float *p1, const float *p2,
                                size_t n) {
#ifdef WITH_AVX2
  if (n % 8 != 0) {
    printf("The size n must be a multiple of 8.");
    return res;
  }
  __m256 a, b;
  size_t THREAD_NUM = 8;
#pragma omp parallel for num_threads(THREAD_NUM) private(a, b)
  for (size_t i = 0; i < n; i += 8) {
    a = _mm256_load_ps(p1 + i);
    b = _mm256_load_ps(p2 + i);
    _mm256_store_ps(res + i, _mm256_div_ps(a, b));
  }
  return res;
#else
  printf("AVX2 is not supported");
  return res;
#endif
}

bool float_matAdd_avx2_omp(Mat *result, const Mat *mat1, const Mat *mat2) {
#ifdef WITH_AVX2
  if (mat1->row != mat2->row || mat1->col != mat2->col) {
    printf("The size of two matrix must be the same.");
    return false;
  }

  float *p1 = mat1->data;
  float *p2 = mat2->data;
  float *res = result->data;
  size_t n = mat1->row * mat1->col;

  float_vectorAdd_avx2_omp(res, p1, p2, n);
  return true;
#else
  printf("AVX2 is not supported");
  return false;
#endif
}

bool float_matSub_avx2_omp(Mat *result, const Mat *mat1, const Mat *mat2) {
#ifdef WITH_AVX2
  if (mat1->row != mat2->row || mat1->col != mat2->col) {
    printf("The size of two matrix must be the same.");
    return false;
  }

  float *p1 = mat1->data;
  float *p2 = mat2->data;
  float *res = result->data;
  size_t n = mat1->row * mat1->col;

  float_vectorSub_avx2_omp(res, p1, p2, n);
  return true;
#else
  printf("AVX2 is not supported");
  return false;
#endif
}

// 矩阵乘法,Strassen算法,使用AVX2指令集，使用OPENMP并行
bool matmul_strassen_avx2_omp(Mat *result, const Mat *mat1, const Mat *mat2) {
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

  int row = mat1->row;
  int col = mat1->col;
  int size = mat1->row / 2;
  size_t THREAD_NUM = 16;

  if (size <= 64) {
    matmul_ikj_avx2_omp(result, mat1, mat2);
    return true;
  }

  Mat *A11 = (Mat *)malloc(sizeof(Mat));
  Mat *A12 = (Mat *)malloc(sizeof(Mat));
  Mat *A21 = (Mat *)malloc(sizeof(Mat));
  Mat *A22 = (Mat *)malloc(sizeof(Mat));

  initMat(A11, size, size);
  initMat(A12, size, size);
  initMat(A21, size, size);
  initMat(A22, size, size);

  Mat *B11 = (Mat *)malloc(sizeof(Mat));
  Mat *B12 = (Mat *)malloc(sizeof(Mat));
  Mat *B21 = (Mat *)malloc(sizeof(Mat));
  Mat *B22 = (Mat *)malloc(sizeof(Mat));

  initMat(B11, size, size);
  initMat(B12, size, size);
  initMat(B21, size, size);
  initMat(B22, size, size);

  Mat *C11 = (Mat *)malloc(sizeof(Mat));
  Mat *C12 = (Mat *)malloc(sizeof(Mat));
  Mat *C21 = (Mat *)malloc(sizeof(Mat));
  Mat *C22 = (Mat *)malloc(sizeof(Mat));

  initMat(C11, size, size);
  initMat(C12, size, size);
  initMat(C21, size, size);
  initMat(C22, size, size);

  Mat *P1 = (Mat *)malloc(sizeof(Mat));
  Mat *P2 = (Mat *)malloc(sizeof(Mat));
  Mat *P3 = (Mat *)malloc(sizeof(Mat));
  Mat *P4 = (Mat *)malloc(sizeof(Mat));
  Mat *P5 = (Mat *)malloc(sizeof(Mat));
  Mat *P6 = (Mat *)malloc(sizeof(Mat));
  Mat *P7 = (Mat *)malloc(sizeof(Mat));

  initMat(P1, size, size);
  initMat(P2, size, size);
  initMat(P3, size, size);
  initMat(P4, size, size);
  initMat(P5, size, size);
  initMat(P6, size, size);
  initMat(P7, size, size);

  Mat *S1 = (Mat *)malloc(sizeof(Mat));
  Mat *S2 = (Mat *)malloc(sizeof(Mat));
  Mat *S3 = (Mat *)malloc(sizeof(Mat));
  Mat *S4 = (Mat *)malloc(sizeof(Mat));
  Mat *S5 = (Mat *)malloc(sizeof(Mat));
  Mat *S6 = (Mat *)malloc(sizeof(Mat));
  Mat *S7 = (Mat *)malloc(sizeof(Mat));
  Mat *S8 = (Mat *)malloc(sizeof(Mat));
  Mat *S9 = (Mat *)malloc(sizeof(Mat));
  Mat *S10 = (Mat *)malloc(sizeof(Mat));

  initMat(S1, size, size);
  initMat(S2, size, size);
  initMat(S3, size, size);
  initMat(S4, size, size);
  initMat(S5, size, size);
  initMat(S6, size, size);
  initMat(S7, size, size);
  initMat(S8, size, size);
  initMat(S9, size, size);
  initMat(S10, size, size);

  Mat *resA = (Mat *)malloc(sizeof(Mat));
  Mat *resB = (Mat *)malloc(sizeof(Mat));

  initMat(resA, size, size);
  initMat(resB, size, size);

// 对矩阵进行分块
#pragma omp parallel for collapse(2) num_threads(THREAD_NUM)
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      A11->data[i * size + j] = mat1->data[i * row + j];
      A12->data[i * size + j] = mat1->data[i * row + j + size];
      A21->data[i * size + j] = mat1->data[(i + size) * row + j];
      A22->data[i * size + j] = mat1->data[(i + size) * row + j + size];

      B11->data[i * size + j] = mat2->data[i * col + j];
      B12->data[i * size + j] = mat2->data[i * col + j + size];
      B21->data[i * size + j] = mat2->data[(i + size) * col + j];
      B22->data[i * size + j] = mat2->data[(i + size) * col + j + size];
    }
  }

  // 计算S1~S10，P1~P7，C11~C22
  float_matSub_avx2_omp(S1, B12, B22);
  float_matAdd_avx2_omp(S2, A11, A12);
  float_matAdd_avx2_omp(S3, A21, A22);
  float_matSub_avx2_omp(S4, B21, B11);
  float_matAdd_avx2_omp(S5, A11, A22);
  float_matAdd_avx2_omp(S6, B11, B22);
  float_matSub_avx2_omp(S7, A12, A22);
  float_matAdd_avx2_omp(S8, B21, B22);
  float_matSub_avx2_omp(S9, A11, A21);
  float_matAdd_avx2_omp(S10, B11, B12);

  matmul_strassen_avx2_omp(P1, A11, S1);
  matmul_strassen_avx2_omp(P2, S2, B22);
  matmul_strassen_avx2_omp(P3, S3, B11);
  matmul_strassen_avx2_omp(P4, A22, S4);
  matmul_strassen_avx2_omp(P5, S5, S6);
  matmul_strassen_avx2_omp(P6, S7, S8);
  matmul_strassen_avx2_omp(P7, S9, S10);

  float_matAdd_avx2_omp(resA, P4, P5);
  float_matSub_avx2_omp(resB, resA, P2);
  float_matAdd_avx2_omp(C11, resB, P6);
  float_matAdd_avx2_omp(C12, P1, P2);
  float_matAdd_avx2_omp(C21, P3, P4);
  float_matAdd_avx2_omp(resA, P5, P1);
  float_matSub_avx2_omp(resB, resA, P3);
  float_matSub_avx2_omp(C22, resB, P7);

// 将C11~C22合并到C中
#pragma omp parallel for collapse(2) num_threads(THREAD_NUM)
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      result->data[i * col + j] = C11->data[i * size + j];
      result->data[i * col + j + size] = C12->data[i * size + j];
      result->data[(i + size) * col + j] = C21->data[i * size + j];
      result->data[(i + size) * col + j + size] = C22->data[i * size + j];
    }
  }

  // 释放内存
  freeMat(A11);
  freeMat(A12);
  freeMat(A21);
  freeMat(A22);
  freeMat(B11);
  freeMat(B12);
  freeMat(B21);
  freeMat(B22);
  freeMat(C11);
  freeMat(C12);
  freeMat(C21);
  freeMat(C22);
  freeMat(P1);
  freeMat(P2);
  freeMat(P3);
  freeMat(P4);
  freeMat(P5);
  freeMat(P6);
  freeMat(P7);
  freeMat(S1);
  freeMat(S2);
  freeMat(S3);
  freeMat(S4);
  freeMat(S5);
  freeMat(S6);
  freeMat(S7);
  freeMat(S8);
  freeMat(S9);
  freeMat(S10);
  freeMat(resA);
  freeMat(resB);

  return true;
}

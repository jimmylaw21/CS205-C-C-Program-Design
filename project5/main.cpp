#include "Mat.hpp"
#include <cstddef>
#include <iostream>
#include <vector>
#include <sys/time.h>
#include <omp.h>
#pragma GCC optimize(3, "Ofast", "inline")
using namespace std;

int main(int argc, const char **argv) {

  size_t size = 4096;
  Mat<float> mat1(size, size, 3);
  Mat<float> mat2(size, size, 3);


  struct timeval start, end;
  gettimeofday(&start, NULL);
  mat2 = mat1;
  gettimeofday(&end, NULL);
  printf("soft copy Time cost: %ld us\n",
         (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));

  gettimeofday(&start, NULL);
  mat2 = mat1.clone();    
  gettimeofday(&end, NULL);
  printf("hard copy Time cost: %ld us\n",
           (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec));



//   Mat<float> mat3(3, 3, 3);
//   Mat<float> mat4(3, 3, 3);
//   Mat<float> mat5(3, 3, 3);
//   Mat<float> mat6(3, 3, 3);
  
//   for (size_t i = 0; i < mat3.getRows(); i++) {
//       for (size_t j = 0; j < mat3.getCols(); j++) {
//           mat3(i, j, 0) = 1.0;
//           mat3(i, j, 1) = 1.0;
//           mat3(i, j, 2) = 1.0;
//       }
//   }
//   for (size_t i = 0; i < mat4.getRows(); i++) {
//       for (size_t j = 0; j < mat4.getCols(); j++) {
//           mat4(i, j, 0) = 2.0;
//           mat4(i, j, 1) = 2.0;
//           mat4(i, j, 2) = 2.0;
//       }
//   }

//   //验证矩阵通道分离
//     cout << mat3 << endl;
//     vector<Mat<float>> mv;
//     mat3.split(mv);
//     cout << mv[0] << endl << endl;
//     cout << mv[1] << endl << endl;
//     cout << mv[2] << endl << endl;
//     mat3.merge(mv);
//     cout << mat3 << endl;
  
  // for (size_t i = 0; i < mat4.getRows(); i++) {
  //     for (size_t j = 0; j < mat4.getCols(); j++) {
  //         mat4(i, j, 0) = 2.0;
  //         mat4(i, j, 1) = 2.0;
  //         mat4(i, j, 2) = 2.0;
  //     }
  // }
  // // 打印矩阵
  // cout << mat3 << endl;
  // cout << mat4 << endl;
  // // 矩阵相加
  // cout << mat3 + mat4 << endl;
  // // 矩阵相减
  // cout << mat3 - mat4 << endl;
  // // 矩阵相乘
  // cout << mat3 * mat4 << endl;
  // // 矩阵相除
  // // cout << mat3 / mat4 << endl;

  // Mat <float> mat5(3, 1, 3);
  // Mat <float> mat6(1, 3, 3);
  // for (size_t i = 0; i < mat5.getRows(); i++) {
  //     for (size_t j = 0; j < mat5.getCols(); j++) {
  //         mat5(i, j, 0) = 1.0;
  //         mat5(i, j, 1) = 1.0;
  //         mat5(i, j, 2) = 1.0;
  //     }
  // }
  // for (size_t i = 0; i < mat6.getRows(); i++) {
  //     for (size_t j = 0; j < mat6.getCols(); j++) {
  //         mat6(i, j, 0) = 2.0;
  //         mat6(i, j, 1) = 2.0;
  //         mat6(i, j, 2) = 2.0;
  //     }
  // }
  // cout << mat5 << endl;
  // cout << mat6 << endl;
  // cout << mat5 * mat6 << endl;

  // mat4 = mat3.ROI(1, 1, 1, 1);
  

//   mat1.loadFromBMP("Lena3");
//   mat1.rgbtogray().saveAsBMP("mat1");
// //   mat1.ROI(mat2, 128, 128, 256, 256);
//   mat2 = mat1;
//   mat2.saveAsBMP("mat2");

  // Mat<int> mat1(3, 3, 1);
  // for(size_t i = 0; i < mat1.getRows(); i++) {
  //   for(size_t j = 0; j < mat1.getCols(); j++) {
  //     mat1(i, j, 0) = i + j;
  //   }
  // }
  // cout << mat1 << endl;
  // cout <<"mat1.getRefCount() = " << mat1.getRefCount() << endl;
  // cout <<"mat1.getNumMatrices() = " << mat1.getNumMatrices() << endl;
  // Mat<int> mat2(3, 3, 1);
  // mat2 = mat1;
  // cout << mat2 << endl;
  // cout <<"mat2.getRefCount() = "  << mat2.getRefCount() << endl;
  // cout <<"mat2.getNumMatrices() = "  << mat2.getNumMatrices() << endl;
  // Mat<int> mat3(mat1);
  // cout << mat3 << endl;
  // cout <<"mat3.getRefCount() = "  << mat3.getRefCount() << endl;
  // cout <<"mat3.getNumMatrices() = "  << mat3.getNumMatrices() << endl;
  // Mat<int> mat4(3, 3, 1);
  // cout << mat4 << endl;
  // cout <<"mat4.getRefCount() = "  << mat4.getRefCount() << endl;
  // cout <<"mat4.getNumMatrices() = "  << mat4.getNumMatrices() << endl;
  return 0;
}
// g++ -o main main.cpp Mat.hpp && ./main
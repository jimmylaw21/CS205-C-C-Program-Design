#include <cstddef>
#include<iostream>
#include<cstring>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include <omp.h>
#include <immintrin.h>
#pragma GCC optimize(3, "Ofast", "inline")

#define UNROLL 4
#define BLOCKSIZE 32

using namespace std;

//cv::Mat
template <class T>
class Mat {
private:
    size_t rows;
    size_t cols;
    size_t channels;
    size_t span;
    // dataStart指向的是矩阵的第一个元素，
    // dataEnd指向的是矩阵的最后一个元素的下一个元素, 
    // data指向的是ROI矩阵的第一个元素, 
    // data+ rows * channels * cols指向的是ROI矩阵的最后一个元素的下一个元素
    T* data;
    T* dataEnd;
    T* dataStart;
    static size_t num_matrices;
    size_t* ref_count;

public:
    // constructors
    Mat(size_t rows, size_t cols, size_t channels);
    Mat(size_t rows, size_t cols, size_t channels, size_t span, T* data = nullptr, T* dataEnd = nullptr, T* dataStart = nullptr, size_t* ref_count = nullptr);
    Mat(const Mat& other);
    //根据传入文件流来创建矩阵的构造函数
    Mat(size_t row, size_t col, size_t channel, ifstream& file);
    ~Mat();

    // accessors
    size_t getRows() const;
    size_t getCols() const;
    size_t getChannels() const;
    size_t getSpan() const;
    T* getData() const;
    T* getDataEnd() const;
    T* getDataStart() const;
    size_t getRefCount() const;
    static size_t getNumMatrices();

    // assignment operators
    Mat& operator=(const Mat& other);

    // array operators
    T& operator()(size_t row, size_t col, size_t channel) const;
    T& operator()(size_t row, size_t col) const;
    T& operator[](size_t index) const;

    //arithmetic functions
    // mat transpose
    Mat& transpose(Mat& dst) const;
    Mat& transpose();
    // mat inverse
    Mat& inv(Mat& dst) const;
    Mat& inv();
    //matdet
    T det() const;

    // Mat& add(Mat& dst, const Mat& other) const;
    // Mat& sub(Mat& dst, const Mat& other) const;
    // Mat& mul(Mat& dst, const Mat& other) const;
    // Mat& dot(Mat& dst, const Mat& other) const;
    // Mat& div(Mat& dst, const Mat& other) const;

    // Mat& add(Mat& dst, const T& other) const;
    // Mat& sub(Mat& dst, const T& other) const;
    // Mat& dot(Mat& dst, const T& other) const;
    // Mat& div(Mat& dst, const T& other) const;

    // Mat& add(const Mat& other);
    // Mat& sub(const Mat& other);
    // Mat& mul(const Mat& other);
    // Mat& dot(const Mat& other);
    // Mat& div(const Mat& other);

    // Mat& add(const T& other);
    // Mat& sub(const T& other);
    // Mat& mul(const T& other);
    // Mat& dot(const T& other);
    // Mat& div(const T& other);

    //arithmetic operaters
    Mat& operator+=(const Mat& other);
    Mat& operator-=(const Mat& other);
    Mat& operator*=(const Mat& other);//dot
    Mat& operator/=(const Mat& other);//dot

    Mat& operator+=(const T& value);
    Mat& operator-=(const T& value);
    Mat& operator*=(const T& value);
    Mat& operator/=(const T& value);

    Mat operator+(const Mat& other) const;
    Mat operator-(const Mat& other) const;
    Mat operator*(const Mat& other) const;
    Mat operator/(const Mat& other) const;

    Mat operator+(const T& value) const;
    Mat operator-(const T& value) const;
    Mat operator*(const T& value) const;
    Mat operator/(const T& value) const;

    // unary operaters
    Mat operator!() const;
    Mat operator~() const;
    Mat& operator++();
    Mat& operator--();
    Mat operator++(int);
    Mat operator--(int);

    // friend operaters
    template<class TT>
    friend Mat operator+(TT num, const Mat<TT>& mat){
        return mat + num;
    }
    template<class TT>
    friend Mat operator-(TT num, const Mat<TT>& mat){
        return mat - num;
    }
    template<class TT>
    friend Mat operator*(TT num, const Mat<TT>& mat){
        return mat * num;
    }
    template<class TT>
    friend Mat operator/(TT num, const Mat<TT>& mat){
        return mat / num;
    }
    template<class TT>
    friend ostream& operator<<(ostream& os, const Mat<TT>& mat){
        size_t r = mat.rows;
        size_t c = mat.cols;
        size_t ch = mat.channels;
        for(size_t i = 0; i < r; i++){
            for(size_t j = 0; j < c; j++){
                os << "{";
                for(size_t k = 0; k < ch; k++){
                    os << mat(i, j, k);
                    if(k != ch - 1){
                        os << " ";
                    };
                }
                os << "}";
                if(j != c - 1){
                    os << ", ";
                }
            }
            if(i != r - 1){
                os << "\n";
            } 
        }
        return os;
    }

    template<class TT>
    friend ostream& operator>>(ostream& os, const Mat<TT>& mat){
        for(size_t i = 0; i < mat.rows; i++){
            for(size_t j = 0; j < mat.cols; j++){
                for(size_t k = 0; k < mat.channels; k++){
                    os >> mat.data[i * mat.span + j * mat.channels + k];
                }
            }
        }
        return os;
    }

    // boolean operaters
    bool operator==(const Mat& other) const;
    bool operator!=(const Mat& other) const;

    // functions
    size_t size(){
        return rows * cols;
    }
    bool empty(){
        return rows == 0 || cols == 0 || channels == 0;
    }

    //ROI
    Mat& ROI(Mat& dst, size_t row, size_t col, size_t height, size_t width) const;
    Mat ROI(size_t row, size_t col, size_t height, size_t width) const;     
              

    //IO
    void save(const string& filename) const;
    void load(const string& filename);
    void saveAsBMP(const std::string& filename) const;
    void loadFromBMP(const std::string& filename);

    //Mat functions
    Mat<T> clone() const;
    Mat<T> copyTo(Mat<T>& dst) const;
    Mat<T> reshape(size_t new_cn, size_t new_rows, size_t new_cols) const;
    Mat& convertTo(Mat& dst, int rtype, double alpha = 1, double beta = 0) const;
    static Mat& zeros(size_t rows, size_t cols, size_t channels);
    static Mat& ones(size_t rows, size_t cols, size_t channels);
    static Mat& eye(size_t rows, size_t cols, size_t channels);
    Mat& append(T *data, int row, int col, int ch);
    Mat& subtract(int row);
    void split(vector<Mat>& mv) const;
    Mat& merge(const vector<Mat>& mv);
    Mat merge_vertical(const Mat& mat) const;
    Mat merge_horizontal(const Mat& mat) const;

    Mat rgbtogray() const;
};
// 静态成员变量初始化
template <class T>
size_t Mat<T>::num_matrices = 0;

// 默认构造函数1
template <typename T>
Mat<T>::Mat(size_t rows, size_t cols, size_t channels){
    //参数检查合法性
    if(rows == 0 || cols == 0 || channels == 0){
        std::cerr << "In default constructor1" << std::endl;
        throw invalid_argument("Mat::Mat: rows, cols and channels must be greater than 0");
    }

    //初始化
    this->rows = rows;
    this->cols = cols;
    this->channels = channels;
    this->data = new T[rows * cols * channels];
    //data数组全置0
    memset(this->data, 0, rows * cols * channels * sizeof(T));
    this->dataEnd = data + rows * cols * channels;
    this->dataStart = data;
    this->span = cols * channels;
    this->ref_count = new size_t(1);

    num_matrices++;
}
// 默认构造函数2
template <typename T>
Mat<T>::Mat(size_t rows, size_t cols, size_t channels, size_t span, T* data, T* dataEnd, T* dataStart, size_t* ref_count){
    //参数检查合法性
    if(rows == 0 || cols == 0 || channels == 0){
        std::cerr << "In default constructor2" << std::endl;
        throw invalid_argument("Mat::Mat: rows, cols and channels must be greater than 0");
    }
    //对ROI的边界的合法性进行检查
    if(data < dataStart || dataEnd < data + rows * cols * channels || dataEnd > dataStart + rows * cols * channels){
        std::cerr << "In default constructor2" << std::endl;
        throw invalid_argument("The region of interest should not exceed the region of the matrix.");
    }
        
    if (ref_count == nullptr){
        this->ref_count = new size_t(1);
    }else{
        this->ref_count = ref_count;
        // (*ref_count)++;
    }

    this->rows = rows;
    this->cols = cols;
    this->channels = channels;
    this->dataEnd = dataEnd;
    this->dataStart = dataStart;

    if (span == 0){
        this->span = cols * channels;
    }else{
        this->span = span;
    }

    this->data = new T[rows * cols * channels];
    if (data != nullptr){
        memcpy(this->data, data, rows * cols * channels * sizeof(T));
    }
    num_matrices++;
}
// 复制构造函数
template <typename T>
Mat<T>::Mat(const Mat& other) : rows(other.rows), cols(other.cols), channels(other.channels), span(other.span), data(other.data),dataEnd(other.dataEnd), dataStart(other.dataStart), ref_count(other.ref_count){
    num_matrices++;
    (*ref_count)++;
}
// 文件构造函数
template <typename T>
Mat<T>::Mat(size_t row, size_t col, size_t channel, ifstream& file){
    //参数检查合法性
    if(row == 0 || col == 0 || channel == 0){
        std::cerr << "In file constructor" << std::endl;
        throw invalid_argument("Mat::Mat: rows, cols and channels must be greater than 0");
    }

    //初始化
    this->rows = row;
    this->cols = col;
    this->channels = channel;
    this->data = new T[row * col * channel];
    this->dataEnd = data + row * col * channel;
    this->dataStart = data;
    this->span = col * channel;
    this->ref_count = new size_t(1);

    //读取文件
    for (size_t i = 0; i < row * col * channel; i++){
        file >> data[i];
    }
    num_matrices++;
}
// 析构函数
template <typename T>
Mat<T>::~Mat(){
    if(--(*ref_count) == 0){
        delete[] data;
        this -> data = nullptr;
        delete[] ref_count;
        this -> ref_count = nullptr;
        num_matrices--;
        // cout << "The matrix is free" << endl;
    }
    // cout << "Destructor called" << endl;
}

// accessors
template <typename T>
size_t Mat<T>::getRows() const{
    return rows;
}
template <typename T>
size_t Mat<T>::getCols() const{
    return cols;
}
template <typename T>
size_t Mat<T>::getChannels() const{
    return channels;
}
template <typename T>
size_t Mat<T>::getSpan() const{
    return span;
}
template <typename T>
T* Mat<T>::getData() const{
    return data;
}
template <typename T>
T* Mat<T>::getDataEnd() const{
    return dataEnd;
}
template <typename T>
T* Mat<T>::getDataStart() const{
    return dataStart;
}
template <typename T>
size_t Mat<T>::getRefCount() const{
    return *ref_count;
}
template <typename T>
size_t Mat<T>::getNumMatrices() {
    return num_matrices;
}

// assignment operator
template <typename T>
Mat<T>& Mat<T>::operator=(const Mat<T>& other){
    if (this != &other){
        if(--(*ref_count) == 0){
            delete[] data;
            this -> data = nullptr;
            delete[] ref_count;
            this -> ref_count = nullptr;
            num_matrices--;
        }
        this->rows = other.rows;
        this->cols = other.cols;
        this->channels = other.channels;
        this->data = other.data;
        this->dataEnd = other.dataEnd;
        this->dataStart = other.dataStart;
        this->span = other.span;
        this->ref_count = other.ref_count;
        (*ref_count)++;
        num_matrices++;
    }
    return *this;
}

// array operaters
// T& operator()(size_t row, size_t col, size_t channel) const;
template <typename T>
T& Mat<T>::operator()(size_t row, size_t col, size_t channel) const{
    //参数检查合法性
    if(row >= rows || col >= cols || channel >= channels){
        std::cerr << "In ()" << std::endl;
        throw invalid_argument("The index is out of range.");
    }
    return data[row * span + col * channels + channel];
}
// T& operator()(size_t row, size_t col) const;
template <typename T>
T& Mat<T>::operator()(size_t row, size_t col) const{
    //参数检查合法性
    if(row >= rows || col >= cols){
        std::cerr << "In ()" << std::endl;
        throw invalid_argument("The index is out of range.");
    }
    return data[row * span + col * channels];
}
// T& operator[](size_t index) const;
template <typename T>
T& Mat<T>::operator[](size_t index) const{
    //参数检查合法性
    if(index >= rows * cols * channels){
        std::cerr << "In []" << std::endl;
        throw invalid_argument("The index is out of range.");
    }
    return data[index];
}

//transpose
//Mat& transpose();
template <typename T>
Mat<T>& Mat<T>::transpose(){
    //判断参数合法性
    if(rows != cols){
        std::cerr << "In transpose function..." << std::endl;
        std::cerr << "The matrix should be square." << std::endl;
        exit(EXIT_FAILURE);
    }
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        T* temp = new T[size];
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++){
                for (size_t ch = 0; ch < channels; ch++)
                {
                    temp[j * rows * channels + i * channels + ch] = (*this)(i, j, ch);
                }
            }
        }
        //复制temp到data
        memcpy(data, temp, size * sizeof(T));
        delete[] temp;
    }else{
        T* temp = new T[size];
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++){
                for (size_t ch = 0; ch < channels; ch++)
                {
                    temp[i * cols * channels + j * channels + ch] = (*this)(i, j, ch);
                }
            }
        }
        delete[] data;
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

//Mat& transpose(Mat& dst) const;
template <typename T>
Mat<T>& Mat<T>::transpose(Mat<T>& dst) const{
    //判断参数合法性
    if(rows != cols){
        std::cerr << "In transpose function..." << std::endl;
        std::cerr << "The matrix should be square." << std::endl;
        exit(EXIT_FAILURE);
    }
    size_t size = rows * cols * channels;
    T* temp = new T[size];
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            for (size_t ch = 0; ch < channels; ch++)
            {
                temp[j * rows * channels + i * channels + ch] = (*this)(i, j, ch);
            }
        }
    }
    dst = Mat<T>(rows, cols, channels, temp);
    return dst;
}

//det
//T det();
template <typename T>
T Mat<T>::det() const{
    //判断参数合法性
    if(rows != cols){
        std::cerr << "In det function..." << std::endl;
        std::cerr << "The matrix should be square." << std::endl;
        exit(EXIT_FAILURE);
    }

    T result = 0;
    if(rows == 1){
        result = (*this)(0, 0, 0);
    }else if(rows == 2){
        result = (*this)(0, 0, 0) * (*this)(1, 1, 0) - (*this)(0, 1, 0) * (*this)(1, 0, 0);
    }else{
        for(int i = 0; i < rows; i ++){
            Mat<T> temp(rows - 1, cols - 1, channels);
            temp.dataEnd = temp.data + (rows - 1) * (cols - 1) * channels;
            temp.dataStart = temp.data;
            for(int j = 0; j < rows - 1; j ++){
                for(int k = 0; k < cols - 1; k ++){
                    for(int ch = 0; ch < channels; ch ++){
                        if(k < i){
                            temp(j, k, ch) = (*this)(j + 1, k, ch);
                        }else{
                            temp(j, k, ch) = (*this)(j + 1, k + 1, ch);
                        }
                    }
                }
            }
            result += pow(-1, i) * (*this)(0, i, 0) * temp.det();
        }
    }

    return result;
}

//inverse
//Mat& inv();
template <typename T>
Mat<T>& Mat<T>::inv(){
    //判断参数合法性
    if(rows != cols){
        std::cerr << "In inv function..." << std::endl;
        std::cerr << "The matrix should be square." << std::endl;
        exit(EXIT_FAILURE);
    }

    T det = this->det();
    if(det == 0){
        std::cerr << "In inv function..." << std::endl;
        std::cerr << "The matrix is not invertible." << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        T* temp = new T[size];
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++){
                Mat<T> temp(rows - 1, cols - 1, channels);
                temp.dataEnd = temp.data + (rows - 1) * (cols - 1) * channels;
                temp.dataStart = temp.data;
                for(int k = 0; k < rows - 1; k ++){
                    for(int l = 0; l < cols - 1; l ++){
                        for(int ch = 0; ch < channels; ch ++){
                            if(k < i){
                                if(l < j){
                                    temp(k, l, ch) = (*this)(k, l, ch);
                                }else{
                                    temp(k, l, ch) = (*this)(k, l + 1, ch);
                                }
                            }else{
                                if(l < j){
                                    temp(k, l, ch) = (*this)(k + 1, l, ch);
                                }else{
                                    temp(k, l, ch) = (*this)(k + 1, l + 1, ch);
                                }
                            }
                        }
                    }
                }
                temp(j, i, 0) = pow(-1, i + j) * temp.det();
            }
        }
        //复制temp到data
        memcpy(data, temp, size * sizeof(T));
        delete[] temp;
    }else{
        T* temp = new T[size];
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < cols; j++){
                Mat<T> temp(rows - 1, cols - 1, channels);
                temp.dataEnd = temp.data + (rows - 1) * (cols - 1) * channels;
                temp.dataStart = temp.data;
                for(int k = 0; k < rows - 1; k ++){
                    for(int l = 0; l < cols - 1; l ++){
                        for(int ch = 0; ch < channels; ch ++){
                            if(k < i){
                                if(l < j){
                                    temp(k, l, ch) = (*this)(k, l, ch);
                                }else{
                                    temp(k, l, ch) = (*this)(k, l + 1, ch);
                                }
                            }else{
                                if(l < j){
                                    temp(k, l, ch) = (*this)(k + 1, l, ch);
                                }else{
                                    temp(k, l, ch) = (*this)(k + 1, l + 1, ch);
                                }
                            }
                        }
                    }
                }
                temp(j, i, 0) = pow(-1, i + j) * temp.det();
            }
        }
        delete[] data;
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

//Mat& inv(Mat& dst) const;
template <typename T>
Mat<T>& Mat<T>::inv(Mat<T>& dst) const{
    //判断参数合法性
    if(rows != cols){
        std::cerr << "In inv function..." << std::endl;
        std::cerr << "The matrix should be square." << std::endl;
        exit(EXIT_FAILURE);
    }

    T det = this->det();
    if(det == 0){
        std::cerr << "In inv function..." << std::endl;
        std::cerr << "The matrix is not invertible." << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t size = rows * cols * channels;
    T* temp = new T[size];
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            Mat<T> temp(rows - 1, cols - 1, channels);
            temp.dataEnd = temp.data + (rows - 1) * (cols - 1) * channels;
            temp.dataStart = temp.data;
            for(int k = 0; k < rows - 1; k ++){
                for(int l = 0; l < cols - 1; l ++){
                    for(int ch = 0; ch < channels; ch ++){
                        if(k < i){
                            if(l < j){
                                temp(k, l, ch) = (*this)(k, l, ch);
                            }else{
                                temp(k, l, ch) = (*this)(k, l + 1, ch);
                            }
                        }else{
                            if(l < j){
                                temp(k, l, ch) = (*this)(k + 1, l, ch);
                            }else{
                                temp(k, l, ch) = (*this)(k + 1, l + 1, ch);
                            }
                        }
                    }
                }
            }
            temp(j, i, 0) = pow(-1, i + j) * temp.det();
        }
    }
    //复制temp到data
    memcpy(dst.data, temp, size * sizeof(T));
    delete[] temp;
    return dst;
}

//arithmetic operaters

template <typename T>
Mat<T>& Mat<T>::operator+=(const Mat<T>& other){
    //参数检查合法性
    if(rows != other.rows || cols != other.cols || channels != other.channels){
        std::cerr << "In +=" << std::endl;
        throw invalid_argument("The size of two matrices must be the same.");
    }
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] += other.data[i];
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] + other.data[i];
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator-=(const Mat<T>& other){
    //参数检查合法性
    if(rows != other.rows || cols != other.cols || channels != other.channels){
        std::cerr << "In -=" << std::endl;
        throw invalid_argument("The size of two matrices must be the same.");
    }
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] -= other.data[i];
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] - other.data[i];
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator*=(const Mat<T>& other){
    //参数检查合法性
    if(rows != other.rows || cols != other.cols || channels != other.channels){
        std::cerr << "In *=" << std::endl;
        throw invalid_argument("The size of two matrices must be the same.");
    }
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] *= other.data[i];
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] * other.data[i];
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator/=(const Mat<T>& other){
    //参数检查合法性
    if(rows != other.rows || cols != other.cols || channels != other.channels){
        std::cerr << "In /=" << std::endl;
        throw invalid_argument("The size of two matrices must be the same.");
    }
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] /= other.data[i];
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] / other.data[i];
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator+=(const T& value){
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] += value;
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] + value;
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator-=(const T& value){
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] -= value;
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] - value;
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator*=(const T& value){
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] *= value;
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] * value;
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T>& Mat<T>::operator/=(const T& value){
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i] /= value;
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] / value;
        }
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

template <typename T>
Mat<T> Mat<T>::operator+(const Mat<T>& other) const{
    Mat<T> result(*this);
    result += other;
    return result;
}

template <typename T>
Mat<T> Mat<T>::operator-(const Mat<T>& other) const{
    Mat<T> result(*this);
    result -= other;
    return result;
}

template <typename T>
Mat<T> Mat<T>::operator*(const Mat<T>& other) const{
    //判断参数合法性
    if(other.rows != this -> cols || other.channels != channels){
        std::cerr << "In operator * friend functoin..." << std::endl;
        std::cerr << "The two matrices' sizes and channel number should be the same." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    Mat<T> result(rows, other.cols, channels);
    result.dataEnd = result.data + rows * other.cols * channels;
    result.dataStart = result.data;

    //矩阵乘法
    #pragma omp parallel for
    for (size_t k = 0; k < other.rows; k++){
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < other.cols; j++){
                for (size_t ch = 0; ch < channels; ch++){
                    result(i, j, ch) += (*this)(i, k, ch) * other(k, j, ch);
                }
            }
        }
    }
    num_matrices++;

    return result;
}

template <typename T>
Mat<T> Mat<T>::operator/(const Mat<T>& other) const{
    //判断参数合法性
    if(other.rows != this -> cols || other.channels != channels){
        std::cerr << "In operator * friend functoin..." << std::endl;
        std::cerr << "The two matrices' sizes and channel number should be the same." << std::endl;
        exit(EXIT_FAILURE);
    }

    Mat<T> result(rows, other.cols, channels);
    result.dataEnd = result.data + rows * other.cols * channels;
    result.dataStart = result.data;

    //求逆矩阵
    Mat<T> other_inv(other.rows, other.cols, other.channels);
    this -> inv(other_inv);
    #pragma omp parallel for
    for (size_t k = 0; k < other_inv.rows; k++){
        for (size_t i = 0; i < rows; i++){
            for (size_t j = 0; j < other_inv.cols; j++){
                for (size_t ch = 0; ch < channels; ch++){
                    result(i, j, ch) += (*this)(i, k, ch) * other_inv(k, j, ch);
                }
            }
        }
    }

    return result;
}

template <typename T>
Mat<T> Mat<T>::operator+(const T& value) const{
    Mat<T> result(*this);
    result += value;
    return result;
}

template <typename T>
Mat<T> Mat<T>::operator-(const T& value) const{
    Mat<T> result(*this);
    result -= value;
    return result;
}

template <typename T>
Mat<T> Mat<T>::operator*(const T& value) const{
    Mat<T> result(*this);
    result *= value;
    return result;
}

template <typename T>
Mat<T> Mat<T>::operator/(const T& value) const{
    Mat<T> result(*this);
    result /= value;
    return result;
}

//Mat operator!() const;
template <typename T>
Mat<T> Mat<T>::operator!() const{
    //判断参数合法性
    if(rows != cols){
        std::cerr << "In operator! function..." << std::endl;
        std::cerr << "The matrix should be square." << std::endl;
        exit(EXIT_FAILURE);
    }

    Mat<T> result(rows, cols, channels);
    result.dataEnd = result.data + rows * cols * channels;
    result.dataStart = result.data;

    //转置
    this -> transpose(result);

    return result;
}

//Mat operator~() const;
template <typename T>
Mat<T> Mat<T>::operator~() const{
    //判断参数合法性
    if(rows != cols){
        std::cerr << "In operator~ function..." << std::endl;
        std::cerr << "The matrix should be square." << std::endl;
        exit(EXIT_FAILURE);
    }

    Mat<T> result(rows, cols, channels);
    result.dataEnd = result.data + rows * cols * channels;
    result.dataStart = result.data;

    //求逆矩阵
    this -> inv(result);

    return result;
}

// Mat& operator++();
template <typename T>
Mat<T>& Mat<T>::operator++(){
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i]++;
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] + 1;
        }
        delete[] data;
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

// Mat& operator--();
template <typename T>
Mat<T>& Mat<T>::operator--(){
    size_t size = rows * cols * channels;
    if (*ref_count == 1){
        for (size_t i = 0; i < size; i++){
            data[i]--;
        }
    }
    else{
        T* temp = new T[size];
        for (size_t i = 0; i < size; i++){
            temp[i] = data[i] - 1;
        }
        delete[] data;
        data = temp;
        dataEnd = data + size;
        dataStart = data;
        (*ref_count)--;
        ref_count = new size_t(1);
    }
    return *this;
}

// Mat operator++(int);
template <typename T>
Mat<T> Mat<T>::operator++(int){
    Mat<T> result(*this);
    ++(*this);
    return result;
}

// Mat operator--(int);
template <typename T>
Mat<T> Mat<T>::operator--(int){
    Mat<T> result(*this);
    --(*this);
    return result;
}

// bool operator==(const Mat& other) const;
template <typename T>
bool Mat<T>::operator==(const Mat<T>& other) const{
    if (rows != other.rows || cols != other.cols || channels != other.channels){
        return false;
    }
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        if (data[i] != other.data[i]){
            return false;
        }
    }
    return true;
}

// bool operator!=(const Mat& other) const;
template <typename T>
bool Mat<T>::operator!=(const Mat<T>& other) const{
    return !(*this == other);
}

//ROI
//Mat& ROI(Mat& dst, size_t row, size_t col, size_t height, size_t width) const;
//ROI区域必须在原图像范围内，与原图像共享内存，dataStart和dataEnd不变，data指向ROI区域的起始位置，即row行col列的位置，height和width变为ROI区域的行数和列数
template <typename T>
Mat<T>& Mat<T>::ROI(Mat<T>& dst, size_t row, size_t col, size_t height, size_t width) const{
    //判断参数合法性
    if(row + height > rows || col + width > cols){
        std::cerr << "In ROI function..." << std::endl;
        std::cerr << "The ROI region should be in the original image." << std::endl;
        exit(EXIT_FAILURE);
    }

    dst.rows = height;
    dst.cols = width;
    dst.channels = channels;
    dst.span = cols * channels;
    //free
    if (dst.ref_count != NULL){
        (*dst.ref_count)--;
        if (*dst.ref_count == 0){
            delete[] dst.data;
            delete dst.ref_count;
        }
    }
    dst.data = data + row * cols * channels + col * channels;
    // cout << "dst.data = " << dst.data << endl;
    // cout << "data = " << data << endl;
    dst.dataStart = dataStart;
    dst.dataEnd = dataEnd;
    dst.ref_count = ref_count;
    (*ref_count)++;

    return dst;
}
//Mat ROI(size_t row, size_t col, size_t height, size_t width) const;
template <typename T>
Mat<T> Mat<T>::ROI(size_t row, size_t col, size_t height, size_t width) const{
    Mat<T> result(rows, cols, channels);
    num_matrices++;
    ROI(result, row, col, height, width);
    return result;
}

//IO
//void save(const std::string& filename) const;
template <typename T>
void Mat<T>::save(const std::string& filename) const{
    std::ofstream ofs(filename  + ".txt");
    if (!ofs){
        std::cerr << "In save function..." << std::endl;
        std::cerr << "Can't open the file." << std::endl;
        exit(EXIT_FAILURE);
    }
    ofs << rows << " " << cols << " " << channels << std::endl;
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        ofs << data[i] << " ";
    }
    ofs.close();
}

//void load(const std::string& filename);
template <typename T>
void Mat<T>::load(const std::string& filename){
    std::ifstream ifs(filename + ".txt");
    if (!ifs){
        std::cerr << "In load function..." << std::endl;
        std::cerr << "Can't open the file." << std::endl;
        exit(EXIT_FAILURE);
    }
    ifs >> rows >> cols >> channels;
    size_t size = rows * cols * channels;
    data = new T[size];
    dataStart = data;
    dataEnd = data + size;
    ref_count = new size_t(1);
    for (size_t i = 0; i < size; i++){
        ifs >> data[i];
    }
    ifs.close();
}

//将三通道矩阵转换为RGB图像，保存为bmp格式，每个通道的值必须在0~255之间，否则会出错，channels必须为3，第一位为R，第二位为G，第三位为B
//void saveAsBMP(const std::string& filename) const;
template <typename T>
void Mat<T>::saveAsBMP(const std::string& filename) const{
    if (channels != 3){
        std::cerr << "In saveAsBMP function..." << std::endl;
        std::cerr << "The channels of the image should be 3." << std::endl;
        exit(EXIT_FAILURE);
    }
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        if (data[i] < -128 || data[i] > 127){
            std::cerr << "In saveAsBMP function..." << std::endl;
            std::cerr << "The value of each channel should be in [0, 255]." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    std::ofstream ofs(filename + ".bmp", std::ios::binary);
    if (!ofs){
        std::cerr << "In saveAsBMP function..." << std::endl;
        std::cerr << "Can't open the file." << std::endl;
        exit(EXIT_FAILURE);
    }
    //文件头
    char fileHeader[14] = {0};
    fileHeader[0] = 'B';
    fileHeader[1] = 'M';
    size_t fileSize = 14 + 40 + rows * cols * 3;
    fileHeader[2] = fileSize & 0xff;
    fileHeader[3] = (fileSize >> 8) & 0xff;
    fileHeader[4] = (fileSize >> 16) & 0xff;
    fileHeader[5] = (fileSize >> 24) & 0xff;
    fileHeader[10] = 54;
    ofs.write(fileHeader, 14);
    //信息头
    char infoHeader[40] = {0};
    infoHeader[0] = 40;
    infoHeader[4] = cols & 0xff;
    infoHeader[5] = (cols >> 8) & 0xff;
    infoHeader[6] = (cols >> 16) & 0xff;
    infoHeader[7] = (cols >> 24) & 0xff;
    infoHeader[8] = rows & 0xff;
    infoHeader[9] = (rows >> 8) & 0xff;
    infoHeader[10] = (rows >> 16) & 0xff;
    infoHeader[11] = (rows >> 24) & 0xff;
    infoHeader[12] = 1;
    infoHeader[14] = 24;
    ofs.write(infoHeader, 40);
    //数据
    char* data = new char[rows * cols * 3];
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            data[(i * cols + j) * 3] = this->data[(i * cols + j) * 3];
            data[(i * cols + j) * 3 + 1] = this->data[(i * cols + j) * 3 + 1];
            data[(i * cols + j) * 3 + 2] = this->data[(i * cols + j) * 3 + 2];
        }
    }
    ofs.write(data, rows * cols * 3);
    ofs.close();
    delete[] data;
}

//将RGB图像转换为三通道矩阵，channels必须为3，第一位为R，第二位为G，第三位为B
//void loadFromBMP(const std::string& filename);
template <typename T>
void Mat<T>::loadFromBMP(const std::string& filename){
    std::ifstream ifs(filename + ".bmp", std::ios::binary);
    if (!ifs){
        std::cerr << "In loadFromBMP function..." << std::endl;
        std::cerr << "Can't open the file." << std::endl;
        exit(EXIT_FAILURE);
    }
    //文件头
    char fileHeader[14] = {0};
    ifs.read(fileHeader, 14);
    if (fileHeader[0] != 'B' || fileHeader[1] != 'M'){
        std::cerr << "In loadFromBMP function..." << std::endl;
        std::cerr << "The file is not a bmp file." << std::endl;
        exit(EXIT_FAILURE);
    }
    //信息头
    char infoHeader[40] = {0};
    ifs.read(infoHeader, 40);
    if (infoHeader[0] != 40){
        std::cerr << "In loadFromBMP function..." << std::endl;
        std::cerr << "The file is not a bmp file." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (infoHeader[14] != 24){
        std::cerr << "In loadFromBMP function..." << std::endl;
        std::cerr << "The file is not a 24-bit bmp file." << std::endl;
        exit(EXIT_FAILURE);
    }
    rows = *(int*)(infoHeader + 8);
    cols = *(int*)(infoHeader + 4);
    channels = 3;
    size_t size = rows * cols * channels;
    data = new T[size];
    dataStart = data;
    dataEnd = data + size;
    ref_count = new size_t(1);
    //数据
    char* data = new char[rows * cols * 3];
    ifs.read(data, rows * cols * 3);
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            this->data[(i * cols + j) * 3] = data[(i * cols + j) * 3];
            this->data[(i * cols + j) * 3 + 1] = data[(i * cols + j) * 3 + 1];
            this->data[(i * cols + j) * 3 + 2] = data[(i * cols + j) * 3 + 2];
        }
    }
    ifs.close();
    delete[] data;
}

//Mat functions
//Mat<T> clone() const;
template <typename T>
Mat<T> Mat<T>::clone() const{
    Mat<T> mat(rows, cols, channels);
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        mat.data[i] = data[i];
    }
    return mat;
}

//Mat<T> reshape(size_t rows, size_t cols) const;
// Mat<T> copyTo(Mat<T>& dst) const;
template <typename T>
Mat<T> Mat<T>::copyTo(Mat<T>& dst) const{
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        dst.data[i] = data[i];
    }
    return dst;
}

template <typename T>
Mat<T> Mat<T>::reshape(size_t new_cn, size_t new_rows, size_t new_cols) const{
    Mat<T> mat;
    mat.rows = new_rows;
    mat.cols = new_cols;
    mat.channels = new_cn;
    mat.data = data;
    mat.dataStart = dataStart;
    mat.dataEnd = dataEnd;
    mat.ref_count = ref_count;
    (*ref_count)++;
    return mat;
}

//Mat& convertTo(Mat& dst, int rtype, double alpha = 1, double beta = 0) const;
// 将矩阵转换为另一种类型。
// dst – 输出矩阵，与*this具有相同的大小和通道数。
// rtype – 输出矩阵的类型，可以是CV_8U， CV_8S， CV_16U， CV_16S， CV_32S， CV_32F或CV_64F。
// alpha – 每个元素的乘数因子。
// beta – 每个元素的加数因子。
template <typename T>
Mat<T>& Mat<T>::convertTo(Mat<T>& dst, int rtype, double alpha, double beta) const{
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        dst.data[i] = (T)(data[i] * alpha + beta);
    }
    return dst;
}

// static Mat& zeros(size_t rows, size_t cols, size_t channels);
template <typename T>
Mat<T>& Mat<T>::zeros(size_t rows, size_t cols, size_t channels){
    Mat<T> mat(rows, cols, channels);
    return mat;
}

//static Mat& ones(size_t rows, size_t cols, size_t channels);
template <typename T>
Mat<T>& Mat<T>::ones(size_t rows, size_t cols, size_t channels){
    Mat<T> mat(rows, cols, channels);
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        mat.data[i] = 1;
    }
    return mat;
}

//static Mat& eye(size_t rows, size_t cols, size_t channels);
template <typename T>
Mat<T>& Mat<T>::eye(size_t rows, size_t cols, size_t channels){
    Mat<T> mat(rows, cols, channels);
    size_t size = rows * cols * channels;
    for (size_t i = 0; i < size; i++){
        mat.data[i] = 0;
    }
    for (size_t i = 0; i < rows; i++){
        mat.data[i * cols * channels + i * channels] = 1;
    }
    return mat;
}

// Mat& append(T *nums, int r, int c, int ch);
//向矩阵尾部添加若干行数据
template <typename T>
Mat<T>& Mat<T>::append(T *data, int row, int col, int ch){
    //参数合法性检查
    if(data == NULL ){
        std::cerr << "In function append..." << std::endl;
        std::cerr << "The input array should be valid." << std::endl;
        throw "data is NULL";
    }
    if (row <= 0 || col <= 0 || ch <= 0 ){
        std::cerr << "In function append..." << std::endl;
        std::cerr << "The input row, col and ch should be positive." << std::endl;
        throw "row, col or ch is not positive";
    }
    if (ch != channels || col != cols){
        std::cerr << "In function append..." << std::endl;
        std::cerr << "The input row, col and ch should be equal to the original matrix." << std::endl;
        throw "row, col or ch is not equal to the original matrix";
    }

    int src_size = rows * cols * channels;
    int tar_size = row * col * ch;
    T * newdata = new T[src_size + tar_size];
    
    memcpy(newdata, this->data, src_size * sizeof(T));
    memcpy(newdata + src_size, data, tar_size * sizeof(T));

    if(this -> ref_count[0] == 1){
        delete[] this -> data;
    }else{
        this -> ref_count[0]--;
        this -> ref_count = new size_t(1);
    }

    this -> data = newdata;
    this -> rows = this -> rows + row;
    return *this;
}

// Matrix& subtract(int r);
//向矩阵尾部删除若干行数据
template <typename T>
Mat<T>& Mat<T>::subtract(int row){
    //参数合法性检查
    if (row <= 0){
        std::cerr << "In function subtract..." << std::endl;
        std::cerr << "The input row should be positive." << std::endl;
        throw "row is not positive";
    }
    if (row > rows){
        std::cerr << "In function subtract..." << std::endl;
        std::cerr << "The input row should be less than the original matrix." << std::endl;
        throw "row is not less than the original matrix";
    }
    int tar_size = (rows - row) * cols * channels;
    T * newdata = new T[tar_size];
    
    memcpy(newdata, this->data, tar_size * sizeof(T));

    if(this -> ref_count[0] == 1){
        delete[] this -> data;
    }else{
        this -> ref_count[0]--;
        this -> ref_count = new size_t(1);
    }

    this -> data = newdata;
    this -> rows = this -> rows - row;
    return *this;
}

// void split(vector<Mat>& mv) const;
//将多通道矩阵分离为多个单通道矩阵
template <typename T>
void Mat<T>::split(vector<Mat<T>>& mv) const{
    if (channels == 1){
        mv.push_back(*this);
        return;
    }
    for (size_t i = 0; i < channels; i++){
        Mat<T> mat(rows, cols, 1);
        for (size_t j = 0; j < rows; j++){
            for (size_t k = 0; k < cols; k++){
                mat.data[j * cols + k] = data[j * cols * channels + k * channels + i];
            }
        }
        mv.push_back(mat);
    }
}

// Mat& merge(const vector<Mat>& mv);
//将多个单通道矩阵合并为多通道矩阵
template <typename T>
Mat<T>& Mat<T>::merge(const vector<Mat<T>>& mv){
    if (mv.size() == 1){
        *this = mv[0];
        return *this;
    }
    size_t rows = mv[0].rows;
    size_t cols = mv[0].cols;
    size_t channels = mv.size();
    for (size_t i = 0; i < mv.size(); i++){
        if (mv[i].rows != rows || mv[i].cols != cols){
            std::cerr << "In function merge..." << std::endl;
            std::cerr << "The input matrix should have the same size." << std::endl;
            throw "The input matrix should have the same size";
        }
    }
    Mat<T> mat(rows, cols, channels);
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            for (size_t k = 0; k < channels; k++){
                mat.data[i * cols * channels + j * channels + k] = mv[k].data[i * cols + j];
            }
        }
    }
    *this = mat;
    return *this;
}

//Mat merge_vertical(const Mat& mat) const;
//将两个矩阵垂直合并
template <typename T>
Mat<T> Mat<T>::merge_vertical(const Mat<T>& mat) const{
    if (cols != mat.cols || channels != mat.channels){
        std::cerr << "In function merge_vertical..." << std::endl;
        std::cerr << "The input matrix should have the same size." << std::endl;
        throw "The input matrix should have the same size";
    }
    Mat<T> new_mat(rows + mat.rows, cols, channels);
    //将原矩阵的数据复制到新矩阵中
    memcpy(new_mat.data, data, rows * cols * channels * sizeof(T));
    //将输入矩阵的数据复制到新矩阵中
    memcpy(new_mat.data + rows * cols * channels, mat.data, mat.rows * mat.cols * mat.channels * sizeof(T));
    return new_mat;
}

//Mat merge_horizontal(const Mat& mat) const;
//将两个矩阵水平合并
template <typename T>
Mat<T> Mat<T>::merge_horizontal(const Mat<T>& mat) const{
    if (rows != mat.rows || channels != mat.channels){
        std::cerr << "In function merge_horizontal..." << std::endl;
        std::cerr << "The input matrix should have the same size." << std::endl;
        throw "The input matrix should have the same size";
    }
    Mat<T> new_mat(rows, cols + mat.cols, channels);
    //将原矩阵的数据复制到新矩阵中
    for (size_t i = 0; i < rows; i++){
        memcpy(new_mat.data + i * new_mat.cols * channels, data + i * cols * channels, cols * channels * sizeof(T));
    }
    //将输入矩阵的数据复制到新矩阵中
    for (size_t i = 0; i < mat.rows; i++){
        memcpy(new_mat.data + i * new_mat.cols * channels + cols * channels, mat.data + i * mat.cols * mat.channels, mat.cols * mat.channels * sizeof(T));
    }
    return new_mat;
}

//将彩色图像转换为三通道黑白图像
template <typename T>
Mat<T> Mat<T>::rgbtogray() const{
    if (channels != 3){
        std::cerr << "In function rgbtogray..." << std::endl;
        std::cerr << "The input matrix should be a color image." << std::endl;
        throw "The input matrix should be a color image";
    }
    Mat<T> gray(rows, cols, 3);
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            T val = (data[i * cols * channels + j * channels] * 0.299 + data[i * cols * channels + j * channels + 1] * 0.587 + data[i * cols * channels + j * channels + 2] * 0.114);
            gray.data[i * cols * channels + j * channels] = val;
            gray.data[i * cols * channels + j * channels + 1] = val;
            gray.data[i * cols * channels + j * channels + 2] = val;
        }
    }
    return gray;
}
    



















    








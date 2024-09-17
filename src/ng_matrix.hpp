#ifndef _NG_MATRIX_HPP_
#define _NG_MATRIX_HPP_
#include <ostream>
#include "ng_mixed_buffer.hpp" // vector

namespace ngineer
{
    typedef std::pair<buflen_t, buflen_t> Indices;

    template<class T>
    class Matrix
    {
    private:
        buflen_t rows;
        buflen_t cols;
        MixedBuffer<T> elements;
        void inplace_invert2(void);
        void inplace_invert3(void);
        void inplace_invert4(void);
        void inplace_invertN(void);
    protected:
        void inplace_pivot_check(buflen_t col);
        void inplace_pivot_check(buflen_t col, Matrix<T> &mirror);
    public:
        Matrix();
        Matrix(buflen_t rows, buflen_t cols);
        Matrix(buflen_t cols, std::vector<T> elems);
        ~Matrix();
        static Matrix<T> identity(buflen_t n);
        bool is_square(void);
        std::vector<Indices> get_indices_iter(void) const noexcept;
        void inplace_row_swap(buflen_t r1, buflen_t r2);
        void inplace_row_scale(buflen_t row, T scalar);
        void inplace_row_add(buflen_t r1, buflen_t r2);
        void inplace_scaled_row_add(buflen_t r1, buflen_t r2, T scalar);
        void inplace_invert(void);
        // Matrix<T> lu_decomposition(void);
        void inplace_transpose(void);
        template<class U>
        friend std::ostream &operator<<(std::ostream &out, Matrix<U> &a);
        T& operator[](Indices index);
        Matrix<T> operator|(Matrix<T> rhs);
        Matrix<T> operator^(char rhs);
        Matrix<T> operator+(Matrix<T> rhs);
        Matrix<T> operator-(Matrix<T> rhs);
        Matrix<T> operator*(Matrix<T> rhs);
        Matrix<T> operator*(T rhs);
        Matrix<T> operator/(Matrix<T> rhs);
    }; // class Matrix<T>

    template<class U>
    inline std::ostream &operator<< (std::ostream &out, Matrix<U>& a)
    {
        buflen_t numElems = a.rows * a.cols;
        out << '[';
        for (buflen_t i = 0; i < a.rows; i++)
        {
            for (buflen_t j = 0; j < a.cols; j++)
            {
                out << a[{ i, j }];
                if (a.cols - 1 != j) // for all elements but the last one...
                {
                    out << ", ";
                }
            }
            if (a.rows - 1 != i) // for all rows but the last one...
            {
                out << "; ";
            }
            else
            {
                out << ']';
            }
        }
        return out;
    }
} // namespace ngineer

template class ngineer::Matrix<char>;
template class ngineer::Matrix<unsigned char>;
template class ngineer::Matrix<short>;
template class ngineer::Matrix<unsigned short>;
template class ngineer::Matrix<int>;
template class ngineer::Matrix<unsigned int>;
template class ngineer::Matrix<long>;
template class ngineer::Matrix<unsigned long>;
template class ngineer::Matrix<long long>;
template class ngineer::Matrix<unsigned long long>;
template class ngineer::Matrix<float>;
template class ngineer::Matrix<double>;

#endif // _MATRIX_HPP_
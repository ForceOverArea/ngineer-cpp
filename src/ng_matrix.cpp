#include <sstream>
#include <stdexcept>
#include "ng_matrix.hpp" // ng_mixed_buffer.hpp, ng_typedefs.hpp, functional, vector
namespace ngineer
{
    template<class T>
    Matrix<T>::Matrix()
    {
        rows = 0;
        cols = 0;
        elements = MixedBuffer<T>(0);
    }

    template<class T>
    Matrix<T>::Matrix(buflen_t rows, buflen_t cols)
    : rows(rows), cols(cols), elements({ rows * cols })
    {
    }

    template<class T>
    Matrix<T>::Matrix(buflen_t cols, std::vector<T> elems)
    : rows(elems.size() / cols), cols(cols), elements(elems)
    {
        if (elems.size() % cols != 0)
        {
            throw std::runtime_error { 
                "`elems` vector used to initialize matrix must have a size evenly divisible by `cols`" };
        }
    }

    template<class T>
    Matrix<T>::~Matrix()
    {
    }

    template<class T>
    Matrix<T> Matrix<T>::identity(buflen_t n) 
    {
        Matrix<T> a { n, n };

        for (buflen_t i = 0; i < n; i++)
        {
            a[{ i, i }] = 1;
        }

        return a;
    }

    template<class T>
    bool Matrix<T>::is_square(void)
    {
        return rows == cols;
    }

    template<class T>
    std::vector<Indices> Matrix<T>::get_indices_iter(void) const noexcept
    {
        std::vector<Indices> range;
        for (buflen_t i = 0; i < rows; i++)
        {
            for (buflen_t j = 0; j < cols; j++)
            {
                range.push_back({ i, j });
            }
        }
        return range;
    }

    template<class T>
    void Matrix<T>::inplace_row_swap(buflen_t r1, buflen_t r2)
    {
        T tmp = 0x00;
        for (buflen_t i = 0; i < cols; i++)
        {
            tmp = (*this)[{ r1 , i }];
            (*this)[{ r1, i }] = (*this)[{ r2, i }];
            (*this)[{ r2, i }] = tmp;
        }
    }

    template<class T>
    void Matrix<T>::inplace_row_scale(buflen_t row, T scalar)
    {
        for (buflen_t i = 0; i < cols; i++)
        {
            (*this)[{ row, i }] *= scalar;
        }
    }

    template<class T>
    void Matrix<T>::inplace_row_add(buflen_t r1, buflen_t r2)
    {
        for (buflen_t i = 0; i < cols; i++)
        {
            (*this)[{ r2, i }] += (*this)[{ r1, i }];
        }
    }

    template<class T>
    void Matrix<T>::inplace_scaled_row_add(buflen_t r1, buflen_t r2, T scalar)
    {
        for (buflen_t i = 0; i < cols; i++)
        {
            (*this)[{ r2, i }] += (*this)[{ r1, i }] * scalar;
        }
    }

    template<class T>
    void Matrix<T>::inplace_pivot_check(buflen_t col_n)
    {
        buflen_t ideal_pivot_row = col_n;
        while (0 == (*this)[{ ideal_pivot_row, col_n }])
        {
            if (ideal_pivot_row >= rows)
            {
                std::stringstream ss;
                ss << "no nonzero pivots found while searching in column " << col_n;
                throw std::runtime_error { ss.str() };
            }
            ideal_pivot_row++;
        }
        inplace_row_swap(ideal_pivot_row, col_n);
    }

    template<class T>
    void Matrix<T>::inplace_pivot_check(buflen_t col_n, Matrix<T> &mirror)
    {
        buflen_t ideal_pivot_row = col_n;
        while (0 == (*this)[{ ideal_pivot_row, col_n }])
        {
            if (ideal_pivot_row >= rows)
            {
                std::stringstream ss;
                ss << "no nonzero pivots found while searching in column " << col_n;
                throw std::runtime_error { ss.str() };
            }
            ideal_pivot_row++;
        }
        inplace_row_swap(ideal_pivot_row, col_n);
        mirror.inplace_row_swap(ideal_pivot_row, col_n);
    }

    template<class T>
    void Matrix<T>::inplace_invert2(void)
    {
        T a11 = (*this)[{ 0, 0 }];
        T a12 = (*this)[{ 0, 1 }];
        T a21 = (*this)[{ 1, 0 }];
        T a22 = (*this)[{ 1, 1 }];

        T det = a11*a22 - a12*a21;

        if (det == 0)
        {
            throw std::runtime_error { 
                "the determinant of the matrix was 0" };
        }

        (*this)[{ 0, 0 }] =   a22 / det;
        (*this)[{ 1, 0 }] = - a21 / det;
        (*this)[{ 0, 1 }] = - a12 / det;
        (*this)[{ 1, 1 }] =   a11 / det;
    }

    template<class T>
    void Matrix<T>::inplace_invert3(void)
    {
        T a11 = (*this)[{ 0, 0 }];
        T a12 = (*this)[{ 0, 1 }];
        T a13 = (*this)[{ 0, 2 }];
        T a21 = (*this)[{ 1, 0 }];
        T a22 = (*this)[{ 1, 1 }];
        T a23 = (*this)[{ 1, 2 }];
        T a31 = (*this)[{ 2, 0 }];
        T a32 = (*this)[{ 2, 1 }];
        T a33 = (*this)[{ 2, 2 }];

        T det       = a11*a22*a33 + a21*a32*a13 + a31*a12*a23 
                    - a11*a32*a23 - a31*a22*a13 - a21*a12*a33;

        if (det == 0)
        {
            throw std::runtime_error { 
                "the determinant of the matrix was 0" };
        }

        (*this)[{ 0, 0 }] = (a22 * a33 - a23 * a32) / det;
        (*this)[{ 1, 0 }] = (a23 * a31 - a21 * a33) / det;
        (*this)[{ 2, 0 }] = (a21 * a32 - a22 * a31) / det;
        (*this)[{ 0, 1 }] = (a13 * a32 - a12 * a33) / det;
        (*this)[{ 1, 1 }] = (a11 * a33 - a13 * a31) / det;
        (*this)[{ 2, 1 }] = (a12 * a31 - a11 * a32) / det;
        (*this)[{ 0, 2 }] = (a12 * a23 - a13 * a22) / det;
        (*this)[{ 1, 2 }] = (a13 * a21 - a11 * a23) / det;
        (*this)[{ 2, 2 }] = (a11 * a22 - a12 * a21) / det;
    }

    template<class T>
    void Matrix<T>::inplace_invert4(void)
    {
        T a11 = (*this)[{ 0, 0 }];
        T a12 = (*this)[{ 0, 1 }];
        T a13 = (*this)[{ 0, 2 }];
        T a14 = (*this)[{ 0, 3 }];
        T a21 = (*this)[{ 1, 0 }];
        T a22 = (*this)[{ 1, 1 }];
        T a23 = (*this)[{ 1, 2 }];
        T a24 = (*this)[{ 1, 3 }];
        T a31 = (*this)[{ 2, 0 }];
        T a32 = (*this)[{ 2, 1 }];
        T a33 = (*this)[{ 2, 2 }];
        T a34 = (*this)[{ 2, 3 }];
        T a41 = (*this)[{ 3, 0 }];
        T a42 = (*this)[{ 3, 1 }];
        T a43 = (*this)[{ 3, 2 }];
        T a44 = (*this)[{ 3, 3 }];

        T det     = a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 +
                    a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + 
                    a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + 
                    a14*a21*a33*a42 + a14*a22*a34*a43 + a14*a23*a32*a41 -
                    a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 -
                    a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 -
                    a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 -
                    a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42;

        if (det == 0)
        {
            throw std::runtime_error { 
                "the determinant of the matrix was 0" };
        }

        (*this)[{ 0, 0 }] = (a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42) / det;
        (*this)[{ 1, 0 }] = (a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43) / det;
        (*this)[{ 2, 0 }] = (a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42) / det;
        (*this)[{ 3, 0 }] = (a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33) / det;
        (*this)[{ 0, 1 }] = (a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43) / det;
        (*this)[{ 1, 1 }] = (a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41) / det;
        (*this)[{ 2, 1 }] = (a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43) / det;
        (*this)[{ 3, 1 }] = (a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31) / det;
        (*this)[{ 0, 2 }] = (a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41) / det;
        (*this)[{ 1, 2 }] = (a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42) / det;
        (*this)[{ 2, 2 }] = (a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41) / det;
        (*this)[{ 3, 2 }] = (a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32) / det;
        (*this)[{ 0, 3 }] = (a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42) / det;
        (*this)[{ 1, 3 }] = (a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41) / det;
        (*this)[{ 2, 3 }] = (a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42) / det;
        (*this)[{ 3, 3 }] = (a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31) / det;
    }

    template<class T>
    void Matrix<T>::inplace_invertN(void)
    {
        const buflen_t n = rows;
        auto inv = Matrix<T>::identity(n);

        for (buflen_t j = 0; j < n; j++)
        {
            for (buflen_t i = 0; i < n; i++)
            {
                if (i == j)
                {
                    continue;
                }

                if ((*this)[{ j, j }] == 0)
                {
                    throw std::runtime_error { 
                        "division by zero was required to invert the given matrix" };
                }

                T scalar = (*this)[{ i, j }] / (*this)[{ j, j }];
                this->inplace_scaled_row_add(i, j, -scalar);
                inv.inplace_scaled_row_add(i, j, -scalar);
            }
        }

        for (buflen_t i = 0; i < n; i++)
        {
            T scalar = 1 / (*this)[{ i, i }];
            this->inplace_row_scale(i, scalar);
            inv.inplace_row_scale(i, scalar);
        }

        this->elements = inv.elements;
    }

    template<class T>
    void Matrix<T>::inplace_invert(void)
    {
        if (rows != cols)
        {
            throw std::runtime_error { 
                "matrix must be square to invert" };
        }

        if ((rows == 1) && (elements[0] == 0))
        {
            throw std::runtime_error {
                "the scalar value 0 is non-invertible" };
        }

        switch(rows)
        {
        case 1:
            elements[0] = 1 / elements[0];
            break;
        case 2:
            inplace_invert2();
            break;
        case 3:
            inplace_invert3();
            break;
        case 4:
            inplace_invert4();
            break;
        default:
            inplace_invertN();
            break;
        }
    }

    // template<class T>
    // Matrix<T> Matrix<T>::lu_decomposition(void)
    // {
    //     auto u = Matrix<T>(*this);
    //     auto l = Matrix<T>::identity(rows);
    //     auto p = Matrix<T>::identity(rows);

    //     if (!is_square())
    //     {
    //         throw std::runtime_error {
    //             "non-square matrices are not yet supported" };
    //     }

    //     // set all elements below the diagonal to 0 for u
    //     for (buflen_t n = 0; n < cols; n++)
    //     {
    //         // Make sure there is a nonzero pivot on the diagonal
    //         try
    //         {
    //             u.inplace_pivot_check(rows, p);
    //             for (buflen_t j = 0; j < rows; j++)
    //             {
    //                 l[{ j, n }] = u[{ j, n }];
    //             }
    //         }
    //         catch (std::runtime_error e)
    //         {

    //         }
    //     }

    //     return l * u;
    // }

    template<class T>
    void Matrix<T>::inplace_transpose(void)
    {
        Matrix<T> transpose = { cols, rows };

        for (buflen_t i = 0; i < rows; i++)
        {
            for (buflen_t j = 0; j < cols; j++)
            {
                transpose[{ j, i }] = (*this)[{ i, j }];
            }
        }
    }

    template<class T>
    T& Matrix<T>::operator[](Indices index)
    {   
        if (index.first >= rows || index.second >= cols)
        {
            std::stringstream ss;
            ss << "index is out of bounds (matrix has " 
               << rows << " rows and " << cols << " columns, "
               << "but user tried to access {" 
               << index.first << ',' << index.second << "})";
            throw std::runtime_error { ss.str() };
        }

        return elements[cols * index.first + index.second];
    }

    template<class T>
    Matrix<T> Matrix<T>::operator|(Matrix<T> rhs)
    {
        Matrix<T> augment = { rows, cols + rhs.cols };

        if (rows != rhs.rows)
        {
            throw std::runtime_error {
                "matrices have mismatching numbers of rows" };
        }

        for (buflen_t i = 0; i < rows; i++)
        {
            for (buflen_t j = 0; j < cols; j++)
            {
                if (j < cols)
                {
                    augment[{ i, j }] = (*this)[{ i, j }];
                }
                else
                {
                    augment[{ i, j }] = rhs[{ i, j - cols }];
                }
            }
        }

        return augment;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator^(char rhs)
    {
        Matrix<T> result;

        switch(rhs)
        {
        case -1:
            result = (*this);
            result.inplace_invert();
            break;
        case 'T':
            result = (*this);
            result.inplace_transpose();
            break;
        default:
            throw std::runtime_error {
                "matrices can only use `-1` and `T` with the `^` operator" };
            break;
        }

        return result;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator+(Matrix<T> rhs)
    {
        Matrix<T> sum { rows, cols };

        if (rows != rhs.rows || cols != rhs.cols)
        {
            throw std::runtime_error {
                "matrices do not have the same shape" };
        }

        for (buflen_t i = 0; i < elements.get_length(); i++)
        {
            sum.elements[i] = elements[i] + rhs.elements[i];
        }

        return sum;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator-(Matrix<T> rhs)
    {
        Matrix<T> difference { rows, cols };

        if (rows != rhs.rows || cols != rhs.cols)
        {
            throw std::runtime_error {
                "matrices do not have the same shape" };
        }

        for (buflen_t i = 0; i < elements.get_length(); i++)
        {
            difference.elements[i] = elements[i] - rhs.elements[i];
        }

        return difference;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator*(Matrix<T> rhs)
    {
        Matrix<T> product { rows, cols };
        buflen_t n = cols;

        if (cols != rhs.rows)
        {
            throw std::runtime_error {
                "matrices do not have the correct shapes" };
        }

        for (buflen_t i = 0; i < rows; i++)
        {
            for (buflen_t j = 0; j < cols; j++)
            {
                for (buflen_t x = 0; x < n; x++)
                {
                    product[{ i, j }] += (*this)[{ i, x }] * rhs[{ x, j }];
                }
            }
        }

        return product;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator*(T rhs)
    {
        Matrix a = (*this);

        for (buflen_t i = 0; i < rows; i++)
        {
            a.inplace_row_scale(i, rhs);
        }

        return a;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator/(Matrix<T> rhs)
    {
        rhs.inplace_invert();
        return rhs * (*this);
    }
} // namespace ngineer

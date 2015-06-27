#ifndef DENSE_MATRIX_H
#define DENSE_MATRIX_H
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
namespace igmdk{

template<typename ITEM> struct DenseMatrix
{
    int rows, columns;
    int index(int row, int column)const
    {
        assert(row >= 0 && row < rows && column >= 0 && column < columns);
        return row + column * rows;
    }
    Vector<ITEM> items;
public:
    ITEM& operator()(int row, int column){return items[index(row, column)];}
    ITEM const& operator()(int row, int column)const
        {return items[index(row, column)];}
    DenseMatrix(int theRows, int theColumns): rows(theRows),
        columns(theColumns), items(rows * columns, rows * columns){}

    DenseMatrix operator*=(ITEM const& scalar)
    {
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j) (*this)(i, j) *= scalar;
        return *this;
    }
    friend DenseMatrix operator*(ITEM const& scalar, DenseMatrix const& a)
    {
        DenseMatrix result(a);
        return result *= scalar;
    }
    friend DenseMatrix operator*(DenseMatrix const& a, ITEM const& scalar)
        {return scalar * a;}

    DenseMatrix& operator+=(DenseMatrix const& rhs)
    {
        assert(rows == rhs.rows && columns = rhs.columns);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j) (*this)(i, j) += rhs(i, j);
        return *this;
    }
    friend DenseMatrix operator+(DenseMatrix const& a, DenseMatrix const& b)
    {
        DenseMatrix result(a);
        return result += b;
    }
    DenseMatrix& operator-=(DenseMatrix const& rhs)
    {
        (*this) *= -1;
        (*this) += rhs;
        (*this) *= -1;
    }
    friend DenseMatrix operator-(DenseMatrix const& a, DenseMatrix const& b)
    {
        DenseMatrix result(a);
        return result -= b;
    }

    friend DenseMatrix operator*(DenseMatrix const& a, DenseMatrix const& b)
    {
        assert(a.columns == b.rows);
        DenseMatrix result(a.rows, b.columns);
        for(int i = 0; i < a.rows; ++i)
            for(int j = 0; j < b.columns; ++j)
            {
                ITEM sum(0);
                for(int k = 0; k < b.rows; ++k) sum += a(i, k) * b(k, j);
                result(i, j) += sum;
            }
        return result;
    }
    DenseMatrix operator*=(DenseMatrix const& rhs)
        {return *this = *this * rhs;}
    Vector<ITEM> operator*(Vector<ITEM> const& v)
    {
        assert(columns == v.getSize());
        Vector<ITEM> result(rows, rows);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j)
                result[i] += (*this)(i, j) * v[j];
        return result;
    }
    friend Vector<ITEM> operator*(Vector<ITEM> const& v,
        DenseMatrix const& m){return m.transpose() * v;}

    static DenseMatrix identity(int n)
    {
        DenseMatrix result(n, n);
        for(int i = 0; i < n; ++i) result(i, i) = ITEM(1);
        return result;
    }
    DenseMatrix transpose()const
    {
        DenseMatrix result(columns, rows);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j) result(j, i) = (*this)(i, j);
        return result;
    }

    void debug()
    {
        for(int i = 0; i < rows; ++i)
        {
            for(int j = 0; j < columns; ++j)
            {
                cout << (*this)(i, j) << " ";
            }
            cout << endl;
        }
    }
};

template<typename ITEM> struct DenseLUP
{
    DenseMatrix<ITEM> d;
    Vector<int> permutation;
    bool isSingular;
    DenseLUP(DenseMatrix<ITEM> const& a): d(a), permutation(a.rows),
        isSingular(false)
    {
        assert(d.rows = d.columns);
        for(int i = 0; i < d.rows; ++i) permutation.append(i);
        for(int i = 0; i < d.rows; ++i)
        {
            ITEM p = 0;
            int entering = -1;
            for(int j = i; j < d.rows; ++j)
                if(abs(d(i, j)) > p)
                {
                    p = abs(d(i, j));
                    entering = i;
                }
            if(entering == -1){isSingular = true; return;}
            swap(permutation[i], permutation[entering]);
            for(int j = 0; j < d.rows; ++j) swap(d(i, j), d(entering, j));
            for(int j = i + 1; j < d.rows; ++j)
            {
                d(j, i) /= d(i, i);
                for(int k = i + 1; k < d.rows; ++k)
                    d(j, k) -= d(j, i) * d(i, k);
            }
        }
    }

    ITEM determinant()
    {
        if(isSingular) return 0;
        ITEM result(1);
        for(int i = 0; i < d.rows; ++i)
        {
            result *= d(i, i);
            if(permutation[i] % 2 != i % 2) result *= -1;
        }
        return result;
    }

    Vector<ITEM> solve(Vector<ITEM> const& b)
    {
        assert(!isSingular);
        Vector<ITEM> result(d.rows, d.rows), y(result);
        for(int i = 0; i < d.rows; ++i)
        {
            y[i] = b[permutation[i]];
            for(int j = 0; j < i; ++j) y[i] -= y[j] * d(i, j);
        }
        for(int i = d.rows - 1; i >= 0; --i)
        {
            result[i] = y[i];
            for(int j = i + 1; j < d.rows; ++j)
                result[i] -= result[j] * d(i, j);
            result[i] /= d(i, i);
        }
        return result;
    }

    DenseMatrix<ITEM> inverse()
    {
        assert(!isSingular);
        Vector<ITEM> identityRow(d.rows, d.rows, 0);
        DenseMatrix<ITEM> result(d.rows, d.rows);
        for(int i = 0; i < d.columns; ++i)
        {
            identityRow[i] = 1;
            Vector<ITEM> column = solve(identityRow);
            identityRow[i] = 0;
            for(int j = 0; j < d.rows; ++j) result(j, i) = column[j];
        }
        return result;
    }
};

struct LinearProgrammingSimplex
{
    DenseMatrix<double> B, N;
    Vector<double> b, cB, cN, x;
    Vector<int> p;
    bool isUnbounded;
    bool performIteration()
    {
        DenseMatrix<double> InvB = DenseLUP<double>(B).inverse();
        x = InvB * b;
        //check if x is optimal or find entering variable
        Vector<double> y = cN - cB * InvB * N;
        int entering = 0;
        double bestValue = y[0];
        for(int i = 1; i < y.getSize(); ++i)
            if(y[i] < bestValue)
            {
                bestValue = y[i];
                entering = i;
            }
        if(bestValue >= 0) return false;
        //find leaving variable
        Vector<double> a;
        for(int i = 0; i < N.rows; ++i) a.append(N(i, entering));
        a = InvB * a;
        int leaving = -1;
        double minRatio, maxA = -1;
        for(int i = 0; i < x.getSize(); ++i)
            if(a[i] > 0)
            {
                double newRatio = x[i]/a[i];
                if(leaving == -1 || minRatio > newRatio)
                {
                    leaving = i;
                    maxA = max(maxA, a[i]);
                    minRatio = newRatio;
                }
            }
        if(maxA <= 0){isUnbounded = true; return false;}
        //swap variables
        for(int i = 0; i < N.rows; ++i){swap(B(i, leaving), N(i, entering));}
        swap(p[leaving], p[entering]);
        swap(cB[leaving], cN[entering]);
        return true;
    }
    LinearProgrammingSimplex(DenseMatrix<double>const&B0, DenseMatrix<double>
        const& N0, Vector<double>const& cB0, Vector<double>const& cN0,
        Vector<double> const& b0): isUnbounded(false), B(B0), N(N0), cB(cB0),
        cN(cN0), b(b0), x(b)
        {for(int i = 0; i < cB.getSize() + cN.getSize(); ++i) p.append(i);}
    Vector<pair<int, double> > solve()
    {
        while(performIteration());
        Vector<pair<int, double> > result;
        if(!isUnbounded)
            for(int i = 0; i < x.getSize(); ++i)
                result.append(make_pair(p[i], x[i]));
        return result;
    }
};

}
#endif

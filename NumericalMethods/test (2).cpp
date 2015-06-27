#include "LargeNumber.h"
#include "DenseMatrix.h"
#include "../Utils/Debug.h"
using namespace igmdk;

void testSimplex()
{
    DenseMatrix<double> B = DenseMatrix<double>::identity(3), N(3, 2);
    N(0, 0) = -2;
    N(0, 1) = 1;
    N(1, 0) = -1;
    N(1, 1) = 2;
    N(2, 0) = 1;
    N(2, 1) = 0;
    Vector<double> b, cB(3, 3, 0), cN;
    b.append(2);
    b.append(7);
    b.append(3);
    cN.append(-1);
    cN.append(-2);
    LinearProgrammingSimplex s(B, N, cB, cN, b);
    Vector<pair<int, double> > result = s.solve();
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i].first);
        DEBUG(result[i].second);
    }
}

int main()
{
    int e = 0;

    DEBUG(frexp(10, &e));
    DEBUG(rationalize(10, e));
    DEBUG(e);
    return 0;

    testSimplex();

    assert((-Number(0)-Number(2)) == (Number(0) - Number(2)));
    assert((Number(2)-Number(0)) == (Number(2) - -Number(0)));

    DenseMatrix<double> a(3, 3);
    a(0, 0) = 1;
    a(0, 1) = 2;
    a(0, 2) = 0;
    a(1, 0) = 3;
    a(1, 1) = 4;
    a(1, 2) = 4;
    a(2, 0) = 5;
    a(2, 1) = 6;
    a(2, 2) = 3;
    DenseLUP<double> lup(a);
    Vector<double> b;
    b.append(3);
    b.append(7);
    b.append(8);
    Vector<double> x = lup.solve(b);
    for(int i = 0; i < b.getSize(); ++i) DEBUG(x[i]);

    DenseMatrix<double> b1(3, 3);
    b1(0, 0) = 1;
    b1(0, 1) = -1;
    b1(0, 2) = 2;
    b1(1, 0) = -2;
    b1(1, 1) = 1;
    b1(1, 2) = 1;
    b1(2, 0) = -1;
    b1(2, 1) = 2;
    b1(2, 2) = 1;

    DenseLUP<double> lup2(b1);

    DenseMatrix<double> c = lup2.inverse();

    for(int i = 0; i < c.rows; ++i)
    {
        for(int j = 0; j < c.columns; ++j)
        {
            cout << (c(i, j)) << " ";
        }
        cout << endl;
    }

    DenseMatrix<double> b2(3, 3);
    b2(0, 0) = 1;
    b2(0, 1) = 4;
    b2(0, 2) = 5;
    b2(1, 0) = 4;
    b2(1, 1) = 20;
    b2(1, 2) = 32;
    b2(2, 0) = 5;
    b2(2, 1) = 32;
    b2(2, 2) = 64;

    PrimeTable pt(1000000);
    DEBUG(pt.isPrime(91));
    DEBUG(pt.isPrime(107));
    DEBUG(pt.isPrime(17));
    DEBUG(pt.isPrime(1));
    DEBUG(pt.isPrime(2));
    DEBUG(pt.isPrime(3));
    for(int i = 0; i < 100; ++i)
    {
        DEBUG(i);
        DEBUG(pt.isPrime(i));
    }
    Number(99).sqrt().debug();

    (Number(-11) % Number(103)).debug();

    Number n(2);
    Number m = n.power(128);
    m.debug();
    m >>= 125;
    m.debug();
    m <<= 125;
    m.debug();
    m = gcd(m, Number(3));
    m = modInverse(Number(4), Number(7));
    m.debug();
    DEBUG(isPrime(Number((8616460799ull))));//no
    DEBUG(isPrime(Number((3ull))));//yes

    Vector<unsigned char> igits = m.toDigitVector();
    for(int i = 0; i < igits.getSize(); ++i) DEBUG(int(igits[i]));

	return 0;
}

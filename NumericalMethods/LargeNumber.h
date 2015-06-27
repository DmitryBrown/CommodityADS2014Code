#ifndef NUMBER_H
#define NUMBER_H
#include "../MiscAlgs/Misc.h"
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/Vector.h"
using namespace std;
namespace igmdk{

class Number
{
    bool isMinus;
    typedef unsigned int DIGIT;
    typedef unsigned long long LARGE_DIGIT;
    enum{BASE_RADIX = numeric_limits<DIGIT>::digits};
    Vector<DIGIT> digits;
public:
    typedef DIGIT DIGIT_TYPE;
    DIGIT getDigit(int i)const{return i < nDigits() ? digits[i] : 0;}
    bool isZero()const{return digits.lastItem() == 0;}
    bool isPositive()const{return !isMinus && !isZero();}
    bool isNegative()const{return isMinus && !isZero();}
    int nDigits()const{return digits.getSize();}
    DIGIT& operator[](unsigned int i){return digits[i];}
    DIGIT const& operator[](unsigned int i)const{return digits[i];}
    void trim(){while(nDigits() > 1 && isZero()) digits.removeLast();}
    void negate(){isMinus = !isMinus;}
    Number operator-()const
    {
        Number result = *this;
        result.negate();
        return result;
    }
    Number abs()const{return isMinus ? -*this : *this;}
    bool isOdd()const{return digits[0] % 2;}
    bool isEven()const{return !isOdd();}
    Number(): isMinus(false){digits.append(0);}
    void appendDigit(DIGIT const& digit){digits.append(digit);}

    void constructFrom(LARGE_DIGIT digit, bool theIsNegative)
    {
        assert(BASE_RADIX * 2 <= numeric_limits<LARGE_DIGIT>::digits);
        isMinus = theIsNegative;
        digits.append(digit);
        DIGIT first = digit >> BASE_RADIX;
        if(first != 0) digits.append(first);
    }
    Number(unsigned int x){constructFrom(x, false);}
    Number(unsigned long long x){constructFrom(x, false);}
    Number(int x){constructFrom(x < 0 ? -x : x, x < 0);}
    Number(long long x){constructFrom(x < 0 ? -x : x, x < 0);}
    Number(int size, unsigned long long fill): digits(size, size, fill),
        isMinus(false){}

    bool absLess(Number const& rhs)const
    {
        if(nDigits() != rhs.nDigits()) return nDigits() < rhs.nDigits();
        for(int i = nDigits() - 1; i >= 0; --i)
            if(digits[i] != rhs[i]) return digits[i] < rhs[i];
        return false;
    }
    bool absEqual(Number const& rhs)const
    {
        if(nDigits() != rhs.nDigits()) return false;
        for(int i = 0; i < nDigits(); ++i)
            if(digits[i] != rhs[i]) return false;
        return true;
    }
    bool operator<(Number const& rhs)const
        {return (isMinus && !rhs.isMinus && !isZero()) || absLess(rhs);}
    bool operator==(Number const& rhs)const
        {return (isMinus == rhs.isMinus || isZero()) && absEqual(rhs);}

    static DIGIT fullAdder(DIGIT a, DIGIT b, DIGIT& carry)
    {
        LARGE_DIGIT sum = LARGE_DIGIT(a) + b + carry;
        carry = sum >> BASE_RADIX;
        return sum;
    }
    static Number add(Number const& a, Number const& b)
    {//O(|a|+|b|)
        int n = max(a.nDigits(), b.nDigits());
        Number result(n + 1, 0);
        DIGIT carry = 0;
        for(int i = 0; i < n; ++i)
            result[i] = fullAdder(a.getDigit(i), b.getDigit(i), carry);
        result[n] = carry;
        result.trim();
        return result;
    }

    static void sub(Number& a, Number const& b)
    {//O(|a| + |b|)
        bool carry = 0;
        for(int i = 0; i < a.nDigits(); ++i)
        {
            LARGE_DIGIT digit = LARGE_DIGIT(b.getDigit(i)) + carry;
            a[i] -= digit;
            carry = a[i] < 0;
        }
        a.trim();
    }

    friend Number operator+(Number const& a, Number const& b)
    {
        if(a.isMinus == b.isMinus)
        {
            Number result = add(a, b);
            result.isMinus = a.isMinus;
            return result;
        }
        else return a - -b;
    }
    Number& operator+=(Number const&rhs){return *this = *this + rhs;}
    Number& operator++(){return *this += 1;}
    Number operator++(int)
    {
        Number result = *this;
        ++*this;
        return result;
    }
    friend Number operator-(Number const& a, Number const& b)
    {
        if(a.isMinus == b.isMinus)
        {
            bool less = a.absLess(b);
            Number larger = less ? b : a;
            sub(larger, less ? a : b);
            if(less) larger.isMinus = !larger.isMinus;
            larger.trim();
            return larger;
        }
        else return a + -b;
    }
    Number& operator-=(Number const&rhs){return *this = *this - rhs;}
    Number& operator--(){return *this -= 1;}
    Number operator--(int)
    {
        Number result = *this;
        --*this;
        return result;
    }

    Number& operator>>=(unsigned int k)
    {
        int skipCells = k/BASE_RADIX, last = nDigits() - skipCells,
            skipBits = k % BASE_RADIX;
        if(skipCells > 0)
            for(int i = 0; i < nDigits(); ++i)
                digits[i] = i < last ? digits[i + skipCells] : 0;
        if(skipBits > 0)
        {
            DIGIT carry = 0, tempCarry;
            for(int i = last - 1; i >= 0; --i)
            {
                tempCarry = digits[i] << (BASE_RADIX - skipBits);
                digits[i] = (digits[i] >> skipBits) | carry;
                carry = tempCarry;
            }
        }
        trim();
        return *this;
    }
    Number operator>>(unsigned int k)const{return Number(*this) >>= k;}
    Number operator<<(unsigned int k)const
    {
        int skipCells = k/BASE_RADIX, skipBits = k % BASE_RADIX;
        Number result(nDigits() + 1 + skipCells, 0);
        result.isMinus = isMinus;
        for(int i = 0; i < result.nDigits(); ++i) result[i] = getDigit(i);
        if(skipCells > 0)
            for(int i = result.nDigits() - 1; i >= 0; --i)
                result[i] = i < skipCells ? 0 : result[i - skipCells];
        if(skipBits > 0)
        {
            DIGIT carry = 0;
            for(int i = skipCells; i < result.nDigits(); ++i)
            {
                DIGIT tempCarry = result[i] >> (BASE_RADIX - skipBits);
                result[i] = (result[i] << skipBits) | carry;
                carry = tempCarry;
            }
        }
        result.trim();
        return result;
    }
    Number& operator<<=(unsigned int k){return *this = *this << k;}

    static DIGIT digitMult(DIGIT a, DIGIT b, DIGIT& carry)
    {
        LARGE_DIGIT prod = LARGE_DIGIT(a) * b;
        carry = prod >> BASE_RADIX;
        return prod;
    }
    static Number mult(Number const& a, DIGIT const& b)
    {
        Number result(a.nDigits() + 1, 0);
        DIGIT carry = 0, multCarry = 0, temp;
        for(int i = 0; i < a.nDigits(); ++i)
        {
            result[i] = fullAdder(digitMult(a[i], b, temp), multCarry,
                carry);
            multCarry = temp;
        }
        result[a.nDigits()] = fullAdder(0, multCarry, carry);
        result.trim();
        return result;
    }
    friend Number operator*(Number const&a, Number const& b)
    {//O(|a| * |b|)
        Number product(a.nDigits() + b.nDigits(), 0);
        for(int j = 0; j < b.nDigits(); ++j)
            product += mult(a, b[j]) << BASE_RADIX * j;
        product.isMinus = a.isMinus != b.isMinus;
        product.trim();
        return product;
    }
    Number& operator*=(Number const&rhs){return *this = *this * rhs;}

    Vector<unsigned char> toDigitVector()const
    {
        Vector<unsigned char> result;
        Number r = *this;
        while(!r.isZero())
        {
            Number q(0);
            result.append(divide(r, 10, q)[0]);
            r = q;
        }
        result.reverse();
        return result;
    }
    Number(Vector<unsigned char> const& digitVector,
       bool theIsNegative = false): isMinus(theIsNegative)
    {
        Number result(0);
        for(int i = digitVector.getSize() - 1; i >= 0; --i)
        {
            result *= Number(10);
            result += Number(digitVector[i]);
        }
        result.trim();
        digits = result.digits;
    }

    void debug()const
    {
        DEBUG("begin");
        for(int i = nDigits() - 1; i >= 0; --i) DEBUG(digits[i]);
        DEBUG(isMinus);
        DEBUG("end");
    }

    static DIGIT findK(Number a, Number b)
    {//O(|a|), find k such that 0 <= k < BASE and kb <= a < (k + 1)b;
        DIGIT guess = a.digits.lastItem()/b.digits.lastItem();
        if(a.nDigits() > b.nDigits())
            guess = (a.digits.lastItem() * (1ull << BASE_RADIX) +
                a[a.nDigits() - 2])/b.digits.lastItem();
        while(mult(b, guess) > a) --guess;//executes <= 2 times
        return guess;
    }
    static Number divide(Number const& a, Number const& b1, Number& q)
    {//O(|a| * |b|)
        assert(!b1.isZero());
        q = 0;
        Number b = b1.abs(), r = a.abs();
        int norm = BASE_RADIX - lgFloor(b.digits.lastItem()) - 1;
        r <<= norm;
        b <<= norm;
        for(int i = r.nDigits() - b.nDigits(); i >= 0; --i)
        {
            int shift = i * BASE_RADIX;
            DIGIT k = findK(r >> shift, b);
            q += mult(Number(1) << shift, k);
            r -= mult(b << shift, k);
        }
        q.isMinus = r.isMinus = a.isMinus != b1.isMinus;
        return r >>= norm;
    }
    friend Number operator%(Number const& a, Number const& b)
    {
        Number quotient(0);
        return divide(a, b, quotient);
    }
    Number& operator/=(Number const&rhs){return *this = *this/rhs;}
    friend Number operator/(Number const& a, Number const& b)
    {
        Number quotient(0);
        divide(a, b, quotient);
        return quotient;
    }
    Number& operator%=(Number const&rhs){return *this = *this % rhs;}

    long long lg()const
        {return BASE_RADIX * (nDigits() - 1) + lgFloor(digits.lastItem());}

    Number sqrt()const
    {
        Number x(Number(1) << (1 + lg() / 2));
        for(;;)
        {
            Number y = (x + *this / x) / 2;
            if(y < x) x = y;
            else return x;
        }
    }
    Number power(Number const& p)const
    {
        Number x = *this, y = 1, n = p;
        for(;;)
        {
            if(n.isOdd()) y *= x;
            n >>= 1;
            if(n.isZero()) break;
            x *= x;
        }
        return y;
    }
    Number modPower(Number const& p, Number const& modulus)const
    {
        assert(!modulus.isZero());
        Number x = *this, y = 1, n = p;
        for(;;)
        {
            if(n.isOdd())
            {
                y *= x;
                y %= modulus;
            }
            n >>= 1;
            if(n.isZero()) break;
            x *= x;
            x %= modulus;
        }
        return y;
    }
};

Number extendedGcdR(Number const& a, Number const& b, Number& x, Number& y)
{
    if(!b.isPositive())
    {
        x = Number(1);
        y = Number(0);
        return a;
    }
    Number q, r = Number::divide(a, b, q), gcd = extendedGcdR(b, r, y, x);
    y -= q * x;
    return gcd;
}
Number extendedGcd(Number const& a, Number const& b, Number& x, Number& y)
{
    assert(a.isPositive() && b.isPositive());
    return a < b ? extendedGcdR(b, a, y, x) : extendedGcdR(a, b, x, y);
}
Number gcd(Number const& a, Number const& b)
{
    Number x, y;
    return extendedGcd(a, b, x, y);
}

Number modInverse(Number a, Number const& n)
{
    assert(a.isPositive() && a < n);
    Number x, y;
    extendedGcd(a, n, x, y);
    if(x.isNegative()) x += n;
    return x;
}

bool provenComposite(Number const& a, Number const& n)
{
    Number ONE = Number(1), oddPart = n - ONE;
    int nSquares = 0;
    while(oddPart.isEven())
    {
        oddPart >>= 1;
        ++nSquares;
    }
    Number x = a.modPower(oddPart, n);
    for(int i = 0; i < nSquares; ++i)
    {
        Number x2 = x.modPower(2, n);
        //if x2 is 1 x must have been 1 or -1 if n is prime
        if(x2 == ONE && x != ONE && x != n - ONE) return true;
        x = x2;
    }
    return x != ONE;
}

bool isPrime(Number const& n)
{
    if(n.isEven() || n < Number(2)) return false;
    int smallPrimes[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    for(int i = 0; i < sizeof(smallPrimes)/sizeof(int); ++i)
    {
        Number p = Number(smallPrimes[i]);
        if(n != p && (n % p).isZero()) return false;
    }
    //Miller-Rabin if trial division was inconclusive
    int nTrials = 1;
    int sizes[] = {73,105,132,198,223,242,253,265,335,480,543,627,747,927,
        1233,1854,4096}, nTests[] = {47,42,35,29,23,20,18,17,16,12,8,7,6,5,4,
        3,2};
    for(int i = 0; i < sizeof(sizes)/sizeof(*sizes); ++i)
        if(n.lg() < sizes[i])
        {
            nTrials = nTests[i];
            break;
        }
    while(nTrials--)
    {//use single digit exponents for efficiency
        Number::DIGIT_TYPE max = numeric_limits<Number::DIGIT_TYPE>::max();
        if(provenComposite(GlobalRNG.inRange(2, (n > Number(max) ?
            max : int(n[0])) - 1), n)) return false;
    }
    return true;
}

long long rationalize(double x, int& e)
{
    for(x = frexp(x, &e); x != (long long)x; x *= 2) --e;
    return x;
}

struct Rational
{
    Number numerator, denominator;
    Rational(Number const& theNumerator = Number(0), Number const&
        theDenominator = Number(1)): numerator(theNumerator),
        denominator(theDenominator)
    {
        assert(!denominator.isZero());
        reduce();
    }
    Rational(double x): denominator(1), numerator(1)
    {
        int e;
        long long n = rationalize(x, e);
        numerator = Number(n);
        if(e < 0) denominator <<= -e;
        else numerator <<= e;
    }
    void reduce()
    {
        Number g = gcd(numerator, denominator);
        numerator /= g;
        denominator /= g;
    }
    bool isZero()const{return numerator.isZero();}
    bool isMinus()const
        {return numerator.isNegative() != denominator.isNegative();}

    Rational operator-()const
    {
        Rational result = *this;
        result.numerator.negate();
        return result;
    }
    friend Rational operator+(Rational const& a, Rational const& b)
    {
        Rational result(a.numerator * b.denominator + b.numerator *
            a.denominator,a.denominator * b.denominator);
        result.reduce();
        return result;
    }
    friend Rational operator-(Rational const& a, Rational const& b)
        {return a + -b;}
    Rational& operator+=(Rational const&rhs){return *this = *this + rhs;}
    Rational& operator-=(Rational const&rhs){return *this = *this - rhs;}

    friend Rational operator*(Rational const& a, Rational const& b)
    {
        Rational result(a.numerator * b.numerator,
            a.denominator * b.denominator);
        result.reduce();
        return result;
    }
    friend Rational operator/(Rational const& a, Rational const& b)
    {
        assert(!b.isZero());
        Rational result(a.numerator * b.denominator,
            a.denominator * b.numerator);
        result.reduce();
        return result;
    }
    Rational& operator*=(Rational const&rhs){return *this = *this * rhs;}
    Rational& operator/=(Rational const&rhs){return *this = *this/rhs;}

    unsigned long long log()const{return numerator.lg() - denominator.lg();}
    Number evaluate(Number const& scale = Number(1))
        {return numerator * scale / denominator;}
};

}
#endif

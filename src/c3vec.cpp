#include "c3vec.hpp"

C3Vec& C3Vec::operator=(const C3Vec& rhs)
{
    if (this != &rhs) {
        c1 = rhs.c1;
        c2 = rhs.c2;
        c3 = rhs.c3;
    }
    return *this;
}
C3VecI& C3VecI::operator=(const C3VecI& rhs)
{
    if (this != &rhs) {
        c1 = rhs.c1;
        c2 = rhs.c2;
        c3 = rhs.c3;
    }
    return *this;
}
C3Vec& C3Vec::operator+=(const C3Vec& rhs)
{
    c1 += rhs.c1;
    c2 += rhs.c2;
    c3 += rhs.c3;
    return *this;
}

C3Vec& C3Vec::operator+=(const float& rhs)
{
    c1 += rhs;
    c2 += rhs;
    c3 += rhs;
    return *this;
}

C3Vec& C3Vec::operator-=(const C3Vec& rhs)
{
    c1 -= rhs.c1;
    c2 -= rhs.c2;
    c3 -= rhs.c3;
    return *this;
}

C3Vec& C3Vec::operator-=(const float& rhs)
{
    c1 -= rhs;
    c2 -= rhs;
    c3 -= rhs;
    return *this;
}
C3VecI& C3VecI::operator-=(const C3VecI& rhs)
{
    c1 -= rhs.c1;
    c2 -= rhs.c2;
    c3 -= rhs.c3;
    return *this;
}

C3VecI& C3VecI::operator-=(const float& rhs)
{
    c1 -= rhs;
    c2 -= rhs;
    c3 -= rhs;
    return *this;
}

C3Vec& C3Vec::operator*=(const C3Vec& rhs)
{
    c1 *= rhs.c1;
    c2 *= rhs.c2;
    c3 *= rhs.c3;
    return *this;
}

C3Vec& C3Vec::operator*=(const float& rhs)
{
    c1 *= rhs;
    c2 *= rhs;
    c3 *= rhs;
    return *this;
}

C3Vec& C3Vec::operator/=(const C3Vec& rhs)
{
    c1 /= rhs.c1;
    c2 /= rhs.c2;
    c3 /= rhs.c3;
    return *this;
}

C3Vec& C3Vec::operator/=(const float& rhs)
{
    c1 /= rhs;
    c2 /= rhs;
    c3 /= rhs;
    return *this;
}

C3VecI& C3VecI::operator/=(const C3VecI& rhs)
{
    c1 /= rhs.c1;
    c2 /= rhs.c2;
    c3 /= rhs.c3;
    return *this;
}

C3VecI& C3VecI::operator/=(const float& rhs)
{
    c1 /= rhs;
    c2 /= rhs;
    c3 /= rhs;
    return *this;
}
C3Vec C3Vec::operator+(const C3Vec& other)
{
    return C3Vec(this->c1 + other.c1, this->c2 + other.c2, this->c3 + other.c3);
}

C3Vec C3Vec::operator+(const float& other)
{
    return C3Vec(*this) += other;
}

C3Vec C3Vec::operator-(const C3Vec& other)
{
    return C3Vec(*this) -= other;
}

C3Vec C3Vec::operator-(const float& other)
{
    return C3Vec(*this) -= other;
}

C3VecI C3VecI::operator-(const C3VecI& other)
{
    return C3VecI(*this) -= other;
}

C3VecI C3VecI::operator-(const float& other)
{
    return C3VecI(*this) -= other;
}

C3Vec C3Vec::operator*(const C3Vec& other)
{
    return C3Vec(*this) *= other;
}

C3Vec C3Vec::operator*(const float& other)
{
    return C3Vec(*this) *= other;
}

C3Vec C3Vec::operator/(const C3Vec& other)
{
    return C3Vec(*this) /= other;
}

C3Vec C3Vec::operator/(const float& other)
{
    return C3Vec(*this) /= other;
}

C3VecI C3VecI::operator/(const C3VecI& other)
{
    return C3VecI(*this) /= other;
}

C3VecI C3VecI::operator/(const float& other)
{
    return C3VecI(*this) /= other;
}
// C3VecI

C3VecI& C3VecI::operator+=(const C3VecI& rhs)
{
    c1 += rhs.c1;
    c2 += rhs.c2;
    c3 += rhs.c3;
    return *this;
}

C3VecI& C3VecI::operator+=(const float& rhs)
{
    c1 += rhs;
    c2 += rhs;
    c3 += rhs;
    return *this;
}

C3VecI& C3VecI::operator*=(const C3VecI& rhs)
{
    c1 *= rhs.c1;
    c2 *= rhs.c2;
    c3 *= rhs.c3;
    return *this;
}

C3VecI& C3VecI::operator*=(const float& rhs)
{
    c1 *= rhs;
    c2 *= rhs;
    c3 *= rhs;
    return *this;
}
C3VecI C3VecI::operator+(const C3VecI& other)
{
    return C3VecI(this->c1 + other.c1, this->c2 + other.c2, this->c3 + other.c3);
}

C3VecI C3VecI::operator+(const float& other)
{
    return C3VecI(*this) += other;
}
C3VecI C3VecI::operator*(const C3VecI& other)
{
    return C3VecI(*this) *= other;
}

C3VecI C3VecI::operator*(const float& other)
{
    return C3VecI(*this) *= other;
}

// Global (not member) functions for lhs operators

C3Vec operator*(const float& other, const C3Vec& rhs)
{
    return C3Vec(rhs) *= other;
}

C3VecI operator*(const float& other, const C3VecI& rhs)
{
    return C3VecI(rhs) *= other;
}

C3VecI operator*(const std::complex<float>& other, const C3VecI& rhs)
{
    C3VecI tmp;
    tmp.c1 = other * rhs.c1;
    tmp.c2 = other * rhs.c2;
    tmp.c3 = other * rhs.c3;
    return tmp;
}

C3Vec operator+(const C3Vec& other, const C3Vec& rhs)
{
    return C3Vec(other.c1 + rhs.c1, other.c2 + rhs.c2, other.c3 + rhs.c3);
}
C3VecI operator+(const C3VecI& other, const C3VecI& rhs)
{
    return C3VecI(other.c1 + rhs.c1, other.c2 + rhs.c2, other.c3 + rhs.c3);
}

std::ostream& operator<<(std::ostream& os, const C3Vec& v)
{
    os << v.c1 << ", " << v.c2 << ", " << v.c3;
    return os;
}

std::ostream& operator<<(std::ostream& os, const C3VecI& v)
{
    os << v.c1 << ", " << v.c2 << ", " << v.c3;
    return os;
}

std::vector<C3Vec> operator-(const std::vector<C3Vec>& other, const C3Vec& rhs)
{
    std::vector<C3Vec> out(other.size());
    for (int i = 0; i < other.size(); i++) {
        out[i].c1 = other[i].c1 - rhs.c1;
        out[i].c2 = other[i].c2 - rhs.c2;
        out[i].c3 = other[i].c3 - rhs.c3;
    }
    return out;
}

std::vector<C3Vec> operator+(const std::vector<C3Vec>& other, const C3Vec& rhs)
{
    std::vector<C3Vec> out(other.size());
    for (int i = 0; i < other.size(); i++) {
        out[i].c1 = other[i].c1 + rhs.c1;
        out[i].c2 = other[i].c2 + rhs.c2;
        out[i].c3 = other[i].c3 + rhs.c3;
    }
    return out;
}
std::vector<C3VecI> operator+(const std::vector<C3VecI>& other, const C3VecI& rhs)
{
    std::vector<C3VecI> out(other.size());
    for (int i = 0; i < other.size(); i++) {
        out[i].c1 = other[i].c1 + rhs.c1;
        out[i].c2 = other[i].c2 + rhs.c2;
        out[i].c3 = other[i].c3 + rhs.c3;
    }
    return out;
}
std::vector<C3Vec> operator-(const std::vector<C3Vec>& other, const std::vector<C3Vec>& rhs)
{
    assert(other.size() == rhs.size());
    std::vector<C3Vec> out(other.size());
    for (int i = 0; i < other.size(); i++) {
        out[i].c1 = other[i].c1 - rhs[i].c1;
        out[i].c2 = other[i].c2 - rhs[i].c2;
        out[i].c3 = other[i].c3 - rhs[i].c3;
    }
    return out;
}

std::vector<C3Vec> operator+(const std::vector<C3Vec>& other, const std::vector<C3Vec>& rhs)
{
    assert(other.size() == rhs.size());
    std::vector<C3Vec> out(other.size());
    for (int i = 0; i < other.size(); i++) {
        out[i].c1 = other[i].c1 + rhs[i].c1;
        out[i].c2 = other[i].c2 + rhs[i].c2;
        out[i].c3 = other[i].c3 + rhs[i].c3;
    }
    return out;
}

std::vector<C3Vec> operator*(const std::vector<C3Vec>& other, const std::vector<float>& rhs)
{
    assert(other.size() == rhs.size());
    std::vector<C3Vec> out(other.size());
    for (int i = 0; i < other.size(); i++) {
        out[i].c1 = other[i].c1 * rhs[i];
        out[i].c2 = other[i].c2 * rhs[i];
        out[i].c3 = other[i].c3 * rhs[i];
    }
    return out;
}

float mag(const C3Vec& in)
{
    return sqrt(pow(in.c1, 2) + pow(in.c2, 2) + pow(in.c3, 2));
}

// This is not really a usefule magnitude
// and is only used to reduce a complex valued
// vector to a number for checking for Inf & NaNs
float mag(const C3VecI& in)
{
    float c1 = abs(in.c1);
    float c2 = abs(in.c2);
    float c3 = abs(in.c3);
    return sqrt(pow(c1, 2) + pow(c2, 2) + pow(c3, 2));
}

C3Vec pow(const C3Vec& in, const int arg)
{
    C3Vec out;
    out.c1 = pow(in.c1, arg);
    out.c2 = pow(in.c2, arg);
    out.c3 = pow(in.c3, arg);
    return out;
}

C3Vec sqrt(const C3Vec& in)
{
    C3Vec out;
    out.c1 = sqrt(in.c1);
    out.c2 = sqrt(in.c2);
    out.c3 = sqrt(in.c3);
    return out;
}

float dot(const C3Vec& Y, const C3Vec& X)
{
    return Y.c1 * X.c1 + Y.c2 * X.c2 + Y.c3 * X.c3;
}

std::complex<float> dot(const C3VecI& Y, const C3Vec& X)
{
    return Y.c1 * X.c1 + Y.c2 * X.c2 + Y.c3 * X.c3;
}

C3Vec atan2(const C3Vec& Y, const C3Vec& X)
{
    C3Vec out;
    out.c1 = atan2(Y.c1, X.c1);
    out.c2 = atan2(Y.c2, X.c2);
    out.c3 = atan2(Y.c3, X.c3);
    return out;
}

C3Vec cross(const C3Vec A, const C3Vec B)
{

    C3Vec answer;
    answer.c1 = (A.c2 * B.c3 - A.c3 * B.c2);
    answer.c2 = -(A.c1 * B.c3 - A.c3 * B.c1);
    answer.c3 = (A.c1 * B.c2 - A.c2 * B.c1);
    return answer;
}

int isnan(const C3Vec arg)
{
    int answer = 0;
    if (std::isnan(mag(arg)))
        answer = 1;
    return answer;
}

int isnan(const C3VecI arg)
{
    int answer = 0;
    if (std::isnan(mag(arg)))
        answer = 1;
    return answer;
}

int isinf(const C3Vec arg)
{
    int answer = 0;
    if (std::isinf(mag(arg)))
        answer = 1;
    return answer;
}

int isinf(const C3VecI arg)
{
    int answer = 0;
    if (std::isinf(mag(arg)))
        answer = 1;
    return answer;
}

float maxC3VecAbs(const std::vector<C3Vec>& input)
{

    std::vector<float> inputAbs(input.size());
    for (int i = 0; i < input.size(); i++) {
        inputAbs[i] = sqrt(pow(input[i].c1, 2) + pow(input[i].c2, 2) + pow(input[i].c3, 2));
    }
    return *std::max_element(inputAbs.begin(), inputAbs.end());
}

complex<float> intVecArray(const vector<float>& x, const vector<complex<float> >& f)
{

    complex<float> result;
    float h = x[1] - x[0];
    for (int i = 1; i < f.size(); i++) {
        result += h / 2.0f * (f[i - 1] + f[i]);
    }

    return result;
}

C3Vec intVecArray(const vector<float>& x, const vector<C3Vec>& f)
{

    C3Vec result;
    float h = x[1] - x[0];
    for (int i = 1; i < f.size(); i++) {
        result += h / 2.0 * (f[i - 1] + f[i]);
    }

    return result;
}

C3VecI intVecArray(const vector<float>& x, const vector<C3VecI>& f)
{

    C3VecI result;
    float h = x[1] - x[0];
    for (int i = 1; i < f.size(); i++) {
        result += h / 2.0 * (f[i - 1] + f[i]);
    }

    return result;
}

void kj_print(const C3Vec arg, string name)
{
    cout << name << ".c1: " << arg.c1 << "  " << name << ".c2: " << arg.c2 << "  " << name << ".c3: " << arg.c3 << endl;
    return;
}

void kj_print(const float arg, string name)
{
    cout << name << ": " << arg << endl;
    return;
}

C3Vec XYZ_to_CYL(const C3Vec xyz)
{
    C3Vec cyl;
    cyl.c1 = sqrt(pow(xyz.c1, 2) + pow(xyz.c2, 2));
    cyl.c2 = atan2(xyz.c2, xyz.c1);
    cyl.c3 = xyz.c3;
    return cyl;
}

C3Vec CYL_to_XYZ(const C3Vec cyl)
{
    C3Vec xyz;
    xyz.c1 = cyl.c1 * cos(cyl.c2);
    xyz.c2 = cyl.c1 * sin(cyl.c2);
    xyz.c3 = cyl.c3;
    return xyz;
}

C3Vec operator*(const float A[][3], const C3Vec x)
{
    C3Vec B;
    B.c1 = A[0][0] * x.c1 + A[0][1] * x.c2 + A[0][2] * x.c3;
    B.c2 = A[1][0] * x.c1 + A[1][1] * x.c2 + A[1][2] * x.c3;
    B.c3 = A[2][0] * x.c1 + A[2][1] * x.c2 + A[2][2] * x.c3;
    return B;
}

C3VecI operator*(const float A[][3], const C3VecI x)
{
    C3VecI B;
    B.c1 = A[0][0] * x.c1 + A[0][1] * x.c2 + A[0][2] * x.c3;
    B.c2 = A[1][0] * x.c1 + A[1][1] * x.c2 + A[1][2] * x.c3;
    B.c3 = A[2][0] * x.c1 + A[2][1] * x.c2 + A[2][2] * x.c3;
    return B;
}

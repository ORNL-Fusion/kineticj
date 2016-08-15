
PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator-=(const C3<T>& rhs)
{
    c1 -= rhs.c1;
    c2 -= rhs.c2;
    c3 -= rhs.c3;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator-=(const float& rhs)
{
    c1 -= rhs;
    c2 -= rhs;
    c3 -= rhs;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator/=(const C3<T>& rhs)
{
    c1 /= rhs.c1;
    c2 /= rhs.c2;
    c3 /= rhs.c3;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator/=(const float& rhs)
{
    c1 /= rhs;
    c2 /= rhs;
    c3 /= rhs;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator-(const C3<T>& other)
{
    return C3<T>(*this) -= other;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator-(const float& other)
{
    return C3<T>(*this) -= other;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator/(const C3<T>& other)
{
    return C3<T>(*this) /= other;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator/(const float& other)
{
    return C3<T>(*this) /= other;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator+=(const C3<T>& rhs)
{
    c1 += rhs.c1;
    c2 += rhs.c2;
    c3 += rhs.c3;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator+=(const float& rhs)
{
    c1 += rhs;
    c2 += rhs;
    c3 += rhs;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator*=(const C3<T>& rhs)
{
    c1 *= rhs.c1;
    c2 *= rhs.c2;
    c3 *= rhs.c3;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T>& C3<T>::operator*=(const float& rhs)
{
    c1 *= rhs;
    c2 *= rhs;
    c3 *= rhs;
    return *this;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator+(const C3<T>& other)
{
    return C3<T>(this->c1 + other.c1, this->c2 + other.c2, this->c3 + other.c3);
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator+(const float& other)
{
    return C3<T>(*this) += other;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator*(const C3<T>& other)
{
    return C3<T>(*this) *= other;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> C3<T>::operator*(const float& other)
{
    return C3<T>(*this) *= other;
}

// Global (not member) functions for lhs operators

PRAGMA
template <typename T>
HOST DEVICE
C3<T> operator*(const float& other, const C3<T>& rhs)
{
    return C3<T>(rhs) *= other;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> operator*(const std::complex<float>& other, const C3<T>& rhs)
{
    C3<T> tmp;
    tmp.c1 = other * rhs.c1;
    tmp.c2 = other * rhs.c2;
    tmp.c3 = other * rhs.c3;
    return tmp;
}

#ifdef __CUDACC__
PRAGMA
inline
HOST DEVICE
C3<thrust::complex<float> > operator*(const thrust::complex<float>& other, const C3<thrust::complex<float> >& rhs)
{
    C3<thrust::complex<float> > tmp;
    tmp.c1 = other * rhs.c1;
    tmp.c2 = other * rhs.c2;
    tmp.c3 = other * rhs.c3;
    return tmp;
}
#endif

PRAGMA
template <typename T>
HOST DEVICE
C3<T> operator+(const C3<T>& other, const C3<T>& rhs)
{
    return C3<T>(other.c1 + rhs.c1, other.c2 + rhs.c2, other.c3 + rhs.c3);
}

PRAGMA
template <typename T>
HOST
std::ostream& operator<<(std::ostream& os, const C3<T>& v)
{
    os << v.c1 << ", " << v.c2 << ", " << v.c3;
    return os;
}

PRAGMA
template <typename T>
HOST DEVICE
std::vector<C3<T> > operator+(const std::vector<C3<T> >& other, const C3<T>& rhs)
{
    std::vector<C3<T> > out(other.size());
    for (int i = 0; i < other.size(); i++) {
        out[i].c1 = other[i].c1 + rhs.c1;
        out[i].c2 = other[i].c2 + rhs.c2;
        out[i].c3 = other[i].c3 + rhs.c3;
    }
    return out;
}

PRAGMA
template <typename T>
HOST DEVICE
T mag(const C3<T>& in)
{
    //float c1 = std::abs(in.c1);
    //float c2 = std::abs(in.c2);
    //float c3 = std::abs(in.c3);
    //return sqrt(pow(c1, 2) + pow(c2, 2) + pow(c3, 2));
    T c1s = T(std::pow(in.c1, 2));
    T c2s = T(std::pow(in.c2, 2));
    T c3s = T(std::pow(in.c3, 2));

    return std::sqrt(c1s + c2s + c3s);
}

template <typename T>
HOST DEVICE
T dot(const C3<T>& Y, const C3<T>& X)
{
    return Y.c1 * X.c1 + Y.c2 * X.c2 + Y.c3 * X.c3;
}

inline
HOST 
std::complex<float> dot(const C3<std::complex<float> >& Y, const C3<float>& X)
{
    return Y.c1 * X.c1 + Y.c2 * X.c2 + Y.c3 * X.c3;
}

#ifdef __CUDACC__
inline
HOST DEVICE 
thrust::complex<float> dot(const C3<thrust::complex<float> >& Y, const C3<float>& X)
{
    return Y.c1 * X.c1 + Y.c2 * X.c2 + Y.c3 * X.c3;
}
#endif 

PRAGMA
template <typename T>
HOST DEVICE
C3<T> cross ( const C3<T> A, const C3<T> B ) {

        C3<T> answer;
        answer.c1 =  (A.c2*B.c3 - A.c3*B.c2);
        answer.c2 = -(A.c1*B.c3 - A.c3*B.c1);
        answer.c3 =  (A.c1*B.c2 - A.c2*B.c1);
        return answer;
}

PRAGMA
template <typename T, typename T2>
HOST DEVICE
C3<T2> cross(const C3<T> A, const C3<T2> B)
{

    C3<T2> answer;
    answer.c1 = (A.c2 * B.c3 - A.c3 * B.c2);
    answer.c2 = -(A.c1 * B.c3 - A.c3 * B.c1);
    answer.c3 = (A.c1 * B.c2 - A.c2 * B.c1);
    return answer;
}

PRAGMA
template <typename T>
HOST DEVICE
int isnan(const C3<T> arg)
{
    int answer = 0;
    if (std::isnan(mag(arg)))
        answer = 1;
    return answer;
}

PRAGMA
template <typename T>
HOST DEVICE
int isinf(const C3<T> arg)
{
    int answer = 0;
    if (std::isinf(mag(arg)))
        answer = 1;
    return answer;
}

PRAGMA
template <typename T>
HOST DEVICE
float maxC3VecAbs(const std::vector<C3<T> >& input)
{

    std::vector<float> inputAbs(input.size());
    for (int i = 0; i < input.size(); i++) {
        inputAbs[i] = sqrt(pow(input[i].c1, 2) + pow(input[i].c2, 2) + pow(input[i].c3, 2));
    }
    return *std::max_element(inputAbs.begin(), inputAbs.end());
}

PRAGMA
HOST DEVICE
inline
std::complex<float> intVecArray(const std::vector<float>& x, const std::vector<std::complex<float> >& f)
{

    std::complex<float> result;
    float h = x[1] - x[0];
    for (int i = 1; i < f.size(); i++) {
        result += h / 2.0f * (f[i - 1] + f[i]);
#ifndef __CUDACC__
#if DEBUG_INTVECARRAY > 0
        std::cout << "result: " << result << std::endl;
        std::cout << "f[i]: " << f[i] << std::endl;
#endif
#endif
    }

    return result;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> intVecArray(const std::vector<float>& x, const std::vector<C3<T> >& f)
{

    C3<T> result;
    float h = x[1] - x[0];
    for (int i = 1; i < f.size(); i++) {
        result += h / 2.0 * (f[i - 1] + f[i]);
    }

    return result;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> XYZ_to_CYL(const C3<T> xyz)
{
    C3<T> cyl;
    cyl.c1 = sqrt(pow(xyz.c1, 2) + pow(xyz.c2, 2));
    cyl.c2 = atan2(xyz.c2, xyz.c1);
    cyl.c3 = xyz.c3;
    return cyl;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> CYL_to_XYZ(const C3<T> cyl)
{
    C3<T> xyz;
    xyz.c1 = cyl.c1 * cos(cyl.c2);
    xyz.c2 = cyl.c1 * sin(cyl.c2);
    xyz.c3 = cyl.c3;
    return xyz;
}

PRAGMA
template <typename T>
HOST DEVICE
C3<T> operator*(const float A[][3], const C3<T> x)
{
    C3<T> B;
    B.c1 = A[0][0] * x.c1 + A[0][1] * x.c2 + A[0][2] * x.c3;
    B.c2 = A[1][0] * x.c1 + A[1][1] * x.c2 + A[1][2] * x.c3;
    B.c3 = A[2][0] * x.c1 + A[2][1] * x.c2 + A[2][2] * x.c3;
    return B;
}


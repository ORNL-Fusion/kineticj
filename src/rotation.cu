#include "rotation.hpp"
#include "constants.hpp"

void transpose(float A[][3])
{

    float B[3][3];

    B[0][0] = A[0][0];
    B[1][0] = A[0][1];
    B[2][0] = A[0][2];

    B[0][1] = A[1][0];
    B[1][1] = A[1][1];
    B[2][1] = A[1][2];

    B[0][2] = A[2][0];
    B[1][2] = A[2][1];
    B[2][2] = A[2][2];

    A = B;
}

C3<float> rot_XYZ_to_abp(const C3<float> A_XYZ, const C3<float> bUnit_XYZ, const int direction)
{

    // If direction<1 then the inverse rotation is applied, i.e., abp_to_XYZ

    C3<float> A_abp;

    C3<float> xu_xyz(1, 0, 0);
    C3<float> yu_xyz(0, 1, 0);
    C3<float> zu_xyz(0, 0, 1);

    C3<float> pu_xyz = bUnit_XYZ;

    // alp is mostly in the +/- x / r direction depending on B toroidal direction
    // bet is mostly z direction

    C3<float> a_xyz = cross(zu_xyz, pu_xyz);
    C3<float> au_xyz = a_xyz / mag(a_xyz);
    //C3<float> au_xyz = ( static_cast<float>(1.0) / mag(a_xyz) ) * a_xyz;

    C3<float> b_xyz = cross(pu_xyz, au_xyz);
    C3<float> bu_xyz = b_xyz / mag(b_xyz);

#if DEBUG_ROTATION >= 1
    
    std::cout<< "bUnit_XYZ: "<<bUnit_XYZ<<std::endl;
    std::cout<< "a_xyz: "<<a_xyz<<std::endl;
    std::cout<< "mag(a_xyz): " << mag(a_xyz) << std::endl;

    C3<float> au_xyz2 = au_xyz;
    C3<float> bu_xyz2 = bu_xyz;
    C3<float> pu_xyz2 = pu_xyz;

    std::cout << "au_xyz: " << au_xyz << std::endl;
    std::cout << "bu_xyz: " << bu_xyz << std::endl;
    std::cout << "pu_xyz: " << pu_xyz << std::endl;
#endif

    // Rotation 1

    float theta = acos(dot(xu_xyz, au_xyz));

#if DEBUG_ROTATION >= 1
    std::cout << "xu_xyz: " << xu_xyz << std::endl;
    std::cout << "au_xyz: " << au_xyz << std::endl;
    std::cout << "dot: " << dot(xu_xyz, au_xyz) << std::endl;

    std::cout << "theta [rad]: " << theta << std::endl;
    std::cout << "theta [deg]: " << theta * 180.0 / physConstants::pi << std::endl;
#endif

    float q0 = cos(theta / 2.0);
    float q1 = sin(theta / 2.0) * (-zu_xyz.c1);
    float q2 = sin(theta / 2.0) * (-zu_xyz.c2);
    float q3 = sin(theta / 2.0) * (-zu_xyz.c3);

    // Construct the rotation matrix

    float rot1[3][3];

    rot1[0][0] = pow(q0, 2.0) + pow(q1, 2.0) - pow(q2, 2.0) - pow(q3, 2.0);
    rot1[0][1] = 2 * (q1 * q2 - q0 * q3);
    rot1[0][2] = 2 * (q1 * q3 + q0 * q2);
    rot1[1][0] = 2 * (q2 * q1 + q0 * q3);
    rot1[1][1] = pow(q0, 2.0) - pow(q1, 2.0) + pow(q2, 2.0) - pow(q3, 2.0);
    rot1[1][2] = 2 * (q2 * q3 - q0 * q1);
    rot1[2][0] = 2 * (q3 * q1 - q0 * q2);
    rot1[2][1] = 2 * (q3 * q2 + q0 * q1);
    rot1[2][2] = pow(q0, 2.0) - pow(q1, 2.0) - pow(q2, 2.0) + pow(q3, 2.0);

    if (direction < 0) {
        transpose(rot1);
    }

    au_xyz = rot1 * au_xyz;
    bu_xyz = rot1 * bu_xyz;
    pu_xyz = rot1 * pu_xyz;

#if DEBUG_ROTATION >= 1
    std::cout << "au_rtz 1: " << au_xyz << std::endl;
    std::cout << "bu_rtz 1: " << bu_xyz << std::endl;
    std::cout << "pu_rtz 1: " << pu_xyz << std::endl;
#endif

    // Rotation 2

    theta = acos(dot(zu_xyz, pu_xyz));

#if DEBUG_ROTATION >= 1
    std::cout << "theta: " << theta * 180.0 / physConstants::pi << std::endl;
#endif
    q0 = cos(theta / 2.0);
    q1 = sin(theta / 2.0) * (-xu_xyz.c1);
    q2 = sin(theta / 2.0) * (-xu_xyz.c2);
    q3 = sin(theta / 2.0) * (-xu_xyz.c3);

    // Construct the rotation matrix

    float rot2[3][3];

    rot2[0][0] = pow(q0, 2.0) + pow(q1, 2.0) - pow(q2, 2.0) - pow(q3, 2.0);
    rot2[0][1] = 2 * (q1 * q2 - q0 * q3);
    rot2[0][2] = 2 * (q1 * q3 + q0 * q2);
    rot2[1][0] = 2 * (q2 * q1 + q0 * q3);
    rot2[1][1] = pow(q0, 2.0) - pow(q1, 2.0) + pow(q2, 2.0) - pow(q3, 2.0);
    rot2[1][2] = 2 * (q2 * q3 - q0 * q1);
    rot2[2][0] = 2 * (q3 * q1 - q0 * q2);
    rot2[2][1] = 2 * (q3 * q2 + q0 * q1);
    rot2[2][2] = pow(q0, 2.0) - pow(q1, 2.0) - pow(q2, 2.0) + pow(q3, 2.0);

    if (direction < 0) {
        transpose(rot2);
    }

    au_xyz = rot2 * au_xyz;
    bu_xyz = rot2 * bu_xyz;
    pu_xyz = rot2 * pu_xyz;

#if DEBUG_ROTATION >= 1
    std::cout << "au_xyz 2: " << au_xyz << std::endl;
    std::cout << "bu_xyz 2: " << bu_xyz << std::endl;
    std::cout << "pu_xyz 2: " << pu_xyz << std::endl;
#endif

    A_abp = rot2 * (rot1 * A_XYZ);

    return A_abp;
}

C3<float> rot_axis_angle(const C3<float> v, const C3<float> u, const float th_deg) {

        // See https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions

        C3<float> result;
        float th = th_deg * physConstants::pi / 180.0;
        float cosTh = cos(th);
        float sinTh = sin(th);

        float ux = u.c1;
        float uy = u.c2;
        float uz = u.c3;

        float R11 = cosTh + std::pow(ux,2)*(1-cosTh);
        float R12 = ux*uy*(1-cosTh)-uz*sinTh;
        float R13 = ux*uz*(1-cosTh)+uy*sinTh;

        float R21 = uy*ux*(1-cosTh)+uz*sinTh;
        float R22 = cosTh+std::pow(uy,2)*(1-cosTh);
        float R23 = uy*uz*(1-cosTh)-ux*sinTh;

        float R31 = uz*ux*(1-cosTh)-uy*sinTh;
        float R32 = uz*uy*(1-cosTh)+ux*sinTh;
        float R33 = cosTh+std::pow(uz,2)*(1-cosTh);

        result.c1 = R11 * v.c1 + R12 * v.c2 + R13 * v.c3;
        result.c2 = R21 * v.c1 + R22 * v.c2 + R23 * v.c3;
        result.c3 = R31 * v.c1 + R32 * v.c2 + R33 * v.c3;

        return result;
}

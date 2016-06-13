#include "rotation.hpp"

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

C3Vec rot_XYZ_to_abp(const C3Vec A_XYZ, const C3Vec bUnit_XYZ, const int direction)
{

    // If direction<1 then the inverse rotation is applied, i.e., abp_to_XYZ

    C3Vec A_abp;

    C3Vec xu_xyz(1, 0, 0);
    C3Vec yu_xyz(0, 1, 0);
    C3Vec zu_xyz(0, 0, 1);

    C3Vec pu_xyz = bUnit_XYZ;

    // alp is mostly in the +/- x / r direction depending on B toroidal direction
    // bet is mostly z direction

    C3Vec a_xyz = cross(zu_xyz, pu_xyz);
    C3Vec au_xyz = a_xyz / mag(a_xyz);

    C3Vec b_xyz = cross(pu_xyz, au_xyz);
    C3Vec bu_xyz = b_xyz / mag(b_xyz);

#if DEBUG_ROTATION >= 1
    C3Vec au_xyz2 = au_xyz;
    C3Vec bu_xyz2 = bu_xyz;
    C3Vec pu_xyz2 = pu_xyz;

    cout << "au_xyz: " << au_xyz << endl;
    cout << "bu_xyz: " << bu_xyz << endl;
    cout << "pu_xyz: " << pu_xyz << endl;
#endif

    // Rotation 1

    float theta = acos(dot(xu_xyz, au_xyz));

#if DEBUG_ROTATION >= 1
    cout << "theta: " << theta * 180.0 / _pi << endl;
#endif

    float q0 = cos(theta / 2.0);
    float q1 = sin(theta / 2.0) * (-zu_xyz.c1);
    float q2 = sin(theta / 2.0) * (-zu_xyz.c2);
    float q3 = sin(theta / 2.0) * (-zu_xyz.c3);

    // Construct the rotation matrix

    float rot1[3][3];

    rot1[0][0] = pow(q0, 2) + pow(q1, 2) - pow(q2, 2) - pow(q3, 2);
    rot1[0][1] = 2 * (q1 * q2 - q0 * q3);
    rot1[0][2] = 2 * (q1 * q3 + q0 * q2);
    rot1[1][0] = 2 * (q2 * q1 + q0 * q3);
    rot1[1][1] = pow(q0, 2) - pow(q1, 2) + pow(q2, 2) - pow(q3, 2);
    rot1[1][2] = 2 * (q2 * q3 - q0 * q1);
    rot1[2][0] = 2 * (q3 * q1 - q0 * q2);
    rot1[2][1] = 2 * (q3 * q2 + q0 * q1);
    rot1[2][2] = pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2);

    if (direction < 0) {
        transpose(rot1);
    }

    au_xyz = rot1 * au_xyz;
    bu_xyz = rot1 * bu_xyz;
    pu_xyz = rot1 * pu_xyz;

#if DEBUG_ROTATION >= 1
    cout << "au_rtz 1: " << au_xyz << endl;
    cout << "bu_rtz 1: " << bu_xyz << endl;
    cout << "pu_rtz 1: " << pu_xyz << endl;
#endif

    // Rotation 2

    theta = acos(dot(zu_xyz, pu_xyz));

#if DEBUG_ROTATION >= 1
    cout << "theta: " << theta * 180.0 / _pi << endl;
#endif
    q0 = cos(theta / 2.0);
    q1 = sin(theta / 2.0) * (-xu_xyz.c1);
    q2 = sin(theta / 2.0) * (-xu_xyz.c2);
    q3 = sin(theta / 2.0) * (-xu_xyz.c3);

    // Construct the rotation matrix

    float rot2[3][3];

    rot2[0][0] = pow(q0, 2) + pow(q1, 2) - pow(q2, 2) - pow(q3, 2);
    rot2[0][1] = 2 * (q1 * q2 - q0 * q3);
    rot2[0][2] = 2 * (q1 * q3 + q0 * q2);
    rot2[1][0] = 2 * (q2 * q1 + q0 * q3);
    rot2[1][1] = pow(q0, 2) - pow(q1, 2) + pow(q2, 2) - pow(q3, 2);
    rot2[1][2] = 2 * (q2 * q3 - q0 * q1);
    rot2[2][0] = 2 * (q3 * q1 - q0 * q2);
    rot2[2][1] = 2 * (q3 * q2 + q0 * q1);
    rot2[2][2] = pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2);

    if (direction < 0) {
        transpose(rot2);
    }

    au_xyz = rot2 * au_xyz;
    bu_xyz = rot2 * bu_xyz;
    pu_xyz = rot2 * pu_xyz;

#if DEBUG_ROTATION >= 1
    cout << "au_xyz 2: " << au_xyz << endl;
    cout << "bu_xyz 2: " << bu_xyz << endl;
    cout << "pu_xyz 2: " << pu_xyz << endl;
#endif

    A_abp = rot2 * (rot1 * A_XYZ);

    return A_abp;
}

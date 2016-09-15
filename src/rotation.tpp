PRAGMA
template <typename T>
HOST DEVICE
T rot_CYL_to_XYZ ( const float t, const T vec, const int direction ) {

    // t here is the the cylindrical angle position in rtz (radians)        

    float rot[3][3];

    rot[0][0] =  cos(t);
    rot[0][1] = -sin(t);
    rot[0][2] = 0;

    rot[1][0] =  sin(t);
    rot[1][1] =  cos(t);
    rot[1][2] = 0;

    rot[2][0] = 0;
    rot[2][1] = 0;
    rot[2][2] = 1;

    if(direction<0) {
        transpose(rot);
    }
    
    return rot * vec;
}



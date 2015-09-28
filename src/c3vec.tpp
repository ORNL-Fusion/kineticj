template <typename C3T1, typename C3T2>
C3VecI cross ( const C3T1 A, const C3T2 B ) {

        C3VecI answer;
        answer.c1 =  (A.c2*B.c3 - A.c3*B.c2);
        answer.c2 = -(A.c1*B.c3 - A.c3*B.c1);
        answer.c3 =  (A.c1*B.c2 - A.c2*B.c1);
        return answer;
}

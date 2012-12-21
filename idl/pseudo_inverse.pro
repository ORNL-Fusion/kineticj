FUNCTION pseudo_inverse,A,U,V,S,TOL=tol,SORT=sort,_extra=extra
;+
;B=pseudo_inverse(A,[U,V,S,tol,EXTRA=extra)
;computes the pseudo-inverse using SVDC.
;B=pseudo_inverse(A,U,V,S) returns the SVDC values
;U,V and S such that A=U##S##Transpose(V). This is
;especially useful in examining the sensitivity of
;the solution. Small values along the diagonal of
;S can be zeroed out in computing the pseudo-inverse
;to reduce numerical sensitivity.
;
;PARAMETERS
;Set TOL to a threshold value on the singular values
;to be accepted in computing the inverse. Example
;Ap=pseudo_inverse(A,TOL=1e-5) causes singular values
;smaller than 1e-5 to be excluded from the inverse
;calculation. Default: tol=1e-10
;
;Set SORT to sort the singular values in descending
;order before computing the inverse.
;
;2004-05-13 HR Added SORT and TOL keywords and replaced
;use of the special DIAG() function with DIAG_MATRIX()
;
;Reference:
;Moon and Stirling,
;Mathematical Methods and Algorithms
;Prentice-Hall, 2000
;Section 7.4
;
;H. Rhody  December 30, 2003
;-
type=size(A,/TYPE)
dim=size(A,/DIM)

IF n_elements(dim) EQ 1 THEN BEGIN
	;Case of a row vector
	;Since SVDC does not handle 1D arrays, we
	;can coerce a solution by making a column vector
	;and extracting a solution from it.
	SVDC,transpose(A),S,U1,V,_EXTRA=extra
	U=V
	V=U1
ENDIF ELSE SVDC,A,S,U,V,_extra=extra

IF N_ELEMENTS(tol) LE 0 THEN tol=1e-10

IF KEYWORD_SET(sort) THEN BEGIN
	I=REVERSE(SORT(S))
	S=S[I]
	U=U[I,*]
	V=V[I,*]
END

W=MAKE_ARRAY(SIZE=SIZE(S),VALUE=0)
I=WHERE(ABS(S) GE tol)
IF MIN(I) LT 0 THEN MESSAGE,'No singular values larger than TOL'

W[I]=1.0/S[I]

W=diag_matrix(W)

if type LE 5 then $
	RETURN,V##W##Transpose(u) ELSE $
		Return,V##W##transpose(conj(U))
END

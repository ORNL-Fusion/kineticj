pro kj_test_mpe

common Cg_comm1, A

; Generate A

n = 100

A = dblArr(n,n)

for i=0,n-1 do begin

	if i eq 0 then a[0:3,i] = [5,2,1,1]
	if i eq 1 then a[0:4,i] = [2,6,3,1,1]
	if i eq 2 then a[0:5,i] = [1,3,6,3,1,1]
	if i gt 2 and i lt n-3 then a[i-3:i+3,i] = [1,1,3,6,3,1,1]
	if i eq n-3 then a[n-1-5:n-1,i] = [1,1,3,6,3,1]
	if i eq n-2 then a[n-1-4:n-1,i] = [1,1,3,6,2]
	if i eq n-1 then a[n-1-3:n-1,i] = [1,1,2,5]

endfor

A = 0.06d0 * A

s_exact = dblArr(n)+1

b = s_exact - A ## s_exact

x = transpose(dblArr(n))

C = dblArr(n,n)
for i=0,n-1 do begin
	for j=0,n-1 do begin
		if i eq j then C[i,j] = 1
	endfor
endfor

C = I - A

eigenVals = la_eigenQL(A)

print, min(eigenVals)
print, max(eigenVals)

nk = 10 
nm = 5 
w = 2d0



xAll = !null
s_mpe = x

p=plot(x)
for m = 0, nm-1 do begin

	for k = 0, nk-1 do begin
	
		xAll = [[xAll],[x]]
		x0 = x
		Fx = A # x0 + b
		x = (1d0-w)*x0 + w*Fx
	
		s1=norm(x-x0,LNORM=2)
		s2=norm(x-s_exact,LNORM=2)
	
		;print, k, s1,s2
		;!null=plot(x,/over)
	endfor

	s0 = s_mpe	
	s_mpe = kj_mpe(xAll)

	!null = plot(s_mpe, thick=2,color='b',/over)

	s1=norm(s_mpe-s0,LNORM=2)
	s2=norm(s_mpe-s_exact,LNORM=2)
	
	print, m, s1,s2
	xAll = s_mpe
	x = s_mpe
;stop
endfor


stop

end

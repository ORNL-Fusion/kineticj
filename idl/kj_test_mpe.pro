function amultp, p
	common Cg_comm1, A
	return, A # p
end

pro kj_test_mpe

common Cg_comm1, A

; Generate A

n = 1000

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

b = s_exact - A # s_exact

x = transpose(dblArr(n))

C = dblArr(n,n)
for i=0,n-1 do begin
	for j=0,n-1 do begin
		if i eq j then C[i,j] = 1
	endfor
endfor

C = I - A

eigenVals = la_eigenQL(A)


nIt = 10
w = 1d0

;result=imsl_sp_cg('amultp',transpose(b))

p=plot(x)

for i = 0, nIt-1 do begin

	x0 = x
	Fx = A ## x + b
	x = (1d0-w)*x + w*Fx

	s = sqrt(mean((x-x0)^2))
	print, s
	!null=plot(x,/over)
endfor


stop

end

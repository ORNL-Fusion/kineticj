pro kj_sigma_vs_t

@constants

n = 100
tMin = 0.1
tMax = 10e3
T_eV = 10d0^(findGen(n)/(n-1)*(alog10(tMax)-alog10(tMin))+alog10(tMin)) 

f = 13d6
Z = +1d0
amu =  2.0;_me_amu
B = 1d0
density = 2d19
harmonicNumber = 1
kPar = 100  
kPer = 10 

eps = ComplexArr(3,3,n)
eps_cold = ComplexArr(3,3,n)

for i=0,n-1 do begin

    thisTeV = T_eV[i]
    this_eps = kj_hot_epsilon(f, amu, Z, B, density, harmonicNumber, kPar, kPer, thisTeV, $
                    epsilon_cold = this_eps_cold )
    
    eps[*,*,i] = this_eps 
    eps_cold[*,*,i] = this_eps_cold 

endfor

layout=[3,3]
pos = 1
thick = 2 
style = '--'
transparency = 50

p=plot(T_eV,eps[0,0,*],layout=[[layout],pos],title='eps[0,0]',/xlog)
p=plot(T_eV,imaginary(eps[0,0,*]),color='r',/over)
p=plot(T_eV,eps_cold[0,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[0,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,eps[0,1,*],layout=[[layout],pos],title='eps[0,1]',/current,/xlog)
p=plot(T_eV,imaginary(eps[0,1,*]),color='r',/over)
p=plot(T_eV,eps_cold[0,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[0,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,eps[0,2,*],layout=[[layout],pos],title='eps[0,2]',/current,/xlog)
p=plot(T_eV,imaginary(eps[0,2,*]),color='r',/over)
p=plot(T_eV,eps_cold[0,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[0,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos
p=plot(T_eV,eps[1,0,*],layout=[[layout],pos],title='eps[1,0]',/current,/xlog)
p=plot(T_eV,imaginary(eps[1,0,*]),color='r',/over)
p=plot(T_eV,eps_cold[1,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[1,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,eps[1,1,*],layout=[[layout],pos],title='eps[1,1]',/current,/xlog)
p=plot(T_eV,imaginary(eps[1,1,*]),color='r',/over)
p=plot(T_eV,eps_cold[1,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[1,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,eps[1,2,*],layout=[[layout],pos],title='eps[1,2]',/current,/xlog)
p=plot(T_eV,imaginary(eps[1,2,*]),color='r',/over)
p=plot(T_eV,eps_cold[1,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[1,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos
p=plot(T_eV,eps[2,0,*],layout=[[layout],pos],title='eps[2,0]',/current,/xlog)
p=plot(T_eV,imaginary(eps[2,0,*]),color='r',/over)
p=plot(T_eV,eps_cold[2,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[2,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,eps[2,1,*],layout=[[layout],pos],title='eps[2,1]',/current,/xlog)
p=plot(T_eV,imaginary(eps[2,1,*]),color='r',/over)
p=plot(T_eV,eps_cold[2,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[2,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,eps[2,2,*],layout=[[layout],pos],title='eps[2,2]',/current,/xlog)
p=plot(T_eV,imaginary(eps[2,2,*]),color='r',/over)
p=plot(T_eV,eps_cold[2,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(eps_cold[2,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

; Create kj input to run this case

kj_create_single_k_input, b0=B, kPar=kPar, kPer=kPer, f_Hz=f

stop

end

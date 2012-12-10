pro kj_sigma_of_energy, hopper=hopper

sig33All1 = !null
sig33All2 = !null

energyAll = !null

ee = [100.0,200.0,500.0,800.0,1000.0,1500.0,2000.0,3000.0,4000.0,5000.0,6e3,7e3,1d4]
;nPts = 30.0
;ee = 10^((fIndGen(nPts)/(nPts/2.0-1)))*1e2
;ee = ee[0:-2]

stop
for i = 0,n_elements(ee)-1 do begin

	create_test_particle_f, /weighted_maxwellian_xyz, rsfwc_1d='data/kj_single_1d.nc', energy_keV = ee[i]*1e-3
	;create_test_particle_f, standard_maxwellian_1d=1, rsfwc_1d='data/kj_single_1d.nc', energy_keV = ee[i]*1e-3

    if keyword_set(hopper) then begin
	    spawn, 'aprun -n 1 -d 24 /global/homes/g/greendl1/kineticj/bin/kineticj'
    endif else begin
	    spawn, '~/code/kineticj/bin/kineticj'
    endelse

	kj_plot_current, sig33 = sig33

	sig33All1 = [sig33All1,sig33[0]]
	sig33All2 = [sig33All2,sig33[1]]

	energyAll = [energyAll,ee[i]]

	print, 'Energy [eV]: ', EE[i], ' sig33: ', sig33[0]

endfor

;sig33All1 = conj(sig33All1)
;sig33All2 = conj(sig33All2)

p=plot(energyAll,sig33All1,/xlog,thick=2.0,transparency=50,color='b',buffer=1)
!null=plot(energyAll,imaginary(sig33All1),/over,thick=2.0,transparency=50)

!null=plot(energyAll,imaginary(sig33All2),/over,color='b',thick=2.0,transparency=50)
!null=plot(energyAll,imaginary(sig33All2),/over,thick=2.0,transparency=50)
p.save, 'sigma_of_energy.png'
stop
end

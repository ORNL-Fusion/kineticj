pro kj_sigma_of_energy

sig33All1 = !null
sig33All2 = !null

energyAll = !null

for EE = 200.0,5000.0,100.0 do begin

	create_test_particle_f, /weighted_maxwellian_xyz, rsfwc_1d='data/kj_single_1d.nc', $
			energy_keV = ee*1e-3

	spawn, '~/code/kineticj/bin/kineticj'

	kj_plot_current, sig33 = sig33

	sig33All1 = [sig33All1,sig33[0]]
	sig33All2 = [sig33All2,sig33[1]]

	energyAll = [energyAll,ee]

	print, 'Energy [eV]: ', EE, ' sig33: ', sig33[0]

endfor

stop

end

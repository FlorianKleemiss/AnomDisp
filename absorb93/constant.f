! CONSTANTS.FOR
	real degrees_to_radians, seconds_to_microrads, two_pi, pi
	1   , planck, fine_structure
	1	, h_bar_mev, hc, boltzmann, light_speed, electron_charge
	1	, electron_mass, electron_energy, electron_radius
	1	, avogadro, avogadro_cc_per_aa, gas_constant, f_to_mu
	1	, electron_ev_to_k, barns_to_electrons

	parameter (DEGREES_TO_RADIANS= 0.0174533D0)
	parameter (SECONDS_TO_MICRORADS= 4.8481368D0)
	parameter (TWO_PI= 6.283185D0)
	parameter (PI= 3.14159265D0)
	parameter (PLANCK= 6.6260755D-34)			! Joule-sec
	parameter (H_BAR_MEV= 6.5821220D-22)		! MeV-sec
	parameter (HC= 12398.4244D0)			! ev-A
	parameter (BOLTZMANN= 1.380658D-23)		! Joules/K
	parameter (LIGHT_SPEED= 2.99792458D8)		! meters/sec
	parameter (ELECTRON_CHARGE= 1.60217733D-19)	! coulombs
	parameter (ELECTRON_MASS= 9.1093897D-31)	! kg
	parameter (ELECTRON_ENERGY= 0.51099906D0)		! MeV
	parameter (ELECTRON_RADIUS= 2.81794092D-5)
	parameter (FINE_STRUCTURE= 7.29735308D0)		! 4pi*e_0*e^2/h_bar*c
	parameter (AVOGADRO= 6.0221367D23)			! atoms/mole
	parameter (AVOGADRO_CC_PER_AA= 0.60221367D0) 	! atoms-cc/mole-A^3
	parameter (GAS_CONSTANT= 1.987216D0)		! calories/mole
	parameter (F_TO_MU= 4208.031548D0)			! Mu in microns
	parameter (ELECTRON_EV_TO_K= 0.512316590D0)		! k in A^-1
	parameter (BARNS_TO_ELECTRONS= 1.43110541D-8) 	! barns to f"

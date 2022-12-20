		subroutine abs$struct_fact( h, k, l, energy, ipol, t_cryst
	1			, cd_strf, thetab, gamma, d_dsp, rho_ave, vol)
		implicit none
		include 'constant.f'
		integer*4 ibasis, nsites, natoms, ized, nstiff, iread, j, nedge
		real*4 fract, atomx, atomy, atomz, a_zero, b_zero, c_zero, alp
	1			, bet, gam, t_debye, temp_coef, stiff, t_cryst, amu_tot
	1			, erg, sinthe, emin, thetab, amu_site, pol_fact
	1			, vol, fpp_ave, rho, amu, ray, compt, rho_ave
	1			, abs$sfcoef, debwal, debye_funct, gamma, denom
		common /atom_pos/ ibasis, nsites(20), natoms(20), ized(20,5)
	1			, fract(20,5), atomx(20,100), atomy(20,100), atomz(20,100)
		common /crystal/ a_zero, b_zero, c_zero, alp, bet, gam
	1			, t_debye, temp_coef 		
		common /stiffness/ nstiff, stiff(10)
		common/ cromer_data/ erg, amu(20,5), rho(20,5), fp(20,5)
	1			, fpp(20,5), comp(20,5)
		integer*4 ipol, i
        real*4 form_fact, h, k, l, ksp, energy, fp, fpp, comp
		real*4 abs_edge(24)
		real*8 d_dsp, d_space
        complex*16 cd_strf, cd_form, cd_cryst
! -- executable code

		cd_strf= (0.d0,0.d0)
		d_dsp = d_space(h, k, l, t_cryst, vol)
		if( vol .eq. 0.) stop 'Crystal volume undefined'
		iread= 1
		if( energy .eq. erg) iread=0		! don't re-read if energy the same
		pol_fact = 1.						! Sigma default polarization
		if( d_dsp .gt. 0.) then
			ksp = TWO_PI/d_dsp
			sinthe = HC/(2.0* energy* d_dsp)
			if( sinthe .lt. 1.0) then
				thetab = asind(sngl( sinthe))
			else
				thetab= 0.
			end if
			if( ipol .eq. 2) pol_fact = abs(cosd(2.*thetab)) ! Pi polarization
		else
			ksp = 0.
		end if

		cd_strf= (0.,0.)
		amu_tot= 0.
! Calculate complex F for each reflection
		do i=1,ibasis				! loop over atomic sites
			amu_site= 0.
			form_fact = 0.
			fpp_ave = 0.
! now loop over different atomic species on each atomic site
			do j=1, natoms(i)
				if( iread) then				! have the data been read in ?
					call abs$atomic_data( ized(i,j), amu(i,j)
	1					, rho(i,j), nedge, abs_edge)
					call abs$cromer( ized(i,j), energy, fp(i,j), fpp(i,j))
					call abs$raycomp( ized(i,j), energy, ray, comp(i,j))
					erg = energy		! store current energy for next time
				end if
			    amu_site = amu_site + fract(i,j)*amu(i,j)
				form_fact = form_fact+ 
	1					fract(i,j)*(abs$sfcoef(ized(i,j), ksp, compt)+ 
	1			fp(i,j))
				fpp_ave = fpp_ave+ fract(i,j)*( fpp(i,j)+ comp(i,j))
			end do
			cd_form = dcmplx( form_fact, fpp_ave*pol_fact)
			denom = 4.*amu_site*(t_debye*d_dsp)**2
			if (t_cryst .le. 0.) t_cryst = .000001		!should be cold 'nuf
			if( denom .gt. 0.) then
				debwal = 1.14904e4*t_cryst*
	1					(debye_funct(t_debye/t_cryst) + (t_debye/t_cryst)/4.)
	1							/denom
			else
				debwal = 0.
			end if

			cd_cryst = (0.,0.)
!  Calculate crystal factor for each reflection
			do j=1, nsites(i)
				cd_cryst = cd_cryst+ cdexp((0.d0,1.d0)*TWO_PI
	1			*( atomx(i,j)*h + atomy(i,j)*k + atomz(i,j)*l))
			end do
			cd_strf = cd_strf + cd_form* cd_cryst* exp( -debwal)
			amu_tot= amu_tot+ amu_site* nsites(i)
		end do
		rho_ave= amu_tot/(AVOGADRO_CC_PER_AA* vol)
		gamma = (ELECTRON_RADIUS/(PI*vol))* (HC/energy)**2

		return
        end
!******************************************************************************
        real*8 function d_space( h, k, l, t_cryst, vol)
		implicit none
        real*4 h, k, l, vol, a_zero, b_zero, c_zero, alp, bet, gam
	1			, t_debye, temp_coef, t_cryst
		real*4 trig_vol, dsq_inv, calp, cbet, cgam, sbet, salp, sgam
	1			, a0_temp, b0_temp, c0_temp
		common /crystal/ a_zero, b_zero, c_zero, alp, bet, gam
	1			, t_debye, temp_coef 		
! -- executable code
		
		a0_temp= a_zero*(1.+temp_coef*1.e-6*(t_cryst-293.15))
		b0_temp= b_zero*(1.+temp_coef*1.e-6*(t_cryst-293.15))
		c0_temp= c_zero*(1.+temp_coef*1.e-6*(t_cryst-293.15))
        calp = cosd( alp)
        salp = sind( alp)
        cbet = cosd( bet)
        sbet = sind( bet)
        cgam = cosd( gam)
        sgam = sind( gam)
     
		trig_vol = 1.+ 2.* calp* cbet* cgam- calp**2- cbet**2- cgam**2
		vol = a0_temp*b0_temp*c0_temp* sqrt(trig_vol)
		dsq_inv=( h**2* salp**2/ a0_temp**2
	1			+ k**2* sbet**2/ b0_temp**2
	1			+ l**2* sgam**2/ c0_temp**2
	1			+2.* h* k*( calp* cbet- cgam)/( a0_temp* b0_temp)
	1			+2.* k* l*( cbet* cgam- calp)/( b0_temp* c0_temp)
	1			+2.* l* h*( cgam* calp- cbet)/( c0_temp* a0_temp))/trig_vol

		if( abs( dsq_inv) .gt. 1.e-12) then
			d_space = sqrt(1./dsq_inv)
		else
			d_space = 0.
		end if
        return
        end
!***********************************************************************
		real*4 function debye_funct (x)
		implicit none
		real*4 x, sum, xi_inc, xi, f2
!  This calculates the Debye function, PHI(x)
!  Written by: P.L.Cowan	1	8jan90		
!  Reference: International Tables for Crystallography (1985), 
!  vol 2. pp.241-265

		if (x.le.0.0001) then
			debye_funct = 1.000
		else if (x.ge.7.) then
			debye_funct = 1.642 / x
		else
			xi_inc = x / 100.001		!make xi_inc slightly too small
			sum = 0.						! so the last point will be less than x
			do xi=xi_inc,x,xi_inc
				f2 = xi / (exp(xi) - 1.)
				sum = sum + f2 * xi_inc						!box integration
			end do						!box integ. works somewhat better
			debye_funct = sum / x
		endif
		return
		end

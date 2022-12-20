		subroutine abs$absorb( name, nam_len, energy, mu, rho, ifail)
		implicit none
		include 'constant.f'
		real*4 energy, mu, rho, gamma, vol, h, k, l, temp, thetab
		real*4 fp, fpp, ray, comp, amu, abs_edge(24)
		real*8 d_space
		character*(*) name
		character*2 element_symbol
		integer*4 ipol, ifail, nam_len, n_edge, zed
		complex*16 cd_strf
		data h/0./, k/0./, l/0./, ipol/1/, temp/293./

! executable code
		
		if( nam_len .le. 2) then	! try an element first
		    element_symbol= name(1:2)
		    zed= -1
		    call abs$element( zed, element_symbol)
		    if( zed .gt. 0) then
				call abs$cromer( zed, energy, fp, fpp)
				call abs$raycomp( zed, energy, ray, comp)
				call abs$atomic_data( zed, amu, rho, n_edge, abs_edge)
				mu= (F_TO_MU* rho*( fpp+ ray+ comp))/( amu* energy)
			else
				ifail= 1
			end if
		else
			call abs$rdlattice( name, nam_len, ifail)
			if(.not. ifail) then
			    call abs$struct_fact( h, k, l, energy, ipol, temp
	1			, cd_strf, thetab, gamma, d_space, rho, vol)
			    if( energy .eq. 0.) stop 'Energy undefined'
			    if( vol .eq. 0.) stop 'Crystal volume undefined'
			    mu= ( 1.e4* TWO_PI* energy* dimag( cd_strf)* gamma)/HC
			else
				ifail= 0
				call abs$matl_abs( name, nam_len, energy, mu, rho, ifail)
			end if
		end if	
		return
		end

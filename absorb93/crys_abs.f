		subroutine abs$crys_abs( crys_name, nchrs, energy, mu, rho, noread)
		implicit none
		include 'constant.f'
		real*4 energy, mu, rho, gamma, vol, h, k, l, temp, thetab
		real*8 d_space
		character*(*) crys_name
		integer*4 ipol, noread, nchrs
		complex*16 cd_strf
		data h/0./, k/0./, l/0./, ipol/1/, temp/293./

! executable code
		
		call abs$rdlattice( nchrs, crys_name, noread)
		if( noread) then
			return
		else
		    call abs$struct_fact( h, k, l, energy, ipol, temp
	1			, cd_strf, thetab, gamma, d_space, rho, vol)
		    if( energy .eq. 0.) stop 'Energy undefined'
		    if( vol .eq. 0.) stop 'Crystal volume undefined'
		    mu= ( 1.e4* TWO_PI* energy* dimag( cd_strf)* gamma)/HC
		end if	
		return
		end

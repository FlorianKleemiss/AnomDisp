! ABS$MASS_ABS
! convert cross-section in electrons/atom to 1/cm
! SMB 22-APR-1991

		subroutine abs$mass_abs( zed, energy, mu, rho)
		implicit none
		include 'constant.f'
		integer zed, n_edge
		real*4 mu, energy, rho, amu, fp, fpp, ray, comp
		real*4 abs_edge(24)

! executable code

		call abs$atomic_data( zed, amu, rho, n_edge, abs_edge)
		if( zed .eq. 0) stop 'Element out of range'
		call abs$cromer( zed, energy, fp, fpp)
		call abs$raycomp( zed, energy, ray, comp)

		if( amu .le. 0.) then
				write(*,'('' Atomic weight undefined'')')
				return
		else if( energy .le. 0.) then
				write(*,'('' Energy undefined'')')
				return
		else if( fpp .le. 0.) then
				write(*,'('' Energy too small, make sure its eV'')')
				return
		end if

		mu= (F_TO_MU* rho*( fpp+ ray+ comp))/( amu* energy)
		return
		end

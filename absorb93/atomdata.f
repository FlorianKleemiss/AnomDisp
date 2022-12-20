C ABS$ATOMIC_DATA.FOR
C NEW 3-MAY-1993 SMB
		subroutine abs$atomic_data( zed, r_amu, r_rho, n_edge, abs_edge)
		implicit none
		include 'amrhoin.f'
		include 'bdnrgin.f'
		integer*4 i, zed, n_edge
		real*4 abs_edge(24), r_amu, r_rho

C executable code

		if(zed .lt.1 .or. zed .gt. 92) then
			write(*,'('' Atomic number out of range'')')
			zed= 0
		else
			r_amu = amu( zed)
			r_rho = rho( zed)

C bind_nrg(1) = K in keV
C bind_nrg(2) = LI "  "
C bind_nrg(5) = M1 "  "

			n_edge = n_orb(zed)	
			do i=1,n_edge
				abs_edge(i) = bind_nrg(zed,i)
			end do

		end if
		return
		end

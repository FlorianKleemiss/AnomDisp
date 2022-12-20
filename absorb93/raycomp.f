! ABS$RAYCOMP.FOR
! 20-APR-91 SMB
! Calculate elastic and compton cross-sections (in electrons/atom)
! 4-MAY-1993 SMB Rewritten for internal data

		subroutine abs$raycomp( zed, energy, ray, comp)
		implicit none
		include 'constant.f'
		include 'raycomin.f'

		real*4 energy, p1, p2, p3
		real*4 ray, comp
		integer zed
! BARNS_TO_ELECTRONS is a constant that converts the cross section in 
! barns/atom into a cross-section in electrons/atom which is the same as f''.
! it is calculated as
! BARNS_TO_ELECTRONS= 1/( 10^8 [b/A^2]*2* rsube=2.8179^-5* hc=12398.52)
! cf V, p79,80 for details
! executable code

		if( energy .le. 0.) then

			write(*,'('' Energy undefined'')')
			zed= 0
			return
		else if(zed .lt.1 .or. zed .gt. 92) then

			write(*,'('' Atomic number out of range'')')
			zed= 0
			return
		end if

		p1 = alog(energy/1000.)
		p2 = (p1)**2
		p3 = (p1)**3

		ray= BARNS_TO_ELECTRONS* energy* exp( dray(zed,1)+ dray(zed,2)* p1+ 
	1			dray(zed,3)* p2+ dray(zed,4)* p3)

		comp= BARNS_TO_ELECTRONS* energy* exp( dcmp(zed,1)+ dcmp(zed,2)* p1+
	1			dcmp(zed,3)* p2+ dcmp(zed,4)* p3)

		return
		end

		subroutine abs$matl_abs( name, nam_len, energy, mu, rho_ave, ifail)
		implicit none
		include 'constant.f'
		real*4 energy, mu, fp, abs_edge(24), fpp(20), ray(20), comp(20)
	1		, fract(20), rho(20), amu(20), density, amu_ave, rho_ave
		integer*4 zed(20), i, ncomp, nedge, ifail
		character*(*) name
		character*80 string
		integer*4 nam_len
! executable code
		
		ifail= 0
		mu= 0.
		rho_ave= 0.
		amu_ave= 0.

		open(unit=18, file=name(1:nam_len)//'.mtl', status='old',readonly
	1	, err=999)

		read(18,'(a)') string
		read(18,'(a)') string
		read(18,*,end=998,err=998)
	1		density, ncomp,(zed(i),fract(i),i=1,ncomp)

		close(unit=18)

		do i=1, ncomp
			call abs$cromer( zed(i), energy, fp, fpp(i))
			call abs$raycomp( zed(i), energy, ray(i), comp(i))
			call abs$atomic_data( zed(i), amu(i), rho(i), nedge, abs_edge)
			mu= mu+ fract(i)* ( fpp(i)+ ray(i)+ comp(i))
			amu_ave= amu_ave+ fract(i)* amu(i)
			rho_ave= rho_ave+ fract(i)* rho(i)
		end do

		if( density .gt. 0.) rho_ave= density
		if( energy .eq. 0) stop 'Energy undefined'
		if( amu_ave .eq. 0) stop 'Atomic Weight undefined'

		mu= (F_TO_MU* rho_ave* mu)/( energy* amu_ave)
		return
999		ifail= 1
		return
998		stop 'Error reading Material data file'
		end

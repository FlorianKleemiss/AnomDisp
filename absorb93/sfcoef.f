!ABS$SFCOEF
!MODIFIED 20-APR-91 SMB
		real*4 function abs$sfcoef( zed, k, compton)
		implicit none
		include 'sf_inc.f'
!  this function calculates f as f(sin(th)/lamda) given zed.
!  it uses the parametrized fits of Cromer and Mann, Acta Cryst. A24
!  321 (1968).  The function also calculates the compton component 
!  using the fits of Balyuzi, Acta. Cryst. A31 600 (1975).
!    NOTE: K IS ASSUMED TO BE 4PISIN(THETA)/LAMBDA !!!!!!!!!!!!!!
!    IT IS CONVERTED TO SIN(TH)/LAMBDA WITHIN THE ROUTINE !!!!!!!
!  compton = incoherent portion of scattering factor
		real*4 k, compton, ksqr, sfr
		integer i, zed

! executable code
		if(( zed .lt. 1).or.( zed.gt. 92)) then
			write(*,'('' Atomic # out of range'')')
			zed= 0
			return
		end if

		if( k .gt. 18.5) then
			write(*,'('' K value too large'')')
			zed= 0
			return
		end if

! 0.0063326=(1/4pi)**2 to convert from k to sin(th)/lambda

		ksqr = -0.0063326* k**2
		sfr=find(zed,9)
		do i=1,8,2
			sfr= sfr+ find(zed,i)* exp( find(zed,i+1)* ksqr)
		end do
		abs$sfcoef = sfr

		compton = 0.		! now calculate compton component
		do i=1,10,2
		   compton= compton+ comp(zed,i)* exp( ksqr* comp(zed,i+1))
		end do
		compton= float( zed) - compton
		return
		end

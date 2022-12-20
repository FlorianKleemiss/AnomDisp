C ABS$CROMER.FOR
C INCLUDES:		ABS$CROMER
C				SIGMA0
C				SIGMA1
C				SIGMA2
C				SIGMA3
C				LGNDR
C				GAUSS
C				AKNINT
C				SORT
C				MCM
C This routine reads data for f' and f" according to an
C algorithm by cromer and lieberman, given to fuoss.
C converted to direct access file by brennan
C converted to internal data 3-may-1993 smb

       subroutine abs_cromer( zed, energy, fp, fpp)
		implicit none
		include 'bdnrgin.f'
		include 'relinc.f'
		include 'xnrg_com.f'
		include 'xsec_com.f'

		integer zed, irb, ipr, icount, iopen, jj, i_zero, i_nxsect

		real energy, fp, fpp

		real
	1			 d_fp_orb			,d_fpp_orb
	1			,d_nrg_s(11)		,d_xsect_s(11)
	1			,d_log_nrg(11)		,d_log_xsect(11)
	1			,d_nrg_int(5)		,d_xsect_int(5)
	1			,d_xsect_edge_au	,d_fp_corr
	1			,d_energy			,d_log_energy
	1			,d_energy_au		,d_bind_nrg_au
	1			,d_xsect_barns		,d_var
	1			,d_sum_fp			,d_sum_fpp
	1			,AU   ,KEV_PER_RYD  ,PI  ,INV_FINE_STRUCT
C                FINE_PI= 1/(alpha*2*pi**2)
	1			,FINE_PI
	1			,sigma0		,sigma1		,sigma2		,sigma3		,gauss
	1			,aknint

		DATA AU/2.80022d+7/		,KEV_PER_RYD/0.02721D0/
		DATA FINE_PI/6.942325/		,PI/3.1415927D0/
		DATA INV_FINE_STRUCT/1.37036D2/
		COMMON /GAUS/ d_xsect_barns
	1			,d_bind_nrg_au
	1			,d_xsect_int
	1			,d_energy_au
	1			,icount
		COMMON /EDGE/ d_xsect_edge_au
		EXTERNAL sigma0,sigma1,sigma2,sigma3

C relcor is the relativistic correction term
C executable code
		if(zed .lt.1 .or. zed .gt. 92) then
			write(*,'('' Atomic number out of range'')')
			zed= 0
			return
		else if( energy .le. 0.) then
			write(*,'('' Energy undefined'')')
			zed= 0
			return
		else if( zed .lt. 3) then

			call mcm( zed, energy, fp, fpp)
			iopen= zed
			return
		end if

C       Change to keV
		d_energy = energy/1000.
                write(*,*) 'energy in keV selected: ',d_energy
		d_log_energy = log(d_energy)
		d_energy_au = d_energy/KEV_PER_RYD
		write(*,*) 'energy in a.u.: ',d_energy_au

		d_sum_fp = 0.d0
		d_sum_fpp = 0.d0

C MAIN LOOP THROUGH THE ORBITALS
		do irb=1,n_orb(zed)
			icount= 6
			d_fp_orb = 0.d0
			d_fp_corr = 0.d0
			d_fpp_orb = 0.d0
			d_xsect_barns = 0.d0
			d_bind_nrg_au = bind_nrg(zed,irb)/KEV_PER_RYD

			if(nparms(zed,irb) .eq. 11) then
				d_xsect_edge_au = xsc(zed,irb,11)/AU
			end if

C also copy subset into second array
			do ipr=6,10
				jj = ipr-5
				d_xsect_int(jj) = xsc(zed,irb,ipr)/AU
				d_nrg_int(jj) = nrg(zed,irb,ipr)
			end do

C the sorting routine messes up subsequent calls with same energy
C so copy to second array before sorting.
			do ipr=1,nparms(zed,irb)
				d_xsect_s(ipr) = xsc(zed,irb,ipr)
				d_nrg_s(ipr) = nrg(zed,irb,ipr)
			end do

			call sort(nparms(zed,irb),d_nrg_s,d_xsect_s)
			call sort(5,d_nrg_int,d_xsect_int)

C           convert to log of energy,xsect
			do ipr=1,nparms(zed,irb)
				d_log_nrg(ipr) = log(d_nrg_s(ipr))
				if(d_xsect_s(ipr) .eq. 0.0) then
						d_log_xsect(ipr) = 0.0
				else
						d_log_xsect(ipr) = log(d_xsect_s(ipr))
				end if
			end do

			if(bind_nrg(zed,irb) .le. d_energy) then
				i_zero = 0
				ipr = 1
				do while(d_log_xsect(ipr) .eq. 0.d0)
						i_zero = i_zero + 1
						ipr = ipr + 1
				end do
				i_nxsect = nparms(zed,irb) - i_zero
				i_zero = i_zero + 1
				d_xsect_barns = aknint(d_log_energy,i_nxsect
	1					,d_log_nrg(i_zero),d_log_xsect(i_zero))
				d_xsect_barns = exp(d_xsect_barns)/AU

				d_fpp_orb = INV_FINE_STRUCT*d_xsect_barns
	1							*d_energy_au/(4.d0*PI)
				d_var = d_energy_au-d_bind_nrg_au
				if(d_var .eq. 0.) d_var = 1.
				d_fp_corr = -0.5*d_xsect_barns*d_energy_au
	1					* FINE_PI
	1					* log((d_energy_au+d_bind_nrg_au)
	1							/d_var)
		   end if

		   if(bind_nrg(zed,irb) .gt. d_energy .and. 
	1					funtype(zed,irb) .eq. 0) then
				d_fp_orb = gauss(sigma3) * FINE_PI

				d_fp_corr = 0.5*d_xsect_edge_au*d_bind_nrg_au**2
	1					* log((-d_bind_nrg_au+d_energy_au)
	1							/(-d_bind_nrg_au-d_energy_au))
	1					/ d_energy_au * FINE_PI
		   else
		   		if(funtype(zed,irb) .eq. 0)
	1				 d_fp_orb = gauss(sigma0) * FINE_PI
		   		if(funtype(zed,irb) .eq. 1)
	1				 d_fp_orb = gauss(sigma1) * FINE_PI
		   		if(funtype(zed,irb) .eq. 2)
	1				 d_fp_orb = gauss(sigma2) * FINE_PI
		   end if

		   d_sum_fp = d_sum_fp + d_fp_orb + d_fp_corr
		   d_sum_fpp = d_sum_fpp + d_fpp_orb

		end do
C this is the end of the loop over orbits

		fpp = d_sum_fpp

C
C Note: the Jensen correction to f' was subsequently shown to be incorrect
C (see L. Kissel and R.H. Pratt, Acta Cryst. A46, 170 (1990))
C and that the relativistic correction that Ludwig used is also
C wrong.  This section retained as comments for historical reasons.
C
C		jensen_cor = -0.5*float(zed)
C     1			*(d_energy_au/INV_FINE_STRUCT**2)**2
C
C       Subtract relcor ala ludwig and change back to real*4
C
C		fp = d_sum_fp+jensen_cor-relcor(zed)
C
C Kissel and Pratt give better corrections.  The relativistic correction
C that Ludwig used is (5/3)(E_tot/mc^2).  Kissel and Pratt say that this
C should be simply (E_tot/mc^2), but their correction (KPCOR) apparently
C takes this into account.  So we can use the old RELCOR and simply add
C the (energy independent) KPCOR term:
C
		write(*,*) 'wihtout correction:',d_sum_fp
		fp = d_sum_fp-relcor(zed)+kpcor(zed)

		return
		end
C***********************************************************************
		real function sigma0( x)
		implicit none
		real x, d_xsect_barns, d_bind_nrg_au, d_xsect_int(5)
	1			, d_energy_au, d_prod
		COMMON /GAUS/ d_xsect_barns, d_bind_nrg_au, d_xsect_int
	1			, d_energy_au, icount
		integer icount
C executable code

		icount = icount-1
		d_prod = d_energy_au**2* x**2- d_bind_nrg_au**2

		if( abs( d_prod) .lt. 1.d-30) then
			sigma0 = d_xsect_int(icount)* d_bind_nrg_au/x**2
		else
			sigma0 =( d_xsect_int(icount)* d_bind_nrg_au**3/ x**2
	1			- d_bind_nrg_au* d_xsect_barns* d_energy_au**2)/ d_prod
		end if
		return
		end
C***********************************************************************
		real function sigma1( x)
		implicit none
		real x, d_xsect_barns, d_bind_nrg_au, d_xsect_int(5)
	1			, d_energy_au
		integer icount
		COMMON /GAUS/ d_xsect_barns, d_bind_nrg_au, d_xsect_int
	1			, d_energy_au, icount
C executable code

		icount = icount-1
		sigma1 = 0.5* d_bind_nrg_au**3* d_xsect_int( icount)
	1			/( sqrt( x)*( d_energy_au**2* x**2- d_bind_nrg_au**2* x))

		return
		end
C***********************************************************************
		real function sigma2(x)
		implicit none
		real x
		real d_xsect_barns, d_bind_nrg_au, d_xsect_int(5)
	1			, d_energy_au
		INTEGER icount
		COMMON /GAUS/ d_xsect_barns, d_bind_nrg_au, d_xsect_int
	1			, d_energy_au, icount
		real denom
		icount=icount-1
C Code modes by Chris T. Chantler, May 12-1992
		IF (abs(x).LT.1.0D-31) THEN
C		 WRITE(5,*) 'Denom. overflow'
		 sigma2=0.0
		ELSEIF (d_energy_au.LT.1.0D-31) THEN
		 sigma2=0.0
C		 WRITE (5,*) 'e simeq 0'
		ELSEIF (abs(d_xsect_int(icount)-d_xsect_barns).LT.1.0D-30)
	1	THEN
C		 WRITE(5,*) 'Factor cancels'
		 sigma2=-2.0d00*d_xsect_int(icount)*d_bind_nrg_au/x**3
        ELSE
		 denom= x**3*d_energy_au**2-d_bind_nrg_au**2/ x
		 IF (abs(denom).LT.1.0D-31) THEN
C		  WRITE(5,*) 'Denom. zero'
		  sigma2=-2.0d00*d_xsect_int(icount)*d_bind_nrg_au/x**3
		 ELSE
		  sigma2=2.0d00*(d_xsect_int(icount)*(d_bind_nrg_au/x)**3/x-
	1	d_bind_nrg_au* d_xsect_barns* d_energy_au**2)/denom
		 ENDIF
		ENDIF
		RETURN
		END
C***********************************************************************
		real function sigma3( x)
		implicit none
		real x, d_xsect_barns, d_bind_nrg_au, d_xsect_int(5)
	1			, d_energy_au, d_xsect_edge_au
		integer icount
		COMMON /EDGE/ d_xsect_edge_au
		COMMON /GAUS/ d_xsect_barns,d_bind_nrg_au,d_xsect_int
	1			,d_energy_au,icount
C executable code
C		write(*,*) 'icount: ',icount
		icount = icount-1
		sigma3 = d_bind_nrg_au**3*( d_xsect_int( icount)
	1			- d_xsect_edge_au* x**2)
	1			/( x**2*( x**2* d_energy_au**2- d_bind_nrg_au**2))

		return
		end
C***********************************************************************
		subroutine lgndr (index,d_bb,d_cc)
		implicit none
		integer index, ip
		real d_bb, d_cc, d_const, d_x(2), d_a(3)
		data d_x(1) /.04691007703067D0/
		data d_x(2) /.23076534494716D0/
		data d_a(1) /.11846344252810D0/
		data d_a(2) /.23931433524968D0/
		data d_a(3) /.28444444444444D0/
C executable code

C  WARNING! THIS ROUTINE HAS BEEN STRIPPED SO IT IS ONLY USEFUL
C  WITH ABS$CROMER IN THIS SET OF ROUTINES.

		d_cc = 0.5d0
		d_const=0.d0
		ip= index
C       ip limited to 1,2,3

		if ( ip .gt. 3) then
			 ip= 6- ip
			 d_const= -1.d0
		end if
		d_bb = d_a( ip)
		if( ip .eq. 3) return
		d_cc= -d_const+ sign( d_x( ip), d_const)
		return
		end
C***********************************************************************

		real function gauss (sigma)
		implicit none
		integer i, index
		real aa, bb, cc, sigma, dd
		external sigma
C executable code

		aa = 0.d0
		dd = 0.d0
		do i=1,5
			index = i
			call lgndr( index, bb, cc)
C			write(*,*) aa, dd
			dd = sigma(cc)
			aa = aa + bb * dd
		end do
		gauss = aa
		return
		end
C*************************************************************
C aitken repeated interpolation
C   d_log_energy = abscissa at which interpolation is desired
C   d_log_nrg    = vector of n values of abscissa
C   d_log_xsect    = vector of n values of ordinate
C   t    = temporary storage vector of 4*(m+1) locations)
		real function aknint 
	1			( d_log_energy, n, d_log_nrg, d_log_xsect)
		implicit none
		integer n, i, index, mend, kk, jj
		real t(20), d_log_energy, d_log_nrg( n)
	1			, d_log_xsect( n), del
C executable code
		index= 0
		if(n .le. 2) then
C			 write(*,'('' Too few points, funct=y(1)'')')
			 aknint = d_log_xsect(1)
			 return
		end if
		del = d_log_nrg(2) - d_log_nrg(1)
C		write(*,*) 'energy target: ',d_log_energy
		if(del .gt. 0.) then
			i=1
			do while(d_log_nrg(i) .lt. d_log_energy .and. i .le. n)
				index = i
				i=i+1
			end do
		else
			i=1
			do while(d_log_nrg(i) .gt. d_log_energy .and. i .le. n)
				index = i
				i=i+1
			end do
		end if
C		write(*,*) 'd: ',del,' index before: ',index
		index = index - 1
		index = max0( index,1)
		index = min0( index, n-2)
C		write(*,*) 'index: ',index
		mend = index+2
		do i= index, mend
			kk= i- index+1
			t( kk)= d_log_xsect( i)
			t( kk+3)= d_log_nrg( i)- d_log_energy
		end do
C        write(*,*) t
C		write(*,*) d_log_nrg
		do i=1,2
C			write(*,*) 'i=',i
			kk=i+1
			do jj=kk,3
C				write(*,*) 'j=',jj
C				write(*,*) d_log_nrg( jj+ index-1)
C				write(*,*) d_log_nrg( i+ index-1)
					t( jj) = ( t( i)*t( jj+3)-t( jj)*t( i+3))
	1			/( d_log_nrg( jj+ index-1)- d_log_nrg( i+ index-1))
			end do
		end do
C		write(*,*) t

		aknint= t(3)
		return
		end
C***********************************************
C bubble sort.  largest becomes last
		subroutine sort (n,a,b)
		implicit none
		integer i, n, j
		real a(1), b(1), x, y
C executable code

		do i=1,n-1
			do j=i+1,n
				if(a(j).lt.a(i)) then
					 x=a(j)
					 y=a(i)
					 a(i)=x
					 a(j)=y
					 x=b(j)
					 y=b(i)
					 b(i)=x
					 b(j)=y
				end if
			end do
		end do
		return
		end
C**************************************************************
		subroutine mcm( zed, energy, fp, fpp)
		implicit none
		include 'constant.f'
		real energy, fp, fpp, p1, p2, p3
		integer zed
C barns_to_electrons is a constant that converts the cross section in 
C barns/atom C into a cross-section in electrons/atom which is the same 
C as f''.  C it is calculated as
C  1/( 10^8 [b/A^2]*2* rsube=2.8179^-5* hc=12398.52)
C executable code

		fp = 0.
		P1 = ALOG( energy/1000.)
		P2 = (P1)**2
		P3 = (P1)**3


		if( zed .eq. 1) then
C       Hydrogen
			if( energy .ge. 14.e-3) then
				fpp= BARNS_TO_ELECTRONS* energy* 
	1			exp( 2.44964- 3.34953* p1 - 0.047137*p2 + 0.0070996*p3)
			else
				fpp= 0.
			end if

		else
C       Helium
			if( energy .ge. 25.e-3) then
				fpp= BARNS_TO_ELECTRONS* energy* 
	1			exp( 6.06488- 3.2905* p1 - 0.107256* p2 + 0.0144465* p3)
			else
				fpp= 0.
			end if
		end if
		return
		end

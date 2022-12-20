C ABSORPTION
C INCLUDES
C	TYPE_DATA
C	OUTPUT_DEFAULT
C	TYPE_MU
C	TYPE_EDGES
C   CONVERT
C   WIDTH (FUNCTION)
C 5-MAY-1993 SMB

		program absorption
		implicit none
		integer*4 oncet, once, elem_len, erg_len
	1		, th_len, first_time, ifail, width, len, places
		real*4 energy, thick, mu
		character*80 element_string
		character*80 energy_string, thick_string
        external XNRG_BLOCK,XSECT ! needed for DEC linker to resolve block data     
C executable code

		write
	1(*,'('' Calculations based on Cromer and Liberman and McMasters'')')
		
		do while (.true.)

		ifail= 0
		first_time= 1
		energy= 8000.
		thick= 100.
		
		write(*,'('' Element(e.g. Si): '',$)')
		call abs$get_string( element_string, elem_len)
		if( elem_len .le. 0) stop 'Glad to be of service!'
			
		once= 1
		do while( energy .gt. 0.)
		    write(*,'('' Enter energy('',f<width(energy)>.0
	1			,'')eV: '',$)') energy
			call abs$get_string( energy_string, erg_len)
			if( erg_len .lt. 0) goto 999
			if( erg_len .eq. 0) then
				if( once) then
					once= 0
				else
					goto 999
				end if
			else
				once= 0
				call default_energy( energy_string, erg_len, energy)
			end if

			call parse_name( element_string, elem_len, energy, mu
	1			, first_time, ifail)
			if( ifail) goto 999

			oncet= 2
			
			do while(( thick.gt. 0.).and.( oncet.gt. 0))
				len= width(thick)
				places= 0
				if( len .lt. 3) then
					len= len+ 3
					places= 3
				end if
			    write(*,'('' Enter thickness(''
	1			,f<len>.<places>,'')microns: '',$)') thick
				call abs$get_string( thick_string, th_len)
				call default_thick( thick_string, th_len, thick)
				if( th_len.lt. 0) then
					goto 999
				else if( th_len .eq. 0) then
					if( oncet.gt. 1) then
						oncet= 1
						call type_mu( mu, thick)
					else
						oncet= 0
					end if
				else
					oncet= 1
					call type_mu( mu, thick)
				end if
			end do
		end do

999		continue
		end do
		end
C*********************************************************************
		subroutine type_energy_data( energy, fp, fpp, ray, comp)
		implicit none
		real*4 energy, fp, fpp, ray, comp
		integer*4 width
C executable code

		write(*,'('' fp:'',1pg10.3
	1		,'' fpp:'',g10.3
	1		,'' Compton:'',g10.3
	1		,'' Rayleigh:'',g10.3
	1		,'' At '',0pf<width(energy)>.0,''eV'')') 
	1		fp, fpp, comp, ray, energy

		return
		end

C***********************************************************************
		subroutine parse_name( name, nam_len, energy, mu
	1		, first_time, ifail)
		implicit none
		include 'constant.f'

		character*80 name
		character*2 element_symbol
		integer*4 nam_len, ifail, first_time
		integer*4 zed, n_edge
		real*4 energy, mu, fp, fpp, ray, comp, rho
	1	, amu, abs_edge(24)
C executable code
		
		if( nam_len .le. 2) then
		    element_symbol= name(1:2)
		    zed= -1
		    call abs$element( zed, element_symbol)
		    if( zed .gt. 0) then
				call abs$atomic_data( zed, amu, rho, n_edge, abs_edge)
				if( first_time) then
					first_time= 0
					call type_element_data( element_symbol, zed, amu
	1					, rho, n_edge, abs_edge)
				end if
				call abs$cromer( zed, energy, fp, fpp)
				call abs$raycomp( zed, energy, ray, comp)
				call type_energy_data( energy, fp, fpp, ray, comp)

				if( amu .le. 0.) stop 'Atomic weight undefined'
				if( energy .le. 0.) stop 'Energy undefined'
				mu= (F_TO_MU* rho*( fpp+ ray+ comp))/( amu* energy)
		    end if
		else
		    call abs$absorb( name, nam_len, energy, mu, rho, ifail)
			if(.not. ifail) then
				call type_crysmtl( name, nam_len, energy, rho)
			end if
		end if
		return
		end
C***********************************************************************
		subroutine type_crysmtl( string, len, energy, rho)
		implicit none
		character*(*) string
		integer*4 len
		real*4 rho, energy
C executable code


		write(*,'(x,a,'' has a density of'',1pg10.3,'' gms/cc'')')
	1		string(1:len), rho
		write(*,'('' At '',f8.0,'' eV'')') energy
		return
		end

C***********************************************************************
		subroutine type_mu( mu, thick)
		implicit none
		real*4 mu, mu_t, tran, thick
C executable code

		mu_t = mu* thick
		tran = exp(-mu_t)
		if( mu .eq. 0.) mu = 1.e-10
		write(*,'('' Trans. Int.= '',1pg12.5
	1		,'' At'',g10.3,''microns; Absorp. len.='',g10.3
	1		,'' microns'')') tran, thick, 1./mu
		return
		end
C************************************************************************
		subroutine type_element_data( string, zed, amu, rho
	1		, n_edge, abs_edge)
		implicit none
		character*2 string
		integer*4 n_edge, i, width, zed
		real*4 abs_edge(24), ka1, ka2, kb1, kb3, lb1, lb2, lb3, lb4
	1		, la1, la2, lg1, lg3, li, le, amu, rho
C executable code

		write(*,'('' Element: '',a,'' Atomic #'',i3,'' Atomic Weight= ''
	1		,f6.2,'' Density= '',1pg10.3,'' grams/cc'')') string, zed
	1		, amu, rho

		do i=1,n_edge
			abs_edge(i)= abs_edge(i)*1000.
		end do

		if( n_edge .eq. 1) then
			write(*,'('' K-edge: '',f<width(abs_edge(1))>.0)') abs_edge(1)
		else if( n_edge .eq. 2) then
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0)') abs_edge(1), abs_edge(2)

		else if( n_edge .eq. 3) then
			ka2= abs_edge(1)- abs_edge(3)
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3)

			write(*,'('' Ka2: '',f<width(ka2)>.0)') ka2

		else if( n_edge .eq. 4) then
			ka1= abs_edge(1)- abs_edge(4)
			ka2= abs_edge(1)- abs_edge(3)
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0
	1			,'' L3: '',f<width(abs_edge(4))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3), abs_edge(4)

			write(*,'('' Ka1: '',f<width(ka1)>.0
	1			,'' Ka2: '',f<width(ka2)>.0)') ka1, ka2

		else if( n_edge .eq. 5) then
			ka1= abs_edge(1)- abs_edge(4)
			ka2= abs_edge(1)- abs_edge(3)
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0
	1			,'' L3: '',f<width(abs_edge(4))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3), abs_edge(4)

			write(*,'('' Ka1: '',f<width(ka1)>.0
	1			,'' Ka2: '',f<width(ka2)>.0)') ka1, ka2

			write(*,'('' M1: '',f<width(abs_edge(5))>.0)') abs_edge(5)

		else if( n_edge .eq. 6) then
			kb3= abs_edge(1)- abs_edge(6)
			ka1= abs_edge(1)- abs_edge(4)
			ka2= abs_edge(1)- abs_edge(3)
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0
	1			,'' L3: '',f<width(abs_edge(4))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3), abs_edge(4)

			write(*,'('' Ka1: '',f<width(ka1)>.0
	1			,'' Ka2: '',f<width(ka2)>.0
	1			,'' Kb3: '',f<width(kb3)>.0)') ka1, ka2, kb3

			write(*,'('' M1: '',f<width(abs_edge(5))>.0
	1			,'' M2: '',f<width(abs_edge(6))>.0)')
	1			abs_edge(5), abs_edge(6)

		else if( n_edge .eq. 7) then
			kb1= abs_edge(1)- abs_edge(7)
			kb3= abs_edge(1)- abs_edge(6)
			ka1= abs_edge(1)- abs_edge(4)
			ka2= abs_edge(1)- abs_edge(3)
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0
	1			,'' L3: '',f<width(abs_edge(4))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3), abs_edge(4)

			write(*,'('' Ka1: '',f<width(ka1)>.0
	1			,'' Ka2: '',f<width(ka2)>.0
	1			,'' Kb1: '',f<width(kb1)>.0
	1			,'' Kb3: '',f<width(kb3)>.0)') ka1, ka2, kb1, kb3

			write(*,'('' M1: '',f<width(abs_edge(5))>.0
	1			,'' M2: '',f<width(abs_edge(6))>.0
	1			,'' M3: '',f<width(abs_edge(7))>.0)')
	1			abs_edge(5), abs_edge(6), abs_edge(7)

		else if( n_edge .eq. 8) then
			lb1= abs_edge(2)- abs_edge(8)
			kb1= abs_edge(1)- abs_edge(7)
			kb3= abs_edge(1)- abs_edge(6)
			ka1= abs_edge(1)- abs_edge(4)
			ka2= abs_edge(1)- abs_edge(3)
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0
	1			,'' L3: '',f<width(abs_edge(4))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3), abs_edge(4)

			write(*,'('' Ka1: '',f<width(ka1)>.0
	1			,'' Ka2: '',f<width(ka2)>.0
	1			,'' Kb1: '',f<width(kb1)>.0
	1			,'' Kb3: '',f<width(kb3)>.0)')
	1			, ka1, ka2, kb1, kb3

			write(*,'('' M1: '',f<width(abs_edge(5))>.0
	1			,'' M2: '',f<width(abs_edge(6))>.0
	1			,'' M3: '',f<width(abs_edge(7))>.0
	1			,'' M4: '',f<width(abs_edge(8))>.0
	1			,'' Lb1: '',f<width(lb1)>.0)')
	1			abs_edge(5), abs_edge(6), abs_edge(7), abs_edge(8)
	1			, lb1

		else if(( n_edge .ge. 9).and.(n_edge .lt. 14)) then
			lb1= abs_edge(2)- abs_edge(8)
			la1= abs_edge(3)- abs_edge(9)
			kb1= abs_edge(1)- abs_edge(7)
			kb3= abs_edge(1)- abs_edge(6)
			ka1= abs_edge(1)- abs_edge(4)
			ka2= abs_edge(1)- abs_edge(3)
			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0
	1			,'' L3: '',f<width(abs_edge(4))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3), abs_edge(4)

			write(*,'('' Ka1: '',f<width(ka1)>.0
	1			,'' Ka2: '',f<width(ka2)>.0
	1			,'' Kb1: '',f<width(kb1)>.0
	1			,'' Kb3: '',f<width(kb3)>.0)')
	1			, ka1, ka2, kb1, kb3

			write(*,'('' M1: '',f<width(abs_edge(5))>.0
	1			,'' M2: '',f<width(abs_edge(6))>.0
	1			,'' M3: '',f<width(abs_edge(7))>.0
	1			,'' M4: '',f<width(abs_edge(8))>.0
	1			,'' M5: '',f<width(abs_edge(9))>.0)')
	1		abs_edge(5), abs_edge(6), abs_edge(7), abs_edge(8), abs_edge(9)

			write(*,'('' Lb1: '',f<width(lb1)>.0,'' La1: '',f<width(la1)>.0)')
	1			, lb1, la1

C		 n_edge >13
		else
			ka1= abs_edge(1)- abs_edge(4)
			ka2= abs_edge(1)- abs_edge(3)
			kb1= abs_edge(1)- abs_edge(7)
			kb3= abs_edge(1)- abs_edge(6)
			lb1= abs_edge(3)- abs_edge(8)
			lb2= abs_edge(4)- abs_edge(14)
			lb3= abs_edge(2)- abs_edge(7)
			lb4= abs_edge(2)- abs_edge(6)
			la1= abs_edge(4)- abs_edge(9)
			la2= abs_edge(4)- abs_edge(8)
			lg1= abs_edge(3)- abs_edge(13)
			lg3= abs_edge(2)- abs_edge(12)
			li = abs_edge(4)- abs_edge(5)
			le = abs_edge(3)- abs_edge(5)

			write(*,'('' K: '',f<width(abs_edge(1))>.0
	1			,'' L1: '',f<width(abs_edge(2))>.0
	1			,'' L2: '',f<width(abs_edge(3))>.0
	1			,'' L3: '',f<width(abs_edge(4))>.0)')
	1			abs_edge(1), abs_edge(2), abs_edge(3), abs_edge(4)

			write(*,'('' Ka1: '',f<width(ka1)>.0
	1			,'' Ka2: '',f<width(ka2)>.0
	1			,'' Kb1: '',f<width(kb1)>.0
	1			,'' Kb3: '',f<width(kb3)>.0)')
	1			, ka1, ka2, kb1, kb3

			write(*,'('' M1: '',f<width(abs_edge(5))>.0
	1			,'' M2: '',f<width(abs_edge(6))>.0
	1			,'' M3: '',f<width(abs_edge(7))>.0
	1			,'' M4: '',f<width(abs_edge(8))>.0
	1			,'' M5: '',f<width(abs_edge(9))>.0)')
	1		abs_edge(5), abs_edge(6), abs_edge(7), abs_edge(8), abs_edge(9)

			write(*,'('' Lb1: '',f<width(lb1)>.0
	1		,'' Lb2: '',f<width(lb2)>.0
	1		,'' Lb3: '',f<width(lb3)>.0
	1		,'' Lb4: '',f<width(lb4)>.0
	1		,'' La1: '',f<width(la1)>.0
	1		,'' La2: '',f<width(la2)>.0)') lb1,lb2,lb3,lb4,la1,la2

			write(*,'('' N1: '',f<width(abs_edge(5))>.0
	1			,'' N2: '',f<width(abs_edge(6))>.0
	1			,'' N3: '',f<width(abs_edge(7))>.0
	1			,'' N4: '',f<width(abs_edge(8))>.0
	1			,'' N5: '',f<width(abs_edge(9))>.0)')
	1		abs_edge(10), abs_edge(11), abs_edge(12), abs_edge(13)
	1		, abs_edge(14)

		write(*,'('' Lg1: '',f<width(lg1)>.0
	1		,'' Lg3: '',f<width(lg3)>.0
	1		,'' Li: '',f<width(li)>.0
	1		,'' Le: '',f<width(le)>.0)') lg1,lg3,li,le
		end if

		return
		end
C***********************************************************************
		subroutine default_energy( energy_string, erg_len, energy)
		implicit none
		real*4 rval, energy
		integer*4 erg_len
		character*(*) energy_string
		if( erg_len .gt. 0) then
			call convert( energy_string, erg_len, rval)
			if( erg_len .lt. 0 .or. rval .lt. 1.) then
			    energy= 8000.
			else
			    energy= rval
			end if
		end if
		return
		end

C***********************************************************************
		subroutine default_thick( thick_string, th_len, thick)
		implicit none
		real*4 rval, thick
		integer*4 th_len
		character*(*) thick_string
		if( th_len .gt. 0) then
			call convert( thick_string, th_len, rval)
			if( th_len .lt. 0 .or. rval .lt. 0.) then
			    thick= 100.
			else
			    thick= rval
			end if
		end if
		return
		end

C**************************************************
		subroutine convert( str, ilen, val)
		integer*4 ilen
		character*(*) str
		real*4 val
		read( str,*,err=99) val
		return
99		ilen = -1
		return
		end
C*****************************************************************
		integer*4 function width( var)
		implicit none
		real*4 var
		if( var .ge. 1000000000.) then
			width= 11
		else if( var .ge. 100000000.) then
			width= 10
		else if( var .ge. 10000000.) then
			width= 9
		else if( var .ge. 1000000.) then
			width= 8
		else if( var .ge. 100000.) then
			width= 7
		else if( var .ge. 10000.) then
			width= 6
		else if( var .ge. 1000.) then
			width= 5
		else if( var .ge. 100.) then
			width= 4
		else if( var .ge. 10.) then
			width= 3
		else
			width= 2
		end if
		return
		end

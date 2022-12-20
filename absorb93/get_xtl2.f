		subroutine abs$get_crystal2_data( nchrs, crys_name, ibasis, natoms
	1			, ized, fract, a_zero, b_zero, c_zero, nstiff, stiff, noread)

		implicit none
		
		integer*4 nchrs, i, j
		integer*4 ibasis, nsites(20), natoms(20), ized(20,5), noread, nstiff
		real*4 fract(20,5), atomx(20,100), atomy(20,100), atomz(20,100)
	1			, a_zero, b_zero, c_zero, alp
	1			, bet, gam, t_debye, temp_coef, stiff(10)
		character*80 input
		character*(*) crys_name
!-- executable code

		open(unit=19, file=crys_name(1:nchrs),status='old',readonly,err=101)
		read(19,'(a)',end=999,err=999) input
		read(19,*,end=999,err=999)  a_zero, b_zero, c_zero, alp
	1			, bet, gam, t_debye, temp_coef
		read(19,*,end=999,err=999) nstiff, (stiff(i),i=1,nstiff)
		read(19,*,end=999,err=999) ibasis, (natoms(i),i=1,ibasis)
		do i=1,ibasis
			read(19,'(a)',end=999,err=999) input
			read(19,*,end=999,err=999) 
	1					(ized(i,j), fract(i,j), j=1,natoms(i))
			read(19,*,end=999,err=999) nsites(i)
			do j=1, nsites(i)
				read(19,*,end=999,err=999) atomx(i,j)
	1					, atomy(i,j), atomz(i,j)
			end do
		end do
		close(unit=19)
		noread = 0
		return
101		noread = 1
		return
999		stop 'Error reading lattice file'
		end

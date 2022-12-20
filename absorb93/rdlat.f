		subroutine abs$rdlattice( crys_name, nchrs, noread)
		implicit none
		common /atom_pos/ ibasis, nsites(20), natoms(20), ized(20,5)
	1		, fract(20,5), atomx(20,100), atomy(20,100), atomz(20,100)
		common /crystal/ a_zero, b_zero, c_zero, alp, bet, gam
	1		, t_debye, temp_coef
		common /stiffness/ nstiff, stiff(10)
		common/ cromer_data/ erg
		integer*4 nchrs, i, j
		integer*4 ibasis, nsites, natoms, ized, noread, nstiff
		real*4 fract, atomx, atomy, atomz, a_zero, b_zero, c_zero, alp
	1		, bet, gam, t_debye, temp_coef, stiff, erg
		character*80 input
		character*(*) crys_name
!-- executable code

		open(unit=17,file=crys_name(1:nchrs)//'.lat',status='old'
	1	,readonly,err=101)
		read(17,'(a)',end=999,err=999) input
		read(17,*,end=999,err=999)  a_zero, b_zero, c_zero, alp
	1		, bet, gam, t_debye, temp_coef
		read(17,*,end=999,err=999) nstiff, (stiff(i),i=1,nstiff)
		read(17,*,end=999,err=999) ibasis, (natoms(i),i=1,ibasis)
		do i=1,ibasis
		    read(17,'(a)',end=999,err=999) input
		    read(17,*,end=999,err=999) 
	1			(ized(i,j), fract(i,j), j=1,natoms(i))
		    read(17,*,end=999,err=999) nsites(i)
		    do j=1, nsites(i)
				read(17,*,end=999,err=999) atomx(i,j)
	1			, atomy(i,j), atomz(i,j)
		    end do
		end do
		erg = 0.	! zero energy for new data set in struct_fact
		close(unit=17)
		noread = 0
		return
101		noread = 1
		return
999		stop 'Error reading Crystal data file'
		end

		subroutine abs$get_crystal_data( jbasis, jatoms, jzed, rfract
	1			, ra_zero, rb_zero, rc_zero, jstiff, rstiff)

		implicit none

		common /atom_pos/ ibasis, nsites(20), natoms(20), ized(20,5)
	1			, fract(20,5), atomx(20,100), atomy(20,100), atomz(20,100)
		common /crystal/ a_zero, b_zero, c_zero, alp, bet, gam
	1			, t_debye, temp_coef
		common /stiffness/ nstiff, stiff(10)

		integer*4 i, j, ibasis, nsites, natoms, ized, nstiff
		real*4 fract, atomx, atomy, atomz, a_zero, b_zero, c_zero, alp
	1			, bet, gam, t_debye, temp_coef, stiff

		real*4 rfract(20,5), ra_zero, rb_zero, rc_zero, rstiff(10)
		integer*4 jbasis, jstiff, jatoms(20), jzed(20,5), jsites, ratoms
!-- executable code

		ra_zero= a_zero
		rb_zero= b_zero
		rc_zero= c_zero

		jbasis= ibasis
		do i=1,ibasis
			jatoms(i)= natoms(i)
			do j=1,natoms(i)
				jzed(i,j)= ized(i,j)
				rfract(i,j)= fract(i,j)
			end do
		end do

		do i=1,nstiff
			rstiff(i)= stiff(i)
		end do
		return
		end

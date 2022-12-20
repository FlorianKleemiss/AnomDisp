		subroutine abs$set_crystal_data( jatoms, jzed, rfract
	1			, ra_zero, rb_zero, rc_zero)

		implicit none

		common /atom_pos/ ibasis, nsites(20), natoms(20), ized(20,5)
	1			, fract(20,5), atomx(20,100), atomy(20,100), atomz(20,100)
		common /crystal/ a_zero, b_zero, c_zero, alp, bet, gam
	1			, t_debye, temp_coef
		common /stiffness/ nstiff, stiff(10)

		integer*4 i, j, ibasis, nsites, natoms, ized, nstiff
		real*4 fract, atomx, atomy, atomz, a_zero, b_zero, c_zero, alp
	1			, bet, gam, t_debye, temp_coef, stiff, erg

		real*4 rfract(20,5), ra_zero, rb_zero, rc_zero
		integer*4 jatoms(20), jzed(20,5)
!-- executable code

		a_zero= ra_zero
		b_zero= rb_zero
		c_zero= rc_zero

		do i=1,ibasis
			natoms(i)= jatoms(i)
			do j=1,natoms(i)
				ized(i,j)= jzed(i,j)
				fract(i,j)= rfract(i,j)
			end do
		end do
		return
		end

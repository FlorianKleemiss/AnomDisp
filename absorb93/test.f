
       include 'cromer.f'
       PROGRAM Test
         
         
         real :: energy
         real*4 :: fp, fpp
         integer :: el
		 
		 fp = 0.0
		 fpp = 0.0
         
C		 write(*,*) 'Please enter energy in eV'
C		 read (*,*) energy
C        write(*,*) 'Please enter element number'
C		 read (*,*) el

		 energy = 19999
		 el = 42
         
         call abs_cromer(el,energy,fp,fpp)
         write(*,*) 'fprime: ',fp
         write(*,*) 'fdouble-prime: ',fpp
       end program Test
  

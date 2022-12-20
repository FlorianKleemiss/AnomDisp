C ABS$GET_STRING
C INCLUDES
C    ABS$UP_LOW
C    ABS$LOW_UP
C 5-MAY-1993 SMB

		subroutine abs$get_string( string, nchrs)
		integer*4 nchrs
		character*(*) string
		read(*,'(q,a)',end=100,err=100) nchrs, string
		if( nchrs .gt. 0) call abs$up_low( string, nchrs)
		return
100		nchrs= -1
		return
		end
C************************************************************************
C convert upper case to lower case
		subroutine abs$up_low( string, len)
		implicit none
		character*(*) string
		character*26 lower/'abcdefghijklmnopqrstuvwxyz'/
		character*26 upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
		integer*4 len, j, i

C executable code
        if( len .gt. 0) then
		    do i=1,len
			    j= index( upper,string(i:i))
			    if( j.gt. 0) string(i:i)= lower(j:j)
		    end do
		end if
		return
		end
C************************************************************************
C convert lower case to upper case
		subroutine abs$low_up( string, len)
		implicit none
		character*(*) string
		character*26 lower/'abcdefghijklmnopqrstuvwxyz'/
		character*26 upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
		integer*4 len, j, i

C executable code
        if( len .gt. 0) then
		    do i=1,len
			    j= index( lower,string(i:i))
			    if( j.gt. 0) string(i:i)= upper(j:j)
		    end do
		end if
		return
		end

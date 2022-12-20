C ABS$ELEMENT.FOR
C 5-MAY-1993 SMB
C********************************************************************
C  this routine returns an integer atomic number for a symbol name.
C  if on entry to the program zed<0 then it compares the value of symbol
C  passed through from the main program.  If zed=0 it prompts for
C  a value of symbol to parse. If 1<zed<92 it returns that symbol.
C  Ctrl/z on input returns zed=0
C  
C*********************************************************************
		subroutine abs$element( zed, symbol)
		implicit none
		character*2 element(92), symbol
		integer*4 zed, i, nchr

        data element/'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     1           'na','mg','al','si','p ','s ','cl','ar',
     2  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     3       'ga','ge','as','se','br','kr',
     4  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag','cd',
     5       'in','sn','sb','te','i ','xe',
     6  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     7       'ho','er','tm','yb','lu',
     8       'hf','ta','w ','re','os','ir','pt','au','hg',
     9       'tl','pb','bi','po','at','rn',
     1           'fr','ra','ac','th','pa','u '/
C executable code
C
		if( zed .eq. 0) then
			write(*,'(''$Enter atomic symbol (e.g. "Si")=> '')')
			call abs$get_string( symbol, nchr)
			if( nchr .le. 0) then
				zed= 0
				return
			end if
		else if( zed .ge. 1 .and. zed .le. 92) then
C if it gets to here, number was passed, wants name

			symbol= element( zed)
			return
		else
C they must have passed in a symbol to be tested (zed<0)

		    call abs$up_low( symbol, 2)
		
		end if

		do i=1,92		! loop through for match
			if( symbol.eq. element( i)) then
				zed = i
				call abs$low_up(symbol(1:1),1)
				return
			end if
		end do
		write(*,'('' Symbol not recognized'')')
		zed = 0
		return
		end

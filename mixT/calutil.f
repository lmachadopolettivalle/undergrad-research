      subroutine fxwrite (string, unit)
      implicit none
      character*(*) string
      integer unit
      integer lnblnk
      write (unit,'(a)') string(1:lnblnk(string))
      return
      end

      integer function clenact (string)
      implicit none
      character*(*) string
      clenact = len(string)
      return
      end

      subroutine crmvblk (string)
      implicit none
      character*(*) string

      call rmblanks(string)
      return
      end

      subroutine crmvlbk (string)
      implicit none
      character*(*) string

      call rmblanks(string)
      return
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
*+WTINFO
        subroutine wtinfo(chatter, wtchatter, level, string)

        IMPLICIT NONE
        integer level, chatter, wtchatter
        character*(*) string
c 
c Description:
c  Noddy routine which passes string and level onto the standard wtout 
c subroutine if chatter.GE.wtchatter
c
c Passed parameters
c  CHATTER       i   : (int) actual chatter flag
c  WTCHATTER     i   : (int) chatter flag at & or above which wtout called
c  LEVEL         i   : (int) importance level (see above)
c  STRING        i   : (char) Context string to be appended to standard msg
c
c Called Routines:
c  subroutine WTOUT      : (CALLIB) Standard roslib/callib string writer
c
c Compilation & Linking
c  link with CALLIB 
c
c Origin:
c  Original
c
c Authors/Modification History:
c  Ian M George     (1.0.0: 1995 Nov 29) original
c       character*7 version
c       parameter (version = '1.0.0')
*- 

        if(chatter.ge.wtchatter) then
                call wtout(level, string)
        endif

        return
        end
c -------------------------------------------------------------------

*+WTOUT
        subroutine wtout(level, string)

	IMPLICIT NONE
	integer level
        character*(*) string
c 
c Description:
c  Writes callib/roslib standard message string(s) to STDOUT, with the 
c format controled by the level parameter:
c 	level = 0 	high level message 
c       level = 1       moderate level message
c	level = 2 	low level message
c
c Passed parameters
c  LEVEL         i   : (int) importance level (see above)
c  STRING	 i   : (char) Context string to be appended to standard msg
c
c Called Routines:
c  function CLENACT      : (CALLIB) Returns actual length of string
c  subroutine CRMVLBK    : (CALLIB) Removes leading blanks from a string
c  subroutine FCECHO     : (FTOOLS) Writes to standard o/p device
c
c Compilation & Linking
c  link with FITSIO & CALLIB & FTOOLS
c
c Origin:
c  Original
c
c Authors/Modification History:
c  Ian M George     (1.0.0: 1995 Nov 29) original
c  Ian M George     (1.1.0: 1996 Feb 06) lobbed in the islop stuff
c  Ian M George     (1.2.0: 1996 Sep 17) removed calls to crmvlbk to prevent
c                                        problems under Solaris 2.2
c	character*7 version
c	parameter (version = '1.2.0')
*- 
c Internals 
	integer istart, istop, strlen
	integer clenact, ifront, i, islop
	character*11 front, blank
	character*80 outstr

c Initialize
	blank = '           '
	strlen = 0	
	istart = 1
	istop = 0
	islop = 10

c Remove all leading blanks from i/p string & check-out size
c	call crmvlbk(string)
	strlen = clenact(string)

c Sort out the level
	if(level.eq.0) then
		front = ' '
		ifront = 1
	elseif(level.eq.1) then
		front = ' ... '
		ifront = 5
	elseif(level.eq.2) then
		front = ' ...... '
		ifront = 8
	else 
		front = ' ......... '
		ifront = 11
	endif

c Dump the first line
	istop = MIN(strlen,80 - ifront)
	if((istop.LT.strlen).and.
     &	   (string(istop+1:istop+1).NE.' ')) then
	   do i = istop,istop-islop,-1
	     if(string(i:i).EQ.' ') then
		istop = i
		goto 122
	     endif	
	   enddo
	endif
122	outstr = front(:ifront)//string(1:istop)
	call fcecho(outstr)

c Return if we've finished
123	if(istop.GE.strlen) then
		return
	else
125		istart = MAX(1,istop+1)
	        if((string(istart:istart).EQ.' ').and.
     &		  (istart+1.LT.strlen)) then
		  istop = istop + 1
		  goto 125
	        endif
		istop = MIN(strlen,istop + 80 - 9)
		if((istop.LT.strlen).and.
     &	   	  (string(istop+1:istop+1).NE.' ')) then
	   	  do i = istop,istop-islop,-1
	     	     if(string(i:i).EQ.' ') then
			istop = i
			goto 124
	     	     endif	
	   	  enddo
	        endif
124		outstr = blank(:ifront)//string(istart:istop)
		call fcecho(outstr)
		go to 123
	endif

	end
c -------------------------------------------------------------------
C******************************************************************************
C SUBROUTINE:
C      fcecho
C
C DESCRIPTION:
C      This subroutine provides a single point to send text to the
C      terminal. This routine should be modified when a new host
C      environment is used.
C
C AUTHOR/DATE:
C      Kent Blackburn  11/5/91
C
C MODIFICATION HISTORY:
C
C  11/28/94 EAG call fxwrite, which in turn calls umsput
C
C NOTES:
C      fcecho uses F77/VOS like calls for terminal I/O
C
C USAGE:
C      call fcecho(string)
C
C ARGUMENTS:
C      string - text string sent to terminal
C
C PRIMARY LOCAL VARIABLES:
C
C CALLED ROUTINES:
C      subroutine umsput - put message
C
C******************************************************************************
      subroutine fcecho(string)

      character*(*) string
      integer dest,prio,irafsts

      dest = 1
      prio = 0
      irafsts = 0

c write to STDOUT and logfile
      call fxwrite (string, 6)
c      call logstr(string)
c      call umsput(string,dest,prio,irafsts)
      return
      end

*+WTERRM
        subroutine wterrm(subrout, version, string)

        IMPLICIT NONE
        character*(*) subrout, version, string
c 
c Description:
c  Writes callib/roslib standard error message string(s) to STDOUT.
c  Current format is approximately:
c       ' ERROR - '//subrout//version//string
c  but with a few extra bits of punctuation, plus some attempt to 
c  handle strings greater than 80 characters.
c
c Passed parameters
c  SUBROUT       i   : (char) Name of the subroutine from which wterr called
c  VERSION       i   : (char) Version of SUBROUT
c  STRING        i   : (char) Context string to be appended to standard msg
c
c Called Routines:
c  function CLENACT      : (CALLIB) Returns actual length of string
c  subroutine CRMVLBK    : (CALLIB) Removes leading blanks from a string
c  subroutine FCECHO     : (FTOOLS) Writes to standard o/p device
c
c Compilation & Linking
c  link with FITSIO & CALLIB & FTOOLS
c
c Origin:
c  Original
c
c Authors/Modification History:
c  Ian M George     (1.0.0: 1995 Nov 29) original
c  Keith Arnaud     (1.1.0: 1996 Aug 21) added tsubrout, tversion, tstring
c                                        internal variables to prevent probs
c  					 under Solaris 2.2 when i/p string 
c    					 starts with leading spaces 
c       character*7 version
c       parameter (version = '1.1.0')
*- 
c Internals 
        integer sublen, verlen, strlen
        integer istart, istop, str1len
        integer clenact
        character*80 outstr
        character*255 tsubrout, tversion, tstring

C Copy the passed parameters into temporary strings to avoid any problems
C with modifying constant strings

        tsubrout = subrout(:MIN(len(tsubrout),len(subrout)))
        tversion = version(:MIN(len(tversion),len(version)))
        tstring  = string(:MIN(len(tstring),len(string)))

c Initialize
        sublen = 0
        verlen = 0
        strlen = 0      
        istart = 1
        istop = 0
        str1len = 0

c Remove all leading blanks from i/p strings
        call crmvlbk(tsubrout)
        call crmvlbk(tversion)
        call crmvlbk(tstring)

c Check out the size of each character string
        sublen = clenact(tsubrout)
        verlen = clenact(tversion)
        strlen = clenact(tstring)

c Work out the 1st bit of the first line (up to where string will begin)        
        outstr = ' ERROR - '//tsubrout(:sublen)//
     &          ' '//tversion(:verlen)//': '
        str1len = clenact(outstr) + 1

c Dump the first line
        istop = MIN(strlen,80 - str1len)
        outstr = outstr(:str1len)//tstring(1:istop)
        call fcecho(outstr)

c Return if we've finished
123     if(istop.GE.strlen) then
                return
        else
                istart = MAX(1,istop+1)
                istop = MIN(strlen,istop + 80 - 9)
                outstr = '         '//tstring(istart:istop)
                call fcecho(outstr)
                go to 123
        endif

        end
c -------------------------------------------------------------------



*+WTFERR
        subroutine wtferr(subrout, version, status, string)

	IMPLICIT NONE
        character*(*) subrout, version, string
	integer status
c 
c Description:
c  Writes callib/roslib standard error message string(s) to STDOUT, 
c  appending the FITSIO error message corresponding to the code 
c  passed down via status
c
c Passed parameters
c  SUBROUT       i   : (char) Name of the subroutine from which wterr called
c  VERSION       i   : (char) Version of SUBROUT
c  STATUS        i   : (int) FITSIO status flag
c  STRING	 i   : (char) Context string to be appended to standard msg
c
c Called Routines:
c  subroutine FTGERR     : (FITSIO) Gets FITSIO error string from status
c  subroutine FTVERS     : (FITSIO) Gets FITSIO version number
c  subroutine WTINFO     : (CALLIB) Writes callib/roslib info message
c  subroutine WTERRM     : (CALLIB) Writes callib/roslib error message
c
c Compilation & Linking
c  link with FITSIO & CALLIB & FTOOLS
c
c Origin:
c  Original
c
c Authors/Modification History:
c  Ian M George     (1.0.0: 1995 Nov 29) original
c	character*7 version
c	parameter (version = '1.0.0')
*- 
c Internals 
	character*80 ftsmsg
	character*160 message
	character*6 cftsver	
	real ftsver

c Return if there's no error
	if(status.EQ.0) return

c Dump the user-defined error message
	call wterrm(subrout, version, string)
	
c Get & dump the FITSIO message & version
	call ftgerr(status, ftsmsg)
	call ftvers(ftsver)	
	write(cftsver,'(f6.3)') ftsver
	message = 'fitsio'//cftsver//' error message: '//ftsmsg
	call wtinfo(1,1,1,message)

	return
	end
c -------------------------------------------------------------------
*+WTFWRN
        subroutine wtfwrn(subrout, version, chatter, wtchatter, 
     &		status, string)

	IMPLICIT NONE
        character*(*) subrout, version, string
	integer status, chatter, wtchatter
c 
c Description:
c  Writes callib/roslib standard warning message string(s) to STDOUT, 
c  appending the FITSIO error message corresponding to the code 
c  passed down via status
c
c Passed parameters
c  SUBROUT       i   : (char) Name of the subroutine from which wterr called
c  VERSION       i   : (char) Version of SUBROUT
c  CHATTER       i   : (int) actual chatter flag
c  WTCHATTER     i   : (int) chatter level at or above which messages written
c  STATUS        i   : (int) FITSIO status flag
c  STRING	 i   : (char) Context string to be appended to standard msg
c
c Called Routines:
c  subroutine FTGERR     : (FITSIO) Gets FITSIO error string from status
c  subroutine FTVERS     : (FITSIO) Gets FITSIO version number
c  subroutine WTINFO     : (CALLIB) Writes callib/roslib info message
c  subroutine WTWARM     : (CALLIB) Writes callib/roslib warning message
c
c Compilation & Linking
c  link with FITSIO & CALLIB & FTOOLS
c
c Origin:
c  Original
c
c Authors/Modification History:
c  Ian M George     (1.0.0: 1995 Nov 29) original
c	character*7 version
c	parameter (version = '1.0.0')
*- 
c Internals 
	character*80 ftsmsg
	character*160 message
	character*6 cftsver	
	real ftsver

c Return if there's no error
	if(status.EQ.0) return

c Return if running in silent mode
	if(chatter.LT.wtchatter) return

c Dump the user-defined error message
	call wtwarm(subrout, version, chatter, wtchatter, string)
	
c Get & dump the FITSIO message & version
	call ftgerr(status, ftsmsg)
	call ftvers(ftsver)	
	write(cftsver,'(f6.3)') ftsver
	message = 'fitsio'//cftsver//' error message: '//ftsmsg
	call wtinfo(chatter,wtchatter,1,message)

	return
	end
c -------------------------------------------------------------------
*+WTWARM
        subroutine wtwarm(subrout, version, chatter, wtchatter, string)

	IMPLICIT NONE
	integer chatter, wtchatter
        character*(*) subrout, version, string
c 
c Description:
c  Writes callib/roslib standard warning message string(s) to STDOUT.
c  Current format is approximately:
c  	' WARNING - '//subrout//version//string
c  but with a few extra bits of punctuation, plus some attempt to 
c  handle strings greater than 80 characters.
c
c Passed parameters
c  SUBROUT       i   : (char) Name of the subroutine from which wterr called
c  VERSION       i   : (char) Version of SUBROUT
c  CHATTER       i   : (int) chatter flag (nowt written if chatter = 0)
c  WTCHATTER     i   : (int) chatter flag at or above which string written
c  STRING	 i   : (char) Context string to be appended to standard msg
c
c Called Routines:
c  function CLENACT      : (CALLIB) Returns actual length of string
c  subroutine CRMVLBK    : (CALLIB) Removes leading blanks from a string
c  subroutine FCECHO     : (FTOOLS) Writes to standard o/p device
c
c Compilation & Linking
c  link with FITSIO & CALLIB & FTOOLS
c
c Origin:
c  Original, based on wterrm (1.0.0)
c
c Authors/Modification History:
c  Ian M George     (1.0.0: 1995 Nov 29) original
c  Ian M George     (1.1.0: 1996 Sep 17) removed calls to crmvlbk to prevent 
c					 problems under Solaris 2.2
c	character*7 version
c	parameter (version = '1.1.0')
*- 
c Internals 
	integer sublen, verlen, strlen
	integer istart, istop, str1len
	integer clenact
	character*80 outstr

c Initialize
	sublen = 0
	verlen = 0
	strlen = 0	
	istart = 1
	istop = 0
	str1len = 0

c Return if chatter flag is zero
	if(chatter.LT.wtchatter) return

c Remove all leading blanks from i/p strings
c	call crmvlbk(subrout)
c	call crmvlbk(version)
c	call crmvlbk(string)

c Check out the size of each character string
	sublen = clenact(subrout)
	verlen = clenact(version)
	strlen = clenact(string)

c Work out the 1st bit of the first line (up to where string will begin)	
	outstr = ' WARNING - '//subrout(:sublen)//
     &		' '//version(:verlen)//': '
	str1len = clenact(outstr) + 1

c Dump the first line
	istop = MIN(strlen,80 - str1len)
	outstr = outstr(:str1len)//string(1:istop)
	call fcecho(outstr)

c Return if we've finished
123	if(istop.GE.strlen) then
		return
	else
		istart = MAX(1,istop+1)
		istop = MIN(strlen,istop + 80 - 11)
		outstr = '           '//string(istart:istop)
		call fcecho(outstr)
		go to 123
	endif

	end
c -------------------------------------------------------------------

      subroutine rmblanks (s)
c
c Remove blanks from a string
c
      implicit none
      character s*(*)
      
      integer lnblnk,len,i,n
      
      
      len=lnblnk(s)
      
      n=0
      do i=1,len
        if(s(i:i).ne.' '.and.s(i:i).ne.'\t')then
          n=n+1
          s(n:n)=s(i:i)
        endif
      enddo
      
      do i=n+1,len
        s(i:i)=' '
      enddo
      
      return
      end
      
      subroutine get_parameter_value (parname,parvalue,type)
      implicit none
      character parname*(*)
      integer parvalue
      character type*(*)
      

      if (type.eq.'j'.or.type.eq.'J') then
        call get_parameter_value_j (parname,parvalue)
      else if (type.eq.'e'.or.type.eq.'E') then
        call get_parameter_value_e (parname,parvalue)
      else if (type.eq.'d'.or.type.eq.'D') then
        call get_parameter_value_d (parname,parvalue)
      else
        write (0,*) ' Use type = j, e, or d in get_parameter_value'
        call exit (1)
      endif
      
      return
      end


*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine get_parameter_value_j (parname,parvalue)
      implicit none
      
      character parname*(*)
      integer parvalue
      character arg*80
      logical defined
      integer lnblnk,status

      call get_cl_par (parname,arg)
      if (.not.defined(arg)) then
        arg = parname
        write (0,*) arg(1:lnblnk(arg))//'=?'
        call exit(1)
      endif

      status=0
      read (arg,*,iostat=status) parvalue
      if (status.ne.0) then
        write (0,*) 'Bad value of parameter ',parname(1:lnblnk(parname)),': ',
     ~      arg(1:lnblnk(arg))
        write (0,*) ' (attempt to read as an integer) '
        call exit(1)
      endif

      return
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine get_parameter_value_e (parname,parvalue)
      implicit none
      
      character parname*(*)
      real parvalue
      character arg*80
      logical defined
      integer lnblnk,status

      call get_cl_par (parname,arg)
      if (.not.defined(arg)) then
        arg = parname
        write (0,*) arg(1:lnblnk(arg))//'=?'
        call exit(1)
      endif

      status=0
      read (arg,*,iostat=status) parvalue
      if (status.ne.0) then
        write (0,*) 'Bad value of parameter ',parname(1:lnblnk(parname)),': ',
     ~      arg(1:lnblnk(arg))
        write (0,*) ' (attempt to read as a real) '
        call exit(1)
      endif

      return
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine get_parameter_value_d (parname,parvalue)
      implicit none
      
      character parname*(*)
      double precision parvalue
      character arg*80
      logical defined
      integer lnblnk,status

      call get_cl_par (parname,arg)
      if (.not.defined(arg)) then
        arg = parname
        write (0,*) arg(1:lnblnk(arg))//'=?'
        call exit(1)
      endif

      status=0
      read (arg,*,iostat=status) parvalue
      if (status.ne.0) then
        write (0,*) 'Bad value of parameter ',parname(1:lnblnk(parname)),': ',
     ~      arg(1:lnblnk(arg))
        write (0,*) ' (attempt to read as a double precision) '
        call exit(1)
      endif

      return
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine get_parameter_value_s (parname,parvalue)
      implicit none
      
      character parname*(*)
      character parvalue*(*)
      character arg*80
      logical defined
      integer lnblnk,status

      call get_cl_par (parname,arg)
      if (.not.defined(arg)) then
        arg = parname
        write (0,*) arg(1:lnblnk(arg))//'=?'
        call exit(1)
      endif

      parvalue=arg

      return
      end

c Print a error message to stderr and exit(1)
      subroutine exiterror (message)
      character message*(*)

      write (0,*) message
      call exit(1)
      end

c-
      subroutine exit_fitsio (message,status)
      implicit none
      character message*(*)
      integer status

      call perror_fitsio (message,status)
      call exit(1)
      return
      end
       
      
      subroutine perror_fitsio (message,status)
      implicit none
      character message*(*)
      integer lnblnk
      integer status
      character errtext*80
      
      if (status.ne.0) then
        call ftgerr(status,errtext)
        write(0,'(a,a,a)')message(1:lnblnk(message)),': '
     ~      ,errtext(1:lnblnk(errtext))
      endif
      return
      end

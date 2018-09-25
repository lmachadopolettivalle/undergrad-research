*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
      subroutine lin_interpolate (x,y,n,x0,y0)
      implicit none
      integer n
      real x(n),y(n)
      real x0,y0
      
      
      integer ibin
      save ibin
      
      if(x0.le.x(1))then
        ibin=1
      else if (x0.ge.x(n)) then
        ibin=n-1
      else
        if (ibin.le.1.or.ibin.gt.n-1) then
          ibin=n/2
        endif
        call hunt (x,n-1,x0,ibin)
        if (ibin.lt.1.or.ibin.gt.n-1) then
          pause 'unexpected error in lin_interpolate'
        endif
      endif
      
      if (x(ibin+1).ne.x(ibin)) then
        y0=y(ibin)+(y(ibin+1)-y(ibin))*(x0-x(ibin))/(x(ibin+1)-x(ibin))
      else
        y0=y(ibin)
      endif

      return
      end

      subroutine hunt(xx,n,x,jlo)
c
c      
c Given an array xx(n) and given value of x returns jlo such as x 
c between xx(jlo) and xx(jlo+1). xx(n) must be monotonic, either increasing 
c or decreasing. jlo=0 or jlo=n return to indicate that x out of range. 
c jlo on input is taken as initial guess for jlo on output.
c
c
c
      implicit none
      integer n,jlo
      real xx(n),x
      logical ascnd
      integer jhi,inc,jm
      
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
 1      jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
      else
        jhi=jlo
 2      jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
 3    if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      subroutine readfile (n,
     ~    v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,
     ~    v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,
     ~    v21,v22,v23,v24,v25,v26,v27,v28,v29,v30)
      implicit none
      integer n

      real v1(*),v2(*),v3(*),v4(*),v5(*),v6(*),v7(*),v8(*),v9(*),v10(*)
      real v11(*),v12(*),v13(*),v14(*),v15(*),v16(*),v17(*),v18(*),v19(*),v20(*)
      real v21(*),v22(*),v23(*),v24(*),v25(*),v26(*),v27(*),v28(*),v29(*),v30(*)
      
      integer icolmax
      
      character*80 word(1024)
      character*1024 line

      integer icol(30), ncol, firstraw, lastraw
      character*128 filename
      common /creadfile__/ icol,ncol,filename,firstraw, lastraw

      integer unit,newunit,status
      integer lnblnk
      integer i,nwords
      real value
      integer iline
      logical closeit

      icolmax = 0
      do i=1,ncol
        icolmax = max(icol(i),icolmax)
      enddo

      if (filename.eq.'-') then
        unit = 5                ! stdin
        closeit = .false.
      else
        unit = newunit()
        open (unit,file=filename,status='old')
        closeit = .true.
      endif
      n = 1
      status = 0
      iline = 0
      do while (status.eq.0)
        read (unit,'(a)',iostat=status) line
        if (status.eq.0) then
          iline = iline + 1
          if (iline.ge.firstraw.and.iline.le.lastraw) then
            call splitwords (line,word,nwords)
            if (nwords.lt.icolmax) then
              write (0,*) 'readfile: cannot read line ',iline,' of file ',
     ~            filename(1:lnblnk(filename))
              go to 1
            endif
            do i=1,ncol
              read (word(icol(i)),*,iostat=status) value
              if (status.ne.0) then
                write (0,*) 'readfile: cannot read column',icol(i),' in line'
     ~              ,n,' of file ',filename(1:lnblnk(filename))
                go to 1
              endif
              
              if (i.eq.1) then
                v1(n)=value
              else if (i.eq.2) then
                v2(n)=value
              else if (i.eq.3) then
                v3(n)=value
              else if (i.eq.4) then
                v4(n)=value
              else if (i.eq.5) then
                v5(n)=value
              else if (i.eq.6) then
                v6(n)=value
              else if (i.eq.7) then
                v7(n)=value
              else if (i.eq.8) then
                v8(n)=value
              else if (i.eq.9) then
                v9(n)=value
              else if (i.eq.10) then
                v10(n)=value
              else if (i.eq.11) then
                v11(n)=value
              else if (i.eq.12) then
                v12(n)=value
              else if (i.eq.13) then
                v13(n)=value
              else if (i.eq.14) then
                v14(n)=value
              else if (i.eq.15) then
                v15(n)=value
              else if (i.eq.16) then
                v16(n)=value
              else if (i.eq.17) then
                v17(n)=value
              else if (i.eq.18) then
                v18(n)=value
              else if (i.eq.19) then
                v19(n)=value
              else if (i.eq.20) then
                v20(n)=value
              else if (i.eq.21) then
                v21(n)=value
              else if (i.eq.22) then
                v22(n)=value
              else if (i.eq.23) then
                v23(n)=value
              else if (i.eq.24) then
                v24(n)=value
              else if (i.eq.25) then
                v25(n)=value
              else if (i.eq.26) then
                v26(n)=value
              else if (i.eq.27) then
                v27(n)=value
              else if (i.eq.28) then
                v28(n)=value
              else if (i.eq.29) then
                v29(n)=value
              else if (i.eq.30) then
                v30(n)=value
              endif
            enddo
            n = n + 1
          endif
        endif
      enddo
 1    continue
      n = n - 1
      if (closeit) then
        close(unit)
      endif
      return
      end

      
      subroutine datafile (file,columns)
      implicit none
      character*(*) file, columns
      integer icol(30), ncol, firstraw, lastraw
      character*128 filename
      common /creadfile__/ icol,ncol,filename,firstraw,lastraw
      
      character*1024 line
      character*80 word(100)
      integer i
      integer lnblnk

      filename = file
      line = columns
      ncol = lnblnk(line)
      do i=1,ncol
        if (line(i:i).eq.',') line(i:i)=' '
      enddo
      firstraw = 1
      lastraw = 1 000 000 000

      call splitwords (line,word,ncol)
      if (ncol.gt.30) then
        write (0,*) 'Error: readfile can read at most 30 columns at a time'
        call exit(1)
      endif

      do i=1,ncol
        read (word(i),*) icol(i)
      enddo

      return
      end

      subroutine linesfile (i1,i2)
      implicit none
      integer i1,i2
      integer icol(30), ncol, firstraw, lastraw
      character*128 filename
      common /creadfile__/ icol,ncol,filename,firstraw,lastraw
      
      firstraw = i1
      lastraw = i2

      return
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      function newunit ()
      
      logical opened
      
      opened=.true.
      i=6
      do while (opened)
        i=i+1
        inquire(unit=i,opened=opened)
      enddo
      
      newunit=i
      
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine splitwords (string,words,nwords)
      implicit none
      character string*(*)
      character words(*)*(*)
      integer nwords
      integer len,i,j,lnblnk
      
      len=lnblnk(string)
      
      i=1
      
      nwords=0
      
      do while (i.le.len)
        do while ((string(i:i).eq.' '.or.string(i:i).eq.'\t').and.i.lt.len)
          i=i+1
        enddo
        if (i.eq.len) then
          if (string(i:i).ne.' '.and.string(i:i).ne.'\t') then
            nwords=nwords+1
            words(nwords)=string(i:i)
          endif
          goto 100
        endif
        
        nwords=nwords+1
        j=1
        words(nwords)=' '
        do while ((string(i:i).ne.' '.and.string(i:i).ne.'\t').and.i.le.len)
          words(nwords)(j:j)=string(i:i)
          i=i+1
          j=j+1
        enddo
      enddo
      
 100  continue
      
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

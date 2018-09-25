      implicit none
      real T,z,nH

      character*200 arffile, rmffile, filename
      character*80 arg, colname, matext, hduclas3,rmfversn
      character*80 telescop, instrume, detnam, filter, chantype
      integer nsmax,nchmax, ns,nch
      parameter (nsmax=3072,nchmax=4096)
      real eb1(nsmax),eb2(nsmax),rsp(nchmax,nsmax)
      real e_bin(0:nsmax), rsspec_cont(nsmax),rsspec_line(nsmax),
     ~    effarea(nsmax)
      real e_min(nchmax), e_max(nchmax)
      integer channel (nchmax)
      real respline(nchmax)
      integer ngrpmax, grpmax
      parameter (ngrpmax=100)
      integer ngrp(nsmax), F_chan(nsmax,ngrpmax), N_chan(nsmax,ngrpmax)
      real lo_thresh, areascal
      integer flchan, chatter
      double precision modelspec_justline(nsmax)
      double precision modelspec_line(nsmax), rsspec8line(nsmax)
      double precision modelspec_cont(nsmax), rsspec8cont(nsmax)
      double precision phmodelspec_line(nsmax), emodelspec_line(nsmax)
      double precision phmodelspec_cont(nsmax), emodelspec_cont(nsmax)
      integer cmin,cmax
      real emin,emax
      integer i,j
      real sigism,e
      double precision absorp
      character*80 plmodel, abtable ! Plasma and metal abundance models

      double precision fluxcont,fluxline, fluxlinetot,emean,weight,logemean

      integer unitrmf, unitarf, colnum, nrows, blocksize, status
      logical anyf
      
      logical defined

      real logT

      call get_cl_par ('plasma_model',plmodel)
      if (.not.defined(plmodel)) then
        plmodel = 'mekal'
      endif
        
      call get_cl_par ('abundance_table',abtable)
      if (.not.defined(abtable)) then
        abtable = 'angers'
      endif

      call get_parameter_value ('z',z,'e')
      call get_parameter_value ('nH',nH,'e')
      
      call get_cl_par ('rmf',rmffile)
      if (.not.defined(rmffile)) call exiterror ('rmf=?')
      call get_cl_par ('arf',arffile)
      if (.not.defined(arffile)) call exiterror ('arf=?')
      
      call get_cl_par ('chanmin',arg)
      if (defined(arg)) then
        read (arg,*) cmin
      else
        call get_cl_par ('emin',arg)
        if (.not.defined(arg)) call exiterror ('please set chanmin or emin')
        cmin = -1000
        read (arg,*) emin
      endif

      call get_cl_par ('chanmax',arg)
      if (defined(arg)) then
        read (arg,*) cmax
      else
        call get_cl_par ('emax',arg)
        if (.not.defined(arg)) call exiterror ('please set chanmax or emax')
        cmax = -1000
        read (arg,*) emax
      endif
      
      
c A) Read RMF
      call ftgiou (unitrmf,status)
      filename = rmffile
      call ftopen (unitrmf,filename,0,blocksize,status)
      if (status.ne.0) call exit_fitsio (filename,status)
      chatter = 0
      
      matext = 'SPECRESP MATRIX'
      call ftmnhd (unitrmf,-1,matext,0,status)
      if (status.ne.0) then
        status = 0
        matext='MATRIX'
        call ftmnhd (unitrmf,-1,matext,0,status)
      endif
      if (status.ne.0) 
     ~    call exit_fitsio ('SPECRESP MATRIX, '//filename,status)
      
      call rdrmf3(unitrmf, chatter,matext,
     &    telescop, instrume, detnam, filter, areascal,
     &    chantype, flchan, 
     &    nch, ns, eb1, eb2,
     &    grpmax,ngrp,F_chan, N_chan,
     &    rsp,lo_thresh,nchmax,nsmax,
     &    rmfversn,hduclas3,status)
      if (status.ne.0) call exit_fitsio ('MATRIX, '//filename,status)

      call ftmnhd (unitrmf,-1,'EBOUNDS',0,status)
      if (status.ne.0) call exit_fitsio ('EBOUNDS, '//filename,status)

      call rdebd3(unitrmf,chatter,nchmax, 
     &    telescop,instrume,detnam,filter,areascal, 
     &    chantype, flchan,
     &    nch,channel,e_min,e_max,rmfversn,status)
      
      call ftclos (unitrmf,status)
      call ftfiou (unitrmf,status)
      if (status.ne.0) call exit_fitsio (filename,status)

c B) Read ARF
      if (arffile.eq.'none') then
        do i=1,ns
          effarea(i)=1
        enddo
      else
        call ftgiou (unitarf,status)
        filename=arffile
        call ftopen(unitarf,filename,0,blocksize,status)
        if (status.ne.0) call exit_fitsio (filename,status)
        
        matext = 'SPECRESP'
        call ftmnhd (unitarf,-1,matext,0,status)
        if (status.ne.0) call exit_fitsio (filename,status)
        
        colname='SPECRESP'
        call ftgcno (unitarf,.false.,colname,colnum,status)
        if (status.ne.0) call exit_fitsio (filename,status)
        
        call ftgnrw(unitarf,nrows,status)
        if (status.ne.0) call exit_fitsio (filename,status)
        if (nrows.ne.ns) call exiterror
     ~      ('Different number of energy channels in ARF and RMF')
        
        call ftgcve(unitarf,colnum,1,1,ns,0.0,effarea,anyf,status)
        if (status.ne.0) call exit_fitsio (filename,status)

        call ftclos(unitarf,status)
        call ftfiou(unitarf,status)
        if (status.ne.0) call exit_fitsio (filename,status)
      endif

      if (cmin.lt.0) then
        do i=1,nch
          if (e_max(i).lt.emin) then
            cmin = i+1
          endif
        enddo
        write (0,*) 'min chan = ',cmin
      endif

      if (cmax.lt.0) then
        cmax=nch
        do i=nch,1,-1
          if (e_min(i).gt.emax) then
            cmax = i
          endif
        enddo
        write (0,*) 'max chan = ',cmax
      endif


      e_bin(0)=eb1(1)
      do i=1,ns
        e_bin(i)=eb2(i)
      enddo

      call mekal_set_abund (abtable)

      T = 0.1
      logT = log(T)
      do while ( T .lt. 22.0)
        call mekal (T, 0.0, z, ns, e_bin, rsspec_cont)
        call mekal (T, 3.0, z, ns, e_bin, rsspec_line)
        do i=1,ns
          rsspec8cont(i)=dble(rsspec_cont(i))
          rsspec8line(i)=dble(rsspec_line(i))
        enddo
      
        do i=1,ns
          e = 0.5*(eb1(i)+eb2(i))*1000.0
          absorp = sigism(e)*nH
          if (absorp.lt.50.0) then
            absorp=exp(-absorp)
          else
            absorp=0.0
          endif
          modelspec_cont(i) = absorp*rsspec8cont(i)*effarea(i)
          modelspec_line(i) = absorp*rsspec8line(i)*effarea(i)
          modelspec_justline(i) = (modelspec_line(i)-modelspec_cont(i))
          phmodelspec_cont(i)=rsspec8cont(i)
          phmodelspec_line(i)=rsspec8line(i)
          emodelspec_cont(i)=rsspec8cont(i)*e*1.60219e-12
          emodelspec_line(i)=rsspec8line(i)*e*1.60219e-12
        enddo
        
        do i = 1,nch
          respline(i)=0
        enddo
        do i=1,ns
          do j=1,nch
            respline(j)=respline(j)+modelspec_justline(i)*rsp(j,i)
          enddo
        enddo
        emean=0
        logemean=0
        weight=0
        do i=1,nch
          if (0.5*(e_max(i)+e_min(i)).ge.0.25) then
            emean = emean + respline(i)*0.5*(e_max(i)+e_min(i))
            logemean = logemean + respline(i)*log(0.5*(e_max(i)+e_min(i)))
            weight = weight + respline(i)
          endif
        enddo
        emean = emean/weight
        logemean = logemean/weight
        fluxlinetot = weight / 3.0 ! to convert to a=1

        fluxcont = 0
        fluxline = 0
        do i=1,ns
          do j=cmin,cmax
            fluxcont=fluxcont+modelspec_cont(i)*dble(rsp(j,i))
            fluxline=fluxline+modelspec_line(i)*dble(rsp(j,i))
          enddo
        enddo
        fluxline = (fluxline - fluxcont)/3.0 ! to convert to a=1

        fluxcont = fluxcont / 2e-13
        fluxline = fluxline / 2e-13
        fluxlinetot = fluxlinetot / 2e-13
        print '(f7.4,1x,1pe12.6,1x,1pe12.6,1x,0pf7.5,1x,1pe12.6,1x,0pf8.5)',
     ~      T,fluxcont,fluxline,emean,fluxlinetot,logemean

        logT = logT + 0.1*log(2.0d0)
        T = exp(logT)
      enddo
      

      end

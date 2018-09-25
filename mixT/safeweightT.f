      implicit none
      integer ncal
      real tcal(10000), fcont_cal(10000), fline_cal(10000), emean_cal(10000)

      !real t(N), a(N), n2(N)
      real t(3), a(3), n2(3)
      integer nt
      real Tmean,Tline,Tcont

      character*128 datfile

      call get_parameter_value_s ('cal',datfile)


      call datafile (datfile,'1 2 3 4')
      call readfile (ncal,tcal,fcont_cal,fline_cal,emean_cal)

      t(1)=4                  ! keV
      a(1)=0.1                  ! Solar
      n2(1)=0.9                 ! \int n^2 dV
      
      t(2)=1                  ! keV
      a(2)=0.1                  ! Solar
      n2(2)=0.1                 ! \int n^2 dV

      nt = 2

      call calc_Txspec (t,a,n2,nt,
     ~    ncal,tcal,fcont_cal,fline_cal,emean_cal,
     ~    Tmean,Tcont,Tline)

      print*,Tmean

      end


*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine  calc_Txspec (t,a,n2,nt,
     ~    ncal,tcal,fcont_cal,fline_cal,emean_cal,Tmean,Tcont,Tline)
c
c  real t(nt)    --- array of input temperatures
c  real a(nt)    --- array of input metallicities
c  real n2(nt)   --- array of input emission measures
c
c  real tcal(nt), fcont_cal(nt), fline_cal(nt), emean_cal(nt) -- calibration
c                                                                arrays
c  Tmean, Tcont, Tline --- "spectro-T", "continuum-T", "line-T" (see PAPER)
c
      implicit none
      integer nt,ncal
      real t(nt), a(nt), n2(nt)
      real tcal(ncal), fcont_cal(ncal), fline_cal(ncal), emean_cal(ncal)
      real Tmean

      real Tcont, Tline

      integer i

      real wline, wcont, fline, fcont, fluxline, fluxcont, emean
      real x

      real acont, delta1, delta2, beta

      character*80 mode

      call get_cl_par ('mode',mode)
      if (mode.eq.'xmm/pn'.or.mode.eq.'XMM/PN') then
        acont = 0.79
        beta = 0.75
        delta1 = 0.270
        delta2 = 0.225
      elseif (mode.eq.'xmm/mos'.or.mode.eq.'XMM/MOS') then
        acont = 0.90
        beta = 1
        delta1 = 0.19
        delta2 = 0.22
      elseif (mode.eq.'xmm/mos+pn'.or.mode.eq.'XMM/MOS+PN') then
        acont = 0.91
        beta = 0.90
        delta1 = 0.19
        delta2 = 0.21
      elseif (mode.eq.'asca/sis'.or.mode.eq.'ASCA/SIS') then
        acont = 0.875
        beta = 0.80
        delta1 = 0.20
        delta2 = 0.22
      elseif (mode.eq.'asca/gis'.or.mode.eq.'ASCA/GIS') then
        acont = 0.79
        beta = 0.75
        delta1 = 0.26
        delta2 = 0.30
      else ! chandra by default
        acont = 0.875
        beta = 1
        delta1 = 0.19
        delta2 = 0.25
      endif

      

      wline = 0
      wcont = 0
      Tline = 0
      Tcont = 0
      fluxline = 0
      fluxcont = 0
      do i=1,nt

c     Find continuum and line flux, and average energy of the line emission
        call lin_interpolate (tcal,fcont_cal,ncal,t(i),fcont)
        call lin_interpolate (tcal,fline_cal,ncal,t(i),fline)
        call lin_interpolate (tcal,emean_cal,ncal,t(i),emean)

c     Multiply fluxes by emission measure and metallicity
        fcont = fcont*n2(i)
        fline = fline*n2(i)*a(i)
        
c     eq.[6] in the paper
        Tcont = Tcont + t(i)*fcont/t(i)**acont
        wcont = wcont + fcont/t(i)**acont
        fluxcont = fluxcont + fcont
        
c     eq.[2] in the paper
        Tline = Tline + emean*fline
        wline = wline + fline
      enddo

c     eq.[4] in the paper
      Tcont = Tcont/wcont

c     eq.[3] in the paper
      emean = Tline/wline
      call lin_interpolate (emean_cal,tcal,ncal,emean,Tline)
      fluxline = wline
        
c     eq.[7,8,12] in the paper
      x = fluxcont/(fluxcont+fluxline)
      wcont =
     ~    exp(-(((x-1)**2)/delta1**2)**beta) *
     ~    exp(-(((x-1)**2)/delta2**2)**4.0)
      Tmean = wcont*Tcont + (1-wcont)*Tline

      return
      end
      
      

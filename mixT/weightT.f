      implicit none
      integer*4 length
      integer ncal
      real tcal(10000), fcont_cal(10000), fline_cal(10000), emean_cal(10000)

      !real t(nt), a(nt), n2(nt)
      !real t(1000), a(1000), n2(1000)
      real t(21567005), a(21567005), n2(21567005)
      integer*4 L(32)
      integer T_cut(32)
      integer nt,first,haloid
      character(len=2) strhaloid
      character(len=1) half
      real Tmean,Tline,Tcont
      real x,y,z
      real Tcut

      character*128 datfile

      call get_parameter_value_s ('cal',datfile)


      call datafile (datfile,'1 2 3 4')
      call readfile (ncal,tcal,fcont_cal,fline_cal,emean_cal)

      !change both! haloid ranges from 0 to 31
      !half refers to whether we use up to 1*R500 or 0.5*R500
      haloid       = 31
      strhaloid =   "31"
      half = 'Y'



      !first L for high cut; second L for half cut
      if (half .eq. 'N') then
          L (1:32) = (/ 13848376, 13470777, 6134598, 4796732, 4526639, 3537781,
     ~2325429, 2494396, 2044658, 2294799, 1701807, 1462543, 677832,
     ~1177737, 776693, 724039, 725305, 1167112, 842166, 583692,
     ~322131, 278142, 411551, 408215, 256655, 131663, 80799,
     ~89587, 238905, 159285, 92468, 78547 /)
      endif
      if (half .eq. 'Y') then
          L (1:32) = (/ 21567005, 3215756, 3172973, 2567207, 2202625, 1718538,
     ~1247439, 1276157, 1113104, 1133731, 1014163, 798373, 382891,
     ~669145, 630423, 476006, 438317, 489002, 615958, 397880,
     ~273114, 192218, 254677, 306655, 170322, 142511, 113358,
     ~91815, 231426, 91874, 51558, 44029 /)
      endif


      nt = 1
      Tcut = 10.0
      length = L(haloid+1)-10

      first = 1
      if (half .eq. 'N') then
          open (1, file = "./HighCuttxt/dataforweightT_"//strhaloid//
     ~".txt") 
      endif
      if (half .eq. 'Y') then
          open (1, file = "./dataforweightT_"//strhaloid//
     ~".txt") 
      endif
      do
          read (1,100,end=30) x,y,z
100       format (f6.5,X,f7.6,X,f11.10)
          if (x > Tcut) then
              !print*,x
              length = length - 1
          else
              if (nt >= first .and. nt <= length) then
                  t(nt-first+1) = x
                  a(nt-first+1) = y
                  n2(nt-first+1) = z
              end if
20            nt = nt + 1
          end if
      end do
30    close (1) 



      call calc_Txspec (t,a,n2,length,
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
      
      

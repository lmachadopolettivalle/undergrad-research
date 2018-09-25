*      implicit none
*      real T, a, z
*      real ebin(0:10000), spec(10000), gov(10000), spec0(10000)
*      integer ne,i
*
*      T=1.0
*      a=0.0
*      z=0.0
*
*      call datafile ('mekal.test.2','1 2')
*      call linesfile (4,10000)
*      call readfile (ne,ebin(1),gov)
*      do i=1,ne
*        ebin(i-1)=ebin(i)-gov(i)
*      enddo
*      ebin(ne)=ebin(ne)+gov(ne)
*
*      call mekal (T, a, z, ne, ebin, spec)
*      call set_abunds ('xspec')
*      call raym  (T, a, z, ne, ebin, spec0)
**      call mekal_set_abund ('aneb')
**      call mekal (T, a, z, ne, ebin, spec0)
*      do i=1,ne
*        print*,i,ebin(i),spec(i),spec0(i),spec(i)/spec0(i)
*      enddo      
*
*      end


      subroutine mekal (T,abund,z,ne,ebin,spec)
c
c Alexey Vikhlinin <alexey@head-cfa.harvard.edu>  Jul 11 2001
c
c Interpolates over the XSPEC pre-calculated MekaL model grid. The meaning
c of input and output parameters is identical to those in J.Hughes' Raymond
c -Smith code:
c
c
c Input:
c 
c real T -- temperature, keV;
c real abund -- relative heavy element abundance
c real z -- redshift;
c integer ne -- number of the spectrum bins, must be Ns.le.NSMAX;
c real ebin(0:ne) -- energy bin boundaries (note the dimension), receiver
c frame, keV
c
c Output:
c 
c real spec(ne) -- spectrum in units of photon*cm**3/s (into 4pi angle) in a
c given ebin interval (energy boundaries are in the receiver frame). Flux is
c given in emitter frame. It should be multiplied by EM=n_e*n_p*V  to get the
c luminosity.               NOTE THE DEFINITION OF EM~~~~~~~~~~~~
c
c The program requires CFITSIO library to read the XSPEC model file and 
c ReadPar library (included in the ZHTOOLS library) to specify the location
c of the model file [call get_cl_par ('mekal_model_file',mekal_model_file)]
c
c The model file should be specified by setting the mekal_model_file=filename
c parameter in the command line.
c
c If you prefer, the model file can be hardcoded by setting the
c -D'MEKALFILE="/path/to/mekal/model/file"' compiler option; in this case the
c ReadPar library and command line arguments will be not required.
c
c
c Entries mekal_set_abund and mekal_set_abund_table can be used to reset the
c abundance table. Usage
c
c  call mekal_set_abund ('allen')
c  call mekal_set_abund ('angers') [or 'xspec' or 'angr' - Anders E. &
c            Grevesse N. (1989,  Geochimica et Cosmochimica Acta 53, 197)]
c  call mekal_set_abund ('feld') [Feldman U. (1992, Physica Scripta 46, 202)]
c  call mekal_set_abund ('aneb') [Anders E. & Ebihara (1982, Geochimica et
c                                  Cosmochimica Acta 46, 2363)]
c
c  real atable(14) ! abundances of He, C, N, O, Ne, Na, Mg, Al, Si, S, 
c                  ! Ar, Ca, Fe, and Ni relative to hydrogen
c  call mekal_set_abund_table (atable)
c
c Fri Sep 17 02:04:58 2004: implement interpolation of log continuum over 1/T
c
      implicit none
      real T, abund, z
      integer ne
      real ebin(0:ne)
      real spec(ne)

      logical firstcall /.true./
      save firstcall
      integer nemax, ntmax, ne0, nt0
      parameter (nemax=3000, ntmax=50)
      real elow(nemax), ehigh(nemax),
     ~    mlspec(nemax,0:14,ntmax),t0(ntmax)
      real mlspec_lines(nemax,ntmax), mlspec_cont(nemax,ntmax)
      save ne0,nt0,elow,ehigh,mlspec,mlspec_lines,mlspec_cont, t0
      character*200 mekal_model_file
      logical defined

      real abund0(14) ! Solar abundances relative to Angers'
      data abund0 /1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
     ~    1.0, 1.0/
      real abund_angers(14)
      data abund_angers /
     ~    9.77e-2,              ! He
     ~    3.63e-4,              ! C
     ~    1.12e-4,              ! N
     ~    8.51e-4,              ! O
     ~    1.23e-4,              ! Ne
     ~    2.14e-6,              ! Na
     ~    3.80e-5,              ! Mg
     ~    2.95e-6,              ! Al
     ~    3.55e-5,              ! Si
     ~    1.62e-5,              ! S 
     ~    3.63e-6,              ! Ar
     ~    2.29e-6,              ! Ca
     ~    4.68e-5,              ! Fe
     ~    1.78e-6               ! Ni
     ~    /
      save abund0, abund_angers
      
      integer i,it1,it2,j,it,iemin,iemax,ie
      real wt1,wt2,wtc1,wtc2,sum,emin,emax
      real atable(14)

      if (firstcall) then
c Load the mekal model file
#ifdef MEKALFILE
c        mekal_model_file = MEKALFILE
        mekal_model_file = '/home/fas/nagai/etl28/programs/Xrays/\
heasoft-6.21/spectral/modelData/mekal.mod'
#else
         mekal_model_file = '/home/fas/nagai/etl28/programs/Xrays/\
heasoft-6.21/spectral/modelData/mekal.mod'
#endif
        call loadmekalfile (mekal_model_file,elow,ehigh,mlspec,t0
     ~      ,nemax,ntmax,ne0,nt0)
        do i=1,ne0
          do it=1,nt0
            mlspec_cont(i,it) = mlspec(i,0,it)+abund0(1)*mlspec(i,1,it)
            mlspec_cont(i,it) = log(mlspec_cont(i,it))
            sum = 0
            do j=2,14
              sum = sum + abund0(j)*mlspec(i,j,it)
            enddo
            mlspec_lines(i,it) = sum
          enddo
        enddo
        firstcall = .false.
      endif


c Find weights for T-interpolation
      if (T.lt.t0(1)) then
        it1=1
        it2=2
        write (0,*) 'warning: T<',t0(1),', the minimum value for MekaL grid'
      else if (T.ge.t0(nt0)) then
        it1=nt0-1
        it2=nt0
        write (0,*) 'warning: T>',t0(nt0),', the maximum value for MekaL grid'
      else
        do i=1,nt0-1
          if (T.ge.t0(i).and.T.lt.t0(i+1)) then
            it1=i
            it2=i+1
          endif
        enddo
      endif
      wt1 = (t0(it2)-T)/(t0(it2)-t0(it1))
      wt2 = (T-t0(it1))/(t0(it2)-t0(it1))
      wtc1 = (1/T - 1/t0(it2))/(1/t0(it1)-1/t0(it2))
      wtc2 = (1/t0(it1) - 1/T)/(1/t0(it1)-1/t0(it2))

c Main interpolation
      iemin=1
      iemax=0
      do i=1,ne
        emin = ebin(i-1)*(1+z)
        emax = ebin(i)*(1+z)
        do ie=iemin,ne0
          if (emin.ge.elow(ie).and.emin.le.ehigh(ie)) then
            iemin = ie
            goto 101
          endif
        enddo
        ! we are here if the energy bin for emin is not found
        write (0,*) emin,' is outside the MekaL energy grid'
        call exit(1)
 101    continue

        do ie=iemin,ne0
          if (emax.ge.elow(ie).and.emax.le.ehigh(ie)) then
            iemax = ie
            goto 102
          endif
        enddo
        ! we are here if the energy bin for emax is not found
        write (0,*) emax,' is outside the MekaL energy grid'
        call exit(1)
 102    continue

        if (iemin.eq.iemax) then
c        emin and emax are within the same energy bin
          spec(i) = ((emax-emin)/(ehigh(iemin)-elow(iemin)))*
     ~        (
     ~        exp(wtc1*mlspec_cont(iemin,it1)+wtc2*mlspec_cont(iemin,it2)) +
     ~        abund*(wt1*mlspec_lines(iemin,it1)+wt2*mlspec_lines(iemin,it2))
     ~        )
        else
          sum = ((ehigh(iemin)-emin)/(ehigh(iemin)-elow(iemin)))*
     ~        (
     ~        exp(wtc1*mlspec_cont(iemin,it1)+wtc2*mlspec_cont(iemin,it2)) +
     ~        abund*(wt1*mlspec_lines(iemin,it1)+wt2*mlspec_lines(iemin,it2))
     ~        )
          do ie=iemin+1,iemax-1
            sum = sum +
     ~          exp(wtc1*mlspec_cont(ie,it1)+wtc2*mlspec_cont(ie,it2)) +
     ~          abund*(wt1*mlspec_lines(ie,it1)+wt2*mlspec_lines(ie,it2))
          enddo
          sum = sum + 
     ~        ((emax-elow(iemax))/(ehigh(iemax)-elow(iemax)))*
     ~        (
     ~        exp(wtc1*mlspec_cont(iemax,it1)+wtc2*mlspec_cont(iemax,it2)) +
     ~        abund*(wt1*mlspec_lines(iemax,it1)+wt2*mlspec_lines(iemax,it2))
     ~        )
          spec(i)=sum
        endif
      enddo

      do i=1,ne
        spec(i)=spec(i)*1e-14
      enddo

      return

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      entry mekal_set_abund_table (atable)
      do i=1,14
        abund0(i)=atable(i)/abund_angers(i)
      enddo
      
      do i=1,ne0
        do it=1,nt0
          mlspec_cont(i,it) = mlspec(i,0,it)+abund0(1)*mlspec(i,1,it)
          sum = 0
          do j=2,14
            sum = sum + abund0(j)*mlspec(i,j,it)
          enddo
          mlspec_lines(i,it) = sum
        enddo
      enddo
      return

      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine mekal_set_abund (table)
      implicit none
      character*(*) table
      real atable(14)
      real abund_angers(14), abund_allen(14), abund_feld(14), abund_aneb(14)
      data abund_angers /
     ~    9.77e-2,              ! He
     ~    3.63e-4,              ! C
     ~    1.12e-4,              ! N
     ~    8.51e-4,              ! O
     ~    1.23e-4,              ! Ne
     ~    2.14e-6,              ! Na
     ~    3.80e-5,              ! Mg
     ~    2.95e-6,              ! Al
     ~    3.55e-5,              ! Si
     ~    1.62e-5,              ! S 
     ~    3.63e-6,              ! Ar
     ~    2.29e-6,              ! Ca
     ~    4.68e-5,              ! Fe
     ~    1.78e-6               ! Ni
     ~    /
      data abund_allen /
     ~    8.51e-02,             ! He
     ~    3.31e-04,             ! C
     ~    9.12e-05,             ! N
     ~    6.61e-04,             ! O
     ~    8.32e-05,             ! Ne
     ~    2.14e-6,              ! Na - taken from Angers table
     ~    2.63e-05,             ! Mg
     ~    2.95e-6,              ! Al - taken from Angers table
     ~    3.31e-05,             ! Si
     ~    1.58e-05,             ! S
     ~    6.31e-06,             ! Ar
     ~    2.00e-06,             ! Ca
     ~    3.98e-05,             ! Fe
     ~    2.00e-06              ! Ni
     ~    /
      data abund_feld /
     ~    9.77e-2,              ! He
     ~    3.98e-4,              ! C 
     ~    1.00e-4,              ! N 
     ~    8.51e-4,              ! O 
     ~    1.29e-4,              ! Ne
     ~    2.14e-6,              ! Na
     ~    3.80e-5,              ! Mg
     ~    2.95e-6,              ! Al
     ~    3.55e-5,              ! Si
     ~    1.62e-5,              ! S 
     ~    4.47e-6,              ! Ar
     ~    2.29e-6,              ! Ca
     ~    3.24e-5,              ! Fe
     ~    1.78e-6               ! Ni
     ~    /
      data abund_aneb /
     ~    8.01e-2,              ! He
     ~    4.45e-4,              ! C 
     ~    9.12e-5,              ! N 
     ~    7.39e-4,              ! O 
     ~    1.38e-4,              ! Ne
     ~    2.10e-6,              ! Na
     ~    3.95e-5,              ! Mg
     ~    3.12e-6,              ! Al
     ~    3.68e-5,              ! Si
     ~    1.89e-5,              ! S 
     ~    3.82e-6,              ! Ar
     ~    2.25e-6,              ! Ca
     ~    3.31e-5,              ! Fe
     ~    8.27e-5               ! Ni
     ~        /
      save abund_angers, abund_allen, abund_feld, abund_aneb
      integer i

      if (table.eq.'allen') then
        do i=1,14
          atable(i)=abund_allen(i)
        enddo
        call mekal_set_abund_table (atable)
      else if (table.eq.'angers'.or.table.eq.'xspec'.or.table.eq.'angr') then
        do i=1,14
          atable(i)=abund_angers(i)
        enddo
        call mekal_set_abund_table (atable)
      elseif (table.eq.'feld') then
        do i=1,14
          atable(i)=abund_feld(i)
        enddo
        call mekal_set_abund_table (atable)
      elseif (table.eq.'aneb') then
        do i=1,14
          atable(i)=abund_aneb(i)
        enddo
        call mekal_set_abund_table (atable)
      else
        write (0,*) 'warning: unknown abundance table: ',table
        write (0,*) 'warning: abundances not reset'
      endif

      return
      end

*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine loadmekalfile (file, elow, ehigh, spec, t0, nemax, ntmax, ne,
     ~    nt)
      implicit none
      character*(*) file
      integer unit, status, icol, coln(0:14)
      logical anyf
      integer ne, nt
      integer nemax, ntmax
      real elow(nemax), ehigh(nemax), spec(nemax,0:14,ntmax), t0(ntmax)
      integer i
      character*80 name
      character*80 errtext
      integer lnblnk


      status = 0
      call ftgiou (unit,status)
      call ftnopn (unit,file,0,status)

c Read energy extension
      call ftmnhd (unit,2,'ENERGIES',0,status)
      call ftgnrw (unit,ne,status)
      if (ne.gt.nemax) then
        write (0,*) 'increase nemax in loadmekalfile'
        call exit(1)
      endif

      call ftgcno (unit,.false.,'ENERG_LO',icol,status)
      call ftgcve (unit,icol,1,1,ne,0.0,elow,anyf,status)

      call ftgcno (unit,.false.,'ENERG_HI',icol,status)
      call ftgcve (unit,icol,1,1,ne,0.0,ehigh,anyf,status)

c Read spectra extension
      call ftmnhd (unit,2,'SPECTRA',0,status)
      call ftgnrw (unit,nt,status)
      if (nt.gt.ntmax) then
        write (0,*) 'increase ntmax in loadmekalfile'
        call exit(1)
      endif
      
      call ftgcno (unit,.false.,'PARAMVAL',icol,status)
      call ftgcve (unit,icol,1,1,nt,0.0,t0,anyf,status)

      call ftgcno (unit,.false.,'INTPSPEC',coln(0),status)
      do icol=1,14
        write (name,'(''ADDSP'',i3.3)') icol
        call ftgcno (unit,.false.,name,coln(icol),status)
      enddo

      do i=1,nt
        do icol=0,14
          call ftgcve (unit,coln(icol),i,1,ne,0.0,spec(1,icol,i),anyf,status)
        enddo
      enddo
      
c
c The 0-th component seems to represent the continuum, while the 14 others - 
c the contribution of 14 elements:
c       1.........He abundance
c       2.........C     "
c       3.........N     "
c       4.........O     "
c       5.........Ne    "
c       6.........Na    "
c       7.........Mg    "
c       8.........Al    "
c       9.........Si    "
c      10.........S     "
c      11.........Ar    "
c      12.........Ca    "
c      13.........Fe    "
c      14.........Ni    "
c abundances seem to be relative to Solar, from Anders E. & Grevesse N. (1989,
c Geochimica et Cosmochimica Acta 53, 197)
c
c Spectrum is given by the number of photons in the given energy interval,
c with the normalization that of the XSPEC model:
c K      = 10**-14 / (4 pi (D_A*(1+z))**2) Int n_e n_H dV, where D_A is 
c          the angular size distance to the source (cm), n_e is the electron 
c          density (cm**-3), and n_H is the hydrogen density (cm**-3)

      call ftclos(unit,status)
      call ftfiou(unit,status)
      
      if (status.ne.0) then
        call ftgerr(status,errtext)
        write(0,'(a,a,a)') file(1:lnblnk(file)),': '
     ~      ,errtext(1:lnblnk(errtext))
        call exit(1)
      endif

      return
      end

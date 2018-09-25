*From:	47413::MCCAMMON     11-MAY-1992 16:33:40.35
*To:	NSSDC::STPMODELS
*CC:	
*Subj:	sigism.for
*
C----------------------------------------------------------------------
C
      FUNCTION SIGISM(ENERGY)  
C
C     Reference:
C      Robert Morrison and Dan McCammon
C      Ap.J., vol. 270,  p. 119 (1983).
C
C
C     Description: 
C
C     This function implements the approximation of Morrison and
C     McCammon (1983) to the interstellar photoelectric absorbtion
C     cross-section.  ENERGY is in eV and the resultant cross-section
C     is in cm**2/hydrogen atom.  Abundances of other elements relative to
C     hydrogen are appropriate for the interstellar medium in the solar
C     neighborhood (see reference for discussion).
C
C
C     Deficiencies:
C     Works only in the range of energy from 30 eV to 10,000 eV.
C       No bounds checking on energy range.
C
C     Bugs:
C     None known -- please report any problems to authors
C  
C     Authors:
C     Dan McCammon               (47413::MCCAMMON)
C     Dick Edgar                 (47413::EDGAR)
C
C     History:
C     25.1.85 - original (By Dick Edgar)
C     19.1.92  - modified format and output units (McCammon)
C
C     Parameters:
C     ENERGY - photon energy in eV
C
C     Type Definitions:
C      IMPLICIT NONE
C
C     Local variables:
      INTEGER I
C          (index for energy interval)
      REAL E
C          (photon energy in keV)

C     Import:
      REAL ENERGY
C          (energy in eV)

C     Export:
      REAL SIGISM
C          (effective cross section in cm**2/H atom)

C     Local constants:


      REAL EMAX(14)
C          (edge energies where polynomial changes -- in keV)
      REAL C0(14)
C          (zero-order polynomial coefficients for 14 energy intervals)
      REAL C1(14)
C          (1st-order polynomial coefficients)
      REAL C2(14)
C          (2nd-order polynomial coefficients)
C
      data emax/.100,.284,.400,.532,.707,.867,1.303,
     #    1.840,2.471,3.210,4.038,7.111,8.331,10.0/ 
      data c0/17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,
     #    202.7,342.7,352.2,433.9,629.0,701.2/  
      data c1/608.1,267.9,18.8,66.8,145.8,-380.6,169.3, 
     #    146.8,104.7,18.7,18.7,-2.4,30.9,25.2/ 
      data c2/-2150.,-476.1,4.3,-51.4,-61.1,294.0,-47.7,
     #    -31.5,-17.0,0.0,0.0,0.75,0.0,0.0/ 
C


C     Start:
C
      E=energy/1.e3    
C          (convert to keV)
      do 100 i=1,14  
        if (E .lt. Emax(i)) goto 200
100   continue  
      i=14
200   sigism=(c0(i)+c1(i)*E+c2(i)*E*E)/E**3 * 1.E-24
      return
      end

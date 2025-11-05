      MODULE parameters
      IMPLICIT NONE 
      DOUBLE PRECISION, PARAMETER :: pi = 4.d0*datan(1.d0)
      DOUBLE PRECISION, PARAMETER :: pi2 = pi**2.d0
      DOUBLE PRECISION, PARAMETER :: alfem = 1.d0/137.d0
      DOUBLE PRECISION, PARAMETER :: sin2 = 0.23d0
      DOUBLE PRECISION, PARAMETER :: aw = DASIN(DSQRT(sin2))
      DOUBLE PRECISION, PARAMETER :: Mz = 91.2d0 !gauge boson mass -> Z0
      !DOUBLE PRECISION, PARAMETER :: M = 80.4d0 !gauge boson mass -> W
      DOUBLE PRECISION, PARAMETER :: rs = 13000.d0   !center of mass energy \sqrt{s}
      END MODULE 


      program pp_Z0_dimuon
      USE parameters
c     -----------------------------------------------------------------
c     The bins
c     -----------------------------------------------------------------
      ny = 60                      
      y_min =  2.0d0
      y_max =  4.5d0
      dy = (y_max-y_min)/ny
c     -----------------------------------------------------------------
      npt = 60             
      pt_min = 0.d0
      pt_max = 150.d0
      dpt = (pt_max-pt_min)/npt
c     -----------------------------------------------------------------      
      nm = 60 
      m_min = 60.d0
      m_max = 120.d0
      dm = (m_max-m_min)/nm


c     Defining the loops

      phiPmax = pi
      phiPmin = 0.d0
      phiMmax = pi
      phiMmin = 0.0

      kPmax = 200.d0
      kPmin = 20.d0
      kMmax = 200.d0
      KMmin = 20.d0

      yPmax = 4.5d0
      yPmin = 2.d0
      yMmax = 4.5d0 
      yMmin = 2.d0

      dphip = (phiPmax - phiPmin)/(nPoints)
      dphim = (phiMmax - phiMmin)/(nPoints)
      dkP = (kPmax - kPmin)/(nPoints)
      dkM = (kMmax - kMmin)/(nPoints)
      dyP = (yPmax - yPmin)/(nPoints)
      dyM = (yMmax - yMmin)/(nPoints)

      do iphiP = 0, nPoints
         phiP = phiPmin + (iphiP - 1)*dphip
      do iphiM = 0, nPoints
         phiM = phiMmin + (iphiM - 1)*dphim
      do ikP = 0, nPoints
         kP = kPmin + (ikP - 1)*dkP
      do ikM = 0, nPoints
         kM = kMmin + (ikM - 1)*dkM
      do iyP = 0, nPoints
         yP = yPmin + (iyP - 1)*dyP
      do iyM = 0, nPoints
         yM = yMmin + (iyM - 1)*dyM
            call CrossSectionDimuon(phiP, phiM, kP, kM, yP, yM, dsigma)
      end do
      end do
      end do
      end do
      end do
      end do

c     =================================================================
c     FILES WITH RESULTS
c     -----------------------------------------------------------------
      open(unit=41,file='Output/dsig_dy.dat',status='unknown')
      open(unit=42,file='Output/dsig_dpt.dat',status='unknown')
      open(unit=43,file='Output/dsig_dm.dat',status='unknown')
      open(unit=21,file='Output/dsig_dydpt_y2p25.dat',status='unknown')
      open(unit=22,file='Output/dsig_dydpt_y2p75.dat',status='unknown')
      open(unit=23,file='Output/dsig_dydpt_y3p25.dat',status='unknown')
      open(unit=24,file='Output/dsig_dydpt_y3p75.dat',status='unknown')
      open(unit=25,file='Output/dsig_dydpt_y4p25.dat',status='unknown')

c     =================================================================
c     building spectra
c     =================================================================
      sum_y = 0.d0
      do iy = 1, ny
      ! We are acessing the y bin middle point
            y = y_min + iy*dy - dy/2.d0
            write(41,101) y, sig_y(iy)
            
            sum_y = sum_y + sig_y(iy)*dy
ctest            write(*,*) y, sig_y(iy), dy
      enddo

      write(*,*) 
      write(*,*) 'sum_y: ', sum_y
c     -----------------------------------------------------------------

      sum_pt = 0.d0
      do ipt = 1, npt
      ! We are acessing the pT bin middle point
            pt = pt_min + ipt*dpt - dpt/2.d0

            write(42,101) pt, sig_pt(ipt)
            write(21,101) pt, sig_ypt1(ipt)
            write(22,101) pt, sig_ypt2(ipt)
            write(23,101) pt, sig_ypt3(ipt)
            write(24,101) pt, sig_ypt4(ipt)
            write(25,101) pt, sig_ypt5(ipt)

            sum_pt = sum_pt + sig_pt(ipt)*dpt
ctest            write(*,*) pt, sig_pt(ipt), dpt
      enddo

      write(*,*) 
      write(*,*) 'sum_pt: ', sum_pt
      
      close(42)
      close(41)


c     -----------------------------------------------------------------
      do iM = 1,nM 
      ! We are acessing the M bin middle point
            M = m_min + im*dm - dm/2.d0
            write(43,101) m, sig_m(im)
      enddo



      end program pp_Z0_dimuon

c     How to construct the weighting function for the cross section
c     calculation for the dimuon production in pp collisions at LHC?


c     In this subroutine we have to calculate the cross section for the 
c     dimuon production, however we must construct the bining process
c     therefore, we have to set the interval with if statements and
c     so on. Moreover, I'm still worried about how this gonna work and
c     if I'll be able to achive an aggriable result.

c     Furthermore, in the future the main ideia is to transpose this code 
c     to c++ in order to use the ROOT framework to handle the data
c     and plot the histograms.
      subroutine CrossSectionDimuon(phiP, phiM, kP, kM, yP, yM, dsigma)
      USE parameters

      


      end subroutine
      
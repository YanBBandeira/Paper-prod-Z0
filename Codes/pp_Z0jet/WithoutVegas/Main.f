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

      end program pp_Z0_dimuon

      subroutine CrossSectionDimuon(phiP, phiM, kP, kM, yP, yM, dsigma)
      USE parameters

      


      end subroutine
      
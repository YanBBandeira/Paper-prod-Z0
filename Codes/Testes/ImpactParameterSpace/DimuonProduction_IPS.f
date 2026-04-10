      MODULE parameters
      IMPLICIT NONE
       
       ! Physical constants
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846D0
      DOUBLE PRECISION, PARAMETER :: alpha_em = 1.0D0/137.035999084D0  
      DOUBLE PRECISION, PARAMETER :: SIN2 = 0.23D0 
      DOUBLE PRECISION, PARAMETER :: MZ = 91.1876D0  ! Z boson mass in GeV
      DOUBLE PRECISION, PARAMETER :: AW = DASIN(DSQRT(SIN2))  ! Weinberg angle
      DOUBLE PRECISION, PARAMETER :: RS = 13000.0D0 !GeV
      END MODULE 

      MODULE globals
      IMPLICIT NONE
      SAVE 
      DOUBLE PRECISION :: pt,Y, x1, m, x2
      END MODULE

      MODULE HankelIntegrals
      IMPLICIT NONE
      SAVE 
      DOUBLE PRECISION :: pt, eps, z
      END MODULE
c     ------------------------------------------------------------------
      PROGRAM DimuonProduction_IPS
      USE parameters
      USE globals
      IMPLICIT NONE
      INTEGER, PARAMETER :: nPoints = 20
      INTEGER :: ipt 

      DOUBLE PRECISION :: pt_min, pt_max, dpt
      DOUBLE PRECISION :: GeVtopb
      DOUBLE PRECISION :: IntegralResult, result
      DOUBLE PRECISION :: TotalIntegrand,DGAUSS
      DOUBLE PRECISION A,B,WK,EPS,RELERR
      INTEGER :: N,IWK,IMINPTS,IMAXPTS,IFFAIL,NFNEVL
      DIMENSION A(2), B(2), WK(1100000)
      EXTERNAL TotalIntegrand, DGAUSS

      pt_min = DLOG10(1.0D0)  ! Minimum pt in GeV
      pt_max = DLOG10(300.0D0) ! Maximum pt in GeV
      dpt = (pt_max - pt_min)/nPoints

      CALL InitPDFsetByName("CT10")

c     ---------------------Loop over pt-----------------------------

      open(10, file='pTdist_ips_results.dat')   
      GeVtopb = 0.398d9 ! Conversion factor from GeV^-2 to pb
      DO ipt = 1, nPoints
          pt = 10.0d0**(pt_min + (ipt-1)*dpt) ! - 0.9do
          
cOpt      call DADMUL(TotalIntegrand,N,A,B,
cOpt     * IMINPTS,IMAXPTS,EPS,WK,IWK,IntegralResult,
cOpt     * RELERR,NFNEVL,IFFAIL)

      IntegralResult = DGAUSS(TotalIntegrand, 2.0d0, 4.5d0, 1.0d-4)  
            result = IntegralResult * GeVtopb
          WRITE(*,*) 'pt = ', pt, ' Result = ', result
          WRITE(10,'(F12.6,1X,E12.5)') pt, result
      END DO

      END PROGRAM DimuonProduction_IPS      

c     ---------------------Integrand over Y----------------------
      DOUBLE PRECISION FUNCTION TotalIntegrand(yVar)
      USE parameters
      USE globals
      IMPLICIT NONE
      DOUBLE PRECISION :: yVar, DGAUSS
      EXTERNAL MIntegrand, DGAUSS 
      y = yVar
      TotalIntegrand = DGAUSS(MIntegrand,60.0d0,120.0d0,1.0d-4)
      RETURN 
      END 
c     ---------------------Integrand over M----------------------
      DOUBLE PRECISION FUNCTION MIntegrand(MVar)
      USE parameters
      USE globals
      IMPLICIT NONE

      DOUBLE PRECISION :: MVar, DileptonDecay
      DOUBLE PRECISION :: HadronicCrossSection,FuncPartonLevelSigma

      EXTERNAL DileptonDecay, FuncPartonLevelSigma
      HadronicCrossSection = DileptonDecay(M)
     &                 *FuncPartonLevelSigma(y,pt,MZ) 
      M = MVar
      MIntegrand = 4.0d0*pi*M*pt*HadronicCrossSection
ctest      write(*,*)'Integrand:',MIntegrand,pt,y,M
      RETURN 
      END  
c     ---------------------Integrand over M and Y----------------------
!       DOUBLE PRECISION FUNCTION TotalIntegrand(N, X)
!       USE parameters
!       USE globals
!       IMPLICIT NONE 
!       INTEGER :: N 
!       DOUBLE PRECISION :: X, y, jac
!       DOUBLE PRECISION :: HadronicCrossSection
!       DOUBLE PRECISION :: DileptonDecay, FuncPartonLevelSigma
!       DIMENSION :: X(2)
!       EXTERNAL DileptonDecay, FuncPartonLevelSigma

!       y = 2.d0 + (4.5d0 - 2.0d0)*X(1)  ! Rapidity from 2 to 4.5
!       M = 60.d0 + (20.d0 - 1.d0)*X(2)  ! Invariant mass from 1 to 20 GeV

!       jac = (4.5d0 - 2.0d0)*(120.d0 - 60.d0)

!       HadronicCrossSection = jac*DileptonDecay(M)
!      &                 *FuncPartonLevelSigma(y,pt,MZ) 

!       TotalIntegrand = 4.0d0*pi*M*pt*HadronicCrossSection
!       write(*,*)'Integrand: ',TotalIntegrand,y,M,jac
!       RETURN 
!       END


c     ---------------------Decay function-----------------------------

      FUNCTION DileptonDecay(MVar) 
      USE parameters
      IMPLICIT NONE
      DOUBLE PRECISION DileptonDecay,MVar, M, M2, MZ2 
      DOUBLE PRECISION DecayWidth, Branch, InvariantMassDist, Result

      M = MVar
      M2 = M*M
      Mz2 = Mz*Mz

      DecayWidth = ((alpha_em*M)/(6.d0*(dsin(2.d0*aw)**2.d0)))*(
     & (160.d0/3.d0)*(dsin(aw)**4.d0) - 40.d0*(dsin(aw)**2.d0) + 21.d0)

      Branch = 3.3d0/100.d0
     
      InvariantMassDist = (1.d0/pi)*
     & ((M*DecayWidth)/((M2 - Mz2)**2.d0 + (M*DecayWidth)**2.d0))

      Result = InvariantMassDist*Branch
      DileptonDecay = Result

ctest      write(*,*)'Decay: ',Result,DecayWidth,Branch,InvariantMassDist
      RETURN 
      END 

c     ---------------------Integrantion over z-------------------------
      FUNCTION FuncPartonLevelSigma(yVar,ptVar,mVar)
      USE globals
      USE parameters
      IMPLICIT NONE 

      DOUBLE PRECISION :: FuncPartonLevelSigma,yVar,ptVar,mVar
      DOUBLE PRECISION :: sqrt_M2pT2, preIntegral, varJacobian
      DOUBLE PRECISION :: IntegralResult, DGAUSS
      DOUBLE PRECISION :: PartonLevelSigma
      DOUBLE PRECISION :: IntegrandHadronicCrossSection

      EXTERNAL IntegrandHadronicCrossSection, DGAUSS

      sqrt_M2pT2 = DSQRT(mVar*mVar + ptVar*ptVar)
      x1 = (sqrt_M2pT2/rs)*DEXP(yVar)
      x2 = (sqrt_M2pT2/rs)*DEXP(-yVar)

      varJacobian = (2.0d0/rs)*DSQRT(mVar*mVar + ptVar*ptVar)
     &              *DCOSH(yVar)

      preIntegral = x1/(x1 + x2)

      IntegralResult = DGAUSS(IntegrandHadronicCrossSection,
     &                        x1,1.0d0,1.0d-4)

      FuncPartonLevelSigma = preIntegral*varJacobian*IntegralResult

ctest      WRITE(*,*) 'PartonLevelSigma: ',
ctest     * FuncPartonLevelSigma,x1,x2,yVar,ptVar,mVar

      RETURN 
      END

c     ---------------------z Integrand---------------------------------
      FUNCTION IntegrandHadronicCrossSection(zVar)
      USE parameters
      USE globals
      IMPLICIT NONE 

      DOUBLE PRECISION :: IntegrandHadronicCrossSection, zVar, xf
      DOUBLE PRECISION :: hs, Q, f( -6:6 )
      DOUBLE PRECISION :: u, d, s, c, b
      DOUBLE PRECISION :: uBar, dBar, sBar, cBar, bBar
      DOUBLE PRECISION :: mLight, mc, mb
      DOUBLE PRECISION :: gfgaup, gfgadw, gfgvup, gfgvdw
      DOUBLE PRECISION :: mf, gfv, gfa, z
      DOUBLE PRECISION :: upQuarkCS, downStrangeQuarksCS
      DOUBLE PRECISION :: charmQuarkCS, bottomQuarkCS
      DOUBLE PRECISION :: upQuarkFunc, dSQuarksFunc
      DOUBLE PRECISION :: charmQuarkFunc, bottomQuarkFunc
      DOUBLE PRECISION :: Result

      EXTERNAL InitPDFsetByName, evolvePDF


      z = zVar
      xf = x1/z
      hs = DSQRT(pt*pt + (1.0d0 - x1)*m*m)

      IF(hs.LE.1.3d0) THEN 
            Q = 0.0d0
      ELSE
            Q = hs
      END IF
      IF(xf.LE.1.d0)THEN
     
ctest      write(*,*) 'xf: ', xf, Q
     
      CALL evolvePDF(xf, Q, f)

      u = f(2)        !u
      d = f(1)        !d
      s = f(3)        !s
      c = f(4)        !c
      b = f(5)        !b
      uBar = f(-2)    !u_bar
      dBar = f(-1)    !d_bar
      sBar = f(-3)    !s_bar 
      cBar = f(-4)    !c_bar
      bBar = f(-5)    !b_bar

      mLight = 0.14d0 

      mc = 1.4d0 
      mb = 4.5d0

ctest      write(*,*) 'Masses: ', mLight, mc, mb
ctest      write(*,*) 'PDFs: ',u, d, s, c, b, uBar, dBar, sBar, cBar, bBar

      gfgaup = 0.5D0          ! FOR u,c,t quarks
      gfgadw = - 0.5d0          ! FOR d,s,b quarks
      gfgvup = 0.5D0 - (4.D0/3.D0)*SIN2   ! FOR u,c,t quarks
      gfgvdw = (2.D0/3.D0)*SIN2 - 0.5D0   ! FOR d,s,b quarks

c     quark up
      mf = mLight
      gfv = gfgvup
      gfa = gfgaup
      call PartonTargetCrossSection(upQuarkFunc,pt,z,M,mf,gfv,gfa)
      upQuarkCS = upQuarkFunc*(u + uBar)

c     quark down and quark strange because both are down type quarks 
c     with similar mass
      mf = mLight
      gfv = gfgvdw
      gfa = gfgadw
      call PartonTargetCrossSection(dSQuarksFunc,pt,z,M,mf,gfv,gfa)
      downStrangeQuarksCS = dSQuarksFunc*(d + s + dBar +sBar)

c     quark charm
      mf = MC
      gfv = gfgvup
      gfa = gfgaup
      call PartonTargetCrossSection(charmQuarkFunc,pt,z,M,mf,gfv,gfa)
      charmQuarkCS = charmQuarkFunc*(c + cBar)

c     quark bottom
      mf = MB
      gfv = gfgvdw
      gfa = gfgadw
      call PartonTargetCrossSection(bottomQuarkFunc,pt,z,M,mf,gfv,gfa)
      bottomQuarkCS = bottomQuarkFunc*(b + bBar)

ctest      write(*,*) 'Up quark cross section: ', upQuarkCS
ctest      write(*,*) 'd and s quarks cross section: ', downStrangeQuarksCS
ctest      write(*,*) 'Charm quark cross section: ', charmQuarkCS
ctest      write(*,*) 'Bottom quark cross section: ', bottomQuarkCS      
c     ------------------------------------------------------------------      
    

      Result = (upQuarkCS + downStrangeQuarksCS + 
     &         charmQuarkCS + bottomQuarkCS )/(z**2.d0)

      IntegrandHadronicCrossSection = Result

ctest      write(*,*)'Hadronic cross section:',IntegrandHadronicCrossSection
      ELSE
ctest      write(*,*) 'xf greater than 1, xf = ', xf
          IntegrandHadronicCrossSection = 0.0d0
      END IF

      RETURN 
      END


c     ---------------------Parton-Target Cross Section------------------
      SUBROUTINE PartonTargetCrossSection(Fvar,ptVar,zVar,mVar,
     *       mfVar,gfv,gfa)
      USE parameters
      USE HankelIntegrals
      IMPLICIT NONE 

      DOUBLE PRECISION :: Fvar, ptVar, zVar, mVar, mfVar,gfv,gfa
      DOUBLE PRECISION :: Eps2, pt2
      DOUBLE PRECISION :: z2, z4, M2, mf2, gfv2, gfa2
      DOUBLE PRECISION :: SigmaT, SigmaL, XiT, XiL
      DOUBLE PRECISION :: I1, I2, I3
      DOUBLE PRECISION :: D1Func, D2Func, Result
      DOUBLE PRECISION :: I1Integrand, I2Integrand, I3Integrand,DGAUSS

      EXTERNAL I1Integrand, I2Integrand, I3Integrand, DGAUSS

ctest      WRITE(*,*) 'Entry variables: ', ptVar, zVar, mVar, mfVar,gfv,gfa

      Eps2 = (1-zVar)*(mVar*mVar) + (zVar*zVar)*(mfVar*mfVar)
      
      Eps = DSQRT(Eps2)
      pt = ptVar
      z = zVar

ctest      WRITE(*,*) 'Eps2, Eps, pt, z: ', Eps2, Eps, pt, z

      I1 = DGAUSS(I1Integrand, 0.0d0, 1.0d0, 1.0d-8)
      I2 = DGAUSS(I2Integrand, 0.0d0, 1.0d0, 1.0d-8)
      I3 = DGAUSS(I1Integrand, 0.0d0, 1.0d0, 1.0d-8)

ctest      WRITE(*,*) 'I1, I2, I3: ', I1, I2, I3

      pt2 = pt*pt

      D1Func = (1.0d0/(pt2 + eps2))*I1 - (1.0d0/(4.0d0*eps))*I2 
      D2Func = (pt/(eps*(pt2 + eps2)))*I1 - (1.0d0/(2.0d0*eps2))*I1 
     &          + (1.0d0/(4.0d0*eps))*I2 

ctest      WRITE(*,*) 'D1Func, D2Func: ', D1Func, D2Func

      z2 = z*z
      z4 = z2*z2
      M2 = mVar*mVar
      mf2 = mfVar*mfVar
      gfv2 = gfv*gfv
      gfa2 = gfa*gfa


      SigmaT  = gfv2*z4*mf2 + gfa2*z2*mf2*((2.d0-z)**2.d0)
      SigmaL  =  gfv2*((1.d0-z)**2.d0)*M2 
     &        + gfa2*(((z2*mf2 + (1.d0-z)*M2)**2.d0)/M2)

      XiT = eps2*(1.d0 + ((1.d0 - z)**2.d0))*(gfv2 + gfa2)
      XiL = gfa2*eps2*((z2*mf2)/M2)

ctest      WRITE(*,*) 'SigmaT, SigmaL, XiT, XiL: ', SigmaT, SigmaL, XiT, XiL

      Result = (SigmaT + 2.0d0*SigmaL)*D1Func 
     &       + (XiT + 2.0d0*XiL)*D2Func
      Fvar = (alpha_em/(2.0d0*(pi*pi)*(DSIN(2.d0*aW)**2.d0)))*Result

ctest      WRITE(*,*) 'Parton-Target Cross Section: ',Fvar

      END  SUBROUTINE

c     --------------------------------Hankel Integrands-----------------
      FUNCTION I1Integrand(u)
      USE HankelIntegrals
      IMPLICIT NONE
      DOUBLE PRECISION :: u, I1Integrand
      DOUBLE PRECISION :: r, SigDip, Result, DipoleModel
      DOUBLE PRECISION :: bessj0, bessk0, bessj1, bessk1
      EXTERNAL bessj0, bessk0, bessj1, bessk1,DipoleModel
      
      r = u/(1.0d0 - u)
      SigDip = DipoleModel(z*r)
      Result = r*BESSEL_J0(pt*r)*bessk0(eps*r)*SigDip
      I1Integrand = Result/(1.0d0 - u)**2.d0
ctest      WRITE(*,*) 'J0, J1: ', BESSEL_J0(pt*r), BESSEL_J1(pt*r)
ctest      WRITE(*,*) 'K0, K1: ', BESSEL_K0(eps*r), BESSEL_K1(eps*r)
ctest      WRITE(*,*) 'I1 Integrand: ', u, r, SigDip, Result, I1Integrand

      RETURN
      END 

      FUNCTION I2Integrand(u)
      USE HankelIntegrals
      IMPLICIT NONE
      DOUBLE PRECISION :: u, I2Integrand
      DOUBLE PRECISION :: r, SigDip, Result, DipoleModel
      DOUBLE PRECISION :: bessj0, bessk0, bessj1, bessk1
      EXTERNAL bessj0, bessk0, bessj1, bessk1,DipoleModel
      
      r = u/(1.0d0 - u)
      SigDip = DipoleModel(z*r)
      Result = (r**2.0d0)*BESSEL_J0(pt*r)*bessk1(eps*r)*SigDip
      I2Integrand = Result/(1.0d0 - u)**2.d0
ctest      WRITE(*,*) 'J0, J1: ', bessj0(pt*r), bessj1(pt*r)
ctest      WRITE(*,*) 'K0, K1: ', bessk0(eps*r), bessk1(eps*r)
ctest      WRITE(*,*) 'I2 Integrand: ', u, r, SigDip, Result, I2Integrand

      RETURN
      END 

      FUNCTION I3Integrand(u)
      USE HankelIntegrals
      IMPLICIT NONE
      DOUBLE PRECISION :: u, I3Integrand
      DOUBLE PRECISION :: r, SigDip, Result, DipoleModel
      DOUBLE PRECISION :: bessj0, bessk0, bessj1, bessk1
      EXTERNAL bessj0, bessk0, bessj1, bessk1,DipoleModel
      
      r = u/(1.0d0 - u)
      SigDip = DipoleModel(z*r)
      Result = r*BESSEL_J1(pt*r)*bessk1(eps*r)*SigDip
      I3Integrand = Result/(1.0d0 - u)**2.d0

ctest      WRITE(*,*) 'J0, J1: ', bessj0(pt*r), bessj1(pt*r)
ctest      WRITE(*,*) 'K0, K1: ', bessk0(eps*r), bessk1(eps*r)
ctest      WRITE(*,*) 'I3 Integrand: ', u, r, SigDip, Result, I3Integrand

      RETURN
      END 

c     --------------------Dipole Model----------------------------------
      FUNCTION DipoleModel(rVar)
      USE globals
      USE parameters
      IMPLICIT NONE
      DOUBLE PRECISION :: DipoleModel,rVar
      DOUBLE PRECISION :: sigma0, lambdaDip, x0
      DOUBLE PRECISION :: xbj,Qs2, Result

      sigma0 = 23.03d0 / 0.389d0 ! in GeV^-2
      lambdaDip = 0.288d0
      x0 = 3.04d-4
      xbj = x2
      
ctest      WRITE(*,*) 'Dipole Model xbj: ', xbj

      Qs2 = (x0/xbj)**lambdaDip  ! in GeV^2



      Result = 1.0d0 - DEXP(- ((rVar*rVar)*Qs2/4.0d0))

      DipoleModel = sigma0 * Result

ctest      WRITE(*,*) 'Dipole Model: ', DipoleModel, rVar, Qs2

      RETURN 
      END

c     --------------------Bessel Functions------------------------------
      FUNCTION bessk0(x)
      DOUBLE PRECISION bessk0,x
      DOUBLE PRECISION bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     *0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     *-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      external bessi0
      if (x.le.2.d0) then
        y=x*x/4.d0
        bessk0=(-log(x/2.d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*
     *(p6+y*p7))))))
      else
        y=(2.d0/x)
        bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END


      FUNCTION bessk1(x)
      DOUBLE PRECISION bessk1,x
      DOUBLE PRECISION bessi1
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     *-0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     *0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      external bessi1
      if (x.le.2.d0) then
        y=x*x/4.d0
        bessk1=(log(x/2.d0)*bessi1(x))+(1.d0/x)*(p1+y*(p2+y*(p3+y*(p4+y*
     *(p5+y*(p6+y*p7))))))
      else
        y=2.d0/x
        bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*
     *q7))))))
      endif
      return
      END


      DOUBLE PRECISION FUNCTION bessi0(X)
      DOUBLE PRECISION X,Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     *    1.2067492D0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75d0)**2.d0
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75d0/AX
        BESSI0=(EXP(AX)/dsqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END



      DOUBLE PRECISION FUNCTION bessi1(X)
      DOUBLE PRECISION X,Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,
     *    0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2.D0
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI1=(EXP(AX)/dsqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
        IF(X.LT.0.D0)BESSI1=-BESSI1
      ENDIF
      RETURN
      END



      FUNCTION bessj0(x)
      DOUBLE PRECISION bessj0,x
      DOUBLE PRECISION y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     *s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END

      FUNCTION bessj1(x)
      DOUBLE PRECISION bessj1,x
      DOUBLE PRECISION y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
      endif
      return
      END


	DOUBLE PRECISION FUNCTION Ii1(PTQ,Z, ETTA) 
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI1/ptt, EPS, alfa
	EXTERNAL I1int, DGAUSS
	PARAMETER(kbes=600)
	DIMENSION zb(kbes)
	PTT  = PTQ 
	PTT2 = PTT*PTT
	ALFA = Z 
	EPS  = ETTA
	!WRITE(*,*) PTT, ALFA, EPS
C	file with zeros of bessel function J0
	OPEN(unit=32,status='old',file='bzero.dat')
		DO k=1,kbes
			READ(32,*) zb(k)
		ENDDO
	CLOSE(32)
	li    = 0.d0
	soma1 = 0.d0
C	loop for the integral between zeros
	DO 10 j =1,kbes
		ls         = zb(j)
		!WRITE(*,*)li,ls
		step_I1int = DGAUSS(I1int,li,ls,1.D-3)!/PTT2
		soma1      = soma1 + step_I1int 
		li 	     = ls
		!WRITE(*,*) soma1, step_I1int
10    CONTINUE 
	Ii1 = soma1 
	!WRITE(*,*)'I1', Ii1 
	!WRITE(*,*) I1
	RETURN 
	END

	DOUBLE PRECISION FUNCTION Ii2(PTQ,Z, ETTA)
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI2/ptt, EPS, alfa
	EXTERNAL I2int, DGAUSS 
	PARAMETER(kbes=600)
	DIMENSION zb(kbes)
	PTT  = PTQ 
	PTT2 = PTT*PTT
	PTT3 = PTT2*PTT
	ALFA = Z 
	EPS  = ETTA
	!WRITE(*,*) PTT, ALFA, EPS
C	file with zeros of bessel function J0
	OPEN(unit=33,status='old',file='bzero.dat')
		DO k=1,kbes
			READ(33,*) zb(k)
		ENDDO
	CLOSE(33)
	li     = 0.D0 
	soma2  = 0.D0
C	loop for the integral between zeros
	DO 10 j=1,kbes
		ls = zb(j)
		!WRITE(*,*)
		step_I2int = DGAUSS(I2int,li,ls,1.D-3)!/PTT3
		soma2 = soma2 + step_I2int 
		li = ls
		!WRITE(*,*) soma2, setp_I2int
10    CONTINUE
	Ii2 = soma2 
	!WRITE(*,*)'I2', Ii2 
	RETURN 
	END

	DOUBLE PRECISION FUNCTION Ii3(PTQ,Z, ETTA)
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI3/PTT,EPS,ALFA
	EXTERNAL I3int, DGAUSS
	PARAMETER(kbes=600)
	DIMENSION zb(kbes)
	PTT   = PTQ 
	PTT2  = PTT*PTT
	ALFA  = Z 
	EPS   = ETTA 
	!WRITE(*,*) PTT, ALFA, EPS 
c     files with zeros of bessel functions J1
	OPEN(unit=42,status='old',file='b1zeros.dat')
		DO k=1,kbes
			READ(42,*) zb(k)
		ENDDO
	CLOSE(42) 	
	li    = 0.D0
	SOMA3 = 0.D0
C	loop for the integral between zeros
	DO 10 j=1,kbes
		ls = zb(j)
		!WRITE(*,*)li,ls
		step_I3int = DGAUSS(I3int,li,ls,1.D-3)!/PTT
		soma3  = soma3 + step_I3int 
		li = ls 
		!WRITE(*,*) soma3, step_I3int 
10    CONTINUE 
	Ii3 = SOMA3 
	!WRITE(*,*) I3,ptt,alfa,eps
	RETURN
	END


	DOUBLE PRECISION FUNCTION I1int(uu) 
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI1/ptt, EPS, ALFA
	EXTERNAL bessj0,bessk0
	EXTERNAL bessj1,bessk1
	EXTERNAL sigdip
	u     = uu
	pt    = PTT 
	pt2   = pt*pt
	eps2  = eps*eps
	arg1  = u
	arg2  = (eps*u)/pt
	sdip  = sigdip(alfa*u/pt)
	I1int = u*bessj0(arg1)*bessk0(arg2)*sdip/pt2
	!WRITE(*,*)'I1int', bessj0(arg1),bessk1(arg2), sdip, u/pt, pt
	RETURN
	END 

	DOUBLE PRECISION FUNCTION I2int(uu)
        IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI2/ptt, EPS, ALFA
	EXTERNAL bessj0,bessk0
	EXTERNAL bessj1,bessk1
	EXTERNAL sigdip
	u     = uu
	u2    = u*u
	pt    = PTT 
	pt2   = pt*pt			
	pt3   = pt2*pt
	eps2  = eps*eps
	arg1  = u
	arg2  = (eps*u)/pt
	sdip  = sigdip(alfa*u/pt)
	I2int = u2*bessj0(arg1)*bessk1(arg2)*sdip/pt3
	!WRITE(*,*)'I2int', bessj0(arg1),bessk1(arg2), sdip, u/pt, pt
	RETURN 
	END

	DOUBLE PRECISION FUNCTION I3int(uu)
        IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI3/PTT,EPS,ALFA
	EXTERNAL bessj0,bessk0
	EXTERNAL bessj1,bessk1
	EXTERNAL sigdip
	u     = uu
	pt    = PTT 
	pt2   = pt*pt
	eps2  = eps*eps
	arg1  = u
	arg2  = EPS*u/pt
	sdip  = sigdip(alfa*u/pt)
	I3int = u*bessj1(arg1)*bessk1(arg2)*sdip/pt2
	!WRITE(*,*)'I3int', bessj1(arg1),bessk1(arg2), sdip, u/pt, I3int
	RETURN
	END

      MODULE parameters
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, PARAMETER :: pi = 4.d0*datan(1.d0)
      DOUBLE PRECISION, PARAMETER :: pi2 = pi**2.d0
      DOUBLE PRECISION, PARAMETER :: alfem = 1.d0/137.d0
      DOUBLE PRECISION, PARAMETER :: sin2 = 0.23d0
      DOUBLE PRECISION, PARAMETER :: aw = DASIN(DSQRT(sin2))
      DOUBLE PRECISION, PARAMETER :: Mz = 91.2d0 !gauge boson mass -> Z0
      !DOUBLE PRECISION, PARAMETER :: M = 80.4d0 !gauge boson mass -> W
      DOUBLE PRECISION, PARAMETER :: rs = 13000.d0   !center of mass energy \sqrt{s}
      END MODULE 
      
      MODULE globals
            IMPLICIT NONE
            SAVE
            DOUBLE PRECISION :: pt, x1, m
            DOUBLE PRECISION :: x2
      END MODULE

      program PL_grids_kslinear
      USE globals
      use parameters
      implicit none
      integer, parameter :: nPoints = 15
c     ------------------------------------------------------------------
      double precision y(nPoints), y_min, y_max, dy 
      double precision pt(nPoints), pt_min, pt_max, dpt
      double precision m(nPoints), m_min, m_max, dm
      double precision partonLevelSigma,FuncPartonLevelSigma 

      integer iy, ipt,im, iset
      external FuncPartonLevelSigma
      
      iset = 400001 !KS-2013-linear 
      call TMDinit(iset)
      call TMDset(iset)

      call InitPDFsetByName('CT10nlo')

      y_min = 1.d0
      y_max = 5.d0
      pt_min = dlog10(1.d0)
      pt_max = dlog10(500.d0)
      m_min = 50.d0
      m_max = 200.d0

      dy = (y_max - y_min)/(nPoints)
      dpt = (pt_max - pt_min)/(nPoints)
      dm = (m_max - m_min)/(nPoints)

c     ==================================================================
c     Output files
c     ==================================================================

      open(unit=11,file="DatFiles/tst_grid.dat",status='replace', 
     * action='write')

c     ==================================================================
c     Grid loop
c     ==================================================================
C     Usar as flags -O3 -fopenmp -march=native -ffast-math -funroll-loops
      !$OMP PARALLEL DO PRIVATE(iy, ipt,im, partonLevelSigma) SCHEDULE(DYNAMIC)
      do iy=1,nPoints
      y(iy) = y_min + (iy-1)*dy
      
            do ipt = 1,nPoints
      pt(ipt) = 10.d0**(pt_min + (ipt-1)*dpt) - 0.9d0
                  do im = 1,nPoints
      m(im) = m_min + (im-1)*dm
      partonLevelSigma = FuncPartonLevelSigma(y(iy),pt(ipt),m(im))


      write(*,*) '-------------------------------------'
      write(*,*) 'Computing grid point number ', iy, ipt, im,
     *  ' of ', nPoints, nPoints, nPoints
      write(*,*) 'Computing point: y = ', y(iy), ' pt = ', pt(ipt),
     *  ' m = ', m(im)  
      write(*,*)'Parton level cross section (pb/GeV): ',partonLevelSigma
      write(11,*) y(iy), pt(ipt), m(im), partonLevelSigma
                  end do
            end do 
      end do
      !$OMP END PARELLEL DO
100   format(2x,6(E10.4,2x))
      close(11)
      end program PL_grids_kslinear

c     ==================================================================
c     parton level cross section function
c     ==================================================================
      function FuncPartonLevelSigma(yVar,ptVar,mVar)
      use globals
      use parameters
      double precision FuncPartonLevelSigma, ptVar, yVar, mVar
      double precision y, pt, pt2, M2, x1, x2, M
      double precision result, units,IntegrandHadronicCrossSection
      double precision dgauss, sqrt_M2pT2
      external IntegrandHadronicCrossSection, dgauss 
      
      pt = ptVar
      y  = yVar

      M  = mVar
      M2 = M**2.d0
      pt2 = pt**2.d0

      sqrt_M2pT2 = DSQRT(M2 + pt2) 
      x1 = (sqrt_M2pT2/RS)*DEXP(y)
      x2 = (sqrt_M2pT2/RS)*DEXP(-y)

      result = dgauss(IntegrandHadronicCrossSection,x1,1.d0,1.d-4) 
              
    
      units = 0.389d9 !GeV-2 to pb
      FuncPartonLevelSigma = result*units
ctest      write(*,*) 'Variables: ', pt, y, x1, x2,M
      return 
      end 





c     =================================================================
c     =================================================================
c     =================================================================

      FUNCTION IntegrandHadronicCrossSection(alf)
      use globals

      USE parameters
      IMPLICIT NONE 
      DOUBLE PRECISION IntegrandHadronicCrossSection, alf, Result
      CHARACTER name*64
      DOUBLE PRECISION f(-6:6)
      DOUBLE PRECISION pt, x1, M, z, xf, hs, Q
      DOUBLE PRECISION pt2, M2 
      DOUBLE PRECISION u,d,s,c,b,uBar,dBar,sBar,cBar,bBar
      DOUBLE PRECISION MU,MD,MS,MC,MB
      DOUBLE PRECISION CFG, gfgaup, gfgadw, gfgvup, gfgvdw
      DOUBLE PRECISION upQuarkFunc, dSQuarksFunc
      DOUBLE PRECISION charmQuarkFunc, bottomQuarkFunc
      DOUBLE PRECISION upQuarkCS, downStrangeQuarksCS
      DOUBLE PRECISION charmQuarkCS, bottomQuarkCS
      DOUBLE PRECISION mf, gfv, gfa
      


      EXTERNAL InitPDFsetByName, evolvePDF

      
      z = alf 
      xf = x1/z 
      pt2 = pt*pt
      M2 = M*M

      hs = dsqrt(pt2 + (1.d0 - x1)*M2)

ctest      write(*,*) 'Hadronic variables: ', pt, x1, M, z, xf, hs

      
      if(hs.le.1.3d0) then
            Q = 0.d0 
            write(*,*) 'Q less than 1.3 GeV, Q = ', Q
      else
            Q = hs 
      end if
      if(xf.le.1.d0)then 
c     ------------------------------------------------------------------
c     This is the parton distribution function (PDF) initialization
c     ------------------------------------------------------------------   
      !name = 'CT10nlo' ! This is the PDF set name
      !call InitPDFsetByName(name)
      call evolvePDF(xf,q,f)

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
c     ------------------------------------------------------------------
c     light quarks
      MU    = 0.14D0   
      MD    = 0.14D0 
      MS    = 0.14D0
c     heavy quarks      
      MC    = 1.4D0 
      MB    = 4.5d0
c     ------------------------------------------------------------------
ctest      write(*,*) 'Masses', MU, MD, MS, MC, MB
ctest      write(*,*) 'Pdfs:', u,d,s,c,b,uBar,dBar,sBar,cBar,bBar

c     ------------------------------------------------------------------
c     Charge and coupling
c     ------------------------------------------------------------------
      CFG   = DSQRT(alfem)/dsin(2.d0*aW)
      gfgaup = 0.5D0          ! FOR u,c,t quarks
      gfgadw = - 0.5d0          ! FOR d,s,b quarks
      gfgvup = 0.5D0 - (4.D0/3.D0)*SIN2   ! FOR u,c,t quarks
      gfgvdw = (2.D0/3.D0)*SIN2 - 0.5D0   ! FOR d,s,b quarks
      
ctest      write(*,*) 'Charge: ', CFG
ctest      write(*,*) 'Axial couplings (up & down): ', gfgaup, gfgadw
ctest      write(*,*) 'Vector couplings (up & down): ', gfgvup, gfgvdw


c     ------------------------------------------------------------------
c     Evaluating each quark flavour
c     ------------------------------------------------------------------
c     quark up
      mf = MU
      gfv = gfgvup
      gfa = gfgaup
      call PartonTargetCrossSection(upQuarkFunc,pt,z,M,mf,gfv,gfa)
      upQuarkCS = upQuarkFunc*(u + uBar)

c     quark down and quark strange because both are down type quarks 
c     with similar mass
      mf = md
      gfv = gfgvdw
      gfa = gfgadw
      call PartonTargetCrossSection(dSQuarksFunc,pt,z,M,mf,gfv,gfa)
      upQuarkCS = dSQuarksFunc*(d + s + dBar +sBar)

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
      IntegrandHadronicCrossSection = 0.d0 
      END IF 
      RETURN
      END 

c     =================================================================
c     =================================================================
c     =================================================================

      SUBROUTINE PartonTargetCrossSection(Fvar,ptVar,zVar,mVar,
     *       mfVar,gfv,gfa)
      USE parameters
      IMPLICIT NONE
      DOUBLE PRECISION Fvar, ptVar, zVar, MVar, mfVar, gfv, gfa
      DOUBLE PRECISION gfv2,gfa2,mf,z,M2,pt
      DOUBLE PRECISION TransverseMomentumIntegralResult, preTerms
      DOUBLE PRECISION A,B,WK,EPS,RELERR
      INTEGER          N,IWK,IMINPTS,IMAXPTS,IFFAIL,NFNEVL
      DIMENSION A(2), B(2), WK(1100000)

      COMMON/PartonLeverCommons/gfv2,gfa2,mf,M2,pt,z
      EXTERNAL TransverseMomentumIntegral

      gfv2 = gfv*gfv
      gfa2 = gfa*gfa
      mf = mfVar
      pt = ptVar
      z  = zVar
      M2 = mVar*mVar

      preTerms = DSQRT(alfem)/(2.d0*pi2*DSIN(2.d0*aW))
ctest      write(*,*) 'Pre-terms: ', preTerms, gfv2, gfa2, mf, M2, pt, z
      !DADMUL routine parameters:
      N = 2               !Dimension
      IMINPTS = 500      !Min number of points
      IMAXPTS = 5000   !Max number of points
      EPS = 1.d-3          !Numerical precision    
      IWK = 110000        !Work array dimension
      A(1) = 0.d0       !Lower limit for variable 1 
      A(2) = 0.d0       !Lower limit for variable 2
      B(1) = 1.d0       !Upper limit for variable 1
      B(2) = 1.d0       !Upper limit for variable 2
    

      CALL DADMUL(TransverseMomentumIntegral,N,A,B,
     * IMINPTS,IMAXPTS,EPS,WK,IWK,TransverseMomentumIntegralResult,
     * RELERR,NFNEVL,IFFAIL)

      Fvar = preTerms*TransverseMomentumIntegralResult
ctest      write(*,*) 'Parton-target cross section: ', Fvar 
      END SUBROUTINE


c     =================================================================
c     =================================================================
c     =================================================================


      FUNCTION TransverseMomentumIntegral(N,X)
      IMPLICIT NONE 
      INTEGER N 
      DOUBLE PRECISION Int, akt, theta, kv, pi
      DOUBLE PRECISION TransverseMomentumIntegral, X, 
     &                 TransverseMomentumIntegrand 
      DIMENSION X(2)
      EXTERNAL TransverseMomentumIntegrand
      pi = 4.d0*datan(1.d0)

      !Integration variables 
      kv    = X(1)/(1.d0 - X(1)) !kv goes from zero to inf
      theta = X(2)*2.d0*pi       !theta goes from 0 to 2pi
ctest      write(*,*) 'int var ',x(2), x(1) 

      Int = TransverseMomentumIntegrand(kv,theta)
     &        *((1.d0 + kv)**2.d0)*(2.d0*pi) 
      ! both 2pi and (1+kv)^2 are jacobian factors
      
      TransverseMomentumIntegral = Int
      
ctest      write(*,*) 'Transverse Momentum Integral: ', Int, kv, theta
      RETURN
      END  

c     =================================================================
c     =================================================================


      FUNCTION TransverseMomentumIntegrand(akt,atheta)
      IMPLICIT NONE 
      
      DOUBLE PRECISION TransverseMomentumIntegrand, akt, atheta
      DOUBLE PRECISION Result, z2, z4, mf2, epps2, epps,
     &                 GammaL, GammaT, LambdaL, LambdaT, UGDF,
     &                 Epsilon1Var, Epsilon2Var
      DOUBLE PRECISION gfv2,gfa2,mf,M2,pt,z
      DOUBLE PRECISION ugd, Epsilon1, Epsilon2
      COMMON/PartonLeverCommons/gfv2,gfa2,mf,M2,pt,z
      EXTERNAL ugd, Epsilon1, Epsilon2

      z2 = z*z 
      z4 = z2*z2
      mf2 = mf*mf 

      GammaL  =  gfv2*((1.d0-z)**2.d0)*M2 
     &        + gfa2*(((z2*mf2 + (1.d0-z)*M2)**2.d0)/M2)
      GammaT  = gfv2*z4*mf2 + gfa2*z2*mf2*((2.d0-z)**2.d0)
      LambdaL = gfa2*((z2*mf2)/M2)
      LambdaT = (1.d0 + ((1.d0 - z)**2.d0))*(gfv2 + gfa2)

      epps2 = (1.d0 - z)*M2 + z2*mf2
      epps  = DSQRT(epps2)
      Epsilon1Var = Epsilon1(akt,atheta,pt,z,epps) 
      Epsilon2Var = Epsilon2(akt,atheta,pt,z,epps)

      UGDF = ugd(akt)

      Result = akt*UGDF*((GammaT + 2.d0*GammaL)*Epsilon1Var
     &        + (LambdaT + 2.d0*LambdaL)*Epsilon2Var)  
      TransverseMomentumIntegrand = Result

ctest      write(*,*)'Trasnverse integrand:', Epsilon1Var, 
ctest     &         Epsilon2Var, UGDF, Result
      RETURN
      END        

      
c     ==================================================================
c                       Epsilon's definitions
c     ==================================================================

      FUNCTION Epsilon1(akt,atheta,pt,z,epps)
      IMPLICIT NONE 
      DOUBLE PRECISION Epsilon1,akt,atheta,pt,z,epps
      DOUBLE PRECISION akt2,pt2,z2,epps2,tau2,eta
      DOUBLE PRECISION term1,term2,term3,Result

      akt2  = akt*akt
      pt2   = pt*pt
      z2    = z*z
      epps2 = epps*epps 
      
      tau2  = pt2 + z2*akt2 - 2.d0*z*pt*akt*dsin(atheta) 

      term1 = 1.d0/((tau2 + epps2)**2.d0 )
      term2 = 2.d0/((pt2+epps2)*(tau2 + epps2))
      term3 = 1.d0/((pt2 + epps2)**2.d0)
      Result = term1 - term2 + term3

      Epsilon1 = 0.5*Result 
      
ctest      write(*,*) 'Epsilon_1', Epsilon1, akt, atheta, pt, z, epps
      RETURN 
      END 

      FUNCTION Epsilon2(akt,atheta,pt,z,epps)
      IMPLICIT NONE 
      DOUBLE PRECISION Epsilon2,akt,atheta,pt,z,epps
      DOUBLE PRECISION akt2,pt2,z2,epps2,tau2,eta
      DOUBLE PRECISION term1,term2,term3,Result
      
      akt2  = akt*akt
      pt2   = pt*pt
      z2    = z*z
      epps2 = epps*epps 
      
      tau2  = pt2 + z2*akt2 - 2.d0*z*pt*akt*dsin(atheta) 
      eta   = pt2 - z*pt*akt*dsin(atheta)

      term1 = tau2/((tau2 + epps2)**2.d0 )
      term2 = (2.d0*eta)/((pt2+epps2)*(tau2 + epps2))
      term3 = pt2/((pt2 + epps2)**2.d0)
      Result = term1 - term2 + term3

      Epsilon2 = 0.5*Result 
      
ctest      write(*,*) 'Epsilon_2', Epsilon2, akt, atheta, pt, z, epps
      RETURN 
      END
c     ==================================================================
c                       Calling TMDlib
c     ==================================================================      


      FUNCTION ugd(akt) !WW
      use globals

      use parameters
      IMPLICIT NONE
      DOUBLE PRECISION ugd, akt, alphas, arg, akt2
      DOUBLE PRECISION pre_ugd, x2, F_KS, sc

      


      EXTERNAL F_KS,sc
      
c      alphas = 0.12d0
      arg =  ( 4.d0 /akt**2.d0)  + 1.9d0 
      alphas = sc(arg)
c      alphas = 0.2d0
      akt2 = akt**2.d0     

      pre_ugd = 4.d0*pi*alphas/3.d0
      ugd = (pre_ugd*F_KS(x2,akt))/akt2

ctest      write(*,*) 'UGD', ugd, pre_ugd, F_KS(x2,akt), akt2, alphas
      RETURN 
      END 

      FUNCTION F_KS(x2,akt)
      use parameters
      IMPLICIT NONE
      DOUBLE PRECISION F_KS,akt,x2 
      DOUBLE PRECISION x,xbar,kt,mu
      double precision up,ubar,dn,dbar,
     &                 strange,sbar,charm,cbar,bottom,bbar,glu

      INTEGER iset
      Integer kf
      Integer TMDnumberPDF

      external TMDnumberPDF
      
      if(akt.lt.399.d0.or.akt.gt.0.011d0)then
      if(x2.gt.1.d-2.or.x2.lt.1.d-8) then
      F_KS = 0.d0 

c      write(*,*) 'not ok', x2, akt

      else 
      x = x2
      kt = akt
      xbar = 0.d0
      mu = 100.d0
      call TMDpdf(kf,x,xbar,kt,mu,up,ubar,dn,dbar,
     & strange,sbar,charm,cbar,bottom,bbar,glu)
      F_KS = glu
c      write(*,*) 'ok', x, kt
ctest      write(*,*) iset, glu, akt    
c      write(*,*) 'ok'
      end if
      else
      F_KS = 0.d0 
      end if
      return
      END

      FUNCTION sc(Q2) !strong coupling
      use parameters
      DOUBLE PRECISION sc, Q2, b0,lambqcd,lambqcd2
      
      lambqcd = 0.35d0 !GeV 
      lambqcd2 = lambqcd**2.d0
      b0 = 11.d0 - (10.d0/3.d0)

      sc = (4.d0*pi)/(b0*dlog(Q2/lambqcd2))
ctest      write(*,*) 'sc', sc, Q2, lambqcd2,pi
      return 
      end 
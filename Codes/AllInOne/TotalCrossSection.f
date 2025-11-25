

      MODULE hadglobals
            IMPLICIT NONE
            SAVE
            DOUBLE PRECISION  pt, x1, m
            DOUBLE PRECISION  x2
      END MODULE

      MODULE parameters
            IMPLICIT NONE 
            SAVE
            DOUBLE PRECISION, PARAMETER :: pi = 4.d0*datan(1.d0)
            DOUBLE PRECISION, PARAMETER :: pi2 = pi**2.d0
            DOUBLE PRECISION, PARAMETER :: alfem = 1.d0/137.d0
            DOUBLE PRECISION, PARAMETER :: sin2 = 0.231
            DOUBLE PRECISION, PARAMETER :: aw = DASIN(DSQRT(sin2))
            DOUBLE PRECISION, PARAMETER :: Mz = 91.2d0 !gauge boson mass -> Z0
            !DOUBLE PRECISION, PARAMETER :: M = 80.4d0 !gauge boson mass -> W
            DOUBLE PRECISION, PARAMETER :: rs = 13000.d0   !center of mass energy \sqrt{s}
      END MODULE 

      MODULE globals
            IMPLICIT NONE
            SAVE
            INTEGER :: ny, npt, nm, iDoHist
            DOUBLE PRECISION :: y_min, y_max, dy
            DOUBLE PRECISION :: pt_min, pt_max, dpt
            DOUBLE PRECISION :: m_min, m_max, dm
            DOUBLE PRECISION, DIMENSION(1000) :: sig_y, sig_pt, sig_m
            DOUBLE PRECISION, DIMENSION(1000) :: sig_ypt1, sig_ypt2
            DOUBLE PRECISION, DIMENSION(1000) :: sig_ypt3, sig_ypt4
            DOUBLE PRECISION, DIMENSION(1000) :: sig_ypt5
      END MODULE globals

      MODULE vegasParams
            IMPLICIT NONE 
            SAVE 
            double precision avgi,sd,chi2a
            double precision alph 
            double precision xl(11),xu(11),acc,si,swgt,schi,xi(50,11) 
            integer ncall,itmx,nprn,ndev,it,ndo 
            integer ndmx, mds 
            integer ncall1, ncall2, itmx1, itmx2 
      END MODULE

      MODULE gridsParams
            IMPLICIT NONE 
            SAVE 
            INTEGER, parameter :: nPoints = 20
            DOUBLE PRECISION :: ent(nPoints + nPoints + nPoints)
            DOUBLE PRECISION :: ptGrid(nPoints)
            DOUBLE PRECISION :: yGrid(nPoints)
            DOUBLE PRECISION :: mGrid(nPoints)
            DOUBLE PRECISION :: PartonLevelGrid(nPoints,nPoints,nPoints)
      END MODULE

c---  Total cross section for pp -> (Z0 -> l l_bar) X 
      PROGRAM sigtot_pp_Z0   
      USE globals
      use gridsParams
      IMPLICIT NONE

      double precision vegasIntegrand 


c---     VEGAS definitions
C      double precision avgi,sd,chi2a
C      double precision alph 
C      double precision xl,xu,acc,si,swgt,schi,xi 
C      integer ncall,itmx,nprn,ndev,it,ndo 
C      integer ndmx, mds 
c      integer ncall1, ncall2, itmx1, itmx2 
c     -----------------------------------------------------------------
C      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
C      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
C      common/bveg3/alph,ndmx,mds
c     -----------------------------------------------------------------
      integer nprn, ncall, itmx, ndev, ndo, ndmx, mds, it, i
      double precision xl(11), xu(11), acc, alph, si, swgt, schi
      double precision xi(50,11), avgi,sd,chi2a
      common /bveg1/ ncall, itmx, nprn, ndev, xl, xu, acc
      common /bveg2/ it, ndo, si, swgt, schi, xi
      common /bveg3/ alph, ndmx, mds



c     =================================================================
c     common blocks (for Bining process)
c     -----------------------------------------------------------------
cold      INTEGER :: iy, ipt, im, ny, npt, nm, iDoHist
      INTEGER :: iy, ipt, im, iset
cold      DOUBLE PRECISION :: y_min, y_max, dy, pt_min, pt_max, y,pt,m
cold      DOUBLE PRECISION :: dpt, m_min, m_max, dm
      DOUBLE PRECISION :: sum_y, sum_pt, y , pt , m
cold      DOUBLE PRECISION, DIMENSION(1000) :: sig_y, sig_pt, sig_m
cold      DOUBLE PRECISION, DIMENSION(1000) :: sig_ypt1, sig_ypt2, sig_ypt3
cold      DOUBLE PRECISION, DIMENSION(1000) :: sig_ypt4, sig_ypt5

c     -----------------------------------------------------------------
cold      common/hist/iDoHist,sig_y,sig_pt,sig_m,
cold     &           sig_ypt1,sig_ypt2,sig_ypt3,sig_ypt4,sig_ypt5 
c     -----------------------------------------------------------------
cold	common/bin/ny,y_min,y_max,npt,pt_min,pt_max,dy,dpt,
cold     &           nm,m_min,m_max,dm
    
c     =================================================================

      EXTERNAL vegasIntegrand


      iset = 400001 !KS-2013-linear 
      call TMDinit(iset)
      call TMDset(iset)
      call InitPDFsetByName("CT10nlo")
c     =================================================================
c     BINNING PARAMETERS (stored in common block)
c     =================================================================
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

c     =================================================================
c     probing of the phase space
c     =================================================================
      nprn=0
      ncall=10000 !0
      itmx=10 
      iDoHist=0
      call VEGAS(3,vegasIntegrand,avgi,sd,chi2a)

c     =================================================================
c     INTEGRATION  ! VEGAS(n,sigma,avgi,sd,chi2a) n=7 - dimensions
c     =================================================================

      print*, avgi, "+-", sd            
      ncall=10000 !0             
      itmx=5 !10                                                  
      iDoHist=1
      call VEGAS1(3,vegasIntegrand,avgi,sd,chi2a)

      print*, avgi, "+-", sd
c     =================================================================
c     FILES WITH RESULTS
c     -----------------------------------------------------------------
      open(unit=41,
     & file='../pp_Z0jet/Output/dsig_dy.dat',status='unknown')
      open(unit=42,
     & file='../pp_Z0jet/Output/dsig_dpt.dat',status='unknown')
      open(unit=43,
     & file='../pp_Z0jet/Output/dsig_dm.dat',status='unknown')
      open(unit=21,
     & file='../pp_Z0jet/Output/dsig_dydpt_y2p25.dat',status='unknown')
      open(unit=22,
     & file='../pp_Z0jet/Output/dsig_dydpt_y2p75.dat',status='unknown')
      open(unit=23,
     & file='../pp_Z0jet/Output/dsig_dydpt_y3p25.dat',status='unknown')
      open(unit=24,
     & file='../pp_Z0jet/Output/dsig_dydpt_y3p75.dat',status='unknown')
      open(unit=25,
     & file='../pp_Z0jet/Output/dsig_dydpt_y4p25.dat',status='unknown')

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

101   FORMAT(1x,f8.4,1x,e12.4)
      END PROGRAM
       
 
c     =================================================================
c     =================================================================
c     =================================================================



      FUNCTION vegasIntegrand(x,vegasWgt)
      IMPLICIT NONE

      DOUBLE PRECISION vegasIntegrand, x(6), vegasWgt
      DOUBLE PRECISION ypVar, ymVar, phip, phim
      DOUBLE PRECISION ypVar_min, ypVar_max
      DOUBLE PRECISION ymVar_min, ymVar_max, jac
      DOUBLE PRECISION yp_min, yp_max, ym_min, ym_max
      DOUBLE PRECISION sigTot, physicalWgt, pi
      DOUBLE PRECISION ktp, ktm, ktp_min, ktp_max
      DOUBLE PRECISION ktm_min, ktm_max
      DOUBLE PRECISION kp_min, kp_max, km_min, km_max
      DOUBLE PRECISION pt2Var, pt2Var_min, pt2Var_max
      DOUBLE PRECISION m2Var, m2Var_min, m2Var_max
      DOUBLE PRECISION yVar, yVar_min, yVar_max
c     -----------------------------------------------------------------
      integer ncall,itmx,nprn,ndev,it,ndo
      integer ndmx,mds
      integer ncall1,ncall2,itmx1,itmx2,ihist
      double precision xl,xu,acc,si,swgt,schi,xi
      double precision alph
      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/bveg3/alph,ndmx,mds

      pi = 4.d0*datan(1.d0)

      ! ----- PHASE SPACE LIMITS -----
      yp_min = 2.0d0
      yp_max = 4.5d0
      ym_min = 2.0d0
      ym_max = 4.5d0
      kp_min = 20.d0
      kp_max = 120.d0
      km_min = 20.d0
      km_max = 120.d0

      ! ----- TRANSFORMATION x_i -> physical variables -----
      ypVar = yp_min + (yp_max - yp_min)*x(1)
      ymVar = ym_min + (ym_max - ym_min)*x(2)
      ktp   = kp_min + (kp_max - kp_min)*x(3)
      ktm   = km_min + (km_max - km_min)*x(4)
      phip  = 2.d0*pi*x(5)
      phim  = 2.d0*pi*x(6)

      ! ----- JACOBIAN -----
      ! Incluindo dφ e dkT para cada múon
      jac = (yp_max - yp_min)*(ym_max - ym_min) 
     &  * ( kp_max - kp_min ) * ( km_max - km_min )
     &  * ( 2.d0 * pi * ktp ) * ( 2.d0 * pi * ktm )

      ! TODO - Isso provavelmente está errado, estou misturando integração
      ! em pt com integração em kt. Verificar limites corretos.
      ! Além disso, ao integrar usando pt estou fazendo pp -> Z0 + X
      ! e no caso onde levo em consideração os kts estou fazendo
      ! pp -> mu + muBar + X. São processos diferentes.
      ! Logo, não sei como fazer as conexões e afins!
      
      physicalWgt = vegasWgt*jac/itmx
      CALL IntegrandSigma(sigTot,ypVar,ymVar,
     &      ktp,ktm,phip,phim,physicalWgt)

ctest      write(*,*) 'Integrand: ', vegasIntegrand, sigTot, jac
      RETURN 
      END 

c     =================================================================
c     =================================================================
c     =================================================================


      SUBROUTINE IntegrandSigma(sigTot,ypVar,ymVar,ktpVar,ktmVar,
     & phipVar,phimVar,physicalWgt) 
      USE globals
      USE parameters
      IMPLICIT NONE


      DOUBLE PRECISION InterpolateGrid
      DOUBLE PRECISION sigTot,ypVar,ymVar,phipVar,phimVar,
     &                 ktpVar,ktmVar,physicalWgt
      DOUBLE PRECISION varJacobian, preIntegral,Result,x1,x2
      DOUBLE PRECISION HadronicCrossSection, DileptonDecay
      DOUBLE PRECISION M2, DGAUSS, deltaY
      DOUBLE PRECISION ktp, ktm, ktp2, ktm2
      DOUBLE PRECISION ktpx, ktpy, ktmx, ktmy
      DOUBLE PRECISION ptx, pty, pt2
      DOUBLE PRECISION phip, phim, yp, ym
      DOUBLE PRECISION mp, mm
      DOUBLE PRECISION mperp_p, mperp_m, mperp_p2, mperp_m2
      DOUBLE PRECISION xp, xm, xf, mVar
      DOUBLE PRECISION pt, m, y, FuncPartonLevelSigma
      

      integer ipt,iy,im

      EXTERNAL DileptonDecay,DGAUSS, FuncPartonLevelSigma
             
      mp = 105.6d-3 !GeV
      mm = 105.6d-3 !GeV

      phip = phipVar
      phim = phimVar

      ktp = ktpVar
      ktm = ktmVar

      yp = ypVar
      ym = ymVar

      ktp2 = ktp**2.d0 
      ktm2 = ktm**2.d0
      ktpx = ktp*DCOS(phip)
      ktpy = ktp*DSIN(phip)
      ktmx = ktm*DCOS(phim)
      ktmy = ktm*DSIN(phim)

      mperp_p2 = ktpx**2.d0 + ktpy**2.d0 + mp**2.d0
      mperp_m2 = ktmx**2.d0 + ktmy**2.d0 + mm**2.d0
     
      mperp_p = dsqrt(mperp_p2)
      mperp_m = dsqrt(mperp_m2)

      m2 = mperp_p**2.d0 + mperp_m**2.d0 
     &       + 2.d0*mperp_p*mperp_m*DCOSH(yp - ym) - pt2
     
c     -----------------------------------------------------------------
      xp = (ktp/rs)*DEXP(yp)
      xm = (ktm/rs)*DEXP(ym)
      xf = xp + xm

ctest      write(*,*) 'phip, phim: ', phip, phim
ctest      write(*,*) 'yp, ym: ', yp, ym
ctest      write(*,*) 'xp, xm, xf: ', xp, xm, xf

ctest      write(*,*) 'ktp, ktm, phip, phim, pt: ', ktp,ktm,phip,phim,pt


c     ==================================================================
c     Boson variables
c     -----------------------------------------------------------------
      ptx = ktpx + ktmx
      pty = ktpy - ktmy
      pt  = DSQRT(ptx**2.d0 + pty**2.d0)
      pt2 = pt**2.d0
      pt = DSQRT(pt2)


      y  = DLOG(xf*(rs/DSQRT(pt2 + M2)))


      x1 = (DSQRT(M2 + pt**2.d0)/RS)*DEXP(y)
      x2 = (DSQRT(M2 + pt**2.d0)/RS)*DEXP(-y)



      

      
     
      mVar = DSQRT(m2) 

ctest      WRITE(*,*) 'Kinematics: ', y, pt, M, x1, x2
ctest      write(*,*) 'M2, M, x1, x2: ', m2,mVar, x1, x2

      varJacobian = (2.d0/rs)*DSQRT(M2 + pt2)*DCOSH(y)
copt      preIntegral = (x1/(x1 + x2))*varJacobian
      
ctest      write(*,*) 'Pre-integral: ', varJacobian, preIntegral


c     =================================================================
c     Delimiting boundaries
c     -----------------------------------------------------------------
      if(y.lt.2.d0.or.y.gt.4.5d0)then
            go to 101
      endif

      if(x1.gt.1.d0.or.x2.gt.1.d0) then
            SigTot = 0.d0
            go to 101
      endif
c     -----------------------------------------------------------------
      if(mVar.lt.60.d0.or.mVar.gt.120.d0) then
ctest            write(*,*) 'fora do range de M: ', mVar
            SigTot = 0.d0
            go to 101
      endif
      

      HadronicCrossSection =  FuncPartonLevelSigma(y,pt,mVar)  
copt      Result = preIntegral*DileptonDecay(Mvar)*HadronicCrossSection
      Result = varJacobian * 0.0336d0 * HadronicCrossSection   
copt      SigTot = xp*xm*ktp*ktm*Result
      SigTot = Result
      write(*,*)'Result: ',Result*physicalWgt*10.d0
ctest      write(*,*) 'Sigma total: ', SigTot, HadronicCrossSection
ctest      write(*,*) sigTot, Result

c     =================================================================
c     indices
c     -----------------------------------------------------------------
      iy  = idint((y-y_min)/dy) + 1
c     -----------------------------------------------------------------
      ipt = idint((pt-pt_min)/dpt) + 1
c     -----------------------------------------------------------------
      im  = idint((mVar-m_min)/dm) + 1
c     -----------------------------------------------------------------
ctest      write(*,*) iy, ipt, iDoHist, dy,dpt
c      write(*,*)
ctest      write(*,*) y, y_min, pt,pt_min
c     =================================================================
c     collecting spectra
c     =================================================================

      if(iDoHist.eq.1) then
c     -----------------------------------------------------------------     
      if(y.gt.y_min.and.y.lt.y_max) then 
            if(iy.gt.0.and.iy.le.ny) then
            sig_y(iy) = sig_y(iy) + SigTot*physicalWgt/dy
ctest            write(*,*) iy, ny
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
            sig_pt(ipt) = sig_pt(ipt) +  SigTot*physicalWgt/dpt
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.2.0d0.and.y.lt.2.5d0) then
                        deltaY = 2.5d0 - 2.0d0
      sig_ypt1(ipt) = sig_ypt1(ipt) +  SigTot*physicalWgt/(dpt*deltaY)
                  endif
            endif 
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.2.5d0.and.y.lt.3.0d0) then
                        deltaY = 3.0d0 - 2.5d0
      sig_ypt2(ipt) = sig_ypt2(ipt) + SigTot*physicalWgt/(dpt*deltaY)
                  endif
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.3.0d0.and.y.lt.3.5d0) then
                        deltaY = 3.5d0 - 3.0d0 
      sig_ypt3(ipt) = sig_ypt3(ipt) + SigTot*physicalWgt/(dpt*deltaY)
                  endif
            endif
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.3.5d0.and.y.lt.4.0d0) then
                        deltaY = 4.0d0 - 3.5d0 
      sig_ypt4(ipt) = sig_ypt4(ipt) + SigTot*physicalWgt/(dpt*deltaY)
                  endif
            endif
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.4.0d0.and.y.lt.4.5d0) then
                        deltaY = 4.5d0 - 4.0d0 
      sig_ypt5(ipt) = sig_ypt5(ipt) + SigTot*physicalWgt/(dpt*deltaY)
                  endif
            endif
      endif
      

c     -----------------------------------------------------------------
      if(mVar.gt.m_min.and.mVar.lt.m_max) then
            if(im.gt.0.and.im.le.nm) then
                  sig_m(im) = sig_m(im) + SigTot*physicalWgt/dm
            endif 
      endif
      

      endif
c     =================================================================
101   continue
      END   

c     =================================================================
c     =================================================================
c     =================================================================

      FUNCTION DileptonDecay(MVar) 
      USE parameters
      IMPLICIT NONE
      DOUBLE PRECISION DileptonDecay,MVar, M, M2, MZ2 
      DOUBLE PRECISION DecayWidth, Branch, InvariantMassDist, Result

      M = MVar
      M2 = M*M
      Mz2 = Mz*Mz

      DecayWidth = ((alfem*M)/(6.d0*(dsin(2.d0*aw)**2.d0)))*(
     & (160.d0/3.d0)*(dsin(aw)**4.d0) - 40.d0*(dsin(aw)**2.d0) + 21.d0)

      Branch = 3.3d0/100.d0
     
      InvariantMassDist = (1.d0/pi)*
     & ((M*DecayWidth)/((M2 - Mz2)**2.d0 + (M*DecayWidth)**2.d0))

      Result = InvariantMassDist*Branch
      DileptonDecay = Result

ctest      write(*,*)'Decay: ',Result,DecayWidth,Branch,InvariantMassDist
      RETURN 
      END 

  

      


c     ==================================================================
c     parton level cross section function
c     ==================================================================
      function FuncPartonLevelSigma(yVar,ptVar,mVar)
      use hadglobals
      use parameters
      double precision FuncPartonLevelSigma, ptVar, yVar, mVar
      double precision pt2, M2
      double precision result, units,IntegrandHadronicCrossSection
      double precision dgauss, sqrt_M2pT2
      external IntegrandHadronicCrossSection, dgauss
      external InitPDFsetByName, evolvePDF
      
      pt = ptVar
      y  = yVar

      M  = mVar
      M2 = M**2.d0
      pt2 = pt**2.d0

      sqrt_M2pT2 = DSQRT(M2 + pt2) 
      x1 = (sqrt_M2pT2/RS)*DEXP(yVar)
      x2 = (sqrt_M2pT2/RS)*DEXP(-yVar)

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
      use hadglobals

      USE parameters
      IMPLICIT NONE 
      DOUBLE PRECISION IntegrandHadronicCrossSection, alf, Result
      CHARACTER name*64
      DOUBLE PRECISION f(-6:6)
      DOUBLE PRECISION  z, xf, hs, Q
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
      
      call evolvePDF(xf,q**2.d0,f)
      
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
c      write(*,*) 'Masses', MU, MD, MS, MC, MB
c      write(*,*) 'Pdfs:', u,d,s,c,b,uBar,dBar,sBar,cBar,bBar

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
     &         charmQuarkCS + bottomQuarkCS )/z!(z**2.d0)

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
      use hadglobals

      use parameters
      IMPLICIT NONE
      DOUBLE PRECISION ugd, akt, alphas, arg, akt2
      DOUBLE PRECISION pre_ugd, F_KS, sc

      


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








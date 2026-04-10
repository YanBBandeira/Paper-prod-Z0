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
            INTEGER, parameter :: nPoints = 600
            DOUBLE PRECISION :: ent(nPoints + nPoints)
            DOUBLE PRECISION :: ptGrid(nPoints)
            DOUBLE PRECISION :: yGrid(nPoints)
            DOUBLE PRECISION :: mGrid(nPoints)
            DOUBLE PRECISION :: PartonLevelGrid(nPoints,nPoints)
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
      INTEGER :: iy, ipt, im
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

      character*2000 :: File
      

      File = "Grids/DatFiles/ksnonlinear_grid.dat"

      call read_grid(File)

      ! Inicialização dos bins
      ny = 100
      y_min = 2.0d0
      y_max = 4.5d0
      dy = (y_max - y_min) / ny

      npt = 100
      pt_min = 0.d0
      pt_max = 200.d0
      dpt = (pt_max - pt_min) / npt

      nm = 100
      m_min = 60.d0
      m_max = 120.d0
      dm = (m_max - m_min) / nm

      ! Inicializa arrays
      !sig_y = 0.d0
      !sig_pt = 0.d0
      !sig_m = 0.d0
      !sig_ypt1 = 0.d0
      !sig_ypt2 = 0.d0
      !sig_ypt3 = 0.d0
      !sig_ypt4 = 0.d0
      !sig_ypt5 = 0.d0

c     =================================================================
c     probing of the phase space
c     =================================================================
      nprn=0
      ncall=10000000 
      itmx= 10 
      iDoHist=0
      do i=1,6
      xl(i) = 0.d0
      xu(i) = 1.d0
      end do
      call VEGAS(6,vegasIntegrand,avgi,sd,chi2a)
      print*, avgi, "+-", sd  

c     =================================================================
c     INTEGRATION  ! VEGAS(n,sigma,avgi,sd,chi2a) n=7 - dimensions
c     =================================================================
      ncall= 900000 !0
      itmx=15 
      iDoHist=0
      call VEGAS(6,vegasIntegrand,avgi,sd,chi2a)
      print*, avgi, "+-", sd          

      ncall= 10000000 !0             
      itmx=30!10                                                  
      iDoHist=1
      call VEGAS(6,vegasIntegrand,avgi,sd,chi2a)

      print*, avgi, "+-", sd
c     =================================================================
c     FILES WITH RESULTS
c     -----------------------------------------------------------------
      open(unit=41,file='Output/ksnonlinear_dsig_dy.dat'
     * ,status='unknown')
      open(unit=42,file='Output/ksnonlinear_dsig_dpt.dat'
     * ,status='unknown')
      open(unit=43,file='Output/ksnonlinear_dsig_dm.dat'
     * ,status='unknown')
      open(unit=21,file='Output/ksnonlinear_dsig_dydpt_y2p25.dat'
     * ,status='unknown')
      open(unit=22,file='Output/ksnonlinear_dsig_dydpt_y2p75.dat'
     * ,status='unknown')
      open(unit=23,file='Output/ksnonlinear_dsig_dydpt_y3p25.dat'
     * ,status='unknown')
      open(unit=24,file='Output/ksnonlinear_dsig_dydpt_y3p75.dat'
     * ,status='unknown')
      open(unit=25,file='Output/ksnonlinear_dsig_dydpt_y4p25.dat'
     * ,status='unknown')

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

      !write(*,*) 
      !write(*,*) 'sum_y: ', sum_y
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

      !write(*,*) 
      !write(*,*) 'sum_pt: ', sum_pt


      
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
      DOUBLE PRECISION ktp, ktm
      DOUBLE PRECISION kp_min, kp_max, km_min, km_max
      integer nprn, ncall, itmx, ndev, ndo, ndmx, mds, it
      double precision xl(11), xu(11), acc, alph, si, swgt, schi
      double precision xi(50,11)
      common /bveg1/ ncall, itmx, nprn, ndev, xl, xu, acc
      common /bveg2/ it, ndo, si, swgt, schi, xi
      common /bveg3/ alph, ndmx, mds
c     -----------------------------------------------------------------


c     =================================================================
c     PHASE SPACE
c     =================================================================
ctest      write(*,*) x(1), x(2), x(3), x(4), x(5), x(6)
      pi  = 4.d0*datan(1.d0)  
      kp_max = 500000.d0
      kp_min = 20.d0
      km_max = 500.d0 
      km_min = 20.d0
      ypVar_max = 4.5d0
      ypVar_min = 2.0d0
      ymVar_max = 4.5d0 
      ymVar_min = 2.d0

ctest      write(*,*) pi, kp_max, kp_min, x(1), x(2), x(3), x(4), x(5), x(6)
c     -----------------------------------------------------------------
      ypVar = ypVar_min + (ypVar_max - ypVar_min)*x(1)
      ymVar = ymVar_min + (ymVar_max - ymVar_min)*x(2)
      ktp   = X(3)/(1.d0 - X(3))!kp_min + (kp_max - kp_min)*x(3)
      ktm   = X(4)/(1.d0 - X(4))!km_min + (km_max - km_min)*x(4)
      phip  = 2.d0*pi*x(5)
      phim  = 2.d0*pi*x(6)
        
     
      
ctest      write(*,*) ypVar, ymVar, ktp, ktm, phip, phim
c     =================================================================
c     jacobian: x(n) ----> phase space
c     =================================================================
      jac = (ypVar_max - ypVar_min)*(ymVar_max - ymVar_min) 
c     &  * ( kp_max - kp_min ) * ( km_max - km_min )
     &  * ( 2.d0 * pi * ktp ) * ( 2.d0 * pi * ktm )
     &  / ( (1.d0 - X(3))**2.d0 * (1.d0 - X(4))**2.d0 )
      ! TODO - Isso provavelmente está errado, estou misturando integração
      ! em pt com integração em kt. Verificar limites corretos.
      ! Além disso, ao integrar usando pt estou fazendo pp -> Z0 + X
      ! e no caso onde levo em consideração os kts estou fazendo
      ! pp -> mu + muBar + X. São processos diferentes.
      ! Logo, não sei como fazer as conexões e afins!
      
      physicalWgt = vegasWgt*jac/itmx
      CALL IntegrandSigma(sigTot,ypVar,ymVar,
     &      ktp,ktm,phip,phim,physicalWgt)
ctest      write(*,*) vegasWgt, jac, itmx
ctest      write(*,*) 'Integrand: ', vegasIntegrand, sigTot, jac
      vegasIntegrand = sigTot*jac
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
      DOUBLE PRECISION pt, m, y

c     =================================================================
c     common blocks (for Bining process)
c     -----------------------------------------------------------------
cold      double precision y,y_min,y_max,dy,sum_y
cold      double precision pt,pt_min,pt_max,dpt,sum_pt
cold      double precision m,m_min,m_max,dm
cold      double precision sig_y(1000),sig_pt(1000), sig_m(1000)
cold      double precision sig_ypt1(1000),sig_ypt2(1000),sig_ypt3(1000)
cold      double precision sig_ypt4(1000),sig_ypt5(1000)
c     -----------------------------------------------------------------
      integer ipt,iy,im
cold      integer ny,npt,nm
cold      integer iDoHist 
c     -----------------------------------------------------------------
cold      common/hist/iDoHist,sig_y,sig_pt,sig_m,
cold     &           sig_ypt1,sig_ypt2,sig_ypt3,sig_ypt4,sig_ypt5 
c     -----------------------------------------------------------------
cold	common/bin/ny,y_min,y_max,npt,pt_min,pt_max,dy,dpt,
cold     &           nm,m_min,m_max,dm
    
c     =================================================================


      EXTERNAL DileptonDecay,DGAUSS, InterpolateGrid
             
c---- Inicializando o result
      sigTot = 0.d0

      mp = 105.6d-3 !GeV
      mm = 105.6d-3 !GeV

      phip = phipVar
      phim = phimVar

      ktp = ktpVar
      ktm = ktmVar

      if(ktpVar.lt.20.d0.or.ktmVar.lt.20.d0) then
            sigTot = 0.d0
            go to 101
      endif

      yp = ypVar
      ym = ymVar
ctest      write(*,*) sigTot,ypVar,ymVar,ktpVar,ktmVar,
ctest     & phipVar,phimVar,physicalWgt
c--- Delimiting boundaries
      if(yp.lt.2.d0.or.yp.gt.4.5d0)then
            sigTot = 0.d0
            go to 101
      endif

      if(ym.lt.2.d0.or.ym.gt.4.5d0)then
            sigTot = 0.d0
            go to 101
      endif


      ktp2 = ktp**2.d0 
      ktm2 = ktm**2.d0
      ktpx = ktp*DCOS(phip)
      ktpy = ktp*DSIN(phip)
      ktmx = ktm*DCOS(phim)
      ktmy = ktm*DSIN(phim)

      ptx = ktpx + ktmx
      pty = ktpy + ktmy
      pt  = DSQRT(ptx**2.d0 + pty**2.d0)
      pt2 = pt**2.d0

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
      



      y  = DLOG(xf*(rs/DSQRT(pt2 + M2)))
      
      !y = (yp + ym)/2.d0
      x1 = (DSQRT(M2 + pt**2.d0)/RS)*DEXP(y)
      x2 = (DSQRT(M2 + pt**2.d0)/RS)*DEXP(-y)
      !x1 = (DSQRT(ktp2)/rs)*DEXP(yp) + (DSQRT(ktm2)/rs)*DEXP(ym) 
      !x2 = (DSQRT(ktp2)/rs)*DEXP(-yp) + (DSQRT(ktm2)/rs)*DEXP(-ym)
ctest      write(*,*) 'rapidity value: ', y, x1, x2  
c--- Delimiting boundaries
      if(y.lt.2.d0.or.y.gt.4.5d0)then
            sigTot = 0.d0
            go to 101
      endif

      if(x1.gt.1.d0.or.x2.gt.1.d0) then
            SigTot = 0.d0
            go to 101
      endif
      mVar = DSQRT(m2) 
      if(mVar.lt.60.d0.or.mVar.gt.120.d0) then
ctest            write(*,*) 'fora do range de M: ', mVar
            SigTot = 0.d0
            go to 101
      endif
      if(pt.gt.300.d0) then
            SigTot = 0.d0
            go to 101
      endif 
      

      
     


ctest      WRITE(*,*) 'Kinematics: ', y, pt, M, x1, x2
ctest      write(*,*) 'M2, M, x1, x2: ', m2,mVar, x1, x2

      varJacobian = (2.d0/rs)*DSQRT(M2 + pt2)*DCOSH(y)
      preIntegral = (x1/(x1 + x2))*varJacobian
      
ctest      write(*,*) 'Pre-integral: ', varJacobian, preIntegral
    

      HadronicCrossSection =  InterpolateGrid(y,pt,mVar)  
      Result = HadronicCrossSection*DileptonDecay(Mvar)
ctest      write(*,*) DileptonDecay(Mvar), HadronicCrossSection   
ctest      print *, ktp, ktm
      SigTot = Result*0.389d9!*1.0d9!*0.389d9*pi ! to pb
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

      FUNCTION DileptonDecay(MVarr) 
      USE parameters
      IMPLICIT NONE
      DOUBLE PRECISION DileptonDecay,MVarr, M, M2, MZ2 
      DOUBLE PRECISION DecayWidth, Branch, InvariantMassDist, Result

      M = MVarr
      M2 = M*M
      Mz2 = Mz*Mz

      DecayWidth = ((alfem*M)/(6.d0*(dsin(2.d0*aw)**2.d0)))*(
     & (160.d0/3.d0)*(dsin(aw)**4.d0) - 40.d0*(dsin(aw)**2.d0) + 21.d0)

      Branch = 3.366d0/100.d0
     
      InvariantMassDist = (1.d0/pi)*
     & ((M*DecayWidth)/((M2 - Mz2)**2.d0 + (M*DecayWidth)**2.d0))

      Result = InvariantMassDist*Branch
      DileptonDecay = Result

ctest      write(*,*)'Decay: ',Result,DecayWidth,Branch,InvariantMassDist
      RETURN 
      END 



      function InterpolateGrid(yVar,ptVar,mVar)
      use gridsParams
      implicit none

      integer, parameter :: narg=2
      
      integer :: nent(2)


      double precision arg(narg)
      double precision y, pt, yVar, ptVar,InterpolateGrid
      double precision m, mVar
      double precision DFINT
     
      
      y = yVar
      pt = ptVar
      m = mVar


      nent(1) = nPoints
      nent(2) = nPoints
copt      nent(3) = nPoints

      arg(1) = y
      arg(2) = pt
copt      arg(3) = m

      InterpolateGrid = DFINT(narg,arg,nent,ent,PartonLevelGrid)
ctest      write(*,*) 'Interpolated value: ', InterpolateGrid, y, pt, m
      return
      end 


      subroutine read_grid(OutputPath)
      use gridsParams
      implicit none
      integer, parameter :: iFile = 11
      integer :: i, j, k
      double precision yVar, ptVar, mVar, PartonLevelVar

      character*2000 :: OutputPath
      

      open(iFile,file=trim(OutputPath),status='old')
c      open(iFile,
c     *  file="Grids/DatFiles/kslinear_grid.dat",status="old")
      
      do i = 1, nPoints
            do j = 1, nPoints
copt                  do k = 1, nPoints
      read(iFile,*) yVar, ptVar, PartonLevelVar
      yGrid(i) = yVar
      ptGrid(j) = ptVar
      PartonLevelGrid(i,j) = PartonLevelVar
ctest            write(*,*) 'Reading grid point:', i, j, k
ctest            write(*,*) 'Grid', yVar, ptVar, mVar, PartonLevelVar
copt                  end do
            end do
      end do
      close(iFile)
      do i = 1, nPoints
         ent(i) = yGrid(i)
      end do 

      do j = 1, nPoints
         ent(nPoints + j) = ptGrid(j)
      end do

copt      do k = 1, nPoints
copt         ent(nPoints + nPoints + k) = mGrid(k)
copt      end do
      
      end subroutine








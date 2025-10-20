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
c     =================================================================
c            total cross section for pp -> (Z0 -> l l_bar) X 
c     =================================================================
      PROGRAM sigtot_pp_Z0   
      IMPLICIT NONE

      double precision vegasIntegrand 
      integer iset


c     =================================================================
c     VEGAS definitions
c     -----------------------------------------------------------------
      double precision avgi,sd,chi2a
      double precision x(3),wgt 
      double precision alph, pi 
      double precision xl,xu,acc,si,swgt,schi,xi 
c     -----------------------------------------------------------------
      integer ncall,itmx,nprn,ndev,it,ndo 
      integer ndmx, mds 
      integer ncall1, ncall2, itmx1, itmx2,ihist 
c     -----------------------------------------------------------------
      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/bveg3/alph,ndmx,mds
c     =================================================================




c     =================================================================
c     common blocks (for Bining process)
c     -----------------------------------------------------------------
      double precision y,y_min,y_max,dy,sum_y
      double precision pt,pt_min,pt_max,dpt,sum_pt
      double precision m,m_min,m_max,dm
      double precision sig_y(1000),sig_pt(1000), sig_m(1000)
      double precision sig_ypt1(1000),sig_ypt2(1000),sig_ypt3(1000)
      double precision sig_ypt4(1000),sig_ypt5(1000)
c     -----------------------------------------------------------------
      integer ipt,iy,im
      integer ny,npt,nm
      integer iDoHist 
c     -----------------------------------------------------------------
      common/hist/iDoHist,sig_y,sig_pt,sig_m,
     &           sig_ypt1,sig_ypt2,sig_ypt3,sig_ypt4,sig_ypt5 
c     -----------------------------------------------------------------
	common/bin/ny,y_min,y_max,npt,pt_min,pt_max,dy,dpt,
     &           nm,m_min,m_max,dm
    
c     =================================================================



      EXTERNAL vegasIntegrand

copt      iset = 400001 !KS-2013-linear 
copt      call TMDinit(iset)
copt      call TMDset(iset)


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

101   FORMAT(1x,f8.4,1x,e12.4)
      END PROGRAM
       
 
c     =================================================================
c     =================================================================
c     =================================================================



      FUNCTION vegasIntegrand(x,vegasWgt)
      IMPLICIT NONE

      DOUBLE PRECISION vegasIntegrand, x(3), vegasWgt
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

c     =================================================================
c     PHASE SPACE
c     =================================================================
      pi  = 4.d0*datan(1.d0)  
!       kp_max = 200.d0
!       kp_min = 20.d0
!       km_max = 200.d0 
!       km_min = 20.d0
!       ypVar_max = 4.5d0
!       ypVar_min = 2.0d0
!       ymVar_max = 4.5d0 
!       ymVar_min = 2.d0

!       kp_max = 120.d0
!       kp_min = 20.d0
!       km_max = 120.d0 
!       km_min = 20.d0
      
! c     -----------------------------------------------------------------
!       ypVar   = yp_min + (yp_max - yp_min)*x(1)
!       ymVar   = ym_min + (ym_max - ym_min)*x(2)
!       ktp   = kp_min + (kp_max - kp_min)*x(3)
!       ktm   = km_min + (km_max - km_min)*x(4)
!       phip = 2.d0*pi*x(5)
!       phim = 2.d0*pi*x(6)
! c     =================================================================
! c     jacobian: x(n) ----> phase space
! c     =================================================================
!       jac = (yp_max - yp_min)*(ym_max - ym_min)*(kp_max - kp_min)
!      &      *(km_max - km_min)*((2.d0*pi)**2.d0)

      ! TODO - Isso provavelmente está errado, estou misturando integração
      ! em pt com integração em kt. Verificar limites corretos.
      ! Além disso, ao integrar usando pt estou fazendo pp -> Z0 + X
      ! e no caso onde levo em consideração os kts estou fazendo
      ! pp -> mu + muBar + X. São processos diferentes.
      ! Logo, não sei como fazer as conexões e afins!
      
c     -----------------------------------------------------------------
      pt2Var_min = 0.d0**2.d0
      pt2Var_max = 200.d0**2.d0
      m2Var_min = 60.d0**2.d0
      m2Var_max = 120.d0**2.d0
      yVar_min = 2.0d0
      yVar_max = 4.5d0

      yVar  = yVar_min + (yVar_max - yVar_min)*x(1)
      pt2Var = pt2Var_min + (pt2Var_max - pt2Var_min)*x(2)
      m2Var = m2Var_min + (m2Var_max - m2Var_min)*x(3)
      
c     =================================================================
c     jacobian: x(n) ----> phase space
c     =================================================================
      jac = (yVar_max - yVar_min)*(pt2Var_max - pt2Var_min)
     &      *(m2Var_max - m2Var_min)

      physicalWgt = vegasWgt*jac/itmx
copt      CALL IntegrandSigma(sigTot,ypVar,ymVar,
copt     &      ktp,ktm,phip,phim,physicalWgt)
      CALL IntegrandSigma(sigTot,yVar,pt2Var,m2Var,physicalWgt)
      vegasIntegrand = jac*sigTot
ctest      write(*,*) 'Integrand: ', vegasIntegrand, sigTot, jac
      RETURN 
      END 

c     =================================================================
c     =================================================================
c     =================================================================


copt      SUBROUTINE IntegrandSigma(sigTot,ypVar,ymVar,ktpVar,ktmVar,
copt     & phipVar,phimVar,physicalWgt) 
      SUBROUTINE IntegrandSigma(sigTot,yVar,pt2Var,m2Var,physicalWgt)
      USE parameters
      IMPLICIT NONE

      double precision yVar, pt2Var, m2Var
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
      

c     =================================================================
c     common blocks (for Bining process)
c     -----------------------------------------------------------------
      double precision y,y_min,y_max,dy,sum_y
      double precision pt,pt_min,pt_max,dpt,sum_pt
      double precision m,m_min,m_max,dm
      double precision sig_y(1000),sig_pt(1000), sig_m(1000)
      double precision sig_ypt1(1000),sig_ypt2(1000),sig_ypt3(1000)
      double precision sig_ypt4(1000),sig_ypt5(1000)
c     -----------------------------------------------------------------
      integer ipt,iy,im
      integer ny,npt,nm
      integer iDoHist 
c     -----------------------------------------------------------------
      common/hist/iDoHist,sig_y,sig_pt,sig_m,
     &           sig_ypt1,sig_ypt2,sig_ypt3,sig_ypt4,sig_ypt5 
c     -----------------------------------------------------------------
	common/bin/ny,y_min,y_max,npt,pt_min,pt_max,dy,dpt,
     &           nm,m_min,m_max,dm
    
c     =================================================================

      
      COMMON/xbj/x2
      COMMON/hadronicVariables/pt, x1, M

      EXTERNAL DileptonDecay,DGAUSS, InterpolateGrid
             
copt      mp = 105.6d-3 !GeV
copt      mm = 105.6d-3 !GeV

copt      phip = phipVar
copt      phim = phimVar

copt      ktp = ktpVar
copt      ktm = ktmVar

copt      yp = ypVar
ctop      ym = ymVar

copt      ktp2 = ktp**2.d0 
copt      ktm2 = ktm**2.d0DileptonDecay
copt      ktpx = ktp*DCOS(phip)
copt      ktpy = ktp*DSIN(phip)
copt      ktmx = ktm*DCOS(phim)
copt      ktmy = ktm*DSIN(phim)

copt      mperp_p2 = ktpx**2.d0 + ktpy**2.d0 + mp**2.d0
copt      mperp_m2 = ktmx**2.d0 + ktmy**2.d0 + mm**2.d0
     
copt      mperp_p = dsqrt(mperp_p2)
copt      mperp_m = dsqrt(mperp_m2)
     
c     -----------------------------------------------------------------
copt      xp = (ktp/rs)*DEXP(yp)
copt      xm = (ktm/rs)*DEXP(ym)
copt      xf = xp + xm

ctest      write(*,*) 'phip, phim: ', phip, phim
ctest      write(*,*) 'yp, ym: ', yp, ym
ctest      write(*,*) 'xp, xm, xf: ', xp, xm, xf

ctest      write(*,*) 'ktp, ktm, phip, phim, pt: ', ktp,ktm,phip,phim,pt


c     ==================================================================
c     Boson variables
c     -----------------------------------------------------------------
copt      ptx = ktpx + ktmx
copt      pty = ktpy + ktmy
copt      pt  = DSQRT(ptx**2.d0 + pty**2.d0)
copt      pt2 = pt**2.d0
      pt = DSQRT(pt2Var)
      pt2 = pt2Var
      y = yVar
      m2 = m2Var
copt      y  = DLOG(xf*(rs/DSQRT(pt2 + M2)))


      x1 = (DSQRT(M2 + pt**2.d0)/RS)*DEXP(y)
      x2 = (DSQRT(M2 + pt**2.d0)/RS)*DEXP(-y)

copt      x1 = (ktp/rs)*DEXP(yp) + (ktm/rs)*DEXP(ym)
copt      x2 = (ktp/rs)*DEXP(-yp) + (ktm/rs)*DEXP(-ym)


copt      m2 = mperp_p**2.d0 + mperp_m**2.d0 
copt     &       + 2.d0*mperp_p*mperp_m*DCOSH(yp - ym) - pt2

ctest      WRITE(*,*) 'Kinematics: ', y, pt, M, x1, x2
ctest      write(*,*) 'M2, M, x1, x2: ', m2,mVar, x1, x2
     
      mVar = DSQRT(m2) 

      varJacobian = (2.d0/rs)*DSQRT(M2 + pt2)*DCOSH(y)
      preIntegral = (x1/(x1 + x2))*varJacobian
      
ctest      write(*,*) 'Pre-integral: ', varJacobian, preIntegral


c     =================================================================
c     Delimiting boundaries
c     -----------------------------------------------------------------f1
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
      

      HadronicCrossSection =  InterpolateGrid(y,pt,mVar)  
      Result = preIntegral*DileptonDecay(Mvar)*HadronicCrossSection
      
copt      SigTot = xp*xm*ktp*ktm*Result
      SigTot = Result
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

      Branch = 3.3d0/100.d0
     
      InvariantMassDist = (1.d0/pi)*
     & ((M*DecayWidth)/((M2 - Mz2)**2.d0 + (M*DecayWidth)**2.d0))

      Result = InvariantMassDist*Branch
      DileptonDecay = Result

ctest      write(*,*)'Decay: ',Result,DecayWidth,Branch,InvariantMassDist
      RETURN 
      END 



      function InterpolateGrid(yVar,ptVar,mVar)
      implicit none
      integer, parameter :: nPoints = 15
      integer, parameter :: narg=3
      
      integer :: nent(3)
      double precision PartonL     evelGrid(nPoints,nPoints,nPoints)
      double precision ptGrid(nPoints), yGrid(nPoints), mGrid(nPoints)
      double precision ent(nPoints + nPoints + nPoints), arg(narg)
      double precision y, pt, yVar, ptVar,InterpolateGrid
      double precision m, mVar
      double precision DFINT
      character*2000 :: File
      
      common/GridArrays/ent,yGrid,ptGrid,mGrid,PartonLevelGrid


      File = "Grids/DatFiles/tst_grid.dat"

      call read_grid(File)
      
      y = yVar
      pt = ptVar
      m = mVar


      nent(1) = nPoints
      nent(2) = nPoints
      nent(3) = nPoints

      arg(1) = y
      arg(2) = pt
      arg(3) = m

      InterpolateGrid = DFINT(narg,arg,nent,ent,PartonLevelGrid)
ctest      write(*,*) 'Interpolated value: ', InterpolateGrid, y, pt, m
      return
      end 


      subroutine read_grid(OutputPath)
      implicit none
      integer, parameter :: iFile = 11
      integer, parameter :: nPoints = 15
      integer :: i, j, k
      double precision yVar, ptVar, mVar, PartonLevelVar
      double precision PartonLevelGrid(nPoints,nPoints,nPoints)
      double precision ptGrid(nPoints), yGrid(nPoints), mGrid(nPoints)
      double precision ent(nPoints + nPoints + nPoints)
      character*2000 :: OutputPath
      
      common/GridArrays/ent,yGrid,ptGrid,mGrid,PartonLevelGrid

      open(iFile,file=trim(OutputPath),status='old')
c      open(iFile,
c     *  file="Grids/DatFiles/kslinear_grid.dat",status="old")
      
      do i = 1, nPoints
            do j = 1, nPoints
                  do k = 1, nPoints
      read(iFile,*) yVar, ptVar, mVar, PartonLevelVar
      yGrid(i) = yVar
      ptGrid(j) = ptVar
      mGrid(k) = mVar
      PartonLevelGrid(i,j,k) = PartonLevelVar
ctest            write(*,*) 'Reading grid point:', i, j, k
ctest            write(*,*) 'Grid', yVar, ptVar, mVar, PartonLevelVar
                  end do
            end do
      end do
      close(iFile)
      do i = 1, nPoints
         ent(i) = yGrid(i)
      end do 

      do j = 1, nPoints
         ent(nPoints + j) = ptGrid(j)
      end do

      do k = 1, nPoints
         ent(nPoints + nPoints + k) = mGrid(k)
      end do
      
      end subroutine








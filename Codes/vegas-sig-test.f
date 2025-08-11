c     TO COMPILE: gfortran -fdefault-real-8 -fdefault-double-8 -fno-align-commons vegas-sig-test.f 
      program VEGAS_DY_CS_test
      implicit none

      
      double precision sigma 
c     -----------------------------------------------------------------
c     =================================================================
    


c     =================================================================
c     VEGAS definitions
c     -----------------------------------------------------------------
      double precision avgi,sd,chi2a
      double precision x(6),wgt 
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
c     common blocks (for DY_CS code)
c     -----------------------------------------------------------------
      double precision y,y_min,y_max,dy,sum_y
      double precision pt,pt_min,pt_max,dpt,sum_pt
      double precision m,m_min,m_max,dm
      double precision sig_y(5000),sig_pt(5000),sig_m(5000)
      double precision sig_ypt1(5000),sig_ypt2(5000),sig_ypt3(5000)
      double precision sig_ypt4(5000),sig_ypt5(5000)
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



c     =================================================================
c     common blocks (for the grid)
c     -----------------------------------------------------------------
      double precision Lfy_grid(400,400), Ly_grid(400)
	double precision Lpt_grid(400), LENT(800)
	double precision Tfy_grid(400,400), Ty_grid(400)
	double precision Tpt_grid(400), TENT(800)
c	common/Lgrid/LENT,Ly_grid,Lpt_grid,Lfy_grid
c	common/Tgrid/TENT,Ty_grid,Tpt_grid,Tfy_grid
	common/gdT/TENT,Ty_grid,Tpt_grid,Tfy_grid
      common/gdL/LENT,Ly_grid,Lpt_grid,Lfy_grid
      

      external sigma




      call read_grid 




c     =================================================================
c     BINNING PARAMETERS (stored in common block)
c     =================================================================
      ny = 10                      
      y_min =  2.0d0
      y_max =  4.5d0
      dy = (y_max-y_min)/ny
c     -----------------------------------------------------------------
      npt = 100             
      pt_min = 0.d0
      pt_max = 150.d0
      dpt = (pt_max-pt_min)/npt
c     -----------------------------------------------------------------      
      nm = 60 
      m_min = 60.d0
      m_max = 120.d0
      dm = (m_max-m_min)/nm

      
      nprn  = 0
      ncall = 10000 
      itmx  = 5
      iDoHist = 0
      call VEGAS(6,sigma,avgi,sd,chi2a) 



      ncall = 1000000
      itmx  = 20
      iDoHist = 1
      call VEGAS1(6,sigma,avgi,sd,chi2a)
      print*, avgi, "+-", sd

c     =================================================================
c     FILES WITH RESULTS
c     -----------------------------------------------------------------
      open(unit=41,file='test/dsig_dy.dat',status='unknown')
      open(unit=42,file='test/dsig_dpt.dat',status='unknown')
      open(unit=43,file='test/dsig_dm.dat',status='unknown')
      open(unit=21,file='test/dsig_dydpt_y2p25.dat',status='unknown')
      open(unit=22,file='test/dsig_dydpt_y2p75.dat',status='unknown')
      open(unit=23,file='test/dsig_dydpt_y3p25.dat',status='unknown')
      open(unit=24,file='test/dsig_dydpt_y3p75.dat',status='unknown')
      open(unit=25,file='test/dsig_dydpt_y4p25.dat',status='unknown')

c     =================================================================
c     building spectra
c     =================================================================
      sum_y = 0.d0
      do iy = 1, ny

      y = y_min + iy*dy - dy/2.d0
      write(41,1001) y, sig_y(iy)
      sum_y = sum_y + sig_y(iy)*dy
      enddo

      write(*,*) 
      write(*,*) 'sum_y'
      write(*,*) sum_y
c     -----------------------------------------------------------------

      sum_pt = 0.d0
      do ipt = 1, npt
      pt = pt_min + ipt*dpt - dpt/2.d0

      write(42,1001) pt, sig_pt(ipt)
      write(21,1001) pt, sig_ypt1(ipt)
      write(22,1001) pt, sig_ypt2(ipt)
      write(23,1001) pt, sig_ypt3(ipt)
      write(24,1001) pt, sig_ypt4(ipt)
      write(25,1001) pt, sig_ypt5(ipt)


      sum_pt = sum_pt + sig_pt(ipt)*dpt
      enddo

      write(*,*) 
      write(*,*) 'sum_pt'
      write(*,*) sum_pt
      close(42)
      close(41)


c     -----------------------------------------------------------------
      do im = 1,nm 

      m = m_min + im*dm - dm/2.d0

      write(43,1001) m, sig_m(im)
      enddo
      close(43)




1001  format(1x,f8.4,1x,e12.4)
      end 


c     =================================================================
c     Integration function
c     -----------------------------------------------------------------  
      function sigma(x,wgt)
      implicit none
      
      double precision sigma, x(6), wgt, jac
      double precision yp,ym,ktp,ktm,phip,phim
      double precision pi,kp_min,kp_max,km_min,km_max
      double precision yp_min,yp_max,ym_min,ym_max
      double precision aintegrand, tmp
c     -----------------------------------------------------------------
      integer ncall,itmx,nprn,ndev,it,ndo
      integer ndmx,mds
      integer ncall1,ncall2,itmx1,itmx2,ihist
      double precision xl,xu,acc,si,swgt,schi,xi
      double precision alph
c     -----------------------------------------------------------------
      logical isnan 
      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/bveg3/alph,ndmx,mds
c     =================================================================
c     PHASE SPACE
c     =================================================================
      pi  = 4.d0*datan(1.d0) 
      kp_max = 200.d0
      kp_min = 20.d0
      km_max = 200.d0 
      km_min = 20.d0
      yp_max = 4.5d0
      yp_min = 2.0d0
      ym_max = 4.5d0 
      ym_min = 2.d0
c     -----------------------------------------------------------------
      yp   = yp_min + (yp_max - yp_min)*x(1)
      ym   = ym_min + (ym_max - ym_min)*x(2)
      ktp   = kp_min + (kp_max - kp_min)*x(3)
      ktm   = km_min + (km_max - km_min)*x(4)
      phip = 2.d0*pi*x(5)
      phim = 2.d0*pi*x(6)
c     =================================================================
c     jacobian: x(n) ----> phase space
c     =================================================================
      jac = (yp_max - yp_min)*(ym_max - ym_min)*(kp_max - kp_min)
     &      *(km_max - km_min)*((2.d0*pi)**2.d0)

      tmp=(wgt*jac)/(itmx)
      
      
      call DY_CS(aintegrand,yp,ym,ktp,ktm,phip,phim,tmp)
      
      sigma = jac*aintegrand
c      write(*,*) yp,ktp,phip
      return
      end 
      

      subroutine DY_CS(aintegrand,yp,ym,ktp,ktm,phim,phip,wgt)
      double precision aintegrand, yp, ym, ktp, ktm,phim,phip,wgt
      double precision rs,s, ktpx,ktpy,ktmx,ktmy
      double precision ax,ay,ax2,ay2
      double precision pi, alfem,Ml,M,Ml2,M2,AW,SIN2,SIN4
      double precision ga,gv,ga2,gv2
      double precision Gamma,Gamma2,propsqrd,CPL
      double precision ktp2,ktm2,xp,xm,xf, U,U2, DL,DT
      double precision y,pt, pt2, ptx,pty
      double precision Ldsdydpt,Tdsdydpt,L_dsdydpt, T_dsdydpt
      double precision L_aintegrand , T_aintegrand
      double precision aintegral, units
      double precision mperp_p,mperp_m,mp,mm,mperp_p2,mperp_m2,ml_perp
      double precision yd, Inv_M, branch, mom,x1


c     =================================================================
c     common blocks (for DY_CS code)
c     -----------------------------------------------------------------
      double precision y_min,y_max,dy
      double precision pt_min,pt_max,dpt
      double precision m_min,m_max,dm
      double precision sig_y(5000),sig_pt(5000),sig_m(5000)
      double precision sig_ypt1(5000),sig_ypt2(5000),sig_ypt3(5000)
      double precision sig_ypt4(5000),sig_ypt5(5000)
      
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

c     =================================================================
c     common blocks (for the grid)
c     -----------------------------------------------------------------
	double precision Lfy_grid(400,400), Ly_grid(400)
	double precision Lpt_grid(400), LENT(800)
	double precision Tfy_grid(400,400), Ty_grid(400)
	double precision Tpt_grid(400), TENT(800)
c     -----------------------------------------------------------------
	common/gdT/TENT,Ty_grid,Tpt_grid,Tfy_grid
      common/gdL/LENT,Ly_grid,Lpt_grid,Lfy_grid

      
      external Ldsdydpt,Tdsdydpt
      
     
c     =================================================================
c     energy in CM in GeV
c     -----------------------------------------------------------------
      rs = 13000.d0      
      s = rs**2.d0
c     =================================================================


c     =================================================================
c     lepton pair variables
c     -----------------------------------------------------------------
      mp = 105.6d-3 !GeV
      mm = 105.6d-3 !GeV



      
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
     

   

      ax = xp*ktmx - xm*ktpx 
      ay = xp*ktmy - xm*ktpy 

      ax2 = ax**2.d0
      ay2 = ay**2.d0
c     -----------------------------------------------------------------
      xp = DSQRT(ktp2/s)*DEXP(yp)
      xm = DSQRT(ktm2/s)*DEXP(ym)
      xf = xp + xm
ctest      write(*,*) xp,xm,xf
      
      
c     ==================================================================
c     Boson variables
c     -----------------------------------------------------------------
      ptx = ktpx + ktmx
      pty = ktpy + ktmy
      pt  = DSQRT(ptx**2.d0 + pty**2.d0)
      pt2 = pt**2.d0
      y   = DLOG(xf*DSQRT(s/(pt2 + M2)))
      mom = dsqrt(pt2 + m2)
      x1  = (mom / rs)*DEXP(Y)
c      write(*,*) y,pt,pt2, yp-ym
      
c     =================================================================

c     =================================================================
c     Definitions
c     -----------------------------------------------------------------
      pi = 4.d0*datan(1.d0)
      alfem = 1.D0/137.D0  
      yd = dcosh(yp-ym)
      ml_perp = mperp_p2 + mperp_m2 + 2.d0*mperp_p*mperp_m*yd
      Ml2  = ml_perp - pt2      
c      Ml2   =  (ax2 + ay2)/(xm*xp)
      M    = 91.1896D0       ! Z mass 
      Ml   = DSQRT(Ml2)
      M2   = M**2.d0 
     
c     -----------------------------------------------------------------
      aw    = 0.5061D0            !WEINBERG ANGLE
      SIN2  = DSIN(AW)**2.D0
      SIN4   = DSIN(AW)**4.d0
c     -----------------------------------------------------------------
      ga  = - 0.5d0
      gv  = 2.d0*sin2 - 0.5d0
      ga2 = ga**2.d0
      gv2 = gv**2.d0
c     -----------------------------------------------------------------
      Gamma  = 2.4952d0 !GeV
      Gamma2 = Gamma**2.d0
      propsqrd = 1.d0/(((Ml2-M2)**2.d0) + M2*Gamma2) !GeV-4
      CPL = (((alfem*Ml)/(dsin(2.d0*aw)**2.d0))*(gv2+ga2)*propsqrd) 
     &       /(4.d0*(pi**2.d0))            !Gev-3
c      !previously, it was ML2, now I changed to ML
c      branch = 0.0336d0
c      Gamma = ((alfem*Ml)/(6.d0*((dsin(2.d0*aw))**2.d0)))*(
c     & (160.d0/3.d0)*(sin2**2.d0) - 40.d0*sin2 + 21.d0)      
c      Gamma2 = Gamma**2.d0
c      Inv_M  = (Ml/pi)*((Ml*Gamma)/(((Ml2-M2)**2.d0) + Ml2*Gamma2))
      
c     -----------------------------------------------------------------   
      U  = xp/xf
      U2 = U**2.D0 
c     -----------------------------------------------------------------
      DT = 4.D0*(U2 + (1.D0-U)**2.D0)
      DL = 8.D0*U*(1.D0-U)
c     =================================================================
ctest      write(*,*) xp, xm, ktp, ktm
     
      
      

c     =================================================================
c     Delimiting boundaries
c     -----------------------------------------------------------------
      
      if(xp.gt.1.d0.or.xm.gt.1.d0) then
            aintegrand = 0.d0
ctest            write(*,*) "acima do x"
            go to 101
      endif
c     -----------------------------------------------------------------
      if(ktp.lt.20.d0.or.ktm.lt.20.d0) then
            aintegrand = 0.d0
ctest            write(*,*) "abaixo do Kt"
            go to 101
      endif

      if(Ml.lt.60.d0.or.Ml.gt.120.d0) then
            aintegrand = 0.d0
ctst            write(*,*) "fora da massa"
            go to 101
      endif


c     ==================================================================
c     Boson differential cross section
c     -----------------------------------------------------------------
      L_dsdydpt = Ldsdydpt(y,pt)
      T_dsdydpt = Tdsdydpt(y,pt)

     
      L_aintegrand = CPL*DL*L_dsdydpt
c      L_aintegrand = 0.d0
      T_aintegrand = CPL*DT*T_dsdydpt
c      T_aintegrand = 0.d0
      aintegral = (xf/x1)*(L_aintegrand + T_aintegrand)
ctest      write(*,*) L_aintegrand, T_aintegrand
ctest      write(*,*) 'aintegrand=', aintegrand, xf, DL, DT
      
c     =================================================================
ctest      write(*,*) CPL*DT, aintegral
ctest      write(*,*) pt, y, aintegrand, wgt, dy
ctest      write(*,*) xp, xm, ktp, ktm
ctest      write(*,*) Inv_M, L_aintegrand, T_aintegrand

     
c     =================================================================
c     The grids values are given in GeV-2
c     -----------------------------------------------------------------c
      units = 0.389d9 !GeV-2 to pb
      if(iDoHist.eq.1.or.iDoHist.eq.0) then 
            aintegrand = aintegral*units*ktp*ktm!output in pb
      endif






c     =================================================================
c     indices
c     -----------------------------------------------------------------
      iy = idint((y-y_min)/dy) + 1
c     -----------------------------------------------------------------
      ipt = idint((pt-pt_min)/dpt) + 1
c     -----------------------------------------------------------------
      im = idint((Ml-m_min)/dm) + 1
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
            sig_y(iy) = sig_y(iy) + aintegrand*wgt/dy
c            write(*,*) sig_y(iy), aintegrand*wgt/dy
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
            sig_pt(ipt) = sig_pt(ipt) + aintegrand*wgt/dpt
            endif
      endif

c     -----------------------------------------------------------------
      if(Ml.gt.m_min.and.Ml.lt.m_max) then
            if(im.gt.0.and.im.le.nm) then
                  sig_m(im) = sig_m(im) + aintegrand*wgt/dm
            endif 
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.y_min.and.y.lt.y_max) then 
                  if(y.gt.2.0d0.and.y.lt.2.5d0) then
            sig_ypt1(ipt) = sig_ypt1(ipt) + aintegrand*wgt/(dpt*dy)
                  endif
                  endif
            endif 
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.2.5d0.and.y.lt.3.0d0) then
            sig_ypt2(ipt) = sig_ypt2(ipt) + aintegrand*wgt/(dpt*dy)!*wgt/dpt
                  endif
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.3.0d0.and.y.lt.3.5d0) then
            sig_ypt3(ipt) = sig_ypt3(ipt) + aintegrand*wgt/(dpt*dy)!*wgt/dpt
                  endif
            endif
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.3.5d0.and.y.lt.4.0d0) then
            sig_ypt4(ipt) = sig_ypt4(ipt) + aintegrand*wgt/(dpt*dy)!*wgt/dpt
                  endif
            endif
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.4.0d0.and.y.lt.4.5d0) then
            sig_ypt5(ipt) = sig_ypt5(ipt) + aintegrand*wgt/(dpt*dy)!*wgt/dpt
                  endif
            endif
      endif

      endif
c     =================================================================



101   continue
      return
      end 

      


























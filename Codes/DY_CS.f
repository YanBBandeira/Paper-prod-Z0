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
      double precision yd


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
	 double precision Lfy_grid(10,100), Ly_grid(10)
	double precision Lpt_grid(100), LENT(110)
	double precision Tfy_grid(10,100), Ty_grid(10)
	double precision Tpt_grid(100), TENT(110)
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
      
      
c     ==================================================================
c     Boson variables
c     -----------------------------------------------------------------
      ptx = ktpx + ktmx
      pty = ktpy + ktmy
      pt  = DSQRT(ptx**2.d0 + pty**2.d0)
      pt2 = pt**2.d0
      y   = DLOG(xf*DSQRT(s/(pt2 + M2)))
      
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
      Ml2   =  (ax2 + ay2)/(xm*xp)
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
      !previously, it was ML2, now I changed to ML
ctest      write(*,*) CPL, Ml
c     -----------------------------------------------------------------   
      U  = xp/xf
      U2 = U**2.D0 
c     -----------------------------------------------------------------
      DT = 4.D0*(U2 + (1.D0-U)**2.D0)
      DL = 8.D0*U*(1.D0-U)
c     =================================================================
ctest      write(*,*) xp, xm, ktp, ktm
     
      
      

c     ==================================================================


c     ==================================================================
c     Boson differential cross section
c     -----------------------------------------------------------------
      L_dsdydpt = Ldsdydpt(y,pt)
      T_dsdydpt = Tdsdydpt(y,pt)

      
      L_aintegrand = CPL*DL*L_dsdydpt
c      L_aintegrand = 0.d0
      T_aintegrand = CPL*DT*T_dsdydpt
c      T_aintegrand = 0.d0
      aintegral = L_aintegrand + T_aintegrand 
ctest      write(*,*) L_aintegrand, T_aintegrand
ctest      write(*,*) 'aintegrand=', aintegrand, xf, DL, DT
      
c     =================================================================
ctest      write(*,*) CPL*DT, aintegral
ctest      write(*,*) pt, y, aintegrand, wgt, dy
ctest      write(*,*) xp, xm, ktp, ktm
c     =================================================================
c     Delimiting boundaries
c     -----------------------------------------------------------------
      if(xp.gt.1.d0.or.xm.gt.1.d0) then
            aintegrand = 0.d0
            go to 101
      endif
c     -----------------------------------------------------------------
      if(ktp.lt.20.d0.or.ktm.lt.20.d0) then
            aintegrand = 0.d0
            go to 101
      endif

      if(Ml.lt.60.d0.or.Ml.gt.120.d0) then
            aintegrand = 0.d0
            go to 101
      endif
     
c     =================================================================
c     The grids values are given in GeV-2
c     -----------------------------------------------------------------c
      units = 0.389d9 !GeV-2 to pb
      if(iDoHist.eq.1.or.iDoHist.eq.0) then 
            aintegrand = (aintegral*units*ktp*ktm)/(2.d0*pi*pt)
c            !! I'm not 100% sure that should be divide by 2*pi*pt
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
ctest            write(*,*) iy, ny
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
            sig_pt(ipt) = sig_pt(ipt) + aintegrand*wgt/dpt
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.2.0d0.and.y.lt.2.5d0) then
            sig_ypt1(ipt) = sig_ypt1(ipt) + aintegrand*wgt/dpt
                  endif
            endif 
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.2.5d0.and.y.lt.3.0d0) then
            sig_ypt2(ipt) = sig_ypt2(ipt) + aintegrand*wgt/dpt
                  endif
            endif
      endif
c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.3.0d0.and.y.lt.3.5d0) then
            sig_ypt3(ipt) = sig_ypt3(ipt) + aintegrand*wgt/dpt
                  endif
            endif
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.3.5d0.and.y.lt.4.0d0) then
            sig_ypt4(ipt) = sig_ypt4(ipt) + aintegrand*wgt/dpt
                  endif
            endif
      endif

c     -----------------------------------------------------------------
      if(pt.gt.pt_min.and.pt.lt.pt_max) then 
            if(ipt.gt.0.and.ipt.le.npt) then
                  if(y.gt.4.0d0.and.y.lt.4.5d0) then
            sig_ypt5(ipt) = sig_ypt5(ipt) + aintegrand*wgt/dpt
                  endif
            endif
      endif
      

c     -----------------------------------------------------------------
      if(Ml.gt.m_min.and.Ml.lt.m_max) then
            if(im.gt.0.and.im.le.nm) then
                  sig_m(im) = sig_m(im) + aintegrand*wgt/dm
            endif 
      endif
      

      endif
c     =================================================================



101   continue
      return
      end 

      



c     CODE TO EVALUATE PP -> Z0 + X      
c     TO COMPILE: gfortran -fdefault-real-8 -fdefault-double-8 -fno-align-commons PL-grids-ms.f dgauss.f CT14Pdf.f
      program PLgrids_MS 
      implicit none 
c     ------------------------------------------------------------------
      double precision y_min, y_max, pt_min, pt_max, dy, dpt
      double precision y(5000), pt(5000), fgrid_T, fgrid_L
      double precision dsig
c     ------------------------------------------------------------------
      integer  ny, npt,iy,ipt
c     ------------------------------------------------------------------
      external dsig 



c     ==================================================================
c     grid parameters
c     ------------------------------------------------------------------
      ny = 400
      npt = 400
      y_min  = 1.d0 
      y_max  = 5.d0 
      pt_min = dlog10(1.0d0)
      pt_max = dlog10(150.d0) 

c     ------------------------------------------------------------------
      dy = (y_max - y_min) / (ny - 1)
      dpt = (pt_max - pt_min) / (npt - 1)

c     ==================================================================
c     output files 
c     ------------------------------------------------------------------
      open(unit=21,file='L_z0_grid_MS_GBW.dat',status='unknown')
      open(unit=22,file='T_z0_grid_MS_GBW.dat',status='unknown')


c     ==================================================================
c     grid loop
c     ------------------------------------------------------------------
      do iy = 1, ny
      y(iy) = y_min + (iy - 1) * dy
            do ipt = 1, npt
                  pt(ipt) = 10.D0**(pt_min + (ipt - 1) * dpt) - 0.9d0
                  fgrid_T = dsig(y(iy),pt(ipt),0) !0 for T 
                  fgrid_L = dsig(y(iy),pt(ipt),1) !1 for L
                  write(*,*) y(iy), pt(ipt), fgrid_T, fgrid_L
                  write(21,100) y(iy), pt(ipt), fgrid_L
                  write(22,100) y(iy), pt(ipt), fgrid_T
            end do 
      end do 
100   format(2x,6(E10.4,2x))
      end program PLgrids_MS


c     ==================================================================
c     parton level cross section function
c     ------------------------------------------------------------------
      function dsig(y,ptq,ifl) 
      implicit none
c     ------------------------------------------------------------------
      double precision dsig,y,ptq,mom,xf,x2,rs,ujac
      double precision alfem,aw,sin2,M,pt,pt2,M2,pi,pi2
      double precision li,ls,pre, dgauss, pl_integrand
c     ------------------------------------------------------------------
      integer iflag, ifl
c     ------------------------------------------------------------------
      common/flag/iflag
      common/ctes/alfem,aw,sin2,M,xf,pt,pt2,M2,pi,pi2
      common/xbj/x2

c     ------------------------------------------------------------------
      external pl_integrand,dgauss 
c     ------------------------------------------------------------------
      iflag = ifl
      pi  = 4.d0*datan(1.d0)
      pi2 = pi**2.d0 
      alfem = 1.d0/137.d0
      sin2 = 0.23d0
      aw   = dasin(dsqrt(sin2))
c     ------------------------------------------------------------------
      rs = 13000.d0
      m  = 91.186d0 

c     ------------------------------------------------------------------
      m2  = m**2.d0 
      pt  = ptq
      pt2 = pt*2.d0 

c     ------------------------------------------------------------------
      mom = dsqrt(pt2 + m2)
      xf = (mom / rs)*DEXP(Y)
      x2 = (mom / rs)*DEXP(-Y)
c     ------------------------------------------------------------------

      li = xf 
      ls = 1.d0 

c     ------------------------------------------------------------------
      ujac = (2.d0/rs)*mom*DCOSH(y)
      pre  = 2.d0*pi*ptq!*(xf/(xf+x2))
      dsig = pre*dgauss(pl_integrand,li,ls,1.d-3)!*ujac
ctest      write(*,*) dsig,pre
      return 
      end 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 INTEGRAL OVER Z 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function pl_integrand(alf)
      implicit none
c     ------------------------------------------------------------------
      double precision alf,pl_integrand,z,x1,xf,hs,pt,pt2,M2,q
      double precision CT14Pdf,u,d,s,c,b,u_bar,d_bar,s_bar,c_bar,b_bar
      double precision mu,md,ms,mc,mb,mu2,md2,ms2,mc2,mb2
      double precision alfem,aw,sin2,M,pi,pi2
      double precision cfg,gfgaup,gfgadw,gfgvup,gfgvdw
      double precision mz,mq,cf,gfa,gfv
      double precision L_func,T_func
      double precision dsu_L,dsu_T,dsds_L,dsds_T,dsc_L,dsc_T,dsb_L,dsb_T
c     ------------------------------------------------------------------
      integer iflag
c     ------------------------------------------------------------------
      character Tablefile*40
      data Tablefile/'CT14LL.pds'/
c     ------------------------------------------------------------------
      common/flag/iflag
      common/ctes/alfem,aw,sin2,M,xf,pt,pt2,M2,pi,pi2
c     ------------------------------------------------------------------
      z = alf  
      x1 = xf/z 
      

      hs = dsqrt(pt2 + (1.d0-xf)*M2) 
c     ------------------------------------------------------------------
      if(hs.le.1.3d0) then
            q = 1.3d0 
      else 
            q = hs
      end if
c     ------------------------------------------------------------------
      if(x1.le.1.d0)then 


c     ==================================================================
c     Calling the ct14 pdf distribution 
c     ------------------------------------------------------------------
      call SetCT14(Tablefile)

      u     = CT14Pdf(1,x1,Q)     !u
      d     = CT14Pdf(2,x1,Q)     !d
      s     = CT14Pdf(3,x1,Q)     !s
      c     = CT14Pdf(4,x1,Q)     !c
      b     = CT14Pdf(5,x1,Q)     !b
      u_bar    = CT14Pdf(-1,x1,Q)     !u_bar
      d_bar    = CT14Pdf(-2,x1,Q)     !d_bar
      s_bar    = CT14Pdf(-3,x1,Q)     !s_bar 
      c_bar    = CT14Pdf(-4,x1,Q)     !c_bar
      b_bar    = CT14Pdf(-5,x1,Q)     !b_bar
c     ------------------------------------------------------------------
c     light quarks
      MU    = 0.14D0   
      MD    = 0.14D0 
      MS    = 0.14D0
c     heavy quarks      
      MC    = 1.4D0 
      MB    = 4.5d0
c     ------------------------------------------------------------------
      MU2 = MU**2.D0
      MD2 = MD**2.D0
      MS2 = MS**2.D0
      MC2 = MC**2.D0
      MB2 = MB**2.d0
c     ------------------------------------------------------------------
ctest      write(*,*) u,d,s,c,b,u_bar,d_bar,s_bar,c_bar,b_bar




c     ------------------------------------------------------------------
c     Charge and coupling
c     ------------------------------------------------------------------
      CFG   = DSQRT(alfem)/dsin(2.d0*aW)
      gfgaup = 0.5D0          ! FOR u,c,t quarks
      gfgadw = - 0.5d0          ! FOR d,s,b quarks
      gfgvup = 0.5D0 - (4.D0/3.D0)*SIN2   ! FOR u,c,t quarks
      gfgvdw = (2.D0/3.D0)*SIN2 - 0.5D0   ! FOR d,s,b quarks



c     ------------------------------------------------------------------
c     Evaluating each quark flavour
c     ------------------------------------------------------------------
      mz = M 
      mq = mu 
      cf = cfg
      gfa = gfgaup
      gfv = gfgvup
      call pl_func(L_func,T_func,pt,z,mz,mq,cf,gfa,gfv)
      dsu_L = (u+u_bar)*L_func  
      dsu_T = (u+u_bar)*T_func 
      
      
c     ------------------------------------------------------------------
c     quark d + quark s beacause both are down quarks with similar mass
      mz = m 
      mq = md 
      cf = cfg
      gfa = gfgadw
      gfv = gfgvdw
      call pl_func(L_func,T_func,pt,z,mz,mq,cf,gfa,gfv)
      dsds_L = (s+s_bar+d+d_bar)*L_func  
      dsds_T = (s+s_bar+d+d_bar)*T_func 

c     ------------------------------------------------------------------
      mz = m 
      mq = mc 
      cf = cfg
      gfa = gfgaup
      gfv = gfgvup
      call pl_func(L_func,T_func,pt,z,mz,mq,cf,gfa,gfv)
      dsc_L = (c+c_bar)*L_func  
      dsc_T = (c+c_bar)*T_func 
      

c     ------------------------------------------------------------------
      mz = m 
      mq = mb 
      cf = cfg
      gfa = gfgadw
      gfv = gfgvdw
      call pl_func(L_func,T_func,pt,z,mz,mq,cf,gfa,gfv)
      dsb_L = (b+b_bar)*L_func  
      dsb_T = (b+b_bar)*T_func 

ctest      write(*,*) dsu_L,dsu_T,dsds_L,dsds_T,dsc_L,dsc_T,dsb_L,dsb_T
ctest      write(*,*) L_func, T_func,m


c     ------------------------------------------------------------------
c     Flag and the integrand result 0 = T and 1 = L
c     ------------------------------------------------------------------

      if(iflag.eq.0) then 
            pl_integrand = (dsu_T + dsds_T + dsc_T + dsb_T)/(z**2.d0)
      end if 
      if(iflag.eq.1) then
            pl_integrand = (dsu_L +dsds_L + dsc_L + dsb_L)/(z**2.d0)
      end if 
c     ------------------------------------------------------------------

      else 
      pl_integrand = 0.d0
      end if
      return 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 PARTON LEVEL CROSS SECTIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine pl_func(L_func,T_func,pt,z,mz,mq,cf,gfa,gfv)
      implicit none
c     ------------------------------------------------------------------
      double precision L_func,T_func 
      double precision tv,ta,lv,la,pt,z,mz,mq,cf,gfa,gfv,eps,eps2
      double precision pi,pi2,li,ls,dgauss
      double precision e1,e2, ptq, epps,zz
      double precision pre_tv,pre_lv,pre_ta,pre_la,ee1,ee2,cf2,gfv2,
     &                 gfa2,z2,z4,mf2,M2
      
      
c     ------------------------------------------------------------------
      common/ints/pre_tv,pre_lv,pre_ta,pre_la,cf2,gfv2,gfa2,z2,z4,mf2,M2             
      common/ees/  ptq, epps,zz
c     ------------------------------------------------------------------
      external e1,e2,tv_int,ta_int,lv_int,la_int 
c     ------------------------------------------------------------------

      pi = 4.d0*datan(1.d0)
      pi2 = pi**2.d0

c     ------------------------------------------------------------------

      cf2 = cf**2.d0
      gfv2= gfv**2.d0
      gfa2= gfa**2.d0
      z2  = z**2.d0
      z4  = z**4.d0
      mf2 = mq**2.d0
      M2  = mz**2.d0
c     ------------------------------------------------------------------

      eps2 = (1.d0 - z)*M2 + z2*mf2
      eps  = dsqrt(eps2)
      ptq  = pt 
      epps = eps 
      zz   = z

ctest      write(*,*) ee1, ee2
c     ------------------------------------------------------------------

      pre_tv = (cf2*gfv2)/(2.d0*pi2) 
      pre_ta = (cf2*gfa2)/(2.d0*pi2) 
      pre_lv = (cf2*gfv2)/(4.d0*pi2)
      pre_la = (cf2*gfa2)/(4.d0*pi2)
c     ------------------------------------------------------------------
      li = 0.d0
      ls = 1.d0
c     ------------------------------------------------------------------
c     UGD integration
c     ------------------------------------------------------------------
      tv    = dgauss(tv_int,li,ls,1.d-4) 
      ta    = dgauss(ta_int,li,ls,1.d-4)
      lv    = dgauss(lv_int,li,ls,1.d-4)
      la    = dgauss(la_int,li,ls,1.d-4)


ctest      write(*,*) tv,ta,lv,la
c     ------------------------------------------------------------------
      T_func = tv + ta 
      L_func = lv + la


      end 



c     ==================================================================
c     Tranverse Vector contribution
c     ==================================================================
      function tv_int(u)
      implicit none
c     ------------------------------------------------------------------
      double precision tv_int,tv_t1,tv_t2,akt2,ugd,x,z,fkt, u
      double precision e1,e2, ptq, epps,zz, akt
      double precision pre_tv,pre_lv,pre_ta,pre_la,ee1,ee2,cf2,gfv2,
     &                 gfa2,z2,z4,mf2,M2
c     ------------------------------------------------------------------
      common/ints/pre_tv,pre_lv,pre_ta,pre_la,cf2,gfv2,gfa2,z2,z4,mf2,M2
      common/ees/  ptq, epps,zz
c     ------------------------------------------------------------------
      external ugd
c     ------------------------------------------------------------------
      akt2 = u/(1.D0-u)
      akt  = dsqrt(akt2)
      z = dsqrt(z2)
      fkt = ugd(akt)
c     ------------------------------------------------------------------
      akt = dsqrt(akt2)
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      tv_t1 = z4*mf2*ee1
      tv_t2 = (1.d0+((1.d0-z)**2.d0))*ee2
      tv_int= 0.5d0*fkt*pre_tv*(tv_t1 + tv_t2)*((1.d0 + akt2)**2.d0)
ctest      write(*,*) ee1,ee2
      return 
      end



c     ==================================================================
c     Tranverse Axial contribution
c     ========================================ee1==========================
      function ta_int(u)
      implicit none
c     ------------------------------------------------------------------
      double precision ta_int,ta_t1,ta_t2,akt2,ugd,x,z, fkt, u
      double precision e1,e2, ptq, epps,zz, akt
      double precision pre_tv,pre_lv,pre_ta,pre_la,ee1,ee2,cf2,gfv2,
     &                 gfa2,z2,z4,mf2,M2
c     ------------------------------------------------------------------
      common/ints/pre_tv,pre_lv,pre_ta,pre_la,cf2,gfv2,gfa2,z2,z4,mf2,M2
      common/ees/  ptq, epps,zz
c     ------------------------------------------------------------------
      external ugd
c     ------------------------------------------------------------------
      akt2 = u/(1.D0-u)     
      akt  = dsqrt(akt2) 
      z = dsqrt(z2)
      fkt = ugd(akt)
c     ------------------------------------------------------------------
      akt = dsqrt(z2*akt2)
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      ta_t1 = z2*mf2*((2.d0-z)**2.d0)*ee1 
      ta_t2 = (1.d0+((1.d0-z)**2.d0))*ee2
      ta_int= 0.5d0*fkt*pre_ta*(ta_t1 + ta_t2)*((1.d0 + akt2)**2.d0)
      return
      end 



c     ==================================================================
c     Longitudinal Vector contribution
c     ==================================================================
      function lv_int(u)
      implicit none
c     ------------------------------------------------------------------
      double precision lv_int,akt2,ugd,x,z,fkt, u
      double precision e1,e2, ptq, epps,zz, akt
      double precision pre_tv,pre_lv,pre_ta,pre_la,ee1,ee2,cf2,gfv2,
     &                 gfa2,z2,z4,mf2,M2, lv_t1
c     ------------------------------------------------------------------
      common/ints/pre_tv,pre_lv,pre_ta,pre_la,cf2,gfv2,gfa2,z2,z4,mf2,M2
      common/ees/  ptq, epps,zz
c     ------------------------------------------------------------------
      external ugd 
c     ------------------------------------------------------------------
      akt2 = u/(1.D0-u)
      akt  = dsqrt(akt2)
      z = dsqrt(z2)
      fkt = ugd(akt)
c     ------------------------------------------------------------------
      akt = dsqrt(akt2)
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      lv_t1 = 4.d0*((1.d0-z)**2.d0)*M2*ee1
      lv_int= 0.5d0*fkt*pre_lv*lv_t1*((1.d0 + akt2)**2.d0)
      return 
      end


c     ==================================================================
c     Tranverse Axial contribution
c     ==================================================================

      function la_int(u)
      implicit none
c     ------------------------------------------------------------------
      double precision la_int,la_t1,la_t2,akt2,ugd,x,z,fkt, u
      double precision e1,e2, ptq, epps,zz, akt
      double precision pre_tv,pre_lv,pre_ta,pre_la,ee1,ee2,cf2,gfv2,
     &                 gfa2,z2,z4,mf2,M2
c     ------------------------------------------------------------------
      common/ints/pre_tv,pre_lv,pre_ta,pre_la,cf2,gfv2,gfa2,z2,z4,mf2,M2
      common/ees/  ptq, epps,zz
c     ------------------------------------------------------------------
      external ugd 
c     ------------------------------------------------------------------
      akt2 = u/(1.D0-u)
      akt  = dsqrt(akt2)
      z = dsqrt(z2)
      fkt = ugd(akt)
c     ------------------------------------------------------------------
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      la_t1 = 4.d0*((z2*mf2 + ((1.d0-z)**2.d0)*M2)**2.d0)*ee1 
      la_t2 = 4.d0*z2*mf2*ee2 
      la_int= 0.5d0*fkt*pre_la*((la_t1/M2) 
     *        + (la_t2/M2))*((1.d0 + akt2)**2.d0)
      return 
      end 

      
     
c     ==================================================================
c     The aux functions
c     ==================================================================c     ------------------------------------------------------------------
c     In what follows I'll define:
c     -->      tau2 =  (p - zk)^2 !for p and k as vectors
c     -->      eta  =  p . (p - zk) !for p and k as vectors and . is the dot product
c     moreover, it'll have the double precision definition in front of the 
c     function definitions to stand out this part of the code.
c     ------------------------------------------------------------------

      
      function e1(kt,pt,zz,epps)
      implicit none 
c     ------------------------------------------------------------------
      double precision e1,kt,pt,zz,epps 
      double precision k,p,z,eps 
      double precision dgauss,li,ls,e1_int,pi
c     ------------------------------------------------------------------
      common/e1var/k,p,z,eps
c     ------------------------------------------------------------------
      external e1_int,dgauss
c     ------------------------------------------------------------------
      pi = 4.d0*datan(1.d0)
      li = 0.d0
      ls = 2.d0*pi
      k = kt
      p = pt
      z = zz
      eps = epps
      e1 = dgauss(e1_int,li,ls,1.d-4)
      return 
      end 

c     ==================================================================
c     ==================================================================

      function e2(kt,pt,zz,epps)
      implicit none
c     ------------------------------------------------------------------
      double precision e2,kt,pt,zz,epps 
      double precision k,p,z,eps 
      double precision dgauss,li,ls,e2_int,pi
c     ------------------------------------------------------------------
      common/e2var/k,p,z,eps
c     ------------------------------------------------------------------
      external e2_int,dgauss 
c     ------------------------------------------------------------------
      pi = 4.d0*datan(1.d0)
      li = 0.d0
      ls = 2.d0*pi
      k = kt
      p = pt
      z = zz
      eps = epps
      e2 = dgauss(e2_int,li,ls,1.d-4)
      return 
      end 

c     ==================================================================
c     ==================================================================

      function e1_int(atheta)
      implicit none 
c     ------------------------------------------------------------------
      double precision e1_int,atheta
      double precision k,p,z,eps 
      double precision k2,z2,eps2,p2,tau2
      double precision e1_int_t1,e1_int_t2,e1_int_t3
c     ------------------------------------------------------------------
      common/e1var/k,p,z,eps
c     ------------------------------------------------------------------

      k2    = k**2.d0
      z2    = z**2.d0 
      eps2  = eps**2.d0
      p2    = p**2.d0
c     ------------------------------------------------------------------

      tau2 = p2 + z2*k2 - 2.d0*z*p*k*dcos(atheta)

c     ------------------------------------------------------------------

      e1_int_t1 = 1.d0/((tau2+eps2)**2.d0)
      e1_int_t2 = 2.d0/((p2+eps2)*(tau2+eps2))
      e1_int_t3 = 1.d0/((p2+eps2)**2.d0)

c     ------------------------------------------------------------------

      e1_int = 0.5d0*(e1_int_t1 - e1_int_t2 + e1_int_t3)
ctest      write(*,*) k2,z2,eps2,p2,tau
      return 
      end 

c     ==================================================================
c     ==================================================================

      function e2_int(atheta) 
      implicit none
c     ------------------------------------------------------------------
      double precision e2_int,atheta
      double precision k,p,z,eps
      double precision k2,z2,eps2,p2,tau2, eta
      double precision e2_int_t1,e2_int_t2,e2_int_t3
c     ------------------------------------------------------------------
      common/e2var/k,p,z,eps
c     ------------------------------------------------------------------

      k2    = k**2.d0
      z2    = z**2.d0 
      eps2  = eps**2.d0
      p2    = p**2.d0

c     ------------------------------------------------------------------

      tau2 = p2 + z2*k2 - 2.d0*z*k*p*dcos(atheta)
      eta  = p2 - z*p*k*dcos(atheta)

c     ------------------------------------------------------------------

      e2_int_t1 = tau2/((tau2+eps2)**2.d0)
      e2_int_t2 = eta/((p2+eps2)*(tau2+eps2))
      e2_int_t3 = p2/((p2+eps2)**2.d0)

c     ------------------------------------------------------------------

      e2_int = 0.5d0*(e2_int_t1 - 2.d0*e2_int_t2 + e2_int_t3)
      return 
      end



c     ==================================================================
c     Call the UGD model - GBW
c     ==================================================================
      !function ugd(akt)
      !implicit none
c     ------------------------------------------------------------------
      !double precision ugd,x,akt,akt2,x2,pi
      !double precision x0,lamb,sig0,R02,pre_ugd,exp_ugd
c     ------------------------------------------------------------------
      !common/xbj/x2
c     ------------------------------------------------------------------
      !pi  = 4.d0*datan(1.d0)
      !x0    = 3.04d-4
      !lamb  = 0.288d0
      !sig0  = 23.03/0.389D0  !mB to GeV^-2
      
      !x   = x2 
      !akt2  = akt**2.d0
      
      !R02     = (x/x0)**lamb
      !exp_ugd = akt2*R02
      !pre_ugd = (sig0*R02)/(4.d0*(pi))
   
      !ugd = pre_ugd*dexp(-exp_ugd)
ctest      write(*,*) ugd,x,akt     
      !return 
      !end


      FUNCTION ugd(akt) !WW
      implicit none

      DOUBLE PRECISION ugd,akt,akt2,alphas,F_WW,xbj,x2,pre_ugd,pi

      COMMON/xbj/x2

      EXTERNAL F_WW
      pi  = 4.d0*datan(1.d0)
      alphas = 0.12d0
      akt2 = akt**2.d0
      xbj = x2 

      pre_ugd = 4.d0*pi*alphas/3.d0
      ugd = (pre_ugd*F_WW(xbj,akt))/akt2
ctest      write(*,*) ugd, pre_ugd, F_WW(xbj,akt), akt4
      RETURN 
      END 

      FUNCTION F_WW(x,kt) !\mathcal{F} for WW 
      implicit none

      DOUBLE PRECISION F_WW,kt,kt2,k0,k02,
     &                 N1, lamb, b, x 

      N1 = 0.889d0 
      lamb = 0.29d0 
      k0 = 1.d0 
      b = 1.d0 

      kt2 = kt**2.d0
      k02 = k0**2.d0

      if (kt.ge.k0) then 
         F_WW = (N1/k02)*((1.d0-x)**7.d0)*(((x**lamb)*(kt2/k02))**(-b))
      end if 
      if (kt.lt.k0) then
         F_WW = (N1/k02)*((1.d0-x)**7.d0)*(x**(-b*lamb))
      end if
ctest      write(*,*) F_WW, kt,k0,x
      RETURN 
      END 


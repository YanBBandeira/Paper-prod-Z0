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
      tv    = dgauss(tv_int,li,ls,1.d-8) 
      ta    = dgauss(ta_int,li,ls,1.d-8)
      lv    = dgauss(lv_int,li,ls,1.d-8)
      la    = dgauss(la_int,li,ls,1.d-8)


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
      fkt = ugd(z2*akt2)
c     ------------------------------------------------------------------
      akt = dsqrt(akt2)
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      tv_t1 = z4*mf2*ee1
      tv_t2 = (1.d0+((1.d0-z)**2.d0))*ee2
      tv_int= fkt*pre_tv*(tv_t1 + tv_t2)*((1.d0 + akt2)**2.d0)
ctest      write(*,*) ee1,ee2
      return 
      end



c     ==================================================================
c     Tranverse Axial contribution
c     ==================================================================
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
      fkt = ugd(z2*akt2)
c     ------------------------------------------------------------------
      akt = dsqrt(z2*akt2)
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      ta_t1 = z2*mf2*((2.d0-z)**2.d0)*ee1 
      ta_t2 = (1.d0+((1.d0-z)**2.d0))*ee2
      ta_int= fkt*pre_ta*(ta_t1 + ta_t2)*((1.d0 + akt2)**2.d0)
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
      fkt = ugd(z2*akt2)
c     ------------------------------------------------------------------
      akt = dsqrt(akt2)
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      lv_t1 = 4.d0*((1.d0-z)**2.d0)*M2*ee1
      lv_int= fkt*pre_lv*lv_t1*((1.d0 + akt2)**2.d0)
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
      fkt = ugd(z2*akt2)
c     ------------------------------------------------------------------
      akt = dsqrt(akt2)
      ee1 = e1(akt,ptq,zz,epps)
      ee2 = e2(akt,ptq,zz,epps)
c     ------------------------------------------------------------------
      la_t1 = 4.d0*((z2*mf2 + ((1.d0-z)**2.d0)*M2)**2.d0)*ee1 
      la_t2 = 4.d0*z2*mf2*ee2 
      la_int= fkt*pre_la*((la_t1/M2) + (la_t2/M2))*((1.d0 + akt2)**2.d0)
      return 
      end 

      


c     ==================================================================
c     Call the UGD model - GBW
c     ==================================================================
      function ugd(akt2)
      implicit none
c     ------------------------------------------------------------------
      double precision ugd,xbj,akt2,akt4,x2,pi
      double precision x0,lamb,sig0,R02,pre_ugd,exp_ugd
c     ------------------------------------------------------------------
      common/xbj/x2
c     ------------------------------------------------------------------
      
      pi = 4.d0*datan(1.d0)
      x0    = 3.04d-4
      lamb  = 0.288d0
      sig0  = 23.03/0.389D0  !mB to GeV^-2
      
c     ------------------------------------------------------------------

      xbj   = x2 
      akt4  = akt2**2.d0
c     ------------------------------------------------------------------
      
      R02     = (xbj/x0)**lamb
      exp_ugd = akt2*R02
      pre_ugd = (sig0*R02)/(4.d0*(pi))
      

c     ------------------------------------------------------------------
      ugd = pre_ugd*dexp(-exp_ugd)
!      ugd = 0.5d0*f_kutak(imode,x,akt2,iread)
ctest      write(*,*) ugd,x,akt2
      return 
      end

      
c     ==================================================================
c     The aux functions
c     ==================================================================
      
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
      double precision k2,z2,eps2,p2,tau
      double precision e1_int_t1,e1_int_t2,e1_int_t3
c     ------------------------------------------------------------------
      common/e1var/k,p,z,eps
c     ------------------------------------------------------------------

      k2    = k**2.d0
      z2    = z**2.d0 
      eps2  = eps**2.d0
      p2    = p**2.d0
c     ------------------------------------------------------------------

      tau = p2 + z2*k2 - 2.d0*z*p*k*dcos(atheta)

c     ------------------------------------------------------------------

      e1_int_t1 = 1.d0/((tau+eps2)**2.d0)
      e1_int_t2 = 2.d0/((p2+eps2)*(tau+eps2))
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
      double precision k2,z2,eps2,p2,tau
      double precision e2_int_t1,e2_int_t2,e2_int_t3
c     ------------------------------------------------------------------
      common/e2var/k,p,z,eps
c     ------------------------------------------------------------------

      k2    = k**2.d0
      z2    = z**2.d0 
      eps2  = eps**2.d0
      p2    = p**2.d0

c     ------------------------------------------------------------------

      tau = p2 + z2*k2 - 2.d0*z*p*k*dcos(atheta)

c     ------------------------------------------------------------------

      e2_int_t1 = tau/((tau+eps2)**2.d0)
      e2_int_t2 = (p2-z*p*k*dcos(atheta))/((p2+eps2)*(tau+eps2))
      e2_int_t3 = p2/((p2+eps2)**2.d0)

c     ------------------------------------------------------------------

      e2_int = 0.5d0*(e2_int_t1 - 2.d0*e2_int_t2 + e2_int_t3)
      return 
      end

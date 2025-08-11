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
      

      hs = dsqrt(pt2+ (1.d0-xf)*M2) 
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

      

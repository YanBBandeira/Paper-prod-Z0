      program ugd_gbw_tmd
      implicit none 

      double precision kt,x,mu,ugd
      double precision kt_min,kt_max,dkt,log_kt_min,log_kt_max
      real f
      integer iset, ikt, nkt 

      external ugd
c     ------------------------------------------------------------------
      iset = 200001 !gbw light

      kt_min = 0.01d0
      kt_max = 50.d0
      nkt = 100 
      dkt = (kt_max - kt_min) / (nkt - 1)

      log_kt_min = dlog(kt_min)  
      log_kt_max = dlog(kt_max)

      call TMDinit(iset)
      call TMDset(iset)

c      open(unit=10,file='Output/ugd_gbw-x4-mu10.dat',
c     *     status='unknown')
c      open(unit=10, file='Output/ugd_gbw-x4-mu100.dat', 
c     *     status='unknown')
c      open(unit=10, file='Output/ugd_gbw-x4-mu300.dat', 
c     *     status='unknown')
c      open(unit=10, file='Output/ugd_gbw-x6-mu10.dat', 
c     *     status='unknown')
c      open(unit=10, file='Output/ugd_gbw-x6-mu100.dat', 
c     *     status='unknown')
c      open(unit=10, file='Output/ugd_gbw-x6-mu300.dat', 
c     *     status='unknown')
c      open(unit=10, file='Output/ugd_gbw-x2-mu100.dat', 
c     *     status='unknown')
      open(unit=10, file='Output/ugd_gbw_tmd.dat', 
     *     status='unknown')

c      x = 1.0d-4
c      x = 1.0d-6
      x = 1.0d-2
c      mu = 10.d0
c      mu = 100.d0
      mu = 300.d0

      do ikt = 1, nkt
            kt = dexp(log_kt_min + (ikt - 1) 
     &       * (log_kt_max - log_kt_min) / (nkt - 1))
            f  = ugd(kt,x,mu) 
            write(*,*) kt, f, x ,mu
            write(10,100) kt, f, x, mu
      end do 
100   format(2x,9(E10.4,2x))
      end



      FUNCTION ugd(akt,xbj,mu) !WW
      IMPLICIT none

      DOUBLE PRECISION ugd,akt,akt2,alphas,F_KS,xbj,x2,pre_ugd,
     &                 sc,pi,arg,mu

      EXTERNAL F_KS,sc  

      pi  = 4.d0*atan(1.d0)
      arg =  ( 4.d0 /akt**2.d0)  + 1.9d0 
      alphas = 0.3d0 !GBW fixed coupling
      akt2 = akt**2.d0

      x2 = xbj

      pre_ugd = 4.d0*pi*alphas/3.d0
      ugd = (pre_ugd*F_KS(x2,akt,mu))/akt2
ctest      write(*,*) ugd, pre_ugd, F_WW(xbj,akt), akt4
ctest      write(*,*) ugd, pre_ugd, x2, alphas
ctest      write(*,*) sc(M**2.d0)
      RETURN 
      END 

      FUNCTION F_KS(x2,akt,mmu)
      IMPLICIT none

      DOUBLE PRECISION F_KS,akt,x2 
      DOUBLE PRECISION x,xbar,kt,mu,mmu
      double precision up,ubar,dn,dbar,
     &                 strange,sbar,charm,cbar,bottom,bbar,glu

      Integer kf

      x = x2
      kt = akt
      xbar = 0.d0
      mu = mmu

      call TMDpdf(kf,x,xbar,kt,mu,up,ubar,dn,dbar,
     & strange,sbar,charm,cbar,bottom,bbar,glu)

      F_KS = glu
      return
      END

      FUNCTION sc(Q2) !strong coupling
      IMPLICIT none
      DOUBLE PRECISION sc, Q2, b0,lambqcd,lambqcd2,pi
      
      pi = 4.d0*atan(1.d0)
      lambqcd = 0.35d0 !GeV 
      lambqcd2 = lambqcd**2.d0
      b0 = 11.d0 - (10.d0/3.d0)

      sc = (4.d0*pi)/(b0*dlog(Q2/lambqcd2))
ctest      write(*,*) 'sc', sc, Q2, lambqcd2,pi
      return 
      end 

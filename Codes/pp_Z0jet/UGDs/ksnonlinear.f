      PROGRAM ksnonlinear_UGD
      implicit none 

      double precision f, x, kt, mu, ugd
      double precision kt_min, kt_max, dkt, log_kt_min 
      double precision log_kt_max, log_dkt
      integer nkt, ikt, iset

      external f

      iset = 400002 !kslinear

      kt_min = 0.01d0
      kt_max = 350.d0
      nkt = 200
      dkt = (kt_max - kt_min) / (nkt - 1)

      log_kt_min = dlog(kt_min)
      log_kt_max = dlog(kt_max)
      log_dkt = (log_kt_max - log_kt_min)/dble(nkt)

      call TMDinit(iset)
      call TMDset(iset)

      open(unit=10, file='ksnonlinear_UGD.dat', status='unknown')

      x = 1.0d-4

      mu = 10.d0 

      do ikt = 0, nkt

         kt = dexp(log_kt_min + dble(ikt)*log_dkt)

         ugd = f(x, kt, mu)

         write(10,100) kt, ugd
         write(*,100)  kt, ugd

      end do
100   format(2x,9(E10.4 , 2x ))     
      end 

c     --------------Unintegrated gluon distribution--------------------
      FUNCTION f(x, kt, mu)
      implicit none
      double precision f
      double precision x, kt, mu, tst, scLO, scNLO
      double precision kt2, alphas, pre_ugd, pi, x2
      double precision TMDgluon

      external TMDgluon, scLO, scNLO
      x2 = x
      pi = 4.d0*atan(1.d0)
      kt2 = kt**2.d0    
      tst = scLO(kt2)
      if(tst.gt.0.3d0) then
      alphas = tst
      else
      alphas = 0.3d0
      endif   
cold      alphas = 0.1d0!0.3d0
      pre_ugd = 4.d0*pi*alphas/3.d0
cOPT      f = (pre_ugd*TMDgluon(x2,kt, mu))/kt2
      f = TMDgluon(x2,kt, mu)/kt2
ctest      write(*,*) 'UGD', ugd, pre_ugd, TMDgluon(x2,kt,mu), kt2, alphas
      RETURN 
      END 


c     --------------Transverse gluon distribution--------------------
      FUNCTION TMDgluon(x2,akt, mu)
      IMPLICIT NONE
      DOUBLE PRECISION TMDgluon,akt,x2 , mu
      DOUBLE PRECISION x,xbar,kt
      double precision up,ubar,dn,dbar,
     &                 strange,sbar,charm,cbar,bottom,bbar,glu

cold     INTEGER iset
      Integer kf
      Integer TMDnumberPDF
      external TMDnumberPDF
      
      if(akt.lt.499.d0.or.akt.gt.0.011d0)then
      if(x2.gt.1.d-2.or.x2.lt.1.d-8) then
      TMDgluon = 0.d0 

c      write(*,*) 'not ok', x2, akt

      else 
      x = x2
      kt = akt
      xbar = 0.d0
      call TMDpdf(kf,x,xbar,kt,mu,up,ubar,dn,dbar,
     & strange,sbar,charm,cbar,bottom,bbar,glu)
      TMDgluon = glu
c      write(*,*) 'ok', x, kt
ctest      write(*,*) 'Transvere', glu, akt    
c      write(*,*) 'ok'
      end if
      else
      TMDgluon = 0.d0 
      end if
      return
      END

      function scLO(Q2) !strong coupling

      DOUBLE PRECISION scLO, Q2, b0,lambqcd,lambqcd2
      double precision t,  ca, nf, cf, pi

      pi = 4.d0*atan(1.d0)

      ca = 3.d0
      nf = 3.d0
      cf = 4.d0/3.d0

      lambqcd = 0.3d0 !GeV 
      lambqcd2 = lambqcd**2.d0

      b0 = (11.d0*ca - 2.d0*nf)/3

      t = dlog(Q2/lambqcd2)

      scLO = (4.d0*pi)/(b0*t)
ctest      write(*,*) 'sc', sc, Q2, lambqcd2,pi
      return 
      end 

      function scNLO(Q2) !NLO strong coupling
      double precision scNLO, Q2, b0,b1,lambqcd,lambqcd2
      double precision t,  ca, nf, cf, pi 

      pi = 4.d0*atan(1.d0)

      ca = 3.d0
      nf = 5.d0
      cf = 4.d0/3.d0

      lambqcd = 0.1d0 !GeV 
      lambqcd2 = lambqcd**2.d0

      b0 = (11.d0*ca - 2.d0*nf)/3
      b1 = (34.d0*ca*ca - 10.d0*ca*nf - 6.d0*cf*nf)/3.d0

      t = dlog(Q2/lambqcd2)

      scNLO = (4.d0*pi)/(b0*t)*(1.d0 - (b1/b0**2.d0)*dlog(t)/t)

      return 
      end
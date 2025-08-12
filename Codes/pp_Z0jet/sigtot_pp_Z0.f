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

      END PROGRAM  




      
      ! TODO - Aqui, será todo o integrando em termos da rapidez e momentos e angulos
      SUBROUTINE VegasIntegrand(SigTot) 
      IMPLICIT NONE

      COMMON/xbj/x2
      EXTERNAL DileptonDecay, HadronicCrossSection
      ! TODO - Aqui devo fornecer os dados de x2
      
      Result = DileptonDecay()*HadronicCrossSection()
      SigTot = Result
      END SUBROUTINE  

      FUNCTION DileptonDecay(M) 
      IMPLICIT NONE
      DOUBLE PRECISION DileptonDecay, M,M2,MZ2 
      DOUBLE PRECISION DecayWidth, Branch, InvariantMassDist, Result
      use parameters

      M2 = M*M
      Mz2 = Mz*Mz

      DecayWidth = ((alfem*M)/(6.d0*(dsin(2.d0*aw)**2.d0)))*(
     & (160.d0/3.d0)*(dsin(aw)**4.d0) - 40.d0*(dsin(aw)**2.d0) + 21.d0)

      Branch = 3.3662d0/100.d0
     
      InvariantMassDist = (1.d0/pi)*
     & ((M*DecayWidth)/((M2 - Mz2)**2.d0 + (M*DecayWidth)**2.d0))

      Result = InvariantMassDist*Branch
      DileptonDecay = Result

ctest      write(*,*)'Decay: ',Result,DecayWidth,Branch,InvariantMassDist
      RETURN 
      END 

      FUNCTION HadronicCrossSection()
      IMPLICIT NONE 
      ! TODO - Tem que ver como chama a LHAPDF
      

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
ctest      write(*,*) 'Masses', MU, MD, MS, MC, MB
ctest      write(*,*) 'Pdfs:' u,d,s,c,b,u_bar,d_bar,s_bar,c_bar,b_bar

      HadronicCrossSection = Result
      RETURN
      END 

      SUBROUTINE PartonTargetCrossSection(pt,z,mf,gfv,gfa)
            ! TODO - Aqui tem que passar todos os parâmetros necessários para baixo
      IMPLICIT NONE
      
      FUNCTION TransverseMomentumIntegral(N,X)
      IMPLICIT NONE 
      use parameters
      INTEGER N 
      DOUBLE PRECISION TransverseMomentumIntegral, X, preIntegrand
      DOUBLE PRECISION Int, akt, theta, kv
      DIMENSION X(2)
      EXTERNAL TransverseMomentumIntegrand

      !Integration variables 
      kv    = X(1)/(1.d0 - X(1)) !kv goes from zero to inf
      theta = X(2)*2.d0*pi       !theta goes from 0 to 2pi

      pi2 = pi*pi

      preIntegrand = DSQRT(alfem)/(2.d0*pi2*DSIN(2.d0*aW))

      Int = preIntegrand*TransverseMomentumIntegrand(kv,theta)
     &        *((1.d0 + kv)**2.d0)*(2.d0*pi) 
      ! both 2pi and (1+kv)^2 are jacobian factors
      
      TransverseMomentumIntegral = Int
      
ctest      write(*,*) 'int var ',x(2), x(1) 
      RETURN
      END  

      FUNCTION TransverseMomentumIntegrand(akt,atheta)
      IMPLICIT NONE 
      ! TODO - Falta definir o tipo das variáveis
      COMMON/PL_inf/gfv2,gfa2,mf,z,M2
      EXTERNAL ugd, Epsilon_1, Epsilon_2

      z2 = z*z 
      z4 = z2*z2
      mf = mf*mf 

      GammaL  =  gfv2*((1.d0-z)**2.d0)*M2 
     &        + gfa2*(((z2*mf2 + (1.d0-z)*M2)**2.d0)/M2)
      GammaT  = gfv2*z4*mf2 + gfa2*z2*mf2*((2.d0-z)**2.d0)
      LambdaL = gfa2*((z2*mf2)/M2)
      LambdaT = (1.d0 + ((1.d0 - z)**2.d0))*(gfv2 + gfa2)

      epps2 = (1.d0 - z)*M2 + z2*mf2
      epps  = DSQRT(epps2)
      Epsilon_1 = Epsilon_1(akt,atheta,pt,z,epps) 
      Epsilon_2 = Epsilon_1(akt,atheta,pt,z,epps)

      UGDF = ugd(akt)

      Result = akt*UGDF*((GammaT + 2.d0*GammaL)*Epsilon_1
     &        + (LambdaT + 2.d0*LambdaL)*Epsilon_2    )  
      TransverseMomentumIntegrand = Result

ctest write(*,*) Epsilon_1, Epsilon_2, UGDF, Result
      RETURN
      END        

      END SUBROUTINE

c     ==================================================================
c                       Epsilon's definitions
c     ==================================================================

      FUNCTION Epsilon_1(akt,atheta,pt,z,epps)
      IMPLICIT NONE 
      DOUBLE PRECISION Epsilon_1,akt,atheta,pt,z,epps
      ! TODO - Falta definir o tipo das variáveis
      akt2  = akt*akt
      pt2   = pt*pt
      z2    = z*z
      epps2 = epps*epps 
      
      tau2  = pt2 + z2*akt2 - 2.d0*z*pt*akt*dsin(atheta) 

      term1 = 1.d0/((tau2 + epps2)**2.d0 )
      term2 = 2.d0/((pt2+epps2)*(tau2 + epps2))
      term3 = 1.d0/((pt2 + epps2)**2.d0)
      Result = term1 - term2 + term3

      Epsilon_1 = 0.5*Result 
      
ctest      write(*,*) 'Epsilon_1', Epsilon_1, akt, atheta, pt, z, epps
      RETURN 
      END 

      FUNCTION Epsilon_2(akt,atheta,pt,z,epps)
      IMPLICIT NONE 
      DOUBLE PRECISION Epsilon_2,akt,atheta,pt,z,epps
      ! TODO - Falta definir o tipo das variáveis
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

      Epsilon_2 = 0.5*Result 
      
ctest      write(*,*) 'Epsilon_2', Epsilon_2, akt, atheta, pt, z, epps
      RETURN 
      END
c     ==================================================================
c                       Calling TMDlib
c     ==================================================================      


      FUNCTION ugd(akt) !WW
      use parameters

      DOUBLE PRECISION ugd,akt,akt2,alphas,F_KS,xbj,x2,pre_ugd,
     &                 M2,qv2,ma2,kv2,x1,kv,z,qv,sc,arg,
     &                 MG,MG2,ma,mb,ptq,alf,epps  

      COMMON/masses/MG,MG2,ma,mb
      COMMON/var/ptq,alf,epps
      COMMON/xx1/x1

      EXTERNAL F_KS,sc
      
      z = alf 
      qv = ptq


c      alphas = 0.12d0
      arg =  ( 4.d0 /akt**2.d0)  + 1.9d0 
      alphas = sc(arg)
c      alphas = 0.2d0
      akt2 = akt**2.d0

      M2 = MG2
      qv2 = qv**2.d0
      ma2 = ma**2.d0
      kv2 = akt2 
      ! TODO - checar a definição de x2!
      x2 = (((M2 + qv2)/z)+((ma2+((kv-qv)**2.d0))/(1.d0-z)))
     &     /(x1*(rs**2.d0))

      pre_ugd = 4.d0*pi*alphas/3.d0
      ugd = (pre_ugd*F_KS(x2,akt))/akt2
ctest      write(*,*) ugd, pre_ugd, F_WW(xbj,akt), akt4
ctest      write(*,*) ugd, pre_ugd, x2
      RETURN 
      END 

      FUNCTION F_KS(x2,akt)
      use parameters

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
c      write(*,*) iset, glu, ugd, akt    
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






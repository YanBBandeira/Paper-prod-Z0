CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC           		HANKEL FUNCTIONS INTEGRATION                 CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c	THE INTEGRATION IS MADE IN UU=PT*R 


	DOUBLE PRECISION FUNCTION Ii1(PTQ,Z, ETTA) 
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI1/ptt, EPS, alfa
	EXTERNAL I1int, DGAUSS
	PARAMETER(kbes=600)
	DIMENSION zb(kbes)
	PTT  = PTQ 
	PTT2 = PTT*PTT
	ALFA = Z 
	EPS  = ETTA
	!WRITE(*,*) PTT, ALFA, EPS
C	file with zeros of bessel function J0
	OPEN(unit=32,status='old',file='bzero.dat')
		DO k=1,kbes
			READ(32,*) zb(k)
		ENDDO
	CLOSE(32)
	li    = 0.d0
	soma1 = 0.d0
C	loop for the integral between zeros
	DO 10 j =1,kbes
		ls         = zb(j)
		!WRITE(*,*)li,ls
		step_I1int = DGAUSS(I1int,li,ls,1.D-3)!/PTT2
		soma1      = soma1 + step_I1int 
		li 	       = ls
		!WRITE(*,*) soma1, step_I1int
10    CONTINUE 
	Ii1 = soma1 
	!WRITE(*,*)'I1', Ii1 
	!WRITE(*,*) I1
	RETURN 
	END

	DOUBLE PRECISION FUNCTION Ii2(PTQ,Z, ETTA)
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI2/ptt, EPS, alfa
	EXTERNAL I2int, DGAUSS 
	PARAMETER(kbes=600)
	DIMENSION zb(kbes)
	PTT  = PTQ 
	PTT2 = PTT*PTT
	PTT3 = PTT2*PTT
	ALFA = Z 
	EPS  = ETTA
	!WRITE(*,*) PTT, ALFA, EPS
C	file with zeros of bessel function J0
	OPEN(unit=33,status='old',file='bzero.dat')
		DO k=1,kbes
			READ(33,*) zb(k)
		ENDDO
	CLOSE(33)
	li     = 0.D0 
	soma2  = 0.D0
C	loop for the integral between zeros
	DO 10 j=1,kbes
		ls = zb(j)
		!WRITE(*,*)
		step_I2int = DGAUSS(I2int,li,ls,1.D-3)!/PTT3
		soma2 = soma2 + step_I2int 
		li = ls
		!WRITE(*,*) soma2, setp_I2int
10    CONTINUE
	Ii2 = soma2 
	!WRITE(*,*)'I2', Ii2 
	RETURN 
	END

	DOUBLE PRECISION FUNCTION Ii3(PTQ,Z, ETTA)
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI3/PTT,EPS,ALFA
	EXTERNAL I3int, DGAUSS
	PARAMETER(kbes=600)
	DIMENSION zb(kbes)
	PTT   = PTQ 
	PTT2  = PTT*PTT
	ALFA  = Z 
	EPS   = ETTA 
	!WRITE(*,*) PTT, ALFA, EPS 
c     files with zeros of bessel functions J1
	OPEN(unit=42,status='old',file='b1zeros.dat')
		DO k=1,kbes
			READ(42,*) zb(k)
		ENDDO
	CLOSE(42) 	
	li    = 0.D0
	SOMA3 = 0.D0
C	loop for the integral between zeros
	DO 10 j=1,kbes
		ls = zb(j)
		!WRITE(*,*)li,ls
		step_I3int = DGAUSS(I3int,li,ls,1.D-3)!/PTT
		soma3  = soma3 + step_I3int 
		li = ls 
		!WRITE(*,*) soma3, step_I3int 
10    CONTINUE 
	Ii3 = SOMA3 
	!WRITE(*,*) I3,ptt,alfa,eps
	RETURN
	END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC           		HANKEL FUNCTIONS INTEGRANDS                  CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




	DOUBLE PRECISION FUNCTION I1int(uu) 
	IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI1/ptt, EPS, ALFA
	EXTERNAL bessj0,bessk0
	EXTERNAL bessj1,bessk1
	EXTERNAL sigdip
	u     = uu
	pt    = PTT 
	pt2   = pt*pt
	eps2  = eps*eps
	arg1  = u
	arg2  = (eps*u)/pt
	sdip  = sigdip(alfa*u/pt)
	I1int = u*bessj0(arg1)*bessk0(arg2)*sdip/pt2
	!WRITE(*,*)'I1int', bessj0(arg1),bessk1(arg2), sdip, u/pt, pt
	RETURN
	END 

	DOUBLE PRECISION FUNCTION I2int(uu)
        IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI2/ptt, EPS, ALFA
	EXTERNAL bessj0,bessk0
	EXTERNAL bessj1,bessk1
	EXTERNAL sigdip
	u     = uu
	u2    = u*u
	pt    = PTT 
	pt2   = pt*pt			
	pt3   = pt2*pt
	eps2  = eps*eps
	arg1  = u
	arg2  = (eps*u)/pt
	sdip  = sigdip(alfa*u/pt)
	I2int = u2*bessj0(arg1)*bessk1(arg2)*sdip/pt3
	!WRITE(*,*)'I2int', bessj0(arg1),bessk1(arg2), sdip, u/pt, pt
	RETURN 
	END

	DOUBLE PRECISION FUNCTION I3int(uu)
        IMPLICIT DOUBLE PRECISION(A-H,I,L-Z)
	COMMON/INTI3/PTT,EPS,ALFA
	EXTERNAL bessj0,bessk0
	EXTERNAL bessj1,bessk1
	EXTERNAL sigdip
	u     = uu
	pt    = PTT 
	pt2   = pt*pt
	eps2  = eps*eps
	arg1  = u
	arg2  = EPS*u/pt
	sdip  = sigdip(alfa*u/pt)
	I3int = u*bessj1(arg1)*bessk1(arg2)*sdip/pt2
	!WRITE(*,*)'I3int', bessj1(arg1),bessk1(arg2), sdip, u/pt, I3int
	RETURN
	END

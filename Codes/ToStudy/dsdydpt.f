    
      DOUBLE PRECISION function Ldsdydpt(y,pt)
c     -----------------------------------------------------------------
      double precision y,pt
c     -----------------------------------------------------------------
      integer ny,npt,NARG,NENT(2)
c     =================================================================
c     common blocks (for the grid)
c     -----------------------------------------------------------------
      double precision Lfy_grid(400,400), Ly_grid(400)
	double precision Lpt_grid(400), LENT(800)
	double precision Tfy_grid(400,400), Ty_grid(400)
	double precision Tpt_grid(400), TENT(800)
      double precision ARG(2)
      double precision DFINT  
c     -----------------------------------------------------------------
      common/gdL/LENT,Ly_grid,Lpt_grid,Lfy_grid  
c     -----------------------------------------------------------------
      external DFINT
ctest      write(*,*) Lfy_grid(10,10)
      ny=400
      npt=400
      NARG=2
      NENT(1) = ny
      NENT(2) = npt
c     -----------------------------------------------------------------
      ARG(1) = y
      ARG(2) = pt
c     -----------------------------------------------------------------
      Ldsdydpt = DFINT(NARG,ARG,NENT,LENT,Lfy_grid)
      !write(*,*)       Ldsdydpt
	return
	end
	
    	DOUBLE PRECISION function Tdsdydpt(y,pt)
c     -----------------------------------------------------------------
      double precision y,pt
c     -----------------------------------------------------------------
      integer ny,npt,NARG,NENT(2)
c     =================================================================
c     common blocks (for the grid)
c     -----------------------------------------------------------------
      double precision Lfy_grid(400,400), Ly_grid(400)
	double precision Lpt_grid(400), LENT(800)
	double precision Tfy_grid(400,400), Ty_grid(400)
	double precision Tpt_grid(400), TENT(800)
      double precision ARG(2)
      double precision DFINT  
c     -----------------------------------------------------------------
      common/gdT/TENT,Ty_grid,Tpt_grid,Tfy_grid 
c     -----------------------------------------------------------------
      external DFINT
      ny=400
      npt=400
      NARG=2

      NENT(1) = ny
      NENT(2) = npt
c     -----------------------------------------------------------------
      ARG(1) = y
      ARG(2) = pt
c     -----------------------------------------------------------------
      Tdsdydpt = DFINT(NARG,ARG,NENT,TENT,Tfy_grid)
      !write(*,*)       Tdsdydpt
	return
	end

      

                

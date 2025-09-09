      function InterpolateGrid(y,pt)
      integer::parameter:: nPoints = 5000
      integer::parameter:: narg=2
      dimension nent(2), arg(narg)
      dimension ent(223)
      common/GridArrays/ent,yGrid,ptGrid,PartonLevelGrid

      double precision PartonLevelGrid(nPoints,nPoints)
      double precision ptGrid(nPoints), yGrid(nPoints)
      double precision ent(nPoints + nPoints)
      char*200 :: File
      
      File = 'kslinear_grid.dat'

      call read_grid(nPoints,File)
      
      nent(1) = nPoints
      nent(2) = nPoints

      arg(1) = y
      arg(2) = pt

      InterpolateGrid = DFINT(narg,arg,nent,ent,PartonLevelGrid)
      return
      end 


      subroutine read_grid(nPoints,OutputPath)
      integer::parameter :: iFile = 11
      integer:: nPoints
      double precision PartonLevelGrid(nPoints,nPoints)
      double precision ptGrid(nPoints), yGrid(nPoints)
      double precision ent(nPoints + nPoints)
      char*200 :: OutputPath
      
      common/GridArrays/ent,yGrid,ptGrid,PartonLevelGrid

      open(iFile,file=OutputPath,status='old')

      read(iFile,*)
      do i = 1, nPoints
            do j = 1, nPoints
      read(iFile,*) yVar, ptVar, PartonLevelVar
      yGrid(i) = yVar
      ptGrid(j) = ptVar
      PartonLevelGrid(i,j) = PartonLevelVar
            end do
      end do

      do i = 1, nPoints
         ent(i) = yGrid(i)
      end do 

      do j = 1, nPoints
         ent(nPoints + j) = ptGrid(j)
      end do

      end 
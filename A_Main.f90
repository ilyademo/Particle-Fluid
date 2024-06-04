Program Main
  use Functions
  implicit none
  include 'omp_lib.h'
  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt', SolutionFile='solution_for_particle_dimension.plt' ! names of input and output files
  character MeshFile*30        ! name of file with computational mesh
  character GradType*30, ctmp
  character gradT_calc*30
  integer NI,NJ,i,j, Niter_GG, mode, iter, Niter, scheme, ip, jp, ip1, jp1, m, St
  integer, parameter:: IO = 12, I5=16 ! input-output unit
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume, T, Res, T_flos, Vmod! scalar arrays
  real,allocatable,dimension(:,:):: divV, divVExact, divVError, LapP, LapPExact, LapPError, rotV, rotVExact, rotVError
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,GradP, GradPExact, GradPError, V, Vel_eq, GradT ! vector arrays
  real rtmp, RO1, RO2, MU, dp, x0, y0, v0, u0, w0, nt, vref, lref, sk, x_m, y_m, u_m, v_m, w_m, fx, fy, mz, xm, ym, x_m1, y_m1, u_m1, v_m1, w_m1, dt, Vx_p, Vy_p, wg
  real u_old, v_old
!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  !read(IO,*) GradType
  !read(IO,*) mode
  NIter_GG = 5
  CLOSE(IO)
!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(Vmod(0:NI,0:NJ))
  allocate(GradP(0:NI, 0:NJ, 2)) !Gradient of pressure
  allocate(GradPExact(0:NI, 0:NJ, 2))
  allocate(GradPError(0:NI, 0:NJ, 2))
  allocate(V(0:NI,0:NJ,2))
  allocate(Vel_eq(0:NI,0:NJ,2))
  allocate(divV(0:NI,0:NJ))
  allocate(divVExact(0:NI,0:NJ))
  allocate(divVError(0:NI,0:NJ))
  allocate(rotV(0:NI,0:NJ))
  allocate(rotVExact(0:NI,0:NJ))
  allocate(rotVError(0:NI,0:NJ))
  allocate(LapP(0:NI,0:NJ))
  allocate(LapPExact(0:NI,0:NJ))
  allocate(LapPError(0:NI,0:NJ))
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces  
  
!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J), I=1,NI),J=1,NJ)
  CLOSE(IO)

!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
      P(I,J) = Pressure(CellCenter(i,j,1),CellCenter(i,j,2))
      Call CalcGradExact(CellCenter(i,j,1), CellCenter(i,j,2), GradPExact(i,j,:))
      Call Velocity(CellCenter(i,j,1), CellCenter(i,j,2), V(i,j,:))
      divVExact(i,j) = divVelocityPExact(CellCenter(i,j,1), CellCenter(i,j,2))
      LapPExact(i,j) = LapPressureExact(CellCenter(i,j,1), CellCenter(i,j,2))
      rotVExact(i,j) = RotVelExact(CellCenter(i,j,1), CellCenter(i,j,2))
    ENDDO
  ENDDO
!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate derivatives'
  if (GradType .EQ. "GG") then
  GradP = 0.0
  !do i = 1, NIter_GG
  Call B_CalcGradient(GradP, NI, NJ, P, CellVolume, CellCenter, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  GradPError = Abs(GradPExact-GradP)/GradPExact
  write(*,*) 'Iteration: ', i, 'maximum error ', maxval(GradPError(1:NI-1, 1:NJ-1, :))
  
  !end do
  
  else if (GradType .EQ. "LS") then
  GradP = 0.0
  Call LeastSquareGradCalc(GradP, NI, NJ, P, CellCenter)
  GradPError = Abs(GradPExact-GradP)/GradPExact
  write(*,*) 'Maximum error LeastSquareCalc:', maxval(GradPError(1:NI-1, 1:NJ-1, :))
  end if
  write(*,*) 'Maximum Gradx-error', maxval(GradPError(1:NI-1, 1:NJ-1, 1))
  write(*,*) 'Maximum Grady-error', maxval(GradPError(1:NI-1, 1:NJ-1, 2))

  !=== CALCULATE DIVERGENCE ===
  write(*,*) "Calculate Divergence"
  divV = 0.0
  call B_CalcDivergence(mode, NI, NJ, V, divV, P, GradP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  divVError = abs(divVExact - divV)/divVExact
  write (*,*) 'Maximum error of divV', maxval(divVError(1:NI-1, 1:NJ-1))
  
  !=== CALCULATE LAPLACIAN ====
  write (*,*) "Calculate Laplacian"
  LapP = 0.0
  Call B_CalcLaplacian(NI, NJ, P, GradP, LapP, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  LapPError = abs(LapPExact - LapP)/LapPExact
  write (*,*) 'Maximum error of Laplacian', maxval(LapPError(1:NI-1, 1:NJ-1))
  
  
  WRITE(*,*) 'Read parameters for particles'
  OPEN(34, FILE = 'input_particle.txt')
  READ(34,*) RO1, RO2
  READ(34,*) MU
  READ(34,*) dp
  READ(34,*) X0, Y0
  READ(34,*) U0, V0
  READ(34,*) W0
  READ(34,*) dt
  READ(34,*) NT
  CLOSE(34)
  
  OPEN(IO, FILE = SOLUTIONFILE)
  READ(IO,*) CTMP
  READ(IO,*) CTMP
  WRITE(*,*) 'Read solution data from file: ', SOLUTIONFILE
  READ(IO,*) ((RTMP,RTMP,V(I,J,1),V(I,J,2),Vmod(I,J),P(I,J),RTMP,RTMP,I=0,NI),J=0,NJ)
  CLOSE(IO)
  
  !=== CALCULATE ROTOR ====
  write (*,*) "Calculate Rotor"
  rotV = 0.0
  Call BCalcRotor(NI, NJ, V, rotV, CellCenter, CellVolume, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  
  Lref = 2e-2
  Vref = 1.82e-5*100/(Lref*ro1)
  
  
  Sk=RO2*dp**2*Vref/(18*MU*Lref)
  
  WRITE(*,*) 'Particle: ', 'd = ',dp, 'ro = ',ro2
  WRITE(*,*) 'Stokes number Sk = ', Sk
  WRITE(*,*) 'Iteration process: ', 'NT = ',NT, 'dt = ',dt
  
  pause
  x_m=x0
  y_m=y0
  u_m=u0
  v_m=v0
  w_m=w0
  u_old=u0
  v_old=v0
  OPEN(35,FILE = 'OutputFile_particle.plt', status = 'replace')
  WRITE(35,*) 'it,t,x,y,u,v,w,Fx,Fy,Mz'
  OPEN(I5,FILE = 'Particle_Force.plt', status = 'replace')
  WRITE(I5,*) 'it,t,F_D,F_VM,F_Saf'
  
  OPEN(IO,FILE = 'Rotor.plt', status = 'replace')
  call B_OutputFields(IO,NI,NJ,X,Y, rotV)
  close (IO)
  DO m=1,NT
      IP = -1
      JP = -1
      call C_Location(x_m, y_m, NI, NJ, x, y, CellVolume, ip, jp)
      
      Vx_p=V(IP,JP,1)
      Vy_p=V(IP,JP,2)
      wg = 0.5*rotV(IP,JP)
      Fx = 0.0
      Fy = 0.0
      Mz = 0.0
      call C_Force(RO1,RO2,MU,dp,u_m,v_m,w_m,Vx_p,Vy_p,wg,Fx,Fy,Mz,I5,m,dt, u_old, v_old)
      WRITE (35, "(i5, 9(1x,f15.10))") m, m*dt, x_m, y_m, u_m, v_m, w_m, Fx, Fy, Mz
      
      x_m1 = x_m+dt*u_m
      y_m1 = y_m+dt*v_m
      u_m1 = u_m+dt*Fx
      v_m1 = v_m+dt*Fy
      w_m1 = w_m+dt*Mz
      call C_Location(x_m1, y_m1, NI, NJ, x, y, CellVolume, ip1, jp1)
      
      if (ip1 .eq. -1 .and. jp1 .eq. -1) then
          call C_Boundary(x, y, iFaceVector, jFaceVector, NI, NJ, x_m, y_m, u_m, v_m, w_m, x_m1, y_m1, u_m1, v_m1, w_m1, ip1, jp1, St, dp)
          if ((St.eq.3).or.(St.eq.4)) then
              write (*,*) "Particle left channel through boundary ", St
              exit
          end if
      end if
      x_m = x_m1
      y_m = y_m1
      u_old = u_m
      v_old = v_m
      u_m = u_m1
      v_m = v_m1
      w_m = w_m1
  enddo
  CLOSE(35)
END PROGRAM Main  


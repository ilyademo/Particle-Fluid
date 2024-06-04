Subroutine B_OutputFields(IO,NI,NJ,X,Y, rotV)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P, divV, divVError, LapP, LapPError, rotV, rotVError, T, T_error
  Real,Dimension(0:NI,0:NJ,2)::GradP, GradPError, V, Vel_eq

  Write(IO, '(a190)') 'VARIABLES = "X", "Y", "RotV"' 
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F25.7)') rotV(1:NI-1,1:NJ-1)
End Subroutine 

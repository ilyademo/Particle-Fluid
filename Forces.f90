Subroutine C_Force(RO1,RO2,MU,dp,u_m,v_m,w_m,Vx_p,Vy_p,wg,Fx,Fy,Mz,I5,m, dt, u_old, v_old)
    implicit none
    INTEGER:: I5, m
    REAL :: RO1,RO2,MU,dp,u_m,v_m,w_m,Vx_p,Vy_p,Fx,Fy,Mz,wg,dt,C_m,I, u_old, v_old
    REAL :: CD, Vr, FD(2),FG(2), Rew,Rep,PI,f,mp, FS(2), Fvm(2)
    
    PI = 3.141592
    mp = 4.0/3.0*PI*(0.5*dp)**3*ro2
    
    Vr = sqrt((u_m-Vx_p)**2+(v_m-Vy_p)**2)
    Rep = dp*ro1*Vr/mu
    Rew = ro1*(dp**2)*(wg-w_m)/mu
    C_m = 64.0/Rew
    I = pi*(dp**5.0)*ro2/60.0
    !Stokes
    f = 1.0 + 0.179*sqrt(Rep) +0.013*Rep
    FD(1) = f*mu*18/(dp**2)/ro2*(Vx_p-u_m)
    FD(2) = f*mu*18/(dp**2)/ro2*(Vy_p-v_m)
    !Saffman
    FS(1) = ((1.414*6.46)/4.0)*dp**2*(mu*ro1)**(0.5)*(abs(wg))**(-0.5)*wg*(Vy_p-v_m)/mp
    FS(2) = ((1.414*6.46)/4.0)*dp**2*(mu*ro1)**(0.5)*(abs(wg))**(-0.5)*wg*(u_m-Vx_p)/mp
    
    !VirtualMass
    Fvm(1) = -(2./3.*PI*(dp/2)**3*ro1)*(u_m-u_old)/dt/mp
    Fvm(2) = -(2./3.*PI*(dp/2)**3*ro1)*(v_m-v_old)/dt/mp
    !print *, sqrt(Fvm(1)**2 + Fvm(2)**2)
    write(I5,"(i5, 9(1x,f25.10))") m, m*dt, sqrt((FD(1))**2+(FD(2))**2), sqrt((Fvm(1))**2+(Fvm(2))**2), sqrt((FS(1))**2+(FS(2))**2)
    !Gravity
    FG(1)=0.0
    FG(2)=0.0 !9.8
    
    Fx = FD(1)+Fvm(1)+FS(1)
    Fy = FD(2)+Fvm(2)+FS(2)
    Mz = 60.0*mu*(wg-w_m)/ro2/dp/dp
    
   
    End Subroutine
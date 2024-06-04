subroutine C_Boundary(x, y, iFaceVector, jFaceVector, NI, NJ, x_m, y_m, u_m, v_m, w_m, x_m1, y_m1, u_m1, v_m1, w_m1, ip1, jp1, St, dp)
    implicit none
    integer Ni, Nj, ip1, jp1, St, i, j, icr, mode
    real x(ni, nj), y(ni, nj)
    real jFaceVector(ni-1, nj, 2), iFaceVector(ni, nj-1, 2)
    real x_m1, y_m1, u_m1, v_m1, w_m1, x_m, y_m, u_m, v_m, w_m, x_m1_new, y_m1_new
    real xcros, ycros
    real D, NV(2), TV(2), Vn, Vtau, dp
    St = 0
    mode = 2
    ! j = 1 - wall St = 1
    ! j = nj - wall St = 2
    ! i = 1 - inlet St = 3
    ! i = ni - out St = 4
    
    do i = 1, ni - 1
        call Cross_edges(x_m, y_m, x_m1, y_m1, x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1), xcros, ycros, icr)
        if (icr.eq.1) then
            St = 1
            call Reflect(x_m, y_m, x_m1, y_m1, x(i,1), y(i,1), x(i+1,1), y(i+1, 1), x_m1_new, y_m1_new)
            write(*,*) "Particle reflect boundary ", St, " in coordinates ", xcros, ycros, " to ", x_m1_new, y_m1_new
            x_m1 = x_m1_new
            y_m1 = y_m1_new
            
            D = norm2(jFaceVector(i,1,:))
            
            NV(1) = -JFaceVector(i, 1, 1)/D
            NV(2) = -jFaceVector(i, 1, 2)/D
            TV(1) = -NV(2)
            TV(2) = NV(1)
            
            Vn = u_m*NV(1) + v_m*NV(2)
            Vtau = u_m*TV(1) + v_m*TV(2)
            call GetNewVeloicitiesAfterReflect(mode, u_m1, v_m1, w_m1, Vtau, Vn, TV, NV, u_m, v_m, w_m, dp)
            return
        endif
        
        call Cross_edges(x_m, y_m, x_m1, y_m1, x(i, NJ), y(i, NJ), x(i+1, NJ), y(i+1, NJ), xcros, ycros, icr)
        if (icr.eq.1) then
            St = 2
            call Reflect(x_m, y_m, x_m1, y_m1, x(i,NJ), y(i,NJ), x(i+1,NJ), y(i+1, NJ), x_m1_new, y_m1_new)
            write(*,*) 'Particle reflect boundary ', St, ' in coordinates ', xcros, ycros, " to ", x_m1_new, y_m1_new
            x_m1 = x_m1_new
            y_m1 = y_m1_new
            D = norm2(jFaceVector(i,NJ,:))
            
            NV(1) = -JFaceVector(i, NJ, 1)/D
            NV(2) = -jFaceVector(i, NJ, 2)/D
            TV(1) = -NV(2)
            TV(2) = NV(1)
            
            Vn = u_m*NV(1) + v_m*NV(2)
            Vtau = u_m*TV(1) + v_m*TV(2)
            call GetNewVeloicitiesAfterReflect(mode, u_m1, v_m1, w_m1, Vtau, Vn, TV, NV, u_m, v_m, w_m, dp)
            return
        endif
    end do
    
    do j = 1, nj -1
        call Cross_edges(x_m, y_m, x_m1, y_m1, x(1, j), y(1, j), x(1, j+1), y(1, j+1), xcros, ycros, icr)
        if (icr.eq.1) then
            St = 3
            call Reflect(x_m, y_m, x_m1, y_m1, x(1,j), y(1,j), x(1,j+1), y(1, j+1), x_m1_new, y_m1_new)
            write(*,*) 'Particle cross boundary', St, 'in coordinates', xcros, ycros
            return
        endif
        
        call Cross_edges(x_m, y_m, x_m1, y_m1, x(NI, j), y(NI, j), x(NI, j+1), y(NI, j+1), xcros, ycros, icr)
        if (icr.eq.1) then
            St = 4
            call Reflect(x_m, y_m, x_m1, y_m1, x(NI,j), y(NI,j), x(NI,j+1), y(NI, j+1), x_m1_new, y_m1_new)
            write(*,*) "Particle cross boundary ", St, " in coordinates ", xcros, ycros
            return
        endif
    end do
    
    end subroutine
    
    subroutine Reflect(x1,y1,x2,y2,x3,y3,x4,y4,x2_new,y2_new)
    
    REAL :: x1,y1,x2,y2,x3,y3,x4,y4
    REAL :: x2_new,y2_new
    REAL :: a,b,c,denom
    
    a = y4-y3
    b = x3-x4
    c = -(a*x4+b*y4)
    denom = a*a+b*b
    
    x2_new = ((b**2-a**2)*x2-2.0*a*b*y2-2.0*c*a)/denom
    y2_new = (-2.0*a*b*x2+(a**2-b**2)*y2-2.0*c*b)/denom
    
    end subroutine
    
    subroutine GetNewVeloicitiesAfterReflect(mode, u_m1, v_m1, w_m1, Vtau, Vn, TV, NV, u_m, v_m, w_m, dp)
    implicit none
    integer:: mode
    real u_m1, v_m1, w_m1, Vtau, Vn, TV(2), NV(2), u_m, v_m, w_m, dp
    if (mode == 1) then
        u_m1 = Vtau*TV(1) - Vn*NV(1)
        v_m1 = Vtau*TV(2) - Vn*NV(2)
        w_m1 = w_m
    else if (mode == 2) then
        u_m1 = u_m - ((1+0.7)*(u_m*NV(1)+v_m*NV(2))*NV(1) + (2.0/7.0)*TV(1)*ABS(Vtau))
        v_m1 = v_m - ((1+0.7)*(u_m*NV(1)+v_m*NV(2))*NV(2) + (2.0/7.0)*TV(2)*ABS(Vtau))
        w_m1 = w_m - (5.0/(7.0*(dp/2.0)))*ABS(Vtau)*(NV(1)*TV(2)-NV(2)*TV(1))
    endif
    end subroutine
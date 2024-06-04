subroutine C_Location(x_n1, y_n1, NI, NJ, X, Y, CellVolume, ip, jp)
    integer NI, NJ, ip, jp, i, j
    real x_n1, y_n1, x(ni, nj), y(ni, nj), CellVolume(NI-1, NJ-1), S, eps
    ip = -1
    jp = -1
    eps = 1e-20
    do j=1,nj-1
        do i=1,ni-1
            s = TS(x(i,j), y(i,j), x(i, j+1), y(i,j+1), x_n1, y_n1) + TS(x(i+1, j+1), y(i+1,j+1), x(i, j+1), y(i,j+1), x_n1, y_n1) + &
                TS(x(i,j), y(i,j), x(i+1,j), y(i+1,j), x_n1, y_n1) + TS(x(i+1, j+1), y(i+1, j+1), x(i+1, j), y(i+1,j), x_n1, y_n1)
        if (abs(s-CellVolume(i,j)).LE.eps) then
            ip = i
            jp = j
            exit
        endif
        
        end do
    end do
    
end subroutine
    
real Function TS(x1, y1, x2, y2, x3, y3) !Triangular square
real x1, y1, x2, y2, x3, y3, a, b, c, p
a = sqrt((x1 - x2)**2 + (y1-y2)**2)
b = sqrt((x3 - x2)**2 + (y3 - y2)**2)
c = sqrt((x1 - x3)**2 + (y1 - y3)**2)
p = (a+b+c)/2
TS = sqrt((p - a) * (p - b) * (p - c) * p)
end function



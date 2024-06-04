subroutine Cross_edges(x1, y1, x2, y2, x3, y3, x4, y4, xcros, ycros, icr)
    implicit none
    real x1, x2, x3, x4, y1, y2, y3, y4
    real xcros, ycros
    real p(2), q(2), r(2), s(2), t, u, denom, num_u, num_t
    integer icr
    real:: eps = 1e-16
    xcros = -1e6
    ycros = -1e6
    p(1) = x1
    p(2) = y1
    q(1) = x3
    q(2) = y3
    
    r(1) = x2-x1
    r(2) = y2-y1
    s(1) = x4-x3
    s(2) = y4-y3
    
    icr = 0 ! no intesction
    
    denom = r(1)*s(2) - r(2)*s(1)
    num_t = ((q(1) - p(1))*s(2) - (q(2)-p(2))*s(1))
    num_u = -((p(1)-q(1))*r(2) - (p(2) - q(2))*r(1))
    
    if (abs(denom).lt.eps) then
        if (abs(num_u).gt.eps) then
            return ! parallel lines
        else
            icr = 1
            xcros = x2
            ycros = y2
            return
        end if
    end if
    
    t = num_t/denom
    u = num_u/denom
    
    if ((t.gt.0.0.and.t.lt.1.0).and.(u.gt.0.0.and.u.lt.1.0)) then
        icr = 1
        xcros = p(1) + t*r(1)
        ycros = p(2) + t*r(2)
    end if
end subroutine
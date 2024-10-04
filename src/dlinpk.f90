subroutine dgefa(a, lda, n, ipvt, info)
    integer :: lda, n, ipvt(*), info
    double precision :: a(lda, *)
    double precision :: t
    integer :: idamax, j, k, kp1, l, nm1
    info = 0
    nm1 = n - 1
    if (nm1 .lt. 1) then
        go to 70
    end if
    do k = 1, nm1
        kp1 = k + 1
        l = idamax(n - k + 1, a(k, k), 1) + k - 1
        ipvt(k) = l
        if (a(l, k) .eq. 0.0d0) then
            go to 40
        end if
        if (l .eq. k) then
            go to 10
        end if
        t = a(l, k)
        a(l, k) = a(k, k)
        a(k, k) = t
10      continue
        t = (-1.0d0)/a(k, k)
        call dscal(n - k, t, a(k + 1, k), 1)
        do j = kp1, n
            t = a(l, j)
            if (l .eq. k) then
                go to 20
            end if
            a(l, j) = a(k, j)
            a(k, j) = t
20          continue
            call daxpy(n - k, t, a(k + 1, k), 1, a(k + 1, j), 1)
30          continue
        end do
        go to 50
40      continue
        info = k
50      continue
60      continue
    end do
70  continue
    ipvt(n) = n
    if (a(n, n) .eq. 0.0d0) then
        info = n
    end if
    return
end subroutine dgefa

subroutine dgesl(a, lda, n, ipvt, b, job)
    integer :: lda, n, ipvt(*), job
    double precision :: a(lda, *), b(*)
    double precision :: ddot, t
    integer :: k, kb, l, nm1
    nm1 = n - 1
    if (job .ne. 0) then
        go to 50
    end if
    if (nm1 .lt. 1) then
        go to 30
    end if
    do k = 1, nm1
        l = ipvt(k)
        t = b(l)
        if (l .eq. k) then
            go to 10
        end if
        b(l) = b(k)
        b(k) = t
10      continue
        call daxpy(n - k, t, a(k + 1, k), 1, b(k + 1), 1)
20      continue
    end do
30  continue
    do kb = 1, n
        k = n + 1 - kb
        b(k) = b(k)/a(k, k)
        t = -b(k)
        call daxpy(k - 1, t, a(1, k), 1, b(1), 1)
40      continue
    end do
    go to 100
50  continue
    do k = 1, n
        t = ddot(k - 1, a(1, k), 1, b(1), 1)
        b(k) = (b(k) - t)/a(k, k)
60      continue
    end do
    if (nm1 .lt. 1) then
        go to 90
    end if
    do kb = 1, nm1
        k = n - kb
        b(k) = b(k) + ddot(n - k, a(k + 1, k), 1, b(k + 1), 1)
        l = ipvt(k)
        if (l .eq. k) then
            go to 70
        end if
        t = b(l)
        b(l) = b(k)
        b(k) = t
70      continue
80      continue
    end do
90  continue
100 continue
    return
end subroutine dgesl

subroutine dgbfa(abd, lda, n, ml, mu, ipvt, info)
    integer :: lda, n, ml, mu, ipvt(*), info
    double precision :: abd(lda, *)
    double precision :: t
    integer :: i, idamax, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm, nm1
    m = ml + mu + 1
    info = 0
    j0 = mu + 2
    j1 = min0(n, m) - 1
    if (j1 .lt. j0) then
        go to 30
    end if
    do jz = j0, j1
        i0 = m + 1 - jz
        do i = i0, ml
            abd(i, jz) = 0.0d0
10          continue
        end do
20      continue
    end do
30  continue
    jz = j1
    ju = 0
    nm1 = n - 1
    if (nm1 .lt. 1) then
        go to 130
    end if
    do k = 1, nm1
        kp1 = k + 1
        jz = jz + 1
        if (jz .gt. n) then
            go to 50
        end if
        if (ml .lt. 1) then
            go to 50
        end if
        do i = 1, ml
            abd(i, jz) = 0.0d0
40          continue
        end do
50      continue
        lm = min0(ml, n - k)
        l = idamax(lm + 1, abd(m, k), 1) + m - 1
        ipvt(k) = l + k - m
        if (abd(l, k) .eq. 0.0d0) then
            go to 100
        end if
        if (l .eq. m) then
            go to 60
        end if
        t = abd(l, k)
        abd(l, k) = abd(m, k)
        abd(m, k) = t
60      continue
        t = (-1.0d0)/abd(m, k)
        call dscal(lm, t, abd(m + 1, k), 1)
        ju = min0(max0(ju, mu + ipvt(k)), n)
        mm = m
        if (ju .lt. kp1) then
            go to 90
        end if
        do j = kp1, ju
            l = l - 1
            mm = mm - 1
            t = abd(l, j)
            if (l .eq. mm) then
                go to 70
            end if
            abd(l, j) = abd(mm, j)
            abd(mm, j) = t
70          continue
            call daxpy(lm, t, abd(m + 1, k), 1, abd(mm + 1, j), 1)
80          continue
        end do
90      continue
        go to 110
100     continue
        info = k
110     continue
120     continue
    end do
130 continue
    ipvt(n) = n
    if (abd(m, n) .eq. 0.0d0) then
        info = n
    end if
    return
end subroutine dgbfa

subroutine dgbsl(abd, lda, n, ml, mu, ipvt, b, job)
    integer :: lda, n, ml, mu, ipvt(*), job
    double precision :: abd(lda, *), b(*)
    double precision :: ddot, t
    integer :: k, kb, l, la, lb, lm, m, nm1
    m = mu + ml + 1
    nm1 = n - 1
    if (job .ne. 0) then
        go to 50
    end if
    if (ml .eq. 0) then
        go to 30
    end if
    if (nm1 .lt. 1) then
        go to 30
    end if
    do k = 1, nm1
        lm = min0(ml, n - k)
        l = ipvt(k)
        t = b(l)
        if (l .eq. k) then
            go to 10
        end if
        b(l) = b(k)
        b(k) = t
10      continue
        call daxpy(lm, t, abd(m + 1, k), 1, b(k + 1), 1)
20      continue
    end do
30  continue
    do kb = 1, n
        k = n + 1 - kb
        b(k) = b(k)/abd(m, k)
        lm = min0(k, m) - 1
        la = m - lm
        lb = k - lm
        t = -b(k)
        call daxpy(lm, t, abd(la, k), 1, b(lb), 1)
40      continue
    end do
    go to 100
50  continue
    do k = 1, n
        lm = min0(k, m) - 1
        la = m - lm
        lb = k - lm
        t = ddot(lm, abd(la, k), 1, b(lb), 1)
        b(k) = (b(k) - t)/abd(m, k)
60      continue
    end do
    if (ml .eq. 0) then
        go to 90
    end if
    if (nm1 .lt. 1) then
        go to 90
    end if
    do kb = 1, nm1
        k = n - kb
        lm = min0(ml, n - k)
        b(k) = b(k) + ddot(lm, abd(m + 1, k), 1, b(k + 1), 1)
        l = ipvt(k)
        if (l .eq. k) then
            go to 70
        end if
        t = b(l)
        b(l) = b(k)
        b(k) = t
70      continue
80      continue
    end do
90  continue
100 continue
    return
end subroutine dgbsl

subroutine daxpy(n, da, dx, incx, dy, incy)
    double precision :: dx(*), dy(*), da
    integer :: i, incx, incy, ix, iy, m, mp1, n
    if (n .le. 0) then
        return
    end if
    if (da .eq. 0.0d0) then
        return
    end if
    if (incx .eq. 1 .and. incy .eq. 1) then
        go to 20
    end if
    ix = 1
    iy = 1
    if (incx .lt. 0) then
        ix = ((-n) + 1)*incx + 1
    end if
    if (incy .lt. 0) then
        iy = ((-n) + 1)*incy + 1
    end if
    do i = 1, n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
10      continue
    end do
    return
20  m = mod(n, 4)
    if (m .eq. 0) then
        go to 40
    end if
    do i = 1, m
        dy(i) = dy(i) + da*dx(i)
30      continue
    end do
    if (n .lt. 4) then
        return
    end if
40  mp1 = m + 1
    do i = mp1, n, 4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
50      continue
    end do
    return
end subroutine daxpy

subroutine dcopy(n, sx, incx, sy, incy)
    double precision :: sx(*), sy(*)
    integer :: i, incx, incy, ix, iy, m, mp1, n
    if (n .le. 0) then
        return
    end if
    if (incx .eq. 1 .and. incy .eq. 1) then
        go to 20
    end if
    ix = 1
    iy = 1
    if (incx .lt. 0) then
        ix = ((-n) + 1)*incx + 1
    end if
    if (incy .lt. 0) then
        iy = ((-n) + 1)*incy + 1
    end if
    do i = 1, n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
10      continue
    end do
    return
20  m = mod(n, 7)
    if (m .eq. 0) then
        go to 40
    end if
    do i = 1, m
        sy(i) = sx(i)
30      continue
    end do
    if (n .lt. 7) then
        return
    end if
40  mp1 = m + 1
    do i = mp1, n, 7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
50      continue
    end do
    return
end subroutine dcopy

subroutine dscal(n, da, dx, incx)
    double precision :: da, dx(*)
    integer :: i, incx, m, mp1, n, nincx
    if (n .le. 0) then
        return
    end if
    if (incx .eq. 1) then
        go to 20
    end if
    nincx = n*incx
    do i = 1, nincx, incx
        dx(i) = da*dx(i)
10      continue
    end do
    return
20  m = mod(n, 5)
    if (m .eq. 0) then
        go to 40
    end if
    do i = 1, m
        dx(i) = da*dx(i)
30      continue
    end do
    if (n .lt. 5) then
        return
    end if
40  mp1 = m + 1
    do i = mp1, n, 5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
50      continue
    end do
    return
end subroutine dscal

double precision function ddot(n, dx, incx, dy, incy)
    double precision :: dx(*), dy(*), dtemp
    integer :: i, incx, incy, ix, iy, m, mp1, n
    ddot = 0.0d0
    dtemp = 0.0d0
    if (n .le. 0) then
        return
    end if
    if (incx .eq. 1 .and. incy .eq. 1) then
        go to 20
    end if
    ix = 1
    iy = 1
    if (incx .lt. 0) then
        ix = ((-n) + 1)*incx + 1
    end if
    if (incy .lt. 0) then
        iy = ((-n) + 1)*incy + 1
    end if
    do i = 1, n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
10      continue
    end do
    ddot = dtemp
    return
20  m = mod(n, 5)
    if (m .eq. 0) then
        go to 40
    end if
    do i = 1, m
        dtemp = dtemp + dx(i)*dy(i)
30      continue
    end do
    if (n .lt. 5) then
        go to 60
    end if
40  mp1 = m + 1
    do i = mp1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
50      continue
    end do
60  ddot = dtemp
    return
end function ddot

double precision function dnrm2(n, dx, incx)
    integer :: next
    double precision :: dx(*), cutlo, cuthi, hitest, sum, xmax, zero, one
    data zero, one/0.0d0, 1.0d0/
    data cutlo, cuthi/8.232d-11, 1.304d19/
    if (n .gt. 0) then
        go to 10
    end if
    dnrm2 = zero
    go to 300
10  assign 30 to next
    sum = zero
    nn = n*incx
    i = 1
20  go to next, (30, 50, 70, 110)
30  if (dabs(dx(i)) .gt. cutlo) then
        go to 85
    end if
    assign 50 to next
    xmax = zero
50  if (dx(i) .eq. zero) then
        go to 200
    end if
    if (dabs(dx(i)) .gt. cutlo) then
        go to 85
    end if
    assign 70 to next
    go to 105
100 i = j
    assign 110 to next
    sum = sum/dx(i)/dx(i)
105 xmax = dabs(dx(i))
    go to 115
70  if (dabs(dx(i)) .gt. cutlo) then
        go to 75
    end if
110 if (dabs(dx(i)) .le. xmax) then
        go to 115
    end if
    sum = one + sum*(xmax/dx(i))**2
    xmax = dabs(dx(i))
    go to 200
115 sum = sum + (dx(i)/xmax)**2
    go to 200
75  sum = sum*xmax*xmax
85  hitest = cuthi/float(n)
    do j = i, nn, incx
        if (dabs(dx(j)) .ge. hitest) then
            go to 100
        end if
95      sum = sum + dx(j)**2
    end do
    dnrm2 = dsqrt(sum)
    go to 300
200 continue
    i = i + incx
    if (i .le. nn) then
        go to 20
    end if
    dnrm2 = xmax*dsqrt(sum)
300 continue
    return
end function dnrm2

integer function idamax(n, dx, incx)
    double precision :: dx(*), dmax
    integer :: i, incx, ix, n
    idamax = 0
    if (n .lt. 1) then
        return
    end if
    idamax = 1
    if (n .eq. 1) then
        return
    end if
    if (incx .eq. 1) then
        go to 20
    end if
    ix = 1
    dmax = dabs(dx(1))
    ix = ix + incx
    do i = 2, n
        if (dabs(dx(ix)) .le. dmax) then
            go to 5
        end if
        idamax = i
        dmax = dabs(dx(ix))
5       ix = ix + incx
10      continue
    end do
    return
20  dmax = dabs(dx(1))
    do i = 2, n
        if (dabs(dx(i)) .le. dmax) then
            go to 30
        end if
        idamax = i
        dmax = dabs(dx(i))
30      continue
    end do
    return
end function idamax

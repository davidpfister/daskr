program dem
    implicit double precision(a - h, o - z)
    external res1, rt1, res2, jac2, rt2, jdum, psdum 
    integer :: idid, iout, iwork, jroot, info, jtype, kroot, leniw, lenrw, liw, lrw, lun, neq, nerr, nrt, nre, nrea, nrte, nje, nst
    double precision :: atol, er, ero, errt, rtol, rwork, t, tout, tzero, y, yt
    dimension :: y(2), yprime(2), rtol(2), atol(2), rwork(100), iwork(100), jroot(2), info(20), rpar(1), ipar(1)
    common/local/neq
    lun = 6
    kprint = 3
    ipass = 1
    nerr = 0
    do i = 1, 20
10      info(i) = 0
    end do
    neq = 1
    t = 1.0d0
    y(1) = 1.0d0
    tout = 2.0d0
    rtol(1) = 0.0d0
    atol(1) = 1.0d-6
    lrw = 100
    liw = 100
    idid = 0
    info(11) = 0
    yprime(1) = 3.0d0
    jtype = 2
    info(5) = 2 - jtype
    nrt = 2
    if (kprint .ge. 2) then
        write (lun, 110) rtol(1), atol(1), jtype
    end if
110 format(' DKRDEM: Demonstration Program for DDASKR'/                                 &
           ' Problem 1.' /                                                              &
           ' Problem is  dY/dT = ((2*LOG(Y)+8)/T - 5)*Y,  Y(1) = 1'/                    &
           ' Solution is  Y(T) = EXP(-T**2 + 5*T - 4)'/                                 &
           ' Root functions are'/                                                       &
           ' R1 = dY/dT  (root at T = 2.5)'/                                            & 
           ' R2 = LOG(Y) - 2.2491  (roots at T = 2.47 and T = 2.53)'/                   & 
           ' RTOL =',e10.1,'   ATOL =',e10.1,'   JTYPE =',i3)
    ero = 0.0d0
    do iout = 1, 5
        call ddaskr(res1, neq, t, y, yprime, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, jdum, psdum, rt1, nrt, jroot)
        yt = exp((-t)*t + 5.0d0*t - 4.0d0)
        er = y(1) - yt
        if (kprint .gt. 2) then
            write (lun, 130) t, y(1), er
        end if
130     format(' At t =', e15.7, 5x, 'y =', e15.7, 5x, 'error =', e12.4)
        if (idid .lt. 0) then
            go to 185
        end if
        er = abs(er)/atol(1)
        ero = max(ero, er)
        if (er .lt. 1000.0d0) then
            go to 140
        end if
        if (kprint .ge. 2) then
            ipass = 0
            write (lun, 135)
        end if
135     format(' WARNING. Error exceeds 1000 * tolerance')
        nerr = nerr + 1
140     continue
        if (idid .ne. 5) then
            go to 175
        end if
        if (kprint .gt. 2) then
            write (lun, 150) t, jroot(1), jroot(2)
        end if
150     format('      Root found at t =', e15.7, 5x, 'JROOT =', 2i5)
        if (jroot(1) .ne. 0) then
            errt = t - 2.5d0
        end if
        if (jroot(2) .ne. 0 .and. t .le. 2.5d0) then
            errt = t - 2.47d0
        end if
        if (jroot(2) .ne. 0 .and. t .gt. 2.5d0) then
            errt = t - 2.53d0
        end if
        if (kprint .gt. 2) then
            write (lun, 160) errt
        end if
160     format('      Error in t location of root is', e12.4)
        if (abs(errt) .lt. 1.0d-3) then
            go to 170
        end if
        if (kprint .ge. 2) then
            ipass = 0
            write (lun, 165)
        end if
165     format(' WARNING.. Root error exceeds 1.0D-3')
        nerr = nerr + 1
170     continue
        cycle
175     tout = tout + 1.0d0
180     continue
    end do
185 continue
    if (idid .lt. 0) then
        nerr = nerr + 1
    end if
    nst = iwork(11)
    nre = iwork(12)
    nje = iwork(13)
    nrte = iwork(36)
    lenrw = 0
    leniw = 0
    nrea = nre
    if (jtype .eq. 2) then
        nre = nre + neq*nje
    end if
    if (kprint .gt. 2) then
        write (lun, 190) nst, nre, nrea, nje, nrte, ero
    end if
190 format(' Final statistics for this run.'/ & 
           ' number of steps =',i5/             &
           ' number of Gs    =',i5/             &
           ' (excluding Js)  =',i5/             & 
           ' number of Js    =',i5/             & 
           ' number of Rs    =',i5/             & 
           ' error overrun   =',e10.2)
    do i = 1, 20
195     info(i) = 0
    end do
    info(2) = 1
    rtol(1) = 1.0d-6
    rtol(2) = 1.0d-6
    atol(1) = 1.0d-6
    atol(2) = 1.0d-4
    if (kprint .ge. 2) then
        write (lun, 200) rtol(1), atol(1), atol(2)
    end if
200 format(80('-')//' Problem 2. Van Der Pol oscillator'/               &
           ' Problem is dY1/dT = Y2,  dY2/dT = 100*(1-Y1**2)*Y2 - Y1'/  &
           '            Y1(0) = 2,  Y2(0) = 0'/                         & 
           ' Root function is  R(T,Y,YP) = Y1'/                         & 
           ' RTOL =',e10.1,'   ATOL =',2e10.1)
    do jtype = 1, 2
        info(1) = 0
        info(5) = 2 - jtype
        neq = 2
        t = 0.0d0
        y(1) = 2.0d0
        y(2) = 0.0d0
        yprime(1) = 0.d0
        yprime(2) = -2.0d0
        tout = 20.0d0
        nrt = 1
        if (kprint .gt. 2) then
            write (lun, 210) jtype
        end if
210     format(80('.')//' Solution with JTYPE =', i2)
        do iout = 1, 10
        call ddaskr(res2, neq, t, y, yprime, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, jac2, psdum, rt2, nrt, jroot)
            if (kprint .gt. 2) then
                write (lun, 230) t, y(1), y(2)
            end if
230         format(' At t =', e15.7, 5x, 'y1 =', e15.7, 5x, 'y2 =', e15.7)
            if (idid .lt. 0) then
                go to 275
            end if
            if (idid .ne. 5) then
                go to 265
            end if
            if (kprint .gt. 2) then
                write (lun, 240) t, jroot(1)
            end if
240         format('      Root found at t =', e15.7, '  JROOT =', i3)
            kroot = int(t/81.2d0 + 0.5d0)
            tzero = 81.17237787055d0 + dfloat(kroot - 1)*81.41853556212d0
            errt = t - tzero
            if (kprint .gt. 2) then
                write (lun, 250) errt
            end if
250         format('      Error in t location of root is', e12.4)
            if (errt .lt. 1.0d0) then
                go to 260
            end if
            if (kprint .ge. 2) then
                ipass = 0
                write (lun, 255)
            end if
255         format(' WARNING. Root error exceeds 1.0')
            nerr = nerr + 1
260         continue
            cycle
265         tout = tout + 20.0d0
270         continue
        end do
275     continue
        if (idid .lt. 0) then
            nerr = nerr + 1
        end if
        nst = iwork(11)
        nre = iwork(12)
        nje = iwork(13)
        nrte = iwork(36)
        lenrw = 0
        leniw = 0
        nrea = nre
        if (jtype .eq. 2) then
            nre = nre + neq*nje
        end if
        if (kprint .ge. 2) then
            write (lun, 280) nst, nre, nrea, nje, nrte
        end if
280 format(' Final statistics for this run.'/   & 
           ' number of steps =',i5/             &
           ' number of Gs    =',i5/             & 
           ' (excluding Js)  =',i5/             & 
           ' number of Js    =',i5/             & 
           ' number of Rs    =',i5)
290     continue
    end do
    if (kprint .ge. 2) then
        write (lun, 300) nerr
    end if
300 format(80('-')//' Number of errors encountered =', i3)
    if (nerr .gt. 0) then
        ipass = 0
    end if
    if (ipass .eq. 1 .and. kprint .gt. 1) then
        write (lun, 700)
    end if
    if (ipass .eq. 0 .and. kprint .ne. 0) then
        write (lun, 800)
    end if
700 format(' **********DDASKR passed all tests**********')
800 format(' **********DDASKR failed some tests*********')
    stop
end program

subroutine res1(t, y, yprime, cj, delta, ires, rpar, ipar)
    implicit double precision(a - h, o - z)
    common/local/neq
    dimension :: y(neq), yprime(neq), delta(neq)

    if (y(1) .le. 0.0d0) then
        ires = -1
        return
    else
        call f1(neq, t, y, delta)
        do i = 1, neq
            delta(i) = yprime(i) - delta(i)
10          continue
        end do
    end if
    return
end subroutine res1

subroutine f1(neq, t, y, ydot)
    implicit double precision(a - h, o - z)
    integer :: neq
    double precision :: t, y, ydot
    dimension :: y(1), ydot(1)
    ydot(1) = ((2.0d0*log(y(1)) + 8.0d0)/t - 5.0d0)*y(1)
    return
end subroutine f1

subroutine rt1(neq, t, y, yp, nrt, rval, rpar, ipar)
    implicit double precision(a - h, o - z)
    integer :: neq, nrt
    double precision :: t, y, rval
    dimension :: rval(2), y(1)
    rval(1) = yp
    rval(2) = log(y(1)) - 2.2491d0
    return
end subroutine rt1

subroutine res2(t, y, yprime, cj, delta, ires, rpar, ipar)
    implicit double precision(a - h, o - z)
    common/local/neq
    dimension :: y(neq), yprime(neq), delta(neq)
    call f2(neq, t, y, delta)
    do i = 1, neq
10      delta(i) = yprime(i) - delta(i)
    end do
    return
end subroutine res2

subroutine f2(neq, t, y, ydot)
    implicit double precision(a - h, o - z)
    integer :: neq
    double precision :: t, y, ydot
    dimension :: y(2), ydot(2)
    ydot(1) = y(2)
    ydot(2) = 100.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
    return
end subroutine f2

subroutine jac2(t, y, yprime, pd, cj, rpar, ipar)
    implicit double precision(a - h, o - z)
    integer :: neq, nrowpd
    double precision :: t, y, pd
    parameter(nrowpd=2)
    dimension :: y(2), pd(nrowpd, 2)
    common/local/neq
    pd(1, 1) = 0.0d0
    pd(1, 2) = 1.0d0
    pd(2, 1) = (-200.0d0)*y(1)*y(2) - 1.0d0
    pd(2, 2) = 100.0d0*(1.0d0 - y(1)*y(1))
    pd(1, 1) = cj - pd(1, 1)
    pd(1, 2) = -pd(1, 2)
    pd(2, 1) = -pd(2, 1)
    pd(2, 2) = cj - pd(2, 2)
    return
end subroutine jac2

subroutine rt2(neq, t, y, yp, nrt, rval, rpar, ipar)
    implicit double precision(a - h, o - z)
    integer :: neq, nrt
    double precision :: t, y, rval
    dimension :: y(2), yp(2), rval(1)
    rval(1) = y(1)
    return
end subroutine rt2

subroutine jdum()
end subroutine

subroutine psdum()
end subroutine

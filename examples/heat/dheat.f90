program heat
    implicit double precision(a - h, o - z)
    external :: resh, rtheat, dbanja, dbanps
    parameter(maxm=10, maxm2=maxm + 2, mxneq=maxm2*maxm2)
    parameter(lenrw=107 + 18*mxneq, leniw=40, lenwp=4*mxneq + 2*(mxneq/3 + 1), leniwp=mxneq)
    dimension :: u(mxneq), uprime(mxneq), rwork(lenrw + lenwp), iwork(leniw + leniwp)
    dimension :: info(20), jroot(2), rpar(2), ipar(4)
    m = maxm
    dx = 1.0d0/(m + 1)
    neq = (m + 2)*(m + 2)
    coeff = 1.0d0/(dx*dx)
    ipar(3) = neq
    ipar(4) = m
    rpar(1) = dx
    rpar(2) = coeff
    nrt = 2
    ml = 1
    mu = 1
    ipar(1) = ml
    ipar(2) = mu
    lenpd = (2*ml + mu + 1)*neq
    mband = ml + mu + 1
    msave = neq/mband + 1
    lwp = lenpd + 2*msave
    liwp = neq
    iwork(27) = lwp
    iwork(28) = liwp
    lrw = lenrw + lwp
    liw = leniw + liwp
    call uinit(u, uprime, rpar, ipar)
    do i = 1, 20
10      info(i) = 0
    end do
    info(12) = 1
    info(15) = 1
    rtol = 0.0d0
    atol = 1.0d-5
    lout = 6
    write (lout, 30) m, neq, info(12), ml, mu, rtol, atol
    30   format(' DHEAT: Heat Equation Example Program for DDASKR'/         &
                '    M+2 by M+2 mesh, M =',i3,',  System size NEQ =',i4/    &
                '    Root functions are: R1 = max(u) - 0.1'/                &
                ' and R2 = max(u) - 0.01'/                                  &
                '    Linear solver method flag INFO(12) =',i3/              &
                '    (0 = direct, 1 = Krylov)'/                             &         
                '    Preconditioner is a banded approximation with ML =', i3,'  MU =',i3/ &                       
                '    Tolerances are RTOL =',e10.1,'   ATOL =',e10.1)
    write (lout, 40)
40  format(5x, 't', 12x, 'UMAX', 8x, 'NQ', 8x, 'H', 8x, 'STEPS', 5x, 'NNI', 5x, 'NLI')
    nout = 11
    t = 0.0d0
    tout = 0.01d0
    do iout = 1, nout
        call ddaskr(resh, neq, t, u, uprime, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, dbanja, dbanps, rtheat, nrt, jroot)
        umax = 0.0d0
        do i = 1, neq
50          umax = max(umax, abs(u(i)))
        end do
        hu = rwork(7)
        nqu = iwork(8)
        nst = iwork(11)
        nni = iwork(19)
        nli = iwork(20)
        write (lout, 60) t, umax, nqu, hu, nst, nni, nli
60      format(e15.5, e12.4, i5, e14.3, i7, i9, i8)
        if (idid .eq. 5) then
            write (6, 61) jroot(1), jroot(2)
61          format(20x, '*****   Root found, JROOT =', 2i3)
            cycle
        end if
        if (idid .lt. 0) then
            write (lout, 65) t
65          format(//' Final time reached =', e12.4//)
            go to 80
        end if
        tout = tout*2.0d0
70      continue
    end do
80  continue
    nst = iwork(11)
    npe = iwork(13)
    nre = iwork(12) + npe*mband
    liw = iwork(17)
    lrw = iwork(18)
    nni = iwork(19)
    nli = iwork(20)
    nps = iwork(21)
    if (nni .ne. 0) then
        avdim = real(nli)/real(nni)
    end if
    ncfn = iwork(15)
    ncfl = iwork(16)
    nrte = iwork(36)
    write (lout, 90) lrw, liw, nst, nre, nrte, npe, nps, nni, nli, avdim, ncfn, ncfl
90  format(' Final statistics for this run.'/ &
           '   RWORK size =',i5/ &
           '   IWORK size =',i4/& 
           '   Number of time steps ................ =',i5/&
           '   Number of residual evaluations ...... =',i5/&
           '   Number of root function evaluations . =',i5/&
           '   Number of preconditioner evaluations  =',i5/&
           '   Number of preconditioner solves ..... =',i5/&
           '   Number of nonlinear iterations ...... =',i5/&
           '   Number of linear iterations ......... =',i5/&
           '   Average Krylov subspace dimension =',f8.4,i5/&
           ' nonlinear conv. failures,',i5/&
           ' linear conv. failures')
    stop
end program

subroutine uinit(u, uprime, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: u(*), uprime(*), rpar(2), ipar(4)
    neq = ipar(3)
    m = ipar(4)
    dx = rpar(1)
    do k = 0, m + 1
        yk = k*dx
        ioff = (m + 2)*k
        do j = 0, m + 1
            xj = j*dx
            i = ioff + j + 1
            u(i) = 16.0d0*xj*(1.0d0 - xj)*yk*(1.0d0 - yk)
10          continue
        end do
20      continue
    end do
    do i = 1, neq
30      uprime(i) = 0.0d0
    end do
    return
end subroutine uinit

subroutine resh(t, u, uprime, cj, delta, ires, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: u(*), uprime(*), delta(*), rpar(2), ipar(4)
    neq = ipar(3)
    m = ipar(4)
    coeff = rpar(2)
    m2 = m + 2
    do i = 1, neq
10      delta(i) = u(i)
    end do
    do k = 1, m
        ioff = m2*k
        do j = 1, m
            i = ioff + j + 1
            temx = u(i - 1) + u(i + 1)
            temy = u(i - m2) + u(i + m2)
            delta(i) = uprime(i) - (temx + temy - 4.0d0*u(i))*coeff
20          continue
        end do
30      continue
    end do
    return
end subroutine resh

subroutine rtheat(neq, t, u, up, nrt, rval, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: u(neq), rval(2)
    umax = 0.0d0
    do i = 1, neq
10      umax = max(umax, u(i))
    end do
    rval(1) = umax - 0.1d0
    rval(2) = umax - 0.01d0
    return
end subroutine rtheat

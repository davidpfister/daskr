program web
    implicit double precision(a - h, o - z)
    external :: resweb, rtweb, jacrs, psolrs
    parameter(maxs=2, maxm=20, maxn=800, lrw=63 + (3*maxs*maxm + 11)*maxn + maxm, liw=40 + 2*maxn)
    dimension :: cc(maxn), ccprime(maxn), rwork(lrw), iwork(liw), info(20), rpar(maxn), ipar(2)
    common/ppar1/aa, ee, gg, bb, dprey, dpred
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    data lout/9/, lcout/10/
    open (unit=lout, file="wdout", status="unknown")
    open (unit=lcout, file="wccout", status="unknown")
    call setpar()
    neq = ns*mx*my
    mxns = mx*ns
    dx = ax/real(mx - 1)
    dy = ay/real(my - 1)
    do i = 1, ns
        cox(i) = diff(i)/dx**2
10      coy(i) = diff(i)/dy**2
    end do
    nrt = 1
    write (lout, 20) ns
20  format(' DWEB: Example program for DDASKR package'/&
          ' Food web problem with NS species, NS =',i4/&
          ' Predator-prey interaction and diffusion on a 2-D square')
    write (lout, 25) aa, ee, gg, bb, dprey, dpred, alph, beta
25  format(' Matrix parameters.  a =',e12.4/&
           '   e =',e12.4/&
           '   g =',e12.4/21x/&
           ' b parameter =',e12.4/&
           ' Diffusion coefficients: dprey =',e12.4,'   dpred =',e12.4/&
           ' Rate parameters: alpha =',e12.4,' and beta =',e12.4)
    write (lout, 30) mx, my, neq
30  format(' Mesh dimensions (MX,MY) =', 2i4, 5x, ' Total system size is NEQ =', i7)
    write (lout, 35)
35  format(' Root function is R(Y) = average(c1) - 20')
    predic = 1.0d5
    do i = 1, 20
50      info(i) = 0
    end do
    info(11) = 1
    info(14) = 1
    info(15) = 1
    info(16) = 1
    rtol = 1.0d-5
    atol = rtol
    write (lout, 70) rtol, atol, info(11), predic, info(16)
70  format(' Tolerance parameters: RTOL =',e10.2,'   ATOL =',e10.2/&
           ' Internal I.C. calculation flag INFO(11) =',i2,'   (0 = off, 1 = on)'/&
           ' Predator I.C. guess =',e10.2/&
           ' Alternate error test flag INFO(16) =',i2,'  (0 = off, 1 = on)')
    nout = 18
    do meth = 0, 2
        info(12) = min(meth, 1)
        info(1) = 0
        info(11) = 1
        jbg = meth - 1
        ipar(2) = jbg
        write (lout, 80) info(12)
80      format(//80('.')//' Linear solver method flag INFO(12) =', i2, '   (0 = direct, 1 = Krylov)')
        if (info(12) .eq. 0) then
            info(6) = 1
            iwork(1) = mxns
            iwork(2) = mxns
            call setid(mx, my, ns, np, 40, iwork)
            write (lout, 90) mxns
90          format(' Difference-quotient banded Jacobian,', ' half-bandwidths =', i4)
        end if
        if (info(12) .eq. 1) then
            jpre = 3
            ipar(1) = jpre
            write (lout, 100) jpre
100 format(' Preconditioner flag is JPRE =', i3/'  (1 = reaction factor A_R, 2 = spatial factor A_S,', ' 3 = A_S*A_R, 4 = A_R*A_S )')
            if (jbg .eq. 0) then
                call dmset2(mx, my, ns, np, 40, iwork)
                write (lout, 110)
110             format(' No block-grouping in reaction factor')
            end if
            if (jbg .eq. 1) then
                nxg = 5
                nyg = 5
                ng = nxg*nyg
                call dgset2(mx, my, ns, np, nxg, nyg, 40, iwork)
                write (lout, 120) ng, nxg, nyg
120         format(' Block-grouping in reaction factor'/' Number of groups =', i5, '   (NGX by NGY, NGX =', i3, ',  NGY =', i3, ')')
            end if
        end if
        t = 0.0d0
        tout = 1.0d-8
        call cinit(cc, ccprime, predic, rpar)
        nli = 0
        nni = 0
        write (lout, 140)
140     format(//'   t', 12x, 'Ave.c1  NSTEP  NRE  NNI  NLI  NPE  NQ', 4x, 'H', 10x, 'AVLIN')
        do iout = 0, nout
    call ddaskr(resweb, neq, t, cc, ccprime, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, jacrs, psolrs, rtweb, nrt, jroot)
            nst = iwork(11)
            nre = iwork(12)
            npe = iwork(13)
            nnidif = iwork(19) - nni
            nni = iwork(19)
            nlidif = iwork(20) - nli
            nli = iwork(20)
            nqu = iwork(8)
            hu = rwork(7)
            avlin = 0.0d0
            if (nnidif .gt. 0) then
                avlin = real(nlidif)/real(nnidif)
            end if
            if (meth .eq. 0) then
                imod3 = iout - 3*iout/3
                if (imod3 .eq. 0) then
                    call outweb(t, cc, ns, mx, my, lcout)
                end if
            end if
            call avc1(cc, c1ave)
            write (lout, 160) t, c1ave, nst, nre, nni, nli, npe, nqu, hu, avlin
160         format(e13.5, f10.5, i5, i6, 3i5, i4, e11.2, f9.4)
            if (idid .eq. 5) then
                write (lout, 165) jroot
165             format(15x, '*****   Root found, JROOT =', i3)
                cycle
            end if
            if (idid .lt. 0) then
                write (lout, 170) t
170             format(//' Final time reached =', e12.4//)
                go to 210
            end if
            if (tout .gt. 0.9d0) then
                tout = tout + 1.0d0
            end if
            if (tout .lt. 0.9d0) then
                tout = tout*10.0d0
            end if
            if (iout .eq. 0) then
                info(11) = 0
                tout = 1.0d-8
                nli = 0
                nni = 0
            end if
200         continue
        end do
210     continue
        lenrw = iwork(18)
        leniw = iwork(17)
        nst = iwork(11)
        nre = iwork(12)
        npe = iwork(13)
        nni = iwork(19)
        nli = iwork(20)
        nps = iwork(21)
        if (nni .gt. 0) then
            avlin = real(nli)/real(nni)
        end if
        ncfn = iwork(15)
        ncfl = iwork(16)
        nrte = iwork(36)
        write (lout, 220) lenrw, leniw, nst, nre, nrte, npe, nps, nni, nli, avlin, ncfn, ncfl
220     format(' Final statistics for this run:'/&
           ' RWORK size =',i8,'   IWORK size =',i6/&
           ' Number of time steps            =',i5/&
           ' Number of residual evaluations  =',i5/&
           ' Number of root fn. evaluations  =',i5/&
           ' Number of Jac. or prec. evals.  =',i5/&
           ' Number of preconditioner solves =',i5/&
           ' Number of nonlinear iterations  =',i5/&
           ' Number of linear iterations     =',i5/&
           ' Average Krylov subspace dimension =',f8.4/i3/&
           ' nonlinear conv. failures,',i5/&
           ' linear conv. failures')
300     continue
    end do
    stop
end program

subroutine setpar()
    implicit double precision(a - h, o - z)
    parameter(maxs=2)
    common/ppar1/aa, ee, gg, bb, dprey, dpred
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    ax = 1.0d0
    ay = 1.0d0
    np = 1
    mx = 20
    my = 20
    aa = 1.0d0
    ee = 1.0d4
    gg = 0.5d-6
    bb = 1.0d0
    dprey = 1.0d0
    dpred = 0.05d0
    alph = 50.0d0
    beta = 100.0d0
    ns = 2*np
    do j = 1, np
        do i = 1, np
            acoef(np + i, j) = ee
            acoef(i, np + j) = -gg
10          continue
        end do
        acoef(j, j) = -aa
        acoef(np + j, np + j) = -aa
        bcoef(j) = bb
        bcoef(np + j) = -bb
        diff(j) = dprey
        diff(np + j) = dpred
20      continue
    end do
    pi = 3.141592653589793d0
    fpi = 4.0d0*pi
    return
end subroutine setpar

subroutine setid(mx, my, ns, nsd, lid, iwork)
    dimension :: iwork(*)
    nsdp1 = nsd + 1
    do jy = 1, my
        i00 = mx*ns*(jy - 1) + lid
        do jx = 1, mx
            i0 = i00 + ns*(jx - 1)
            do i = 1, nsd
10              iwork(i0 + i) = 1
            end do
            do i = nsdp1, ns
20              iwork(i0 + i) = -1
            end do
30          continue
        end do
40      continue
    end do
    return
end subroutine setid

subroutine cinit(cc, ccprime, predic, rpar)
    implicit double precision(a - h, o - z)
    dimension :: cc(*), ccprime(*), rpar(*)
    parameter(maxs=2)
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    npp1 = np + 1
    do jy = 1, my
        y = real(jy - 1)*dy
        argy = 16.0d0*y*y*(ay - y)*(ay - y)
        iyoff = mxns*(jy - 1)
        do jx = 1, mx
            x = real(jx - 1)*dx
            argx = 16.0d0*x*x*(ax - x)*(ax - x)
            ioff = iyoff + ns*(jx - 1)
            fac = 1.0d0 + alph*x*y + beta*sin(fpi*x)*sin(fpi*y)
            do i = 1, np
10              cc(ioff + i) = 10.0d0 + real(i)*argx*argy
            end do
            do i = npp1, ns
15              cc(ioff + i) = predic
            end do
20          continue
        end do
30      continue
    end do
    t = 0.0d0
    call fweb(t, cc, ccprime, rpar)
    do jy = 1, my
        iyoff = mxns*(jy - 1)
        do jx = 1, mx
            ioff = iyoff + ns*(jx - 1)
            do i = npp1, ns
40              ccprime(ioff + i) = 0.0d0
            end do
50          continue
        end do
60      continue
    end do
    return
end subroutine cinit

subroutine outweb(t, c, ns, mx, my, lun)
    implicit double precision(a - h, o - z)
    dimension :: c(ns, mx, my)
    write (lun, 10) t
10  format(1x, 79('-')/30x, 'At time t = ', e16.8/1x, 79('-'))
    do i = 1, ns
        write (lun, 20) i
20      format(' the species c(', i2, ') values are ')
        do jy = my, 1, -1
            write (lun, 25) (c(i, jx, jy), jx=1, mx)
25          format(6(1x, g12.6))
30          continue
        end do
        write (lun, 35)
35      format(1x, 79('-'))
40      continue
    end do
    return
end subroutine outweb

subroutine resweb(t, u, uprime, cj, delta, ires, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: u(*), uprime(*), delta(*), rpar(*), ipar(*)
    parameter(maxs=2)
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    call fweb(t, u, delta, rpar)
    do jy = 1, my
        iyoff = mxns*(jy - 1)
        do jx = 1, mx
            ic0 = iyoff + ns*(jx - 1)
            do i = 1, ns
                ici = ic0 + i
                if (i .gt. np) then
                    delta(ici) = -delta(ici)
                else
                    delta(ici) = uprime(ici) - delta(ici)
                end if
10              continue
            end do
20          continue
        end do
30      continue
    end do
    return
end subroutine resweb

subroutine fweb(t, cc, crate, rpar)
    implicit double precision(a - h, o - z)
    dimension :: cc(*), crate(*), rpar(*)
    parameter(maxs=2)
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    do jy = 1, my
        iyoff = mxns*(jy - 1)
        idyu = mxns
        if (jy .eq. my) then
            idyu = -mxns
        end if
        idyl = mxns
        if (jy .eq. 1) then
            idyl = -mxns
        end if
        do jx = 1, mx
            ic = iyoff + ns*(jx - 1) + 1
            call webr(t, jx, jy, cc(ic), rpar(ic))
            idxu = ns
            if (jx .eq. mx) then
                idxu = -ns
            end if
            idxl = ns
            if (jx .eq. 1) then
                idxl = -ns
            end if
            do i = 1, ns
                ici = ic + i - 1
                dcyli = cc(ici) - cc(ici - idyl)
                dcyui = cc(ici + idyu) - cc(ici)
                dcxli = cc(ici) - cc(ici - idxl)
                dcxui = cc(ici + idxu) - cc(ici)
                crate(ici) = coy(i)*(dcyui - dcyli) + cox(i)*(dcxui - dcxli) + rpar(ici)
20              continue
            end do
40          continue
        end do
60      continue
    end do
    return
end subroutine fweb

subroutine webr(t, jx, jy, c, crate)
    implicit double precision(a - h, o - z)
    dimension :: c(*), crate(*)
    parameter(maxs=2)
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    y = real(jy - 1)*dy
    x = real(jx - 1)*dx
    do i = 1, ns
10      crate(i) = 0.0d0
    end do
    do j = 1, ns
        call daxpy(ns, c(j), acoef(1, j), 1, crate, 1)
15      continue
    end do
    fac = 1.0d0 + alph*x*y + beta*sin(fpi*x)*sin(fpi*y)
    do i = 1, ns
20      crate(i) = c(i)*(bcoef(i)*fac + crate(i))
    end do
    return
end subroutine webr

subroutine jacrs(res, ires, neq, t, cc, ccprime, rewt, savr, wk, h, cj, wp, iwp, ier, rpar, ipar)
    implicit double precision(a - h, o - z)
    external :: webr
    dimension :: cc(*), ccprime(*), rewt(*), savr(*), wk(*), wp(*), iwp(*), rpar(*), ipar(*)
    jbg = ipar(2)
    if (jbg .eq. 0) then
        call drbdja(t, cc, rpar, webr, wk, rewt, cj, wp, iwp, ier)
    else
        call drbgja(t, cc, rpar, webr, wk, rewt, cj, wp, iwp, ier)
    end if
    return
end subroutine jacrs

subroutine psolrs(neq, t, cc, ccprime, savr, wk, cj, wt, wp, iwp, b, eplin, ier, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: cc(*), ccprime(*), savr(*), wk(*), wp(*), iwp(*), b(*), rpar(*), ipar(*)
    jpre = ipar(1)
    ier = 0
    hl0 = 1.0d0/cj
    jbg = ipar(2)
    if (jpre .eq. 2 .or. jpre .eq. 3) then
        call gs(neq, hl0, b, wk)
    end if
    if (jpre .ne. 2) then
        if (jbg .eq. 0) then
            call drbdps(b, wp, iwp)
        end if
        if (jbg .eq. 1) then
            call drbgps(b, wp, iwp)
        end if
    end if
    if (jpre .eq. 4) then
        call gs(neq, hl0, b, wk)
    end if
    return
end subroutine psolrs

subroutine gs(n, hl0, z, x)
    implicit double precision(a - h, o - z)
    dimension :: z(*), x(*)
    parameter(maxs=2)
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    dimension :: beta1(maxs), gamma1(maxs), beta2(maxs), gamma2(maxs), dinv(maxs)
    data itmax/5/
    do i = 1, ns
        elamda = 1.0d0/(1.0d0 + 2.0d0*hl0*(cox(i) + coy(i)))
        beta1(i) = hl0*cox(i)*elamda
        beta2(i) = 2.0d0*beta1(i)
        gamma1(i) = hl0*coy(i)*elamda
        gamma2(i) = 2.0d0*gamma1(i)
        dinv(i) = elamda
10      continue
    end do
    iter = 1
    do ii = 1, n
15      x(ii) = 0.0d0
    end do
    do jy = 1, my
        iyoff = mxns*(jy - 1)
        do jx = 1, mx
            ic = iyoff + ns*(jx - 1)
            do i = 1, ns
                ici = ic + i
                x(ici) = dinv(i)*z(ici)
                z(ici) = 0.0d0
30              continue
            end do
40          continue
        end do
50      continue
    end do
    go to 160
70  continue
    iter = iter + 1
    jy = 1
    jx = 1
    ic = ns*(jx - 1)
    do i = 1, ns
        ici = ic + i
75      x(ici) = beta2(i)*x(ici + ns) + gamma2(i)*x(ici + mxns)
    end do
    do jx = 2, mx - 1
        ic = ns*(jx - 1)
        do i = 1, ns
            ici = ic + i
80          x(ici) = beta1(i)*x(ici + ns) + gamma2(i)*x(ici + mxns)
        end do
85      continue
    end do
    jx = mx
    ic = ns*(jx - 1)
    do i = 1, ns
        ici = ic + i
90      x(ici) = gamma2(i)*x(ici + mxns)
    end do
    do jy = 2, my - 1
        iyoff = mxns*(jy - 1)
        jx = 1
        ic = iyoff
        do i = 1, ns
            ici = ic + i
95          x(ici) = beta2(i)*x(ici + ns) + gamma1(i)*x(ici + mxns)
        end do
        do jx = 2, mx - 1
            ic = iyoff + ns*(jx - 1)
            do i = 1, ns
                ici = ic + i
100             x(ici) = beta1(i)*x(ici + ns) + gamma1(i)*x(ici + mxns)
            end do
105         continue
        end do
        jx = mx
        ic = iyoff + ns*(jx - 1)
        do i = 1, ns
            ici = ic + i
110         x(ici) = gamma1(i)*x(ici + mxns)
        end do
115     continue
    end do
    jy = my
    iyoff = mxns*(jy - 1)
    jx = 1
    ic = iyoff
    do i = 1, ns
        ici = ic + i
120     x(ici) = beta2(i)*x(ici + ns)
    end do
    do jx = 2, mx - 1
        ic = iyoff + ns*(jx - 1)
        do i = 1, ns
            ici = ic + i
125         x(ici) = beta1(i)*x(ici + ns)
        end do
130     continue
    end do
    jx = mx
    ic = iyoff + ns*(jx - 1)
    do i = 1, ns
        ici = ic + i
135     x(ici) = 0.0d0
    end do
160 continue
    jy = 1
    do jx = 2, mx - 1
        ic = ns*(jx - 1)
        do i = 1, ns
            ici = ic + i
170         x(ici) = x(ici) + beta1(i)*x(ici - ns)
        end do
175     continue
    end do
    jx = mx
    ic = ns*(jx - 1)
    do i = 1, ns
        ici = ic + i
180     x(ici) = x(ici) + beta2(i)*x(ici - ns)
    end do
    do jy = 2, my - 1
        iyoff = mxns*(jy - 1)
        jx = 1
        ic = iyoff
        do i = 1, ns
            ici = ic + i
185         x(ici) = x(ici) + gamma1(i)*x(ici - mxns)
        end do
        do jx = 2, mx - 1
            ic = iyoff + ns*(jx - 1)
            do i = 1, ns
                ici = ic + i
                x(ici) = x(ici) + beta1(i)*x(ici - ns) + gamma1(i)*x(ici - mxns)
195             continue
            end do
200         continue
        end do
        jx = mx
        ic = iyoff + ns*(jx - 1)
        do i = 1, ns
            ici = ic + i
            x(ici) = x(ici) + beta2(i)*x(ici - ns) + gamma1(i)*x(ici - mxns)
205         continue
        end do
210     continue
    end do
    jy = my
    iyoff = mxns*(jy - 1)
    jx = 1
    ic = iyoff
    do i = 1, ns
        ici = ic + i
215     x(ici) = x(ici) + gamma2(i)*x(ici - mxns)
    end do
    do jx = 2, mx - 1
        ic = iyoff + ns*(jx - 1)
        do i = 1, ns
            ici = ic + i
            x(ici) = x(ici) + beta1(i)*x(ici - ns) + gamma2(i)*x(ici - mxns)
220         continue
        end do
225     continue
    end do
    jx = mx
    ic = iyoff + ns*(jx - 1)
    do i = 1, ns
        ici = ic + i
        x(ici) = x(ici) + beta2(i)*x(ici - ns) + gamma2(i)*x(ici - mxns)
230     continue
    end do
    do i = 1, n
300     z(i) = z(i) + x(i)
    end do
    if (iter .lt. itmax) then
        go to 70
    end if
    return
end subroutine gs

subroutine avc1(cc, c1ave)
    implicit double precision(a - h, o - z)
    dimension :: cc(*)
    parameter(maxs=2)
common/ppar2/np, ns, ax, ay, acoef(maxs, maxs), bcoef(maxs), dx, dy, alph, beta, fpi, diff(maxs), cox(maxs), coy(maxs), mx, my, mxns
    sum = 0.0d0
    npp1 = np + 1
    do jy = 1, my
        iyoff = mxns*(jy - 1)
        do jx = 1, mx
            ioff = iyoff + ns*(jx - 1)
            sum = sum + cc(ioff + 1)
20          continue
        end do
30      continue
    end do
    c1ave = sum/(mx*my)
    return
end subroutine avc1

subroutine rtweb(neq, t, cc, cp, nrt, rval, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: cc(neq), cp(neq)
    call avc1(cc, c1ave)
    rval = c1ave - 20.0d0
    return
end subroutine rtweb

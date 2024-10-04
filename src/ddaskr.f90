subroutine ddaskr(res, neq, t, y, yprime, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, rpar, ipar, jac, psol, rt, nrt, jroot)
    implicit double precision(a - h, o - z)
    logical :: done, lavl, lcfn, lcfl, lwarn
    dimension :: y(*), yprime(*)
    dimension :: info(20)
    dimension :: rwork(lrw), iwork(liw)
    dimension :: rtol(*), atol(*)
    dimension :: rpar(*), ipar(*)
    dimension :: jroot(*)
    character :: msg*80
    external :: res, jac, psol, rt, ddasid, ddasik, dnedd, dnedk
parameter(lml = 1, lmu = 2, lmtype = 4, liwm = 1, lmxord = 3, ljcalc = 5, lphase = 6, lk = 7, lkold = 8, lns = 9, lnstl = 10, lnst = 11, lnre = 12, lnje = 13, letf = 14, lncfn = 15, lncfl = 16, lniw = 17, lnrw = 18, lnni = 19, lnli = 20, lnps = 21, lnpd = 22, lmiter = 23, lmaxl = 24, lkmp = 25, lnrmax = 26, llnwp = 27, llniwp = 28, llocwp = 29, llciwp = 30, lkprin = 31, lmxnit = 32, lmxnj = 33, lmxnh = 34, llsoff = 35, lnrte = 36, lirfnd = 37, licns = 41)
parameter(ltstop = 1, lhmax = 2, lh = 3, ltn = 4, lcj = 5, lcjold = 6, lhold = 7, ls = 8, lround = 9, lepli = 10, lsqrn = 11, lrsqrn = 12, lepcon = 13, lstol = 14, lepin = 15, lalpha = 21, lbeta = 27, lgamma = 33, lpsi = 39, lsigma = 45, lt0 = 51, ltlast = 52, ldelta = 61)
    save :: lid, lenid, nonneg, ncphi
    if (info(1) .ne. 0) then
        go to 100
    end if
    do i = 2, 9
        itemp = i
        if (info(i) .ne. 0 .and. info(i) .ne. 1) then
            go to 701
        end if
10      continue
    end do
    itemp = 10
    if (info(10) .lt. 0 .or. info(10) .gt. 3) then
        go to 701
    end if
    itemp = 11
    if (info(11) .lt. 0 .or. info(11) .gt. 2) then
        go to 701
    end if
    do i = 12, 17
        itemp = i
        if (info(i) .ne. 0 .and. info(i) .ne. 1) then
            go to 701
        end if
15      continue
    end do
    itemp = 18
    if (info(18) .lt. 0 .or. info(18) .gt. 2) then
        go to 701
    end if
    if (neq .le. 0) then
        go to 702
    end if
    mxord = 5
    if (info(9) .ne. 0) then
        mxord = iwork(lmxord)
        if (mxord .lt. 1 .or. mxord .gt. 5) then
            go to 703
        end if
    end if
    iwork(lmxord) = mxord
    icnflg = 0
    nonneg = 0
    lid = licns
    if (info(10) .eq. 0) then
        go to 20
    end if
    if (info(10) .eq. 1) then
        icnflg = 1
        nonneg = 0
        lid = licns + neq
    else
        if (info(10) .eq. 2) then
            icnflg = 0
            nonneg = 1
        else
            icnflg = 1
            nonneg = 1
            lid = licns + neq
        end if
    end if
20  continue
    if (info(12) .eq. 0) then
        go to 25
    end if
    iwork(lmiter) = info(12)
    if (info(13) .eq. 0) then
        iwork(lmaxl) = min(5, neq)
        iwork(lkmp) = iwork(lmaxl)
        iwork(lnrmax) = 5
        rwork(lepli) = 0.05d0
    else
        if (iwork(lmaxl) .lt. 1 .or. iwork(lmaxl) .gt. neq) then
            go to 720
        end if
        if (iwork(lkmp) .lt. 1 .or. iwork(lkmp) .gt. iwork(lmaxl)) then
            go to 721
        end if
        if (iwork(lnrmax) .lt. 0) then
            go to 722
        end if
        if (rwork(lepli) .le. 0.0d0 .or. rwork(lepli) .ge. 1.0d0) then
            go to 723
        end if
    end if
25  continue
    if (info(11) .eq. 0) then
        go to 30
    end if
    if (info(17) .eq. 0) then
        iwork(lmxnit) = 5
        if (info(12) .gt. 0) then
            iwork(lmxnit) = 15
        end if
        iwork(lmxnj) = 6
        if (info(12) .gt. 0) then
            iwork(lmxnj) = 2
        end if
        iwork(lmxnh) = 5
        iwork(llsoff) = 0
        rwork(lepin) = 0.01d0
    else
        if (iwork(lmxnit) .le. 0) then
            go to 725
        end if
        if (iwork(lmxnj) .le. 0) then
            go to 725
        end if
        if (iwork(lmxnh) .le. 0) then
            go to 725
        end if
        lsoff = iwork(llsoff)
        if (lsoff .lt. 0 .or. lsoff .gt. 1) then
            go to 725
        end if
        if (rwork(lepin) .le. 0.0d0) then
            go to 725
        end if
    end if
30  continue
    lenic = 0
    if (info(10) .eq. 1 .or. info(10) .eq. 3) then
        lenic = neq
    end if
    lenid = 0
    if (info(11) .eq. 1 .or. info(16) .eq. 1) then
        lenid = neq
    end if
    if (info(12) .eq. 0) then
        ncphi = max(mxord + 1, 4)
        if (info(6) .eq. 0) then
            lenpd = neq**2
            lenrw = 60 + 3*nrt + (ncphi + 3)*neq + lenpd
            if (info(5) .eq. 0) then
                iwork(lmtype) = 2
            else
                iwork(lmtype) = 1
            end if
        else
            if (iwork(lml) .lt. 0 .or. iwork(lml) .ge. neq) then
                go to 717
            end if
            if (iwork(lmu) .lt. 0 .or. iwork(lmu) .ge. neq) then
                go to 718
            end if
            lenpd = (2*iwork(lml) + iwork(lmu) + 1)*neq
            if (info(5) .eq. 0) then
                iwork(lmtype) = 5
                mband = iwork(lml) + iwork(lmu) + 1
                msave = neq/mband + 1
                lenrw = 60 + 3*nrt + (ncphi + 3)*neq + lenpd + 2*msave
            else
                iwork(lmtype) = 4
                lenrw = 60 + 3*nrt + (ncphi + 3)*neq + lenpd
            end if
        end if
        leniw = 40 + lenic + lenid + neq
        lenwp = 0
        leniwp = 0
    else
        if (info(12) .eq. 1) then
            ncphi = mxord + 1
            maxl = iwork(lmaxl)
            lenwp = iwork(llnwp)
            leniwp = iwork(llniwp)
            lenpd = (maxl + 3 + min0(1, maxl - iwork(lkmp)))*neq + (maxl + 3)*maxl + 1 + lenwp
            lenrw = 60 + 3*nrt + (mxord + 5)*neq + lenpd
            leniw = 40 + lenic + lenid + leniwp
        end if
    end if
    if (info(16) .ne. 0) then
        lenrw = lenrw + neq
    end if
    iwork(lniw) = leniw
    iwork(lnrw) = lenrw
    iwork(lnpd) = lenpd
    iwork(llocwp) = lenpd - lenwp + 1
    if (lrw .lt. lenrw) then
        go to 704
    end if
    if (liw .lt. leniw) then
        go to 705
    end if
    if (lenic .gt. 0) then
        do i = 1, neq
            ici = iwork(licns - 1 + i)
            if (ici .lt. -2 .or. ici .gt. 2) then
                go to 726
            end if
40          continue
        end do
    end if
    if (lenic .gt. 0) then
        call dcnst0(neq, y, iwork(licns), iret)
        if (iret .ne. 0) then
            go to 727
        end if
    end if
    index = 1
    if (lenid .gt. 0) then
        index = 0
        do i = 1, neq
            idi = iwork(lid - 1 + i)
            if (idi .ne. 1 .and. idi .ne. -1) then
                go to 724
            end if
            if (idi .eq. -1) then
                index = 1
            end if
50          continue
        end do
    end if
    if (tout .eq. t) then
        go to 719
    end if
    if (nrt .lt. 0) then
        go to 730
    end if
    if (info(7) .ne. 0) then
        hmax = rwork(lhmax)
        if (hmax .le. 0.0d0) then
            go to 710
        end if
    end if
    iwork(lnst) = 0
    iwork(lnre) = 0
    iwork(lnje) = 0
    iwork(letf) = 0
    iwork(lncfn) = 0
    iwork(lnni) = 0
    iwork(lnli) = 0
    iwork(lnps) = 0
    iwork(lncfl) = 0
    iwork(lnrte) = 0
    iwork(lkprin) = info(18)
    idid = 1
    go to 200
100 continue
    if (info(1) .eq. 1) then
        go to 110
    end if
    itemp = 1
    if (info(1) .ne. -1) then
        go to 701
    end if
    msg = "DASKR--  THE LAST STEP TERMINATED WITH A NEGATIVE"
    call xerrwd(msg, 49, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    msg = "DASKR--  VALUE (=I1) OF IDID AND NO APPROPRIATE"
    call xerrwd(msg, 47, 202, 0, 1, idid, 0, 0, 0.0d0, 0.0d0)
    msg = "DASKR--  ACTION WAS TAKEN. RUN TERMINATED"
    call xerrwd(msg, 41, 203, 1, 0, 0, 0, 0, 0.0d0, 0.0d0)
    return
110 continue
200 continue
    iwork(lnstl) = iwork(lnst)
    nli0 = iwork(lnli)
    nni0 = iwork(lnni)
    ncfn0 = iwork(lncfn)
    ncfl0 = iwork(lncfl)
    nwarn = 0
    nzflg = 0
    rtoli = rtol(1)
    atoli = atol(1)
    do i = 1, neq
        if (info(2) .eq. 1) then
            rtoli = rtol(i)
        end if
        if (info(2) .eq. 1) then
            atoli = atol(i)
        end if
        if (rtoli .gt. 0.0d0 .or. atoli .gt. 0.0d0) then
            nzflg = 1
        end if
        if (rtoli .lt. 0.0d0) then
            go to 706
        end if
        if (atoli .lt. 0.0d0) then
            go to 707
        end if
210     continue
    end do
    if (nzflg .eq. 0) then
        go to 708
    end if
    iwork(llciwp) = lid + lenid
    lsavr = ldelta
    if (info(12) .ne. 0) then
        lsavr = ldelta + neq
    end if
    le = lsavr + neq
    lwt = le + neq
    lvt = lwt
    if (info(16) .ne. 0) then
        lvt = lwt + neq
    end if
    lphi = lvt + neq
    lr0 = lphi + ncphi*neq
    lr1 = lr0 + nrt
    lrx = lr1 + nrt
    lwm = lrx + nrt
    if (info(1) .eq. 1) then
        go to 400
    end if
300 continue
    tn = t
    idid = 1
    call ddawts(neq, info(2), rtol, atol, y, rwork(lwt), rpar, ipar)
    call dinvwt(neq, rwork(lwt), ier)
    if (ier .ne. 0) then
        go to 713
    end if
    if (info(16) .ne. 0) then
        do i = 1, neq
305         rwork(lvt + i - 1) = max(iwork(lid + i - 1), 0)*rwork(lwt + i - 1)
        end do
    end if
    uround = d1mach(4)
    rwork(lround) = uround
    hmin = 4.0d0*uround*max(abs(t), abs(tout))
    if (info(11) .ne. 0) then
        if (info(17) .eq. 0) then
            rwork(lstol) = uround**.6667d0
        else
            if (rwork(lstol) .le. 0.0d0) then
                go to 725
            end if
        end if
    end if
    rwork(lepcon) = 0.33d0
    floatn = neq
    rwork(lsqrn) = sqrt(floatn)
    rwork(lrsqrn) = 1.d0/rwork(lsqrn)
    tdist = abs(tout - t)
    if (tdist .lt. hmin) then
        go to 714
    end if
    if (info(8) .eq. 0) then
        go to 310
    end if
    h0 = rwork(lh)
    if ((tout - t)*h0 .lt. 0.0d0) then
        go to 711
    end if
    if (h0 .eq. 0.0d0) then
        go to 712
    end if
    go to 320
310 continue
    h0 = 0.001d0*tdist
    ypnorm = ddwnrm(neq, yprime, rwork(lvt), rpar, ipar)
    if (ypnorm .gt. 0.5d0/h0) then
        h0 = 0.5d0/ypnorm
    end if
    h0 = sign(h0, tout - t)
320 if (info(7) .eq. 0) then
        go to 330
    end if
    rh = abs(h0)/rwork(lhmax)
    if (rh .gt. 1.0d0) then
        h0 = h0/rh
    end if
330 if (info(4) .eq. 0) then
        go to 340
    end if
    tstop = rwork(ltstop)
    if ((tstop - t)*h0 .lt. 0.0d0) then
        go to 715
    end if
    if ((t + h0 - tstop)*h0 .gt. 0.0d0) then
        h0 = tstop - t
    end if
    if ((tstop - tout)*h0 .lt. 0.0d0) then
        go to 709
    end if
340 if (info(11) .eq. 0) then
        go to 370
    end if
    nwt = 1
    epconi = rwork(lepin)*rwork(lepcon)
    tscale = 0.0d0
    if (index .eq. 0) then
        tscale = tdist
    end if
350 if (info(12) .eq. 0) then
        lyic = lphi + 2*neq
        lypic = lyic + neq
        lpwk = lypic
    call ddasic(tn, y, yprime, neq, info(11), iwork(lid), res, jac, psol, h0, tscale, rwork(lwt), nwt, idid, rpar, ipar, rwork(lphi), rwork(lsavr), rwork(ldelta), rwork(le), rwork(lyic), rwork(lypic), rwork(lpwk), rwork(lwm), iwork(liwm), rwork(lround), rwork(lepli), rwork(lsqrn), rwork(lrsqrn), epconi, rwork(lstol), info(15), icnflg, iwork(licns), ddasid)
    else
        if (info(12) .eq. 1) then
            lyic = lwm
            lypic = lyic + neq
            lpwk = lypic + neq
        call ddasic(tn, y, yprime, neq, info(11), iwork(lid), res, jac, psol, h0, tscale, rwork(lwt), nwt, idid, rpar, ipar, rwork(lphi), rwork(lsavr), rwork(ldelta), rwork(le), rwork(lyic), rwork(lypic), rwork(lpwk), rwork(lwm), iwork(liwm), rwork(lround), rwork(lepli), rwork(lsqrn), rwork(lrsqrn), epconi, rwork(lstol), info(15), icnflg, iwork(licns), ddasik)
        end if
    end if
    if (idid .lt. 0) then
        go to 600
    end if
    if (nwt .eq. 2) then
        go to 355
    end if
    nwt = 2
    call ddawts(neq, info(2), rtol, atol, y, rwork(lwt), rpar, ipar)
    call dinvwt(neq, rwork(lwt), ier)
    if (ier .ne. 0) then
        go to 713
    end if
    go to 350
355 if (info(14) .eq. 1) then
        idid = 4
        h = h0
        if (info(11) .eq. 1) then
            rwork(lhold) = h0
        end if
        go to 590
    end if
    call ddawts(neq, info(2), rtol, atol, y, rwork(lwt), rpar, ipar)
    call dinvwt(neq, rwork(lwt), ier)
    if (ier .ne. 0) then
        go to 713
    end if
    if (info(16) .ne. 0) then
        do i = 1, neq
357         rwork(lvt + i - 1) = max(iwork(lid + i - 1), 0)*rwork(lwt + i - 1)
        end do
    end if
    if (info(8) .ne. 0) then
        h0 = rwork(lh)
        go to 360
    end if
    h0 = 0.001d0*tdist
    ypnorm = ddwnrm(neq, yprime, rwork(lvt), rpar, ipar)
    if (ypnorm .gt. 0.5d0/h0) then
        h0 = 0.5d0/ypnorm
    end if
    h0 = sign(h0, tout - t)
360 if (info(7) .ne. 0) then
        rh = abs(h0)/rwork(lhmax)
        if (rh .gt. 1.0d0) then
            h0 = h0/rh
        end if
    end if
    if (info(4) .ne. 0) then
        tstop = rwork(ltstop)
        if ((t + h0 - tstop)*h0 .gt. 0.0d0) then
            h0 = tstop - t
        end if
    end if
370 h = h0
    rwork(lh) = h
    itemp = lphi + neq
    do i = 1, neq
        rwork(lphi + i - 1) = y(i)
380     rwork(itemp + i - 1) = h*yprime(i)
    end do
    rwork(lt0) = t
    iwork(lirfnd) = 0
    rwork(lpsi) = h
    rwork(lpsi + 1) = 2.0d0*h
    iwork(lkold) = 1
    if (nrt .eq. 0) then
        go to 390
    end if
call drchek(1, rt, nrt, neq, t, tout, y, yprime, rwork(lphi), rwork(lpsi), iwork(lkold), rwork(lr0), rwork(lr1), rwork(lrx), jroot, irt, rwork(lround), info(3), rwork, iwork, rpar, ipar)
    if (irt .lt. 0) then
        go to 731
    end if
390 go to 500
400 continue
    uround = rwork(lround)
    done = .false.
    tn = rwork(ltn)
    h = rwork(lh)
    if (nrt .eq. 0) then
        go to 405
    end if
call drchek(2, rt, nrt, neq, tn, tout, y, yprime, rwork(lphi), rwork(lpsi), iwork(lkold), rwork(lr0), rwork(lr1), rwork(lrx), jroot, irt, rwork(lround), info(3), rwork, iwork, rpar, ipar)
    if (irt .lt. 0) then
        go to 731
    end if
    if (irt .ne. 1) then
        go to 405
    end if
    iwork(lirfnd) = 1
    idid = 5
    t = rwork(lt0)
    done = .true.
    go to 490
405 continue
    if (info(7) .eq. 0) then
        go to 410
    end if
    rh = abs(h)/rwork(lhmax)
    if (rh .gt. 1.0d0) then
        h = h/rh
    end if
410 continue
    if (t .eq. tout) then
        go to 719
    end if
    if ((t - tout)*h .gt. 0.0d0) then
        go to 711
    end if
    if (info(4) .eq. 1) then
        go to 430
    end if
    if (info(3) .eq. 1) then
        go to 420
    end if
    if ((tn - tout)*h .lt. 0.0d0) then
        go to 490
    end if
    call ddatrp(tn, tout, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
    t = tout
    idid = 3
    done = .true.
    go to 490
420 if ((tn - t)*h .le. 0.0d0) then
        go to 490
    end if
    if ((tn - tout)*h .ge. 0.0d0) then
        go to 425
    end if
    call ddatrp(tn, tn, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
    t = tn
    idid = 1
    done = .true.
    go to 490
425 continue
    call ddatrp(tn, tout, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
    t = tout
    idid = 3
    done = .true.
    go to 490
430 if (info(3) .eq. 1) then
        go to 440
    end if
    tstop = rwork(ltstop)
    if ((tn - tstop)*h .gt. 0.0d0) then
        go to 715
    end if
    if ((tstop - tout)*h .lt. 0.0d0) then
        go to 709
    end if
    if ((tn - tout)*h .lt. 0.0d0) then
        go to 450
    end if
    call ddatrp(tn, tout, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
    t = tout
    idid = 3
    done = .true.
    go to 490
440 tstop = rwork(ltstop)
    if ((tn - tstop)*h .gt. 0.0d0) then
        go to 715
    end if
    if ((tstop - tout)*h .lt. 0.0d0) then
        go to 709
    end if
    if ((tn - t)*h .le. 0.0d0) then
        go to 450
    end if
    if ((tn - tout)*h .ge. 0.0d0) then
        go to 445
    end if
    call ddatrp(tn, tn, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
    t = tn
    idid = 1
    done = .true.
    go to 490
445 continue
    call ddatrp(tn, tout, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
    t = tout
    idid = 3
    done = .true.
    go to 490
450 continue
    if (abs(tn - tstop) .gt. 100.0d0*uround*(abs(tn) + abs(h))) then
        go to 460
    end if
    call ddatrp(tn, tstop, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
    idid = 2
    t = tstop
    done = .true.
    go to 490
460 tnext = tn + h
    if ((tnext - tstop)*h .le. 0.0d0) then
        go to 490
    end if
    h = tstop - tn
    rwork(lh) = h
490 if (done) then
        go to 590
    end if
500 continue
    if (iwork(lnst) - iwork(lnstl) .lt. 500) then
        go to 505
    end if
    idid = -1
    go to 527
505 if (info(12) .eq. 0) then
        go to 510
    end if
    nstd = iwork(lnst) - iwork(lnstl)
    nnid = iwork(lnni) - nni0
    if (nstd .lt. 10 .or. nnid .eq. 0) then
        go to 510
    end if
    avlin = real(iwork(lnli) - nli0)/real(nnid)
    rcfn = real(iwork(lncfn) - ncfn0)/real(nstd)
    rcfl = real(iwork(lncfl) - ncfl0)/real(nnid)
    fmaxl = iwork(lmaxl)
    lavl = avlin .gt. fmaxl
    lcfn = rcfn .gt. 0.9d0
    lcfl = rcfl .gt. 0.9d0
    lwarn = lavl .or. lcfn .or. lcfl
    if (.not. lwarn) then
        go to 510
    end if
    nwarn = nwarn + 1
    if (nwarn .gt. 10) then
        go to 510
    end if
    if (lavl) then
        msg = "DASKR-- Warning. Poor iterative algorithm performance   "
        call xerrwd(msg, 56, 501, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
        msg = "      at T = R1. Average no. of linear iterations = R2  "
        call xerrwd(msg, 56, 501, 0, 0, 0, 0, 2, tn, avlin)
    end if
    if (lcfn) then
        msg = "DASKR-- Warning. Poor iterative algorithm performance   "
        call xerrwd(msg, 56, 502, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
        msg = "      at T = R1. Nonlinear convergence failure rate = R2"
        call xerrwd(msg, 56, 502, 0, 0, 0, 0, 2, tn, rcfn)
    end if
    if (lcfl) then
        msg = "DASKR-- Warning. Poor iterative algorithm performance   "
        call xerrwd(msg, 56, 503, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
        msg = "      at T = R1. Linear convergence failure rate = R2   "
        call xerrwd(msg, 56, 503, 0, 0, 0, 0, 2, tn, rcfl)
    end if
510 call ddawts(neq, info(2), rtol, atol, rwork(lphi), rwork(lwt), rpar, ipar)
    call dinvwt(neq, rwork(lwt), ier)
    if (ier .ne. 0) then
        idid = -3
        go to 527
    end if
    if (info(16) .ne. 0) then
        do i = 1, neq
515         rwork(lvt + i - 1) = max(iwork(lid + i - 1), 0)*rwork(lwt + i - 1)
        end do
    end if
    r = ddwnrm(neq, rwork(lphi), rwork(lwt), rpar, ipar)*100.0d0*uround
    if (r .le. 1.0d0) then
        go to 525
    end if
    if (info(2) .eq. 1) then
        go to 523
    end if
    rtol(1) = r*rtol(1)
    atol(1) = r*atol(1)
    idid = -2
    go to 527
523 do i = 1, neq
        rtol(i) = r*rtol(i)
524     atol(i) = r*atol(i)
    end do
    idid = -2
    go to 527
525 continue
    hmin = 4.0d0*uround*max(abs(tn), abs(tout))
    if (info(7) .ne. 0) then
        rh = abs(h)/rwork(lhmax)
        if (rh .gt. 1.0d0) then
            h = h/rh
        end if
    end if
    if (info(12) .eq. 0) then
    call ddstp(tn, y, yprime, neq, res, jac, psol, h, rwork(lwt), rwork(lvt), info(1), idid, rpar, ipar, rwork(lphi), rwork(lsavr), rwork(ldelta), rwork(le), rwork(lwm), iwork(liwm), rwork(lalpha), rwork(lbeta), rwork(lgamma), rwork(lpsi), rwork(lsigma), rwork(lcj), rwork(lcjold), rwork(lhold), rwork(ls), hmin, rwork(lround), rwork(lepli), rwork(lsqrn), rwork(lrsqrn), rwork(lepcon), iwork(lphase), iwork(ljcalc), info(15), iwork(lk), iwork(lkold), iwork(lns), nonneg, info(12), dnedd)
    else
        if (info(12) .eq. 1) then
        call ddstp(tn, y, yprime, neq, res, jac, psol, h, rwork(lwt), rwork(lvt), info(1), idid, rpar, ipar, rwork(lphi), rwork(lsavr), rwork(ldelta), rwork(le), rwork(lwm), iwork(liwm), rwork(lalpha), rwork(lbeta), rwork(lgamma), rwork(lpsi), rwork(lsigma), rwork(lcj), rwork(lcjold), rwork(lhold), rwork(ls), hmin, rwork(lround), rwork(lepli), rwork(lsqrn), rwork(lrsqrn), rwork(lepcon), iwork(lphase), iwork(ljcalc), info(15), iwork(lk), iwork(lkold), iwork(lns), nonneg, info(12), dnedk)
        end if
    end if
527 if (idid .lt. 0) then
        go to 600
    end if
    if (nrt .eq. 0) then
        go to 530
    end if
call drchek(3, rt, nrt, neq, tn, tout, y, yprime, rwork(lphi), rwork(lpsi), iwork(lkold), rwork(lr0), rwork(lr1), rwork(lrx), jroot, irt, rwork(lround), info(3), rwork, iwork, rpar, ipar)
    if (irt .ne. 1) then
        go to 530
    end if
    iwork(lirfnd) = 1
    idid = 5
    t = rwork(lt0)
    go to 580
530 if (info(4) .eq. 0) then
        if ((tn - tout)*h .ge. 0.0d0) then
            call ddatrp(tn, tout, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
            t = tout
            idid = 3
            go to 580
        end if
        if (info(3) .eq. 0) then
            go to 500
        end if
        t = tn
        idid = 1
        go to 580
    end if
540 if (info(3) .ne. 0) then
        go to 550
    end if
    if (abs(tn - tstop) .le. 100.0d0*uround*(abs(tn) + abs(h))) then
        call ddatrp(tn, tstop, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
        t = tstop
        idid = 2
        go to 580
    end if
    if ((tn - tout)*h .ge. 0.0d0) then
        call ddatrp(tn, tout, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
        t = tout
        idid = 3
        go to 580
    end if
    tnext = tn + h
    if ((tnext - tstop)*h .le. 0.0d0) then
        go to 500
    end if
    h = tstop - tn
    go to 500
550 continue
    if (abs(tn - tstop) .le. 100.0d0*uround*(abs(tn) + abs(h))) then
        call ddatrp(tn, tstop, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
        t = tstop
        idid = 2
        go to 580
    end if
    if ((tn - tout)*h .ge. 0.0d0) then
        call ddatrp(tn, tout, y, yprime, neq, iwork(lkold), rwork(lphi), rwork(lpsi))
        t = tout
        idid = 3
        go to 580
    end if
    t = tn
    idid = 1
580 continue
590 continue
    rwork(ltn) = tn
    rwork(ltlast) = t
    rwork(lh) = h
    return
600 continue
    itemp = -idid
    go to(610, 620, 630, 700, 655, 640, 650, 660, 670, 675, 680, 685, 690, 695), itemp
610 msg = "DASKR--  AT CURRENT T (=R1)  500 STEPS"
    call xerrwd(msg, 38, 610, 0, 0, 0, 0, 1, tn, 0.0d0)
    msg = "DASKR--  TAKEN ON THIS CALL BEFORE REACHING TOUT"
    call xerrwd(msg, 48, 611, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
620 msg = "DASKR--  AT T (=R1) TOO MUCH ACCURACY REQUESTED"
    call xerrwd(msg, 47, 620, 0, 0, 0, 0, 1, tn, 0.0d0)
    msg = "DASKR--  FOR PRECISION OF MACHINE. RTOL AND ATOL"
    call xerrwd(msg, 48, 621, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    msg = "DASKR--  WERE INCREASED BY A FACTOR R (=R1)"
    call xerrwd(msg, 43, 622, 0, 0, 0, 0, 1, r, 0.0d0)
    go to 700
630 msg = "DASKR--  AT T (=R1) SOME ELEMENT OF WT"
    call xerrwd(msg, 38, 630, 0, 0, 0, 0, 1, tn, 0.0d0)
    msg = "DASKR--  HAS BECOME .LE. 0.0"
    call xerrwd(msg, 28, 631, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
640 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 640, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN"
    call xerrwd(msg, 57, 641, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
650 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 650, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  NONLINEAR SOLVER FAILED TO CONVERGE"
    call xerrwd(msg, 44, 651, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    msg = "DASKR--  REPEATEDLY OR WITH ABS(H)=HMIN"
    call xerrwd(msg, 40, 652, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
655 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 655, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  PRECONDITIONER HAD REPEATED FAILURES."
    call xerrwd(msg, 46, 656, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
660 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 660, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  ITERATION MATRIX IS SINGULAR."
    call xerrwd(msg, 38, 661, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
670 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 670, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  NONLINEAR SOLVER COULD NOT CONVERGE."
    call xerrwd(msg, 45, 671, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    msg = "DASKR--  ALSO, THE ERROR TEST FAILED REPEATEDLY."
    call xerrwd(msg, 49, 672, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
675 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 675, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  NONLINEAR SYSTEM SOLVER COULD NOT CONVERGE"
    call xerrwd(msg, 51, 676, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    msg = "DASKR--  BECAUSE IRES WAS EQUAL TO MINUS ONE"
    call xerrwd(msg, 44, 677, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
680 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2)"
    call xerrwd(msg, 40, 680, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  IRES WAS EQUAL TO MINUS TWO"
    call xerrwd(msg, 36, 681, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
685 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 685, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    msg = "DASKR--  INITIAL (Y,YPRIME) COULD NOT BE COMPUTED"
    call xerrwd(msg, 49, 686, 0, 0, 0, 0, 2, tn, h0)
    go to 700
690 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2)"
    call xerrwd(msg, 40, 690, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  IER WAS NEGATIVE FROM PSOL"
    call xerrwd(msg, 35, 691, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
695 msg = "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE"
    call xerrwd(msg, 44, 695, 0, 0, 0, 0, 2, tn, h)
    msg = "DASKR--  LINEAR SYSTEM SOLVER COULD NOT CONVERGE."
    call xerrwd(msg, 50, 696, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 700
700 continue
    info(1) = -1
    t = tn
    rwork(ltn) = tn
    rwork(lh) = h
    return
701 msg = "DASKR--  ELEMENT (=I1) OF INFO VECTOR IS NOT VALID"
    call xerrwd(msg, 50, 1, 0, 1, itemp, 0, 0, 0.0d0, 0.0d0)
    go to 750
702 msg = "DASKR--  NEQ (=I1) .LE. 0"
    call xerrwd(msg, 25, 2, 0, 1, neq, 0, 0, 0.0d0, 0.0d0)
    go to 750
703 msg = "DASKR--  MAXORD (=I1) NOT IN RANGE"
    call xerrwd(msg, 34, 3, 0, 1, mxord, 0, 0, 0.0d0, 0.0d0)
    go to 750
704 msg = "DASKR--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)"
    call xerrwd(msg, 60, 4, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
    go to 750
705 msg = "DASKR--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)"
    call xerrwd(msg, 60, 5, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
    go to 750
706 msg = "DASKR--  SOME ELEMENT OF RTOL IS .LT. 0"
    call xerrwd(msg, 39, 6, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
707 msg = "DASKR--  SOME ELEMENT OF ATOL IS .LT. 0"
    call xerrwd(msg, 39, 7, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
708 msg = "DASKR--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO"
    call xerrwd(msg, 47, 8, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
709 msg = "DASKR--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2)"
    call xerrwd(msg, 54, 9, 0, 0, 0, 0, 2, tstop, tout)
    go to 750
710 msg = "DASKR--  HMAX (=R1) .LT. 0.0"
    call xerrwd(msg, 28, 10, 0, 0, 0, 0, 1, hmax, 0.0d0)
    go to 750
711 msg = "DASKR--  TOUT (=R1) BEHIND T (=R2)"
    call xerrwd(msg, 34, 11, 0, 0, 0, 0, 2, tout, t)
    go to 750
712 msg = "DASKR--  INFO(8)=1 AND H0=0.0"
    call xerrwd(msg, 29, 12, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
713 msg = "DASKR--  SOME ELEMENT OF WT IS .LE. 0.0"
    call xerrwd(msg, 39, 13, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
714 msg = "DASKR-- TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION"
    call xerrwd(msg, 60, 14, 0, 0, 0, 0, 2, tout, t)
    go to 750
715 msg = "DASKR--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2)"
    call xerrwd(msg, 49, 15, 0, 0, 0, 0, 2, tstop, t)
    go to 750
717 msg = "DASKR--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ"
    call xerrwd(msg, 52, 17, 0, 1, iwork(lml), 0, 0, 0.0d0, 0.0d0)
    go to 750
718 msg = "DASKR--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ"
    call xerrwd(msg, 52, 18, 0, 1, iwork(lmu), 0, 0, 0.0d0, 0.0d0)
    go to 750
719 msg = "DASKR--  TOUT (=R1) IS EQUAL TO T (=R2)"
    call xerrwd(msg, 39, 19, 0, 0, 0, 0, 2, tout, t)
    go to 750
720 msg = "DASKR--  MAXL (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. NEQ"
    call xerrwd(msg, 54, 20, 0, 1, iwork(lmaxl), 0, 0, 0.0d0, 0.0d0)
    go to 750
721 msg = "DASKR--  KMP (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. MAXL"
    call xerrwd(msg, 54, 21, 0, 1, iwork(lkmp), 0, 0, 0.0d0, 0.0d0)
    go to 750
722 msg = "DASKR--  NRMAX (=I1) ILLEGAL. .LT. 0"
    call xerrwd(msg, 36, 22, 0, 1, iwork(lnrmax), 0, 0, 0.0d0, 0.0d0)
    go to 750
723 msg = "DASKR--  EPLI (=R1) ILLEGAL. EITHER .LE. 0.D0 OR .GE. 1.D0"
    call xerrwd(msg, 58, 23, 0, 0, 0, 0, 1, rwork(lepli), 0.0d0)
    go to 750
724 msg = "DASKR--  ILLEGAL IWORK VALUE FOR INFO(11) .NE. 0"
    call xerrwd(msg, 48, 24, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
725 msg = "DASKR--  ONE OF THE INPUTS FOR INFO(17) = 1 IS ILLEGAL"
    call xerrwd(msg, 54, 25, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
726 msg = "DASKR--  ILLEGAL IWORK VALUE FOR INFO(10) .NE. 0"
    call xerrwd(msg, 48, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
    go to 750
727 msg = "DASKR--  Y(I) AND IWORK(40+I) (I=I1) INCONSISTENT"
    call xerrwd(msg, 49, 27, 0, 1, iret, 0, 0, 0.0d0, 0.0d0)
    go to 750
730 msg = "DASKR--  NRT (=I1) .LT. 0"
    call xerrwd(msg, 25, 30, 1, 1, nrt, 0, 0, 0.0d0, 0.0d0)
    go to 750
731 msg = "DASKR--  R IS ILL-DEFINED.  ZERO VALUES WERE FOUND AT TWO"
    call xerrwd(msg, 57, 31, 1, 0, 0, 0, 0, 0.0d0, 0.0d0)
    msg = "         VERY CLOSE T VALUES, AT T = R1"
    call xerrwd(msg, 39, 31, 1, 0, 0, 0, 1, rwork(lt0), 0.0d0)
750 if (info(1) .eq. -1) then
        go to 760
    end if
    info(1) = -1
    idid = -33
    return
760 msg = "DASKR--  REPEATED OCCURRENCES OF ILLEGAL INPUT"
    call xerrwd(msg, 46, 701, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
770 msg = "DASKR--  RUN TERMINATED. APPARENT INFINITE LOOP"
    call xerrwd(msg, 47, 702, 1, 0, 0, 0, 0, 0.0d0, 0.0d0)
    return
end subroutine ddaskr

subroutine drchek(job, rt, nrt, neq, tn, tout, y, yp, phi, psi, kold, r0, r1, rx, jroot, irt, uround, info3, rwork, iwork, rpar, ipar)
    implicit double precision(a - h, o - z)
    parameter(lnrte=36, lirfnd=37)
    parameter(lt0=51, ltlast=52)
    external :: rt
    integer :: job, nrt, neq, kold, jroot, irt, info3, iwork, ipar
    double precision :: tn, tout, y, yp, phi, psi, r0, r1, rx, uround, rwork, rpar
    dimension :: y(*), yp(*), phi(neq, *), psi(*), r0(*), r1(*), rx(*), jroot(*), rwork(*), iwork(*), rpar(*), ipar(*)
    integer :: i, jflag
    double precision :: h
    double precision :: hminr, t1, temp1, temp2, x, zero
    logical :: zroot
    data zero/0.0d0/
    h = psi(1)
    irt = 0
    do i = 1, nrt
10      jroot(i) = 0
    end do
    hminr = (abs(tn) + abs(h))*uround*100.0d0
    go to(100, 200, 300), job
100 continue
    call ddatrp(tn, rwork(lt0), y, yp, neq, kold, phi, psi)
    call rt(neq, rwork(lt0), y, yp, nrt, r0, rpar, ipar)
    iwork(lnrte) = 1
    zroot = .false.
    do i = 1, nrt
110     if (abs(r0(i)) .eq. zero) then
            zroot = .true.
        end if
    end do
    if (.not. zroot) then
        go to 190
    end if
    temp2 = max(hminr/abs(h), 0.1d0)
    temp1 = temp2*h
    rwork(lt0) = rwork(lt0) + temp1
    do i = 1, neq
120     y(i) = y(i) + temp2*phi(i, 2)
    end do
    call rt(neq, rwork(lt0), y, yp, nrt, r0, rpar, ipar)
    iwork(lnrte) = iwork(lnrte) + 1
    zroot = .false.
    do i = 1, nrt
130     if (abs(r0(i)) .eq. zero) then
            zroot = .true.
        end if
    end do
    if (.not. zroot) then
        go to 190
    end if
    irt = -1
    return
190 continue
    return
200 continue
    if (iwork(lirfnd) .eq. 0) then
        go to 260
    end if
    call ddatrp(tn, rwork(lt0), y, yp, neq, kold, phi, psi)
    call rt(neq, rwork(lt0), y, yp, nrt, r0, rpar, ipar)
    iwork(lnrte) = iwork(lnrte) + 1
    zroot = .false.
    do i = 1, nrt
        if (abs(r0(i)) .eq. zero) then
            zroot = .true.
            jroot(i) = 1
        end if
210     continue
    end do
    if (.not. zroot) then
        go to 260
    end if
    temp1 = sign(hminr, h)
    rwork(lt0) = rwork(lt0) + temp1
    if ((rwork(lt0) - tn)*h .lt. zero) then
        go to 230
    end if
    temp2 = temp1/h
    do i = 1, neq
220     y(i) = y(i) + temp2*phi(i, 2)
    end do
    go to 240
230 call ddatrp(tn, rwork(lt0), y, yp, neq, kold, phi, psi)
240 call rt(neq, rwork(lt0), y, yp, nrt, r0, rpar, ipar)
    iwork(lnrte) = iwork(lnrte) + 1
    do i = 1, nrt
        if (abs(r0(i)) .gt. zero) then
            go to 250
        end if
        if (jroot(i) .eq. 1) then
            irt = -2
            return
        else
            jroot(i) = -sign(1.0d0, r0(i))
            irt = 1
        end if
250     continue
    end do
    if (irt .eq. 1) then
        return
    end if
260 if (tn .eq. rwork(ltlast)) then
        return
    end if
300 continue
    if ((tout - tn)*h .ge. zero) then
        t1 = tn
        go to 330
    end if
    t1 = tout
    if ((t1 - rwork(lt0))*h .le. zero) then
        go to 390
    end if
330 call ddatrp(tn, t1, y, yp, neq, kold, phi, psi)
    call rt(neq, t1, y, yp, nrt, r1, rpar, ipar)
    iwork(lnrte) = iwork(lnrte) + 1
    jflag = 0
350 continue
    call droots(nrt, hminr, jflag, rwork(lt0), t1, r0, r1, rx, x, jroot)
    if (jflag .gt. 1) then
        go to 360
    end if
    call ddatrp(tn, x, y, yp, neq, kold, phi, psi)
    call rt(neq, x, y, yp, nrt, rx, rpar, ipar)
    iwork(lnrte) = iwork(lnrte) + 1
    go to 350
360 rwork(lt0) = x
    call dcopy(nrt, rx, 1, r0, 1)
    if (jflag .eq. 4) then
        go to 390
    end if
    call ddatrp(tn, x, y, yp, neq, kold, phi, psi)
    irt = 1
    return
390 continue
    return
end subroutine drchek

subroutine droots(nrt, hmin, jflag, x0, x1, r0, r1, rx, x, jroot)
    integer :: nrt, jflag, jroot
    double precision :: hmin, x0, x1, r0, r1, rx, x
    dimension :: r0(nrt), r1(nrt), rx(nrt), jroot(nrt)
    integer :: i, imax, imxold, last, nxlast
    double precision :: alpha, t2, tmax, x2, fracint, fracsub, zero, tenth, half, five
    logical :: zroot, sgnchg, xroot
    save :: alpha, x2, imax, last
    data zero/0.0d0/, tenth/0.1d0/, half/0.5d0/, five/5.0d0/
    if (jflag .eq. 1) then
        go to 200
    end if
    imax = 0
    tmax = zero
    zroot = .false.
    do i = 1, nrt
        if (abs(r1(i)) .gt. zero) then
            go to 110
        end if
        zroot = .true.
        go to 120
110     if (sign(1.0d0, r0(i)) .eq. sign(1.0d0, r1(i))) then
            go to 120
        end if
        t2 = abs(r1(i)/(r1(i) - r0(i)))
        if (t2 .le. tmax) then
            go to 120
        end if
        tmax = t2
        imax = i
120     continue
    end do
    if (imax .gt. 0) then
        go to 130
    end if
    sgnchg = .false.
    go to 140
130 sgnchg = .true.
140 if (.not. sgnchg) then
        go to 400
    end if
    xroot = .false.
    nxlast = 0
    last = 1
150 continue
    if (xroot) then
        go to 300
    end if
    if (nxlast .eq. last) then
        go to 160
    end if
    alpha = 1.0d0
    go to 180
160 if (last .eq. 0) then
        go to 170
    end if
    alpha = 0.5d0*alpha
    go to 180
170 alpha = 2.0d0*alpha
180 x2 = x1 - (x1 - x0)*r1(imax)/(r1(imax) - alpha*r0(imax))
    if (abs(x2 - x0) .lt. half*hmin) then
        fracint = abs(x1 - x0)/hmin
        if (fracint .gt. five) then
            fracsub = tenth
        else
            fracsub = half/fracint
        end if
        x2 = x0 + fracsub*(x1 - x0)
    end if
    if (abs(x1 - x2) .lt. half*hmin) then
        fracint = abs(x1 - x0)/hmin
        if (fracint .gt. five) then
            fracsub = tenth
        else
            fracsub = half/fracint
        end if
        x2 = x1 - fracsub*(x1 - x0)
    end if
    jflag = 1
    x = x2
    return
200 imxold = imax
    imax = 0
    tmax = zero
    zroot = .false.
    do i = 1, nrt
        if (abs(rx(i)) .gt. zero) then
            go to 210
        end if
        zroot = .true.
        go to 220
210     if (sign(1.0d0, r0(i)) .eq. sign(1.0d0, rx(i))) then
            go to 220
        end if
        t2 = abs(rx(i)/(rx(i) - r0(i)))
        if (t2 .le. tmax) then
            go to 220
        end if
        tmax = t2
        imax = i
220     continue
    end do
    if (imax .gt. 0) then
        go to 230
    end if
    sgnchg = .false.
    imax = imxold
    go to 240
230 sgnchg = .true.
240 nxlast = last
    if (.not. sgnchg) then
        go to 250
    end if
    x1 = x2
    call dcopy(nrt, rx, 1, r1, 1)
    last = 1
    xroot = .false.
    go to 270
250 if (.not. zroot) then
        go to 260
    end if
    x1 = x2
    call dcopy(nrt, rx, 1, r1, 1)
    xroot = .true.
    go to 270
260 continue
    call dcopy(nrt, rx, 1, r0, 1)
    x0 = x2
    last = 0
    xroot = .false.
270 if (abs(x1 - x0) .le. hmin) then
        xroot = .true.
    end if
    go to 150
300 jflag = 2
    x = x1
    call dcopy(nrt, r1, 1, rx, 1)
    do i = 1, nrt
        jroot(i) = 0
        if (abs(r1(i)) .eq. zero) then
            jroot(i) = -sign(1.0d0, r0(i))
            go to 320
        end if
        if (sign(1.0d0, r0(i)) .ne. sign(1.0d0, r1(i))) then
            jroot(i) = sign(1.0d0, r1(i) - r0(i))
        end if
320     continue
    end do
    return
400 if (.not. zroot) then
        go to 420
    end if
    x = x1
    call dcopy(nrt, r1, 1, rx, 1)
    do i = 1, nrt
        jroot(i) = 0
        if (abs(r1(i)) .eq. zero) then
            jroot(i) = -sign(1.0d0, r0(i))
        end if
410     continue
    end do
    jflag = 3
    return
420 call dcopy(nrt, r1, 1, rx, 1)
    x = x1
    jflag = 4
    return
end subroutine droots

subroutine ddasic(x, y, yprime, neq, icopt, id, res, jac, psol, h, tscale, wt, nic, idid, rpar, ipar, phi, savr, delta, e, yic, ypic, pwk, wm, iwm, uround, epli, sqrtn, rsqrtn, epconi, stptol, jflg, icnflg, icnstr, nlsic)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), id(*), wt(*), phi(neq, *)
    dimension :: savr(*), delta(*), e(*), yic(*), ypic(*), pwk(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*), icnstr(*)
    external :: res, jac, psol, nlsic
    parameter(lcfn=15)
    parameter(lmxnh=34)
    save :: rhcut, ratemx
    data rhcut/0.1d0/, ratemx/0.8d0/
    mxnh = iwm(lmxnh)
    idid = 1
    nh = 1
    jskip = 0
    if (nic .eq. 2) then
        jskip = 1
    end if
    call dcopy(neq, y, 1, phi(1, 1), 1)
    call dcopy(neq, yprime, 1, phi(1, 2), 1)
    if (icopt .eq. 2) then
        cj = 0.0d0
    else
        cj = 1.0d0/h
    end if
200 continue
call nlsic(x, y, yprime, neq, icopt, id, res, jac, psol, h, tscale, wt, jskip, rpar, ipar, savr, delta, e, yic, ypic, pwk, wm, iwm, cj, uround, epli, sqrtn, rsqrtn, epconi, ratemx, stptol, jflg, icnflg, icnstr, iernls)
    if (iernls .eq. 0) then
        return
    end if
    iwm(lcfn) = iwm(lcfn) + 1
    jskip = 0
    if (iernls .eq. -1) then
        go to 350
    end if
    if (icopt .eq. 2) then
        go to 350
    end if
    if (nh .eq. mxnh) then
        go to 350
    end if
    nh = nh + 1
    h = h*rhcut
    cj = 1.0d0/h
    if (iernls .eq. 1) then
        go to 200
    end if
    call dcopy(neq, phi(1, 1), 1, y, 1)
    call dcopy(neq, phi(1, 2), 1, yprime, 1)
    go to 200
350 idid = -12
    return
end subroutine ddasic

subroutine dyypnw(neq, y, yprime, cj, rl, p, icopt, id, ynew, ypnew)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), ynew(*), ypnew(*), id(*), p(*)
    if (icopt .eq. 1) then
        do i = 1, neq
            if (id(i) .lt. 0) then
                ynew(i) = y(i) - rl*p(i)
                ypnew(i) = yprime(i)
            else
                ynew(i) = y(i)
                ypnew(i) = yprime(i) - rl*cj*p(i)
            end if
10          continue
        end do
    else
        do i = 1, neq
            ynew(i) = y(i) - rl*p(i)
            ypnew(i) = yprime(i)
20          continue
        end do
    end if
    return
end subroutine dyypnw

subroutine ddstp(x, y, yprime, neq, res, jac, psol, h, wt, vt, jstart, idid, rpar, ipar, phi, savr, delta, e, wm, iwm, alpha, beta, gamma, psi, sigma, cj, cjold, hold, s, hmin, uround, epli, sqrtn, rsqrtn, epcon, iphase, jcalc, jflg, k, kold, ns, nonneg, ntype, nls)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), wt(*), vt(*)
    dimension :: phi(neq, *), savr(*), delta(*), e(*)
    dimension :: wm(*), iwm(*)
    dimension :: psi(*), alpha(*), beta(*), gamma(*), sigma(*)
    dimension :: rpar(*), ipar(*)
    external :: res, jac, psol, nls
    parameter(lmxord=3)
    parameter(lnst=11, letf=14, lcfn=15)
    xold = x
    ncf = 0
    nef = 0
    if (jstart .ne. 0) then
        go to 120
    end if
    k = 1
    kold = 0
    hold = 0.0d0
    psi(1) = h
    cj = 1.d0/h
    iphase = 0
    ns = 0
120 continue
200 continue
    kp1 = k + 1
    kp2 = k + 2
    km1 = k - 1
    if (h .ne. hold .or. k .ne. kold) then
        ns = 0
    end if
    ns = min0(ns + 1, kold + 2)
    nsp1 = ns + 1
    if (kp1 .lt. ns) then
        go to 230
    end if
    beta(1) = 1.0d0
    alpha(1) = 1.0d0
    temp1 = h
    gamma(1) = 0.0d0
    sigma(1) = 1.0d0
    do i = 2, kp1
        temp2 = psi(i - 1)
        psi(i - 1) = temp1
        beta(i) = beta(i - 1)*psi(i - 1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        sigma(i) = (i - 1)*sigma(i - 1)*alpha(i)
        gamma(i) = gamma(i - 1) + alpha(i - 1)/h
210     continue
    end do
    psi(kp1) = temp1
230 continue
    alphas = 0.0d0
    alpha0 = 0.0d0
    do i = 1, k
        alphas = alphas - 1.0d0/i
        alpha0 = alpha0 - alpha(i)
240     continue
    end do
    cjlast = cj
    cj = (-alphas)/h
    ck = abs(alpha(kp1) + alphas - alpha0)
    ck = max(ck, alpha(kp1))
    if (kp1 .lt. nsp1) then
        go to 280
    end if
    do j = nsp1, kp1
        do i = 1, neq
260         phi(i, j) = beta(j)*phi(i, j)
        end do
270     continue
    end do
280 continue
    x = x + h
    idid = 1
call nls(x, y, yprime, neq, res, jac, psol, h, wt, jstart, idid, rpar, ipar, phi, gamma, savr, delta, e, wm, iwm, cj, cjold, cjlast, s, uround, epli, sqrtn, rsqrtn, epcon, jcalc, jflg, kp1, nonneg, ntype, iernls)
    if (iernls .ne. 0) then
        go to 600
    end if
    enorm = ddwnrm(neq, e, vt, rpar, ipar)
    erk = sigma(k + 1)*enorm
    terk = (k + 1)*erk
    est = erk
    knew = k
    if (k .eq. 1) then
        go to 430
    end if
    do i = 1, neq
405     delta(i) = phi(i, kp1) + e(i)
    end do
    erkm1 = sigma(k)*ddwnrm(neq, delta, vt, rpar, ipar)
    terkm1 = k*erkm1
    if (k .gt. 2) then
        go to 410
    end if
    if (terkm1 .le. 0.5*terk) then
        go to 420
    end if
    go to 430
410 continue
    do i = 1, neq
415     delta(i) = phi(i, k) + delta(i)
    end do
    erkm2 = sigma(k - 1)*ddwnrm(neq, delta, vt, rpar, ipar)
    terkm2 = (k - 1)*erkm2
    if (max(terkm1, terkm2) .gt. terk) then
        go to 430
    end if
420 continue
    knew = k - 1
    est = erkm1
430 continue
    err = ck*enorm
    if (err .gt. 1.0d0) then
        go to 600
    end if
    idid = 1
    iwm(lnst) = iwm(lnst) + 1
    kdiff = k - kold
    kold = k
    hold = h
    if (knew .eq. km1 .or. k .eq. iwm(lmxord)) then
        iphase = 1
    end if
    if (iphase .eq. 0) then
        go to 545
    end if
    if (knew .eq. km1) then
        go to 540
    end if
    if (k .eq. iwm(lmxord)) then
        go to 550
    end if
    if (kp1 .ge. ns .or. kdiff .eq. 1) then
        go to 550
    end if
    do i = 1, neq
510     delta(i) = e(i) - phi(i, kp2)
    end do
    erkp1 = 1.0d0/(k + 2)*ddwnrm(neq, delta, vt, rpar, ipar)
    terkp1 = (k + 2)*erkp1
    if (k .gt. 1) then
        go to 520
    end if
    if (terkp1 .ge. 0.5d0*terk) then
        go to 550
    end if
    go to 530
520 if (terkm1 .le. min(terk, terkp1)) then
        go to 540
    end if
    if (terkp1 .ge. terk .or. k .eq. iwm(lmxord)) then
        go to 550
    end if
530 k = kp1
    est = erkp1
    go to 550
540 k = km1
    est = erkm1
    go to 550
545 k = kp1
    hnew = h*2.0d0
    h = hnew
    go to 575
550 hnew = h
    temp2 = k + 1
    r = (2.0d0*est + 0.0001d0)**((-1.0d0)/temp2)
    if (r .lt. 2.0d0) then
        go to 555
    end if
    hnew = 2.0d0*h
    go to 560
555 if (r .gt. 1.0d0) then
        go to 560
    end if
    r = max(0.5d0, min(0.9d0, r))
    hnew = h*r
560 h = hnew
575 continue
    if (kold .eq. iwm(lmxord)) then
        go to 585
    end if
    do i = 1, neq
580     phi(i, kp2) = e(i)
    end do
585 continue
    do i = 1, neq
590     phi(i, kp1) = phi(i, kp1) + e(i)
    end do
    do j1 = 2, kp1
        j = kp1 - j1 + 1
        do i = 1, neq
595         phi(i, j) = phi(i, j) + phi(i, j + 1)
        end do
    end do
    jstart = 1
    return
600 iphase = 1
    x = xold
    if (kp1 .lt. nsp1) then
        go to 630
    end if
    do j = nsp1, kp1
        temp1 = 1.0d0/beta(j)
        do i = 1, neq
610         phi(i, j) = temp1*phi(i, j)
        end do
620     continue
    end do
630 continue
    do i = 2, kp1
640     psi(i - 1) = psi(i) - h
    end do
    if (iernls .eq. 0) then
        go to 660
    end if
    iwm(lcfn) = iwm(lcfn) + 1
    if (iernls .lt. 0) then
        go to 675
    end if
    ncf = ncf + 1
    r = 0.25d0
    h = h*r
    if (ncf .lt. 10 .and. abs(h) .ge. hmin) then
        go to 690
    end if
    if (idid .eq. 1) then
        idid = -7
    end if
    if (nef .ge. 3) then
        idid = -9
    end if
    go to 675
660 nef = nef + 1
    iwm(letf) = iwm(letf) + 1
    if (nef .gt. 1) then
        go to 665
    end if
    k = knew
    temp2 = k + 1
    r = 0.90d0*(2.0d0*est + 0.0001d0)**((-1.0d0)/temp2)
    r = max(0.25d0, min(0.9d0, r))
    h = h*r
    if (abs(h) .ge. hmin) then
        go to 690
    end if
    idid = -6
    go to 675
665 if (nef .gt. 2) then
        go to 670
    end if
    k = knew
    r = 0.25d0
    h = r*h
    if (abs(h) .ge. hmin) then
        go to 690
    end if
    idid = -6
    go to 675
670 k = 1
    r = 0.25d0
    h = r*h
    if (abs(h) .ge. hmin) then
        go to 690
    end if
    idid = -6
    go to 675
675 continue
    call ddatrp(x, x, y, yprime, neq, k, phi, psi)
    jstart = 1
    if (idid .ge. 0) then
        idid = -7
    end if
    return
690 if (kold .eq. 0) then
        psi(1) = h
        do i = 1, neq
695         phi(i, 2) = r*phi(i, 2)
        end do
    end if
    go to 200
end subroutine ddstp

subroutine dcnstr(neq, y, ynew, icnstr, tau, rlx, iret, ivar)
    implicit double precision(a - h, o - z)
    dimension :: y(neq), ynew(neq), icnstr(neq)
    save :: fac, fac2, zero
    data fac/0.6d0/, fac2/0.9d0/, zero/0.0d0/
    iret = 0
    rdymx = zero
    ivar = 0
    do i = 1, neq
        if (icnstr(i) .eq. 2) then
            rdy = abs((ynew(i) - y(i))/y(i))
            if (rdy .gt. rdymx) then
                rdymx = rdy
                ivar = i
            end if
            if (ynew(i) .le. zero) then
                tau = fac*tau
                ivar = i
                iret = 1
                return
            end if
        else
            if (icnstr(i) .eq. 1) then
                if (ynew(i) .lt. zero) then
                    tau = fac*tau
                    ivar = i
                    iret = 1
                    return
                end if
            else
                if (icnstr(i) .eq. -1) then
                    if (ynew(i) .gt. zero) then
                        tau = fac*tau
                        ivar = i
                        iret = 1
                        return
                    end if
                else
                    if (icnstr(i) .eq. -2) then
                        rdy = abs((ynew(i) - y(i))/y(i))
                        if (rdy .gt. rdymx) then
                            rdymx = rdy
                            ivar = i
                        end if
                        if (ynew(i) .ge. zero) then
                            tau = fac*tau
                            ivar = i
                            iret = 1
                            return
                        end if
                    end if
                end if
            end if
        end if
100     continue
    end do
    if (rdymx .ge. rlx) then
        tau = fac2*tau*rlx/rdymx
        iret = 1
    end if
    return
end subroutine dcnstr

subroutine dcnst0(neq, y, icnstr, iret)
    implicit double precision(a - h, o - z)
    dimension :: y(neq), icnstr(neq)
    save :: zero
    data zero/0.d0/
    iret = 0
    do i = 1, neq
        if (icnstr(i) .eq. 2) then
            if (y(i) .le. zero) then
                iret = i
                return
            end if
        else
            if (icnstr(i) .eq. 1) then
                if (y(i) .lt. zero) then
                    iret = i
                    return
                end if
            else
                if (icnstr(i) .eq. -1) then
                    if (y(i) .gt. zero) then
                        iret = i
                        return
                    end if
                else
                    if (icnstr(i) .eq. -2) then
                        if (y(i) .ge. zero) then
                            iret = i
                            return
                        end if
                    end if
                end if
            end if
        end if
100     continue
    end do
    return
end subroutine dcnst0

subroutine ddawts(neq, iwt, rtol, atol, y, wt, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: rtol(*), atol(*), y(*), wt(*)
    dimension :: rpar(*), ipar(*)
    rtoli = rtol(1)
    atoli = atol(1)
    do i = 1, neq
        if (iwt .eq. 0) then
            go to 10
        end if
        rtoli = rtol(i)
        atoli = atol(i)
10      wt(i) = rtoli*abs(y(i)) + atoli
20      continue
    end do
    return
end subroutine ddawts

subroutine dinvwt(neq, wt, ier)
    implicit double precision(a - h, o - z)
    dimension :: wt(*)
    do i = 1, neq
        if (wt(i) .le. 0.0d0) then
            go to 30
        end if
10      continue
    end do
    do i = 1, neq
20      wt(i) = 1.0d0/wt(i)
    end do
    ier = 0
    return
30  ier = i
    return
end subroutine dinvwt

subroutine ddatrp(x, xout, yout, ypout, neq, kold, phi, psi)
    implicit double precision(a - h, o - z)
    dimension :: yout(*), ypout(*)
    dimension :: phi(neq, *), psi(*)
    koldp1 = kold + 1
    temp1 = xout - x
    do i = 1, neq
        yout(i) = phi(i, 1)
10      ypout(i) = 0.0d0
    end do
    c = 1.0d0
    d = 0.0d0
    gamma = temp1/psi(1)
    do j = 2, koldp1
        d = d*gamma + c/psi(j - 1)
        c = c*gamma
        gamma = (temp1 + psi(j - 1))/psi(j)
        do i = 1, neq
            yout(i) = yout(i) + c*phi(i, j)
20          ypout(i) = ypout(i) + d*phi(i, j)
        end do
30      continue
    end do
    return
end subroutine ddatrp

double precision function ddwnrm(neq, v, rwt, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: v(*), rwt(*)
    dimension :: rpar(*), ipar(*)
    ddwnrm = 0.0d0
    vmax = 0.0d0
    do i = 1, neq
        if (abs(v(i)*rwt(i)) .gt. vmax) then
            vmax = abs(v(i)*rwt(i))
        end if
10      continue
    end do
    if (vmax .le. 0.0d0) then
        go to 30
    end if
    sum = 0.0d0
    do i = 1, neq
20      sum = sum + (v(i)*rwt(i)/vmax)**2
    end do
    ddwnrm = vmax*sqrt(sum/neq)
30  continue
    return
end function ddwnrm

subroutine ddasid(x, y, yprime, neq, icopt, id, res, jacd, pdum, h, tscale, wt, jsdum, rpar, ipar, dumsvr, delta, r, yic, ypic, dumpwk, wm, iwm, cj, uround, dume, dums, dumr, epcon, ratemx, stptol, jfdum, icnflg, icnstr, iernls)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), id(*), wt(*), icnstr(*)
    dimension :: delta(*), r(*), yic(*), ypic(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    external :: res, jacd
    parameter(lnre=12, lnje=13, lmxnit=32, lmxnj=33)
    mxnit = iwm(lmxnit)
    mxnj = iwm(lmxnj)
    iernls = 0
    nj = 0
    ires = 0
    iwm(lnre) = iwm(lnre) + 1
    call res(x, y, yprime, cj, delta, ires, rpar, ipar)
    if (ires .lt. 0) then
        go to 370
    end if
300 continue
    ierj = 0
    ires = 0
    iernew = 0
    nj = nj + 1
    iwm(lnje) = iwm(lnje) + 1
    call dmatd(neq, x, y, yprime, delta, cj, h, ierj, wt, r, wm, iwm, res, ires, uround, jacd, rpar, ipar)
    if (ires .lt. 0 .or. ierj .ne. 0) then
        go to 370
    end if
call dnsid(x, y, yprime, neq, icopt, id, res, wt, rpar, ipar, delta, r, yic, ypic, wm, iwm, cj, tscale, epcon, ratemx, mxnit, stptol, icnflg, icnstr, iernew)
    if (iernew .eq. 1 .and. nj .lt. mxnj) then
        iwm(lnre) = iwm(lnre) + 1
        call res(x, y, yprime, cj, delta, ires, rpar, ipar)
        if (ires .lt. 0) then
            go to 370
        end if
        go to 300
    end if
    if (iernew .ne. 0) then
        go to 380
    end if
    return
370 iernls = 2
    if (ires .le. -2) then
        iernls = -1
    end if
    return
380 iernls = min(iernew, 2)
    return
end subroutine ddasid

subroutine dnsid(x, y, yprime, neq, icopt, id, res, wt, rpar, ipar, delta, r, yic, ypic, wm, iwm, cj, tscale, epcon, ratemx, maxit, stptol, icnflg, icnstr, iernew)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), wt(*), r(*)
    dimension :: id(*), delta(*), yic(*), ypic(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    dimension :: icnstr(*)
    external :: res
    parameter(lnni=19, llsoff=35)
    lsoff = iwm(llsoff)
    m = 0
    rate = 1.0d0
    rlx = 0.4d0
    call dslvd(neq, delta, wm, iwm)
    delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
    fnrm = delnrm
    if (tscale .gt. 0.0d0) then
        fnrm = fnrm*tscale*abs(cj)
    end if
    if (fnrm .le. epcon) then
        return
    end if
300 continue
    iwm(lnni) = iwm(lnni) + 1
    oldfnm = fnrm
call dlinsd(neq, y, x, yprime, cj, tscale, delta, delnrm, wt, lsoff, stptol, iret, res, ires, wm, iwm, fnrm, icopt, id, r, yic, ypic, icnflg, icnstr, rlx, rpar, ipar)
    rate = fnrm/oldfnm
    if (iret .ne. 0) then
        go to 390
    end if
    if (fnrm .le. epcon) then
        return
    end if
    m = m + 1
    if (m .ge. maxit) then
        go to 380
    end if
    call dcopy(neq, r, 1, delta, 1)
    delnrm = fnrm
    go to 300
380 if (rate .le. ratemx) then
        iernew = 1
    else
        iernew = 2
    end if
    return
390 if (ires .le. -2) then
        iernew = -1
    else
        iernew = 3
    end if
    return
end subroutine dnsid

subroutine dlinsd(neq, y, t, yprime, cj, tscale, p, pnrm, wt, lsoff, stptol, iret, res, ires, wm, iwm, fnrm, icopt, id, r, ynew, ypnew, icnflg, icnstr, rlx, rpar, ipar)
    implicit double precision(a - h, o - z)
    external :: res
    dimension :: y(*), yprime(*), wt(*), r(*), id(*)
    dimension :: wm(*), iwm(*)
    dimension :: ynew(*), ypnew(*), p(*), icnstr(*)
    dimension :: rpar(*), ipar(*)
    character :: msg*80
    parameter(lnre=12, lkprin=31)
    save :: alpha, one, two
    data alpha/1.0d-4/, one/1.0d0/, two/2.0d0/
    kprin = iwm(lkprin)
    f1nrm = fnrm*fnrm/two
    ratio = one
    if (kprin .ge. 2) then
        msg = "------ IN ROUTINE DLINSD-- PNRM = (R1)"
        call xerrwd(msg, 38, 901, 0, 0, 0, 0, 1, pnrm, 0.0d0)
    end if
    tau = pnrm
    rl = one
    if (icnflg .ne. 0) then
10      continue
        call dyypnw(neq, y, yprime, cj, rl, p, icopt, id, ynew, ypnew)
        call dcnstr(neq, y, ynew, icnstr, tau, rlx, iret, ivar)
        if (iret .eq. 1) then
            ratio1 = tau/pnrm
            ratio = ratio*ratio1
            do i = 1, neq
20              p(i) = p(i)*ratio1
            end do
            pnrm = tau
            if (kprin .ge. 2) then
                msg = "------ CONSTRAINT VIOL., PNRM = (R1), INDEX = (I1)"
                call xerrwd(msg, 50, 902, 0, 1, ivar, 0, 1, pnrm, 0.0d0)
            end if
            if (pnrm .le. stptol) then
                iret = 1
                return
            end if
            go to 10
        end if
    end if
    slpi = (-two)*f1nrm*ratio
    rlmin = stptol/pnrm
    if (lsoff .eq. 0 .and. kprin .ge. 2) then
        msg = "------ MIN. LAMBDA = (R1)"
        call xerrwd(msg, 25, 903, 0, 0, 0, 0, 1, rlmin, 0.0d0)
    end if
100 continue
    call dyypnw(neq, y, yprime, cj, rl, p, icopt, id, ynew, ypnew)
    call dfnrmd(neq, ynew, t, ypnew, r, cj, tscale, wt, res, ires, fnrmp, wm, iwm, rpar, ipar)
    iwm(lnre) = iwm(lnre) + 1
    if (ires .ne. 0) then
        iret = 2
        return
    end if
    if (lsoff .eq. 1) then
        go to 150
    end if
    f1nrmp = fnrmp*fnrmp/two
    if (kprin .ge. 2) then
        msg = "------ LAMBDA = (R1)"
        call xerrwd(msg, 20, 904, 0, 0, 0, 0, 1, rl, 0.0d0)
        msg = "------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)"
        call xerrwd(msg, 43, 905, 0, 0, 0, 0, 2, f1nrm, f1nrmp)
    end if
    if (f1nrmp .gt. f1nrm + alpha*slpi*rl) then
        go to 200
    end if
150 iret = 0
    call dcopy(neq, ynew, 1, y, 1)
    call dcopy(neq, ypnew, 1, yprime, 1)
    fnrm = fnrmp
    if (kprin .ge. 1) then
        msg = "------ LEAVING ROUTINE DLINSD, FNRM = (R1)"
        call xerrwd(msg, 42, 906, 0, 0, 0, 0, 1, fnrm, 0.0d0)
    end if
    return
200 continue
    if (rl .lt. rlmin) then
        iret = 1
        return
    end if
    rl = rl/two
    go to 100
end subroutine dlinsd

subroutine dfnrmd(neq, y, t, yprime, r, cj, tscale, wt, res, ires, fnorm, wm, iwm, rpar, ipar)
    implicit double precision(a - h, o - z)
    external :: res
    dimension :: y(*), yprime(*), wt(*), r(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    ires = 0
    call res(t, y, yprime, cj, r, ires, rpar, ipar)
    if (ires .lt. 0) then
        return
    end if
    call dslvd(neq, r, wm, iwm)
    fnorm = ddwnrm(neq, r, wt, rpar, ipar)
    if (tscale .gt. 0.0d0) then
        fnorm = fnorm*tscale*abs(cj)
    end if
    return
end subroutine dfnrmd

subroutine dnedd(x, y, yprime, neq, res, jacd, pdum, h, wt, jstart, idid, rpar, ipar, phi, gamma, dumsvr, delta, e, wm, iwm, cj, cjold, cjlast, s, uround, dume, dums, dumr, epcon, jcalc, jfdum, kp1, nonneg, ntype, iernls)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), wt(*)
    dimension :: delta(*), e(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    dimension :: phi(neq, *), gamma(*)
    external :: res, jacd
    parameter(lnre=12, lnje=13)
    save :: muldel, maxit, xrate
    data muldel/1/, maxit/4/, xrate/0.25d0/
    iertyp = 0
    if (ntype .ne. 0) then
        iertyp = 1
        go to 380
    end if
    if (jstart .eq. 0) then
        cjold = cj
        jcalc = -1
    end if
    iernls = 0
    temp1 = (1.0d0 - xrate)/(1.0d0 + xrate)
    temp2 = 1.0d0/temp1
    if (cj/cjold .lt. temp1 .or. cj/cjold .gt. temp2) then
        jcalc = -1
    end if
    if (cj .ne. cjlast) then
        s = 100.d0
    end if
300 continue
    ierj = 0
    ires = 0
    iernew = 0
    do i = 1, neq
        y(i) = phi(i, 1)
310     yprime(i) = 0.0d0
    end do
    do j = 2, kp1
        do i = 1, neq
            y(i) = y(i) + phi(i, j)
320         yprime(i) = yprime(i) + gamma(j)*phi(i, j)
        end do
330     continue
    end do
    pnorm = ddwnrm(neq, y, wt, rpar, ipar)
    tolnew = 100.d0*uround*pnorm
    iwm(lnre) = iwm(lnre) + 1
    call res(x, y, yprime, cj, delta, ires, rpar, ipar)
    if (ires .lt. 0) then
        go to 380
    end if
    if (jcalc .eq. -1) then
        iwm(lnje) = iwm(lnje) + 1
        jcalc = 0
        call dmatd(neq, x, y, yprime, delta, cj, h, ierj, wt, e, wm, iwm, res, ires, uround, jacd, rpar, ipar)
        cjold = cj
        s = 100.d0
        if (ires .lt. 0) then
            go to 380
        end if
        if (ierj .ne. 0) then
            go to 380
        end if
    end if
    temp1 = 2.0d0/(1.0d0 + cj/cjold)
call dnsd(x, y, yprime, neq, res, pdum, wt, rpar, ipar, dumsvr, delta, e, wm, iwm, cj, dums, dumr, dume, epcon, s, temp1, tolnew, muldel, maxit, ires, idum, iernew)
    if (iernew .gt. 0 .and. jcalc .ne. 0) then
        jcalc = -1
        go to 300
    end if
    if (iernew .ne. 0) then
        go to 380
    end if
375 if (nonneg .eq. 0) then
        go to 390
    end if
    do i = 1, neq
377     delta(i) = min(y(i), 0.0d0)
    end do
    delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
    if (delnrm .gt. epcon) then
        go to 380
    end if
    do i = 1, neq
378     e(i) = e(i) - delta(i)
    end do
    go to 390
380 continue
    if (ires .le. -2 .or. iertyp .ne. 0) then
        iernls = -1
        if (ires .le. -2) then
            idid = -11
        end if
        if (iertyp .ne. 0) then
            idid = -15
        end if
    else
        iernls = 1
        if (ires .lt. 0) then
            idid = -10
        end if
        if (ierj .ne. 0) then
            idid = -8
        end if
    end if
390 jcalc = 1
    return
end subroutine dnedd

subroutine dnsd(x, y, yprime, neq, res, pdum, wt, rpar, ipar, dumsvr, delta, e, wm, iwm, cj, dums, dumr, dume, epcon, s, confac, tolnew, muldel, maxit, ires, idum, iernew)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), wt(*), delta(*), e(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    external :: res
    parameter(lnre=12, lnni=19)
    m = 0
    do i = 1, neq
100     e(i) = 0.0d0
    end do
300 continue
    iwm(lnni) = iwm(lnni) + 1
    if (muldel .eq. 1) then
        do i = 1, neq
320         delta(i) = delta(i)*confac
        end do
    end if
    call dslvd(neq, delta, wm, iwm)
    do i = 1, neq
        y(i) = y(i) - delta(i)
        e(i) = e(i) - delta(i)
340     yprime(i) = yprime(i) - cj*delta(i)
    end do
    delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
    if (m .eq. 0) then
        oldnrm = delnrm
        if (delnrm .le. tolnew) then
            go to 370
        end if
    else
        rate = (delnrm/oldnrm)**(1.0d0/m)
        if (rate .gt. 0.9d0) then
            go to 380
        end if
        s = rate/(1.0d0 - rate)
    end if
    if (s*delnrm .le. epcon) then
        go to 370
    end if
    m = m + 1
    if (m .ge. maxit) then
        go to 380
    end if
    iwm(lnre) = iwm(lnre) + 1
    call res(x, y, yprime, cj, delta, ires, rpar, ipar)
    if (ires .lt. 0) then
        go to 380
    end if
    go to 300
370 return
380 continue
    if (ires .le. -2) then
        iernew = -1
    else
        iernew = 1
    end if
    return
end subroutine dnsd

subroutine dmatd(neq, x, y, yprime, delta, cj, h, ier, ewt, e, wm, iwm, res, ires, uround, jacd, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), delta(*), ewt(*), e(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    external :: res, jacd
    parameter(lml=1, lmu=2, lmtype=4, lnre=12, lnpd=22, llciwp=30)
    lipvt = iwm(llciwp)
    ier = 0
    mtype = iwm(lmtype)
    go to(100, 200, 300, 400, 500), mtype
100 lenpd = iwm(lnpd)
    do i = 1, lenpd
110     wm(i) = 0.0d0
    end do
    call jacd(x, y, yprime, wm, cj, rpar, ipar)
    go to 230
200 ires = 0
    nrow = 0
    squr = sqrt(uround)
    do i = 1, neq
        del = max(squr*max(abs(y(i)), abs(h*yprime(i))), 1.0d0/ewt(i))
        del = sign(del, h*yprime(i))
        del = y(i) + del - y(i)
        ysave = y(i)
        ypsave = yprime(i)
        y(i) = y(i) + del
        yprime(i) = yprime(i) + cj*del
        iwm(lnre) = iwm(lnre) + 1
        call res(x, y, yprime, cj, e, ires, rpar, ipar)
        if (ires .lt. 0) then
            return
        end if
        delinv = 1.0d0/del
        do l = 1, neq
220         wm(nrow + l) = (e(l) - delta(l))*delinv
        end do
        nrow = nrow + neq
        y(i) = ysave
        yprime(i) = ypsave
210     continue
    end do
230 call dgefa(wm, neq, neq, iwm(lipvt), ier)
    return
300 return
400 lenpd = iwm(lnpd)
    do i = 1, lenpd
410     wm(i) = 0.0d0
    end do
    call jacd(x, y, yprime, wm, cj, rpar, ipar)
    meband = 2*iwm(lml) + iwm(lmu) + 1
    go to 550
500 mband = iwm(lml) + iwm(lmu) + 1
    mba = min0(mband, neq)
    meband = mband + iwm(lml)
    meb1 = meband - 1
    msave = neq/mband + 1
    isave = iwm(lnpd)
    ipsave = isave + msave
    ires = 0
    squr = sqrt(uround)
    do j = 1, mba
        do n = j, neq, mband
            k = (n - j)/mband + 1
            wm(isave + k) = y(n)
            wm(ipsave + k) = yprime(n)
            del = max(squr*max(abs(y(n)), abs(h*yprime(n))), 1.0d0/ewt(n))
            del = sign(del, h*yprime(n))
            del = y(n) + del - y(n)
            y(n) = y(n) + del
510         yprime(n) = yprime(n) + cj*del
        end do
        iwm(lnre) = iwm(lnre) + 1
        call res(x, y, yprime, cj, e, ires, rpar, ipar)
        if (ires .lt. 0) then
            return
        end if
        do n = j, neq, mband
            k = (n - j)/mband + 1
            y(n) = wm(isave + k)
            yprime(n) = wm(ipsave + k)
            del = max(squr*max(abs(y(n)), abs(h*yprime(n))), 1.0d0/ewt(n))
            del = sign(del, h*yprime(n))
            del = y(n) + del - y(n)
            delinv = 1.0d0/del
            i1 = max0(1, n - iwm(lmu))
            i2 = min0(neq, n + iwm(lml))
            ii = n*meb1 - iwm(lml)
            do i = i1, i2
520             wm(ii + i) = (e(i) - delta(i))*delinv
            end do
530         continue
        end do
540     continue
    end do
550 call dgbfa(wm, meband, neq, iwm(lml), iwm(lmu), iwm(lipvt), ier)
    return
end subroutine dmatd

subroutine dslvd(neq, delta, wm, iwm)
    implicit double precision(a - h, o - z)
    dimension :: delta(*), wm(*), iwm(*)
    parameter(lml=1, lmu=2, lmtype=4, llciwp=30)
    lipvt = iwm(llciwp)
    mtype = iwm(lmtype)
    go to(100, 100, 300, 400, 400), mtype
100 call dgesl(wm, neq, neq, iwm(lipvt), delta, 0)
    return
300 continue
    return
400 meband = 2*iwm(lml) + iwm(lmu) + 1
    call dgbsl(wm, meband, neq, iwm(lml), iwm(lmu), iwm(lipvt), delta, 0)
    return
end subroutine dslvd

subroutine ddasik(x, y, yprime, neq, icopt, id, res, jack, psol, h, tscale, wt, jskip, rpar, ipar, savr, delta, r, yic, ypic, pwk, wm, iwm, cj, uround, epli, sqrtn, rsqrtn, epcon, ratemx, stptol, jflg, icnflg, icnstr, iernls)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), id(*), wt(*), icnstr(*)
    dimension :: savr(*), delta(*), r(*), yic(*), ypic(*), pwk(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    external :: res, jack, psol
    parameter(lnre=12, lnje=13, llocwp=29, llciwp=30)
    parameter(lmxnit=32, lmxnj=33)
    lwp = iwm(llocwp)
    liwp = iwm(llciwp)
    mxnit = iwm(lmxnit)
    mxnj = iwm(lmxnj)
    iernls = 0
    nj = 0
    eplin = epli*epcon
    ires = 0
    iwm(lnre) = iwm(lnre) + 1
    call res(x, y, yprime, cj, delta, ires, rpar, ipar)
    if (ires .lt. 0) then
        go to 370
    end if
300 continue
    ierpj = 0
    ires = 0
    iernew = 0
    if (jflg .eq. 1 .and. jskip .eq. 0) then
        nj = nj + 1
        iwm(lnje) = iwm(lnje) + 1
        call jack(res, ires, neq, x, y, yprime, wt, delta, r, h, cj, wm(lwp), iwm(liwp), ierpj, rpar, ipar)
        if (ires .lt. 0 .or. ierpj .ne. 0) then
            go to 370
        end if
    end if
    jskip = 0
call dnsik(x, y, yprime, neq, icopt, id, res, psol, wt, rpar, ipar, savr, delta, r, yic, ypic, pwk, wm, iwm, cj, tscale, sqrtn, rsqrtn, eplin, epcon, ratemx, mxnit, stptol, icnflg, icnstr, iernew)
    if (iernew .eq. 1 .and. nj .lt. mxnj .and. jflg .eq. 1) then
        call dcopy(neq, savr, 1, delta, 1)
        go to 300
    end if
    if (iernew .ne. 0) then
        go to 380
    end if
    return
370 iernls = 2
    if (ires .le. -2) then
        iernls = -1
    end if
    return
380 iernls = min(iernew, 2)
    return
end subroutine ddasik

subroutine dnsik(x, y, yprime, neq, icopt, id, res, psol, wt, rpar, ipar, savr, delta, r, yic, ypic, pwk, wm, iwm, cj, tscale, sqrtn, rsqrtn, eplin, epcon, ratemx, maxit, stptol, icnflg, icnstr, iernew)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), wt(*), id(*), delta(*), r(*), savr(*)
    dimension :: yic(*), ypic(*), pwk(*), wm(*), iwm(*), rpar(*), ipar(*)
    dimension :: icnstr(*)
    external :: res, psol
    parameter(lnni=19, lnps=21, llocwp=29, llciwp=30)
    parameter(llsoff=35, lstol=14)
    lsoff = iwm(llsoff)
    m = 0
    rate = 1.0d0
    lwp = iwm(llocwp)
    liwp = iwm(llciwp)
    rlx = 0.4d0
    call dcopy(neq, delta, 1, savr, 1)
call dfnrmk(neq, y, x, yprime, savr, r, cj, tscale, wt, sqrtn, rsqrtn, res, ires, psol, 1, ier, fnrm, eplin, wm(lwp), iwm(liwp), pwk, rpar, ipar)
    iwm(lnps) = iwm(lnps) + 1
    if (ier .ne. 0) then
        iernew = 3
        return
    end if
    if (fnrm .le. epcon) then
        return
    end if
    fnrm0 = fnrm
300 continue
    iwm(lnni) = iwm(lnni) + 1
    call dslvk(neq, y, x, yprime, savr, delta, wt, wm, iwm, res, ires, psol, iersl, cj, eplin, sqrtn, rsqrtn, rhok, rpar, ipar)
    if (ires .ne. 0 .or. iersl .ne. 0) then
        go to 390
    end if
    delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
    if (delnrm .eq. 0.0d0) then
        return
    end if
    oldfnm = fnrm
call dlinsk(neq, y, x, yprime, savr, cj, tscale, delta, delnrm, wt, sqrtn, rsqrtn, lsoff, stptol, iret, res, ires, psol, wm, iwm, rhok, fnrm, icopt, id, wm(lwp), iwm(liwp), r, eplin, yic, ypic, pwk, icnflg, icnstr, rlx, rpar, ipar)
    rate = fnrm/oldfnm
    if (iret .ne. 0) then
        go to 390
    end if
    if (fnrm .le. epcon) then
        return
    end if
    m = m + 1
    if (m .ge. maxit) then
        go to 380
    end if
    call dcopy(neq, savr, 1, delta, 1)
    go to 300
380 if (rate .le. ratemx .or. fnrm .le. 0.1d0*fnrm0) then
        iernew = 1
    else
        iernew = 2
    end if
    return
390 if (ires .le. -2 .or. iersl .lt. 0) then
        iernew = -1
    else
        iernew = 3
        if (ires .eq. 0 .and. iersl .eq. 1 .and. m .ge. 2 .and. rate .lt. 1.0d0) then
            iernew = 1
        end if
    end if
    return
end subroutine dnsik

subroutine dlinsk(neq, y, t, yprime, savr, cj, tscale, p, pnrm, wt, sqrtn, rsqrtn, lsoff, stptol, iret, res, ires, psol, wm, iwm, rhok, fnrm, icopt, id, wp, iwp, r, eplin, ynew, ypnew, pwk, icnflg, icnstr, rlx, rpar, ipar)
    implicit double precision(a - h, o - z)
    external :: res, psol
    dimension :: y(*), yprime(*), p(*), wt(*), savr(*), r(*), id(*)
    dimension :: wm(*), iwm(*), ynew(*), ypnew(*), pwk(*), icnstr(*)
    dimension :: wp(*), iwp(*), rpar(*), ipar(*)
    character :: msg*80
    parameter(lnre=12, lnps=21, lkprin=31)
    save :: alpha, one, two
    data alpha/1.0d-4/, one/1.0d0/, two/2.0d0/
    kprin = iwm(lkprin)
    f1nrm = fnrm*fnrm/two
    ratio = one
    if (kprin .ge. 2) then
        msg = "------ IN ROUTINE DLINSK-- PNRM = (R1)"
        call xerrwd(msg, 38, 921, 0, 0, 0, 0, 1, pnrm, 0.0d0)
    end if
    tau = pnrm
    rl = one
    if (icnflg .ne. 0) then
10      continue
        call dyypnw(neq, y, yprime, cj, rl, p, icopt, id, ynew, ypnew)
        call dcnstr(neq, y, ynew, icnstr, tau, rlx, iret, ivar)
        if (iret .eq. 1) then
            ratio1 = tau/pnrm
            ratio = ratio*ratio1
            do i = 1, neq
20              p(i) = p(i)*ratio1
            end do
            pnrm = tau
            if (kprin .ge. 2) then
                msg = "------ CONSTRAINT VIOL., PNRM = (R1), INDEX = (I1)"
                call xerrwd(msg, 50, 922, 0, 1, ivar, 0, 1, pnrm, 0.0d0)
            end if
            if (pnrm .le. stptol) then
                iret = 1
                return
            end if
            go to 10
        end if
    end if
    slpi = (-two)*f1nrm*ratio
    rlmin = stptol/pnrm
    if (lsoff .eq. 0 .and. kprin .ge. 2) then
        msg = "------ MIN. LAMBDA = (R1)"
        call xerrwd(msg, 25, 923, 0, 0, 0, 0, 1, rlmin, 0.0d0)
    end if
100 continue
    call dyypnw(neq, y, yprime, cj, rl, p, icopt, id, ynew, ypnew)
call dfnrmk(neq, ynew, t, ypnew, savr, r, cj, tscale, wt, sqrtn, rsqrtn, res, ires, psol, 0, ier, fnrmp, eplin, wp, iwp, pwk, rpar, ipar)
    iwm(lnre) = iwm(lnre) + 1
    if (ires .ge. 0) then
        iwm(lnps) = iwm(lnps) + 1
    end if
    if (ires .ne. 0 .or. ier .ne. 0) then
        iret = 2
        return
    end if
    if (lsoff .eq. 1) then
        go to 150
    end if
    f1nrmp = fnrmp*fnrmp/two
    if (kprin .ge. 2) then
        msg = "------ LAMBDA = (R1)"
        call xerrwd(msg, 20, 924, 0, 0, 0, 0, 1, rl, 0.0d0)
        msg = "------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)"
        call xerrwd(msg, 43, 925, 0, 0, 0, 0, 2, f1nrm, f1nrmp)
    end if
    if (f1nrmp .gt. f1nrm + alpha*slpi*rl) then
        go to 200
    end if
150 iret = 0
    call dcopy(neq, ynew, 1, y, 1)
    call dcopy(neq, ypnew, 1, yprime, 1)
    fnrm = fnrmp
    if (kprin .ge. 1) then
        msg = "------ LEAVING ROUTINE DLINSK, FNRM = (R1)"
        call xerrwd(msg, 42, 926, 0, 0, 0, 0, 1, fnrm, 0.0d0)
    end if
    return
200 continue
    if (rl .lt. rlmin) then
        iret = 1
        return
    end if
    rl = rl/two
    go to 100
end subroutine dlinsk

subroutine dfnrmk(neq, y, t, yprime, savr, r, cj, tscale, wt, sqrtn, rsqrtn, res, ires, psol, irin, ier, fnorm, eplin, wp, iwp, pwk, rpar, ipar)
    implicit double precision(a - h, o - z)
    external :: res, psol
    dimension :: y(*), yprime(*), wt(*), savr(*), r(*), pwk(*)
    dimension :: wp(*), iwp(*), rpar(*), ipar(*)
    if (irin .eq. 0) then
        ires = 0
        call res(t, y, yprime, cj, savr, ires, rpar, ipar)
        if (ires .lt. 0) then
            return
        end if
    end if
    call dcopy(neq, savr, 1, r, 1)
    call dscal(neq, rsqrtn, wt, 1)
    ier = 0
    call psol(neq, t, y, yprime, savr, pwk, cj, wt, wp, iwp, r, eplin, ier, rpar, ipar)
    call dscal(neq, sqrtn, wt, 1)
    if (ier .ne. 0) then
        return
    end if
    fnorm = ddwnrm(neq, r, wt, rpar, ipar)
    if (tscale .gt. 0.0d0) then
        fnorm = fnorm*tscale*abs(cj)
    end if
    return
end subroutine dfnrmk

subroutine dnedk(x, y, yprime, neq, res, jack, psol, h, wt, jstart, idid, rpar, ipar, phi, gamma, savr, delta, e, wm, iwm, cj, cjold, cjlast, s, uround, epli, sqrtn, rsqrtn, epcon, jcalc, jflg, kp1, nonneg, ntype, iernls)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), wt(*)
    dimension :: phi(neq, *), savr(*), delta(*), e(*)
    dimension :: wm(*), iwm(*)
    dimension :: gamma(*), rpar(*), ipar(*)
    external :: res, jack, psol
    parameter(lnre=12, lnje=13, llocwp=29, llciwp=30)
    save :: muldel, maxit, xrate
    data muldel/0/, maxit/4/, xrate/0.25d0/
    iertyp = 0
    if (ntype .ne. 1) then
        iertyp = 1
        go to 380
    end if
    if (jstart .eq. 0) then
        cjold = cj
        jcalc = -1
        s = 100.d0
    end if
    iernls = 0
    lwp = iwm(llocwp)
    liwp = iwm(llciwp)
    if (jflg .ne. 0) then
        temp1 = (1.0d0 - xrate)/(1.0d0 + xrate)
        temp2 = 1.0d0/temp1
        if (cj/cjold .lt. temp1 .or. cj/cjold .gt. temp2) then
            jcalc = -1
        end if
        if (cj .ne. cjlast) then
            s = 100.d0
        end if
    else
        jcalc = 0
    end if
300 continue
    ierpj = 0
    ires = 0
    iersl = 0
    iernew = 0
    do i = 1, neq
        y(i) = phi(i, 1)
310     yprime(i) = 0.0d0
    end do
    do j = 2, kp1
        do i = 1, neq
            y(i) = y(i) + phi(i, j)
320         yprime(i) = yprime(i) + gamma(j)*phi(i, j)
        end do
330     continue
    end do
    eplin = epli*epcon
    tolnew = eplin
    iwm(lnre) = iwm(lnre) + 1
    call res(x, y, yprime, cj, delta, ires, rpar, ipar)
    if (ires .lt. 0) then
        go to 380
    end if
    if (jcalc .eq. -1) then
        iwm(lnje) = iwm(lnje) + 1
        jcalc = 0
        call jack(res, ires, neq, x, y, yprime, wt, delta, e, h, cj, wm(lwp), iwm(liwp), ierpj, rpar, ipar)
        cjold = cj
        s = 100.d0
        if (ires .lt. 0) then
            go to 380
        end if
        if (ierpj .ne. 0) then
            go to 380
        end if
    end if
call dnsk(x, y, yprime, neq, res, psol, wt, rpar, ipar, savr, delta, e, wm, iwm, cj, sqrtn, rsqrtn, eplin, epcon, s, temp1, tolnew, muldel, maxit, ires, iersl, iernew)
    if (iernew .gt. 0 .and. jcalc .ne. 0) then
        jcalc = -1
        go to 300
    end if
    if (iernew .ne. 0) then
        go to 380
    end if
    if (nonneg .eq. 0) then
        go to 390
    end if
    do i = 1, neq
360     delta(i) = min(y(i), 0.0d0)
    end do
    delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
    if (delnrm .gt. epcon) then
        go to 380
    end if
    do i = 1, neq
370     e(i) = e(i) - delta(i)
    end do
    go to 390
380 continue
    if (ires .le. -2 .or. iersl .lt. 0 .or. iertyp .ne. 0) then
        iernls = -1
        if (ires .le. -2) then
            idid = -11
        end if
        if (iersl .lt. 0) then
            idid = -13
        end if
        if (iertyp .ne. 0) then
            idid = -15
        end if
    else
        iernls = 1
        if (ires .eq. -1) then
            idid = -10
        end if
        if (ierpj .ne. 0) then
            idid = -5
        end if
        if (iersl .gt. 0) then
            idid = -14
        end if
    end if
390 jcalc = 1
    return
end subroutine dnedk

subroutine dnsk(x, y, yprime, neq, res, psol, wt, rpar, ipar, savr, delta, e, wm, iwm, cj, sqrtn, rsqrtn, eplin, epcon, s, confac, tolnew, muldel, maxit, ires, iersl, iernew)
    implicit double precision(a - h, o - z)
    dimension :: y(*), yprime(*), wt(*), delta(*), e(*), savr(*)
    dimension :: wm(*), iwm(*), rpar(*), ipar(*)
    external :: res, psol
    parameter(lnni=19, lnre=12)
    m = 0
    do i = 1, neq
100     e(i) = 0.0d0
    end do
300 continue
    iwm(lnni) = iwm(lnni) + 1
    if (muldel .eq. 1) then
        do i = 1, neq
320         delta(i) = delta(i)*confac
        end do
    end if
    do i = 1, neq
340     savr(i) = delta(i)
    end do
    call dslvk(neq, y, x, yprime, savr, delta, wt, wm, iwm, res, ires, psol, iersl, cj, eplin, sqrtn, rsqrtn, rhok, rpar, ipar)
    if (ires .ne. 0 .or. iersl .ne. 0) then
        go to 380
    end if
    do i = 1, neq
        y(i) = y(i) - delta(i)
        e(i) = e(i) - delta(i)
360     yprime(i) = yprime(i) - cj*delta(i)
    end do
    delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
    if (m .eq. 0) then
        oldnrm = delnrm
        if (delnrm .le. tolnew) then
            go to 370
        end if
    else
        rate = (delnrm/oldnrm)**(1.0d0/m)
        if (rate .gt. 0.9d0) then
            go to 380
        end if
        s = rate/(1.0d0 - rate)
    end if
    if (s*delnrm .le. epcon) then
        go to 370
    end if
    m = m + 1
    if (m .ge. maxit) then
        go to 380
    end if
    iwm(lnre) = iwm(lnre) + 1
    call res(x, y, yprime, cj, delta, ires, rpar, ipar)
    if (ires .lt. 0) then
        go to 380
    end if
    go to 300
370 return
380 continue
    if (ires .le. -2 .or. iersl .lt. 0) then
        iernew = -1
    else
        iernew = 1
    end if
    return
end subroutine dnsk

subroutine dslvk(neq, y, tn, yprime, savr, x, ewt, wm, iwm, res, ires, psol, iersl, cj, eplin, sqrtn, rsqrtn, rhok, rpar, ipar)
    integer :: neq, iwm, ires, iersl, ipar
    double precision :: y, tn, yprime, savr, x, ewt, wm, cj, eplin, sqrtn, rsqrtn, rhok, rpar
    dimension :: y(*), yprime(*), savr(*), x(*), ewt(*), wm(*), iwm(*), rpar(*), ipar(*)
    integer :: iflag, irst, nrsts, nrmax, lr, ldl, lhes, lgmr, lq, lv, lwk, lz, maxlp1, npsl
    integer :: nli, nps, ncfl, nre, maxl, kmp, miter
    external :: res, psol
    parameter(lnre=12, lncfl=16, lnli=20, lnps=21)
    parameter(llocwp=29, llciwp=30)
    parameter(lmiter=23, lmaxl=24, lkmp=25, lnrmax=26)
    data irst/1/
    liwp = iwm(llciwp)
    nli = iwm(lnli)
    nps = iwm(lnps)
    ncfl = iwm(lncfl)
    nre = iwm(lnre)
    lwp = iwm(llocwp)
    maxl = iwm(lmaxl)
    kmp = iwm(lkmp)
    nrmax = iwm(lnrmax)
    miter = iwm(lmiter)
    iersl = 0
    ires = 0
    maxlp1 = maxl + 1
    lv = 1
    lr = lv + neq*maxl
    lhes = lr + neq + 1
    lq = lhes + maxl*maxlp1
    lwk = lq + 2*maxl
    ldl = lwk + min0(1, maxl - kmp)*neq
    lz = ldl + neq
    call dscal(neq, rsqrtn, ewt, 1)
    call dcopy(neq, x, 1, wm(lr), 1)
    do i = 1, neq
110     x(i) = 0.d0
    end do
    nrsts = -1
115 continue
    nrsts = nrsts + 1
    if (nrsts .gt. 0) then
        call dcopy(neq, wm(ldl), 1, wm(lr), 1)
    end if
call dspigm(neq, tn, y, yprime, savr, wm(lr), ewt, maxl, maxlp1, kmp, eplin, cj, res, ires, nres, psol, npsl, wm(lz), wm(lv), wm(lhes), wm(lq), lgmr, wm(lwp), iwm(liwp), wm(lwk), wm(ldl), rhok, iflag, irst, nrsts, rpar, ipar)
    nli = nli + lgmr
    nps = nps + npsl
    nre = nre + nres
    do i = 1, neq
120     x(i) = x(i) + wm(lz + i - 1)
    end do
    if (iflag .eq. 1 .and. nrsts .lt. nrmax .and. ires .eq. 0) then
        go to 115
    end if
    if (ires .lt. 0) then
        ncfl = ncfl + 1
    else
        if (iflag .ne. 0) then
            ncfl = ncfl + 1
            if (iflag .gt. 0) then
                iersl = 1
            end if
            if (iflag .lt. 0) then
                iersl = -1
            end if
        end if
    end if
    iwm(lnli) = nli
    iwm(lnps) = nps
    iwm(lncfl) = ncfl
    iwm(lnre) = nre
    call dscal(neq, sqrtn, ewt, 1)
    return
end subroutine dslvk

subroutine dspigm(neq, tn, y, yprime, savr, r, wght, maxl, maxlp1, kmp, eplin, cj, res, ires, nre, psol, npsl, z, v, hes, q, lgmr, wp, iwp, wk, dl, rhok, iflag, irst, nrsts, rpar, ipar)
    integer :: neq, maxl, maxlp1, kmp, ires, nre, npsl, lgmr, iwp, iflag, irst, nrsts, ipar
    double precision :: tn, y, yprime, savr, r, wght, eplin, cj, z, v, hes, q, wp, wk, dl, rhok, rpar
dimension :: y(*), yprime(*), savr(*), r(*), wght(*), z(*), v(neq,*), hes(maxlp1,*), q(*), wp(*), iwp(*), wk(*), dl(*), rpar(*), ipar(*)
    integer :: i, ier, info, ip1, i2, j, k, ll, llp1
    double precision :: rnrm, c, dlnrm, prod, rho, s, snormw, dnrm2, tem
    external :: res, psol
    ier = 0
    iflag = 0
    lgmr = 0
    npsl = 0
    nre = 0
    do i = 1, neq
10      z(i) = 0.0d0
    end do
    if (nrsts .eq. 0) then
        call psol(neq, tn, y, yprime, savr, wk, cj, wght, wp, iwp, r, eplin, ier, rpar, ipar)
        npsl = 1
        if (ier .ne. 0) then
            go to 300
        end if
        do i = 1, neq
30          v(i, 1) = r(i)*wght(i)
        end do
    else
        do i = 1, neq
35          v(i, 1) = r(i)
        end do
    end if
    rnrm = dnrm2(neq, v, 1)
    if (rnrm .le. eplin) then
        rhok = rnrm
        return
    end if
    tem = 1.0d0/rnrm
    call dscal(neq, tem, v(1, 1), 1)
    do j = 1, maxl
        do i = 1, maxlp1
60          hes(i, j) = 0.0d0
        end do
65      continue
    end do
    prod = 1.0d0
    do ll = 1, maxl
        lgmr = ll
    call datv(neq, y, tn, yprime, savr, v(1, ll), wght, z, res, ires, psol, v(1, ll + 1), wk, wp, iwp, cj, eplin, ier, nre, npsl, rpar, ipar)
        if (ires .lt. 0) then
            return
        end if
        if (ier .ne. 0) then
            go to 300
        end if
        call dorth(v(1, ll + 1), v, hes, neq, ll, maxlp1, kmp, snormw)
        hes(ll + 1, ll) = snormw
        call dheqr(hes, maxlp1, ll, q, info, ll)
        if (info .eq. ll) then
            go to 120
        end if
        prod = prod*q(2*ll)
        rho = abs(prod*rnrm)
        if (ll .gt. kmp .and. kmp .lt. maxl) then
            if (ll .eq. kmp + 1) then
                call dcopy(neq, v(1, 1), 1, dl, 1)
                do i = 1, kmp
                    ip1 = i + 1
                    i2 = i*2
                    s = q(i2)
                    c = q(i2 - 1)
                    do k = 1, neq
70                      dl(k) = s*dl(k) + c*v(k, ip1)
                    end do
75                  continue
                end do
            end if
            s = q(2*ll)
            c = q(2*ll - 1)/snormw
            llp1 = ll + 1
            do k = 1, neq
80              dl(k) = s*dl(k) + c*v(k, llp1)
            end do
            dlnrm = dnrm2(neq, dl, 1)
            rho = rho*dlnrm
        end if
        if (rho .le. eplin) then
            go to 200
        end if
        if (ll .eq. maxl) then
            go to 100
        end if
        tem = 1.0d0/snormw
        call dscal(neq, tem, v(1, ll + 1), 1)
90      continue
    end do
100 continue
    if (rho .lt. rnrm) then
        go to 150
    end if
120 continue
    iflag = 2
    do i = 1, neq
130     z(i) = 0.d0
    end do
    return
150 iflag = 1
    if (irst .gt. 0) then
        if (kmp .eq. maxl) then
            call dcopy(neq, v(1, 1), 1, dl, 1)
            maxlm1 = maxl - 1
            do i = 1, maxlm1
                ip1 = i + 1
                i2 = i*2
                s = q(i2)
                c = q(i2 - 1)
                do k = 1, neq
170                 dl(k) = s*dl(k) + c*v(k, ip1)
                end do
175             continue
            end do
            s = q(2*maxl)
            c = q(2*maxl - 1)/snormw
            do k = 1, neq
180             dl(k) = s*dl(k) + c*v(k, maxlp1)
            end do
        end if
        tem = rnrm*prod
        call dscal(neq, tem, dl, 1)
    end if
200 continue
    ll = lgmr
    llp1 = ll + 1
    do k = 1, llp1
210     r(k) = 0.0d0
    end do
    r(1) = rnrm
    call dhels(hes, maxlp1, ll, q, r)
    do k = 1, neq
220     z(k) = 0.0d0
    end do
    do i = 1, ll
        call daxpy(neq, r(i), v(1, i), 1, z, 1)
230     continue
    end do
    do i = 1, neq
240     z(i) = z(i)/wght(i)
    end do
    rhok = rho
    return
300 continue
    if (ier .lt. 0) then
        iflag = -1
    end if
    if (ier .gt. 0) then
        iflag = 3
    end if
    return
end subroutine dspigm

subroutine datv(neq, y, tn, yprime, savr, v, wght, yptem, res, ires, psol, z, vtem, wp, iwp, cj, eplin, ier, nre, npsl, rpar, ipar)
    integer :: neq, ires, iwp, ier, nre, npsl, ipar
    double precision :: y, tn, yprime, savr, v, wght, yptem, z, vtem, wp, cj, rpar
    dimension :: y(*), yprime(*), savr(*), v(*), wght(*), yptem(*), z(*), vtem(*), wp(*), iwp(*), rpar(*), ipar(*)
    integer :: i
    double precision :: eplin
    external :: res, psol
    ires = 0
    do i = 1, neq
10      vtem(i) = v(i)/wght(i)
    end do
    ier = 0
    do i = 1, neq
        yptem(i) = yprime(i) + vtem(i)*cj
20      z(i) = y(i) + vtem(i)
    end do
    continue
    call res(tn, z, yptem, cj, vtem, ires, rpar, ipar)
    nre = nre + 1
    if (ires .lt. 0) then
        return
    end if
    do i = 1, neq
70      z(i) = vtem(i) - savr(i)
    end do
    call psol(neq, tn, y, yprime, savr, yptem, cj, wght, wp, iwp, z, eplin, ier, rpar, ipar)
    npsl = npsl + 1
    if (ier .ne. 0) then
        return
    end if
    do i = 1, neq
90      z(i) = z(i)*wght(i)
    end do
    return
end subroutine datv

subroutine dorth(vnew, v, hes, n, ll, ldhes, kmp, snormw)
    integer :: n, ll, ldhes, kmp
    double precision :: vnew, v, hes, snormw
    dimension :: vnew(*), v(n, *), hes(ldhes, *)
    integer :: i, i0
    double precision :: arg, ddot, dnrm2, sumdsq, tem, vnrm
    vnrm = dnrm2(n, vnew, 1)
    i0 = max0(1, ll - kmp + 1)
    do i = i0, ll
        hes(i, ll) = ddot(n, v(1, i), 1, vnew, 1)
        tem = -hes(i, ll)
        call daxpy(n, tem, v(1, i), 1, vnew, 1)
10      continue
    end do
    snormw = dnrm2(n, vnew, 1)
    if (vnrm + 0.001d0*snormw .ne. vnrm) then
        return
    end if
    sumdsq = 0.0d0
    do i = i0, ll
        tem = -ddot(n, v(1, i), 1, vnew, 1)
        if (hes(i, ll) + 0.001d0*tem .eq. hes(i, ll)) then
            go to 30
        end if
        hes(i, ll) = hes(i, ll) - tem
        call daxpy(n, tem, v(1, i), 1, vnew, 1)
        sumdsq = sumdsq + tem**2
30      continue
    end do
    if (sumdsq .eq. 0.0d0) then
        return
    end if
    arg = max(0.0d0, snormw**2 - sumdsq)
    snormw = sqrt(arg)
    return
end subroutine dorth

subroutine dheqr(a, lda, n, q, info, ijob)
    integer :: lda, n, info, ijob
    double precision :: a(lda, *), q(*)
    integer :: i, iq, j, k, km1, kp1, nm1
    double precision :: c, s, t, t1, t2
    if (ijob .gt. 1) then
        go to 70
    end if
    info = 0
    do k = 1, n
        km1 = k - 1
        kp1 = k + 1
        if (km1 .lt. 1) then
            go to 20
        end if
        do j = 1, km1
            i = 2*(j - 1) + 1
            t1 = a(j, k)
            t2 = a(j + 1, k)
            c = q(i)
            s = q(i + 1)
            a(j, k) = c*t1 - s*t2
            a(j + 1, k) = s*t1 + c*t2
10          continue
        end do
20      continue
        iq = 2*km1 + 1
        t1 = a(k, k)
        t2 = a(kp1, k)
        if (t2 .ne. 0.0d0) then
            go to 30
        end if
        c = 1.0d0
        s = 0.0d0
        go to 50
30      continue
        if (abs(t2) .lt. abs(t1)) then
            go to 40
        end if
        t = t1/t2
        s = (-1.0d0)/sqrt(1.0d0 + t*t)
        c = (-s)*t
        go to 50
40      continue
        t = t2/t1
        c = 1.0d0/sqrt(1.0d0 + t*t)
        s = (-c)*t
50      continue
        q(iq) = c
        q(iq + 1) = s
        a(k, k) = c*t1 - s*t2
        if (a(k, k) .eq. 0.0d0) then
            info = k
        end if
60      continue
    end do
    return
70  continue
    nm1 = n - 1
    do k = 1, nm1
        i = 2*(k - 1) + 1
        t1 = a(k, n)
        t2 = a(k + 1, n)
        c = q(i)
        s = q(i + 1)
        a(k, n) = c*t1 - s*t2
        a(k + 1, n) = s*t1 + c*t2
100     continue
    end do
    info = 0
    t1 = a(n, n)
    t2 = a(n + 1, n)
    if (t2 .ne. 0.0d0) then
        go to 110
    end if
    c = 1.0d0
    s = 0.0d0
    go to 130
110 continue
    if (abs(t2) .lt. abs(t1)) then
        go to 120
    end if
    t = t1/t2
    s = (-1.0d0)/sqrt(1.0d0 + t*t)
    c = (-s)*t
    go to 130
120 continue
    t = t2/t1
    c = 1.0d0/sqrt(1.0d0 + t*t)
    s = (-c)*t
130 continue
    iq = 2*n - 1
    q(iq) = c
    q(iq + 1) = s
    a(n, n) = c*t1 - s*t2
    if (a(n, n) .eq. 0.0d0) then
        info = n
    end if
    return
end subroutine dheqr

subroutine dhels(a, lda, n, q, b)
    integer :: lda, n
    double precision :: a(lda, *), b(*), q(*)
    integer :: iq, k, kb, kp1
    double precision :: c, s, t, t1, t2
    do k = 1, n
        kp1 = k + 1
        iq = 2*(k - 1) + 1
        c = q(iq)
        s = q(iq + 1)
        t1 = b(k)
        t2 = b(kp1)
        b(k) = c*t1 - s*t2
        b(kp1) = s*t1 + c*t2
20      continue
    end do
    do kb = 1, n
        k = n + 1 - kb
        b(k) = b(k)/a(k, k)
        t = -b(k)
        call daxpy(k - 1, t, a(1, k), 1, b(1), 1)
40      continue
    end do
    return
end subroutine dhels

subroutine dspsetup(neq, lwp, liwp, rpar, ipar, ierr, lwp_min, liwp_min)
    implicit none
    integer :: neq
    real(8) :: rpar(*)
    integer :: ipar(*)
    integer :: lwp
    integer :: liwp
    integer :: ierr
    integer :: lwp_min
    integer :: liwp_min
    integer :: lbw, ubw, lenplumx, ljac, ljaci, ljacj, lrownrms, lrwk1, liwk1
    integer :: lenpfac, lenplufac, lfililut, ipremeth, neqp1, nnzmx
    integer :: lplu, lju, ljlu, lperm, lqperm, llevels, lmask
    integer :: isrnorm
    integer :: normtype
    integer :: ireorder
    integer :: jacout
    integer :: jscalcol
    real(8) :: tolilut, permtol
    lbw = ipar(1)
    if (lbw .le. 0) then
        ierr = 1
        return
    end if
    ubw = ipar(2)
    if (ubw .le. 0) then
        ierr = 2
        return
    end if
    lenpfac = ipar(3)
    if (lenpfac .le. 1) then
        ierr = 3
        return
    end if
    lenplufac = ipar(4)
    if (lenplufac .le. 1) then
        ierr = 4
        return
    end if
    ipremeth = ipar(5)
    if (ipremeth .ne. 1 .and. ipremeth .ne. 2) then
        ierr = 5
        return
    end if
    lfililut = ipar(6)
    if (lfililut .lt. 0) then
        ierr = 6
        return
    end if
    ireorder = ipar(7)
    if (ireorder .lt. 0 .or. ireorder .gt. 1) then
        ierr = 7
        return
    end if
    isrnorm = ipar(8)
    if (isrnorm .lt. 0 .or. isrnorm .gt. 1) then
        ierr = 8
        return
    end if
    normtype = ipar(9)
    if (normtype .lt. 0 .or. normtype .gt. 2) then
        ierr = 9
        return
    end if
    jacout = ipar(10)
    if (jacout .lt. 0 .or. jacout .gt. 1) then
        ierr = 10
        return
    end if
    jscalcol = ipar(11)
    if (jscalcol .lt. 0 .or. jscalcol .gt. 1) then
        ierr = 11
        return
    end if
    if (jacout .eq. 1) then
        if (ipar(29) .le. 0) then
            ierr = 12
            return
        end if
    end if
    tolilut = rpar(1)
    if (tolilut .lt. 0.) then
        ierr = 21
        return
    end if
    if (ipremeth .eq. 2) then
        permtol = rpar(2)
        if (permtol .lt. 0.) then
            ierr = 22
            return
        end if
    end if
    neqp1 = neq + 1
    nnzmx = lenpfac*neq
    lenplumx = nnzmx + lenplufac*neq
    ljac = 1
    lrownrms = nnzmx + ljac
    if (isrnorm .eq. 1) then
        lplu = lrownrms + neq
    else
        lplu = lrownrms
    end if
    lrwk1 = lplu + lenplumx
    lwp_min = lrwk1 + neq - 1
    if (lwp .lt. lwp_min) then
        ierr = 30
        return
    end if
    ljaci = 1
    ljacj = ljaci + neqp1
    lju = ljacj + nnzmx
    ljlu = lju + max(lenplumx, neqp1)
    if (ireorder .ne. 0) then
        lperm = ljlu + lenplumx
        lqperm = lperm + neq
        liwk1 = lqperm + neq
        llevels = ljlu + nnzmx
        lmask = llevels + neq
    else
        lperm = 0
        lqperm = 0
        llevels = 0
        lmask = 0
        liwk1 = ljlu + lenplumx
    end if
    liwp_min = liwk1 + 2*neq - 1
    if (ipremeth .eq. 2) then
        liwp_min = liwp_min + 2*neq
    end if
    if (liwp .lt. liwp_min) then
        ierr = 31
        return
    end if
    ierr = 0
    return
end subroutine dspsetup

subroutine djacilu(res, ires, neq, t, y, yprime, rewt, savr, wk, h, cj, wp, iwp, ierr, rpar, ipar)
    implicit none
    integer :: neq
    real(8) :: t
    real(8) :: y(neq)
    real(8) :: yprime(neq)
    real(8) :: savr(neq)
    real(8) :: rewt(neq)
    external :: res
    integer :: ires
    real(8) :: wk(neq)
    real(8) :: h
    real(8) :: cj
    real(8) :: rpar(*)
    integer :: ipar(*)
    real(8) :: wp(*)
    integer :: iwp(*)
    integer :: nre
    integer :: ierr
    real(8) :: tolilut, permtol, sqrtn
    integer :: i, lbw, ubw, lenplumx, ljac, ljaci, ljacj, liperm, lrownrms, lrwk1, liwk1, ifmt
    integer :: lenpfac, lenplufac, lfililut, ipremeth, neqp1, nnzmx
    integer :: lplu, lju, ljlu, lperm, lqperm, llevels, lmask
    integer :: isrnorm
    integer :: normtype
    integer :: ireorder
    integer :: jacout
    integer :: iunit
    integer :: jscalcol
    character(len=8) :: pmeth(4), premeth
    character(len=72) :: title
    character(len=80) :: msg
    save
    data pmeth/"ILUT", "ILUTP", "ILU0", "MILU0"/
    nre = 0
    lbw = ipar(1)
    ubw = ipar(2)
    lenpfac = ipar(3)
    lenplufac = ipar(4)
    ipremeth = ipar(5)
    lfililut = ipar(6)
    ireorder = ipar(7)
    isrnorm = ipar(8)
    normtype = ipar(9)
    jacout = ipar(10)
    jscalcol = ipar(11)
    tolilut = rpar(1)
    permtol = rpar(2)
    premeth = pmeth(ipremeth)
    neqp1 = neq + 1
    nnzmx = lenpfac*neq
    lenplumx = nnzmx + lenplufac*neq
    ljac = 1
    lrownrms = nnzmx + ljac
    if (isrnorm .eq. 1) then
        lplu = lrownrms + neq
    else
        lplu = lrownrms
    end if
    lrwk1 = lplu + lenplumx
    ljaci = 1
    ljacj = ljaci + neqp1
    lju = ljacj + nnzmx
    ljlu = lju + lenplumx
    ierr = 0
call djcalc(neq, t, y, yprime, savr, lbw, ubw, wk, rewt, res, h, cj, nnzmx, wp(ljac), iwp(ljacj), iwp(ljaci), wp(lplu), iwp(ljlu), iwp(lju), ipar, rpar, ires, nre, ierr)
    if (ires .lt. 0) then
        return
    end if
    if (ierr .ne. 0) then
        return
    end if
    ipar(30) = ipar(30) + nre
    ljlu = lju + neqp1
    if (ireorder .ne. 0) then
        lperm = ljlu + lenplumx
        lqperm = lperm + neq
        liwk1 = lqperm + neq
        llevels = ljlu + nnzmx
        lmask = llevels + neq
    else
        lperm = 0
        lqperm = 0
        llevels = 0
        lmask = 0
        liwk1 = ljlu + lenplumx
    end if
    if (premeth .eq. "ILUTP") then
        liperm = liwk1 + 2*neq
    else
        liperm = liwk1
    end if
    if (jscalcol .eq. 1) then
        sqrtn = sqrt(real(neq))
        do i = 1, neq
            wk(i) = sqrtn/rewt(i)
10          continue
        end do
        call amudia(neq, 0, wp(ljac), iwp(ljacj), iwp(ljaci), wk, wp(ljac), iwp(ljacj), iwp(ljaci))
    end if
    if (isrnorm .eq. 1) then
        call roscal(neq, 0, normtype, wp(ljac), iwp(ljacj), iwp(ljaci), wp(lrownrms), wp(ljac), iwp(ljacj), iwp(ljaci), ierr)
        if (ierr .ne. 0) then
            return
        end if
    end if
    if (ireorder .ne. 0) then
    call djreord(neq, neqp1, nnzmx, premeth, wp(ljac), iwp(ljacj), iwp(ljaci), wp(lplu), iwp(ljlu), iwp(lju), iwp(lperm), iwp(lqperm), iwp(llevels), iwp(lmask), ireorder)
    end if
    if (jacout .eq. 1) then
        iunit = ipar(29)
        if (isrnorm .eq. 1) then
            do i = 1, neq
                savr(i) = savr(i)*wp(lrownrms + i - 1)
20              continue
            end do
        end if
        if (ireorder .ne. 0) then
            call dvperm(neq, savr, iwp(lperm))
        end if
        title = " DDASPK Test Matrix "
        ifmt = 15
        call prtmt(neq, neq, wp(ljac), iwp(ljacj), iwp(ljaci), savr, "NN", title, "SPARSKIT", "RUA", ifmt, 3, iunit)
        msg = "DJACILU -- Jacobian Matrix written to file."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        ierr = 1
        ires = -2
        return
    end if
call djilu(neq, neq + 1, nnzmx, wp(ljac), iwp(ljacj), iwp(ljaci), iwp(lju), wp(lplu), iwp(ljlu), wp(lrwk1), iwp(liwk1), lenplumx, tolilut, lfililut, permtol, premeth, iwp(liperm), ierr)
    if (ierr .eq. -2 .or. ierr .eq. -3) then
        ires = -2
    end if
    ipar(21) = lplu
    ipar(22) = lju
    ipar(23) = ljlu
    ipar(24) = lrownrms
    ipar(25) = lperm
    ipar(26) = lqperm
    return
end subroutine djacilu

subroutine djcalc(neq, t, y, yprime, r0, ml, mu, r1, rewt, res, h, cj, nnzmx, jac, ja, ia, rcoo, jcoo, icoo, ipar, rpar, ires, nre, ierr)
    implicit none
    integer :: neq
    real(8) :: t
    real(8) :: y(neq)
    real(8) :: yprime(neq)
    real(8) :: r0(neq)
    real(8) :: rewt(neq)
    external :: res
    integer :: ires
    integer :: ml, mu
    integer :: nnzmx
    real(8) :: h
    real(8) :: cj
    real(8) :: rpar(*)
    integer :: ipar(*)
    real(8) :: r1(neq)
    real(8) :: jac(nnzmx)
    integer :: ja(nnzmx)
    integer :: ia(neq + 1)
    integer :: ierr
    real(8) :: rcoo(nnzmx)
    integer :: jcoo(nnzmx)
    integer :: icoo(nnzmx)
    integer :: nnz, i, i1, i2, j, jj, mba, meband, meb1, mband, nre
    real(8) :: jacelem, uround, d1mach, squr, del, delinv
    character(len=80) :: msg
    nnz = 1
    mband = ml + mu + 1
    mba = min0(mband, neq)
    meband = mband + ml
    meb1 = meband - 1
    uround = d1mach(4)
    squr = sqrt(uround)
    ierr = 0
    ires = 0
    do j = 1, mba
        do jj = j, neq, mband
            jac(jj) = y(jj)
            jac(jj + neq) = yprime(jj)
            del = squr*max(abs(y(jj)), abs(h*yprime(jj)), abs(1.0/rewt(jj)))
            del = sign(del, h*yprime(jj))
            del = y(jj) + del - y(jj)
            y(jj) = y(jj) + del
            yprime(jj) = yprime(jj) + cj*del
10          continue
        end do
        call res(t, y, yprime, cj, r1, ires, rpar, ipar)
        if (ires .lt. 0) then
            return
        end if
        nre = nre + 1
        do jj = j, neq, mband
            y(jj) = jac(jj)
            yprime(jj) = jac(jj + neq)
            del = squr*max(abs(y(jj)), abs(h*yprime(jj)), abs(1.0/rewt(jj)))
            del = sign(del, h*yprime(jj))
            del = y(jj) + del - y(jj)
            delinv = 1.0/del
            i1 = max(1, jj - mu)
            i2 = min(neq, jj + ml)
            do i = i1, i2
                jacelem = (r1(i) - r0(i))*delinv
                if (jacelem .ne. 0.) then
                    if (nnz .gt. nnzmx) then
                        msg = "DJCALC -- More storage needed for Jacobian."
                        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
                        msg = "DJCALC -- Increase LENPFAC."
                        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
                        msg = "DJCALC -- Storage exceeded at (I,J) = (I1,I2)"
                        call xerrwd(msg, 80, 0, 0, 2, i, jj, 0, 0.0, 0.0)
                        ierr = 1
                        ires = -2
                        return
                    end if
                    rcoo(nnz) = jacelem
                    jcoo(nnz) = jj
                    icoo(nnz) = i
                    nnz = nnz + 1
                end if
20              continue
            end do
30          continue
        end do
40      continue
    end do
    nnz = nnz - 1
    call coocsr(neq, nnz, rcoo, icoo, jcoo, jac, ja, ia)
    return
end subroutine djcalc

subroutine dpsolilu(neq, t, y, yprime, r0, wk, cj, wght, wp, iwp, bl, eplin, ierr, rpar, ipar)
    implicit none
    integer :: neq
    real(8) :: t
    real(8) :: y(neq)
    real(8) :: yprime(neq)
    real(8) :: r0(neq)
    real(8) :: wght(neq)
    real(8) :: cj
    real(8) :: eplin
    real(8) :: wp(*)
    integer :: iwp(*)
    integer :: ipar(*)
    real(8) :: rpar(*)
    real(8) :: wk(neq)
    real(8) :: bl(neq)
    integer :: ierr
    integer :: i, lplu, lju, ljlu, lrownrms, lperm, lqperm, ireorder, isrnorm, ipremeth, jscalcol
    ipremeth = ipar(5)
    ireorder = ipar(7)
    isrnorm = ipar(8)
    jscalcol = ipar(11)
    lplu = ipar(21)
    lju = ipar(22)
    ljlu = ipar(23)
    lrownrms = ipar(24)
    lperm = ipar(25)
    lqperm = ipar(26)
    if (isrnorm .eq. 1) then
        do i = 1, neq
            bl(i) = bl(i)*wp(lrownrms + i - 1)
10          continue
        end do
    end if
    if (ipremeth .eq. 1 .or. ipremeth .eq. 2) then
        if (ireorder .eq. 1) then
            call dvperm(neq, bl, iwp(lperm))
        end if
        call lusol(neq, bl, wk, wp(lplu), iwp(ljlu), iwp(lju))
        if (ireorder .eq. 1) then
            call dvperm(neq, wk, iwp(lqperm))
        end if
    end if
    if (jscalcol .eq. 1) then
        do i = 1, neq
40          bl(i) = wk(i)/wght(i)
        end do
    else
        do i = 1, neq
50          bl(i) = wk(i)
        end do
    end if
    ierr = 0
    return
end subroutine dpsolilu

subroutine djilu(neq, neqp1, nnzmx, jac, ja, ia, ju, plu, jlu, rwk1, iwk1, lenplumx, tolilut, lfililut, permtol, premeth, iperm, ierr)
    implicit none
    integer :: neq
    integer :: neqp1
    integer :: nnzmx
    real(8) :: jac(nnzmx)
    integer :: ja(nnzmx)
    integer :: ia(neqp1)
    character(len=8) :: premeth
    real(8) :: tolilut
    real(8) :: permtol
    integer :: lenplumx
    integer :: lfililut
    real(8) :: rwk1(neq)
    integer :: iwk1(2*neq), iperm(2*neq)
    real(8) :: plu(lenplumx)
    integer :: jlu(lenplumx)
    integer :: ju(neq)
    integer :: ierr
    character(len=80) :: msg
    logical :: error
    error = .false.
    if (premeth .eq. "ILUT") then
        call ilut(neq, jac, ja, ia, lfililut, tolilut, plu, jlu, ju, lenplumx, rwk1, iwk1, ierr)
        if (ierr .ne. 0) then
            msg = "DJILU -- Error return from ILUT: IERR = (I1)"
            call xerrwd(msg, 80, 0, 0, 1, ierr, 0, 0, 0.0, 0.0)
            error = .true.
        end if
    else
        if (premeth .eq. "ILUTP") then
            call ilutp(neq, jac, ja, ia, lfililut, tolilut, permtol, neq, plu, jlu, ju, lenplumx, rwk1, iwk1, iperm, ierr)
            if (ierr .ne. 0) then
                msg = "DJILU -- Error return from ILUTP: IERR = (I1)"
                call xerrwd(msg, 80, 0, 0, 1, ierr, 0, 0, 0.0, 0.0)
                error = .true.
            end if
        end if
    end if
    if (error) then
        msg = "DJILU -- IERR .NE. 0 means one of the following has occurred:"
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    IERR >  0   --> Zero pivot encountered at step number IERR."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    IERR = -1   --> Error. input matrix may be wrong."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "                     (The elimination process has generated a"
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "                     row in L or U with length > NEQ.)"
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    IERR = -2   --> Matrix L overflows."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    IERR = -3   --> Matrix U overflows."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    IERR = -4   --> Illegal value for LFILILUT."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    IERR = -5   --> Zero row encountered."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    "
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    For IERR = -2 or -3, increase the value of LENPLUFAC or"
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    decrease the value of LFILILUT if LENPLUFAC cannot be"
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        msg = "    increased."
        call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
    end if
    return
end subroutine djilu

subroutine djreord(neq, neqp1, nnzmx, premeth, jac, ja, ia, awk, jwk, iwk, perm, qperm, levels, mask, ireorder)
    implicit none
    integer :: neq
    integer :: neqp1
    integer :: nnzmx
    real(8) :: jac(nnzmx)
    integer :: ja(nnzmx)
    integer :: ia(neqp1)
    character(len=8) :: premeth
    real(8) :: awk(nnzmx)
    integer :: jwk(nnzmx)
    integer :: iwk(neqp1)
    integer :: perm(neq)
    integer :: qperm(neq)
    integer :: levels(neq)
    integer :: mask(neq)
    integer :: ireorder
    integer :: nlev
    integer :: maskval
    integer :: i, nfirst
    if (ireorder .eq. 1) then
        call atob(neq, jac, ja, ia, awk, jwk, iwk)
        nfirst = 1
        perm(1) = 0
        do i = 1, neq
            mask(i) = 1
        end do
        maskval = 1
        qperm(1) = 1
        call bfs(neq, jwk, iwk, nfirst, perm, mask, maskval, qperm, levels, nlev)
        call rversp(neq, qperm)
        do i = 1, neq
            perm(qperm(i)) = i
        end do
        call dperm(neq, awk, jwk, iwk, jac, ja, ia, perm, perm, 1)
    end if
    return
end subroutine djreord

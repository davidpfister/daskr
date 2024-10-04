subroutine dmset2(mx, my, ns, nsd, lid, iwork)
    implicit double precision(a - h, o - z)
    dimension :: iwork(*)
    common/drpre1/srur, mp, mpd, mpsq, meshx, meshy, mxmp
    uround = d1mach(4)
    srur = sqrt(uround)
    mp = ns
    mpd = nsd
    mpsq = ns*ns
    meshx = mx
    meshy = my
    mxmp = meshx*mp
    iwork(27) = mpsq*meshx*meshy
    iwork(28) = mp*meshx*meshy
    if (lid .eq. 0) then
        return
    end if
    i0 = lid
    do jy = 1, my
        do jx = 1, mx
            do i = 1, mpd
10              iwork(i0 + i) = 1
            end do
            do i = mpd + 1, mp
20              iwork(i0 + i) = -1
            end do
            i0 = i0 + mp
30          continue
        end do
40      continue
    end do
    return
end subroutine dmset2

subroutine drbdja(t, u, r0, rblock, r1, rewt, cj, bd, ipbd, ier)
    implicit double precision(a - h, o - z)
    external :: rblock
    dimension :: u(*), r0(*), r1(*), rewt(*), bd(*), ipbd(*)
    common/drpre1/srur, mp, mpd, mpsq, meshx, meshy, mxmp
    dfac = 1.0d-2
    ibd = 0
    j0 = 0
    do jy = 1, meshy
        do jx = 1, meshx
            do js = 1, mp
                j = j0 + js
                uj = u(j)
                del = max(srur*abs(uj), dfac/rewt(j))
                u(j) = u(j) + del
                fac = (-1.0d0)/del
                call rblock(t, jx, jy, u(j0 + 1), r1)
                do i = 1, mp
10                  bd(ibd + i) = (r1(i) - r0(j0 + i))*fac
                end do
                u(j) = uj
                ibd = ibd + mp
20              continue
            end do
            j0 = j0 + mp
30          continue
        end do
40      continue
    end do
    ibd = 1
    iip = 1
    do j = 1, meshx*meshy
        idiag = ibd
        do i = 1, mp
            if (i .le. mpd) then
                bd(idiag) = bd(idiag) + cj
            end if
70          idiag = idiag + mp + 1
        end do
        call dgefa(bd(ibd), mp, mp, ipbd(iip), ier)
        if (ier .ne. 0) then
            go to 90
        end if
        ibd = ibd + mpsq
        iip = iip + mp
80      continue
    end do
90  return
end subroutine drbdja

subroutine drbdps(b, bd, ipbd)
    implicit double precision(a - h, o - z)
    dimension :: b(*), bd(*), ipbd(*)
    common/drpre1/srur, mp, mpd, mpsq, meshx, meshy, mxmp
    ier = 0
    ib = 1
    ibd = 1
    do jy = 1, meshy
        do jx = 1, meshx
            call dgesl(bd(ibd), mp, mp, ipbd(ib), b(ib), 0)
            ib = ib + mp
            ibd = ibd + mpsq
10          continue
        end do
20      continue
    end do
    return
end subroutine drbdps

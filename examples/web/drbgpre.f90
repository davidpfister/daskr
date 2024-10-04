subroutine dgset2(mx, my, ns, nsd, nxg, nyg, lid, iwork)
    implicit double precision(a - h, o - z)
    dimension :: iwork(*)
    parameter(maxm=50)
    common/drpre1/srur, mp, mpd, mpsq, meshx, meshy, mxmp
    common/rpre2/ngx, ngy, ngrp, jgx(maxm + 1), jgy(maxm + 1), jigx(maxm), jigy(maxm), jxr(maxm), jyr(maxm)
    uround = d1mach(4)
    srur = sqrt(uround)
    mp = ns
    mpd = nsd
    mpsq = ns*ns
    meshx = mx
    meshy = my
    mxmp = meshx*mp
    ngx = nxg
    ngy = nyg
    ngrp = ngx*ngy
    call gset1(meshx, ngx, jgx, jigx, jxr)
    call gset1(meshy, ngy, jgy, jigy, jyr)
    iwork(27) = mpsq*ngrp
    iwork(28) = mp*ngrp
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
end subroutine dgset2

subroutine gset1(m, ng, jg, jig, jr)
    dimension :: jg(*), jig(*), jr(*)
    mper = m/ng
    do ig = 1, ng
10      jg(ig) = 1 + (ig - 1)*mper
    end do
    jg(ng + 1) = m + 1
    ngm1 = ng - 1
    len1 = ngm1*mper
    do j = 1, len1
20      jig(j) = 1 + (j - 1)/mper
    end do
    len1 = len1 + 1
    do j = len1, m
25      jig(j) = ng
    end do
    do ig = 1, ngm1
30      jr(ig) = 0.5d0 + (real(ig) - 0.5d0)*real(mper)
    end do
    jr(ng) = 0.5d0*real(1 + ngm1*mper + m)
    return
end subroutine gset1

subroutine drbgja(t, u, r0, rblock, r1, rewt, cj, bd, ipbd, ier)
    implicit double precision(a - h, o - z)
    external :: rblock
    dimension :: u(*), r0(*), r1(*), rewt(*), bd(*), ipbd(*)
    parameter(maxm=50)
    common/drpre1/srur, mp, mpd, mpsq, meshx, meshy, mxmp
    common/rpre2/ngx, ngy, ngrp, jgx(maxm + 1), jgy(maxm + 1), jigx(maxm), jigy(maxm), jxr(maxm), jyr(maxm)
    dfac = 1.0d-2
    ibd = 0
    do igy = 1, ngy
        jy = jyr(igy)
        j00 = (jy - 1)*mxmp
        do igx = 1, ngx
            jx = jxr(igx)
            j0 = j00 + (jx - 1)*mp
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
30          continue
        end do
40      continue
    end do
    ibd = 1
    iip = 1
    do ig = 1, ngrp
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
end subroutine drbgja

subroutine drbgps(b, bd, ipbd)
    implicit double precision(a - h, o - z)
    dimension :: b(*), bd(*), ipbd(*)
    parameter(maxm=50)
    common/drpre1/srur, mp, mpd, mpsq, meshx, meshy, mxmp
    common/rpre2/ngx, ngy, ngrp, jgx(maxm + 1), jgy(maxm + 1), jigx(maxm), jigy(maxm), jxr(maxm), jyr(maxm)
    ier = 0
    ib = 1
    do jy = 1, meshy
        igy = jigy(jy)
        ig0 = (igy - 1)*ngx
        do jx = 1, meshx
            igx = jigx(jx)
            igm1 = igx - 1 + ig0
            ibd = 1 + igm1*mpsq
            iip = 1 + igm1*mp
            call dgesl(bd(ibd), mp, mp, ipbd(iip), b(ib), 0)
            ib = ib + mp
10          continue
        end do
20      continue
    end do
    return
end subroutine drbgps

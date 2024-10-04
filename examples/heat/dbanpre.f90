subroutine dbanja(res, ires, neq, t, y, yprime, rewt, savr, wk, h, cj, wp, iwp, ier, rpar, ipar)
    implicit double precision(a - h, o - z)
    external :: res
    dimension :: y(*), yprime(*), savr(*), rewt(*), wk(*)
    dimension :: wp(*), iwp(*), rpar(*), ipar(*)
    ml = ipar(1)
    mu = ipar(2)
    mband = ml + mu + 1
    mba = min(mband, neq)
    meband = mband + ml
    meb1 = meband - 1
    uround = d1mach(4)
    squr = sqrt(uround)
    lenp = (2*ml + mu + 1)*neq
    msave = neq/mband + 1
    isave = lenp
    ipsave = isave + msave
    ier = 0
    ires = 0
    do j = 1, mba
        do n = j, neq, mband
            k = (n - j)/mband + 1
            wp(isave + k) = y(n)
            wp(ipsave + k) = yprime(n)
            del = squr*max(abs(y(n)), abs(h*yprime(n)), abs(1.d0/rewt(n)))
            del = sign(del, h*yprime(n))
            del = y(n) + del - y(n)
            y(n) = y(n) + del
            yprime(n) = yprime(n) + cj*del
10          continue
        end do
        call res(t, y, yprime, cj, wk, ires, rpar, ipar)
        if (ires .lt. 0) then
            return
        end if
        do n = j, neq, mband
            k = (n - j)/mband + 1
            y(n) = wp(isave + k)
            yprime(n) = wp(ipsave + k)
            del = squr*max(abs(y(n)), abs(h*yprime(n)), abs(1.d0/rewt(n)))
            del = sign(del, h*yprime(n))
            del = y(n) + del - y(n)
            delinv = 1.0d0/del
            i1 = max(1, n - mu)
            i2 = min(neq, n + ml)
            ii = n*meb1 - ml
            do i = i1, i2
20              wp(ii + i) = (wk(i) - savr(i))*delinv
            end do
30          continue
        end do
40      continue
    end do
    call dgbfa(wp, meband, neq, ml, mu, iwp, ier)
    return
end subroutine dbanja

subroutine dbanps(neq, t, y, yprime, savr, wk, cj, wght, wp, iwp, b, eplin, ier, rpar, ipar)
    implicit double precision(a - h, o - z)
    dimension :: b(*), wp(*), iwp(*), rpar(*), ipar(*)
    ml = ipar(1)
    mu = ipar(2)
    meband = 2*ml + mu + 1
    call dgbsl(wp, meband, neq, ml, mu, iwp, b, 0)
    return
end subroutine dbanps

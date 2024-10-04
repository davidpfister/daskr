double precision function d1mach(idummy)
    integer :: idummy
    double precision :: u, comp
    u = 1.0d0
10  u = u*0.5d0
    call dumsum(1.0d0, u, comp)
    if (comp .ne. 1.0d0) then
        go to 10
    end if
    d1mach = u*2.0d0
    return
end function d1mach

subroutine dumsum(a, b, c)
    double precision :: a, b, c
    c = a + b
    return
end subroutine dumsum

subroutine xerrwd(msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
    double precision :: r1, r2
    integer :: nmes, nerr, level, ni, i1, i2, nr
    character*(*) :: msg
    integer :: lunit, ixsav, mesflg
    lunit = ixsav(1, 0, .false.)
    mesflg = ixsav(2, 0, .false.)
    if (mesflg .eq. 0) then
        go to 100
    end if
    write (lunit, 10) msg
10  format(1x, a)
    if (ni .eq. 1) then
        write (lunit, 20) i1
    end if
20  format(6x, 'In above message,  I1 =', i10)
    if (ni .eq. 2) then
        write (lunit, 30) i1, i2
    end if
30  format(6x, 'In above message,  I1 =', i10, 3x, 'I2 =', i10)
    if (nr .eq. 1) then
        write (lunit, 40) r1
    end if
40  format(6x, 'In above message,  R1 =', d21.13)
    if (nr .eq. 2) then
        write (lunit, 50) r1, r2
    end if
50  format(6x, 'In above,  R1 =', d21.13, 3x, 'R2 =', d21.13)
100 if (level .ne. 2) then
        return
    end if
    stop
end subroutine xerrwd

subroutine xsetf(mflag)
    integer :: mflag, junk, ixsav
    if (mflag .eq. 0 .or. mflag .eq. 1) then
        junk = ixsav(2, mflag, .true.)
    end if
    return
end subroutine xsetf

subroutine xsetun(lun)
    integer :: lun, junk, ixsav
    if (lun .gt. 0) then
        junk = ixsav(1, lun, .true.)
    end if
    return
end subroutine xsetun

integer function ixsav(ipar, ivalue, iset)
    logical :: iset
    integer :: ipar, ivalue
    integer :: lunit, lundef, mesflg
    save :: lunit, lundef, mesflg
    data lunit/-1/, lundef/6/, mesflg/1/
    if (ipar .eq. 1) then
        if (lunit .eq. -1) then
            lunit = lundef
        end if
        ixsav = lunit
        if (iset) then
            lunit = ivalue
        end if
    end if
    if (ipar .eq. 2) then
        ixsav = mesflg
        if (iset) then
            mesflg = ivalue
        end if
    end if
    return
end function ixsav

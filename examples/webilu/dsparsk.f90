subroutine prtmt(nrow, ncol, a, ja, ia, rhs, guesol, title, key, type, &
                 ifmt, job, iounit)

    !*****************************************************************************80
    !
   !! PRTMT writes a matrix in Harwell-Boeing format into a file.
    !
    !  Discussion:
    !
    !    This routine assumes that the matrix is stored in CSC format
    !    (Compressed Sparse Column format).
    !    There is some limited functionality for right hand sides.
    !
    !    This code attempts to pack as many elements as possible per
    !    80-character line.
    !
    !    This code attempts to avoid as much as possible to put
    !    blanks in the formats that are written in the 4-line header
    !    This is done for purely esthetical reasons since blanks
    !    are ignored in format descriptors.
    !
    !    sparse formats for right hand sides and guesses not supported.
    !
    !  Modified:
    !
    !    07 January 2004
    !
    !  Author:
    !
    !    Youcef Saad
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
    !
    !    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
    !
    !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NCOL+1), the matrix in CSC
    !    Compressed Sparse Column format.
    !
    !    Input, real RHS(*), contains the right hand sides and optionally
    !    the associated initial guesses and/or exact solutions
    !    in this order.  See also GUESOL for details.   RHS will
    !    be used only if 2 < JOB.  Only full storage for the right hand
    !    sides is supported.
    !
    ! guesol = a 2-character string indicating whether an initial guess
    !          (1-st character) and / or the exact solution (2-nd)
    !          character) is provided with the right hand side.
    !         if the first character of guesol is 'G' it means that an
    !          an intial guess is provided for each right hand sides.
    !          These are assumed to be appended to the right hand sides in
    !          the array rhs.
    !         if the second character of guesol is 'X' it means that an
    !          exact solution is provided for each right hand side.
    !          These are assumed to be appended to the right hand sides
    !          and the initial guesses (if any) in the array rhs.
    !
    ! title  = character*71 = title of matrix test ( character a*71 ).
    ! key    = character*8  = key of matrix
    ! type   = charatcer*3  = type of matrix.
    !
    ! ifmt       = integer ( kind = 4 ) specifying the format chosen for the real values
    !         to be output (i.e., for a, and for rhs-guess-sol if
    !          applicable). the meaning of ifmt is as follows.
    !        * if (ifmt < 100) then the E descriptor is used,
    !           format Ed.m, in which the length (m) of the mantissa is
    !           precisely the integer ( kind = 4 ) ifmt (and d = ifmt+6)
    !        * if (ifmt > 100) then prtmt will use the
    !           F- descriptor (format Fd.m) in which the length of the
    !           mantissa (m) is the integer ( kind = 4 ) mod(ifmt,100) and the length
    !           of the integer ( kind = 4 ) part is k = ifmt/100 (and d = k+m+2)
    !          Thus  ifmt= 4   means  E10.4  +.xxxxD+ee    while
    !                ifmt=104  means  F7.4   +x.xxxx
    !                ifmt=205  means  F9.5   +xx.xxxxx
    !          Note: formats for ja, and ia are internally computed.
    !
    ! job       = integer ( kind = 4 ) to indicate whether matrix values and
    !         a right hand side is available to be written
    !          job = 1   write srtucture only, i.e., the arrays ja and ia.
    !          job = 2   write matrix including values, i.e., a, ja, ia
    !          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
    !         job = nrhs+2 write matrix and nrhs successive right hand sides
    !         Note that there cannot be any right hand side if the matrix
    !         has no values. Also the initial guess and exact solutions when
    !          provided are for each right hand side. For example if nrhs=2
    !          and guesol='GX' there are 6 vectors to write.
    !
    !
    ! iounit = logical unit number where to write the matrix into.
    !
    ! on return:
    !
    ! the matrix a, ja, ia will be written in output unit iounit
    ! in the Harwell-Boeing format. Noe of the inputs is modofied.
    !
    implicit none

    integer(kind=4) ncol

    real(kind=8) a(*)
    character(len=2) guesol
    integer(kind=4) i
    integer(kind=4) ia(ncol + 1)
    integer(kind=4) iend
    integer(kind=4) ifmt
    integer(kind=4) ihead
    integer(kind=4) indcrd
    character(len=16) indfmt
    integer(kind=4) iounit
    integer(kind=4) ix
    integer(kind=4) ja(*)
    integer(kind=4) job
    character(len=8) key
    integer(kind=4) len
    integer(kind=4) next
    integer(kind=4) nnz
    integer(kind=4) nperli
    integer(kind=4) nrhs
    integer(kind=4) nrow
    integer(kind=4) ptrcrd
    character(len=16) ptrfmt
    real(kind=8) rhs(*)
    integer(kind=4) rhscrd
    character(len=3) rhstyp
    character(len=72) title
    integer(kind=4) totcrd
    character(len=3) type
    integer(kind=4) valcrd
    character(len=20) valfmt
    !
    !  Compute pointer format.
    !
    nnz = ia(ncol + 1) - 1
    len = int(log10(0.1d+00 + real(nnz + 1, kind=8))) + 1
    nperli = 80/len
    ptrcrd = ncol/nperli + 1

    if (9 .lt. len) then
        assign 101 to ix
    else
        assign 100 to ix
    end if

    write (ptrfmt, ix) nperli, len
100 format(1h(, i2, 1hi, i1, 1h))
101 format(1h(, i2, 1hi, i2, 1h))
    !
    !  Compute the ROW index format.
    !
    len = int(log10(0.1d+00 + real(nrow, kind=8))) + 1
    nperli = min(80/len, nnz)
    indcrd = (nnz - 1)/nperli + 1
    write (indfmt, 100) nperli, len
    !
    !  Compute values and RHS format (using the same for both).
    !
    valcrd = 0
    rhscrd = 0
    !
    !  Skip this part if no values provided.
    !
    if (job .le. 1) then
        go to 20
    end if

    if (100 .le. ifmt) then
        ihead = ifmt/100
        ifmt = ifmt - 100*ihead
        len = ihead + ifmt + 2
        nperli = 80/len

        if (len .le. 9) then
            assign 102 to ix
        elseif (ifmt .le. 9) then
            assign 103 to ix
        else
            assign 104 to ix
        end if

        write (valfmt, ix) nperli, len, ifmt
102     format(1h(, i2, 1hf, i1, 1h., i1, 1h))
103     format(1h(, i2, 1hf, i2, 1h., i1, 1h))
104     format(1h(, i2, 1hf, i2, 1h., i2, 1h))

    else
        len = ifmt + 6
        nperli = 80/len
        !
        !  Try to minimize the blanks in the format strings.
        !
        if (nperli .le. 9) then
            if (len .le. 9) then
                assign 105 to ix
            else if (ifmt .le. 9) then
                assign 106 to ix
            else
                assign 107 to ix
            end if
        else
            if (len .le. 9) then
                assign 108 to ix
            else if (ifmt .le. 9) then
                assign 109 to ix
            else
                assign 110 to ix
            end if
        end if

        write (valfmt, ix) nperli, len, ifmt
105     format(1h(, i1, 1he, i1, 1h., i1, 1h))
106     format(1h(, i1, 1he, i2, 1h., i1, 1h))
107     format(1h(, i1, 1he, i2, 1h., i2, 1h))
108     format(1h(, i2, 1he, i1, 1h., i1, 1h))
109     format(1h(, i2, 1he, i2, 1h., i1, 1h))
110     format(1h(, i2, 1he, i2, 1h., i2, 1h))

    end if
    valcrd = (nnz - 1)/nperli + 1
    nrhs = job - 2

    if (1 .le. nrhs) then
        i = (nrhs*nrow - 1)/nperli + 1
        rhscrd = i
        if (guesol(1:1) .eq. 'G') then
            rhscrd = rhscrd + i
        end if
        if (guesol(2:2) .eq. 'X') then
            rhscrd = rhscrd + i
        end if
        rhstyp = 'F'//guesol
    end if

20  continue

    totcrd = ptrcrd + indcrd + valcrd + rhscrd
    !
    !  Write four line or five line header.
    !
    write (iounit, 10) title, key, totcrd, ptrcrd, indcrd, valcrd, &
        rhscrd, type, nrow, ncol, nnz, nrhs, ptrfmt, indfmt, valfmt, valfmt

    if (1 .le. nrhs) then
        write (iounit, 11) rhstyp, nrhs
    end if

10  format(a72, a8/5i14/a3, 11x, 4i14/2a16, 2a20)
11  format(A3, 11x, i4)

    write (iounit, ptrfmt) ia(1:ncol + 1)
    write (iounit, indfmt) ja(1:nnz)

    if (job .le. 1) then
        return
    end if

    write (iounit, valfmt) (a(i), i=1, nnz)
    if (job .le. 2) then
        return
    end if
    len = nrow*nrhs
    next = 1
    iend = len
    write (iounit, valfmt) (rhs(i), i=next, iend)
    !
    !  Write initial guesses if available
    !
    if (guesol(1:1) .eq. 'G') then
        next = next + len
        iend = iend + len
        write (iounit, valfmt) (rhs(i), i=next, iend)
    end if
    !
    !  Write exact solutions if available.
    !
    if (guesol(2:2) .eq. 'X') then
        next = next + len
        iend = iend + len
        write (iounit, valfmt) (rhs(i), i=next, iend)
    end if

    return
end subroutine

subroutine aplb(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, iw, ierr)
    real(8) :: a(*), b(*), c(*)
    integer :: ja(*), jb(*), jc(*), ia(nrow + 1), ib(nrow + 1), ic(nrow + 1), iw(ncol)
    logical :: values
    values = job .ne. 0
    ierr = 0
    len = 0
    ic(1) = 1
    do j = 1, ncol
        iw(j) = 0
1       continue
    end do
    do ii = 1, nrow
        do ka = ia(ii), ia(ii + 1) - 1
            len = len + 1
            jcol = ja(ka)
            if (len .gt. nzmax) then
                go to 999
            end if
            jc(len) = jcol
            if (values) then
                c(len) = a(ka)
            end if
            iw(jcol) = len
200         continue
        end do
        do kb = ib(ii), ib(ii + 1) - 1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
                len = len + 1
                if (len .gt. nzmax) then
                    go to 999
                end if
                jc(len) = jcol
                if (values) then
                    c(len) = b(kb)
                end if
                iw(jcol) = len
            else
                if (values) then
                    c(jpos) = c(jpos) + b(kb)
                end if
            end if
300         continue
        end do
        do k = ic(ii), len
            iw(jc(k)) = 0
301         continue
        end do
        ic(ii + 1) = len + 1
500     continue
    end do
    return
999 ierr = ii
    return
end subroutine aplb

subroutine aplb1(nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, ierr)
    real(8) :: a(*), b(*), c(*)
    integer :: ja(*), jb(*), jc(*), ia(nrow + 1), ib(nrow + 1), ic(nrow + 1)
    logical :: values
    values = job .ne. 0
    ierr = 0
    kc = 1
    ic(1) = kc
    do i = 1, nrow
        ka = ia(i)
        kb = ib(i)
        kamax = ia(i + 1) - 1
        kbmax = ib(i + 1) - 1
5       continue
        if (ka .le. kamax) then
            j1 = ja(ka)
        else
            j1 = ncol + 1
        end if
        if (kb .le. kbmax) then
            j2 = jb(kb)
        else
            j2 = ncol + 1
        end if
        if (j1 .eq. j2) then
            if (values) then
                c(kc) = a(ka) + b(kb)
            end if
            jc(kc) = j1
            ka = ka + 1
            kb = kb + 1
            kc = kc + 1
        else
            if (j1 .lt. j2) then
                jc(kc) = j1
                if (values) then
                    c(kc) = a(ka)
                end if
                ka = ka + 1
                kc = kc + 1
            else
                if (j1 .gt. j2) then
                    jc(kc) = j2
                    if (values) then
                        c(kc) = b(kb)
                    end if
                    kb = kb + 1
                    kc = kc + 1
                end if
            end if
        end if
        if (kc .gt. nzmax) then
            go to 999
        end if
        if (ka .le. kamax .or. kb .le. kbmax) then
            go to 5
        end if
        ic(i + 1) = kc
6       continue
    end do
    return
999 ierr = i
    return
end subroutine aplb1

subroutine aplsb(nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, nzmax, ierr)
    real(8) :: a(*), b(*), c(*), s
    integer :: ja(*), jb(*), jc(*), ia(nrow + 1), ib(nrow + 1), ic(nrow + 1)
    ierr = 0
    kc = 1
    ic(1) = kc
    do i = 1, nrow
        ka = ia(i)
        kb = ib(i)
        kamax = ia(i + 1) - 1
        kbmax = ib(i + 1) - 1
5       continue
        if (ka .le. kamax .or. kb .le. kbmax) then
            if (ka .le. kamax) then
                j1 = ja(ka)
            else
                j1 = ncol + 1
            end if
            if (kb .le. kbmax) then
                j2 = jb(kb)
            else
                j2 = ncol + 1
            end if
            if (j1 .eq. j2) then
                c(kc) = a(ka) + s*b(kb)
                jc(kc) = j1
                ka = ka + 1
                kb = kb + 1
                kc = kc + 1
            else
                if (j1 .lt. j2) then
                    jc(kc) = j1
                    c(kc) = a(ka)
                    ka = ka + 1
                    kc = kc + 1
                else
                    if (j1 .gt. j2) then
                        jc(kc) = j2
                        c(kc) = s*b(kb)
                        kb = kb + 1
                        kc = kc + 1
                    end if
                end if
            end if
            if (kc .gt. nzmax) then
                go to 999
            end if
            go to 5
        end if
        ic(i + 1) = kc
6       continue
    end do
    return
999 ierr = i
    return
end subroutine aplsb

subroutine diamua(nrow, job, a, ja, ia, diag, b, jb, ib)
    real(8) :: a(*), b(*), diag(nrow), scal
    integer :: ja(*), jb(*), ia(nrow + 1), ib(nrow + 1)
    do ii = 1, nrow
        k1 = ia(ii)
        k2 = ia(ii + 1) - 1
        scal = diag(ii)
        do k = k1, k2
            b(k) = a(k)*scal
2           continue
        end do
1       continue
    end do
    if (job .eq. 0) then
        return
    end if
    do ii = 1, nrow + 1
        ib(ii) = ia(ii)
3       continue
    end do
    do k = ia(1), ia(nrow + 1) - 1
        jb(k) = ja(k)
31      continue
    end do
    return
end subroutine diamua

subroutine amudia(nrow, job, a, ja, ia, diag, b, jb, ib)
    real(8) :: a(*), b(*), diag(nrow)
    integer :: ja(*), jb(*), ia(nrow + 1), ib(nrow + 1)
    do ii = 1, nrow
        k1 = ia(ii)
        k2 = ia(ii + 1) - 1
        do k = k1, k2
            b(k) = a(k)*diag(ja(k))
2           continue
        end do
1       continue
    end do
    if (job .eq. 0) then
        return
    end if
    do ii = 1, nrow + 1
        ib(ii) = ia(ii)
3       continue
    end do
    do k = ia(1), ia(nrow + 1) - 1
        jb(k) = ja(k)
31      continue
    end do
    return
end subroutine amudia

subroutine aplsca(nrow, a, ja, ia, scal, iw)
    real(8) :: a(*), scal
    integer :: ja(*), ia(nrow + 1), iw(*)
    logical :: test
    call diapos(nrow, ja, ia, iw)
    icount = 0
    do j = 1, nrow
        if (iw(j) .eq. 0) then
            icount = icount + 1
        else
            a(iw(j)) = a(iw(j)) + scal
        end if
1       continue
    end do
    if (icount .eq. 0) then
        return
    end if
    ko = ia(nrow + 1) + icount
    do ii = nrow, 1, -1
        k1 = ia(ii)
        k2 = ia(ii + 1) - 1
        ia(ii + 1) = ko
        test = iw(ii) .eq. 0
        do k = k2, k1, -1
            j = ja(k)
            if (test .and. j .lt. ii) then
                test = .false.
                ko = ko - 1
                a(ko) = scal
                ja(ko) = ii
                iw(ii) = ko
            end if
            ko = ko - 1
            a(ko) = a(k)
            ja(ko) = j
4           continue
        end do
        if (test) then
            ko = ko - 1
            a(ko) = scal
            ja(ko) = ii
            iw(ii) = ko
        end if
5       continue
    end do
    ia(1) = ko
    return
end subroutine aplsca

subroutine amux(n, x, y, a, ja, ia)
    real(8) :: x(*), y(*), a(*)
    integer :: n, ja(*), ia(*)
    real(8) :: t
    integer :: i, k
    do i = 1, n
        t = 0.0d0
        do k = ia(i), ia(i + 1) - 1
            t = t + a(k)*x(ja(k))
99          continue
        end do
        y(i) = t
100     continue
    end do
    return
end subroutine amux

subroutine csrdns(nrow, ncol, a, ja, ia, dns, ndns, ierr)
    real(8) :: dns(ndns, *), a(*)
    integer :: ja(*), ia(*)
    ierr = 0
    do i = 1, nrow
        do j = 1, ncol
            dns(i, j) = 0.0d0
2           continue
        end do
1       continue
    end do
    do i = 1, nrow
        do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            if (j .gt. ncol) then
                ierr = i
                return
            end if
            dns(i, j) = a(k)
3           continue
        end do
4       continue
    end do
    return
end subroutine csrdns

subroutine coocsr(nrow, nnz, a, ir, jc, ao, jao, iao)
    real(8) :: a(*), ao(*), x
    integer :: ir(*), jc(*), jao(*), iao(*)
    do k = 1, nrow + 1
        iao(k) = 0
1       continue
    end do
    do k = 1, nnz
        iao(ir(k)) = iao(ir(k)) + 1
2       continue
    end do
    k = 1
    do j = 1, nrow + 1
        k0 = iao(j)
        iao(j) = k
        k = k + k0
3       continue
    end do
    do k = 1, nnz
        i = ir(k)
        j = jc(k)
        x = a(k)
        iad = iao(i)
        ao(iad) = x
        jao(iad) = j
        iao(i) = iad + 1
4       continue
    end do
    do j = nrow, 1, -1
        iao(j + 1) = iao(j)
5       continue
    end do
    iao(1) = 1
    return
end subroutine coocsr

subroutine coicsr(n, nnz, job, a, ja, ia, iwk)
    integer :: ia(nnz), ja(nnz), iwk(n + 1)
    real(8) :: a(*)
    real(8) :: t, tnext
    logical :: values
    values = job .eq. 1
    do i = 1, n + 1
        iwk(i) = 0
35      continue
    end do
    do k = 1, nnz
        i = ia(k)
        iwk(i + 1) = iwk(i + 1) + 1
4       continue
    end do
    iwk(1) = 1
    do i = 2, n
        iwk(i) = iwk(i - 1) + iwk(i)
44      continue
    end do
    init = 1
    k = 0
5   if (values) then
        t = a(init)
    end if
    i = ia(init)
    j = ja(init)
    ia(init) = -1
6   k = k + 1
    ipos = iwk(i)
    if (values) then
        tnext = a(ipos)
    end if
    inext = ia(ipos)
    jnext = ja(ipos)
    if (values) then
        a(ipos) = t
    end if
    ja(ipos) = j
    iwk(i) = ipos + 1
    if (ia(ipos) .lt. 0) then
        go to 65
    end if
    t = tnext
    i = inext
    j = jnext
    ia(ipos) = -1
    if (k .lt. nnz) then
        go to 6
    end if
    go to 70
65  init = init + 1
    if (init .gt. nnz) then
        go to 70
    end if
    if (ia(init) .lt. 0) then
        go to 65
    end if
    go to 5
70  do i = 1, n
        ia(i + 1) = iwk(i)
80      continue
    end do
    ia(1) = 1
    return
end subroutine coicsr

subroutine csrcoo(nrow, job, nzmax, a, ja, ia, nnz, ao, ir, jc, ierr)
    real(8) :: a(*), ao(*)
    integer :: ir(*), jc(*), ja(*), ia(nrow + 1)
    ierr = 0
    nnz = ia(nrow + 1) - 1
    if (nnz .gt. nzmax) then
        ierr = 1
        return
    end if
    go to(3, 2, 1), job
1   do k = 1, nnz
        ao(k) = a(k)
10      continue
    end do
2   do k = 1, nnz
        jc(k) = ja(k)
11      continue
    end do
3   do i = nrow, 1, -1
        k1 = ia(i + 1) - 1
        k2 = ia(i)
        do k = k1, k2, -1
            ir(k) = i
12          continue
        end do
13      continue
    end do
    return
end subroutine csrcoo

subroutine csrcsc(n, job, ipos, a, ja, ia, ao, jao, iao)
    integer :: ia(n + 1), iao(n + 1), ja(*), jao(*)
    real(8) :: a(*), ao(*)
    call csrcsc2(n, n, job, ipos, a, ja, ia, ao, jao, iao)
end subroutine csrcsc

subroutine csrcsc2(n, n2, job, ipos, a, ja, ia, ao, jao, iao)
    integer :: ia(n + 1), iao(n2 + 1), ja(*), jao(*)
    real(8) :: a(*), ao(*)
    do i = 1, n2 + 1
        iao(i) = 0
1       continue
    end do
    do i = 1, n
        do k = ia(i), ia(i + 1) - 1
            j = ja(k) + 1
            iao(j) = iao(j) + 1
2           continue
        end do
3       continue
    end do
    iao(1) = ipos
    do i = 1, n2
        iao(i + 1) = iao(i) + iao(i + 1)
4       continue
    end do
    do i = 1, n
        do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            next = iao(j)
            if (job .eq. 1) then
                ao(next) = a(k)
            end if
            jao(next) = i
            iao(j) = next + 1
62          continue
        end do
6       continue
    end do
    do i = n2, 1, -1
        iao(i + 1) = iao(i)
7       continue
    end do
    iao(1) = ipos
end subroutine csrcsc2

subroutine csrdia(n, idiag, job, a, ja, ia, ndiag, diag, ioff, ao, jao, iao, ind)
    real(8) :: diag(ndiag, idiag), a(*), ao(*)
    integer :: ia(*), ind(*), ja(*), jao(*), iao(*), ioff(*)
    job1 = job/10
    job2 = job - job1*10
    if (job1 .eq. 0) then
        go to 50
    end if
    n2 = n + n - 1
    call infdia(n, ja, ia, ind, idum)
    ii = 0
4   ii = ii + 1
    jmax = 0
    do k = 1, n2
        j = ind(k)
        if (j .le. jmax) then
            go to 41
        end if
        i = k
        jmax = j
41      continue
    end do
    if (jmax .le. 0) then
        ii = ii - 1
        go to 42
    end if
    ioff(ii) = i - n
    ind(i) = -jmax
    if (ii .lt. idiag) then
        go to 4
    end if
42  idiag = ii
50  continue
    do j = 1, idiag
        do i = 1, n
            diag(i, j) = 0.0d0
54          continue
        end do
55      continue
    end do
    ko = 1
    do i = 1, n
        do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            do l = 1, idiag
                if (j - i .ne. ioff(l)) then
                    go to 52
                end if
                diag(i, l) = a(k)
                go to 51
52              continue
            end do
            if (job2 .eq. 0) then
                go to 51
            end if
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko + 1
51          continue
        end do
        if (job2 .ne. 0) then
            ind(i + 1) = ko
        end if
6       continue
    end do
    if (job2 .eq. 0) then
        return
    end if
    iao(1) = 1
    do i = 2, n + 1
        iao(i) = ind(i)
7       continue
    end do
    return
end subroutine csrdia

subroutine csrbnd(n, a, ja, ia, job, abd, nabd, lowd, ml, mu, ierr)
    real(8) :: a(*), abd(nabd, n)
    integer :: ia(n + 1), ja(*)
    ierr = 0
    if (job .eq. 1) then
        call getbwd(n, a, ja, ia, ml, mu)
    end if
    m = ml + mu + 1
    if (lowd .eq. 0) then
        lowd = m
    end if
    if (m .gt. lowd) then
        ierr = -2
    end if
    if (lowd .gt. nabd .or. lowd .lt. 0) then
        ierr = -1
    end if
    if (ierr .lt. 0) then
        return
    end if
    do i = 1, m
        ii = lowd - i + 1
        do j = 1, n
            abd(ii, j) = 0.0d0
10          continue
        end do
15      continue
    end do
    mdiag = lowd - ml
    do i = 1, n
        do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            abd(i - j + mdiag, j) = a(k)
20          continue
        end do
30      continue
    end do
    return
end subroutine csrbnd

subroutine rperm(nrow, a, ja, ia, ao, jao, iao, perm, job)
    integer :: nrow, ja(*), ia(nrow + 1), jao(*), iao(nrow + 1), perm(nrow), job
    real(8) :: a(*), ao(*)
    logical :: values
    values = job .eq. 1
    do j = 1, nrow
        i = perm(j)
        iao(i + 1) = ia(j + 1) - ia(j)
50      continue
    end do
    iao(1) = 1
    do j = 1, nrow
        iao(j + 1) = iao(j + 1) + iao(j)
51      continue
    end do
    do ii = 1, nrow
        ko = iao(perm(ii))
        do k = ia(ii), ia(ii + 1) - 1
            jao(ko) = ja(k)
            if (values) then
                ao(ko) = a(k)
            end if
            ko = ko + 1
60          continue
        end do
100     continue
    end do
    return
end subroutine rperm

subroutine cperm(nrow, a, ja, ia, ao, jao, iao, perm, job)
    integer :: nrow, ja(*), ia(nrow + 1), jao(*), iao(nrow + 1), perm(*), job
    real(8) :: a(*), ao(*)
    integer :: k, i, nnz
    nnz = ia(nrow + 1) - 1
    do k = 1, nnz
        jao(k) = perm(ja(k))
100     continue
    end do
    if (job .ne. 1) then
        return
    end if
    do i = 1, nrow + 1
        iao(i) = ia(i)
1       continue
    end do
    do k = 1, nnz
        ao(k) = a(k)
2       continue
    end do
    return
end subroutine cperm

subroutine dperm(nrow, a, ja, ia, ao, jao, iao, perm, qperm, job)
    integer :: nrow, ja(*), ia(nrow + 1), jao(*), iao(nrow + 1), perm(nrow), qperm(*), job
    real(8) :: a(*), ao(*)
    integer :: locjob, mod
    locjob = mod(job, 2)
    call rperm(nrow, a, ja, ia, ao, jao, iao, perm, locjob)
    locjob = 0
    if (job .le. 2) then
        call cperm(nrow, ao, jao, iao, ao, jao, iao, perm, locjob)
    else
        call cperm(nrow, ao, jao, iao, ao, jao, iao, qperm, locjob)
    end if
    return
end subroutine dperm

subroutine dvperm(n, x, perm)
    integer :: n, perm(n)
    real(8) :: x(n)
    real(8) :: tmp, tmp1
    init = 1
    tmp = x(init)
    ii = perm(init)
    perm(init) = -perm(init)
    k = 0
6   k = k + 1
    tmp1 = x(ii)
    x(ii) = tmp
    next = perm(ii)
    if (next .lt. 0) then
        go to 65
    end if
    if (k .gt. n) then
        go to 101
    end if
    tmp = tmp1
    perm(ii) = -perm(ii)
    ii = next
    go to 6
65  init = init + 1
    if (init .gt. n) then
        go to 101
    end if
    if (perm(init) .lt. 0) then
        go to 65
    end if
    tmp = x(init)
    ii = perm(init)
    perm(init) = -perm(init)
    go to 6
101 continue
    do j = 1, n
        perm(j) = -perm(j)
200     continue
    end do
    return
end subroutine dvperm

subroutine ivperm(n, ix, perm)
    integer :: n, perm(n), ix(n)
    integer :: tmp, tmp1
    init = 1
    tmp = ix(init)
    ii = perm(init)
    perm(init) = -perm(init)
    k = 0
6   k = k + 1
    tmp1 = ix(ii)
    ix(ii) = tmp
    next = perm(ii)
    if (next .lt. 0) then
        go to 65
    end if
    if (k .gt. n) then
        go to 101
    end if
    tmp = tmp1
    perm(ii) = -perm(ii)
    ii = next
    go to 6
65  init = init + 1
    if (init .gt. n) then
        go to 101
    end if
    if (perm(init) .lt. 0) then
        go to 65
    end if
    tmp = ix(init)
    ii = perm(init)
    perm(init) = -perm(init)
    go to 6
101 continue
    do j = 1, n
        perm(j) = -perm(j)
200     continue
    end do
    return
end subroutine ivperm

subroutine diapos(n, ja, ia, idiag)
    integer :: ia(n + 1), ja(*), idiag(n)
    do i = 1, n
        idiag(i) = 0
1       continue
    end do
    do i = 1, n
        do k = ia(i), ia(i + 1) - 1
            if (ja(k) .eq. i) then
                idiag(i) = k
            end if
51          continue
        end do
6       continue
    end do
    return
end subroutine diapos

subroutine getbwd(n, a, ja, ia, ml, mu)
    real(8) :: a(*)
    integer :: ja(*), ia(n + 1), ml, mu, ldist, i, k
    ml = -n
    mu = -n
    do i = 1, n
        do k = ia(i), ia(i + 1) - 1
            ldist = i - ja(k)
            ml = max(ml, ldist)
            mu = max(mu, -ldist)
31          continue
        end do
3       continue
    end do
    return
end subroutine getbwd

subroutine infdia(n, ja, ia, ind, idiag)
    integer :: ia(*), ind(*), ja(*)
    n2 = n + n - 1
    do i = 1, n2
        ind(i) = 0
1       continue
    end do
    do i = 1, n
        do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            ind(n + j - i) = ind(n + j - i) + 1
2           continue
        end do
3       continue
    end do
    idiag = 0
    do k = 1, n2
        if (ind(k) .ne. 0) then
            idiag = idiag + 1
        end if
41      continue
    end do
    return
end subroutine infdia

subroutine rnrms(nrow, nrm, a, ja, ia, diag)
    real(8) :: a(*), diag(nrow), scal
    integer :: ja(*), ia(nrow + 1)
    do ii = 1, nrow
        scal = 0.0d0
        k1 = ia(ii)
        k2 = ia(ii + 1) - 1
        if (nrm .eq. 0) then
            do k = k1, k2
                scal = max(scal, abs(a(k)))
2               continue
            end do
        else
            if (nrm .eq. 1) then
                do k = k1, k2
                    scal = scal + abs(a(k))
3                   continue
                end do
            else
                do k = k1, k2
                    scal = scal + a(k)**2
4                   continue
                end do
            end if
        end if
        if (nrm .eq. 2) then
            scal = sqrt(scal)
        end if
        diag(ii) = scal
1       continue
    end do
    return
end subroutine rnrms

subroutine roscal(nrow, job, nrm, a, ja, ia, diag, b, jb, ib, ierr)
    real(8) :: a(*), b(*), diag(nrow)
    integer :: nrow, job, nrm, ja(*), jb(*), ia(nrow + 1), ib(nrow + 1), ierr
    call rnrms(nrow, nrm, a, ja, ia, diag)
    ierr = 0
    do j = 1, nrow
        if (diag(j) .eq. 0.0d0) then
            ierr = j
            return
        else
            diag(j) = 1.0d0/diag(j)
        end if
1       continue
    end do
    call diamua(nrow, job, a, ja, ia, diag, b, jb, ib)
    return
end subroutine roscal

subroutine ilut(n, a, ja, ia, lfil, droptol, alu, jlu, ju, iwk, w, jw, ierr)
    integer :: n
    real(8) :: a(*), alu(*), w(n), droptol
    integer :: ja(*), ia(n + 1), jlu(*), ju(n), jw(2*n), lfil, iwk, ierr
    integer :: ju0, k, j1, j2, j, ii, i, lenl, lenu, jj, jrow, jpos, len
    real(8) :: tnorm, t, abs, s, fact
    if (lfil .lt. 0) then
        go to 998
    end if
    ju0 = n + 2
    jlu(1) = ju0
    do j = 1, n
        jw(n + j) = 0
1       continue
    end do
    do ii = 1, n
        j1 = ia(ii)
        j2 = ia(ii + 1) - 1
        tnorm = 0.0d0
        do k = j1, j2
            tnorm = tnorm + abs(a(k))
501         continue
        end do
        if (tnorm .eq. 0.0d0) then
            go to 999
        end if
        tnorm = tnorm/real(j2 - j1 + 1)
        lenu = 1
        lenl = 0
        jw(ii) = ii
        w(ii) = 0.0
        jw(n + ii) = ii
        do j = j1, j2
            k = ja(j)
            t = a(j)
            if (k .lt. ii) then
                lenl = lenl + 1
                jw(lenl) = k
                w(lenl) = t
                jw(n + k) = lenl
            else
                if (k .eq. ii) then
                    w(ii) = t
                else
                    lenu = lenu + 1
                    jpos = ii + lenu - 1
                    jw(jpos) = k
                    w(jpos) = t
                    jw(n + k) = jpos
                end if
            end if
170         continue
        end do
        jj = 0
        len = 0
150     jj = jj + 1
        if (jj .gt. lenl) then
            go to 160
        end if
        jrow = jw(jj)
        k = jj
        do j = jj + 1, lenl
            if (jw(j) .lt. jrow) then
                jrow = jw(j)
                k = j
            end if
151         continue
        end do
        if (k .ne. jj) then
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
            jw(n + jrow) = jj
            jw(n + j) = k
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
        end if
        jw(n + jrow) = 0
        fact = w(jj)*alu(jrow)
        if (abs(fact) .le. droptol) then
            go to 150
        end if
        do k = ju(jrow), jlu(jrow + 1) - 1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n + j)
            if (j .ge. ii) then
                if (jpos .eq. 0) then
                    lenu = lenu + 1
                    if (lenu .gt. n) then
                        go to 995
                    end if
                    i = ii + lenu - 1
                    jw(i) = j
                    jw(n + j) = i
                    w(i) = -s
                else
                    w(jpos) = w(jpos) - s
                end if
            else
                if (jpos .eq. 0) then
                    lenl = lenl + 1
                    if (lenl .gt. n) then
                        go to 995
                    end if
                    jw(lenl) = j
                    jw(n + j) = lenl
                    w(lenl) = -s
                else
                    w(jpos) = w(jpos) - s
                end if
            end if
203         continue
        end do
        len = len + 1
        w(len) = fact
        jw(len) = jrow
        go to 150
160     continue
        do k = 1, lenu
            jw(n + jw(ii + k - 1)) = 0
308         continue
        end do
        lenl = len
        len = min0(lenl, lfil)
        call qsplit(w, jw, lenl, len)
        do k = 1, len
            if (ju0 .gt. iwk) then
                go to 996
            end if
            alu(ju0) = w(k)
            jlu(ju0) = jw(k)
            ju0 = ju0 + 1
204         continue
        end do
        ju(ii) = ju0
        len = 0
        do k = 1, lenu - 1
            if (abs(w(ii + k)) .gt. droptol*tnorm) then
                len = len + 1
                w(ii + len) = w(ii + k)
                jw(ii + len) = jw(ii + k)
            end if
        end do
        lenu = len + 1
        len = min0(lenu, lfil)
        if (lenu .gt. 1) then
            call qsplit(w(ii + 1), jw(ii + 1), lenu - 1, len)
        end if
        t = abs(w(ii))
        if (len + ju0 .gt. iwk) then
            go to 997
        end if
        do k = ii + 1, ii + len - 1
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k))
            ju0 = ju0 + 1
302         continue
        end do
        if (w(ii) .eq. 0.0d0) then
            w(ii) = (0.0001 + droptol)*tnorm
        end if
        alu(ii) = 1.0d0/w(ii)
        jlu(ii + 1) = ju0
500     continue
    end do
    ierr = 0
    return
995 ierr = -1
    return
996 ierr = -2
    return
997 ierr = -3
    return
998 ierr = -4
    return
999 ierr = -5
    return
end subroutine ilut

subroutine ilutp(n, a, ja, ia, lfil, droptol, permtol, mbloc, alu, jlu, ju, iwk, w, jw, iperm, ierr)
    integer :: n, ja(*), ia(n + 1), lfil, jlu(*), ju(n), jw(2*n), iwk, iperm(2*n), ierr
    real(8) :: a(*), alu(*), w(n), droptol
    integer :: k, i, j, jrow, ju0, ii, j1, j2, jpos, len, imax, lenu, lenl, jj, mbloc, icut
    real(8) :: s, tmp, tnorm, xmax, xmax0, fact, abs, t, permtol
    if (lfil .lt. 0) then
        go to 998
    end if
    ju0 = n + 2
    jlu(1) = ju0
    do j = 1, n
        jw(n + j) = 0
        iperm(j) = j
        iperm(n + j) = j
1       continue
    end do
    do ii = 1, n
        j1 = ia(ii)
        j2 = ia(ii + 1) - 1
        tnorm = 0.0d0
        do k = j1, j2
            tnorm = tnorm + abs(a(k))
501         continue
        end do
        if (tnorm .eq. 0.0d0) then
            go to 999
        end if
        tnorm = tnorm/(j2 - j1 + 1)
        lenu = 1
        lenl = 0
        jw(ii) = ii
        w(ii) = 0.0
        jw(n + ii) = ii
        do j = j1, j2
            k = iperm(n + ja(j))
            t = a(j)
            if (k .lt. ii) then
                lenl = lenl + 1
                jw(lenl) = k
                w(lenl) = t
                jw(n + k) = lenl
            else
                if (k .eq. ii) then
                    w(ii) = t
                else
                    lenu = lenu + 1
                    jpos = ii + lenu - 1
                    jw(jpos) = k
                    w(jpos) = t
                    jw(n + k) = jpos
                end if
            end if
170         continue
        end do
        jj = 0
        len = 0
150     jj = jj + 1
        if (jj .gt. lenl) then
            go to 160
        end if
        jrow = jw(jj)
        k = jj
        do j = jj + 1, lenl
            if (jw(j) .lt. jrow) then
                jrow = jw(j)
                k = j
            end if
151         continue
        end do
        if (k .ne. jj) then
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
            jw(n + jrow) = jj
            jw(n + j) = k
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
        end if
        jw(n + jrow) = 0
        fact = w(jj)*alu(jrow)
        if (abs(fact) .le. droptol) then
            go to 150
        end if
        do k = ju(jrow), jlu(jrow + 1) - 1
            s = fact*alu(k)
            j = iperm(n + jlu(k))
            jpos = jw(n + j)
            if (j .ge. ii) then
                if (jpos .eq. 0) then
                    lenu = lenu + 1
                    i = ii + lenu - 1
                    if (lenu .gt. n) then
                        go to 995
                    end if
                    jw(i) = j
                    jw(n + j) = i
                    w(i) = -s
                else
                    w(jpos) = w(jpos) - s
                end if
            else
                if (jpos .eq. 0) then
                    lenl = lenl + 1
                    if (lenl .gt. n) then
                        go to 995
                    end if
                    jw(lenl) = j
                    jw(n + j) = lenl
                    w(lenl) = -s
                else
                    w(jpos) = w(jpos) - s
                end if
            end if
203         continue
        end do
        len = len + 1
        w(len) = fact
        jw(len) = jrow
        go to 150
160     continue
        do k = 1, lenu
            jw(n + jw(ii + k - 1)) = 0
308         continue
        end do
        lenl = len
        len = min0(lenl, lfil)
        call qsplit(w, jw, lenl, len)
        do k = 1, len
            if (ju0 .gt. iwk) then
                go to 996
            end if
            alu(ju0) = w(k)
            jlu(ju0) = iperm(jw(k))
            ju0 = ju0 + 1
204         continue
        end do
        ju(ii) = ju0
        len = 0
        do k = 1, lenu - 1
            if (abs(w(ii + k)) .gt. droptol*tnorm) then
                len = len + 1
                w(ii + len) = w(ii + k)
                jw(ii + len) = jw(ii + k)
            end if
        end do
        lenu = len + 1
        len = min0(lenu, lfil)
        if (lenu .gt. 1) then
            call qsplit(w(ii + 1), jw(ii + 1), lenu - 1, len)
        end if
        imax = ii
        xmax = abs(w(imax))
        xmax0 = xmax
        icut = ii - 1 + mbloc - mod(ii - 1, mbloc)
        do k = ii + 1, ii + len - 1
            t = abs(w(k))
            if (t .gt. xmax .and. t*permtol .gt. xmax0 .and. jw(k) .le. icut) then
                imax = k
                xmax = t
            end if
        end do
        tmp = w(ii)
        w(ii) = w(imax)
        w(imax) = tmp
        j = jw(imax)
        i = iperm(ii)
        iperm(ii) = iperm(j)
        iperm(j) = i
        iperm(n + iperm(ii)) = ii
        iperm(n + iperm(j)) = j
        if (len + ju0 .gt. iwk) then
            go to 997
        end if
        do k = ii + 1, ii + len - 1
            jlu(ju0) = iperm(jw(k))
            alu(ju0) = w(k)
            ju0 = ju0 + 1
302         continue
        end do
        if (w(ii) .eq. 0.0d0) then
            w(ii) = (1.0d-4 + droptol)*tnorm
        end if
        alu(ii) = 1.0d0/w(ii)
        jlu(ii + 1) = ju0
500     continue
    end do
    do k = jlu(1), jlu(n + 1) - 1
        jlu(k) = iperm(n + jlu(k))
    end do
    do k = ia(1), ia(n + 1) - 1
        ja(k) = iperm(n + ja(k))
    end do
    ierr = 0
    return
995 ierr = -1
    return
996 ierr = -2
    return
997 ierr = -3
    return
998 ierr = -4
    return
999 ierr = -5
    return
end subroutine ilutp

subroutine qsplit(a, ind, n, ncut)
    real(8) :: a(n)
    integer :: ind(n), n, ncut
    real(8) :: tmp, abskey
    integer :: itmp, first, last
    first = 1
    last = n
    if (ncut .lt. first .or. ncut .gt. last) then
        return
    end if
1   mid = first
    abskey = abs(a(mid))
    do j = first + 1, last
        if (abs(a(j)) .gt. abskey) then
            mid = mid + 1
            tmp = a(mid)
            itmp = ind(mid)
            a(mid) = a(j)
            ind(mid) = ind(j)
            a(j) = tmp
            ind(j) = itmp
        end if
2       continue
    end do
    tmp = a(mid)
    a(mid) = a(first)
    a(first) = tmp
    itmp = ind(mid)
    ind(mid) = ind(first)
    ind(first) = itmp
    if (mid .eq. ncut) then
        return
    end if
    if (mid .gt. ncut) then
        last = mid - 1
    else
        first = mid + 1
    end if
    go to 1
end subroutine qsplit

subroutine lusol(n, y, x, alu, jlu, ju)
    real(8) :: x(n), y(n), alu(*)
    integer :: n, jlu(*), ju(*)
    integer :: i, k
    do i = 1, n
        x(i) = y(i)
        do k = jlu(i), ju(i) - 1
            x(i) = x(i) - alu(k)*x(jlu(k))
41          continue
        end do
40      continue
    end do
    do i = n, 1, -1
        do k = ju(i), jlu(i + 1) - 1
            x(i) = x(i) - alu(k)*x(jlu(k))
91          continue
        end do
        x(i) = alu(i)*x(i)
90      continue
    end do
    return
end subroutine lusol

subroutine dblstr(n, ja, ia, ip1, ip2, nfirst, riord, ndom, map, mapptr, mask, levels, iwk)
    integer :: ndom, ja(*), ia(*), ip1, ip2, nfirst, riord(*), map(*), mapptr(*), mask(*), levels(*), iwk(*), nextdom
    integer :: n, j, idom, kdom, jdom, maskval, k, nlev, init, ndp1, numnod
    maskval = 1
    do j = 1, n
        mask(j) = maskval
    end do
    iwk(1) = 0
    call bfs(n, ja, ia, nfirst, iwk, mask, maskval, riord, levels, nlev)
    call stripes(nlev, riord, levels, ip1, map, mapptr, ndom)
    if (ip2 .eq. 1) then
        return
    end if
    ndp1 = ndom + 1
    do j = 1, ndom + 1
        iwk(j) = ndp1 + mapptr(j)
    end do
    do j = 1, mapptr(ndom + 1) - 1
        iwk(ndp1 + j) = map(j)
    end do
    do idom = 1, ndom
        j = iwk(idom)
        numnod = iwk(idom + 1) - iwk(idom)
        init = iwk(j)
        do k = j, iwk(idom + 1) - 1
        end do
    end do
    do idom = 1, ndom
        do k = mapptr(idom), mapptr(idom + 1) - 1
            mask(map(k)) = idom
        end do
    end do
    nextdom = 1
    jdom = 1
    mapptr(jdom) = 1
    do idom = 1, ndom
        maskval = idom
        nfirst = 1
        numnod = iwk(idom + 1) - iwk(idom)
        j = iwk(idom)
        init = iwk(j)
        nextdom = mapptr(jdom)
        call perphn(numnod, ja, ia, init, mask, maskval, nlev, riord, levels)
        call stripes(nlev, riord, levels, ip2, map(nextdom), mapptr(jdom), kdom)
        mapptr(jdom) = nextdom
        do j = jdom, jdom + kdom - 1
            mapptr(j + 1) = nextdom + mapptr(j + 1) - 1
        end do
        jdom = jdom + kdom
    end do
    ndom = jdom - 1
    return
end subroutine dblstr

subroutine bfs(n, ja, ia, nfirst, iperm, mask, maskval, riord, levels, nlev)
    integer :: n, ja(*), ia(*), nfirst, iperm(n), mask(n), riord(*), levels(*), nlev, maskval
    integer :: j, ii, nod, istart, iend
    logical :: permut
    permut = iperm(1) .ne. 0
    nlev = 0
    istart = 0
    ii = 0
    iend = nfirst
    do j = 1, nfirst
        mask(riord(j)) = 0
12      continue
    end do
13  continue
1   nlev = nlev + 1
    levels(nlev) = istart + 1
    call add_lvst(istart, iend, nlev, riord, ja, ia, mask, maskval)
    if (istart .lt. iend) then
        go to 1
    end if
2   ii = ii + 1
    if (ii .le. n) then
        nod = ii
        if (permut) then
            nod = iperm(nod)
        end if
        if (mask(nod) .eq. maskval) then
            istart = iend
            iend = iend + 1
            riord(iend) = nod
            mask(nod) = 0
            go to 1
        else
            go to 2
        end if
    end if
3   levels(nlev + 1) = iend + 1
    do j = 1, iend
        mask(riord(j)) = maskval
    end do
    return
end subroutine bfs

subroutine add_lvst(istart, iend, nlev, riord, ja, ia, mask, maskval)
    integer :: nlev, nod, riord(*), ja(*), ia(*), mask(*)
    nod = iend
    do ir = istart + 1, iend
        i = riord(ir)
        do k = ia(i), ia(i + 1) - 1
            j = ja(k)
            if (mask(j) .eq. maskval) then
                nod = nod + 1
                mask(j) = 0
                riord(nod) = j
            end if
24          continue
        end do
25      continue
    end do
    istart = iend
    iend = nod
    return
end subroutine add_lvst

subroutine stripes(nlev, riord, levels, ip, map, mapptr, ndom)
    integer :: nlev, riord(*), levels(nlev + 1), ip, map(*), mapptr(*), ndom
    integer :: ib, ktr, ilev, k, nsiz, psiz
    ndom = 1
    ib = 1
    nsiz = levels(nlev + 1) - levels(1)
    psiz = (nsiz - ib)/max(1, ip - ndom + 1) + 1
    mapptr(ndom) = ib
    ktr = 0
    do ilev = 1, nlev
        do k = levels(ilev), levels(ilev + 1) - 1
            map(ib) = riord(k)
            ib = ib + 1
            ktr = ktr + 1
            if (ktr .ge. psiz .or. k .ge. nsiz) then
                ndom = ndom + 1
                mapptr(ndom) = ib
                psiz = (nsiz - ib)/max(1, ip - ndom + 1) + 1
                ktr = 0
            end if
3           continue
        end do
10      continue
    end do
    ndom = ndom - 1
    return
end subroutine stripes

integer function maskdeg(ja, ia, nod, mask, maskval)
    integer :: ja(*), ia(*), nod, mask(*), maskval
    integer :: deg, k
    deg = 0
    do k = ia(nod), ia(nod + 1) - 1
        if (mask(ja(k)) .eq. maskval) then
            deg = deg + 1
        end if
    end do
    maskdeg = deg
    return
end function maskdeg

subroutine perphn(n, ja, ia, init, mask, maskval, nlev, riord, levels)
    integer :: n, ja(*), ia(*), init, mask(*), maskval, nlev, riord(*), levels(*)
    integer :: j, nlevp, deg, nfirst, mindeg, nod, maskdeg
    integer :: iperm(1)
    nlevp = 0
1   continue
    riord(1) = init
    nfirst = 1
    iperm(1) = 0
    call bfs(n, ja, ia, nfirst, iperm, mask, maskval, riord, levels, nlev)
    if (nlev .gt. nlevp) then
        mindeg = n + 1
        do j = levels(nlev), levels(nlev + 1) - 1
            nod = riord(j)
            deg = maskdeg(ja, ia, nod, mask, maskval)
            if (deg .lt. mindeg) then
                init = nod
                mindeg = deg
            end if
        end do
        nlevp = nlev
        go to 1
    end if
    return
end subroutine perphn

subroutine rversp(n, riord)
    integer :: n, riord(n)
    integer :: j, k
    do j = 1, n/2
        k = riord(j)
        riord(j) = riord(n - j + 1)
        riord(n - j + 1) = k
26      continue
    end do
    return
end subroutine rversp

subroutine atob(n, a, ja, ia, b, jb, ib)
    integer :: n
    real(8) :: a(*)
    integer :: ja(*)
    integer :: ia(n + 1)
    real(8) :: b(*)
    integer :: jb(*)
    integer :: ib(n + 1)
    integer :: i
    do i = 1, ia(n + 1) - 1
        b(i) = a(i)
        jb(i) = ja(i)
    end do
    do i = 1, n + 1
        ib(i) = ia(i)
    end do
    return
end subroutine atob

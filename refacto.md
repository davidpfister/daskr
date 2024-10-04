# Refactoring DASKR

## Introduction

I am working on the refactoring of daskr using the historical code on netlib as a starting point. 
It's a legacy fixed-form fortran 77 piece of code with it's implicit typing, implicit procedures, 
goto's, etc. 
It's a great piece of work and I have in mind to turn it into a modern code. 

## Automatic refactoring

### Using lfortran

LFortran has the amazing capability of parsing pretty much any fortran file and output an ast. Reversely, it can also read the ast and output a fortran file, the latter being written in modern fortran this time.

One can use the following command 
```
lfortran dlinpk.f --fixed-form --implicit-typing --show-ast-f90 --no-color > dlinpk.f90
```

Or this script to deal with batches
```powershell
Get-ChildItem "./src" -Filter *.f | Foreach-Object {
    & lfortran $_.FullName --fixed-form --implicit-typing --show-ast-f90 --no-color > $($_.FullName + "90")
}
```

The output is the following

Addition of '::' separator
```
!old
integer lda,n,ipvt(1),info

!new
integer :: lda, n, ipvt(1), info
```

Removal of all comments

Using c-style mathematical expressions

```
!old
if (nm1 .lt. 1) go to 70

!new
if (nm1 < 1) then
```

Expanding if-block
```
!old
if (nm1 .lt. 1) go to 70

!new
if (nm1 < 1) then
    go to 70
end if
```

Removing numbered labelled do-loop
```
!old
do 30 j = kp1, n
            t = a(l,j)
            if (l .eq. k) go to 20
                
            
                a(l,j) = a(k,j)
                a(k,j) = t
20          continue
            call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30       continue

!new
do j = kp1, n
        t = a(l, j)
        if (l == k) then
            go to 20
        end if
        a(l, j) = a(k, j)
        a(k, j) = t
        20 continue
        call daxpy(n - k, t, a(k + 1, k), 1, a(k + 1, j), 1)
        30 continue
    end do
end if
```

Add missing 'end subroutine'
```
return
end

!new
return
end subroutine dgesl
```

Add missing `program` and implicit none statement
```
program implicit_program_lfortran
implicit double precision (a-h,o-z)
```

## f77_to_f90
For converting sparskit

### Using fprettify

```
fprettify .\src\ -r --case 1 1 1 1 -i 4 --strict-indent --enable-replacements --strip-comments
```

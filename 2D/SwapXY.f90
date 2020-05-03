subroutine SwapXY()

implicit none

! Variables (general)
real, pointer :: a(:, :), rhs(:, :)

real, pointer :: x(:), y(:)

real dx, dy
integer lx, ly
real llx, lly							! These are Lx and Ly, respectively, in MatLab.
integer axisPlusX, axisPlusY

real, pointer :: myu(:, :)

real pi 

common /arhs/ a, rhs
common /axes/ x, y
common /steps/ dx, dy
common /limits/ llx, lly
common /sizes/ lx, ly
common /axisPluses/ axisPlusX, axisPlusY
common /myu/ myu
common /consts/ pi

! Variables (partial)
real, pointer :: xNew(:)
integer lxNew
real llxNew

real, pointer :: yNew(:)
integer lyNew
real llyNew							! This is LyNew in MatLab

integer axisPlusXNew, axisPlusYNew

real, pointer :: myuNew(:, :)
real, pointer :: aNew(:, :)

integer i

! Body
axisPlusXNew = 0
axisPlusYNew = 1

llxNew = pi
llyNew = 2.0 * pi

allocate(xNew(int((llxNew - dx) / dx + 1)));					lxNew = size(xNew, 1)
xNew = [0 : lxNew - 1] * dx + dx / 2.0

allocate(yNew(int((llyNew + dy) / dy + 1)));					lyNew = size(yNew, 1)
yNew = -lly + dy / 2.0 + [0 : lyNew - 1] * dy

allocate(myuNew(lxNew, lyNew))
allocate(aNew(lxNew, lyNew))

do i = 1, lxNew
   myuNew(i, :) = [myu(i, :), myu(size(myu, 1) - axisPlusX - i, size(myu, 2) : 1 : -1), myu(i, 1 : 2)]
    
   aNew(i, :) = [a(i, :), a(size(a, 1) - axisPlusX - i, size(a, 2) : 1 : -1), a(i, 1 : 2)]
end do

llx = llxNew
lly = llyNew

deallocate(x)
allocate(x(lxNew))
x = xNew
lx = lxNew

deallocate(y)
allocate(y(lyNew))
y = yNew
ly = lyNew

axisPlusX = axisPlusXNew
axisPlusY = axisPlusYNew

deallocate(myu)
allocate(myu(lxNew, lyNew))
myu = myuNew
!do i = 1, lxNew
!   un(i, :) = unNew(i, :)
!   vn(i, :) = vnNew(i, :)
!   zn(i, :) = znNew(i, :)
!end do
! The assignement "un = unNew" (and the similars) does not work when "lxNew" and "lyNew" are large (say, for lxNew = 360 and lyNew = 180)
! So, we used a cycle

deallocate(a)
allocate(a(lxNew, lyNew))
a = aNew
!do i = 1, lxNew
!   U(i, :) = UNew(i, :)
!   V(i, :) = VNew(i, :)
!   a(i, :) = aNew(i, :)
!   zb(i, :) = zbNew(i, :)
!end do

deallocate(xNew)
deallocate(yNew)
deallocate(myuNew)
deallocate(aNew)

end subroutine SwapXY

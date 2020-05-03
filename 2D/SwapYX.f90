subroutine SwapYX()

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

integer j

! Body
axisPlusXNew = 1
axisPlusYNew = 0

llxNew = 2.0 * pi
llyNew = pi / 2.0
!llxNew = 1.0
!llyNew = 0.25

allocate(xNew(int((llxNew + dx) / dx + 1)));					lxNew = size(xNew, 1)
xNew = [0 : lxNew - 1] * dx + dx / 2.0
	
allocate(yNew(int((2.0 * llyNew - dy) / dy + 1)));				lyNew = size(yNew, 1)
yNew = -llyNew + dy / 2.0 + [0 : lyNew - 1] * dy

allocate(myuNew(lxNew, lyNew))
allocate(aNew(lxNew, lyNew))

do j = 1, (lxNew - axisPlusXNew * 2) / 2
   myuNew(j, :) = myu(j, 1 : (size(myu, 2) - axisPlusXNew * 2) / 2)
   myuNew(size(myuNew, 1) - axisPlusXNew - j, :) = myu(j, size(myu, 2) - axisPlusXNew * 2 : (size(myu, 2) - axisPlusXNew * 2) / 2 + 1 : -1)
   myuNew(size(myuNew, 1) - axisPlusXNew : size(myuNew, 1), :) = myuNew(1 : 2, :)

   aNew(j, :) = a(j, 1 : (size(a, 2) - axisPlusXNew * 2) / 2)
   aNew(size(aNew, 1) - axisPlusXNew - j, :) = a(j, size(a, 2) - axisPlusXNew * 2 : (size(a, 2) - axisPlusXNew * 2) / 2 + 1 : -1)
   aNew(size(aNew, 1) - axisPlusXNew : size(aNew, 1), :) = aNew(1 : 2, :)
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

deallocate(a)
allocate(a(lxNew, lyNew))
a = aNew

deallocate(xNew)
deallocate(yNew)
deallocate(myuNew)
deallocate(aNew)

end subroutine SwapYX

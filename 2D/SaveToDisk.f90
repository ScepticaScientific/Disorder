subroutine SaveToDisk(filenameNo)

implicit none

! Variables (general)
real, pointer :: a(:, :), rhs(:, :)

common /arhs/ a, rhs

! Variables (partial)
integer filenameNo
character*5 the_file

integer i, j

! Body
open (1, file = 'dummy.tx2')
write (1, *) filenameNo
rewind(1)
read (1, *) the_file
close (1)

! Velocity field
! Depth
open (1, file = 'a_' // trim(the_file) // '.tx2')
do i = 1, size(a, 1)
   !write (1, '(D20.10$)') (a(i, j), j = 1, size(a, 2) - 1)
   !write (1, '(D20.10)') a(i, size(a, 2))
   write (1, '(E20.10E3$)') (a(i, j), j = 1, size(a, 2) - 1)
   write (1, '(E20.10E3)') a(i, size(a, 2))
end do 
close (1)

end subroutine SaveToDisk
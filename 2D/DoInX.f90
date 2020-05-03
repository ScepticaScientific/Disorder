subroutine DoInX(tau, res)

!DEC$ DEFINE APPROXIMATION_TYPE_1
! APPROXIMATION_TYPE_1 - Crank-Nicolson
! APPROXIMATION_TYPE_2 - Euler (purely explicit)

!DEC$ IF DEFINED (VF)
use linear_operators
!DEC$ ENDIF

use isDiagPredomFunction

implicit none

!DEC$ IF .NOT. DEFINED (VF)
include 'mpif.h'
!DEC$ ENDIF

! Variables (general)
integer lx, ly
real dx, dy

integer axisPlusX, axisPlusY

logical IsPlaneModel

real, pointer :: a(:, :), rhs(:, :)

real, pointer :: mycosTabX(:), mysinTabX(:), mytanTabX(:)

real, pointer :: myu(:, :)

real R_Earth

integer ApprOrder
logical IsFirstOrderMinusPlus, IsCheck

real, pointer :: Ax(:, :), bxa(:), ipivx(:)

integer, pointer :: indexMainDiagX(:), indexSuperDiagX(:), indexSuperDiagX2(:), indexSuperDiagX3(:), indexSuperDiagX4(:), indexSubDiagX(:), indexSubDiagX2(:), indexSubDiagX3(:), indexSubDiagX4(:)
real, pointer :: workingArrayX(:)

real cflX, cflY

!DEC$ IF .NOT. DEFINED (VF)
integer ierr, processID, NumberOfProcesses
!DEC$ ENDIF

common /arhs/ a, rhs
common /trigsX/ mycosTabX, mysinTabX, mytanTabX
common /steps/ dx, dy
common /sizes/ lx, ly
common /axisPluses/ axisPlusX, axisPlusY
common /myu/ myu
common /Radius/ R_Earth
common /specials/ ApprOrder, IsFirstOrderMinusPlus, IsCheck
common /matvecX/ Ax, bxa, ipivx
common /diagsX/ indexMainDiagX, indexSuperDiagX, indexSuperDiagX2, indexSuperDiagX3, indexSuperDiagX4, indexSubDiagX, indexSubDiagX2, indexSubDiagX3, indexSubDiagX4, workingArrayX
common /CFL/ cflX, cflY
common /theGeometry/ IsPlaneModel

!DEC$ IF .NOT. DEFINED (VF)
common /mpiData/ ierr, processID, NumberOfProcesses
!DEC$ ENDIF

! Variables (partial)
real tau
integer res

real, pointer :: a_next(:, :)

integer j
integer jFrom, jTo

integer, pointer :: sgn(:)

integer RemApprOrder

!DEC$ IF .NOT. DEFINED (VF)
integer info
external dgetrf, dgetrs
!DEC$ ENDIF

! Body
allocate(a_next(lx, ly))
allocate(sgn(lx))

RemApprOrder = ApprOrder
!ApprOrder = 2

jFrom = 1
jTo = ly - axisPlusY

!DEC$ IF .NOT. DEFINED (VF)
jFrom = 1 + (ly - axisPlusY) * 1.0 / NumberOfProcesses * processID
jTo = 1 + (ly - axisPlusY) * 1.0 / NumberOfProcesses * (processID + 1) - 1

!write (*, *) processID, jFrom, jTo
!DEC$ ENDIF

do j = jFrom, jTo
   if (IsCheck .eq. 1) then
      ApprOrder = RemApprOrder
   end if

   Ax(:, :) = 0.0

   Ax(1, 1) = 1.0
   Ax(1, lx - 1) = -1.0
    
   Ax(lx, 2) = 1.0
   Ax(lx, lx) = -1.0
    
   !a_next(:, j) = a(:, j)
    
   ! Computation.
   workingArrayX = reshape(Ax, (/product(shape(Ax))/))
   workingArrayX(indexSuperDiagX4) = 0.0
   workingArrayX(indexSuperDiagX3) = 0.0
   workingArrayX(indexSuperDiagX2) = 0.0
   workingArrayX(indexSubDiagX2) = 0.0
   workingArrayX(indexSubDiagX3) = 0.0
   workingArrayX(indexSubDiagX4) = 0.0

   if (ApprOrder .eq. 2) then
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
      workingArrayX(indexSuperDiagX2) = -myu(3 : lx, j) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
	  workingArrayX(indexMainDiagX) = 1.0 / tau + (myu(3 : lx, j) + myu(1 : lx - 2, j)) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
      workingArrayX(indexSubDiagX2) = -myu(1 : lx - 2, j) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
	  workingArrayX(indexMainDiagX) = 1.0 / tau
!DEC$ ENDIF
      Ax = reshape(workingArrayX, (/lx, lx/))
   else if (ApprOrder .eq. 4) then
	  workingArrayX(indexSuperDiagX4) = -myu([4 : lx, 3], j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
      workingArrayX(indexSuperDiagX3) = (8.0 * myu([4 : lx, 3], j) + 8.0 * myu(3 : lx, j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
      workingArrayX(indexSuperDiagX2) = -64.0 * myu(3 : lx, j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
 	  workingArrayX(indexSuperDiagX) = -(8.0 * myu([4 : lx, 3], j) + 8.0 * myu(1 : lx - 2, j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
	  workingArrayX(indexMainDiagX) = 1.0 / tau + (myu([4 : lx, 3], j) + 64.0 * myu(3 : lx, j) + 64.0 * myu(1 : lx - 2, j) + myu([lx - 2, 1 : lx - 3], j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
	  workingArrayX(indexSubDiagX) = -(8.0 * myu(3 : lx, j) + 8.0 * myu([lx - 2, 1 : lx - 3], j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
      workingArrayX(indexSubDiagX2) = -64.0 * myu(1 : lx - 2, j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
      workingArrayX(indexSubDiagX3) = (8.0 * myu([lx - 2, 1 : lx - 3], j) + 8.0 * myu(1 : lx - 2, j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
	  workingArrayX(indexSubDiagX4) = -myu([lx - 2, 1 : lx - 3], j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
      Ax = reshape(workingArrayX, (/lx, lx/))
   end if
    
   if (isDiagPredom(Ax) .eq. 0) then
      write (*, *) 'The matrix Ax is not diagonally predominant.'
      !res = 0
      !return
   end if
    
   bxa(1) = -a(1, j) + a(lx - 1, j)
   bxa(lx) = -a(2, j) + a(lx, j)
    
   if (ApprOrder .eq. 2) then
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  bxa(2 : lx - 1) = a(2 : lx - 1, j) * (1.0 / tau - (myu(3 : lx, j) + myu(1 : lx - 2, j)) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)) + &
												a([4 : lx, 3], j) * myu(3 : lx, j) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a([lx - 2, 1 : lx - 3], j) * myu(1 : lx - 2, j) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
	  bxa(2 : lx - 1) = a(2 : lx - 1, j) * (1.0 / tau - (myu(3 : lx, j) + myu(1 : lx - 2, j)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)) + &
												a([4 : lx, 3], j) * myu(3 : lx, j) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a([lx - 2, 1 : lx - 3], j) * myu(1 : lx - 2, j) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
!DEC$ ENDIF
   else if (ApprOrder .eq. 4) then
	  bxa(2 : lx - 1) = a(2 : lx - 1, j) * (1.0 / tau - (myu([4 : lx, 3], j) + 64.0 * myu(3 : lx, j) + 64.0 * myu(1 : lx - 2, j) + myu([lx - 2, 1 : lx - 3], j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)) + &
												a([6 : lx, 3, 4, 5], j) * myu([4 : lx, 3], j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) - &
												a([5 : lx, 3, 4], j) * (8.0 * myu([4 : lx, 3], j) + 8.0 * myu(3 : lx, j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a([4 : lx, 3], j) * 64.0 * myu(3 : lx, j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a(3 : lx, j) * (8.0 * myu([4 : lx, 3], j) + 8.0 * myu(1 : lx - 2, j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a(1 : lx - 2, j) * (8.0 * myu(3 : lx, j) + 8.0 * myu([lx - 2, 1 : lx - 3], j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a([lx - 2, 1 : lx - 3], j) * 64.0 * myu(1 : lx - 2, j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) - &
												a([lx - 3, lx - 2, 1 : lx - 4], j) * (8.0 * myu([lx - 2, 1 : lx - 3], j) + 8.0 * myu(1 : lx - 2, j)) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a([lx - 4, lx - 3, lx - 2, 1 : lx - 5], j) * myu([lx - 2, 1 : lx - 3], j) / (288.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
   end if

   cflX = max(maxval(abs(Ax)) * tau, cflX)

   !DEC$ IF DEFINED (VF)
   a_next(:, j) = Ax .ix. bxa
   !DEC$ ELSE
   ipivx(:) = 0
   call dgetrf(lx, lx, Ax, lx, ipivx, info)

   if (info .eq. 0) then
      call dgetrs('N', lx, 1, Ax, lx, ipivx, bxa, lx, info)

      a_next(:, j) = bxa
   else
      write (*, *) 'The matrix Ax is singular.'
   end if
   !DEC$ ENDIF
       
   ! Update
   a(:, j) = a_next(:, j)
end do

deallocate(sgn)
deallocate(a_next)

res = 1

ApprOrder = RemApprOrder

! Periodicity
if (IsPlaneModel .eq. 1) then
!DEC$ IF .NOT. DEFINED (VF)
   if (processID .eq. 0) then
!DEC$ ENDIF
	  a(1 : lx, ly - 1 : ly) = a(1 : lx, 1 : 2)
!DEC$ IF .NOT. DEFINED (VF)
   end if
!DEC$ ENDIF
end if

!DEC$ IF .NOT. DEFINED (VF)
call mpi_gather(a(1, jFrom), lx * (jTo - jFrom + 1), MPI_REAL8, a, lx * (jTo - jFrom + 1), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(a, product(shape(a)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

end subroutine DoInX
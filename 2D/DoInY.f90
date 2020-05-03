subroutine DoInY(tau, res)

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

real, pointer :: mycosTabY(:), mysinTabY(:), mytanTabY(:)

real, pointer :: myu(:, :)

real R_Earth

integer ApprOrder
logical IsFirstOrderMinusPlus, IsCheck

real, pointer :: Ay(:, :), bya(:), ipivy(:)

integer, pointer :: indexMainDiagY(:), indexSuperDiagY(:), indexSuperDiagY2(:), indexSuperDiagY3(:), indexSuperDiagY4(:), indexSubDiagY(:), indexSubDiagY2(:), indexSubDiagY3(:), indexSubDiagY4(:)
real, pointer :: workingArrayY(:)

real cflX, cflY

!DEC$ IF .NOT. DEFINED (VF)
integer ierr, processID, NumberOfProcesses
!DEC$ ENDIF

common /arhs/ a, rhs
common /trigsY/ mycosTabY, mysinTabY, mytanTabY
common /steps/ dx, dy
common /sizes/ lx, ly
common /axisPluses/ axisPlusX, axisPlusY
common /myu/ myu
common /Radius/ R_Earth
common /specials/ ApprOrder, IsFirstOrderMinusPlus, IsCheck
common /matvecY/ Ay, bya, ipivy
common /diagsY/ indexMainDiagY, indexSuperDiagY, indexSuperDiagY2, indexSuperDiagY3, indexSuperDiagY4, indexSubDiagY, indexSubDiagY2, indexSubDiagY3, indexSubDiagY4, workingArrayY
common /CFL/ cflX, cflY
common /theGeometry/ IsPlaneModel

!DEC$ IF .NOT. DEFINED (VF)
common /mpiData/ ierr, processID, NumberOfProcesses
!DEC$ ENDIF

! Variables (partial)
real tau
integer res

real, pointer :: a_next(:, :)

integer i
integer iFrom, iTo

integer, pointer :: sgn(:)

integer RemApprOrder

!DEC$ IF .NOT. DEFINED (VF)
integer info
external dgetrf, dgetrs
!DEC$ ENDIF

!DEC$ IF .NOT. DEFINED (VF)
integer j
!DEC$ ENDIF

! Body
allocate(a_next(lx, ly))
allocate(sgn(ly))

RemApprOrder = ApprOrder
!ApprOrder = 2

iFrom = 1
iTo = lx - 2 * axisPlusX

!DEC$ IF .NOT. DEFINED (VF)
iFrom = 1 + (lx - 2 * axisPlusX) * 1.0 / NumberOfProcesses * processID
iTo = 1 + (lx - 2 * axisPlusX) * 1.0 / NumberOfProcesses * (processID + 1) - 1

!write (*, *) processID, iFrom, iTo
!DEC$ ENDIF

do i = iFrom, iTo
   if (IsCheck .eq. 1) then
      ApprOrder = RemApprOrder
   end if

   Ay(:, :) = 0.0

   Ay(1, 1) = 1.0
   Ay(1, ly - 1) = -1.0
    
   Ay(ly, 2) = 1.0
   Ay(ly, ly) = -1.0
    
   !a_next(i, :) = a(i, :)
    
   ! Computation.
   workingArrayY = reshape(Ay, (/product(shape(Ay))/))
   workingArrayY(indexSuperDiagY4) = 0.0
   workingArrayY(indexSuperDiagY3) = 0.0
   workingArrayY(indexSuperDiagY2) = 0.0
   workingArrayY(indexSubDiagY2) = 0.0
   workingArrayY(indexSubDiagY3) = 0.0
   workingArrayY(indexSubDiagY4) = 0.0

   if (ApprOrder .eq. 2) then
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  workingArrayY(indexSuperDiagY2) = -myu(i, 3 : ly) * mycosTabY(3 : ly) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	  workingArrayY(indexMainDiagY) = 1.0 / tau + (myu(i, 3 : ly) * mycosTabY(3 : ly) + myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2)) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	  workingArrayY(indexSubDiagY2) = -myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
	  workingArrayY(indexMainDiagY) = 1.0 / tau
!DEC$ ENDIF
      Ay = reshape(workingArrayY, (/ly, ly/))
   else if (ApprOrder .eq. 4) then
	  workingArrayY(indexSuperDiagY4) = -myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
      workingArrayY(indexSuperDiagY3) = (8.0 * myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) + 8.0 * myu(i, 3 : ly) * mycosTabY(3 : ly)) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
      workingArrayY(indexSuperDiagY2) = -64.0 * myu(i, 3 : ly) * mycosTabY(3 : ly) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
 	  workingArrayY(indexSuperDiagY) = -(8.0 * myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) + 8.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2)) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	  workingArrayY(indexMainDiagY) = 1.0 / tau + (myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) + 64.0 * myu(i, 3 : ly) * mycosTabY(3 : ly) + 64.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2) + myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3])) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	  workingArrayY(indexSubDiagY) = -(8.0 * myu(i, 3 : ly) * mycosTabY(3 : ly) + 8.0 * myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3])) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
      workingArrayY(indexSubDiagY2) = -64.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
      workingArrayY(indexSubDiagY3) = (8.0 * myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3]) + 8.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2)) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	  workingArrayY(indexSubDiagY4) = -myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3]) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
      Ay = reshape(workingArrayY, (/ly, ly/))
   end if
    
   if (isDiagPredom(Ay) .eq. 0) then
      write (*, *) 'The matrix Ay is not diagonally predominant.'
      !res = 0
      !return
   end if
    
   bya(1) = -a(i, 1) + a(i, ly - 1)
   bya(ly) = -a(i, 2) + a(i, ly)
    
   if (ApprOrder .eq. 2) then
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  bya(2 : ly - 1) = a(i, 2 : ly - 1) * (1.0 / tau - (myu(i, 3 : ly) * mycosTabY(3 : ly) + myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2)) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)) + &
													a(i, [4 : ly, 3]) * myu(i, 3 : ly) * mycosTabY(3 : ly) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
													a(i, [ly - 2, 1 : ly - 3]) * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
	  bya(2 : ly - 1) = a(i, 2 : ly - 1) * (1.0 / tau - (myu(i, 3 : ly) * mycosTabY(3 : ly) + myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)) + &
													a(i, [4 : ly, 3]) * myu(i, 3 : ly) * mycosTabY(3 : ly) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
													a(i, [ly - 2, 1 : ly - 3]) * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
!DEC$ ENDIF
   else if (ApprOrder .eq. 4) then
	  bya(2 : ly - 1) = a(i, 2 : ly - 1) * (1.0 / tau - (myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) + 64.0 * myu(i, 3 : ly) * mycosTabY(3 : ly) + 64.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2) + myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3])) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)) + &
												a(i, [6 : ly, 3, 4, 5]) * myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) - &
												a(i, [5 : ly, 3, 4]) * (8.0 * myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) + 8.0 * myu(i, 3 : ly) * mycosTabY(3 : ly)) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
												a(i, [4 : ly, 3]) * 64.0 * myu(i, 3 : ly) * mycosTabY(3 : ly) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
												a(i, 3 : ly) * (8.0 * myu(i, [4 : ly, 3]) * mycosTabY([4 : ly, 3]) + 8.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2)) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
												a(i, 1 : ly - 2) * (8.0 * myu(i, 3 : ly) * mycosTabY(3 : ly) + 8.0 * myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3])) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
												a(i, [ly - 2, 1 : ly - 3]) * 64.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) - &
												a(i, [ly - 3, ly - 2, 1 : ly - 4]) * (8.0 * myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3]) + 8.0 * myu(i, 1 : ly - 2) * mycosTabY(1 : ly - 2)) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
												a(i, [ly - 4, ly - 3, ly - 2, 1 : ly - 5]) * myu(i, [ly - 2, 1 : ly - 3]) * mycosTabY([ly - 2, 1 : ly - 3]) / (288.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
   end if

   cflY = max(maxval(abs(Ay)) * tau, cflY)

   !DEC$ IF DEFINED (VF)
   a_next(i, :) = Ay .ix. bya
   !DEC$ ELSE
   ipivy(:) = 0
   call dgetrf(ly, ly, Ay, ly, ipivy, info)

   if (info .eq. 0) then
      call dgetrs('N', ly, 1, Ay, ly, ipivy, bya, ly, info)

      a_next(i, :) = bya
   else
      write (*, *) 'The matrix Ay is singular.'
   end if
   !DEC$ ENDIF

   ! Update
   a(i, :) = a_next(i, :)
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
	  a(lx - 1 : lx, 1 : ly) = a(1 : 2, 1 : ly)
!DEC$ IF .NOT. DEFINED (VF)
   end if
!DEC$ ENDIF
end if

!DEC$ IF .NOT. DEFINED (VF)
do j = 1, ly
   call mpi_gather(a(iFrom, j), iTo - iFrom + 1, MPI_REAL8, a(1, j), iTo - iFrom + 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
end do

call mpi_bcast(a, product(shape(a)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

end subroutine DoInY
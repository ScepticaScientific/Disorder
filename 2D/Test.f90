program DiffusionSphere

! To do: Use 'MYU_FACTOR_IS_ARRAY' if myu is a function of (lambda, phi) or use 'NOTMYU_FACTOR_IS_ARRAY' if myu is a constant
!DEC$ DEFINE NOTMYU_FACTOR_IS_ARRAY

! To do: Use 'SOLUTION_ENERGY' if you need to compute the solution's L2-norm or use 'ERROR_ENERGY' if you need to compute the solution relative error's L2-norm (i.e. ||(numerics - analytics) / analytics||_L2)
!DEC$ DEFINE SOLUTION_ENERGY

!DEC$ IF DEFINED (VF)
use linear_operators
!DEC$ ENDIF

use kronFunction
use GetMassSubroutine
use fRHSFunction

implicit none

!DEC$ IF .NOT. DEFINED (VF)
include 'mpif.h'
!DEC$ ENDIF

! Variables (general)
real, pointer :: a(:, :), rhs(:, :)

real, pointer :: x(:), y(:)

real, pointer :: xs(:, :), ys(:, :)

real, pointer :: mycosTabX(:), mysinTabX(:), mytanTabX(:)
real, pointer :: mycosTabY(:), mysinTabY(:), mytanTabY(:)

real, pointer :: coslambda(:, :), sinlambda(:, :), tanlambda(:, :)
real, pointer :: cosphi(:, :), sinphi(:, :), tanphi(:, :)

real dx, dy
integer lx, ly
real llx, lly							! These are Lx and Ly, respectively, in MatLab
integer axisPlusX, axisPlusY

real, pointer :: yNew(:)
integer lyNew
real llyNew							! This is LyNew in MatLab

real, pointer :: myu(:, :)

real weight

real R_Earth

integer ApprOrder
logical IsFirstOrderMinusPlus, IsCheck

real pi

real mass_dummy, energy_aux
real alpha

real, pointer :: Ax(:, :), bxa(:), ipivx(:)
real, pointer :: Ay(:, :), bya(:), ipivy(:)

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
real, pointer :: myu_factor(:, :)
!DEC$ ENDIF
!DEC$ IF .NOT. DEFINED (MYU_FACTOR_IS_ARRAY)
real myu_factor
!DEC$ ENDIF

integer, pointer :: indexMainDiagX(:), indexSuperDiagX(:), indexSuperDiagX2(:), indexSuperDiagX3(:), indexSuperDiagX4(:), indexSubDiagX(:), indexSubDiagX2(:), indexSubDiagX3(:), indexSubDiagX4(:)
integer, pointer :: indexMainDiagY(:), indexSuperDiagY(:), indexSuperDiagY2(:), indexSuperDiagY3(:), indexSuperDiagY4(:), indexSubDiagY(:), indexSubDiagY2(:), indexSubDiagY3(:), indexSubDiagY4(:)

real, pointer :: workingArrayX(:), workingArrayY(:)

real cflX, cflY

logical IsPlaneModel

!DEC$ IF .NOT. DEFINED (VF)
integer ierr, processID, NumberOfProcesses
!DEC$ ENDIF

common /arhs/ a, rhs
common /axes/ x, y
common /axesmat/ xs, ys
common /trigsX/ mycosTabX, mysinTabX, mytanTabX
common /trigsY/ mycosTabY, mysinTabY, mytanTabY
common /trigsXYmat/ coslambda, sinlambda, tanlambda, cosphi, sinphi, tanphi
common /steps/ dx, dy
common /limits/ llx, lly
common /sizes/ lx, ly
common /tau/ tau
common /axisPluses/ axisPlusX, axisPlusY
common /myu/ myu
common /weight/ weight
common /Radius/ R_Earth
common /specials/ ApprOrder, IsFirstOrderMinusPlus, IsCheck
common /consts/ pi
common /matvecX/ Ax, bxa, ipivx
common /diagsX/ indexMainDiagX, indexSuperDiagX, indexSuperDiagX2, indexSuperDiagX3, indexSuperDiagX4, indexSubDiagX, indexSubDiagX2, indexSubDiagX3, indexSubDiagX4, workingArrayX
common /matvecY/ Ay, bya, ipivy
common /diagsY/ indexMainDiagY, indexSuperDiagY, indexSuperDiagY2, indexSuperDiagY3, indexSuperDiagY4, indexSubDiagY, indexSubDiagY2, indexSubDiagY3, indexSubDiagY4, workingArrayY
common /CFL/ cflX, cflY
common /theGeometry/ IsPlaneModel

!DEC$ IF .NOT. DEFINED (VF)
common /mpiData/ ierr, processID, NumberOfProcesses
!DEC$ ENDIF

! Variables (partial)
real T, tau
real a_ampl

integer Nx, Ny

integer fieldNo
!logical IsDepth

integer currTime
real, pointer :: mass(:), energy(:)

real, pointer :: analytic_solution(:, :)

integer outputStep, outputToDiskStep
integer i, j
integer res
logical IsSources

integer cx

! Body
!DEC$ IF .NOT. DEFINED (VF)
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, processID, ierr)
call mpi_comm_size(MPI_COMM_WORLD, NumberOfProcesses, ierr)
!DEC$ ENDIF

pi = 3.141592653589793238D+000

! Definitions
R_Earth = 1.0 / (2.0 * pi)     ! Earth's radius

!DEC$ IF .NOT. DEFINED (VF)
if (processID .eq. 0) then
!DEC$ ENDIF

open (1, file = 'init.txt')
read (1, *) T, tau, outputStep, outputToDiskStep
read (1, *) Nx, Ny
read (1, *) IsSources
read (1, *) IsPlaneModel
read (1, *) weight
read (1, *) IsCheck, ApprOrder, IsFirstOrderMinusPlus
close (1)

!DEC$ IF .NOT. DEFINED (VF)
end if

call mpi_bcast(T, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(tau, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(Nx, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(Ny, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(IsSources, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(IsPlaneModel, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(weight, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(isCheck, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(ApprOrder, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(IsFirstOrderMinusPlus, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

!write (*, *) processID, T, tau, Nx, Ny, fieldNo, isDepth, IsPlaneModel, weight, IsCheck, ApprOrder, IsFirstOrderMinusPlus
!DEC$ ENDIF

if (IsPlaneModel .eq. 0) then
   llx = 2.0 * pi;  dx = llx / (Nx - 1)
   lly = pi / 2.0;  dy = 2.0 * lly / (Ny - 1)
else if (IsPlaneModel .eq. 1) then
   llx = 1.0;   dx = llx / (Nx - 1)
   lly = 0.25;  dy = 2.0 * lly / (Ny - 1)
end if 

!outputStep = max(int(5 * 5.0e-4 / tau), 50)

axisPlusX = 1
axisPlusY = 0

allocate(x(int((llx + dx) / dx + 1)));						lx = size(x, 1)
x = [0 : lx - 1] * dx + dx / 2.0

if (IsPlaneModel .eq. 0) then
   allocate(y(int((2.0 * lly - dy) / dy + 1)));				ly = size(y, 1)
   y = -lly + dy / 2.0 + [0 : ly - 1] * dy
else if (IsPlaneModel .eq. 1) then
   allocate(y(int((2.0 * lly + dy) / dy + 1)));				ly = size(y, 1)
   y = [0 : ly - 1] * dy
   axisPlusY = 2
end if

!
allocate(xs(lx, ly), ys(lx, ly))
xs = spread(x, 2, ly)
ys = spread(y, 1, lx)

allocate(coslambda(lx, ly), sinlambda(lx, ly), tanlambda(lx, ly))
allocate(cosphi(lx, ly), sinphi(lx, ly), tanphi(lx, ly))
coslambda = cos(xs)
sinlambda = sin(xs)
tanlambda = tan(xs)
cosphi = cos(ys)
sinphi = sin(ys)
tanphi = tan(ys)

allocate(mycosTabX(ly), mysinTabX(ly), mytanTabX(ly))
mycosTabX = cos(y)
mysinTabX = sin(y)
mytanTabX = tan(y)

! The parameters "LyNew" and "yNew" must be the same as those used in the function "SwapXY()"
if (IsPlaneModel .eq. 0) then
   llyNew = 2.0 * pi

   allocate(yNew(int((llyNew + dy) / dy + 1)));				lyNew = size(yNew, 1)
   yNew = -lly + dy / 2.0 + [0 : lyNew - 1] * dy
else if (IsPlaneModel .eq. 1) then
   llyNew = lly

   allocate(yNew(size(y, 1)));								lyNew = size(yNew, 1)
   yNew = y
end if 

allocate(mycosTabY(lyNew), mysinTabY(lyNew), mytanTabY(lyNew))
mycosTabY = abs(cos(yNew))
mysinTabY = sin(yNew)
mytanTabY = tan(yNew)

!
allocate(myu(lx, ly))

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
allocate(myu_factor(lx, ly))
!DEC$ ENDIF

!
allocate(a(lx, ly), rhs(lx, ly), analytic_solution(lx, ly))

rhs(:, :) = 0.0 * kron(sin(x) ** 2.0, sin(y) ** 2.0)

!DEC$ IF .NOT. DEFINED (VF)
if (processID .eq. 0) then
!DEC$ ENDIF

! Sources.
if (IsSources .eq. 1) then
   open (1, file = 'SourcesInit.txt')
   do i = 1, size(rhs, 1)
      read (1, *) (rhs(i, j), j = 1, size(rhs, 2))
   end do 
   close (1)
end if

!DEC$ IF .NOT. DEFINED (VF)
end if

call mpi_bcast(rhs, product(shape(rhs)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

if (IsPlaneModel .eq. 1) then
   R_Earth = 1.0

   mycosTabX(1 : ly) = 1.0
   mysinTabX(1 : ly) = 1.0
   mytanTabX(1 : ly) = 0.0

   mycosTabY(1 : lyNew) = 1.0
   mysinTabY(1 : lyNew) = 1.0
   mytanTabY(1 : lyNew) = 0.0
end if

!DEC$ IF .NOT. DEFINED (VF)
if (processID .eq. 0) then
!DEC$ ENDIF

open (1, file = 'aInit.txt')
do i = 1, size(a, 1)
   read (1, *) (a(i, j), j = 1, size(a, 2))
end do 
close (1)

open (1, file = 'myuInit.txt')
do i = 1, size(myu, 1)
   read (1, *) (myu(i, j), j = 1, size(myu, 2))
end do 
close (1)

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
myu_factor = myu
!DEC$ ENDIF

!DEC$ IF .NOT. DEFINED (VF)
end if

call mpi_bcast(a, product(shape(a)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(myu, product(shape(myu)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

!write (*, *) processID, un(lx, ly), product(shape(un)), shape(un), vn(lx, ly), a(lx, ly), U(lx, ly), V(lx, ly)
!DEC$ ENDIF

! These allocations and assignements are equal to those performed in MatLab
allocate(Ax(lx, lx), bxa(lx), ipivx(lx))

allocate(indexSuperDiagX4(size(Ax, 2) - 1 - 1))
indexSuperDiagX4 = (/(size(Ax, 1) * (j - 1) + j - 4, j = 6, size(Ax, 2)), size(Ax, 1) * (4 - 1) - 3, size(Ax, 1) * (5 - 1) - 2, size(Ax, 1) * (6 - 1) - 1/)

allocate(indexSuperDiagX3(size(Ax, 2) - 1 - 1))
indexSuperDiagX3 = (/(size(Ax, 1) * (j - 1) + j - 3, j = 5, size(Ax, 2)), size(Ax, 1) * (4 - 1) - 2, size(Ax, 1) * (5 - 1) - 1/)

allocate(indexSuperDiagX2(size(Ax, 2) - 1 - 1))
indexSuperDiagX2 = (/(size(Ax, 1) * (j - 1) + j - 2, j = 4, size(Ax, 2)), size(Ax, 1) * (4 - 1) - 1/)

allocate(indexSuperDiagX(size(Ax, 2) - 1 - 1))
indexSuperDiagX = (/(size(Ax, 1) * (j - 1) + j - 1, j = 3, size(Ax, 2))/)

allocate(indexMainDiagX(size(Ax, 2) - 1 - 1))
indexMainDiagX = (/(size(Ax, 1) * (j - 1) + j, j = 2, size(Ax, 2) - 1)/)

allocate(indexSubDiagX(size(Ax, 2) - 1 - 1))
indexSubDiagX = (/(size(Ax, 1) * (j - 1) + j + 1, j = 1, size(Ax, 2) - 2)/)

allocate(indexSubDiagX2(size(Ax, 2) - 1 - 1))
indexSubDiagX2 = (/size(Ax, 1) * (size(Ax, 2) - 2 - 1) + 2, (size(Ax, 1) * (j - 1) + j + 2, j = 1, size(Ax, 2) - 3)/)

allocate(indexSubDiagX3(size(Ax, 2) - 1 - 1))
indexSubDiagX3 = (/size(Ax, 1) * (size(Ax, 2) - 2 - 2) + 2, size(Ax, 1) * (size(Ax, 2) - 2 - 1) + 3, (size(Ax, 1) * (j - 1) + j + 3, j = 1, size(Ax, 2) - 4)/)

allocate(indexSubDiagX4(size(Ax, 2) - 1 - 1))
indexSubDiagX4 = (/size(Ax, 1) * (size(Ax, 2) - 2 - 3) + 2, size(Ax, 1) * (size(Ax, 2) - 2 - 2) + 3, size(Ax, 1) * (size(Ax, 2) - 2 - 1) + 4, (size(Ax, 1) * (j - 1) + j + 4, j = 1, size(Ax, 2) - 5)/)

allocate(Ay(lyNew, lyNew), bya(lyNew), ipivy(lyNew))

allocate(indexSuperDiagY4(size(Ay, 2) - 1 - 1))
indexSuperDiagY4 = (/(size(Ay, 1) * (j - 1) + j - 4, j = 6, size(Ay, 2)), size(Ay, 1) * (4 - 1) - 3, size(Ay, 1) * (5 - 1) - 2, size(Ay, 1) * (6 - 1) - 1/)

allocate(indexSuperDiagY3(size(Ay, 2) - 1 - 1))
indexSuperDiagY3 = (/(size(Ay, 1) * (j - 1) + j - 3, j = 5, size(Ay, 2)), size(Ay, 1) * (4 - 1) - 2, size(Ay, 1) * (5 - 1) - 1/)

allocate(indexSuperDiagY2(size(Ay, 2) - 1 - 1))
indexSuperDiagY2 = (/(size(Ay, 1) * (j - 1) + j - 2, j = 4, size(Ay, 2)), size(Ay, 1) * (4 - 1) - 1/)

allocate(indexSuperDiagY(size(Ay, 2) - 1 - 1))
indexSuperDiagY = (/(size(Ay, 1) * (j - 1) + j - 1, j = 3, size(Ay, 2))/)

allocate(indexMainDiagY(size(Ay, 2) - 1 - 1))
indexMainDiagY = (/(size(Ay, 1) * (j - 1) + j, j = 2, size(Ay, 2) - 1)/)

allocate(indexSubDiagY(size(Ay, 2) - 1 - 1))
indexSubDiagY = (/(size(Ay, 1) * (j - 1) + j + 1, j = 1, size(Ay, 2) - 2)/)

allocate(indexSubDiagY2(size(Ay, 2) - 1 - 1))
indexSubDiagY2 = (/size(Ay, 1) * (size(Ay, 2) - 2 - 1) + 2, (size(Ay, 1) * (j - 1) + j + 2, j = 1, size(Ay, 2) - 3)/)

allocate(indexSubDiagY3(size(Ay, 2) - 1 - 1))
indexSubDiagY3 = (/size(Ay, 1) * (size(Ay, 2) - 2 - 2) + 2, size(Ay, 1) * (size(Ay, 2) - 2 - 1) + 3, (size(Ay, 1) * (j - 1) + j + 3, j = 1, size(Ay, 2) - 4)/)

allocate(indexSubDiagY4(size(Ay, 2) - 1 - 1))
indexSubDiagY4 = (/size(Ay, 1) * (size(Ay, 2) - 2 - 3) + 2, size(Ay, 1) * (size(Ay, 2) - 2 - 2) + 3, size(Ay, 1) * (size(Ay, 2) - 2 - 1) + 4, (size(Ay, 1) * (j - 1) + j + 4, j = 1, size(Ay, 2) - 5)/)

allocate(workingArrayX(product(shape(Ax))), workingArrayY(product(shape(Ay)))) 

allocate(mass(int(T / tau)), energy(int(T / tau)))

res = 1

cx = 0

!DEC$ IF .NOT. DEFINED (VF)
if (processID .eq. 0) then
!DEC$ ENDIF

call GetMass(mass(1), energy(1))
call SaveToDisk(cx)

!DEC$ IF .NOT. DEFINED (VF)
end if

call mpi_barrier(MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

alpha = 0.0
!DEC$ IF .NOT. DEFINED (MYU_FACTOR_IS_ARRAY)
myu_factor = 1.0e-04
!DEC$ ENDIF

! Computing
do currTime = 1, int(T / tau)
   cflX = -1.0D+200
   cflY = -1.0D+200

   myu = myu_factor * a ** alpha
   
   Ax = 0.0
   bxa = 0.0

   ! Computing in "x"
   call DoInX(tau / 2.0, res)

   if (res .eq. 0) then
      pause
      stop
   end if

   if (IsPlaneModel .eq. 0) then
      call SwapXY()
   end if 

   Ay = 0.0
   bya = 0.0

   ! Computing in "y"
   call DoInY(tau / 2.0, res)

   if (res .eq. 0) then
      pause
      stop
   end if
    
   if (IsPlaneModel .eq. 0) then
      call SwapYX()
   end if 

   ! Sources
   if (IsSources .eq. 1) then
	  a = a + tau * rhs
   else if (IsSources .eq. 2) then
!DEC$ IF DEFINED (ERROR_ENERGY)
      ! Third experiment for NMPDE (2010-2011)
	  !a = a + tau * fRHS(1.0, 2.0, 0.0, -pi / 3.0, pi, currTime * tau - tau / 2.0)

      ! First experiment for SIMULTECH 2011
	  !a = a + tau * fRHS(1.0, 1.0, 0.0, 0.0, pi, currTime * tau - tau / 2.0)

	  ! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
	  a = a + tau * fRHS(alpha, 9.0, 5.0 * alpha, 3.0, -2.5, 50.0, myu_factor, currTime * tau - tau / 2.0)

	  ! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
	  !a = a + tau * fRHS(alpha, 7.0, 1.0, 0.5, 0.0 * pi / 4.0, 15.0, 0.7, 0.0 * pi / 4.0, 4.0, 10.0, 100.0, myu_factor, currTime * tau - tau / 2.0)
!DEC$ ENDIF
      ! NOVA Publishers, third experiment (combustion); Applied Mathematics and Computation, fifth experiment (combustion); WCE2013, temperature wave; WCE2013, combustion
	  a = a + tau * 1.0e+01 * (a - a ** 3.0)	! for 'myu_factor = 1.0e-04'
	  !a = a + tau * 3.9e+00 * a ** 3.0			! for 'myu_factor = 1.0e-03'
	  !a = a + tau * 4.1e+00 * a ** 4.0			! for 'myu_factor = 1.0e-03'
	  !a = a + tau * 4.5e+00 * a ** 2.0			! for 'myu_factor = 1.0e-03'
   end if

   if (IsPlaneModel .eq. 0) then
      call SwapXY()
   end if 

   Ay = 0.0
   bya = 0.0

   ! Computing in "y"
   call DoInY(tau / 2.0, res)

   if (res .eq. 0) then
      pause
      stop
   end if
    
   if (IsPlaneModel .eq. 0) then
      call SwapYX()
   end if 

   Ax = 0.0
   bxa = 0.0

   ! Computing in "x"
   call DoInX(tau / 2.0, res)

   if (res .eq. 0) then
      pause
      stop
   end if

!DEC$ IF DEFINED (ERROR_ENERGY)
   ! Third experiment for NMPDE (2010-2011)
   !analytic_solution = theT(1.0, 2.0, 0.0, -pi / 3.0, pi, currTime * tau)

   ! First experiment for SIMULTECH 2011
   !analytic_solution = theT(1.0, 1.0, 0.0, 0.0, pi, currTime * tau)

   ! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
   analytic_solution = theT(9.0, 5.0 * alpha, 3.0, -2.5, 50.0, currTime * tau)
      
   ! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
   !analytic_solution = theT(7.0, 1.0, 0.5, 0.0 * pi / 4.0, 15.0, 0.7, 0.0 * pi / 4.0, 4.0, 10.0, 100.0, currTime * tau)
!DEC$ ENDIF

   ! Computing the total mass and energy
   !DEC$ IF .NOT. DEFINED (VF)   
   if (processID .eq. 0) then
   !DEC$ ENDIF

!DEC$ IF DEFINED (ERROR_ENERGY)
   a = a - analytic_solution
!DEC$ ENDIF
   call GetMass(mass_dummy, energy(currTime))
!DEC$ IF DEFINED (ERROR_ENERGY)
   a = a + analytic_solution
!DEC$ ENDIF
   call GetMass(mass(currTime), energy_aux)
!DEC$ IF DEFINED (ERROR_ENERGY)
   energy(currTime) = energy(currTime) / energy_aux	! Actually, this yields (numerics - analytics) / numerics, which is not exactly the relative error, as well. Indeed, it should be (numerics - analytics) / analytics instead.
!DEC$ ENDIF

   if (mod(currTime, outputStep) .eq. 0) then
      write (*, '(F15.10, A, D20.10, D20.10)') currTime * tau, ': ', mass(currTime), energy(currTime)
	  write (*, '(D20.10, A, D20.10)') cflX, ' ', cflY
   end if

   if (mod(currTime, outputToDiskStep) .eq. 0) then
      cx = cx + 1
      call SaveToDisk(cx)
   end if

   !DEC$ IF .NOT. DEFINED (VF)   
   end if

   !DEC$ ENDIF
end do

!DEC$ IF .NOT. DEFINED (VF)   
if (processID .eq. 0) then
!DEC$ ENDIF

! Saving the output
! Axes
open (1, file = 'x.tx2')
write (1, '(D20.10)') x
close (1)

open (1, file = 'y.tx2')
write (1, '(D20.10)') y
close (1)

! Mass
open (1, file = 'mass.tx2')
!write (1, '(D20.10)') mass
write (1, '(E20.10E3)') mass
close (1)

! Energy
open (1, file = 'energy.tx2')
!write (1, '(D20.10)') energy
write (1, '(E20.10E3)') energy
close (1)

pause

!DEC$ IF .NOT. DEFINED (VF)   
end if
!DEC$ ENDIF

deallocate(mass, energy)

deallocate(workingArrayX, workingArrayY) 
deallocate(Ay, bya, ipivy)
deallocate(Ax, bxa, ipivx)

deallocate(indexSuperDiagY4, indexSuperDiagY3, indexSuperDiagY2, indexSuperDiagY, indexMainDiagY, indexSubDiagY, indexSubDiagY2, indexSubDiagY3, indexSubDiagY4)
deallocate(indexSuperDiagX4, indexSuperDiagX3, indexSuperDiagX2, indexSuperDiagX, indexMainDiagX, indexSubDiagX, indexSubDiagX2, indexSubDiagX3, indexSubDiagX4)

deallocate(a, rhs, analytic_solution)

deallocate(myu)

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
deallocate(myu_factor)
!DEC$ ENDIF

deallocate(yNew)

deallocate(coslambda, sinlambda, tanlambda)
deallocate(cosphi, sinphi, tanphi)

deallocate(mycosTabY, mysinTabY, mytanTabY)
deallocate(mycosTabX, mysinTabX, mytanTabX)

deallocate(xs, ys)
deallocate(x, y)

!DEC$ IF .NOT. DEFINED (VF)
call mpi_finalize(ierr)
!DEC$ ENDIF

end program DiffusionSphere


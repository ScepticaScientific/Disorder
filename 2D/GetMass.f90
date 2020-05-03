module GetMassSubroutine

contains

	subroutine GetMass(mass, energy)

	implicit none

	! Variables (general)
	real, pointer :: a(:, :), rhs(:, :)

	real, pointer :: mycosTabX(:), mysinTabX(:), mytanTabX(:)

	real dx, dy
	integer lx, ly
	integer axisPlusX, axisPlusY

	real tau

	real R_Earth

	common /arhs/ a, rhs
	common /trigsX/ mycosTabX, mysinTabX, mytanTabX
	common /steps/ dx, dy
	common /sizes/ lx, ly
	common /tau/ tau
	common /axisPluses/ axisPlusX, axisPlusY
	common /Radius/ R_Earth

	! Variables (partial)
	real mass, energy

	integer j, i

	! Body
	mass = 0.0
	energy = 0.0

	do j = 1, ly - axisPlusY
	   mass = mass + mycosTabX(j) * sum(a(1 : lx - axisPlusX - 1, j))
	   energy = energy + mycosTabX(j) * sum(a(1 : lx - axisPlusX - 1, j) ** 2.0)
	end do

	mass = mass * R_Earth * R_Earth * dx * dy
	energy = sqrt(energy * R_Earth * R_Earth * dx * dy)

	end subroutine GetMass

end module GetMassSubroutine
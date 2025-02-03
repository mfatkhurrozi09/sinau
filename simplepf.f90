program simplepf

! modified

! Simple Potential Flow using Finite Difference Methdod
! governing equation is Laplace Equation , del^2(phi) = 0.
! solution is of explicit type

! update 2025 01 31
! Computational Hydrodynamics, KL3252, class of 2025
! week 2 - 7
	
implicit none

! Parameters
integer, parameter :: nx=20, ny=20
real, parameter :: dx=1., dy=1.

! variabel preparation
real :: phi(nx,ny), phin(nx,ny)
real :: x,y,r,pii
integer :: i,j,m,iter,maxit
real :: tolr,error
character*40 fo,foo
real :: phii(10000)

 ! Convergence parameters
 maxit=10000 
 tolr=1.0e-4

print*, '-- KL3252 Komputasi Hidrodinamika --'
print*
print*, 'SIMFLOW started ...'
call sleep(1)
print*

fo = 'field.txt'
open(15,file=fo,status='unknown')
print*, 'Output file is ',fo

foo = 'obs.txt'
open(16,file=foo,status='unknown')


 ! Initialize potential field
 phi=0.0
 
! Apply boundary conditions at left
do i=1, ny
	phi(i,1)=0.5
end do

! Iterative solver for Laplace equation
do iter = 1, maxit
	phin = phi

    do j=2,ny-1
		do i=2,nx-1
			phin(i,j)=0.25*(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1))
			if (j>(ny/2).and.i>(nx/2)) then
				phin(i,j)=0.6 
			end if 
		end do
    end do


	do j=2,ny-1
		phin(j,nx)=0. !-0.01+phin(j,nx-1)
		if (j>(ny/2)) then
			m=(nx/2)-1
			phin(j,nx/2)=phi(j,m)
		end if
		
	enddo
				
	do i=2,nx-1			
		phin(1,i)=phin(2,i) 
		phin(ny,i)=phin(ny-1,i)
		if (i>(nx/2)) then
			m=(ny/2)-1
			phin(ny/2,i)=phi(m,i)
		end if		
	enddo 
				
	! Compute error
	error=maxval(abs(phin - phi))
	pii=phin(3,18)-phi(3,18)
	write(16,*) pii

	! Update potential field
	phi=phin

	if (error<tolr) then
		!call sleep(1)
		print*, 'Converged after iterations:', iter
		exit
	end if
	
	! Print selected observation point
	!write(16,*) phin(3,18)
	
 end do

if (iter == maxit) then
    print *, 'Maximum iterations reached. Solution may not have converged.'
end if

print *, '-- Simulation completed. Matur suwun --'
  
do i=2,nx-1
	write(15,11) (phi(i,j),j=2,ny-1)
enddo
11	  format(500f10.3)
 
end program simplepf

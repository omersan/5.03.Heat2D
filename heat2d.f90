!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> ADI Method for solving unsteady heat equation in 2D
!     du/dt = a*(d2u/dx2 + d2u/dy2) + q(x,y) 
!     Approximate factorization + Crank-Nicolson
!     Drichlet b.c.
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 13, 2015
!-----------------------------------------------------------------------------!

program heat2d
implicit none
integer::i,j,k,nx,ny,nt,ns,nf
real*8 ::dx,dy,dt,x0,xL,y0,yL,t,Tmax,alpha
real*8,allocatable ::u(:,:),x(:),y(:)

!Domain
x0 =-1.0d0 !left
xL = 1.0d0 !right

y0 =-1.0d0 !bottom
yL = 1.0d0 !up

!number of points
nx = 20
ny = 20

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do


!maximum time desired
Tmax = 1.0d0

!diffusion constant:
alpha = 1.0d0

!time step
dt = 0.05d0

!number of points in time
nt = nint(Tmax/dt)

!number of snapshots to plot
ns = nt

!frequency for plotting
nf = nint(dfloat(nt)/dfloat(ns))

!u: convected variable 
allocate(u(0:nx,0:ny))

!initial condition
t = 0.0d0
do j=0,ny
do i=0,nx
u(i,j) = 0.0d0
end do
end do

!Plot initial condition
open(18,file='u.plt')
write(18,*) 'variables ="x","y","u"'
write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
do j=0,ny
do i=0,nx
write(18,*) x(i),y(j),u(i,j)
end do
end do

!time integration
do k=1,nt
    
    call CN(nx,ny,dx,dy,dt,t,alpha,x,y,u)
      
    !update t
    t = t+dt 

    !plot field
    if (mod(k,nf).eq.0) then
	write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
	do j=0,ny
	do i=0,nx
	write(18,*) x(i),y(j),u(i,j)
	end do
	end do
    end if

    
	print*,k,t,maxval(u)


end do
  
close(18)


100 format(a16,i8,a4,i8,a10,f10.4,a3)
end

!-----------------------------------------------------------------------------!
!Crank-Nicolson(CN) scheme with ADI
!-----------------------------------------------------------------------------!
subroutine CN(nx,ny,dx,dy,dt,t,alpha,x,y,u)
implicit none
integer::nx,ny,i,j
real*8 ::dx,dy,dt,t,alpha,bx,by,f,ua,ub
real*8 ::u(0:nx,0:ny),x(0:nx),y(0:ny)
real*8,allocatable ::a(:),b(:),c(:),r(:),q(:)
real*8 ::s(0:nx,0:ny),p(0:nx,0:ny),z(0:nx,0:ny)

do j=0,ny
do i=0,nx
s(i,j) = 0.0d0
p(i,j) = 0.0d0
z(i,j) = 0.0d0
end do
end do

bx = 0.5d0*alpha*dt/(dx*dx)
by = 0.5d0*alpha*dt/(dy*dy)

!compute source term in known step
do j=1,ny-1
do i=0,nx
s(i,j) = u(i,j) + by*(u(i,j+1)-2.0d0*u(i,j)+u(i,j-1))
end do
end do

do j=1,ny-1
do i=1,nx-1
p(i,j) = s(i,j) + bx*(s(i+1,j)-2.0d0*s(i,j)+s(i-1,j)) &
       + 0.5d0*dt*(f(t,x(i),y(j))+f(t+dt,x(i),y(j))) 
end do
end do


!x-sweep to compute intermediate values:
do j=1,ny-1

  
	!Build coefficient matrix:
	allocate(a(1:nx-1),b(1:nx-1),c(1:nx-1),r(1:nx-1),q(1:nx-1))

	do i=1,nx-1
	a(i) = -bx
	b(i) = (1.0d0+2.0d0*bx)
	c(i) = -bx
	r(i) = p(i,j) 
	end do
    
	!apply boundary conditions
    ua = u(0,j)  - by*(u(0,j+1)-2.0d0*u(0,j)+u(0,j-1))
    ub = u(nx,j) - by*(u(nx,j+1)-2.0d0*u(nx,j)+u(nx,j-1))
    
	r(1)   = r(1) - a(1)*ua        !b.c.
	r(nx-1) = r(nx-1) - c(nx-1)*ub !b.c.
    
	call tdma(a,b,c,r,q,1,nx-1)

	!assign solutions for as z
	do i=1,nx-1
	z(i,j)=q(i)
	end do
    z(0,j) =ua
    z(nx,j)=ub

    deallocate(a,b,c,r,q)

end do


!y-sweep to compute final solution:
do i=1,nx-1

  
	!Build coefficient matrix:
	allocate(a(1:ny-1),b(1:ny-1),c(1:ny-1),r(1:ny-1),q(1:ny-1))

	do j=1,ny-1
	a(j) = -by
	b(j) = (1.0d0+2.0d0*by)
	c(j) = -by
	r(j) = z(i,j) 
	end do
    
	!apply boundary conditions
    ua = u(i,0)  
    ub = u(i,ny) 
    
	r(1)   = r(1) - a(1)*ua        !b.c.
	r(ny-1) = r(ny-1) - c(ny-1)*ub !b.c.
    
	call tdma(a,b,c,r,q,1,ny-1)

	!assign solutions for as z
	do j=1,ny-1
	u(i,j)=q(j)
	end do

    deallocate(a,b,c,r,q)
end do



end 

!------------------------------------------------------------------!
!Tridiagonal matrix algorithm (TDMA)
!Thomas algorithm
!solution tridiagonal systems
!a: lower diagonal
!b: main diagonal
!c: upper diagonal
!r: source vector
!x: solution vector
!   for indices s(start) to e(end)
!   i: s,s+1,s+2, ....,i,....,e 
!
!Note: a(s) and c(e) are dummy coefficients, not used.
!------------------------------------------------------------------!

subroutine tdma(a,b,c,r,x,s,e)
implicit none
integer s,e,i
real*8, dimension(s:e) ::a,b,c,r,x    

! forward elimination phase
do i=s+1,e
b(i) = b(i) - a(i)/b(i-1)*c(i-1)
r(i) = r(i) - a(i)/b(i-1)*r(i-1)
end do
! backward substitution phase 
x(e) = r(e)/b(e)
do i=e-1,s,-1
x(i) = (r(i)-c(i)*x(i+1))/b(i)
end do

return
end


!-----------------------------------------------------------------------------!
!source term (nonhomogenous forcing term in heat equation)
!-----------------------------------------------------------------------------!
real*8 function f(t,x,y)
implicit none
real*8::t,x,y
f = 2.0d0*(2.0d0-x*x-y*y)
end











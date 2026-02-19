!NON-TIME DEPENDENT BILLIARD

module parametros
real*8, parameter::r0=1.0d0
integer, parameter::m=3
real*8, parameter::eps=0.1d0
real*8, parameter::pi=3.141592653589793d0
real*8, parameter::tolerancia=1e-12
integer, parameter::ipasso=1d3
integer, parameter::itpasso=2d2
integer, parameter::bic=1d6
integer, parameter::conto=10
end module parametros


program bilhar
use parametros
implicit none

integer:: i, j, n, k
real*8::alpha, phi, theta, thetai, phii, alphai, r, ri, dr, dri, x, dx, y, dy, xi, dxi, yi, dyi, solution(4), daa, dtt, alpha0, &
& passoa1, passoa2, ai, af, bi, bf
open(1,file="phase_espace.dat", status="unknown")
!_____________________________________________________________________

ai=0.0d0
af=pi
passoa1=dabs(af-ai)/dble(conto)
bi=0.0d0
bf=2*pi
passoa2=dabs(bf-bi)/dble(conto)
solution=0.0d0
!_____________________________________________________________________

do n=0, conto
alpha0=ai+n*passoa1
t=0.0d00
v=0.5d0
print*, n
do k=0, conto
theta=bi+k*passoa2
alpha=alpha0
print*, theta, alpha
do i=1, 1000
  solution=0.0d0
  r=r0/(1.0d0+eps*dcos(m*theta))
  dr=(r0*eps*m*dsin(m*theta))/((1.0d0+eps*dcos(m*theta))**2)
  x=r*dcos(theta)
  y=r*dsin(theta)
  dx=dr*dcos(theta)-y
  dy=dr*dsin(theta)+x
  phi = datan2(dy, dx)

if(phi<0.0d0)then
  phi=phi+2*pi
end if

  call root(x, y, alpha, phi, solution)
!print*, solution
  if(solution(3)==0.0d0) then
    call cup(thetai, x, y, solution)
  else 
    call cake(x, y, thetai, theta, solution)
!print*, solution, thetai
!read*, 
  end if
  
!print*, solution
!print*, thetai
  ri=r0/(1.0d0+eps*dcos(m*thetai))
  dri=(r0*eps*m*dsin(m*thetai))/((1.0d0+eps*dcos(m*thetai))**2)
  xi=ri*dcos(thetai)
  yi=ri*dsin(thetai)
  dxi=dri*dcos(thetai)-yi
  dyi=dri*dsin(thetai)+xi
  phii = datan2(dyi, dxi)

if(phii<0.0d0)then
  phii=phii+2*pi
end if

  alphai=phii-(alpha+phi)
  alphai = mod(alphai, pi)
  
  if (alphai < 0.0d0) then 
    alphai = alphai + pi
  end if
  
    write(1,*)thetai, alphai


  theta=thetai
  alpha=alphai
  !print*, i

end do

end do
end do

close(1)

end program bilhar

!_____________________________________________________________________

subroutine root(x, y, alpha, phi, solution)

use parametros
implicit none
integer::i, j, k
real*8, intent(in)::x, y, alpha, phi
real*8, intent(inout)::solution(4)
real*8::a, b, c, fa, fb, fc, tpasso, rc, ra, rb

k=0

tpasso=(2*pi)/dble(ipasso)
do i=1, ipasso

  a=(i-1)*tpasso
  b=i*tpasso
  ra=r0/(1.0d0+eps*dcos(m*a))
fa = (ra*dsin(a) - y) * dcos(alpha+phi) - (ra*dcos(a) - x) * dsin(alpha+phi)
  !fa=ra*dsin(a)-y-dtan(alpha+phi)*(ra*dcos(a)-x)
  rb=r0/(1.0d0+eps*dcos(m*b))
fb = (rb*dsin(b) - y) * dcos(alpha+phi) - (rb*dcos(b) - x) * dsin(alpha+phi)
  !fb=rb*dsin(b)-y-dtan(alpha+phi)*(rb*dcos(b)-x)

  if(fa*fb<0.0d0)then
    bics: do j=1,bic
      c=(a+b)/2.0d0
      rc=r0/(1.0d0+eps*dcos(m*c))
fc = (rc*dsin(c) - y) * dcos(alpha+phi) - (rc*dcos(c) - x) * dsin(alpha+phi)
      !fc=rc*dsin(c)-y-dtan(alpha+phi)*(rc*dcos(c)-x)

      if(fc*fa<0.0d0)then
        b=c
        rb=rc
        fb=fc
      else 
        a=c
        ra=rc
        fa=fc
      end if


      if(dabs(a-b)<=tolerancia) then
        k=k+1
        solution(k)=a
        exit bics
      else if(k==4)then
      return
      end if

      end do bics

    end if

  end do

  end subroutine

!

subroutine cup(thetai, x, y, solution)

use parametros
implicit none
integer::i
real*8, intent(out)::thetai
real*8, intent(in)::x, y, solution(4)
real*8::theta, xn, yn, r, dist

do i=1, 2

  theta=solution(i)
  r=r0/(1.0d0+eps*dcos(m*theta))
  xn=r*dcos(theta)
  yn=r*dsin(theta)
  dist=dsqrt((x-xn)**2+(y-yn)**2)
  
  if(dist>=tolerancia) then
    thetai=solution(i)
  end if
  
end do

end subroutine

!
subroutine cake(x, y, thetai, theta, solution)
use parametros
implicit none
integer::i, j
real*8::ac, x0, y0, r, v(4, 2), px, py, dist, xr, yr, thetar
real*8, intent(in)::x, y, theta
real*8, intent(out):: thetai
real*8, intent(inout)::solution(4)

v=0.0d0

do i=1,4

if(dabs(theta-solution(i))>tolerancia .and. solution(i) /= 0.0d0) then
  r=r0/(1.0d0+eps*dcos(m*solution(i)))
  x0=r*dcos(solution(i))
  y0=r*dsin(solution(i))
  v(i, 1)=x0-x
  v(i, 2)=y0-y
end if

end do

do i=1,4

if(v(i,1) /= 0.0 .or. v(i,2) /= 0.0d0) then
px=1.0d0/dble(itpasso)
py=1.0d0/dble(itpasso)

fronteira: do j=0, itpasso
xr=x+dble(j)*px*v(i,1)
yr=y+dble(j)*py*v(i,2)
thetar=datan2(yr,xr)
r=r0/(1.0d0+eps*dcos(m*thetar))
dist=dsqrt(xr**2+yr**2)

if(dist>r+tolerancia)then
solution(i)=0.0d0
exit fronteira
end if

end do fronteira
end if
end do

do i=1,4

if(solution(i) /= 0.0d0) then
thetai=solution(i)
end if

end do

end


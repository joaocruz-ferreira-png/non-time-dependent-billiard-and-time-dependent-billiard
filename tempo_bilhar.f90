module parametros
real*8, parameter::r0=1.0d0
integer, parameter::m=3
real*8, parameter::eps=0.3d0
real*8, parameter::eta=0.0d0
real*8, parameter::pi=3.141592653589793d0
real*8, parameter::tolerancia=1e-8
integer, parameter::ipasso=1d3
integer, parameter::itpasso=2d2
integer, parameter::refinamento=1d6
integer, parameter::conto=25
end module parametros


program bilhar
use parametros
implicit none

integer:: i, j, n, k
real*8::alpha, phi, theta, thetai, phii, alphai, r, ri, dr, dri, x, dx, y, dy, xi, dxi, yi, dyi, solution(4), daa, dtt, alpha0, &
& passoa1, passoa2, ai, af, bi, bf, tempo, t, ti, v, vi, xp, yp, tp, tempoi, delx, dely, tpi
open(1,file="phase_espace.dat", status="unknown")
!_____________________________________________________________________

ai =5e-3
af = pi -5e-3
passoa1=dabs(af-ai)/dble(conto)
bi=0.0d0
bf=2*pi
passoa2=dabs(bf-bi)/dble(conto)
!_____________________________________________________________________

do n=0, conto
alpha=ai+n*passoa1
print*, n
do k=0, conto
theta=bi+k*passoa2
t=0.0d0
v=0.9d0
print*, k	
 do i=1, 1000
    tempo=(1.0d0+eta*dcos(t))
    r=(r0*tempo)/(1.0d0+eps*dcos(m*theta))
    dr=(r0*eps*m*tempo*dsin(m*theta))/((1.0d0+eps*dcos(m*theta))**2)
    x=r*dcos(theta)
    y=r*dsin(theta)
    dx=dr*dcos(theta)-y
    dy=dr*dsin(theta)+x
    phi=datan2(dy, dx)
    if(phi<0.0d0)then
      phi=phi+2*pi
    end if
    call root(x, y, r, phi, alpha, t, thetai, v, tpi)
    tempoi=(1.0d0+eta*dcos(tpi))
    ri=(r0*tempoi)/(1.0d0+eps*dcos(m*thetai))
    xi=ri*dcos(thetai)
    yi=ri*dsin(thetai)
    delx=xi-x
    dely=yi-y
    ti=t+(dsqrt(delx**2+dely**2))/(dabs(v))
    tempoi=(1.0d0+eta*dcos(ti))
    dri=(r0*eps*m*tempoi*dsin(m*thetai))/((1.0d0+eps*dcos(m*thetai))**2)
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
    if(thetai<0.0d0) thetai=thetai+2*pi
    write(1,*)thetai, alphai
    theta=thetai
    alpha=alphai
    t=ti
  end do

end do
end do


close(1)

end program bilhar

!_____________________________________________________________________

subroutine root(x, y, r, phi, alpha, t, thetai, v, tpi)

use parametros
implicit none
integer::i, j
real*8, intent(in)::x, y, alpha, phi, r, t, v
real*8, intent(out)::thetai, tpi
real*8::passot, ypa, xpa, thetapa, tpa, ypb, xpb, thetapb, tpb, ra, rb, tempoa, tempob, fa, fb, fc, tpc, ypc, xpc, tempoc, thetapc, &
& rc

!passot=dabs(((r0*(1.0d0+eta))/(1.0d0-eps))-((r0*(1.0d0-eta))/(1.0d0-eps)))/ipasso
passot = (3.0d0 * r0 / dabs(v)) /dble(ipasso)
!passot=3.0d0*r0/dble(ipasso)
do i=1, ipasso
  tpa=t+i*passot
  tpb=t+(i-1)*passot
  tempoa=1.0d0+eta*dcos(tpa)
  tempob=1.0d0+eta*dcos(tpb)
  xpa=x+dabs(v)*dcos(phi+alpha)*(tpa-t)
  ypa=y+dabs(v)*dsin(phi+alpha)*(tpa-t)
  thetapa=datan2(ypa, xpa)
  xpb=x+dabs(v)*dcos(phi+alpha)*(tpb-t)
  ypb=y+dabs(v)*dsin(phi+alpha)*(tpb-t)
  thetapb=datan2(ypb, xpb)
  ra=(r0*tempoa)/(1.0d0+eps*dcos(m*thetapa))
  rb=(r0*tempob)/(1.0d0+eps*dcos(m*thetapb))
  fa=ra-dsqrt(xpa**2+ypa**2)
  fb=rb-dsqrt(xpb**2+ypb**2)
  if(fa*fb <= 0.0d0 .and. (tpb - t) > 1.0d-8) then
    do j=1,refinamento
       tpc=(tpa+tpb)/2.0d0
       tempoc=1.0d0+eta*dcos(tpc)
       xpc=x+dabs(v)*dcos(phi+alpha)*(tpc-t)
       ypc=y+dabs(v)*dsin(phi+alpha)*(tpc-t)
       thetapc=datan2(ypc, xpc)
       rc=(r0*tempoc)/(1.0d0+eps*dcos(m*thetapc))
       fc=rc-dsqrt(xpc**2+ypc**2)
       if(fc*fa<0.0d0) then
         tpb=tpc
         thetapb=thetapc
         fb=fc
       else 
         tpa=tpc
         thetapa=thetapc
         fa=fc
       end if
       if(dabs(tpa-tpb)<tolerancia) then
         thetai=(thetapa+thetapb)/2.0d0
         tpi=(tpa+tpb)/2.0d0
         return
       end if
    end do
  end if
end do
end subroutine
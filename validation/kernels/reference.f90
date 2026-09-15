program mizuroute_kernel_reference
  use iso_fortran_env, only: real64
  implicit none
  integer, parameter :: dp=real64, n=20, nt=120
  real(dp), parameter :: L=4500._dp, S=0.002_dp, b=15._dp, zc=0.5_dp, mann=0.035_dp
  real(dp), parameter :: dt=3600._dp
  real(dp) :: qkw(n), qdw(n), vkw, vdw, vmc, iprev, oprev
  real(dp) :: qin, qlat, okw, odw, omc
  integer :: t

  qkw=0._dp; qdw=0._dp
  vkw=0._dp; vdw=0._dp; vmc=0._dp
  iprev=0._dp; oprev=0._dp

  open(unit=10,file='kernel_reference.csv',status='replace',action='write')
  write(10,'(A)') 't,qin,qlat,kw_qout,kw_vol,mc_qout,mc_vol,dw_qout,dw_vol'
  do t=1,nt
    qin = 18._dp + 85._dp*exp(-((real(t,dp)-30._dp)/9._dp)**2) &
                + 55._dp*exp(-((real(t,dp)-78._dp)/14._dp)**2) &
                + 4._dp*sin(real(t,dp)*0.31_dp)
    qlat = 1.5_dp + 0.8_dp*(1._dp + sin(real(t,dp)*0.17_dp))
    call kwe_step(qkw,vkw,qin,qlat,okw)
    call mc_step(vmc,iprev,oprev,qin,qlat,omc)
    call dw_step(qdw,vdw,qin,qlat,odw)
    write(10,'(I0,8(",",ES24.16E3))') t,qin,qlat,okw,vkw,omc,vmc,odw,vdw
  end do
  close(10)

contains

  pure function btop(y) result(bt)
    real(dp), intent(in) :: y
    real(dp) :: bt
    bt = b + 2._dp*y*zc
  end function

  pure function pwet(y) result(p)
    real(dp), intent(in) :: y
    real(dp) :: p
    p = b + 2._dp*y*sqrt(1._dp+zc*zc)
  end function

  pure function area(y) result(a)
    real(dp), intent(in) :: y
    real(dp) :: a
    a = y*(b+zc*y)
  end function

  function flow_depth(q) result(y)
    real(dp), intent(in) :: q
    real(dp) :: y, y0, coef1, coef2, a, bt, p, h, dhdy, err
    integer :: it
    if (q <= 1.e-50_dp) then
      y=0._dp; return
    end if
    coef1=(sqrt(S)/mann/q)**3
    coef2=2._dp*sqrt(zc*zc+1._dp)
    y0=(1._dp/(coef1*b**3))**(1._dp/5._dp)
    err=100._dp
    do it=1,100
      a=area(y0); bt=btop(y0); p=pwet(y0)
      h=coef1*a**5/p**2-1._dp
      dhdy=coef1*(5._dp*a**4*bt*p-2._dp*coef2*a**5)/p**3
      y=y0-h/dhdy
      err=abs((y-y0)/y)
      y0=y
      if (err <= 0.005_dp) exit
    end do
  end function

  function cel(q,y) result(c)
    real(dp), intent(in) :: q,y
    real(dp) :: c,a,p,bt,sf
    if (y <= 0._dp) then
      c=0._dp; return
    end if
    a=area(y); p=pwet(y); bt=btop(y)
    sf=(q*mann/a/(a/p)**(2._dp/3._dp))**2
    c=(5._dp/3._dp)*sf**0.3_dp*q**0.4_dp/bt**0.4_dp/mann**0.6_dp
  end function

  function diffu(q,y) result(d)
    real(dp), intent(in) :: q,y
    real(dp) :: d,a,p,bt,sf
    if (y <= 0._dp) then
      d=0._dp; return
    end if
    a=area(y); p=pwet(y); bt=btop(y)
    sf=(q*mann/a/(a/p)**(2._dp/3._dp))**2
    d=abs(q)/sf/bt/2._dp
  end function

  subroutine tdma(nn, mat, rhs, x)
    integer, intent(in) :: nn
    real(dp), intent(in) :: mat(nn,3), rhs(nn)
    real(dp), intent(out) :: x(nn)
    real(dp) :: dd(nn), b1(nn), coef
    integer :: i
    dd=mat(:,2); b1=rhs
    do i=2,nn
      coef=mat(i-1,3)/dd(i-1)
      dd(i)=dd(i)-coef*mat(i,1)
      b1(i)=b1(i)-coef*b1(i-1)
    end do
    x(nn)=b1(nn)/dd(nn)
    do i=nn-1,1,-1
      x(i)=(b1(i)-mat(i+1,1)*x(i+1))/dd(i)
    end do
  end subroutine

  subroutine solve_ade(qprev, qup, ck, dk, qnew)
    real(dp), intent(in) :: qprev(n), qup, ck, dk
    real(dp), intent(out) :: qnew(n)
    real(dp) :: mat(n,3), rhs(n), cd, ca, dx, sbc
    integer :: nx
    nx=n-1
    dx=L/real(nx-1,dp)
    cd=dk*dt/(dx*dx)
    ca=ck*dt/dx
    mat=0._dp; rhs=0._dp
    mat(1,2)=1._dp
    mat(2:n-1,2)=2._dp+4._dp*cd
    mat(n,2)=1._dp
    mat(3:n,1)=ca-2._dp*cd
    mat(1:n-2,3)=-ca-2._dp*cd
    mat(n-1,3)=-1._dp
    rhs(1)=qup
    sbc=qprev(n)-qprev(n-1)
    rhs(n)=sbc
    rhs(2:n-1)=2._dp*qprev(2:n-1)
    call tdma(n,mat,rhs,qnew)
  end subroutine

  subroutine kwe_step(profile,vol,qin,qlat,qout)
    real(dp), intent(inout) :: profile(n), vol
    real(dp), intent(in) :: qin, qlat
    real(dp), intent(out) :: qout
    real(dp) :: qbar,y,ck,qnew(n),channel,reduction,vtmp
    qbar=(qin+profile(1)+profile(n-1))/3._dp
    y=flow_depth(abs(qbar)); ck=cel(abs(qbar),y)
    call solve_ade(profile,qin,ck,0._dp,qnew)
    channel=qnew(n-1)
    if (abs(channel)>0._dp) then
      vtmp=max(0._dp,vol)
      reduction=min((vtmp+dt*qin)*0.999_dp/(channel*dt),1._dp)
      qnew(2:n)=qnew(2:n)*reduction
      channel=qnew(n-1)
    end if
    vol=vol+(qin-channel)*dt
    profile=qnew
    qout=channel+qlat
  end subroutine

  subroutine dw_step(profile,vol,qin,qlat,qout)
    real(dp), intent(inout) :: profile(n), vol
    real(dp), intent(in) :: qin, qlat
    real(dp), intent(out) :: qout
    real(dp) :: qbar,y,ck,dk,qnew(n),channel,reduction,vtmp
    qbar=(qin+profile(1)+profile(n-1))/3._dp
    y=flow_depth(abs(qbar)); ck=cel(abs(qbar),y); dk=diffu(abs(qbar),y)
    call solve_ade(profile,qin,ck,dk,qnew)
    channel=qnew(n-1)
    if (abs(channel)>0._dp) then
      vtmp=max(0._dp,vol)
      reduction=min((vtmp+dt*qin)*0.999_dp/(channel*dt),1._dp)
      qnew(2:n)=qnew(2:n)*reduction
      channel=qnew(n-1)
    end if
    vol=vol+(qin-channel)*dt
    profile=qnew
    qout=channel+qlat
  end subroutine

  subroutine mc_step(vol,iprev,oprev,qin,qlat,qout)
    real(dp), intent(inout) :: vol,iprev,oprev
    real(dp), intent(in) :: qin,qlat
    real(dp), intent(out) :: qout
    real(dp) :: qbar0,y,ck0,cn0,dts,oldout,oldin,newin,qbar,bt,ck,x,cn,den,c0,c1,c2,channel,reduction
    real(dp), allocatable :: qouts(:)
    integer :: nsub,k
    qbar0=(iprev+qin+oprev)/3._dp
    channel=0._dp
    if (qbar0>1.e-50_dp) then
      y=flow_depth(abs(qbar0)); ck0=cel(abs(qbar0),y)
      cn0=ck0*dt/L
      if (cn0>1._dp) then
        nsub=ceiling(cn0)
      else
        nsub=1
      end if
      dts=dt/real(nsub,dp)
      allocate(qouts(nsub))
      oldout=oprev; oldin=iprev
      do k=1,nsub
        newin=qin
        qbar=(newin+oldin+oldout)/3._dp
        if (qbar>1.e-50_dp) then
          y=flow_depth(abs(qbar)); bt=btop(y); ck=cel(abs(qbar),y)
          x=0.5_dp*(1._dp-qbar/(bt*S*ck*L))
          cn=ck*dts/L
          den=1._dp-x+0.5_dp*cn
          c0=(-x+0.5_dp*cn)/den
          c1=( x+0.5_dp*cn)/den
          c2=(1._dp-x-0.5_dp*cn)/den
          qouts(k)=max(0._dp,c0*newin+c1*oldin+c2*oldout)
        else
          qouts(k)=0._dp
        end if
        oldout=qouts(k); oldin=newin
      end do
      channel=sum(qouts)/real(nsub,dp)
      if (abs(channel)>0._dp) then
        reduction=min((vol/dt+qin)*0.999_dp/channel,1._dp)
        channel=channel*reduction
      end if
      deallocate(qouts)
    end if
    vol=vol+(qin-channel)*dt
    iprev=qin; oprev=channel
    qout=max(0._dp,channel+qlat)
  end subroutine

end program

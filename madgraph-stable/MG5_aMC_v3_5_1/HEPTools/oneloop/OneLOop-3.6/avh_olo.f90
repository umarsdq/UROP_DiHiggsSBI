!
! Copyright (C) 2015 Andreas van Hameren. 
!
! This file is part of OneLOop-3.6.
!
! OneLOop-3.6 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! OneLOop-3.6 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with OneLOop-3.6.  If not, see <http://www.gnu.org/licenses/>.
!


module avh_olo_version
  implicit none
  private
  public :: olo_version
  logical ,save :: done=.false.
contains
  subroutine olo_version
  if (done) return ;done=.true.
  write(*,'(a72)') '########################################################################'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '#                      You are using OneLOop-3.6                       #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '# for the evaluation of 1-loop scalar 1-, 2-, 3- and 4-point functions #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '# author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>           #'
  write(*,'(a72)') '#   date: 18-02-2015                                                   #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '# Please cite                                                          #'
  write(*,'(a72)') '#    A. van Hameren,                                                   #'
  write(*,'(a72)') '#      Comput.Phys.Commun. 182 (2011) 2427-2438, arXiv:1007.4716       #'
  write(*,'(a72)') '#    A. van Hameren, C.G. Papadopoulos and R. Pittau,                  #'
  write(*,'(a72)') '#      JHEP 0909:106,2009, arXiv:0903.4665                             #'
  write(*,'(a72)') '# in publications with results obtained with the help of this program. #'
  write(*,'(a72)') '#                                                                      #'
  write(*,'(a72)') '########################################################################'
  end subroutine
end module


module avh_olo_units
  implicit none
  integer :: eunit=6
  integer :: wunit=6
  integer :: munit=6
  integer :: punit=0 ! print all
contains
  subroutine set_unit( message ,val )
!***********************************************************************
! message is intended to be one of the following:
! 'printall', 'message' ,'warning' ,'error'
!***********************************************************************
  character(*) ,intent(in) :: message
  integer      ,intent(in) :: val
  if (.false.) then
  elseif (message(1:8).eq.'printall') then ;punit=val
  elseif (message(1:7).eq.'message' ) then ;munit=val
  elseif (message(1:7).eq.'warning' ) then ;wunit=val
  elseif (message(1:5).eq.'error'   ) then ;eunit=val
  else
    eunit=val
    wunit=val
    munit=val
    punit=0
  endif
  end subroutine
end module


module avh_olo_dp_kinds
  integer ,parameter :: kindr2=selected_real_kind(15) 
end module


module avh_olo_dp_arrays
  use avh_olo_units
  use avh_olo_dp_kinds 
  implicit none
  private
  public :: shift1,shift2,shift3,resize,enlarge

! Increase the size of the last dimension by one,
! and move  x(...,n:nsize)  to  x(...,n+1:nsize+1).
  interface shift1 ! for x(:)
    module procedure shift1_r,shift1_i
  end interface
  interface shift2 ! for x(:,:)
    module procedure shift2_r,shift2_i
  end interface
  interface shift3 ! for x(:,:,:)
    module procedure shift3_r,shift3_i
  end interface

! Resize x to the new bounds. Anything that doesn't fit anymore is lost.
  interface resize
    module procedure resize1_r,resize2_r
  end interface

! Resize x to the maximum of the bounds it has and then new bounds.
  interface enlarge
    module procedure enlarge1_r,enlarge2_r
  end interface

contains

  subroutine shift1_r( xx ,nn )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:)
  integer        ,intent(in   ) :: nn
  real(kindr2) &  
    ,allocatable :: tt(:)
  integer ,parameter :: dm=1
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift1_r'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(dm):ub(dm)))
  xx(lb(dm):nn-1) = tt(lb(dm):nn-1)
  xx(nn+1:ub(dm)) = tt(nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift1_i( xx ,nn )
  integer ,allocatable ,intent(inout) :: xx(:)
  integer              ,intent(in   ) :: nn
  integer ,allocatable :: tt(:)
  integer ,parameter :: dm=1
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift1_i'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(dm):ub(dm)))
  xx(lb(dm):nn-1) = tt(lb(dm):nn-1)
  xx(nn+1:ub(dm)) = tt(nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift2_r( xx ,nn )
  real(kindr2) &  
          ,allocatable ,intent(inout) :: xx(:,:)
  integer              ,intent(in   ) :: nn
  real(kindr2) &  
          ,allocatable :: tt(:,:)
  integer ,parameter :: dm=2
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift2_r'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(dm):ub(dm)))
  xx(:,lb(dm):nn-1) = tt(:,lb(dm):nn-1)
  xx(:,nn+1:ub(dm)) = tt(:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift2_i( xx ,nn )
  integer ,allocatable ,intent(inout) :: xx(:,:)
  integer              ,intent(in   ) :: nn
  integer ,allocatable :: tt(:,:)
  integer ,parameter :: dm=2
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift2_i'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(dm):ub(dm)))
  xx(:,lb(dm):nn-1) = tt(:,lb(dm):nn-1)
  xx(:,nn+1:ub(dm)) = tt(:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift3_r( xx ,nn )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:,:,:)
  integer        ,intent(in   ) :: nn
  real(kindr2) &  
    ,allocatable :: tt(:,:,:)
  integer ,parameter :: dm=3
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift3_r'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  xx(:,:,lb(dm):nn-1) = tt(:,:,lb(dm):nn-1)
  xx(:,:,nn+1:ub(dm)) = tt(:,:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift3_i( xx ,nn )
  integer ,allocatable ,intent(inout) :: xx(:,:,:)
  integer              ,intent(in   ) :: nn
  integer ,allocatable :: tt(:,:,:)
  integer ,parameter :: dm=3
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift3_i'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  xx(:,:,lb(dm):nn-1) = tt(:,:,lb(dm):nn-1)
  xx(:,:,nn+1:ub(dm)) = tt(:,:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

 
  subroutine resize1_r( xx ,l1,u1 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:)
  integer        ,intent(in   ) :: l1,u1
  real(kindr2) &  
    ,allocatable :: tt(:)
  integer :: lb(1),ub(1)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1)))
  tt = xx
  deallocate(xx)
  allocate( xx(l1:u1) )
  lb(1)=max(l1,lb(1)) ;ub(1)=min(u1,ub(1))
  xx(lb(1):ub(1)) = tt(lb(1):ub(1))
  deallocate(tt)
  end subroutine 

  subroutine resize2_r( xx ,l1,u1 ,l2,u2 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:,:)
  integer        ,intent(in   ) :: l1,u1,l2,u2
  real(kindr2) &  
    ,allocatable :: tt(:,:)
  integer :: lb(2),ub(2)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1,l2:u2))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(2):ub(2)))
  tt = xx
  deallocate(xx)
  allocate( xx(l1:u1,l2:u2) )
  lb(1)=max(l1,lb(1)) ;ub(1)=min(u1,ub(1))
  lb(2)=max(l2,lb(2)) ;ub(2)=min(u2,ub(2))
  xx(lb(1):ub(1),lb(2):ub(2)) = &
  tt(lb(1):ub(1),lb(2):ub(2))
  deallocate(tt)
  end subroutine 


  subroutine enlarge1_r( xx ,l1,u1 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:)
  integer        ,intent(in   ) :: l1,u1
  real(kindr2) &  
    ,allocatable :: tt(:)
  integer :: lb(1),ub(1)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  if (lb(1).le.l1.and.u1.le.ub(1)) return
  if (lb(1).gt.ub(1)) then
    deallocate( xx )
    allocate( xx(min(l1,lb(1)):max(u1,ub(1))) )
    return
  endif
  allocate(tt(lb(1):ub(1)))
  tt = xx
  deallocate(xx)
  allocate( xx(min(l1,lb(1)):max(u1,ub(1))) )
  xx(lb(1):ub(1)) = tt(lb(1):ub(1))
  deallocate(tt)
  end subroutine 

  subroutine enlarge2_r( xx ,l1,u1 ,l2,u2 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:,:)
  integer        ,intent(in   ) :: l1,u1,l2,u2
  real(kindr2) &  
    ,allocatable :: tt(:,:)
  integer :: lb(2),ub(2)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1,l2:u2))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  if (lb(1).le.l1.and.u1.le.ub(1).and. &
      lb(2).le.l2.and.u2.le.ub(2)      ) return
  if (lb(1).gt.ub(1).or.lb(2).gt.ub(2)) then
    deallocate( xx )
    allocate( xx(min(l1,lb(1)):max(u1,ub(1))  &
                ,min(l2,lb(2)):max(u2,ub(2))) )
    return
  endif
  allocate(tt(lb(1):ub(1),lb(2):ub(2)))
  tt = xx
  deallocate(xx)
  allocate( xx(min(l1,lb(1)):max(u1,ub(1))  &
              ,min(l2,lb(2)):max(u2,ub(2))) )
  xx(lb(1):ub(1),lb(2):ub(2)) = &
  tt(lb(1):ub(1),lb(2):ub(2))
  deallocate(tt)
  end subroutine 

end module


module avh_olo_dp_prec
  use avh_olo_dp_kinds

  implicit none
  public
  private :: IMAG,acmplx_r,acmplx_rr,acmplx_ir,acmplx_ri,acmplx_c

  integer ,save :: prcpar=0
  integer ,save :: ndecim(1)
  real(kindr2) &
          ,save :: epsilo(1),neglig(1)

  real(kindr2) &
    ,save :: RZRO ,RONE ,EPSN ,EPSN2 ,TWOPI ,ONEPI
  complex(kindr2) &
    ,save :: IEPS ,CZRO ,CONE ,IMAG ,PISQo24 ,IPI

  interface acmplx
    module procedure acmplx_r,acmplx_rr,acmplx_ir,acmplx_ri,acmplx_c
  end interface

contains


  subroutine set_precision( newprc )
!***********************************************************************
!***********************************************************************
  use avh_olo_units
  logical ,intent(out) :: newprc
  integer :: ndec                                  
  if (prcpar.eq.1) then                    
    newprc = .false.                             
    return                                       
  endif
  prcpar = 1                                   
  call set_epsn
  newprc = .true.                              
  RZRO=0
  RONE=1
  IMAG=cmplx(0,1,kind=kind(IMAG))
  CZRO=RZRO
  CONE=RONE
  ONEPI=4*atan(RONE)
  TWOPI=2*ONEPI
  PISQo24=CONE*ONEPI*ONEPI/24
  IPI=IMAG*ONEPI
  EPSN2= EPSN*EPSN
  IEPS= EPSN2*IMAG
!
  contains
!
  subroutine set_epsn
  EPSN = epsilon(EPSN)                         
  ndec = -log10(EPSN)                            
  ndecim(prcpar) = ndec                          
  epsilo(prcpar) = EPSN                        
  neglig(prcpar) = EPSN*10**(ndec/7)       
  end subroutine
!
  end subroutine


  function adble(xx) result(rslt)
!***********************************************************************
! Turn real(kindr2) into kind(1d0)
!***********************************************************************
  real(kindr2) ,intent(in) :: xx
  real(kind(1d0)) :: rslt
  rslt = real(xx,kind=kind(rslt))
  end function

  function convert(xx) result(rslt)
!***********************************************************************
! Turn kind(1d0) into real(kindr2)
!***********************************************************************
  real(kind(1d0)) ,intent(in) :: xx
  real(kindr2) :: rslt
  rslt = real(xx,kind=kind(rslt))
  end function

  function areal(zz) result(rslt)
!***********************************************************************
! Get real part of a complex
!***********************************************************************
  complex(kindr2) &
    ,intent(in) :: zz
  real(kindr2) &
    :: rslt
  rslt = zz
  end function

  function acmplx_r(xx) result(rslt)
!***********************************************************************
! Turn a real into a complex
!***********************************************************************
  real(kindr2) &
    ,intent(in) :: xx
  complex(kindr2) &
    :: rslt
  rslt = xx
  end function
  
  function acmplx_rr(xx,yy) result(rslt)
!***********************************************************************
! Turn two reals into one complex
!***********************************************************************
  real(kindr2) &
    ,intent(in) :: xx,yy
  complex(kindr2) &
    :: rslt
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
  function acmplx_ri(xx,yy) result(rslt)
!***********************************************************************
! Turn a real and an integer into one complex
!***********************************************************************
  real(kindr2) &
           ,intent(in) :: xx
  integer  ,intent(in) :: yy
  complex(kindr2) &
    :: rslt
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
  function acmplx_ir(xx,yy) result(rslt)
!***********************************************************************
! Turn an integer and a real into one complex
!***********************************************************************
  integer ,intent(in) :: xx
  real(kindr2) &
          ,intent(in) :: yy
  complex(kindr2) &
    :: rslt
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
  function acmplx_c(zz) result(rslt)
!***********************************************************************
! Replaces the real part of zz by its absolute value
!***********************************************************************
  complex(kindr2) &
    ,intent(in) :: zz
  complex(kindr2) &
    :: rslt
  real(kindr2) &
    :: xx,yy
  xx = zz
  xx = abs(xx)
  yy = aimag(zz)
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
end module


module avh_olo_dp_print
  use avh_olo_dp_prec
  implicit none
  private
  public :: myprint

  integer ,parameter :: novh=10 !maximally 6 decimals for exponent
  integer ,parameter :: nxtr=4  !extra decimals

  interface myprint
    module procedure printr,printc,printi
  end interface

contains

  function printc( zz ,ndec ) result(rslt)
  complex(kindr2) &   
    ,intent(in) :: zz
  integer,optional,intent(in) :: ndec
  character((ndecim(prcpar)+nxtr+novh)*2+3) :: rslt
  if (present(ndec)) then
    rslt = '('//trim(printr(areal(zz),ndec)) &
         //','//trim(printr(aimag(zz),ndec)) &
         //')'
  else
    rslt = '('//trim(printr(areal(zz))) &
         //','//trim(printr(aimag(zz))) &
         //')'
  endif
  rslt = adjustl(rslt)
  end function

  function printr( xx_in ,ndec_in ) result(rslt)
  real(kindr2) &  
                  ,intent(in) :: xx_in
  integer,optional,intent(in) :: ndec_in
  character(ndecim(prcpar)+nxtr+novh  ) :: rslt
  character(ndecim(prcpar)+nxtr+novh+1) :: cc
  character(10) :: aa,bb
  integer :: ndec
  real(kindr2) :: xx     
  xx = xx_in
  if (present(ndec_in)) then ;ndec=ndec_in
                        else ;ndec=ndecim(prcpar)+nxtr
  endif
  write(aa,'(i10)') min(len(cc),ndec+novh+1) ;aa=adjustl(aa)
  write(bb,'(i10)') min(len(cc),ndec       ) ;bb=adjustl(bb)
  aa = '(e'//trim(aa)//'.'//trim(bb)//')'
  write(cc,aa) xx  ;cc=adjustl(cc)
  if (cc(1:2).eq.'-0') then ;rslt = '-'//cc(3:len(cc))
  else                      ;rslt = ' '//cc(2:len(cc))
  endif
  end function

  function printi( ii ) result(rslt)
  integer ,intent(in) :: ii
  character(ndecim(prcpar)) :: rslt
  character(ndecim(prcpar)) :: cc
  character(10) :: aa
  write(aa,'(i10)') ndecim(prcpar) ;aa=adjustl(aa)
  aa = '(i'//trim(aa)//')'
  write(cc,aa) ii ;cc=adjustl(cc)
  if (cc(1:1).ne.'-') then ;rslt=' '//cc
  else                     ;rslt=cc 
  endif
  end function

end module


module avh_olo_dp_auxfun
  use avh_olo_units
  use avh_olo_dp_prec

  implicit none
  private
  public :: mysqrt,eta5,eta3,eta2,sgnIm,sgnRe,kallen
  public :: solabc,rfun,rfun0,solabc_rcc

  interface mysqrt
    module procedure mysqrt_c,mysqrt_cr,mysqrt_ci
  end interface

  interface eta5
    module procedure eta5_0
  end interface
  interface eta3
    module procedure eta3_r,eta3_0
  end interface
  interface eta2
    module procedure eta2_r,eta2_0
  end interface

  interface sgnIm
    module procedure sgnIm_c,sgnIm_ci
  end interface
  interface sgnRe
    module procedure sgnRe_c,sgnRe_r,sgnRe_ri
  end interface

contains


  function mysqrt_c(xx) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! negative imaginary.
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt ,zz
  real(kindr2) &  
    :: xim,xre
  xim = aimag(xx)
  if (xim.eq.RZRO) then
    xre = areal(xx)
    if (xre.ge.RZRO) then
      zz = acmplx(sqrt(xre),0)
    else
      zz = acmplx(0,-sqrt(-xre))
    endif
  else
    zz = sqrt(xx)
  endif
  rslt = zz
  end function

  function mysqrt_cr(xx,sgn) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! imaginary and has the same sign as  sgn .
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  real(kindr2) &  
    ,intent(in) :: sgn
  complex(kindr2) &   
    :: rslt ,zz
  real(kindr2) &  
    :: xim,xre
  xim = aimag(xx)
  if (xim.eq.RZRO) then
    xre = areal(xx)
    if (xre.ge.RZRO) then
      zz = acmplx(sqrt(xre),0)
    else
      zz = acmplx(0,sign(sqrt(-xre),sgn))
    endif
  else
    zz = sqrt(xx)
  endif
  rslt = zz
  end function

  function mysqrt_ci(xx,sgn) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! imaginary and has the same sign as  sgn .
!*******************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: sgn
  complex(kindr2) &   
    :: rslt ,zz
  real(kindr2) &  
    :: xim,xre,hh
  xim = aimag(xx)
  if (xim.eq.RZRO) then
    xre = areal(xx)
    if (xre.ge.RZRO) then
      zz = acmplx(sqrt(xre),0)
    else
      hh = sgn
      zz = acmplx(0,sign(sqrt(-xre),hh))
    endif
  else
    zz = sqrt(xx)
  endif
  rslt = zz
  end function


  subroutine solabc( x1,x2 ,dd ,aa,bb,cc ,imode )
!*******************************************************************
! Returns the solutions  x1,x2  to the equation  aa*x^2+bb*x+cc=0
! Also returns  dd = aa*(x1-x2)
! If  imode=/=0  it uses  dd  as input as value of  sqrt(b^2-4*a*c)
!*******************************************************************
  complex(kindr2) &   
    ,intent(out)   :: x1,x2
  complex(kindr2) &   
    ,intent(inout) :: dd
  complex(kindr2) &   
    ,intent(in) :: aa,bb,cc
  integer         ,intent(in) :: imode
  complex(kindr2) &   
    :: qq,hh
  real(kindr2) &  
    :: r1,r2

  if (aa.eq.CZRO) then
    if (bb.eq.CZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop solabc: ' &
        ,'no solutions, returning 0'
      x1 = 0
      x2 = 0
      dd = 0
    else
      x1 = -cc/bb
      x2 = x1
      dd = bb
    endif
  elseif (cc.eq.CZRO) then
    dd = -bb
    x1 = dd/aa
    x2 = 0
  else
    if (imode.eq.0) dd = sqrt(bb*bb - 4*aa*cc)
    qq = -bb+dd
    hh = -bb-dd
    r1 = abs(qq)
    r2 = abs(hh)
    if (r1.ge.r2) then
      x1 = qq/(2*aa)
      x2 = (2*cc)/qq
    else
      qq = hh
      x2 = qq/(2*aa)
      x1 = (2*cc)/qq
    endif
  endif
  end subroutine


  subroutine solabc_rcc( x1,x2 ,aa,bb,cc )
!*******************************************************************
! Tested
!*******************************************************************
  intent(out) :: x1,x2
  intent(in ) :: aa,bb,cc
  complex(kindr2) &   
    :: x1,x2,bb,cc ,t1,t2
  real(kindr2) &  
    :: aa,xx,yy,pp,qq,uu,vv,pq1,pq2,uv1,uv2,dd,xd1,xd2,yd1,yd2 &
      ,gg,hh,rx1,rx2,ix1,ix2
  if (aa.eq.RZRO) then
    if (bb.eq.CZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop solabc: ' &
        ,'no solutions, returning 0'
      x1 = 0
      x2 = 0
    else
      x1 = -cc/bb
      x2 = x1
    endif
  elseif (cc.eq.CZRO) then
    x1 = -bb/aa
    x2 = 0
  else
    t1 = cc/aa          ;xx= areal(t1) ;yy= aimag(t1)
    t2 = bb/(aa*2)      ;pp=-areal(t2) ;uu=-aimag(t2)
    t2 = sqrt(t2*t2-t1) ;qq= areal(t2) ;vv= aimag(t2)
    pq1=pp+qq ;uv1=uu+vv
    pq2=pp-qq ;uv2=uu-vv
    dd=pq1*pq1+uv1*uv1 ;xd1=xx/dd ;yd1=yy/dd
    dd=pq2*pq2+uv2*uv2 ;xd2=xx/dd ;yd2=yy/dd
    if (abs(pq1).gt.abs(pq2)) then
      rx1 = pq1
      gg=xd1*pq1 ;hh=yd1*uv1
      rx2 = gg+hh
      if (abs(rx2).lt.neglig(prcpar)*max(abs(gg),abs(hh))) rx2 = 0
    elseif (abs(pq2).gt.abs(pq1)) then
      rx2 = pq2
      gg=xd2*pq2 ;hh=yd2*uv2
      rx1 = gg+hh
      if (abs(rx1).lt.neglig(prcpar)*max(abs(gg),abs(hh))) rx1 = 0
    else
      rx1 = pq1
      rx2 = pq2
    endif
    if (abs(uv1).gt.abs(uv2)) then
      ix1 = uv1
      gg=yd1*pq1 ;hh=xd1*uv1
      ix2 = gg-hh
      if (abs(ix2).lt.neglig(prcpar)*max(abs(gg),abs(hh))) ix2 = 0
    elseif (abs(uv2).gt.abs(uv1)) then
      ix2 = uv2
      gg=yd2*pq2 ;hh=xd2*uv2
      ix1 = gg-hh
      if (abs(ix1).lt.neglig(prcpar)*max(abs(gg),abs(hh))) ix1 = 0
    else
      ix1 = uv1
      ix2 = uv2
    endif
    x1 = acmplx(rx1,ix1)
    x2 = acmplx(rx2,ix2)
  endif
  end subroutine


  subroutine rfun(rr,dd ,qq)
!*******************************************************************
! Returns  rr  such that  qq = rr + 1/rr  and  Im(rr)  has the same
! sign as  Im(qq) .
! If  Im(qq)  is zero, then  Im(rr)  is negative or zero.
! If  Im(rr)  is zero, then  |rr| > 1/|rr| .
! Also returns  dd = rr - 1/rr .
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rr,dd
  complex(kindr2) &   
    ,intent(in)  :: qq
  complex(kindr2) &   
    :: r2
  real(kindr2) &  
    :: aa,bb
  integer :: ir,ik
  dd = sqrt(qq*qq-4)
  rr = qq+dd
  r2 = qq-dd
  aa = abs(rr)
  bb = abs(r2)
  if (bb.gt.aa) then
    rr = r2
    dd = -dd
  endif
  aa = aimag(qq)
  bb = aimag(rr)
  if (aa.eq.RZRO) then
    if (bb.le.RZRO) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  else
    ik = sgnRe(aa)
    ir = sgnRe(bb)
    if (ir.eq.ik) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  endif
  end subroutine

  subroutine rfun0(rr ,dd,qq)
!*******************************************************************
! Like rfun, but now  dd  is input, which may get a minus sign
!*******************************************************************
  complex(kindr2) &   
    ,intent(out)   :: rr
  complex(kindr2) &   
    ,intent(inout) :: dd
  complex(kindr2) &   
    ,intent(in)  :: qq
  complex(kindr2) &   
    :: r2
  real(kindr2) &  
    :: aa,bb
  integer :: ir,ik
  rr = qq+dd
  r2 = qq-dd
  aa = abs(rr)
  bb = abs(r2)
  if (bb.gt.aa) then
    rr = r2
    dd = -dd
  endif
  aa = aimag(qq)
  bb = aimag(rr)
  if (aa.eq.RZRO) then
    if (bb.le.RZRO) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  else
    ik = sgnRe(aa)
    ir = sgnRe(bb)
    if (ir.eq.ik) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  endif
  end subroutine


  function eta3_r( aa,sa ,bb,sb ,cc,sc ) result(rslt)
!*******************************************************************
! 2*pi*imag times the result of
!     theta(-Im(a))*theta(-Im(b))*theta( Im(c))
!   - theta( Im(a))*theta( Im(b))*theta(-Im(c))
! where a,b,c are interpreted as a+i|eps|sa, b+i|eps|sb, c+i|eps|sc
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb,cc
  real(kindr2) &  
    ,intent(in) :: sa,sb,sc
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: ima,imb,imc
  ima = aimag(aa)
  imb = aimag(bb)
  imc = aimag(cc)
  if (ima.eq.RZRO) ima = sa
  if (imb.eq.RZRO) imb = sb
  if (imc.eq.RZRO) imc = sc
  ima = sgnRe(ima)
  imb = sgnRe(imb)
  imc = sgnRe(imc)
  if (ima.eq.imb.and.ima.ne.imc) then
    rslt = acmplx(0,imc*TWOPI)
  else
    rslt = 0
  endif
  end function

  function eta3_0( aa ,bb ,cc ) result(rslt)
!*******************************************************************
! 2*pi*imag times the result of
!     theta(-Im(a))*theta(-Im(b))*theta( Im(c))
!   - theta( Im(a))*theta( Im(b))*theta(-Im(c))
! where a,b,c are interpreted as a+i|eps|sa, b+i|eps|sb, c+i|eps|sc
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb,cc
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: ima,imb,imc
  ima = sgnIm(aa)
  imb = sgnIm(bb)
  imc = sgnIm(cc)
  if (ima.eq.imb.and.ima.ne.imc) then
    rslt = acmplx(0,imc*TWOPI)
  else
    rslt = 0
  endif
  end function

  function eta5_0( aa ,b1,c1 ,b2,c2 ) result(rslt)
!*******************************************************************
! eta3(aa,b1,c1) - eta3(aa,b2,c2)
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,b1,c1 ,b2,c2
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: imaa,imb1,imc1,imb2,imc2
  imaa = sgnIm(aa)
  imb1 = sgnIm(b1)
  imb2 = sgnIm(b2)
  imc1 = sgnIm(c1)
  imc2 = sgnIm(c2)
  if (imaa.eq.imb1) then
    if (imaa.eq.imb2) then
      if (imc1.eq.imc2) then
        rslt = 0
      elseif (imaa.ne.imc1) then
        rslt = acmplx(0, imc1*TWOPI)
      else
        rslt = acmplx(0,-imc2*TWOPI)
      endif
    elseif (imaa.ne.imc1) then
      rslt = acmplx(0, imc1*TWOPI)
    else
      rslt = 0
    endif
  elseif (imaa.eq.imb2.and.imaa.ne.imc2) then
    rslt = acmplx(0,-imc2*TWOPI)
  else
    rslt = 0
  endif
  end function

  function eta2_r( aa,sa ,bb,sb ) result(rslt)
!*******************************************************************
! The same as  eta3, but with  c=a*b, so that
!   eta(a,b) = log(a*b) - log(a) - log(b)
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb
  real(kindr2) &  
    ,intent(in) :: sa,sb
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rea,reb,ima,imb,imab
  rea = areal(aa)  ;ima = aimag(aa)
  reb = areal(bb)  ;imb = aimag(bb)
  imab = rea*imb + reb*ima
  if (ima .eq.RZRO) ima = sa
  if (imb .eq.RZRO) imb = sb
  if (imab.eq.RZRO) imab = sign(rea,sb) + sign(reb,sa)
  ima  = sgnRe(ima)
  imb  = sgnRe(imb)
  imab = sgnRe(imab)
  if (ima.eq.imb.and.ima.ne.imab) then
    rslt = acmplx(0,imab*TWOPI)
  else
    rslt = 0
  endif
  end function
 
  function eta2_0( aa ,bb ) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rea,reb,ima,imb,imab
  rea = areal(aa)  ;ima = aimag(aa)
  reb = areal(bb)  ;imb = aimag(bb)
  rea = rea*imb
  reb = reb*ima
  imab = rea+reb
  ima  = sgnRe(ima)
  imb  = sgnRe(imb)
  imab = sgnRe(imab)
  if (ima.eq.imb.and.ima.ne.imab) then
    rslt = acmplx(0,imab*TWOPI)
  else
    rslt = 0
  endif
  end function 


  function kallen( p1,p2,p3 ) result(rslt)
!*******************************************************************
!  p1^2 + p2^2 + p3^2 - 2*p1*p2 - 2*p2*p3 - 2*p3*p1
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3
  complex(kindr2) &   
    :: rslt ,y1,y2,y3
  real(kindr2) &  
    :: b1,b2,b3
  y1=p2*p3 ;b1=areal(y1)
  y2=p3*p1 ;b2=areal(y2)
  y3=p1*p2 ;b3=areal(y3)
      if (b1.le.RZRO) then  ;rslt = (p1-p2-p3)**2 - 4*y1
  elseif (b2.le.RZRO) then  ;rslt = (p2-p3-p1)**2 - 4*y2
  elseif (b3.le.RZRO) then  ;rslt = (p3-p1-p2)**2 - 4*y3
  elseif (b1.le.b2.and.b1.le.b3) then  ;rslt = (p1-p2-p3)**2 - 4*y1
  elseif (b2.le.b3.and.b2.le.b1) then  ;rslt = (p2-p3-p1)**2 - 4*y2
                                 else  ;rslt = (p3-p1-p2)**2 - 4*y3
  endif
  end function


  function sgnIm_c(zz) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: zz
  integer :: rslt
  real(kindr2) &  
    :: imz
  imz = aimag(zz)
  if (imz.ge.RZRO) then ;rslt= 1
                   else ;rslt=-1
  endif
  end function

  function sgnIm_ci(zz,ii) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
          ,intent(in) :: zz
  integer ,intent(in) :: ii
  integer :: rslt
  real(kindr2) &  
    :: imz
  imz = aimag(zz)
  if     (imz.gt.RZRO) then ;rslt= 1
  elseif (imz.lt.RZRO) then ;rslt=-1
                       else ;rslt= sign(1,ii)
  endif
  end function

  function sgnRe_c(zz) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: zz
  integer :: rslt
  real(kindr2) &  
    :: rez
  rez = zz
  if (rez.ge.RZRO) then ;rslt= 1
                   else ;rslt=-1
  endif
  end function

  function sgnRe_r(rez) result(rslt)
!*******************************************************************
!*******************************************************************
  real(kindr2) &  
    ,intent(in) :: rez
  integer :: rslt
  if (rez.ge.RZRO) then ;rslt= 1
                   else ;rslt=-1
  endif
  end function

  function sgnRe_ri(rez,ii) result(rslt)
!*******************************************************************
!*******************************************************************
  real(kindr2) &  
          ,intent(in) :: rez
  integer ,intent(in) :: ii
  integer :: rslt
  if     (rez.gt.RZRO) then ;rslt= 1
  elseif (rez.lt.RZRO) then ;rslt=-1
                       else ;rslt=sign(1,ii)
  endif
  end function

end module


module avh_olo_dp_olog
!***********************************************************************
! Provides the functions
!   olog(x,n) = log(x) + n*pi*imag  
!   olog1(x,n) = olog(x,n)/(x-1)
!   olog2(x,n) = ( olog1(x,n) - 1 )/(x-1)
!   olog3(x,n) = ( olog2(x,n) + 1/2 )/(x-1)
! In the vicinity of x=1,n=0, the logarithm of complex argument is
! evaluated with a series expansion.
!***********************************************************************
  use avh_olo_units
  use avh_olo_dp_prec
  use avh_olo_dp_print
  use avh_olo_dp_auxfun
  implicit none
  private
  public :: update_olog,olog,olog1,olog2,olog3

  real(kindr2) &  
         ,allocatable,save :: thrs(:,:)
  integer,allocatable,save :: ntrm(:,:)
  integer,parameter :: nStp=6

  interface olog
    module procedure log_c,log_r
  end interface
  interface olog1
    module procedure log1_c,log1_r
  end interface
  interface olog2
    module procedure log2_c,log2_r
  end interface
  interface olog3
    module procedure log3_c,log3_r
  end interface

contains

  subroutine update_olog
!***********************************************************************
!***********************************************************************
  use avh_olo_dp_arrays
  real(kindr2) &  
    :: tt
  integer :: nn,mm,ii,jj
!  real(kind(1d0)) :: xx(6) !DEBUG
  if (allocated(thrs)) then
    call shift2( thrs ,prcpar )
    call shift2( ntrm ,prcpar )
  else
    allocate(thrs(1:nStp,1:1))
    allocate(ntrm(1:nStp,1:1))
    if (prcpar.ne.1) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop update_olog'
      stop
    endif
  endif
  if (prcpar.gt.1) then ;nn=ntrm(nStp,prcpar-1)-1
                   else ;nn=1
  endif
  do
    nn = nn+1
    mm = 2*nn-1
    tt = 1
    tt = (EPSN*mm)**(tt/(mm-1))
    tt = 2*tt/(1-tt)
! expansion from x=1+d with |d|=1/1000
    if (1000*tt.gt.RONE) exit
  enddo
  ntrm(nStp,prcpar) = nn
  thrs(nStp,prcpar) = tt
  nn = max(1,nint(nn*1d0/nStp))
  do ii=nStp-1,1,-1
    ntrm(ii,prcpar) = ntrm(ii+1,prcpar)-nn
    if (ntrm(ii,prcpar).le.1) then
      do jj=1,ii
        ntrm(jj,prcpar) = ntrm(ii,prcpar)
        thrs(jj,prcpar) = 0 
      enddo
      exit
    endif
    mm = 2*ntrm(ii,prcpar)-1
    tt = 1
    tt = (EPSN*mm)**(tt/(mm-1))
    thrs(ii,prcpar) = 2*tt/(1-tt)
  enddo
!  do ii=lbound(thrs,2),ubound(thrs,2) !DEBUG
!    do jj=1,nStp                      !DEBUG
!      xx(jj) = thrs(jj,ii)            !DEBUG
!    enddo                             !DEBUG
!    write(*,'(99e10.3)') xx(:)        !DEBUG
!    write(*,'(99i10)'  ) ntrm(:,ii)   !DEBUG
!  enddo                               !DEBUG
  end subroutine


  function log_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,iyy
!
  rex = areal(xx)
  imx = aimag(xx)
  iyy = iph
!
  if (abs(imx).le.EPSN*abs(rex)) then
    if (rex.ge.RZRO) then
      rslt = log_r( rex, iyy )
    else
      rslt = log_r(-rex, iyy+sgnRe(imx) )
    endif
    return
  endif
!
  if (mod(iyy,2).eq.0) then
    yy = acmplx(rex,imx)
  else
    yy = acmplx(-rex,-imx)
    iyy = iyy+sgnRe(imx)
  endif
!
  if (iyy.ne.0) then
    rslt = log(yy) + IPI*iyy
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = log(yy)
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = 2 * ( z + z^3/3 + z^5/5 + z^7/7 + ... )  
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,1,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = zz*rslt
  end function


  function log_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rr
  integer :: jj
!
  if (xx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log_r: ' &
       ,'xx =',trim(myprint(xx)),', returning 0'
    rslt = 0
    return
  elseif (xx.gt.RZRO) then ;rr= xx ;jj= iph
                      else ;rr=-xx ;jj= iph+1 ! log(-1)=i*pi
  endif
!
  rslt = log(rr) + IPI*jj
  end function


  function log1_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,jj
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (abs(imx).le.EPSN*abs(rex)) then
    if (rex.ge.RZRO) then
      rslt = log1_r( rex, iph )
    else
      rslt = log1_r(-rex, iph+sgnRe(imx) )
    endif
    return
  endif
!
  if (mod(iph,2).eq.0) then ;yy= xx ;jj=iph
                       else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  if (jj.ne.0) then
    rslt = ( log(yy) + IPI*jj )/(yy-1)
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = log(yy)/zz
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = 2/(y+1) * ( 1 + z^2/3 + z^4/5 + z^6/7 + ... )  
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,1,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = rslt/(yy+1)
  end function


  function log1_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rr,yy
  integer :: jj
!  include 'avh_olo_dp_real.h90'
!    :: aa,zz,z2
!  integer :: nn,ii
!
  if (xx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log1_r: ' &
       ,'xx =',trim(myprint(xx)),', returning 0'
    rslt = 0
    return
  elseif (xx.gt.RZRO) then ;rr= xx ;jj=iph
                      else ;rr=-xx ;jj=iph+1 ! log(-1)=i*pi
  endif
!
  yy=rr ;if (mod(jj,2).ne.0) yy=-rr
!
  if (abs(yy-1).le.10*EPSN) then
    if (jj.ne.0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log1_r: ' &
        ,'rr,jj =',trim(myprint(rr)),jj,', putting jj to 0'
    endif
    rslt = 1 - (yy-1)/2
    return
  endif
!
  rslt = ( log(rr) + IPI*jj )/(yy-1)
  end function


  function log2_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  rslt = log2_c(xx*CONE,iph)
  end function


  function log2_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,jj
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (rex.eq.RZRO.and.imx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log2_c: ' &
       ,'xx = 0, returning 0'
    rslt = 0
    return
  endif
!
  if (mod(iph,2).eq.0) then ;yy= xx ;jj=iph
                       else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  if (jj.ne.0) then
    rslt = ( olog1(yy,jj) - 1 )/(yy-1)
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = (log(yy)/zz-1)/zz
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = -1/(y+1) + 2/(y+1)^2 * ( z/3 + z^3/5 + z^5/7 + ... )  
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,3,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = ( -1 + zz*rslt/(yy+1) )/(yy+1)
  end function


  function log3_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  rslt = log3_c(xx*CONE,iph)
  end function


  function log3_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2,HLF
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,jj
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (rex.eq.RZRO.and.imx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log2_c: ' &
       ,'xx = 0, returning 0'
    rslt = 0
    return
  endif
!
  HLF = CONE/2
!
  if (mod(iph,2).eq.0) then ;yy= xx ;jj=iph
                       else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  if (jj.ne.0) then
    rslt = ( olog2(xx,jj) + HLF )/(yy-1)
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = ((log(yy)/zz-1)/zz+HLF)/zz
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = 1/(2*(y+1)) + 2/(y+1)^3 * ( 1/3 + z^2/5 + z^4/7 + ... )
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,3,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = ( HLF + rslt/(yy+1)**2 )/(yy+1)
  end function

end module




module avh_olo_dp_dilog
!***********************************************************************
!                     /1    ln(1-zz*t)
!   dilog(xx,iph) = - |  dt ---------- 
!                     /0        t
! with  zz = 1 - xx*exp(imag*pi*iph)  [pi, NOT 2*pi]
!
!   dilog(x1,i1,x2,i2) = ( dilog(x1,i1)-dilog(x2,i2) )/( x1-x2 )
!
! Arguments xx,x1,x2, may be all real or all complex,
! arguments iph,i1,i2 must be all integer.
!***********************************************************************
  use avh_olo_units
  use avh_olo_dp_prec
  use avh_olo_dp_print
  use avh_olo_dp_auxfun
  use avh_olo_dp_arrays
  implicit none
  private
  public :: update_dilog,dilog

  real(kindr2) &  
         ,allocatable,save :: coeff(:)
  real(kindr2) &  
         ,allocatable,save :: thrs(:,:)
  integer,allocatable,save :: ntrm(:,:)
  integer,parameter :: nStp=6

  real(kindr2) &  
         ,allocatable :: bern(:),fact(:)

  interface dilog
    module procedure dilog_c,dilog_r,dilog2_c,dilog2_r
  end interface

contains

  subroutine update_dilog
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
    :: tt
  integer :: nn,ii,jj
  logical :: highestSoFar
!  real(kind(1d0)) :: xx(6) !DEBUG
!
  if (allocated(thrs)) then
    call shift2( thrs ,prcpar )
    call shift2( ntrm ,prcpar )
  else
    allocate(thrs(1:nStp,1:1))
    allocate(ntrm(1:nStp,1:1))
    if (prcpar.ne.1) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop update_dilog'
      stop
    endif
  endif
!
  highestSoFar = prcpar.eq.ubound(ntrm,2)
  if (highestSoFar) then
    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(0:-1)) ! allocate at size=0
  endif
!
  if (prcpar.gt.1) then ;nn=ntrm(nStp,prcpar-1)-1
                   else ;nn=2
  endif
!
  do
    nn = nn+1
    if (nn.gt.ubound(coeff,1)) call update_coeff( 2*nn )
    tt = 1
    tt = (EPSN/abs(coeff(nn)))**(tt/(2*nn))
! expansion parameter is smaller than 1.05
    if (100*tt.gt.105*RONE) exit
  enddo
!
  if (highestSoFar) call resize( coeff ,0,nn )
!
  ntrm(nStp,prcpar) = nn
  thrs(nStp,prcpar) = tt
  nn = max(1,nint(nn*1d0/nStp))
  do ii=nStp-1,1,-1
    ntrm(ii,prcpar) = ntrm(ii+1,prcpar)-nn
    if (ntrm(ii,prcpar).le.2) then
      do jj=1,ii
        ntrm(jj,prcpar) = max(2,ntrm(ii,prcpar))
        thrs(jj,prcpar) = 0 
      enddo
      exit
    endif
    jj = ntrm(ii,prcpar)
    tt = 1
    tt = (EPSN/abs(coeff(jj)))**(tt/(2*jj))
    thrs(ii,prcpar) = tt
  enddo
!
  if (allocated(bern)) deallocate(bern)
  if (allocated(fact)) deallocate(fact)
!
!  do ii=lbound(thrs,2),ubound(thrs,2) !DEBUG
!    do jj=1,nStp                      !DEBUG
!      xx(jj) = thrs(jj,ii)            !DEBUG
!    enddo                             !DEBUG
!    write(*,'(99e10.3)') xx(:)        !DEBUG
!    write(*,'(99i10)'  ) ntrm(:,ii)   !DEBUG
!  enddo                               !DEBUG
  end subroutine


  subroutine update_coeff( ncf )
!*******************************************************************
!   coeff(0)=-1/4
!   coeff(n)=bern(2*n)/(2*n+1)
!    bern(n)=bernoulli(n)/n!
!    fact(n)=n!
! DO NOT SKIP THE ODD bern IN THE RECURSIVE LOOP
! DO NOT PUT THE ODD bern TO ZERO
!*******************************************************************
  integer ,intent(in) :: ncf
  integer :: ii,jj,nbern,nold
!
  if (allocated(bern)) then ;nold=ubound(bern,1)
                       else ;nold=0
  endif
!
  nbern = 2*ncf
!
  call enlarge( bern  ,1,nbern   )
  call enlarge( fact  ,0,nbern+1 )
  call enlarge( coeff ,0,ncf     )
!
  fact(0) = 1
  do ii=nold+1,nbern+1
    fact(ii) = fact(ii-1)*ii
  enddo
!
  do ii=nold+1,nbern
    bern(ii) = -1/fact(ii+1)
    do jj=1,ii-1
      bern(ii) = bern(ii) - bern(jj)/fact(ii+1-jj)
    enddo
  enddo
!
  coeff(0) = 1
  coeff(0) =-coeff(0)/4
  do ii=nold+2,nbern,2
    coeff(ii/2) = bern(ii)/(ii+1)
  enddo
!
  end subroutine


  function dilog_c(xx,iph) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,lyy,loy,zz,z2
  real(kindr2) &  
    :: rex,imx,az
  integer :: ii,jj,ntwo,odd,nn
  logical :: r_gt_1 , y_lt_h
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (abs(imx).le.EPSN*abs(rex)) then
    if (rex.ge.RZRO) then
      rslt = dilog_r( rex, iph )
    else
      rslt = dilog_r(-rex, iph+sgnRe(imx) )
    endif
    return
  endif
!
  if (rex.gt.RZRO) then ;yy= xx ;jj=iph
                   else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  odd = mod(jj,2)
  ntwo = jj-odd
! 
  r_gt_1 = (rex*rex+imx*imx.gt.RONE)
  lyy = log(yy)
  if (odd.ne.0) yy = -yy
!
  if (r_gt_1) then
    yy   = 1/yy
    lyy  =-lyy
    ntwo =-ntwo
    odd  =-odd
  endif
  loy = log(1-yy)
!
  y_lt_h = (2*areal(yy).lt.RONE)
  if (y_lt_h) then ;zz=-loy
              else ;zz=-lyy
  endif
!
  az = abs(zz)
! if (az.gt.thrs(6,prcpar)) ERROR az to big 
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
  z2 = zz*zz
  rslt = coeff(nn)
  do ii=nn,2,-1
    rslt = coeff(ii-1) + z2*rslt
  enddo
  rslt = zz*( 1 + zz*( coeff(0) + zz*rslt ) )
!
  if (y_lt_h) then
    rslt = 4*PISQo24 - rslt - loy*(lyy+IPI*(ntwo+odd))
  else
    rslt = rslt - loy*IPI*ntwo
  endif
!
  if (r_gt_1) rslt = -rslt - (lyy+IPI*(ntwo+odd))**2/2
  end function



  function dilog_r(xx,iph) result(rslt)
!*******************************************************************
!*******************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: yy,lyy,loy,zz,z2,liox,az
  integer :: jj,ii,ntwo,odd,nn
  logical :: r_gt_1 , y_lt_h
!
  if (xx.eq.RZRO) then
    rslt = 4*PISQo24
    return
  elseif (xx.gt.RZRO) then ;yy= xx ;jj=iph
                      else ;yy=-xx ;jj=iph+1 ! log(-1)=i*pi
  endif
!
  odd = mod(jj,2)
  ntwo = jj-odd
! 
  if (yy.eq.RONE.and.odd.eq.0) then
    if (ntwo.ne.0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog_r: ' &
        ,'|x|,iph = ',trim(myprint(yy)),',',jj,', returning 0'
    endif
    rslt = 0
    return
  endif
!
  r_gt_1 = (yy.gt.RONE)
  lyy = log(yy)
  if (odd.ne.0) yy = -yy
!
  if (r_gt_1) then
    yy   = 1/yy
    lyy  =-lyy
    ntwo =-ntwo
    odd  =-odd
  endif
  loy = log(1-yy) ! log(1-yy) is always real
!
  y_lt_h = (2*yy.lt.RONE)
  if (y_lt_h) then
    zz = -loy ! log(1-yy) is real
  else
    zz = -lyy ! yy>0.5 => log(yy) is real
  endif
!
  az = abs(zz)
! if (az.gt.thrs(6,prcpar)) ERROR az to big 
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
  z2 = zz*zz
  liox = coeff(nn)
  do ii=nn,2,-1
    liox = coeff(ii-1) + z2*liox
  enddo
  liox = zz*( 1 + zz*( coeff(0) + zz*liox ) )
!
  rslt = acmplx(liox)
!
  if (y_lt_h) then
    rslt = 4*PISQo24 - rslt - acmplx(loy*lyy,loy*ONEPI*(ntwo+odd))
  else
    rslt = rslt + acmplx( 0 ,-loy*ONEPI*ntwo )
  endif
!
  if (r_gt_1) rslt = -rslt - acmplx(lyy,ONEPI*(ntwo+odd))**2/2
  end function


  function dilog2_c( x1,i1 ,x2,i2 ) result(rslt)
!*******************************************************************
!*******************************************************************
  use avh_olo_dp_olog
  complex(kindr2) &   
          ,intent(in) :: x1,x2
  integer ,intent(in) :: i1,i2
  complex(kindr2) &   
    :: rslt ,y1,y2 ,ff,gg,logr1,logr2,logo1,logo2,r1,r2,rr
  real(kindr2) &  
    :: eps ,re1,im1,re2,im2,a1,a2,aa,ao1,ao2
  integer :: j1,j2,ii,nn,oo
  integer,parameter :: pp(-1:1,-1:1)=&
                      reshape((/-2,-2,2 ,-2,0,2 ,-2,2,2/),(/3,3/))
!
  re1=areal(x1) ;re2=areal(x2)
  im1=aimag(x1) ;im2=aimag(x2)
!
  if (abs(im1).le.EPSN*abs(re1).and.abs(im2).le.EPSN*abs(re2)) then
    if (re1.ge.RZRO) then
      if (re2.ge.RZRO) then
        rslt = dilog2_r( re1,i1 , re2,i2 )
      else
        rslt = dilog2_r( re1,i1 ,-re2,i2+sgnRe(im2) )
      endif
    elseif (re2.ge.RZRO) then
      rslt = dilog2_r(-re1,i1+sgnRe(im1) , re2,i2 )
    else
      rslt = dilog2_r(-re1,i1+sgnRe(im1) ,-re2,i2+sgnRe(im2) )
    endif
    return
  endif
!
  if (re1.ge.RZRO) then ;r1= x1 ;j1=i1
                   else ;r1=-x1 ;j1=i1+sgnRe(im1,1)
  endif
  if (re2.ge.RZRO) then ;r2= x2 ;j2=i2
                   else ;r2=-x2 ;j2=i2+sgnRe(im2,1)
  endif
!
  a1=abs(r1) ;a2=abs(r2)
  if (a1.gt.a2) then
    aa=a1;a1=a2;a2=aa
    rr=r1;r1=r2;r2=rr
    ii=j1;j1=j2;j2=ii
  endif
!
  oo=mod(j1,2) ;nn=j1-oo ;y1=r1 ;if (oo.ne.0) y1=-y1
  oo=mod(j2,2) ;nn=j2-oo ;y2=r2 ;if (oo.ne.0) y2=-y2
!
  eps = 10*EPSN
!
  if (j1.ne.j2) then
    if (r1.eq.r2) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'j1,j2,r1-r2',j1,j2,',',trim(myprint(r1-r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_c j1=/=j2,r1=r2' !DEBUG
      return
    else
      rslt = ( dilog_c(r1,j1)-dilog_c(r2,j2) )/(y1-y2)
!      write(*,*) 'dilog2_c j1=/=j2' !DEBUG
      return
    endif
  endif
!
  if (a1.lt.eps) then
    if (a2.lt.eps) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r1,r2 =',trim(myprint(r1)),',',trim(myprint(r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_c r1<eps,r2<eps' !DEBUG
      return
    else
      rslt = (dilog_c(r2,j2)-4*PISQo24)/y2
!      write(*,*) 'dilog2_c r1<eps' !DEBUG
      return
    endif
  endif
!
  logr1=log(r1) ;logr2=log(r2)
!
  ao1=abs(1-y1) ;ao2=abs(1-y2)
  if (10*ao1.lt.RONE.or.10*ao2.lt.RONE) then
    aa = abs(r1/r2-1)
    if (10*aa.gt.RONE) then
      rslt = (dilog_c(r1,j1)-dilog_c(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_c ||1-y1|/|1-y2|-1|>0.1' !DEBUG
      return
    elseif (oo.eq.0.and.ao1.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r1,oo,nn =',trim(myprint(r1)),',',oo,nn,', putting nn=0'
      if (ao2.lt.eps) then
        rslt = -1
!        write(*,*) 'dilog2_c |1-y1|' !DEBUG
        return
      else
        y1=1-eps ;nn=0 ;logr1=0 ;r1=1-eps
      endif
    elseif (oo.eq.0.and.ao2.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r2,oo,nn =',trim(myprint(r2)),',',oo,nn,', putting nn=0'
      y2=1-eps ;nn=0 ;logr2=0 ;r2=1-eps
    endif
  else
    aa = abs((logr1+oo*IPI)/(logr2+oo*IPI)-1)
    if (10*aa.gt.RONE) then
      rslt = (dilog_c(r1,j1)-dilog_c(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_c |logr1/logr2-1|>0.1',logr1,logr2 !DEBUG
      return
    elseif (aa.lt.eps) then
      ii = 0
      if (a1.gt.RONE) ii = ii + (nn+pp(oo,sgnIm(y2)))
      if (a2.gt.RONE) ii = ii - (nn+pp(oo,sgnIm(y2)))
      ii = nn*ii
      if (ii.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r1,r2,nn =',trim(myprint(r1)),',',trim(myprint(r2)),',',nn &
        ,', putting nn=0'
      rslt = -olog1(y2,0)
!      write(*,*) 'dilog2_c |logr1/lorg2|<eps' !DEBUG
      return
    endif
  endif
!
  if (a1.gt.RONE) then
    y1=1/y1 ;logr1=-logr1
    y2=1/y2 ;logr2=-logr2
    nn=-nn ;oo=-oo
  endif
!
  ff=y1/y2         ;ff=-olog1(ff,0)/y2
  gg=(1-y1)/(1-y2) ;gg=-olog1(gg,0)/(1-y2)
!
  if (2*areal(y1).ge.RONE) then
!    write(*,*) 'dilog2_c re>1/2' !DEBUG
    rslt = ff*sumterms_c(-logr1,-logr2) - nn*IPI*gg
  else
!    write(*,*) 'dilog2_c re<1/2' !DEBUG
    logo1 = log(1-y1)
    logo2 = log(1-y2)
    rslt = gg*( sumterms_c(-logo1,-logo2) - (nn+oo)*IPI - logr2 ) + ff*logo1
  endif
!
  if (a1.gt.RONE) then !implies also r2>1
!    write(*,*) 'dilog2_c r1>1,r2>1' !DEBUG
    rslt = y1*y2*( rslt - ff*((logr1+logr2)/2 + (nn+oo)*IPI) )
  elseif (a2.gt.RONE.and.nn.ne.0) then
!    write(*,*) 'dilog2_c r1<1,r2>1',oo,sgnIm(y2)!DEBUG
    rslt = rslt - 12*nn*( nn + pp(oo,sgnIm(y2)) )*PISQo24/(y1-y2)
  endif
!
  end function


  function dilog2_r( x1,i1 ,x2,i2 ) result(rslt)
!*******************************************************************
!*******************************************************************
  use avh_olo_dp_olog
  real(kindr2) &  
          ,intent(in) :: x1,x2
  integer ,intent(in) :: i1,i2
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: y1,y2 ,ff,gg,logr1,logr2,logo1,logo2
  real(kindr2) &  
    :: eps,r1,r2,rr,ro1,ro2
  integer :: j1,j2,ii,nn,oo
!
  if (x1.ge.RZRO) then ;r1= x1 ;j1=i1
                  else ;r1=-x1 ;j1=i1+1 ! log(-1)=i*pi
  endif
  if (x2.ge.RZRO) then ;r2= x2 ;j2=i2
                  else ;r2=-x2 ;j2=i2+1 ! log(-1)=i*pi
  endif
!
  if (r1.gt.r2) then
    rr=r1;r1=r2;r2=rr
    ii=j1;j1=j2;j2=ii
  endif
!
  oo=mod(j1,2) ;nn=j1-oo ;y1=r1 ;if (oo.ne.0) y1=-y1
  oo=mod(j2,2) ;nn=j2-oo ;y2=r2 ;if (oo.ne.0) y2=-y2
!
  eps = 10*EPSN
!
  if (j1.ne.j2) then
    if (r1.eq.r2) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'j1,j2,r1-r2',j1,j2,',',trim(myprint(r1-r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_r j1=/=j2,r1=r2' !DEBUG
      return
    else
      rslt = ( dilog_r(r1,j1)-dilog_r(r2,j2) )/(y1-y2)
!      write(*,*) 'dilog2_r j1=/=j2' !DEBUG
      return
    endif
  endif
!
  if (r1.lt.eps) then
    if (r2.lt.eps) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r1,r2 =',trim(myprint(r1)),',',trim(myprint(r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_r r1<eps,r2<eps' !DEBUG
      return
    else
      rslt = (dilog_r(r2,j2)-4*PISQo24)/y2
!      write(*,*) 'dilog2_r r1<eps' !DEBUG
      return
    endif
  endif
!
  logr1=log(r1) ;logr2=log(r2)
!
  ro1=abs(1-y1) ;ro2=abs(1-y2)
  if (10*ro1.lt.RONE.or.10*ro2.lt.RONE) then
    rr = abs(r1/r2-1)
    if (10*rr.gt.RONE) then
      rslt = (dilog_r(r1,j1)-dilog_r(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_r ||1-y1|/|1-y2|-1|>0.1' !DEBUG
      return
    elseif (oo.eq.0.and.ro1.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r1,oo,nn =',trim(myprint(r1)),',',oo,nn,', putting nn=0'
      if (ro2.lt.eps) then
        rslt = -1
!        write(*,*) 'dilog2_r |1-y1|' !DEBUG
        return
      else
        y1=1-eps ;nn=0 ;logr1=0 ;r1=1-eps
      endif
    elseif (oo.eq.0.and.ro2.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r2,oo,nn =',trim(myprint(r2)),',',oo,nn,', putting nn=0'
      y2=1-eps ;nn=0 ;logr2=0 ;r2=1-eps
    endif
  else
    rr = abs((logr1+oo*IPI)/(logr2+oo*IPI)-1)
    if (10*rr.gt.RONE) then
      rslt = (dilog_r(r1,j1)-dilog_r(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_r |logr1/logr2-1|>0.1',logr1,logr2 !DEBUG
      return
    elseif (rr.lt.eps) then
      ii = 0
      if (r1.gt.RONE) ii = ii + (nn+2*oo)
      if (r2.gt.RONE) ii = ii - (nn+2*oo)
      ii = nn*ii
      if (ii.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r1,r2,nn =',trim(myprint(r1)),',',trim(myprint(r2)),',',nn &
        ,', putting nn=0'
      rslt = -olog1(y2,0)
!      write(*,*) 'dilog2_r |logr1/lorg2|<eps' !DEBUG
      return
    endif
  endif
!
  if (r1.gt.RONE) then
    y1=1/y1 ;logr1=-logr1
    y2=1/y2 ;logr2=-logr2
    nn=-nn ;oo=-oo
  endif
!
  ff=y1/y2         ;ff=-olog1(ff,0)/y2
  gg=(1-y1)/(1-y2) ;gg=-olog1(gg,0)/(1-y2)
!
  if (2*y1.ge.RONE) then
!    write(*,*) 'dilog2_r re>1/2' !DEBUG
    rslt = ff*sumterms_r(-logr1,-logr2) - nn*IPI*gg
  else
!    write(*,*) 'dilog2_r re<1/2' !DEBUG
    logo1 = log(1-y1)
    logo2 = log(1-y2)
    rslt = gg*( sumterms_r(-logo1,-logo2) - (nn+oo)*IPI - logr2 ) + ff*logo1
  endif
!
  if (r1.gt.RONE) then !implies also r2>1
!    write(*,*) 'dilog2_r r1>1,r2>1' !DEBUG
    rslt = y1*y2*( rslt - ff*((logr1+logr2)/2 + (nn+oo)*IPI) )
  elseif (r2.gt.RONE.and.nn.ne.0) then
!    write(*,*) 'dilog2_r r1<1,r2>1' !DEBUG
    rslt = rslt - 12*nn*PISQo24*(nn+2*oo)/(y1-y2)
  endif
!
  end function


  function sumterms_c( z1,z2 ) result(rslt)
!***********************************************************************
! ( f(z1)-f(z2) )/( z1-z2 ), where
! f(z)= z + c0*z^2 + c1*z^3 + c2*z^5 + c3*z^7 + ...
!***********************************************************************
  complex(kindr2) &   
    ,intent(in) :: z1,z2
  complex(kindr2) &   
    :: rslt,yy,zz
  real(kindr2) &  
    :: az
  integer :: ii,nn
  az = max(abs(z1),abs(z2))
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! calculates all z(i)=(z1^i-z2^i)/(z1-z2) numerically stable
!  zz(1) = 1
!  yy    = 1
!  do ii=2,2*nn+1
!    yy = z2*yy
!    zz(ii) = z1*zz(ii-1) + yy
!  enddo
  zz = 1
  yy = 1
  rslt = zz
  yy = z2*yy
  zz = z1*zz+yy
  rslt = rslt + coeff(0)*zz
  do ii=1,nn
    yy = z2*yy
    zz = z1*zz+yy
    rslt = rslt + coeff(ii)*zz
    yy = z2*yy
    zz = z1*zz+yy
  enddo
  end function  


  function sumterms_r( z1,z2 ) result(rslt)
!***********************************************************************
! ( f(z1)-f(z2) )/( z1-z2 ), where
! f(z)= z + c0*z^2 + c1*z^3 + c2*z^5 + c3*z^7 + ...
!***********************************************************************
  real(kindr2) &  
    ,intent(in) :: z1,z2
  real(kindr2) &  
    :: rslt,yy,zz
  real(kindr2) &  
    :: az
  integer :: ii,nn
  az = max(abs(z1),abs(z2))
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
  zz = 1
  yy = 1
  rslt = zz
  yy = z2*yy
  zz = z1*zz+yy
  rslt = rslt + coeff(0)*zz
  do ii=1,nn
    yy = z2*yy
    zz = z1*zz+yy
    rslt = rslt + coeff(ii)*zz
    yy = z2*yy
    zz = z1*zz+yy
  enddo
  end function  

end module


module avh_olo_dp_bnlog
!***********************************************************************
!                      /1    
!   bnlog(n,x) = (n+1) |  dt t^n ln(1-t/x) 
!                      /0 
!***********************************************************************
  use avh_olo_units
  use avh_olo_dp_prec
  use avh_olo_dp_auxfun
  use avh_olo_dp_arrays
  use avh_olo_dp_olog
  use avh_olo_dp_print
  implicit none
  private
  public :: update_bnlog,bnlog

  real(kindr2) &  
         ,allocatable,save :: coeff(:,:)
  real(kindr2) &  
         ,allocatable,save :: thrs(:,:,:)
  integer,allocatable,save :: ntrm(:,:,:)
  integer,parameter :: nStp=6
  integer,parameter :: rank=4
  integer,parameter :: aCoef(0:rank,0:rank)=reshape((/ &
                         1, 0, 0, 0, 0 & ! 1
                       , 1, 2, 0, 0, 0 & ! 1/2,1
                       , 2, 3, 6, 0, 0 & ! 1/3,1/2,1
                       , 3, 4, 6,12, 0 & ! 1/4,1/3,1/2,1
                       ,12,15,20,30,60 & ! 1/5,1/4,1/3,1/2,1
                       /),(/rank+1,rank+1/))

  interface bnlog
    module procedure bnlog_c,bnlog_r
  end interface

contains


  subroutine update_bnlog
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
    :: tt
  integer :: nn,ii,jj,n1,nmax,irank
  logical :: highestSoFar
!  real(kind(1d0)) :: xx(6) !DEBUG
!
  if (allocated(thrs)) then
    call shift3( thrs ,prcpar )
    call shift3( ntrm ,prcpar )
  else
    allocate(thrs(1:nStp,0:rank,1:1))
    allocate(ntrm(1:nStp,0:rank,1:1))
    if (prcpar.ne.1) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop update_bnlog'
      stop
    endif
  endif
!
  highestSoFar = prcpar.eq.ubound(ntrm,3)
!
  if (highestSoFar) then
    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(0:-1,0:2)) ! allocate at size=0
  endif
!
  nmax = 0
!
  do irank=0,rank
!
    n1 = 2+irank
!
    if (prcpar.gt.1) then ;nn=ntrm(nStp,irank,prcpar-1)-1
                     else ;nn=n1
    endif
!  
    do
      nn = nn+1
      if (highestSoFar.and.nn.gt.ubound(coeff,1)) call update_coeff( 2*nn )
      tt = 1
      tt = (EPSN*abs(coeff(n1,irank)/coeff(nn,irank)))**(tt/(nn-n1))
      if (8*(irank+1)*tt.gt.RONE) exit
    enddo
!
    if (nn.gt.nmax) nmax=nn
!  
    ntrm(nStp,irank,prcpar) = nn
    thrs(nStp,irank,prcpar) = tt
    nn = max(1,nint(nn*1d0/nStp))
    do ii=nStp-1,1,-1
      ntrm(ii,irank,prcpar) = ntrm(ii+1,irank,prcpar)-nn
      if (ntrm(ii,irank,prcpar).le.n1) then
        do jj=1,ii
          ntrm(jj,irank,prcpar) = max(n1,ntrm(ii,irank,prcpar))
          thrs(jj,irank,prcpar) = 0 
        enddo
        exit
      endif
      jj = ntrm(ii,irank,prcpar)
      tt = 1
      tt = (EPSN*abs(coeff(n1,irank)/coeff(jj,irank)))**(tt/(jj-n1))
      thrs(ii,irank,prcpar) = tt
    enddo
!  
  enddo!irank=1,nrank
!  
  if (highestSoFar) call resize( coeff ,2,nmax ,0,rank )
!
!  do ii=lbound(thrs,3),ubound(thrs,3)        !DEBUG
!  do irank=0,rank                            !DEBUG
!    do jj=1,nStp                             !DEBUG
!      xx(jj) = thrs(jj,irank,ii)             !DEBUG
!    enddo                                    !DEBUG
!    write(*,'(i2,99e10.3)') irank,xx(:)      !DEBUG
!    write(*,'(2x,99i10)'  ) ntrm(:,irank,ii) !DEBUG
!  enddo                                      !DEBUG
!  enddo                                      !DEBUG
  end subroutine


  subroutine update_coeff( ncf )
!*******************************************************************
! Coefficients of the expansion of
!   f(n,x) = -int( t^n*log(1-t*x) ,t=0..1 )
! in terms of log(1-x)
!*******************************************************************
  integer ,intent(in) :: ncf
  integer :: ii,jj
  real(kindr2) &  
    :: fact,tt(rank)
!
  call enlarge( coeff ,2,ncf ,0,rank )
!
  do jj=0,rank
  do ii=2,1+jj
    coeff(ii,jj) = 0
  enddo
  enddo
  fact = 1
  do ii=1,rank ;tt(ii)=1 ;enddo
  do ii=2,ncf
    fact = fact*ii
    coeff(ii,0) = (ii-1)/fact
    if (ii.eq.2) cycle
    do jj=1,rank ;tt(jj)=tt(jj)*(jj+1) ;enddo
    coeff(ii,1) = coeff(ii,0)*(1-tt(1))
    if (ii.eq.3) cycle
    coeff(ii,2) = coeff(ii,0)*(1-2*tt(1)+tt(2))
    if (ii.eq.4) cycle
    coeff(ii,3) = coeff(ii,0)*(1-3*tt(1)+3*tt(2)-tt(3))
    if (ii.eq.5) cycle
    coeff(ii,4) = coeff(ii,0)*(1-4*tt(1)+6*tt(2)-4*tt(3)+tt(4))
!   if (ii.eq.n+1) cycle
!   coeff(ii,n) = coeff(ii,0)
!               * ( 1 - binom(n,1)*tt(1) + binom(n,2)*tt(2)...)
  enddo
!
  end subroutine


  function bnlog_c( irank ,xx ) result(rslt)
!*******************************************************************
!*******************************************************************
  integer ,intent(in) :: irank
  complex(kindr2) &   
    ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt,yy,omx
  real(kindr2) &  
    :: aa,rex,imx
  integer :: ii,nn
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (abs(imx).le.EPSN*abs(rex)) then
    rslt = bnlog_r( irank ,rex ,sgnRe(imx,1) )
    return
  endif
!
  if (abs(xx-1).le.EPSN*10) then
    aa = 1
    rslt = -1
    do ii=2,irank+1
      rslt = rslt - aa/ii
    enddo
    return
  endif
!
  yy = olog(1-1/xx,0)
  aa = abs(yy)
  if     (aa.ge.thrs(6,irank,prcpar)) then
     omx = 1
    rslt = aCoef(irank,irank)
    do ii=irank,1,-1
       omx = 1 + xx*omx
      rslt = aCoef(ii-1,irank) + xx*rslt
    enddo
     omx = (1-xx)*omx
    rslt = omx*yy - rslt/aCoef(irank,irank)
!    if     (irank.eq.0) then
!      rslt = (1-xx)*yy - 1
!    elseif (irank.eq.1) then
!      rslt = (1-xx)*(1+xx)*yy - (1+xx*2)/2
!    elseif (irank.eq.2) then
!      rslt = (1-xx)*(1+xx*(1+xx))*yy - (2+xx*(3+xx*6))/6
!    elseif (irank.eq.3) then
!      rslt = (1-xx)*(1+xx*(1+xx*(1+xx)))*yy &
!           - (3+xx*(4+xx*(6+xx*12)))/12
!    elseif (irank.eq.4) then
!      rslt = (1-xx)*(1+xx*(1+xx*(1+xx*(1+xx))))*yy &
!           - (12+xx*(15+xx*(20+xx*(30+xx*60))))/60
!    endif
    return
  elseif (aa.ge.thrs(5,irank,prcpar)) then ;nn=ntrm(6,irank,prcpar)
  elseif (aa.ge.thrs(4,irank,prcpar)) then ;nn=ntrm(5,irank,prcpar)
  elseif (aa.ge.thrs(3,irank,prcpar)) then ;nn=ntrm(4,irank,prcpar)
  elseif (aa.ge.thrs(2,irank,prcpar)) then ;nn=ntrm(3,irank,prcpar)
  elseif (aa.ge.thrs(1,irank,prcpar)) then ;nn=ntrm(2,irank,prcpar)
                                      else ;nn=ntrm(1,irank,prcpar)
  endif
!
  rslt = coeff(nn,irank)
  do ii=nn-1,2+irank,-1
    rslt = coeff(ii,irank) + yy*rslt
  enddo
  rslt = -(irank+1)*rslt*yy*(yy*xx)**(irank+1)
!
  aa = areal(rslt)
  if (abs(aimag(rslt)).le.EPSN*abs(aa)) rslt = acmplx(aa)
!
  end function


  function bnlog_r( irank ,xx ,sgn ) result(rslt)
!*******************************************************************
!*******************************************************************
  integer ,intent(in) :: irank
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: sgn
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: yy,aa,omx
  integer :: ii,nn
  logical :: y_lt_0
!
  if (abs(xx).eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop bnlog_r: ' &
      ,'argument xx=',trim(myprint(xx,8)),', returning 0'
    rslt = 0
    return
  elseif (abs(xx-1).le.EPSN*10) then
    aa = 1
    rslt = -1
    do ii=2,irank+1
      rslt = rslt - aa/ii
    enddo
    return
  endif
!
  yy = 1-1/xx
  y_lt_0 = (yy.lt.RZRO)
  if (y_lt_0) then 
    yy = log(-yy)
    aa = sqrt(yy*yy+ONEPI*ONEPI)
  else
    yy = log( yy)
    aa = abs(yy)
  endif
!
  omx = 1
  do ii=irank,1,-1
    omx = 1+xx*omx
  enddo
  omx = (1-xx)*omx ! (1-x^{rank+1})
!
  if     (aa.ge.thrs(6,irank,prcpar)) then
    rslt = aCoef(irank,irank)
    do ii=irank,1,-1
      rslt = aCoef(ii-1,irank) + xx*rslt
    enddo
    rslt = omx*yy - rslt/aCoef(irank,irank)
!    if     (irank.eq.0) then
!      rslt = omx*yy - 1
!    elseif (irank.eq.1) then
!      rslt = omx*yy - (1+xx*2)/2
!    elseif (irank.eq.2) then
!      rslt = omx*yy - (2+xx*(3+xx*6))/6
!    elseif (irank.eq.3) then
!      rslt = omx*yy - (3+xx*(4+xx*(6+xx*12)))/12
!    elseif (irank.eq.4) then
!      rslt = omx*yy - (12+xx*(15+xx*(20+xx*(30+xx*60))))/60
!    endif
    if (y_lt_0) rslt = rslt + sgn*omx*IPI
    return
  elseif (aa.ge.thrs(5,irank,prcpar)) then ;nn=ntrm(6,irank,prcpar)
  elseif (aa.ge.thrs(4,irank,prcpar)) then ;nn=ntrm(5,irank,prcpar)
  elseif (aa.ge.thrs(3,irank,prcpar)) then ;nn=ntrm(4,irank,prcpar)
  elseif (aa.ge.thrs(2,irank,prcpar)) then ;nn=ntrm(3,irank,prcpar)
  elseif (aa.ge.thrs(1,irank,prcpar)) then ;nn=ntrm(2,irank,prcpar)
                                      else ;nn=ntrm(1,irank,prcpar)
  endif
!
  aa = coeff(nn,irank)
  do ii=nn-1,2+irank,-1
    aa = coeff(ii,irank) + yy*aa
  enddo
  rslt = -(irank+1)*aa*yy*(yy*xx)**(irank+1)
  if (y_lt_0) rslt = rslt + sgn*omx*IPI
!  
  end function

end module


module avh_olo_dp_qmplx
  use avh_olo_units
  use avh_olo_dp_prec
  use avh_olo_dp_auxfun
  use avh_olo_dp_olog
  use avh_olo_dp_dilog

  implicit none
  private
  public :: qmplx_type,qonv,directly,sheet,logc,logc2,logc3,li2c,li2c2
  public :: operator (*) ,operator (/)

  type :: qmplx_type
  complex(kindr2) &   
          :: c
  integer :: p
  end type

  interface qonv
    module procedure qonv_cr,qonv_ci,qonv_c,qonv_i
  end interface

  interface operator (*)
    module procedure prduct_qq,prduct_qr
  end interface
  interface operator (/)
    module procedure ratio_qq,ratio_qr
  end interface

contains


  function qonv_cr(xx,sgn) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz  becomes the
! sign of  sgn .
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  real(kindr2) &  
    ,intent(in) :: sgn
  type(qmplx_type) :: rslt
  real(kindr2) &  
    :: xre,xim
  xre = areal(xx)
  if (xre.ge.RZRO) then
    rslt%c = xx
    rslt%p = 0
  else
    xim = aimag(xx)
    if (xim.eq.RZRO) then
      rslt%c = -xre
      rslt%p = sgnRe(sgn)
    else
      rslt%c = -xx
      rslt%p = sgnRe(xim) ! xim = -Im(rslt%c)
    endif
  endif
  end function

  function qonv_ci(xx,sgn) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz  becomes the
! sign of  sgn .
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  integer         ,intent(in) :: sgn
  type(qmplx_type) :: rslt
  real(kindr2) &  
    :: xre,xim
  xre = areal(xx)
  if (xre.ge.RZRO) then
    rslt%c = xx
    rslt%p = 0
  else
    xim = aimag(xx)
    if (xim.eq.RZRO) then
      rslt%c = -xre
      rslt%p = sign(1,sgn)
    else
      rslt%c = -xx
      rslt%p = sgnRe(xim) ! xim = -Im(rslt%c)
    endif
  endif
  end function

  function qonv_c(xx) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz=1
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  type(qmplx_type) :: rslt
  real(kindr2) &  
    :: xre,xim
  xre = areal(xx)
  if (xre.ge.RZRO) then
    rslt%c = xx
    rslt%p = 0
  else
    xim = aimag(xx)
    if (xim.eq.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop qonv_c: ' &
        ,'negative input with undefined sign for the imaginary part, ' &
        ,'putting +ieps'
      rslt%c = -xre
      rslt%p = 1
    else
      rslt%c = -xx
      rslt%p = sgnRe(xim) ! xim = -Im(rslt%c)
    endif
  endif
  end function

  function qonv_i(xx) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz=1
!*******************************************************************
  integer ,intent(in) :: xx
  type(qmplx_type) :: rslt
  if (xx.ge.0) then
    rslt%c = xx
    rslt%p = 0
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop qonv_i: ' &
      ,'negative input with undefined sign for the imaginary part, ' &
      ,'putting +ieps'
    rslt%c = -xx
    rslt%p = 1
  endif
  end function

  function directly(xx,ix) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  integer         ,intent(in) :: ix
  type(qmplx_type) :: rslt
  rslt%c = xx
  rslt%p = ix
  end function


  function sheet(xx) result(ii)
!*******************************************************************
! Returns the number of the Riemann-sheet (times 2) for the complex
! number  xx*exp(ix*imag*pi) . The real part of xx is assumed to be
! positive or zero. Examples:
! xx=1+imag, ix=-1 -> ii= 0 
! xx=1+imag, ix= 1 -> ii= 2 
! xx=1-imag, ix=-1 -> ii=-2 
! xx=1-imag, ix= 1 -> ii= 0 
! xx=1     , ix= 1 -> ii= 0  convention that log(-1)=pi on
! xx=1     , ix=-1 -> ii=-2  the principal Riemann-sheet
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  integer :: ii,jj
  real(kindr2) &  
    :: xim
  jj = mod(xx%p,2)
  ii = xx%p-jj
  xim = aimag(xx%c)
  if (xim.le.RZRO) then ! also xim=0 <==> log(-1)=pi, not -pi
    if (jj.eq.-1) ii = ii-2
  else
    if (jj.eq. 1) ii = ii+2
  endif
  end function


  function prduct_qq(yy,xx) result(zz)
!*******************************************************************
! Return the product  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy,xx
  type(qmplx_type) :: zz
  zz%c = yy%c*xx%c
  zz%p = yy%p+xx%p
  if (areal(zz%c).lt.RZRO) then
    zz%p = zz%p + sgnIm(xx%c)
    zz%c = -zz%c
  endif
  end function

  function prduct_qr(yy,xx) result(zz)
!*******************************************************************
! Return the product  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy
  real(kindr2) &  
    ,intent(in) :: xx
  type(qmplx_type) :: zz
  zz%c = yy%c*abs(xx)
  zz%p = yy%p
  end function

  function ratio_qq(yy,xx) result(zz)
!*******************************************************************
! Return the ratio  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy,xx
  type(qmplx_type) :: zz
  zz%c = yy%c/xx%c
  zz%p = yy%p-xx%p
  if (areal(zz%c).lt.RZRO) then
    zz%p = zz%p - sgnIm(xx%c)
    zz%c = -zz%c
  endif
  end function

  function ratio_qr(yy,xx) result(zz)
!*******************************************************************
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy
  real(kindr2) &  
    ,intent(in) :: xx
  type(qmplx_type) :: zz
  zz%c = yy%c/abs(xx)
  zz%p = yy%p
  end function


  function logc(xx) result(rslt)
!*******************************************************************
! log(xx)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = olog(acmplx(xx%c),xx%p)
  rslt = olog(xx%c,xx%p)
  end function

  function logc2(xx) result(rslt)
!*******************************************************************
! log(xx)/(1-xx)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = -olog1(acmplx(xx%c),xx%p)
  rslt = -olog1(xx%c,xx%p)
  end function

  function logc3(xx) result(rslt)
!*******************************************************************
!  ( log(xx)/(1-xx) + 1 )/(1-xx)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = olog2(acmplx(xx%c),xx%p)
  rslt = olog2(xx%c,xx%p)
  end function

  function li2c(xx) result(rslt)
!*******************************************************************
!    /1    ln(1-(1-xx)*t)
!  - |  dt -------------- 
!    /0        t
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = dilog(acmplx(xx%c),xx%p)
  rslt = dilog(xx%c,xx%p)
  end function

  function li2c2(xx,yy) result(rslt)
!*******************************************************************
! ( li2(xx) - li2(yy) )/(xx-yy)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx,yy
  complex(kindr2) &   
    :: rslt
!  rslt = dilog( acmplx(xx%c),xx%p ,acmplx(yy%c),yy%p )
!  write(*,*) 'li2c2 x:',xx%c,xx%p !DEBUG
!  write(*,*) 'li2c2 y:',yy%c,yy%p !DEBUG
  rslt = dilog( xx%c,xx%p ,yy%c,yy%p )
!  write(*,*) 'li2c2 out:',rslt !DEBUG
  end function


end module


module avh_olo_dp_bub
  use avh_olo_units
  use avh_olo_dp_prec
  use avh_olo_dp_auxfun
  use avh_olo_dp_bnlog
  use avh_olo_dp_qmplx
  use avh_olo_dp_olog
  implicit none
  private
  public :: tadp ,tadpn ,bub0 ,dbub0 ,bub1 ,bub11 ,bub111 ,bub1111

contains

  subroutine tadp( rslt ,mm ,amm ,rmu2 )
!*******************************************************************
! The 1-loop scalar 1-point function.
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
    ,intent(in)  :: amm,rmu2
!
  rslt(2) = 0
  if (amm.eq.RZRO.or.mm.eq.CZRO) then
    rslt(1) = 0
    rslt(0) = 0
  else
    rslt(1) = mm
    rslt(0) = mm - mm*logc( qonv(mm/rmu2,-1) )
  endif
  end subroutine


  subroutine tadpn( rslt ,rank ,mm ,amm ,rmu2 )
!*******************************************************************
! The 1-loop tensor 1-point functions.
!   rslt(:,0) = A0
!   rslt(:,1) = A00
!   rslt(:,2) = A0000  etc.
! For input  rank  only  rslt(:,0:rank/2)  is filled.
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
    ,intent(in)  :: amm,rmu2
  integer ,intent(in) :: rank
  complex(kindr2) &   
    :: aa
  real(kindr2) &  
    :: bb
  integer :: ii
!
  do ii=0,rank
    rslt(2,ii) = 0
    rslt(1,ii) = 0
    rslt(0,ii) = 0
  enddo
  if (amm.eq.RZRO.or.mm.eq.CZRO) then
    return
  else
    rslt(1,0) = mm
    rslt(0,0) = mm - mm*logc( qonv(mm/rmu2,-1) )
    aa = 1
    bb = 0
    do ii=1,rank/2
      aa = aa*mm/(2*(ii+1))
      bb = bb + RONE/(ii+1)
      rslt(1,ii) = aa*( rslt(1,0) )
      rslt(0,ii) = aa*( rslt(0,0) + mm*bb )
    enddo
  endif
  end subroutine


!*******************************************************************
! Return the Passarino-Veltman functions
!
!      C   /      d^(Dim)q
!   ------ | -------------------- = b0
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!      C   /    d^(Dim)q q^mu
!   ------ | -------------------- = p^mu b1
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!      C   /  d^(Dim)q q^mu q^nu
!   ------ | -------------------- = g^{mu,nu} b00 + p^mu p^nu b11
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!   etc.
!
! Based on the formulas from
! A. Denner, M. Dittmaier, Nucl.Phys. B734 (2006) 62-115
!*******************************************************************

  subroutine bub0( b0 &
                  ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
  endif
!
  b0(2) = 0
  b0(1) = CONE
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
  endif
!
  end subroutine

  subroutine bub1( b1,b0 &
                  ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  end subroutine

  subroutine bub11( b11,b00,b1,b0 &
                   ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      b00(0)=0 ;b00(1)=0 ;b00(2)=0 
      b11(0)=0 ;b11(1)=0 ;b11(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
  b11(2) = 0      
  b11(1) = CONE/3 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
      b11(0) = lna/3 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
      b11(0) = ( lna - bnlog(2,x1) )/3 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
      b11(0) = ( lna*3 + 2 )/9 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
      b11(0) = ( (lna-bnlog(2,x1))*3 + 1 )/9 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    b11(0) = ( lna - bnlog(2,x1) - bnlog(2,x2) )/3 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  rank = 0 
  call tadpn( a0 ,rank ,m1 ,am1 ,rmu2 )
  x1 = (m1-m0)-pp
  x2 = 2*m0
  b00(2) = 0
  b00(1) = ( a0(1,0) - x1*b1(1) + x2*b0(1) )/6
  b00(0) = ( a0(0,0) - x1*b1(0) + x2*b0(0) + 4*b00(1) )/6
  end subroutine

  subroutine bub111( b111,b001,b11,b00,b1,b0 &
                    ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b111(0:2),b001(0:2),b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      b00(0)=0 ;b00(1)=0 ;b00(2)=0 
      b11(0)=0 ;b11(1)=0 ;b11(2)=0 
      b001(0)=0 ;b001(1)=0 ;b001(2)=0 
      b111(0)=0 ;b111(1)=0 ;b111(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
  b11(2) = 0      
  b11(1) = CONE/3 
  b111(2) = 0      
  b111(1) =-CONE/4 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
      b11(0) = lna/3 
      b111(0) =-lna/4 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
      b11(0) = ( lna - bnlog(2,x1) )/3 
      b111(0) =-( lna - bnlog(3,x1) )/4 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
      b11(0) = ( lna*3 + 2 )/9 
      b111(0) =-( lna*4 + 2 )/16 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
      b11(0) = ( (lna-bnlog(2,x1))*3 + 1 )/9 
      b111(0) =-( (lna-bnlog(3,x1))*4 + 1 )/16 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    b11(0) = ( lna - bnlog(2,x1) - bnlog(2,x2) )/3 
    b111(0) =-( lna - bnlog(3,x1) - bnlog(3,x2) )/4 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  rank = 0 
  rank = 1 
  call tadpn( a0 ,rank ,m1 ,am1 ,rmu2 )
  x1 = (m1-m0)-pp
  x2 = 2*m0
  b00(2) = 0
  b00(1) = ( a0(1,0) - x1*b1(1) + x2*b0(1) )/6
  b00(0) = ( a0(0,0) - x1*b1(0) + x2*b0(0) + 4*b00(1) )/6
  b001(2) = 0
  b001(1) = (-a0(1,0) - x1*b11(1) + x2*b1(1) )/8
  b001(0) = (-a0(0,0) - x1*b11(0) + x2*b1(0) + 4*b001(1) )/8
  end subroutine

  subroutine bub1111( b1111,b0011,b0000,b111,b001,b11,b00,b1,b0 &
                    ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b1111(0:2),b0011(0:2),b0000(0:2) &
                   ,b111(0:2),b001(0:2),b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      b00(0)=0 ;b00(1)=0 ;b00(2)=0 
      b11(0)=0 ;b11(1)=0 ;b11(2)=0 
      b001(0)=0 ;b001(1)=0 ;b001(2)=0 
      b111(0)=0 ;b111(1)=0 ;b111(2)=0 
      b0000(0)=0 ;b0000(1)=0 ;b0000(2)=0 
      b0011(0)=0 ;b0011(1)=0 ;b0011(2)=0 
      b1111(0)=0 ;b1111(1)=0 ;b1111(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
  b11(2) = 0      
  b11(1) = CONE/3 
  b111(2) = 0      
  b111(1) =-CONE/4 
  b1111(2) = 0      
  b1111(1) = CONE/5 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
      b11(0) = lna/3 
      b111(0) =-lna/4 
      b1111(0) = lna/5 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
      b11(0) = ( lna - bnlog(2,x1) )/3 
      b111(0) =-( lna - bnlog(3,x1) )/4 
      b1111(0) = ( lna - bnlog(4,x1) )/5 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b1111(0) =-b1111(0)-4*b111(0)-6*b11(0)-4*b1(0)-b0(0) 
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
      b11(0) = ( lna*3 + 2 )/9 
      b111(0) =-( lna*4 + 2 )/16 
      b1111(0) = ( lna*5 + 2 )/25 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
      b11(0) = ( (lna-bnlog(2,x1))*3 + 1 )/9 
      b111(0) =-( (lna-bnlog(3,x1))*4 + 1 )/16 
      b1111(0) = ( (lna-bnlog(4,x1))*5 + 1 )/25 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1111(0) =-b1111(0)-4*b111(0)-6*b11(0)-4*b1(0)-b0(0) 
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    b11(0) = ( lna - bnlog(2,x1) - bnlog(2,x2) )/3 
    b111(0) =-( lna - bnlog(3,x1) - bnlog(3,x2) )/4 
    b1111(0) = ( lna - bnlog(4,x1) - bnlog(4,x2) )/5 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1111(0) =-b1111(0)-4*b111(0)-6*b11(0)-4*b1(0)-b0(0) 
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  rank = 0 
  rank = 1 
  rank = 2 
  call tadpn( a0 ,rank ,m1 ,am1 ,rmu2 )
  x1 = (m1-m0)-pp
  x2 = 2*m0
  b00(2) = 0
  b00(1) = ( a0(1,0) - x1*b1(1) + x2*b0(1) )/6
  b00(0) = ( a0(0,0) - x1*b1(0) + x2*b0(0) + 4*b00(1) )/6
  b001(2) = 0
  b001(1) = (-a0(1,0) - x1*b11(1) + x2*b1(1) )/8
  b001(0) = (-a0(0,0) - x1*b11(0) + x2*b1(0) + 4*b001(1) )/8
  b0000(2) = 0
  b0000(1) = ( a0(1,1) - x1*b001(1) + x2*b00(1) )/10
  b0000(0) = ( a0(0,1) - x1*b001(0) + x2*b00(0) + 4*b0000(1) )/10
  b0011(2) = 0
  b0011(1) = ( a0(1,0) - x1*b111(1) + x2*b11(1) )/10
  b0011(0) = ( a0(0,0) - x1*b111(0) + x2*b11(0) + 4*b0011(1) )/10
  end subroutine


!*******************************************************************
! Derivative of B0
! expects  m0<m1
! only finite case, so input must not be  m0=0 & m1=pp
!*******************************************************************

  subroutine dbub0( rslt &
                   ,pp,m0,m1 ,app,am0,am1 )
  complex(kindr2) &   
    ,intent(out) :: rslt
  complex(kindr2) &   
    ,intent(in)  :: pp,m0,m1
  real(kindr2) &  
    ,intent(in)  :: app,am0,am1
  complex(kindr2) &   
    :: ch,x1,x2,lambda
  real(kindr2) &  
    :: ax1,ax2,ax1x2,maxa
  type(qmplx_type) :: q1,q2,q1o,q2o
  integer :: sgn
!
  if (am1.eq.RZRO) then
    if (app.eq.RZRO) then
      rslt = 0
      return
    endif
  endif
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      rslt = 1/(6*m1)
    else
      ch = m0/m1
      rslt = ( CONE/2 - ch*olog3(ch,0) )/m1 
    endif
  elseif (am1.eq.RZRO) then
    rslt =-1/pp
  else
    call solabc( x1,x2 ,lambda ,pp ,(m0-m1)-pp ,m1 ,0 )
    sgn =-sgnRe(pp)*sgnRe(x2-x1)
    q1  = qonv(x1  , sgn)
    q1o = qonv(x1-1, sgn)
    q2  = qonv(x2  ,-sgn)
    q2o = qonv(x2-1,-sgn)
    ax1 = abs(x1)
    ax2 = abs(x2)
    ax1x2 = abs(x1-x2)
    maxa = max(ax1,ax2)
    if (ax1x2.lt.maxa*EPSN*10) then
      rslt = ( (x1+x2-1)*logc(q2/q2o) - 2 )/pp
    elseif (ax1x2*2.lt.maxa) then
      if     (x1.eq.CZRO.or.x1.eq.CONE) then
        rslt = ( (x1+x2-1)*logc(q2/q2o) - 1 )/pp
      elseif (x2.eq.CZRO.or.x2.eq.CONE) then
        rslt = ( (x1+x2-1)*logc(q1/q1o) - 1 )/pp
      else
        rslt = x1*(x1-1)*( logc2(q1o/q2o)/(x2-1) - logc2(q1/q2)/x2 ) &
             + (x1+x2-1)*logc(q2/q2o) - 1
        rslt = rslt/pp
      endif
    else
      rslt = 0
      if (ax1.ne.RZRO) then
        if (ax1.lt.2*RONE) then
          rslt = rslt - x1
          if (x1.ne.CONE) rslt = rslt - x1*logc2(q1/q1o)
        else
          rslt = rslt + x1/(x1-1)*logc3(q1/q1o)
        endif
      endif
      if (ax2.ne.RZRO) then
        if (ax2.lt.2*RONE) then
          rslt = rslt + x2
          if (x2.ne.CONE) rslt = rslt + x2*logc2(q2/q2o)
        else
          rslt = rslt - x2/(x2-1)*logc3(q2/q2o)
        endif
      endif
      rslt = rslt/lambda
    endif
  endif
!
  end subroutine


end module


module avh_olo_dp_tri
  use avh_olo_units
  use avh_olo_dp_prec
  use avh_olo_dp_auxfun
  use avh_olo_dp_qmplx
  implicit none
  private
  public :: tria0,tria1,tria2,tria3,tria4,trif0,trif1,trif2,trif3 &
           ,trif3HV &
           ,permtable,casetable,base
  integer ,parameter :: permtable(3,0:7)=reshape((/ &
       1,2,3 &! 0, 0 masses non-zero, no permutation
      ,1,2,3 &! 1, 1 mass non-zero,   no permutation
      ,3,1,2 &! 2, 1 mass non-zero,   1 cyclic permutation
      ,1,2,3 &! 3, 2 masses non-zero, no permutation
      ,2,3,1 &! 4, 1 mass non-zero,   2 cyclic permutations
      ,2,3,1 &! 5, 2 masses non-zero, 2 cyclic permutations
      ,3,1,2 &! 6, 2 masses non-zero, 1 cyclic permutation
      ,1,2,3 &! 7, 3 masses non-zero, no permutation
      /) ,(/3,8/))                     ! 0,1,2,3,4,5,6,7
  integer ,parameter :: casetable(0:7)=(/0,1,1,2,1,2,2,3/)
  integer ,parameter :: base(3)=(/4,2,1/)

contains

   subroutine tria4( rslt ,cpp,cm2,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /             d^(Dim)q
!            ------ | ----------------------------------
!            i*pi^2 / q^2 [(q+k1)^2-m2] [(q+k1+k2)^2-m3]
!
! with  k1^2=m2, k2^2=pp, (k1+k2)^2=m3.
! m2,m3 should NOT be identically 0d0.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cm2,cm3,cpp
  real(kindr2) &  
     ,intent(in)  :: rmu2
   type(qmplx_type) :: q23,qm3,q32
  complex(kindr2) &   
     :: sm2,sm3,k23,r23,d23,cc
!
   sm2 = mysqrt(cm2)
   sm3 = mysqrt(cm3)
   k23 = (cm2+cm3-cpp)/(sm2*sm3)
   call rfun( r23,d23, k23 )
   if (r23.eq.-CONE) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop tria4: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   q23 = qonv(r23,-1)
   qm3 = qonv(cm3/rmu2,-1)
   q32 = qonv(sm3)/qonv(sm2)
!
   rslt(2) = 0
   cc = logc2(q23) * r23/(1+r23)/(sm2*sm3)
   rslt(1) = -cc
   rslt(0) = cc*( logc(qm3) - logc(q23) ) &
           - li2c2(q32*q23,q32/q23) / cm2 &
           + li2c2(q23*q23,qonv(1)) * r23/(sm2*sm3)
   end subroutine


   subroutine tria3( rslt ,cp2,cp3,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  p2=k2^2, p3=(k1+k2)^2.
! mm should NOT be identically 0d0,
! and p2 NOR p3 should be identical to mm. 
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp2,cp3,cm3
  real(kindr2) &  
     ,intent(in)  :: rmu2
   type(qmplx_type) :: q13,q23,qm3,x1,x2
  complex(kindr2) &   
     :: r13,r23
!
   r13 = cm3-cp3
   r23 = cm3-cp2
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   qm3 = qonv(cm3,-1)
   x1 = q23/qm3
   x2 = q13/qm3
   rslt(2) = 0
   rslt(1) = -logc2( q23/q13 )/r13
   rslt(0) = -li2c2( x1,x2 )/cm3 &
           - rslt(1)*( logc(x1*x2)+logc(qm3/rmu2) )
   end subroutine


   subroutine tria2( rslt ,cp3,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  k1^2 = 0 , k2^2 = m3  and  (k1+k2)^2 = p3.
! mm should NOT be identically 0d0,
! and pp should NOT be identical to mm. 
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cm3
  real(kindr2) &  
     ,intent(in)  :: rmu2
   type(qmplx_type) :: q13,qm3,qxx
  complex(kindr2) &   
     :: r13,logm,z2,z1,z0,cc
!
   r13 = cm3-cp3
   q13 = qonv(r13,-1)
   qm3 = qonv(cm3,-1)
   logm = logc( qm3/rmu2 )
   qxx = qm3/q13
   z2 = 1 
   z2 = z2/2
   z1 = logc(qxx)
   z0 = PISQo24 + z1*z1/2 - li2c(qxx)
   cc = -1/r13
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine tria1( rslt ,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  k1^2 = (k1+k2)^2 = m3.
! mm should NOT be identically 0d0.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cm3
  real(kindr2) &  
     ,intent(in)  :: rmu2
  complex(kindr2) &   
     :: zm
!
   zm = 1/(2*cm3)
   rslt(2) = 0
   rslt(1) = -zm
   rslt(0) = zm*( 2 + logc(qonv(cm3/rmu2,-1)) )
   end subroutine


   subroutine tria0( rslt ,cp ,ap ,rmu2 )
!*******************************************************************
! calculates
!               C   /         d^(Dim)q
!            ------ | ------------------------
!            i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  p1 = k1^2,  p2 = k2^2,  p3 = k3^2
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! If any of these numbers is IDENTICALLY 0d0, the corresponding
! IR-singular case is returned.
!*******************************************************************
   use avh_olo_dp_olog
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp(3)
  real(kindr2) &  
     ,intent(in)  :: ap(3),rmu2
  real(kindr2) &  
     :: pp(3),rp1,rp2,rp3
  complex(kindr2) &   
     :: log2,log3
   integer :: icase,i1,i2,i3
!
   pp(1)=areal(cp(1))
   pp(2)=areal(cp(2))
   pp(3)=areal(cp(3))
!
   icase = 0
   if (ap(1).gt.RZRO) icase = icase + base(1)
   if (ap(2).gt.RZRO) icase = icase + base(2)
   if (ap(3).gt.RZRO) icase = icase + base(3)
   rp1 = pp(permtable(1,icase))
   rp2 = pp(permtable(2,icase))
   rp3 = pp(permtable(3,icase))
   icase  = casetable(  icase)
!
   i1=0 ;if (-rp1.lt.RZRO) i1=-1
   i2=0 ;if (-rp2.lt.RZRO) i2=-1
   i3=0 ;if (-rp3.lt.RZRO) i3=-1
!
   if     (icase.eq.0) then
! 0 masses non-zero
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop tria0: ' &
       ,'all external masses equal zero, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
   elseif (icase.eq.1) then
! 1 mass non-zero
    log3 = olog( abs(rp3/rmu2) ,i3 )
    rslt(2) = 1/rp3
    rslt(1) = -log3/rp3
    rslt(0) = ( log3**2/2 - 2*PISQo24 )/rp3
  elseif (icase.eq.2) then
! 2 masses non-zero
    log2 = olog( abs(rp2/rmu2) ,i2 )
    log3 = olog( abs(rp3/rmu2) ,i3 )
    rslt(2) = 0
    rslt(1) = -olog1( abs(rp3/rp2) ,i3-i2 )/rp2
    rslt(0) = -rslt(1)*(log3+log2)/2
  elseif (icase.eq.3) then
! 3 masses non-zero
    call trif0( rslt ,cp(1),cp(2),cp(3) )
  endif
  end subroutine


   subroutine trif0( rslt ,p1,p2,p3 )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! equal zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1,p2,p3
   type(qmplx_type) :: q23,q24,q34,qx1,qx2
  complex(kindr2) &   
     :: r23,r24,r34,aa,bb,cc,dd,x1,x2
  real(kindr2) &  
     :: hh
!
   r23 = -p1
   r24 = -p3
   r34 = -p2
!
   aa = r34*r24
   bb = r24 + r34 - r23
   cc = 1
   hh = areal(r23)
   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,1 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1, hh)
   qx2 = qonv(x2,-hh)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   rslt(0) = li2c2( qx1*q34 ,qx2*q34 )*r34 &
           + li2c2( qx1*q24 ,qx2*q24 )*r24 &
           - logc2( qx1/qx2 )*logc( qx1*qx2 )/(x2*2) &
           - logc2( qx1/qx2 )*logc( q23 )/x2
!
   rslt(0) = rslt(0)/aa
   end subroutine


   subroutine trif1( rslt ,p1i,p2i,p3i ,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with one internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1i,p2i,p3i ,m3i 
   type(qmplx_type) :: q23,q24,q34,qm4,qx1,qx2,qss
  complex(kindr2) &   
     :: p2,p3,p4,p12,p23,m4,sm2,sm3,sm4 &
                     ,aa,bb,cc,dd,x1,x2,r23,r24,r34
  real(kindr2) &  
     :: mhh
   logical :: r24Not0,r34Not0
!
!   p1 = nul
   p2 = p1i
   p3 = p2i
   p4 = p3i
   p12 = p1i
   p23 = p3i
!   m1 = infinite
!   m2 = m1i = 0
!   m3 = m2i = 0
   m4 = m3i
!
   sm4 = mysqrt(m4)
   mhh = abs(sm4)
   sm3 = mhh
   sm2 = sm3
!
   r23 = (   -p2 -p2 *IEPS )/(sm2*sm3)
   r24 = ( m4-p23-p23*IEPS )/(sm2*sm4)
   r34 = ( m4-p3 -p3 *IEPS )/(sm3*sm4)
!
   r24Not0 = (abs(areal(r24))+abs(aimag(r24)).ge.neglig(prcpar))
   r34Not0 = (abs(areal(r34))+abs(aimag(r34)).ge.neglig(prcpar))
!
   aa = r34*r24 - r23
!
   if (aa.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif1: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   bb = r24/sm3 + r34/sm2 - r23/sm4
   cc = 1/(sm2*sm3)
!   hh = areal(r23)
!   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,0 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1) ! x2 SHOULD HAVE im. part
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
   qm4 = qonv(sm4,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   rslt(0) = -logc2( qx1/qx2 )*logc( qx1*qx2/(qm4*qm4) )/(x2*2) &
             -li2c2( qx1*qm4 ,qx2*qm4 )*sm4
!
   if (r34Not0) then
     qss = q34*mhh
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34*sm3
   endif
!
   if (r24Not0) then
     qss = q24*mhh
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r24*sm2
   endif
!
   rslt(0) = rslt(0) - logc2( qx1/qx2 )*logc( q23*(mhh*mhh) )/x2
!
   rslt(0) = rslt(0)/(aa*sm2*sm3*sm4)
   end subroutine


   subroutine trif2( rslt ,p1i,p2i,p3i ,m2i,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with two internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1i,p2i,p3i ,m2i,m3i
   type(qmplx_type) :: q23,q34,q24,qm2,qm3,qm4,qx1,qx2,qss,qy1,qy2
  complex(kindr2) &   
     :: p2,p3,p23,m2,m4,sm2,sm3,sm4,aa,bb,cc,dd,x1,x2 &
                     ,r23,k24,r34,r24,d24
   logical :: r23Not0,r34Not0
!
!   p1 = nul
   p2 = p3i
   p3 = p1i
!   p4 = p2i
!   p12 = p3i
   p23 = p2i
!   m1 = infinite
   m2 = m3i
!   m3 = m1i = 0
   m4 = m2i
!
!   sm1 = infinite
   sm2 = mysqrt(m2)
   sm3 = abs(sm2) !mysqrt(m3)
   sm4 = mysqrt(m4)
!
   r23 = (    m2-p2 -p2 *IEPS )/(sm2*sm3) ! p2
   k24 = ( m2+m4-p23-p23*IEPS )/(sm2*sm4) ! p2+p3
   r34 = (    m4-p3 -p3 *IEPS )/(sm3*sm4) ! p3
!
   r23Not0 = (abs(areal(r23))+abs(aimag(r23)).ge.neglig(prcpar))
   r34Not0 = (abs(areal(r34))+abs(aimag(r34)).ge.neglig(prcpar))
!
   call rfun( r24,d24 ,k24 )
!
   aa = r34/r24 - r23
!
   if (aa.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif2: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   bb = -d24/sm3 + r34/sm2 - r23/sm4
   cc = (sm4/sm2 - r24)/(sm3*sm4)
!   hh = areal(r23 - r24*r34)
!   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc(x1,x2,dd ,aa,bb,cc ,0)
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1 ) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1 ) ! x2 SHOULD HAVE im. part
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
   qm2 = qonv(sm2,-1)
   qm3 = qonv(sm3,-1)
   qm4 = qonv(sm4,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   qy1 = qx1/q24
   qy2 = qx2/q24
!
   rslt(0) = li2c2( qy1*qm2 ,qy2*qm2 )/r24*sm2
!
   if (x2.ne.CZRO) then ! better to put a threshold on cc 
     rslt(0) = rslt(0) + ( logc2( qy1/qy2 )*logc( qy1*qy2/(qm2*qm2) ) &
                          -logc2( qx1/qx2 )*logc( qx1*qx2/(qm4*qm4) ) )/(x2*2)
   endif
!
   rslt(0) = rslt(0) - li2c2( qx1*qm4 ,qx2*qm4 )*sm4
!
   if (r23Not0) then
     qss = q23*qm3/q24
     rslt(0) = rslt(0) - li2c2( qx1*qss ,qx2*qss )*r23*sm3/r24
   endif
!
   if (r34Not0) then
     qss = q34*qm3
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34*sm3
   endif
!
   rslt(0) = rslt(0)/(aa*sm2*sm3*sm4)
   end subroutine


   subroutine trif3( rslt ,p1i,p2i,p3i ,m1i,m2i,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1i,p2i,p3i,m1i,m2i,m3i
   type(qmplx_type) :: q12,q13,q23,qm1,qm2,qm3,qx1,qx2,qz1,qz2,qtt
  complex(kindr2) &   
     :: p1,p2,p3,m1,m2,m3,sm1,sm2,sm3,aa,bb,cc,dd,x1,x2 &
                     ,k12,k13,k23,r12,r13,r23,d12,d13,d23 
  real(kindr2) &  
     :: h1,h2,h3
!
   h1 = -aimag(m1i)
   h2 = -aimag(m2i)
   h3 = -aimag(m3i)
   if (h2.ge.h1.and.h2.ge.h3) then
     p1=p3i ;p2=p1i ;p3=p2i ;m1=m3i ;m2=m1i ;m3=m2i
   else
     p1=p1i ;p2=p2i ;p3=p3i ;m1=m1i ;m2=m2i ;m3=m3i
   endif
!
   sm1 = mysqrt(m1)
   sm2 = mysqrt(m2)
   sm3 = mysqrt(m3)
!
   k12 = 0
   k13 = 0
   k23 = 0
   if (m1+m2.ne.p1) k12 = ( m1+m2-p1-p1*IEPS )/(sm1*sm2) ! p1
   if (m1+m3.ne.p3) k13 = ( m1+m3-p3-p3*IEPS )/(sm1*sm3) ! p1+p2 => p12
   if (m2+m3.ne.p2) k23 = ( m2+m3-p2-p2*IEPS )/(sm2*sm3) ! p2
!
   call rfun( r12,d12 ,k12 )
   call rfun( r13,d13 ,k13 )
   call rfun( r23,d23 ,k23 )
!
   aa = sm2/sm3 - k23 + r13*(k12 - sm2/sm1)
!
   if (aa.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   bb = d13/sm2 + k12/sm3 - k23/sm1
   cc = ( sm1/sm3 - 1/r13 )/(sm1*sm2)
!   hh = areal( (r13-sm1/sm3)/(sm1*sm2) )
!   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,0 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1) ! x2 SHOULD HAVE im. part
   q12 = qonv(r12,-1)
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   qm1 = qonv(sm1,-1)
   qm2 = qonv(sm2,-1)
   qm3 = qonv(sm3,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   qz1 = qx1*qm2
   qz2 = qx2*qm2
   rslt(0) = rslt(0) + ( li2c2( qz1*q12 ,qz2*q12 )*r12 &
                        +li2c2( qz1/q12 ,qz2/q12 )/r12 )*sm2
   qtt = q13*qm2
   qz1 = qx1*qtt
   qz2 = qx2*qtt
   rslt(0) = rslt(0) - ( li2c2( qz1*q23 ,qz2*q23 )*r23 &
                        +li2c2( qz1/q23 ,qz2/q23 )/r23 )*r13*sm2
   qz1 = qx1*q13
   qz2 = qx2*q13
   rslt(0) = rslt(0) + li2c2( qz1*qm3 ,qz2*qm3 )*r13*sm3 &
                     - li2c2( qx1*qm1 ,qx2*qm1 )*sm1
   if (x2.ne.CZRO) then
     rslt(0) = rslt(0) + ( logc2( qz1/qz2 )*logc( qz1*qz2/(qm3*qm3) ) &
                          -logc2( qx1/qx2 )*logc( qx1*qx2/(qm1*qm1) ) )/(x2*2)
   endif
!
   rslt(0) = rslt(0)/(aa*sm1*sm2*sm3)
   end subroutine
   

   subroutine trif3HV( rslt ,pp,mm ,ap ,smax ,lam )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! non-zero. Based on the fomula of 't Hooft & Veltman
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: pp(3),mm(3)
  real(kindr2) &  
     ,intent(in)  :: ap(3),smax
  complex(kindr2) &   
     ,optional ,intent(in) :: lam
  complex(kindr2) &   
     :: p1,p2,p3,m1,m2,m3,slam,yy
  complex(kindr2) &   
     :: sm1,sm2,sm3
   type(qmplx_type) :: qm1,qm2,qm3
  real(kindr2) &  
     :: a12,a23,a31,thrs,a1,a2,a3
!
! Order squared momenta, first one smallest
   if     (ap(1).le.ap(2).and.ap(1).le.ap(3)) then
     if (ap(2).le.ap(3)) then
       a1=ap(1) ;a2=ap(2) ;a3=ap(3)
       p1=pp(1) ;p2=pp(2) ;p3=pp(3)
       m1=mm(1) ;m2=mm(2) ;m3=mm(3)
     else
       a1=ap(1) ;a2=ap(3) ;a3=ap(2)
       p1=pp(1) ;p2=pp(3) ;p3=pp(2)
       m1=mm(2) ;m2=mm(1) ;m3=mm(3)
     endif
   elseif (ap(2).le.ap(3).and.ap(2).le.ap(1)) then
     if (ap(3).le.ap(1)) then
       a1=ap(2) ;a2=ap(3) ;a3=ap(1)
       p1=pp(2) ;p2=pp(3) ;p3=pp(1)
       m1=mm(2) ;m2=mm(3) ;m3=mm(1)
     else
       a1=ap(2) ;a2=ap(1) ;a3=ap(3)
       p1=pp(2) ;p2=pp(1) ;p3=pp(3)
       m1=mm(3) ;m2=mm(2) ;m3=mm(1)
     endif
   else
     if (ap(1).le.ap(2)) then
       a1=ap(3) ;a2=ap(1) ;a3=ap(2)
       p1=pp(3) ;p2=pp(1) ;p3=pp(2)
       m1=mm(3) ;m2=mm(1) ;m3=mm(2)
     else
       a1=ap(3) ;a2=ap(2) ;a3=ap(1)
       p1=pp(3) ;p2=pp(2) ;p3=pp(1)
       m1=mm(1) ;m2=mm(3) ;m3=mm(2)
     endif
   endif
!
! Need to cut out negligible squared momenta
   thrs = smax*neglig(prcpar)
!
! Add infinitesimal imaginary parts to masses
   m1 = m1 - abs(areal(m1))*IEPS
   m2 = m2 - abs(areal(m2))*IEPS
   m3 = m3 - abs(areal(m3))*IEPS
!       
   if (a1.gt.thrs) then ! 3 non-zero squared momenta
     if (present(lam)) then ;slam=lam
                       else ;slam=kallen(p1,p2,p3)
     endif
     if (slam.eq.CZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3HV: ' &
         ,'threshold singularity, returning 0'
       rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
       return
     endif
     slam = mysqrt( slam ,1 )
     sm1=mysqrt(m1,-1) ;sm2=mysqrt(m2,-1) ;sm3=mysqrt(m3,-1)
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     rslt(0) = s3fun( p1,sm1,sm2 , (m2-m3)+p2    ,p3-p1-p2 ,p2 ,slam ) &
             - s3fun( p3,sm1,sm3 ,-(m1-m2)+p3-p2 ,p2-p1-p3 ,p1 ,slam ) &
             + s3fun( p2,sm2,sm3 ,-(m1-m2)+p3-p2 ,p1+p2-p3 ,p1 ,slam )
     rslt(0) = -rslt(0)/slam
!
   elseif (a2.gt.thrs) then ! 2 non-zero squared momenta
     if (p2.eq.p3) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3HV: ' &
         ,'threshold singularity, returning 0'
       rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
       return
     endif
     sm1=mysqrt(m1,-1) ;sm2=mysqrt(m2,-1) ;sm3=mysqrt(m3,-1)
     yy = ( (m1-m2)-p3+p2 )/( p2-p3 )
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     rslt(0) = s3fun( p3,sm1,sm3 ,yy ) - s3fun( p2,sm2,sm3 ,yy )
     rslt(0) = rslt(0)/(p2-p3)
!
   elseif (a3.gt.thrs) then ! 1 non-zero squared momentum
     sm1=mysqrt(m1,-1) ;sm3=mysqrt(m3,-1)
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     yy = -( (m1-m2)-p3 )/p3
     rslt(0) = s3fun( p3,sm1,sm3 ,yy ) - s2fun( m2-m3 ,m3 ,yy )
     rslt(0) = -rslt(0)/p3
!
   else ! all squared momenta zero
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     a12=abs(m1-m2) ;a23=abs(m2-m3) ;a31=abs(m3-m1)
     if     (a12.ge.a23.and.a12.ge.a31) then
       if (a12.eq.RZRO) then ;rslt(0)=-1/(2*m3) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm3/qm1) - logc2(qm3/qm2) )/(m1-m2)
       endif
     elseif (a23.ge.a12.and.a23.ge.a31) then
       if (a23.eq.RZRO) then ;rslt(0)=-1/(2*m1) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm1/qm2) - logc2(qm1/qm3) )/(m2-m3)
       endif
     else
       if (a31.eq.RZRO) then ;rslt(0)=-1/(2*m2) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm2/qm3) - logc2(qm2/qm1) )/(m3-m1)
       endif
     endif
   endif
!
   contains
!
     function s3fun( aa,s1,s2 ,t1,t2,t3,t4 ) result(rslt)
!***************************************************************
! int( ( ln(a*y^2+b*y+c) - ln(a*y0^2+b*y0+c) )/(y-y0) ,y=0..1 )
! with  b=s1^2-s2^2-aa  and  c=s2^2
! and with  y0  in terms of t1,t2,t3,t4 defined at the "present"
! function below.
! t4  should be  sqrt(lambda(aa,t2,t3))
!***************************************************************
  complex(kindr2) &   
       ,intent(in) :: aa,s1,s2,t1
  complex(kindr2) &   
       ,optional,intent(in) :: t2,t3
  complex(kindr2) &   
       ,optional,intent(inout) :: t4
  complex(kindr2) &   
       :: rslt ,cc,bb,dd,y0,y1,y2,zz,hh,alpha
  real(kindr2) &  
       :: rez,arez,aimz
     type(qmplx_type) :: q1,q2
!
     bb = (s1+s2)*(s1-s2)-aa
     cc = s2*s2
     dd = (aa-(s1+s2)**2)*(aa-(s1-s2)**2)
     dd = sqrt( dd )!+ sign(abs(dd),areal(aa))*IEPS )
     call solabc( y1,y2 ,dd ,aa,bb,cc ,1 )
!
     if (present(t4)) then
       call solabc( alpha,hh ,t4 ,aa,t2,t3 ,1 )
       y0 = -(t1+bb*alpha)/t4
     else
       y0 = t1
     endif
!
     q1 = qonv(y0-y1)
     q2 = qonv(y0-y2)
     rslt = li2c(qonv(-y1)/q1) - li2c(qonv(1-y1)/q1) &
          + li2c(qonv(-y2)/q2) - li2c(qonv(1-y2)/q2)
! Take some care about the imaginary part of  a*y0^2+b*y0+c=a*(y0-y1)*(y0-y2)
     zz = y0*(aa*y0+bb)
     rez=areal(zz)  ;arez=abs(rez) ;aimz=abs(aimag(zz))
     if (arez*EPSN*EPSN.le.aimz*neglig(prcpar).and.aimz.le.arez*neglig(prcpar)) then
! Here, the value of Imz is just numerical noise due to cancellations.
! Realize that |Imz|~eps^2 indicates there were no such cancellations,
! so the lower limit is needed in in the if-statement!
       zz = (rez + cc)/aa
     else
       zz = (zz + cc)/aa
     endif
     hh = eta3(-y1,-y2,cc/aa) - eta3(y0-y1,y0-y2,zz)
     if (areal(aa).lt.RZRO.and.aimag(zz).lt.RZRO) hh = hh - 2*IPI
     if (hh.ne.CZRO) rslt = rslt + hh*logc(qonv((y0-1)/y0,1))
!
     end function
!
     function s2fun( aa,bb ,y0 ) result(rslt)
!**************************************************
! int( ( ln(a*y+b) - ln(a*y0+b) )/(y-y0) ,y=0..1 )
!**************************************************
  complex(kindr2) &   
       ,intent(in) :: aa,bb,y0
  complex(kindr2) &   
       :: rslt ,y1,hh
     type(qmplx_type) :: q1
     y1 = -bb/aa
     q1 = qonv(y0-y1)
     rslt = li2c(qonv(-y1,-1)/q1) - li2c(qonv(1-y1,-1)/q1)
! aa may have imaginary part, so  theta(-aa)*theta(-Im(y0-y1))  is not
! sufficient and need the following:
     hh = eta5( aa ,-y1,bb ,y0-y1,aa*(y0-y1) )
     if (hh.ne.CZRO) rslt = rslt + hh*logc(qonv((y0-1)/y0,1))
     end function
!
   end subroutine


end module


module avh_olo_dp_box
  use avh_olo_units
  use avh_olo_dp_prec
  use avh_olo_dp_auxfun
  use avh_olo_dp_qmplx
  implicit none
  private
  public :: box00,box03,box05,box06,box07,box08,box09,box10,box11,box12 &
           ,box13,box14,box15,box16,boxf1,boxf2,boxf3,boxf5,boxf4 &
           ,permtable,casetable,base
  integer ,parameter ::  permtable(6,0:15)=reshape((/ &
     1,2,3,4 ,5,6 &! 0, 0 masses non-zero,           no perm
    ,1,2,3,4 ,5,6 &! 1, 1 mass non-zero,             no perm
    ,4,1,2,3 ,6,5 &! 2, 1 mass non-zero,             1 cyclic perm
    ,1,2,3,4 ,5,6 &! 3, 2 neighbour masses non-zero, no perm
    ,3,4,1,2 ,5,6 &! 4, 1 mass   non-zero,           2 cyclic perm's
    ,1,2,3,4 ,5,6 &! 5, 2 opposite masses non-zero,  no perm
    ,4,1,2,3 ,6,5 &! 6, 2 neighbour masses non-zero, 1 cyclic perm
    ,1,2,3,4 ,5,6 &! 7, 3 masses non-zero,           no perm
    ,2,3,4,1 ,6,5 &! 8, 1 mass   non-zero,           3 cyclic perm's
    ,2,3,4,1 ,6,5 &! 9, 2 neighbour masses non-zero, 3 cyclic perm's
    ,4,1,2,3 ,6,5 &!10, 2 opposite masses non-zero,  1 cyclic perm
    ,2,3,4,1 ,6,5 &!11, 3 masses non-zero,           3 cyclic perm's
    ,3,4,1,2 ,5,6 &!12, 2 neighbour masses non-zero, 2 cyclic perm's
    ,3,4,1,2 ,5,6 &!13, 3 masses non-zero,           2 cyclic perm's
    ,4,1,2,3 ,6,5 &!14, 3 masses non-zero,           1 cyclic perm
    ,1,2,3,4 ,5,6 &!15, 4 masses non-zero,           no perm
    /),(/6,16/)) !          0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
  integer ,parameter :: casetable(0:15)= &
                          (/0,1,1,2,1,5,2,3,1,2, 5, 3, 2, 3, 3, 4/)
  integer ,parameter :: base(4)=(/8,4,2,1/)
contains

   subroutine box16( rslt ,p2,p3,p12,p23 ,m2,m3,m4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                     d^(Dim)q
! ------ | ------------------------------------------------------
! i*pi^2 / q^2 [(q+k1)^2-m2] [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=m2, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=m4
! m2,m4 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p12,p23 ,m2,m3,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp12,cp23,cm2,cm3,cm4,sm1,sm2,sm3,sm4 &
                     ,r13,r23,r24,r34,d23,d24,d34,log24,cc
   type(qmplx_type) :: q13,q23,q24,q34,qss,qy1,qy2,qz1,qz2
!
   if (abs(m2).gt.abs(m4)) then
     cm2=m2 ;cm4=m4 ;cp2=p2 ;cp3=p3
   else
     cm2=m4 ;cm4=m2 ;cp2=p3 ;cp3=p2
   endif
   cm3=m3 ;cp12=p12 ;cp23=p23
!
   if (cp12.eq.cm3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box16: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm1 = abs(rmu)
   sm2 = mysqrt(cm2)
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
!
   r13 = (cm3-cp12)/(sm1*sm3)
   call rfun( r23,d23 ,(cm2+cm3-cp2 )/(sm2*sm3) )
   call rfun( r24,d24 ,(cm2+cm4-cp23)/(sm2*sm4) )
   call rfun( r34,d34 ,(cm3+cm4-cp3 )/(sm3*sm4) )
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
!
   if (r24.eq.-CONE) then 
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box16: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   qss = q23*q34
   qy1 = qss*q24
   qy2 = qss/q24
!
   qss = q23/q34
   qz1 = qss*q24
   qz2 = qss/q24
!
   qss = q13*q23
   qss = (qss*qss)/q24
!
   cc = 1/( sm2*sm4*(cp12-cm3) )
   log24 = logc2(q24)*r24/(1+r24)
   rslt(2) = 0
   rslt(1) = -log24
   rslt(0) = log24*logc(qss) + li2c2(q24*q24,qonv(1))*r24 &
           - li2c2(qy1,qy2)*r23*r34 - li2c2(qz1,qz2)*r23/r34
   rslt(1) = cc*rslt(1)
   rslt(0) = cc*rslt(0)
   end subroutine


   subroutine box15( rslt ,p2,p3,p12,p23 ,m2,m4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 [(q+k1)^2-m2] (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=m2, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=m4
! m2,m4 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p12,p23 ,m2,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp12,cp23,cm2,cm4,sm1,sm2,sm3,sm4 &
                     ,r13,r23,r24,r34,d24,log24,cc
   type(qmplx_type) :: q13,q23,q24,q34,qss,qz1,qz2
!
   if (abs(m2-p2).gt.abs(m4-p3)) then
     cm2=m2 ;cm4=m4 ;cp2=p2 ;cp3=p3
   else
     cm2=m4 ;cm4=m2 ;cp2=p3 ;cp3=p2
   endif
   cp12=p12 ;cp23=p23
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box15: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm1 = abs(rmu)
   sm2 = mysqrt(cm2)
   sm4 = mysqrt(cm4)
   sm3 = abs(sm2)
   r13 = (       -cp12)/(sm1*sm3)
   r23 = (cm2    -cp2 )/(sm2*sm3)
   r34 = (    cm4-cp3 )/(sm3*sm4)
   call rfun( r24,d24 ,(cm2+cm4-cp23)/(sm2*sm4) )
!
   if (r24.eq.-CONE) then 
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box15: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
!
   qss = q13/q23
   qss = (qss*qss)/q24
!
   cc = r24/(sm2*sm4*cp12)
   log24 = logc2(q24)/(1+r24)
   rslt(2) = 0
   rslt(1) = -log24
   rslt(0) = log24 * logc(qss) + li2c2(q24*q24,qonv(1))
   if (r34.ne.CZRO) then
     qss = q34/q23
     qz1 = qss*q24
     qz2 = qss/q24
     rslt(0) = rslt(0) - li2c2(qz1,qz2)*r34/(r23*r24)
   endif
   rslt(1) = cc*rslt(1)
   rslt(0) = cc*rslt(0)
   end subroutine


   subroutine box14( rslt ,cp12,cp23 ,cm2,cm4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 [(q+k1)^2-m2] (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=m2, k2^2=m2, k3^2=m4, (k1+k2+k3)^2=m4
! m2,m4 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp12,cp23,cm2,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: sm2,sm4,r24,d24,cc
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box14: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm2 = mysqrt(cm2)
   sm4 = mysqrt(cm4)
   call rfun( r24,d24 ,(cm2+cm4-cp23)/(sm2*sm4) )
!
   if (r24.eq.-CONE) then 
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box14: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   cc = -2*logc2(qonv(r24,-1))*r24/(1+r24)/(sm2*sm4*cp12)
!
   rslt(2) = 0
   rslt(1) = cc
   rslt(0) = -cc*logc(qonv(-cp12/(rmu*rmu),-1))
   end subroutine


   subroutine box13( rslt ,p2,p3,p4,p12,p23 ,m3,m4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=p4
! m3,m4 should NOT be identically 0d0
! p4 should NOT be identical to m4
! p2 should NOT be identical to m3
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p4,p12,p23,m3,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp4,cp12,cp23,cm3,cm4,sm3,sm4,sm1,sm2 &
             ,r13,r14,r23,r24,r34,d34,cc,logd,li2d,loge,li2f,li2b,li2e
   type(qmplx_type) :: q13,q14,q23,q24,q34,qy1,qy2
  real(kindr2) &  
     :: h1,h2
!
   if (p12.eq.m3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box13: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (p23.eq.m4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box13: ' &
       ,'p23=m4, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   h1 = abs((m3-p12)*(m4-p23))
   h2 = abs((m3-p2 )*(m4-p4 ))
   if (h1.ge.h2) then
     cp2=p2  ;cp3=p3 ;cp4=p4  ;cp12=p12 ;cp23=p23 ;cm3=m3 ;cm4=m4
   else
     cp2=p12 ;cp3=p3 ;cp4=p23 ;cp12=p2  ;cp23=p4  ;cm3=m3 ;cm4=m4
   endif
!
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
   sm1 = abs(rmu)
   sm2 = sm1
!
   r13 = (cm3-cp12)/(sm1*sm3)
   r14 = (cm4-cp4 )/(sm1*sm4)
   r23 = (cm3-cp2 )/(sm2*sm3)
   r24 = (cm4-cp23)/(sm2*sm4)
   call rfun( r34,d34 ,(cm3+cm4-cp3)/(sm3*sm4) )
!
   q13 = qonv(r13,-1)
   q14 = qonv(r14,-1)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1) 
!
   qy1 = q14*q23/q13/q24
   logd = logc2(qy1     )/(r13*r24)
   li2d = li2c2(qy1,qonv(1))/(r13*r24)
   loge = logc(q13)
!
   qy1 = q23/q24
   qy2 = q13/q14
   li2f = li2c2( qy1*q34,qy2*q34 )*r34/(r14*r24)
   li2b = li2c2( qy1/q34,qy2/q34 )/(r34*r14*r24)
   li2e = li2c2( q14/q24,q13/q23 )/(r23*r24)
!
   rslt(2) = 0
   rslt(1) = logd
   rslt(0) = li2f + li2b + 2*li2e - 2*li2d - 2*logd*loge
   cc = sm1*sm2*sm3*sm4
   rslt(1) = rslt(1)/cc
   rslt(0) = rslt(0)/cc
   end subroutine


   subroutine box12( rslt ,cp3,cp4,cp12,cp23 ,cm3,cm4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=m3, k3^2=p3, (k1+k2+k3)^2=p4
! m3,m4 should NOT be indentiallcy 0d0
! p4 should NOT be identical to m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cp4,cp12,cp23,cm3,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: sm3,sm4,sm1,sm2,r13,r14,r24,r34,d34,cc &
                     ,log13,log14,log24,log34,li2f,li2b,li2d
   type(qmplx_type) :: q13,q14,q24,q34,qyy
!
   if (cp12.eq.cm3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box12: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box12: ' &
       ,'p23=m4, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
   sm1 = abs(rmu)
   sm2 = sm1
!
   r13 = (cm3-cp12)/(sm1*sm3)
   r14 = (cm4-cp4 )/(sm1*sm4)
   r24 = (cm4-cp23)/(sm2*sm4)
   call rfun( r34,d34 ,(cm3+cm4-cp3)/(sm3*sm4) )
!
   q13 = qonv(r13,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1) 
!
   log13 = logc(q13) 
   log14 = logc(q14) 
   log24 = logc(q24) 
   log34 = logc(q34) 
!
   qyy = q14/q13
   li2f = li2c(qyy*q34)
   li2b = li2c(qyy/q34)
   li2d = li2c(q14/q24)
!
   rslt(2) = 1
   rslt(2) = rslt(2)/2
   rslt(1) = log14 - log24 - log13
   rslt(0) = 2*log13*log24 - log14*log14 - log34*log34 &
           - 2*li2d - li2f - li2b - 3*PISQo24
   cc = (cm3-cp12)*(cm4-cp23) ! = sm1*sm2*sm3*sm4*r13*r24
   rslt(2) = rslt(2)/cc
   rslt(1) = rslt(1)/cc
   rslt(0) = rslt(0)/cc
   end subroutine


   subroutine box11( rslt ,cp3,cp12,cp23 ,cm3,cm4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=m3, k3^2=p3, (k1+k2+k3)^2=m4
! m3,m4 should NOT be indentiallcy 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cp12,cp23,cm3,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: sm3,sm4,sm1,sm2,r13,r24,r34,d34 &
                     ,cc,log13,log24,log34
!
   if (cp12.eq.cm3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box11: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box11: ' &
       ,'p23=m4, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
   sm1 = abs(rmu)
   sm2 = sm1
!
   r13 = (cm3-cp12)/(sm1*sm3)
   r24 = (cm4-cp23)/(sm2*sm4)
   call rfun( r34,d34 ,(cm3+cm4-cp3 )/(sm3*sm4) )
!
   log13 = logc(qonv(r13,-1)) 
   log24 = logc(qonv(r24,-1)) 
   log34 = logc(qonv(r34,-1)) 
!
   rslt(2) = 1
   rslt(1) = -log13-log24
   rslt(0) = 2*log13*log24 - log34*log34 - 14*PISQo24
   cc = (cm3-cp12)*(cm4-cp23) ! = sm1*sm2*sm3*sm4*r13*r24
   rslt(2) = rslt(2)/cc
   rslt(1) = rslt(1)/cc
   rslt(0) = rslt(0)/cc
   end subroutine


   subroutine box10( rslt ,p2,p3,p4,p12,p23 ,m4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=p4
! m4 should NOT be identically 0d0
! p2 should NOT be identically 0d0
! p4 should NOT be identical to m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p4,p12,p23,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp4,cp12,cp23,cm4,r13,r14,r23,r24,r34,z1,z0
   type(qmplx_type) :: q13,q14,q23,q24,q34,qm4,qxx,qx1,qx2
  real(kindr2) &  
     :: h1,h2
!
   if (p12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box10: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (p23.eq.m4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box10: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   h1 = abs(p12*(m4-p23))
   h2 = abs( p2*(m4-p4 ))
   if (h1.ge.h2) then
     cp2=p2  ;cp3=p3 ;cp4=p4  ;cp12=p12 ;cp23=p23 ;cm4=m4
   else
     cp2=p12 ;cp3=p3 ;cp4=p23 ;cp12=p2  ;cp23=p4  ;cm4=m4
   endif
!
   r23 =    -cp2
   r13 =    -cp12
   r34 = cm4-cp3
   r14 = cm4-cp4
   r24 = cm4-cp23
   q23 = qonv(r23,-1)
   q13 = qonv(r13,-1)
   q34 = qonv(r34,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   if (r34.ne.CZRO) then
     qx1 = q34/qm4
     qx2 = qx1*q14/q13
     qx1 = qx1*q24/q23
     z0 = -li2c2(qx1,qx2)*r34/(2*cm4*r23)
   else
     z0 = 0
   endif
!
   qx1 = q23/q13
   qx2 = q24/q14
   qxx = qx1/qx2
   z1 = -logc2(qxx)/r24
   z0 = z0 - li2c2(qx1,qx2)/r14
   z0 = z0 + li2c2(qxx,qonv(1))/r24
   z0 = z0 + z1*( logc(qm4/q24) - logc(qm4/(rmu*rmu))/2 )
!
   rslt(2) = 0
   rslt(1) = -z1/r13
   rslt(0) = -2*z0/r13
   end subroutine


   subroutine box09( rslt ,cp2,cp3,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=m4
! m4 should NOT be identically 0d0
! p2 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp2,cp3,cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: logm,log12,log23,li12,li23,z2,z1,z0,cc &
                     ,r13,r23,r24,r34
   type(qmplx_type) :: q13,q23,q24,q34,qm4,qxx
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box09: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box09: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   r23 =    -cp2
   r13 =    -cp12
   r34 = cm4-cp3
   r24 = cm4-cp23
   q23 = qonv(r23,-1)
   q13 = qonv(r13,-1)
   q34 = qonv(r34,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   logm  = logc(qm4/(rmu*rmu))
   qxx = q13/q23
   log12 = logc(qxx)
   li12  = li2c(qxx)
!
   qxx = q24/qm4
   log23 = logc(qxx)
   li23  = li2c(qxx*q34/q23)
!
   z2 = 1
   z2 = z2/2
   z1 = -log12 - log23
   z0 = li23 + 2*li12 + z1*z1 + PISQo24
   cc = 1/(r13*r24)
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine box08( rslt ,cp3,cp4,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=k2^2=0, k3^2=p3, (k1+k2+k3)^2=p4
! mm should NOT be identically 0d0
! p3 NOR p4 should be identically m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cp4,cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q13,q14,q24,q34,qm4,qxx,qx1,qx2,qx3
  complex(kindr2) &   
     :: r13,r14,r24,r34,z1,z0,cc
  real(kindr2) &  
     :: rmu2
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box08: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box08: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   rmu2 = rmu*rmu
   r13 =    -cp12
   r34 = cm4-cp3
   r14 = cm4-cp4
   r24 = cm4-cp23
   q13 = qonv(r13,-1)
   q34 = qonv(r34,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   qx1 = q34/q24
   qx2 = q14/q24
   qx3 = q13/rmu2
   z1 = logc(qx1*qx2/qx3)
   z0 = 2*( logc(q24/rmu2)*logc(qx3) - (li2c(qx1)+li2c(qx2)) )
!
   qx1 = q34/rmu2
   qx2 = q14/rmu2
   qxx = qx1*qx2/qx3
   z0 = z0 - logc(qx1)**2 - logc(qx2)**2 &
           + logc(qxx)**2/2 + li2c(qm4/qxx/rmu2)
!
   cc = 1/(r13*r24)
   rslt(2) = cc
   rslt(1) = cc*z1
   rslt(0) = cc*( z0 - 6*PISQo24 )
   end subroutine


   subroutine box07( rslt ,cp4,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=k2^2=0, k3^2=m4, (k1+k2+k3)^2=p4
! m3 should NOT be identically 0d0
! p4 should NOT be identically m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp4,cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q13,q14,q24,qm4
  complex(kindr2) &   
     :: r13,r14,r24,logm,log12,log23,log4,li423 &
                     ,z2,z1,z0,cc
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box07: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box07: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   r13 =    -cp12
   r14 = cm4-cp4
   r24 = cm4-cp23
   q13 = qonv(r13,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   logm  = logc(qm4/(rmu*rmu))
   log12 = logc(q13/qm4)
   log23 = logc(q24/qm4)
   log4  = logc(q14/qm4)
   li423 = li2c(q14/q24)
!
   z2 = 3
   z2 = z2/2
   z1 = -2*log23 - log12 + log4
   z0 = 2*(log12*log23 - li423) - log4*log4 - 13*PISQo24
   cc = 1/(r13*r24)
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine box06( rslt ,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=k2^2=0, k3^2=(k1+k2+k3)^2=m4
! m3 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q13,q24,qm4
  complex(kindr2) &   
     :: r13,r24,logm,log1,log2,z2,z1,z0,cc
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box06: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box06: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   r13 =    -cp12
   r24 = cm4-cp23
   q13 = qonv(r13,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   logm = logc(qm4/(rmu*rmu))
   log1 = logc(q13/qm4)
   log2 = logc(q24/qm4)
!
   z2 = 2
   z1 = -2*log2 - log1
   z0 = 2*(log2*log1 - 8*PISQo24)
   cc = 1/(r13*r24)
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine box03( rslt ,p2,p4,p5,p6 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | ---------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 (q+k1+k2+k3)^2
!
! with  k1^2=k3^2=0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p4,p5,p6 
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q2,q4,q5,q6,q26,q54,qy
  complex(kindr2) &   
     :: logy
  real(kindr2) &  
     :: rmu2
!
   rmu2 = rmu*rmu
   q2 = qonv(-p2,-1)
   q4 = qonv(-p4,-1)
   q5 = qonv(-p5,-1)
   q6 = qonv(-p6,-1)
   q26 = q2/q6
   q54 = q5/q4
   qy = q26/q54
   logy = logc2(qy)/(p5*p6)
   rslt(1) = logy
   rslt(0) = li2c2(q6/q4,q2/q5)/(p4*p5) &
           + li2c2(q54,q26)/(p4*p6)     &
           - li2c2(qonv(1),qy)/(p5*p6) &
           - logy*logc(q54*q2*q6/(rmu2*rmu2))/2
   rslt(2) = 0
   rslt(1) = 2*rslt(1)
   rslt(0) = 2*rslt(0)
   end subroutine


   subroutine box05( rslt ,p2,p3,p4,p5,p6 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | ---------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 (q+k1+k2+k3)^2
!
! with  k1^2=0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p4,p5,p6
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) ::q2,q3,q4,q5,q6 ,q25,q64,qy,qz
  complex(kindr2) &   
     :: logy
  real(kindr2) &  
     :: rmu2
!
   rmu2 = rmu*rmu
   q2 = qonv(-p2,-1)
   q3 = qonv(-p3,-1)
   q4 = qonv(-p4,-1)
   q5 = qonv(-p5,-1)
   q6 = qonv(-p6,-1)
   q25 = q2/q5
   q64 = q6/q4
   qy = q25/q64
   qz = q64*q2*q5*q6*q6/q3/q3/(rmu2*rmu2)
!
   logy = logc2(qy)/(p5*p6)
   rslt(2) = 0
   rslt(1) = logy
   rslt(0) = li2c2(q64,q25)/(p4*p5) &
           - li2c2(qonv(1),qy)/(p5*p6) &
           - logy*logc(qz)/4
   rslt(0) = 2*rslt(0)
   end subroutine


   subroutine box00( rslt ,cp ,api ,rmu )
!*******************************************************************
! calculates
!               C   /              d^(Dim)q
!            ------ | ---------------------------------------
!            i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 (q+k1+k2+k3)^2
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  p1 = k1^2,  p2 = k2^2,  p3 = k3^2,  p4 = (k1+k2+k3)^2,
!         p12 = (k1+k2)^2,  p23 = (k2+k3)^2
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! If any of these numbers is IDENTICALLY 0d0, the corresponding
! IR-singular case is returned.
!*******************************************************************
   use avh_olo_dp_olog
   use avh_olo_dp_dilog
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp(6)
  real(kindr2) &  
     ,intent(in)  :: api(6),rmu
  complex(kindr2) &   
     :: log3,log4,log5,log6,li24,li25,li26 &
                     ,li254,li263
  real(kindr2) &  
     :: rp1,rp2,rp3,rp4,rp5,rp6,pp(6),ap(6),gg,ff,hh,arg,rmu2
   integer :: icase,sf,sgn,i3,i4,i5,i6
   integer ,parameter :: base(4)=(/8,4,2,1/)
!
   rmu2 = rmu*rmu
   ff = api(5)*api(6)
   gg = api(2)*api(4)
   hh = api(1)*api(3)
   if     (ff.ge.gg.and.ff.ge.hh) then
     pp(1)=areal(cp(1)) ;ap(1)=api(1)
     pp(2)=areal(cp(2)) ;ap(2)=api(2)
     pp(3)=areal(cp(3)) ;ap(3)=api(3)
     pp(4)=areal(cp(4)) ;ap(4)=api(4)
     pp(5)=areal(cp(5)) ;ap(5)=api(5)
     pp(6)=areal(cp(6)) ;ap(6)=api(6)
   elseif (gg.ge.ff.and.gg.ge.hh) then
     pp(1)=areal(cp(1)) ;ap(1)=api(1)
     pp(2)=areal(cp(6)) ;ap(2)=api(6)
     pp(3)=areal(cp(3)) ;ap(3)=api(3)
     pp(4)=areal(cp(5)) ;ap(4)=api(5)
     pp(5)=areal(cp(4)) ;ap(5)=api(4)
     pp(6)=areal(cp(2)) ;ap(6)=api(2)
   else
     pp(1)=areal(cp(5)) ;ap(1)=api(5)
     pp(2)=areal(cp(2)) ;ap(2)=api(2)
     pp(3)=areal(cp(6)) ;ap(3)=api(6)
     pp(4)=areal(cp(4)) ;ap(4)=api(4)
     pp(5)=areal(cp(1)) ;ap(5)=api(1)
     pp(6)=areal(cp(3)) ;ap(6)=api(3)
   endif
!
   icase = 0
   if (ap(1).gt.RZRO) icase = icase + base(1)
   if (ap(2).gt.RZRO) icase = icase + base(2)
   if (ap(3).gt.RZRO) icase = icase + base(3)
   if (ap(4).gt.RZRO) icase = icase + base(4)
   rp1 = pp(permtable(1,icase))
   rp2 = pp(permtable(2,icase))
   rp3 = pp(permtable(3,icase))
   rp4 = pp(permtable(4,icase))
   rp5 = pp(permtable(5,icase))
   rp6 = pp(permtable(6,icase))
   icase = casetable(   icase)
!
   i3=0 ;if (-rp3.lt.RZRO) i3=-1
   i4=0 ;if (-rp4.lt.RZRO) i4=-1
   i5=0 ;if (-rp5.lt.RZRO) i5=-1
   i6=0 ;if (-rp6.lt.RZRO) i6=-1
!
   if     (icase.eq.0) then
! 0 masses non-zero
     gg = 1/( rp5 * rp6 )
     log5 = olog(abs(rp5/rmu2),i5)
     log6 = olog(abs(rp6/rmu2),i6)
     rslt(2) = gg*( 4 )
     rslt(1) = gg*(-2*(log5 + log6) )
     rslt(0) = gg*( log5**2 + log6**2 - olog(abs(rp5/rp6),i5-i6)**2 - 32*PISQo24 )
   elseif (icase.eq.1) then
! 1 mass non-zero
     gg = 1/( rp5 * rp6 )
     ff =  gg*( rp5 + rp6 - rp4 )
     log4 = olog(abs(rp4/rmu2),i4)
     log5 = olog(abs(rp5/rmu2),i5)
     log6 = olog(abs(rp6/rmu2),i6)
     sf = sgnRe(ff)
     sgn = 0
       arg = rp4*ff 
       if (arg.lt.RZRO) sgn = sf
       li24 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp5*ff 
       if (arg.lt.RZRO) sgn = sf
       li25 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp6*ff 
       if (arg.lt.RZRO) sgn = sf
       li26 = dilog(abs(arg),sgn)
     rslt(2) = gg*( 2 )
     rslt(1) = gg*( 2*(log4-log5-log6) )
     rslt(0) = gg*( log5**2 + log6**2 - log4**2 - 12*PISQo24 &
                   + 2*(li25 + li26 - li24) )
   elseif (icase.eq.2) then
! 2 neighbour masses non-zero
     gg = 1/( rp5 * rp6 )
     ff =  gg*( rp5 + rp6 - rp4 )
     log3 = olog(abs(rp3/rmu2),i3)
     log4 = olog(abs(rp4/rmu2),i4)
     log5 = olog(abs(rp5/rmu2),i5)
     log6 = olog(abs(rp6/rmu2),i6)
     li254 = dilog( abs(rp4/rp5) ,i4-i5 )
     li263 = dilog( abs(rp3/rp6) ,i3-i6 )
     sf = sgnRe(ff)
     sgn = 0
       arg = rp4*ff 
       if (arg.lt.RZRO) sgn = sf
       li24 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp5*ff 
       if (arg.lt.RZRO) sgn = sf
       li25 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp6*ff 
       if (arg.lt.RZRO) sgn = sf
       li26 = dilog(abs(arg),sgn)
     rslt(2) = gg
     rslt(1) = gg*( log4 + log3 - log5 - 2*log6 )
     rslt(0) = gg*( log5**2 + log6**2 - log3**2 - log4**2 &
                   + (log3 + log4 - log5)**2/2 &
                   - 2*PISQo24 + 2*(li254 - li263 + li25 + li26 - li24) )
   elseif (icase.eq.5) then
! 2 opposite masses non-zero
     call box03( rslt ,acmplx(rp2),acmplx(rp4) &
                      ,acmplx(rp5),acmplx(rp6) ,rmu )
   elseif (icase.eq.3) then
! 3 masses non-zero
     call box05( rslt ,acmplx(rp2),acmplx(rp3) &
                      ,acmplx(rp4),acmplx(rp5) &
                      ,acmplx(rp6) ,rmu )
   elseif (icase.eq.4) then
! 4 masses non-zero
     call boxf0( rslt ,acmplx(rp1),acmplx(rp2) &
                      ,acmplx(rp3),acmplx(rp4) &
                      ,acmplx(rp5),acmplx(rp6) )
   endif
   end subroutine

  
  subroutine boxf0( rslt ,p1,p2,p3,p4,p12,p23 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with all internal masses
! equal zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23
  type(qmplx_type) :: q12,q13,q14,q23,q24,q34,qx1,qx2,qss
  complex(kindr2) &   
    :: aa,bb,cc,dd,x1,x2,ss,r12,r13,r14,r23,r24,r34
  real(kindr2) &  
    :: hh
!
  r12 = -p1  !  p1
  r13 = -p12 !  p1+p2
  r14 = -p4  !  p1+p2+p3
  r23 = -p2  !  p2
  r24 = -p23 !  p2+p3
  r34 = -p3  !  p3      
!
  aa = r34*r24
!
  if (r13.eq.CZRO.or.aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf0: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = r13*r24 + r12*r34 - r14*r23
  cc = r12*r13
  hh = areal(r23)
  dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
  call solabc(x1,x2,dd ,aa,bb,cc ,1)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 , hh)
  qx2 = qonv(x2 ,-hh)
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  qss = q34/q13
  rslt(0) = rslt(0) + li2c2(qx1*qss,qx2*qss) * r34/r13
!
  qss = q24/q12
  rslt(0) = rslt(0) + li2c2(qx1*qss,qx2*qss) * r24/r12
!
  ss = -logc2(qx1/qx2) / x2
  rslt(0) = rslt(0) + ss*( logc(qx1*qx2)/2 - logc(q12*q13/q14/q23) )
!
  rslt(0) = -rslt(0) / aa
  end subroutine


  subroutine boxf1( rslt ,p1,p2,p3,p4,p12,p23 ,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with one internal mass
! non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23 ,m4
  type(qmplx_type) :: qx1,qx2,qss,q12,q13,q14,q23,q24,q34
  complex(kindr2) &   
    :: smm,sm4,aa,bb,cc,dd,x1,x2,r12,r13,r14,r23,r24,r34
  logical :: r12zero,r13zero,r14zero
!
  sm4 = mysqrt(m4)
  smm = abs(sm4) 
!
  r12 = ( m4-p4 -p4 *IEPS )/(smm*sm4)
  r13 = ( m4-p23-p23*IEPS )/(smm*sm4)
  r14 = ( m4-p3 -p3 *IEPS )/(smm*sm4)
  r23 = (   -p1 -p1 *IEPS )/(smm*smm)
  r24 = (   -p12-p12*IEPS )/(smm*smm)
  r34 = (   -p2 -p2 *IEPS )/(smm*smm)
!
  r12zero=(abs(areal(r12))+abs(aimag(r12)).lt.neglig(prcpar))
  r13zero=(abs(areal(r13))+abs(aimag(r13)).lt.neglig(prcpar))
  r14zero=(abs(areal(r14))+abs(aimag(r14)).lt.neglig(prcpar))
!
  aa = r34*r24
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf1: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = r13*r24 + r12*r34 - r14*r23
  cc = r12*r13 - r23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 ,1 )
  qx2 = qonv(x2 ,1 )
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  if (r12zero.and.r13zero) then
    qss = qx1*qx2*q34*q24/q23
    qss = qss*qss
    rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
  else
    if (r13zero) then
      qss = q34*q12/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    else
      qss = q34/q13
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34/r13
    endif
    if (r12zero) then
      qss = q24*q13/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    else
      qss = q24/q12
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r24/r12
    endif
    if (.not.r12zero.and..not.r13zero) then
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( q12*q13/q23 )/x2
    endif
  endif
!
  if (.not.r14zero) then
    rslt(0) = rslt(0) - li2c2( qx1*q14 ,qx2*q14 )*r14
  endif
!
  rslt(0) = -rslt(0)/(aa*smm*smm*smm*sm4)
  end subroutine


  subroutine boxf5( rslt ,p1,p2,p3,p4,p12,p23, m2,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with two opposite internal
! masses non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m2,m4
  call boxf2( rslt ,p12,p2,p23,p4,p1,p3 ,m2,m4 )
  end subroutine


  subroutine boxf2( rslt ,p1,p2,p3,p4,p12,p23 ,m3,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with two adjacent internal
! masses non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m3,m4
  type(qmplx_type) :: qx1,qx2,qss,q12,q13,q14,q23,q24,q34
  complex(kindr2) &   
    :: smm,sm3,sm4,aa,bb,cc,dd,x1,x2 &
                    ,r12,r13,r14,r23,r24,r34,d14,k14
  logical :: r12zero,r13zero,r24zero,r34zero
!
  sm3 = mysqrt(m3)
  sm4 = mysqrt(m4)
!
  smm = abs(sm3)
!
  r12 = (    m4-p4 -p4 *IEPS )/(smm*sm4)
  r13 = (    m4-p23-p23*IEPS )/(smm*sm4)
  k14 = ( m3+m4-p3 -p3 *IEPS )/(sm3*sm4)
  r23 = (      -p1 -p1 *IEPS )/(smm*smm)
  r24 = (    m3-p12-p12*IEPS )/(smm*sm3)
  r34 = (    m3-p2 -p2 *IEPS )/(smm*sm3)
!
  r12zero = (abs(areal(r12))+abs(aimag(r12)).lt.neglig(prcpar))
  r13zero = (abs(areal(r13))+abs(aimag(r13)).lt.neglig(prcpar))
  r24zero = (abs(areal(r24))+abs(aimag(r24)).lt.neglig(prcpar))
  r34zero = (abs(areal(r34))+abs(aimag(r34)).lt.neglig(prcpar))
!
  if (r12zero.and.r24zero) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf2: ' &
       ,'m4=p4 and m3=p12, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
  if (r13zero.and.r34zero) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf2: ' &
       ,'m4=p23 and m3=p2, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  call rfun( r14,d14 ,k14 )
!
  aa = r34*r24 - r23
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf2: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = r13*r24 + r12*r34 - k14*r23
  cc = r12*r13 - r23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 ,1 )
  qx2 = qonv(x2 ,1 )
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  rslt(0) = rslt(0) - li2c2( qx1*q14 ,qx2*q14 )*r14
  rslt(0) = rslt(0) - li2c2( qx1/q14 ,qx2/q14 )/r14
!
  if (r12zero.and.r13zero) then
    qss = qx1*qx2*q34*q24/q23
    qss = qss*qss
    rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
  else
    if (r13zero) then
      qss = q34*q12/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    elseif (.not.r34zero) then
      qss = q34/q13
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34/r13
    endif
    if (r12zero) then
      qss = q24*q13/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    elseif (.not.r24zero) then
      qss = q24/q12
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r24/r12
    endif
    if (.not.r12zero.and..not.r13zero) then
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( q12*q13/q23 )/x2 
    endif
  endif
!
  rslt(0) = -rslt(0)/(aa*smm*smm*sm3*sm4)
  end subroutine


  subroutine boxf3( rslt ,pp ,mm )
!*******************************************************************
! Finite 1-loop scalar 4-point function with three internal masses
! non-zero.
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: pp(6),mm(4)
  integer :: j
  integer ,parameter :: ip(6)=(/4,5,2,6,3,1/)
  integer ,parameter :: im(4)=(/4,1,3,2/)
  integer ,parameter :: ic(4,6)=reshape((/1,2,3,4 ,2,3,4,1 ,3,4,1,2 &
                                  ,4,1,2,3 ,5,6,5,6 ,6,5,6,5/),(/4,6/))
!
  if     (mm(1).eq.CZRO) then ;j=3
  elseif (mm(2).eq.CZRO) then ;j=4
  elseif (mm(3).eq.CZRO) then ;j=1
  else                        ;j=2
  endif
  call boxf33( rslt ,pp(ic(j,ip(1))) ,pp(ic(j,ip(2))) ,pp(ic(j,ip(3))) &
                    ,pp(ic(j,ip(4))) ,pp(ic(j,ip(5))) ,pp(ic(j,ip(6))) &
                    ,mm(ic(j,im(1))) ,mm(ic(j,im(2))) ,mm(ic(j,im(4))) )
  end subroutine

  subroutine boxf33( rslt ,p1,p2,p3,p4,p12,p23, m1,m2,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with three internal masses
! non-zero, and m3=0. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m1,m2,m4
  type(qmplx_type) :: qx1,qx2,qss,q12,q13,q14,q23,q24,q34,qy1,qy2
  complex(kindr2) &   
    :: sm1,sm2,sm3,sm4 ,aa,bb,cc,dd,x1,x2 &
                    ,r12,r13,r14,r23,r24,r34,d12,d14,d24,k12,k14,k24
  logical ::r13zero,r23zero,r34zero
!
  sm1 = mysqrt(m1)
  sm2 = mysqrt(m2)
  sm4 = mysqrt(m4)
  sm3 = abs(sm2)
!
  k12 = ( m1+m2-p1 -p1 *IEPS )/(sm1*sm2) ! p1
  r13 = ( m1   -p12-p12*IEPS )/(sm1*sm3) ! p1+p2
  k14 = ( m1+m4-p4 -p4 *IEPS )/(sm1*sm4) ! p1+p2+p3
  r23 = ( m2   -p2 -p2 *IEPS )/(sm2*sm3) ! p2
  k24 = ( m2+m4-p23-p23*IEPS )/(sm2*sm4) ! p2+p3
  r34 = (    m4-p3 -p3 *IEPS )/(sm3*sm4) ! p3
!
  r13zero = (abs(areal(r13))+abs(aimag(r13)).lt.neglig(prcpar))
  r23zero = (abs(areal(r23))+abs(aimag(r23)).lt.neglig(prcpar))
  r34zero = (abs(areal(r34))+abs(aimag(r34)).lt.neglig(prcpar))
!
  if (r13zero) then
    if     (r23zero) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf33: ' &
       ,'m4=p4 and m3=p12, returning 0'
      rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
      return
    elseif (r34zero) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf33: ' &
       ,'m2=p1 and m3=p12, returning 0'
      rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
      return
    endif
  endif
!
  call rfun( r12,d12 ,k12 )
  call rfun( r14,d14 ,k14 )
  call rfun( r24,d24 ,k24 )
!
  aa = r34/r24 - r23
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf33: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = -r13*d24 + k12*r34 - k14*r23
  cc = k12*r13 + r24*r34 - k14*r24*r13 - r23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 ,1 ) ! x1 SHOULD HAVE im. part
  qx2 = qonv(x2 ,1 ) ! x2 SHOULD HAVE im. part
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  qy1 = qx1/q24
  qy2 = qx2/q24
  rslt(0) = rslt(0) + li2c2( qy1*q12 ,qy2*q12 )/r24*r12
  rslt(0) = rslt(0) + li2c2( qy1/q12 ,qy2/q12 )/r24/r12
  rslt(0) = rslt(0) - li2c2( qx1*q14 ,qx2*q14 )*r14
  rslt(0) = rslt(0) - li2c2( qx1/q14 ,qx2/q14 )/r14
!
  if (.not.r13zero) then
    if (.not.r23zero) then
      qss = q23/q13/q24
      rslt(0) = rslt(0) - li2c2( qx1*qss ,qx2*qss )*r23/(r13*r24)
    endif
    if (.not.r34zero) then
      qss = q34/q13
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34/r13
    endif
  else
    rslt(0) = rslt(0) - logc2( qx1/qx2 )*logc( q23/q24/q34 )/x2 
  endif
!
  rslt(0) = -rslt(0)/(aa*sm1*sm2*sm3*sm4)
  end subroutine


  subroutine boxf4( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with all internal masses
! non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m1,m2,m3,m4
  type(qmplx_type) :: q12,q13,q14,q23,q24,q34,qx1,qx2,qy1,qy2,qtt
  complex(kindr2) &   
    :: sm1,sm2,sm3,sm4 ,aa,bb,cc,dd,x1,x2,tt &
                    ,k12,k13,k14,k23,k24,k34 &
                    ,r12,r13,r14,r23,r24,r34 &
                    ,d12,d13,d14,d23,d24,d34
  real(kindr2) &  
    :: h1,h2
!
  sm1 = mysqrt(m1)
  sm2 = mysqrt(m2)
  sm3 = mysqrt(m3)
  sm4 = mysqrt(m4)
!
  k12 = ( m1+m2-p1 -p1 *IEPS)/(sm1*sm2) ! p1
  k13 = ( m1+m3-p12-p12*IEPS)/(sm1*sm3) ! p1+p2
  k14 = ( m1+m4-p4 -p4 *IEPS)/(sm1*sm4) ! p1+p2+p3
  k23 = ( m2+m3-p2 -p2 *IEPS)/(sm2*sm3) ! p2
  k24 = ( m2+m4-p23-p23*IEPS)/(sm2*sm4) ! p2+p3
  k34 = ( m3+m4-p3 -p3 *IEPS)/(sm3*sm4) ! p3
!
  call rfun( r12,d12 ,k12 )
  call rfun( r13,d13 ,k13 )
  call rfun( r14,d14 ,k14 )
  call rfun( r23,d23 ,k23 )
  call rfun( r24,d24 ,k24 )
  call rfun( r34,d34 ,k34 )
!
  aa = k34/r24 + r13*k12 - k14*r13/r24 - k23
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf4: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = d13*d24 + k12*k34 - k14*k23
  cc = k12/r13 + r24*k34 - k14*r24/r13 - k23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
!
  h1 = areal(k23 - r13*k12 - r24*k34 + r13*r24*k14)
  h2 = h1*areal(aa)*areal(x1)
  h1 = h1*areal(aa)*areal(x2)
!
  qx1 = qonv(-x1,-h1) ! x1 should have im. part
  qx2 = qonv(-x2,-h2) ! x2 should have im. part
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  qy1 = qx1/q24
  qy2 = qx2/q24
  rslt(0) = rslt(0) + ( li2c2( qy1*q12 ,qy2*q12 )*r12 &
                      + li2c2( qy1/q12 ,qy2/q12 )/r12 )/r24
  tt = r13/r24
  qtt = qonv(tt,-areal(r24) )
  qy1 = qx1*qtt
  qy2 = qx2*qtt
  rslt(0) = rslt(0) - ( li2c2( qy1*q23 ,qy2*q23 )*r23 &
                      + li2c2( qy1/q23 ,qy2/q23 )/r23 )*tt
  qy1 = qx1*q13
  qy2 = qx2*q13
  rslt(0) = rslt(0) + ( li2c2( qy1*q34 ,qy2*q34 )*r34 &
                      + li2c2( qy1/q34 ,qy2/q34 )/r34 )*r13
!
  rslt(0) = rslt(0) - ( li2c2( qx1*q14 ,qx2*q14 )*r14 &
                      + li2c2( qx1/q14 ,qx2/q14 )/r14 )
!
  rslt(0) = -rslt(0)/(aa*sm1*sm2*sm3*sm4)
  end subroutine

end module


module avh_olo_dp_boxc
   use avh_olo_units
   use avh_olo_dp_prec
   use avh_olo_dp_auxfun
   use avh_olo_dp_qmplx
   implicit none
   private
   public :: boxc

contains

   subroutine boxc( rslt ,pp_in ,mm_in ,ap_in ,smax )
!*******************************************************************
! Finite 1-loop scalar 4-point function for complex internal masses
! Based on the formulas from
!   Dao Thi Nhung and Le Duc Ninh, arXiv:0902.0325 [hep-ph]
!   G. 't Hooft and M.J.G. Veltman, Nucl.Phys.B153:365-401,1979 
!*******************************************************************
   use avh_olo_dp_box ,only: base,casetable,ll=>permtable
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: pp_in(6),mm_in(4)
  real(kindr2) &  
     ,intent(in)  :: ap_in(6),smax
  complex(kindr2) &   
     :: pp(6),mm(4)
  real(kindr2) &  
     :: ap(6),aptmp(6),rem,imm,hh
  complex(kindr2) &   
     :: a,b,c,d,e,f,g,h,j,k,dpe,epk,x1,x2,sdnt,o1,j1,e1 &
       ,dek,dpf,def,dpk,abc,bgj,jph,cph
   integer :: icase,jcase,ii
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   hh = neglig(prcpar)*smax
   do ii=1,6
     if (ap_in(ii).ge.hh) then ;ap(ii)=ap_in(ii)
                          else ;ap(ii)=0
     endif
   enddo
!
   do ii=1,4
     if (ap(ii).eq.RZRO) then ;pp(ii)=0
                         else ;pp(ii)=pp_in(ii)
     endif
   enddo
   if (ap(5).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
       ,' |s| too small, putting it by hand'
     ap(5) = hh
     pp(5) = acmplx(sign(hh,areal(pp_in(5))))
   else
     pp(5) = pp_in(5)
   endif
   if (ap(6).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
       ,' |t| too small, putting it by hand'
     ap(6) = hh
     pp(6) = acmplx(sign(hh,areal(pp_in(6))))
   else
     pp(6) = pp_in(6)
   endif
!
   do ii=1,4
     rem = areal(mm_in(ii))
     imm = aimag(mm_in(ii))
     hh = EPSN*abs(rem)
     if (abs(imm).lt.hh) imm = -hh
     mm(ii) = acmplx(rem,imm)
   enddo
!
   icase = 0
   do ii=1,4
     if (ap(ii).gt.RZRO) icase = icase + base(ii)
   enddo
!
   if (icase.lt.15) then
! at least one exernal mass equal zero
     jcase = casetable(icase)
     if (jcase.eq.0.or.jcase.eq.1.or.jcase.eq.5) then
! two opposite masses equal zero
       a = pp(ll(5,icase)) - pp(ll(1,icase))
       c = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(3,icase))
       g = pp(ll(2,icase))
       h = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
       d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
       e = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
       f = mm(ll(4,icase))
       k = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
       dpe = (mm(ll(1,icase)) - mm(ll(4,icase))) - pp(ll(4,icase))
       dpk = (mm(ll(2,icase)) - mm(ll(4,icase))) - pp(ll(6,icase))
       dpf = mm(ll(3,icase)) - pp(ll(3,icase))
       rslt(0) = t13fun( a,c,g,h ,d,e,f,k ,dpe,dpk,dpf )
     else
       a = pp(ll(3,icase))
       b = pp(ll(2,icase))
       c = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
       h = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(6,icase)) + pp(ll(2,icase))
       j = pp(ll(5,icase)) - pp(ll(1,icase)) - pp(ll(2,icase))
       d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
       e = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
       k = (mm(ll(1,icase)) - mm(ll(2,icase))) + pp(ll(6,icase)) - pp(ll(4,icase))
       f = mm(ll(4,icase))
       cph = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(3,icase))
       dpe = (mm(ll(2,icase)) - mm(ll(4,icase))) - pp(ll(6,icase))
       epk = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
       dek = (mm(ll(1,icase)) - mm(ll(4,icase))) - pp(ll(4,icase))
       dpf = mm(ll(3,icase)) - pp(ll(3,icase))
       rslt(0) = tfun( a,b  ,c  ,h,j ,d,e  ,f ,k ,dpe,dpf ) &
               - tfun( a,b+j,cph,h,j ,d,epk,f ,k ,dek,dpf )
     endif
   else
! no extenal mass equal zero
     if    (areal((pp(5)-pp(1)-pp(2))**2-4*pp(1)*pp(2)).gt.RZRO)then ;icase=0 !12, no permutation
     elseif(areal((pp(6)-pp(2)-pp(3))**2-4*pp(2)*pp(3)).gt.RZRO)then ;icase=8 !23, 1 cyclic permutation
     elseif(areal((pp(4)-pp(5)-pp(3))**2-4*pp(5)*pp(3)).gt.RZRO)then ;icase=4 !34, 2 cyclic permutations
     elseif(areal((pp(4)-pp(1)-pp(6))**2-4*pp(1)*pp(6)).gt.RZRO)then ;icase=2 !41, 3 cyclic permutations
     else
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
         ,'no positive lambda, returning 0'
       return
     endif
     a = pp(ll(3,icase))
     b = pp(ll(2,icase))
     g = pp(ll(1,icase))
     c = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
     h = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(6,icase)) + pp(ll(2,icase))
     j = pp(ll(5,icase)) - pp(ll(1,icase)) - pp(ll(2,icase))
     d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
     e = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
     k = (mm(ll(1,icase)) - mm(ll(2,icase))) + pp(ll(6,icase)) - pp(ll(4,icase))
     f = mm(ll(4,icase))
     abc = pp(ll(6,icase))
     bgj = pp(ll(5,icase))
     jph = pp(ll(4,icase)) - pp(ll(1,icase)) - pp(ll(6,icase))
     cph = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(3,icase))
     dpe = (mm(ll(2,icase)) - mm(ll(4,icase))) - pp(ll(6,icase))
     epk = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
     dek = (mm(ll(1,icase)) - mm(ll(4,icase))) - pp(ll(4,icase))
     dpf = mm(ll(3,icase)) - pp(ll(3,icase))
     def = mm(ll(2,icase)) - pp(ll(6,icase))
     call solabc( x1,x2 ,sdnt ,g,j,b ,0 )
     if (aimag(sdnt).ne.RZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
         ,'no real solution for alpha, returning 0'
       return
     endif
!BAD        if (abs(areal(x1)).gt.abs(areal(x2))) then
     if (abs(areal(x1)).lt.abs(areal(x2))) then !BETTER
       sdnt = x1
       x1 = x2
       x2 = sdnt
     endif
     o1 = 1-x1
     j1 = j+2*g*x1
     e1 = e+k*x1
     rslt(0) =   -tfun( abc,g  ,jph,c+2*b+(h+j)*x1, j1   ,dpe,k  ,f,e1 ,dek,def ) &
             + o1*tfun( a  ,bgj,cph,c+h*x1        , o1*j1,d  ,epk,f,e1 ,dek,dpf ) &
             + x1*tfun( a  ,b  ,c  ,c+h*x1        ,-j1*x1,d  ,e  ,f,e1 ,dpe,dpf )
   endif
   end subroutine


   function t13fun( aa,cc,gg,hh ,dd,ee,ff,jj ,dpe,dpj,dpf ) result(rslt)
!*******************************************************************
! /1   /x                             y
! | dx |  dy -----------------------------------------------------
! /0   /0    (gy^2 + hxy + dx + jy + f)*(ax^2 + cxy + dx + ey + f)
!
! jj should have negative imaginary part
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: aa,cc,gg,hh ,dd,ee,ff,jj ,dpe,dpj,dpf
  complex(kindr2) &   
     :: rslt ,kk,ll,nn,y1,y2,sdnt
!
!
   kk = hh*aa - cc*gg
   ll = aa*dd + hh*ee - dd*gg - cc*jj
   nn = dd*(ee - jj) + (hh - cc)*(ff-IEPS*abs(areal(ff)))
   call solabc( y1,y2 ,sdnt ,kk,ll,nn ,0 )
!
   rslt = - s3fun( y1,y2 ,CZRO,CONE ,aa   ,ee+cc,dpf ) &
          + s3fun( y1,y2 ,CZRO,CONE ,gg   ,jj+hh,dpf ) &
          - s3fun( y1,y2 ,CZRO,CONE ,gg+hh,dpj  ,ff  ) &
          + s3fun( y1,y2 ,CZRO,CONE ,aa+cc,dpe  ,ff  )
!
   rslt = rslt/kk
   end function


   function t1fun( aa,cc,gg,hh ,dd,ee,ff,jj ,dpe ) result(rslt)
!*******************************************************************
! /1   /x                         1
! | dx |  dy ----------------------------------------------
! /0   /0    (g*x + h*x + j)*(a*x^2 + c*xy + d*x + e*y + f)
!
! jj should have negative imaginary part
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: aa,cc,gg,hh ,dd,ee,ff,jj,dpe
  complex(kindr2) &   
     ::rslt ,kk,ll,nn,y1,y2,sdnt
!
!
   kk = hh*aa - cc*gg
   ll = hh*dd - cc*jj - ee*gg
   nn = hh*(ff-IEPS*abs(areal(ff))) - ee*jj
   call solabc( y1,y2 ,sdnt ,kk,ll,nn ,0 )
!
   rslt = - s3fun( y1,y2 ,CZRO,CONE ,aa+cc,dpe  ,ff ) &
          + s3fun( y1,y2 ,CZRO,CONE ,CZRO ,gg+hh,jj ) &
          - s3fun( y1,y2 ,CZRO,CONE ,CZRO ,gg   ,jj ) &
          + s3fun( y1,y2 ,CZRO,CONE ,aa   ,dd   ,ff )
!
   rslt = rslt/kk
   end function


   function tfun( aa,bb,cc ,gin,hin ,dd,ee,ff ,jin ,dpe ,dpf ) result(rslt)
!*******************************************************************
! /1   /x                             1
! | dx |  dy ------------------------------------------------------
! /0   /0    (g*x + h*x + j)*(a*x^2 + b*y^2 + c*xy + d*x + e*y + f)
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: aa,bb,cc ,gin,hin ,dd,ee,ff ,jin ,dpe ,dpf
  complex(kindr2) &   
     :: rslt ,gg,hh,jj,zz(2),beta,tmpa(2),tmpb(2) &
       ,tmpc(2),kiz(2),ll,nn,kk,y1,y2,yy(2,2),sdnt
  real(kindr2) &  
     :: ab1,ab2,ac1,ac2,abab,acac,abac,det,ap1,ap2 &
                  ,apab,apac,x1(2,2),x2(2,2),xmin
   integer :: iz,iy,izmin,sj
   logical :: pp(2,2),p1,p2
!
!
   sj = sgnIm(jin,-1)
   gg = -sj*gin
   hh = -sj*hin
   jj = -sj*jin
!
   if     (bb.eq.CZRO) then
     rslt = -sj*t1fun( aa,cc,gg,hh ,dd,ee,ff,jj ,dpe )
     return
   elseif (aa.eq.CZRO) then
     rslt = -sj*t1fun( bb+cc,-cc,-gg-hh,gg, -dpe-2*(bb+cc),dd+cc &
                      ,dpe+bb+cc+ff,gg+hh+jj ,-ee-2*bb-cc )
     return
   endif
!
   call solabc( zz(1),zz(2) ,sdnt ,bb,cc,aa ,0 )
   if (abs(zz(1)).gt.abs(zz(2))) then
     beta = zz(1)
     zz(1) = zz(2)
     zz(2) = beta
   endif
!
   do iz=1,2
     beta = zz(iz)
     tmpa(iz) = gg + beta*hh
     tmpb(iz) = cc + 2*beta*bb
     tmpc(iz) = dd + beta*ee
     kiz(iz) =        bb*tmpa(iz)               - hh*tmpb(iz)
     ll      =        ee*tmpa(iz) - hh*tmpc(iz) - jj*tmpb(iz)
     nn      = (ff-IEPS*abs(areal(ff)))*tmpa(iz) - jj*tmpc(iz)
     call solabc( yy(iz,1),yy(iz,2) ,sdnt ,kiz(iz),ll,nn ,0 )
     if (abs(aimag(beta)).ne.RZRO) then
       ab1 = areal(-beta)
       ab2 = aimag(-beta)
       ac1 = ab1+1 !areal(1-beta)
       ac2 = ab2   !aimag(1-beta)
       abab = ab1*ab1 + ab2*ab2
       acac = ac1*ac1 + ac2*ac2
       abac = ab1*ac1 + ab2*ac2
       det = abab*acac - abac*abac
       do iy=1,2
         ap1 = areal(yy(iz,iy))
         ap2 = aimag(yy(iz,iy))
         apab = ap1*ab1 + ap2*ab2
         apac = ap1*ac1 + ap2*ac2
         x1(iz,iy) = ( acac*apab - abac*apac )/det
         x2(iz,iy) = (-abac*apab + abab*apac )/det
       enddo
     else
       do iy=1,2
         x1(iz,iy) = -1
         x2(iz,iy) = -1
       enddo
     endif
   enddo
   xmin = 1
   izmin = 2
   do iz=1,2
   do iy=1,2
     if ( x1(iz,iy).ge.RZRO.and.x2(iz,iy).ge.RZRO &
                 .and.x1(iz,iy)+x2(iz,iy).le.RONE ) then
       pp(iz,iy) = .true.
       if (x1(iz,iy).lt.xmin) then
         xmin = x1(iz,iy)
         izmin = iz
       endif
       if (x2(iz,iy).lt.xmin) then
         xmin = x2(iz,iy)
         izmin = iz
       endif
     else
       pp(iz,iy) = .false.
     endif
   enddo
   enddo
   iz = izmin+1
   if (iz.eq.3) iz = 1
!
   beta = zz(iz)
   kk = kiz(iz)
   y1 = yy(iz,1)
   y2 = yy(iz,2)
   p1 = pp(iz,1)
   p2 = pp(iz,2)
!
   rslt =   s3fun( y1,y2 ,beta ,CONE      ,CZRO    ,hh   ,gg+jj  ) &
          - s3fun( y1,y2 ,CZRO ,CONE-beta ,CZRO    ,gg+hh,   jj  ) &
          + s3fun( y1,y2 ,CZRO ,    -beta ,CZRO    ,gg   ,   jj  ) &
          - s3fun( y1,y2 ,beta ,CONE      ,bb      ,cc+ee,aa+dpf ) &
          + s3fun( y1,y2 ,CZRO ,CONE-beta ,aa+bb+cc,dpe  ,ff     ) &
          - s3fun( y1,y2 ,CZRO ,    -beta ,aa      ,dd   ,ff     )
!
   sdnt = plnr( y1,y2 ,p1,p2, tmpa(iz),tmpb(iz),tmpc(iz) )
   if (aimag(beta).le.RZRO) then ;rslt = rslt + sdnt
                            else ;rslt = rslt - sdnt
   endif
!
   rslt = -sj*rslt/kk
   end function


   function s3fun( y1i,y2i ,dd,ee ,aa,bb,cin ) result(rslt)
!*******************************************************************
! Calculate
!            ( S3(y1i) - S3(y2i) )/( y1i - y2i )
! where
!               /1    ee * ln( aa*x^2 + bb*x + cc )
!       S3(y) = |  dx -----------------------------
!               /0           ee*x - y - dd
!
! y1i,y2i should have a non-zero imaginary part
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) ::  y1i,y2i ,dd,ee ,aa,bb,cin
  complex(kindr2) &   
     :: rslt ,y1,y2,fy1y2,z1,z2,tmp,cc
  real(kindr2) &  
     ::rea,reb,rez1,rez2,imz1,imz2,simc,hh
!
!
   if (ee.eq.CZRO) then
     rslt = 0
     return
   endif
!
   cc = cin
   rea = abs(aa)
   reb = abs(bb)
   simc = abs(cc)
   if (simc.lt.10*neglig(prcpar)*min(rea,reb)) cc = 0
!
   simc = aimag(cc)
   if (simc.eq.RZRO) then
     simc = aimag(bb)
     if (simc.eq.RZRO) simc = -1
   endif
   simc = sgnRe(simc)
!
   y1 = (dd+y1i)/ee
   y2 = (dd+y2i)/ee
   if (aimag(y1).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'y1 has zero imaginary part'
   endif
   if (aimag(y2).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'y2 has zero imaginary part'
   endif
   fy1y2 = r0fun( y1,y2 )
!
   if     (aa.ne.CZRO) then
!
!     call solabc( z1,z2 ,tmp ,aa,bb,cc ,0 )
     call solabc_rcc( z1,z2 ,areal(aa),bb,cc )
     rea  = sgnRe(aa)
     rez1 = areal(z1)
     rez2 = areal(z2) 
     imz1 = aimag(z1) ! sign(Im(a*z1*z2)) = simc
     imz2 = aimag(z2)
     hh = abs(EPSN2*rez1)
!     if (abs(imz1).lt.EPSN*hh) imz1 = simc*rea*sgnRe(rez2)*hh
     if (imz1.eq.RZRO) imz1 = simc*rea*sgnRe(rez2)*hh
     hh = abs(EPSN2*rez2)
!     if (abs(imz2).lt.EPSN*hh) imz2 = simc*rea*sgnRe(rez1)*hh
     if (imz2.eq.RZRO) imz2 = simc*rea*sgnRe(rez1)*hh
     z1 = acmplx( rez1,imz1)
     z2 = acmplx( rez2,imz2)
     rslt = fy1y2 * ( logc(qonv(aa,simc)) &
                    + eta3( -z1,-imz1,-z2,-imz2,CZRO,simc*rea ) ) &
          + r1fun( z1,y1,y2,fy1y2 ) &
          + r1fun( z2,y1,y2,fy1y2 )
!
   elseif (bb.ne.CZRO) then
!
     z1 = -cc/bb ! - i|eps|Re(b)
     reb  = areal(bb)
     rez1 = areal(z1)
     imz1 = aimag(z1)
     if (abs(imz1).eq.RZRO) then
       imz1 = -simc*reb*abs(EPSN2*rez1/reb)
       z1 = acmplx( rez1,imz1)
     endif
     rslt = fy1y2 * ( logc(qonv(bb,simc)) &
                    + eta3(bb,simc ,-z1,-imz1 ,cc,simc) ) &
          + r1fun( z1,y1,y2,fy1y2 )
!
   elseif (cc.ne.CZRO) then
!
     rslt = logc( qonv(cc,simc) )*fy1y2
!
   else!if (aa=bb=cc=0)
!
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'cc equal zero, returning 0'
     rslt = 0
!
   endif
!
   rslt = rslt/ee
   end function


   function r1fun( zz,y1,y2,fy1y2 ) result(rslt)
!*******************************************************************
! calculates  ( R1(y1,z) - R1(y2,z) )/( y1 - y2 )
! where
!                          /     / 1-y \       / 1-z \ \
!      R1(y,z) = ln(y-z) * | log |-----| - log |-----| |
!                          \     \ -y  /       \ -z  / / 
!
!                      /    y-z \       /    y-z \
!                - Li2 |1 - ----| + Li2 |1 - ----|
!                      \    -z  /       \    1-z /
!
!                                     / 1-y1 \       / 1-y2 \
!                                 log |------| - log |------| 
! input fy1y2 should be equal to      \  -y1 /       \  -y2 /
!                                 ---------------------------
!                                           y1 - y2
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: y1,y2,zz,fy1y2
  complex(kindr2) &   
     :: rslt ,oz
   type(qmplx_type) :: q1z,q2z,qq
  real(kindr2) &  
     :: h12,hz1,hz2,hzz,hoz
   logical :: zzsmall,ozsmall
!
!
   oz = 1-zz
   h12 = abs(y1-y2)
   hz1 = abs(y1-zz)
   hz2 = abs(y2-zz)
   hzz = abs(zz)
   hoz = abs(oz)
   q1z = qonv(y1-zz)
   q2z = qonv(y2-zz)
!
   zzsmall = .false.
   ozsmall = .false.
   if     (hzz.lt.hz1.and.hzz.lt.hz2.and.hzz.lt.hoz) then ! |z| < |y1-z|,|y2-z|
     zzsmall = .true.
     rslt = fy1y2*logc( q1z ) &
          - ( logc(q1z*q2z)/2 + logc(qonv((y2-1)/y2)) &
                                     - logc(qonv(oz)) )*logc2(q1z/q2z)/(y2-zz)
   elseif (hoz.lt.hz1.and.hoz.lt.hz2) then ! |1-z| < |y1-z|,|y2-z|
     ozsmall = .true.
     rslt = fy1y2*logc( q1z ) &
          - (-logc(q1z*q2z)/2 + logc(qonv((y2-1)/y2)) &
                                    + logc(qonv(-zz)) )*logc2(q1z/q2z)/(y2-zz)
   elseif (h12.le.hz2.and.hz2.le.hz1) then ! |y1-y2| < |y2-z| < |y1-z|
     rslt = fy1y2*logc( q1z ) - r0fun( y2,zz )*logc2( q1z/q2z )        
   elseif (h12.le.hz1.and.hz1.le.hz2) then ! |y1-y2| < |y2-z| < |y1-z|
     rslt = fy1y2*logc( q2z ) - r0fun( y1,zz )*logc2( q2z/q1z )        
   else!if(hz1.lt.h12.or.hz2.lt.h12) then ! |y2-z|,|y1-z| < |y1-y2|
     rslt = 0
     if (hz1.ne.RZRO) rslt = rslt + (y1-zz)*logc( q1z )*r0fun( y1,zz )
     if (hz2.ne.RZRO) rslt = rslt - (y2-zz)*logc( q2z )*r0fun( y2,zz )
     rslt = rslt/(y1-y2)
   endif
!
   if (zzsmall) then ! |z| < |y1-z|,|y2-z|
     qq  = qonv(-zz)
     rslt = rslt + ( li2c( qq/q1z ) - li2c( qq/q2z ) )/(y1-y2)
   else
     qq  = qonv(-zz)
     rslt = rslt + li2c2( q1z/qq ,q2z/qq )/zz
   endif
!
   if (ozsmall) then ! |1-z| < |y1-z|,|y2-z|
     qq  = qonv(oz)
     rslt = rslt - ( li2c( qq/q1z ) - li2c( qq/q2z ) )/(y1-y2)
   else
     qq = qonv(oz)
     rslt = rslt + li2c2( q1z/qq ,q2z/qq )/oz
   endif
   end function


   function r0fun( y1,y2 ) result(rslt)
!*******************************************************************
!      / 1-y1 \       / 1-y2 \
!  log |------| - log |------| 
!      \  -y1 /       \  -y2 /
!  ---------------------------
!            y1 - y2
!
! y1,y2 should have non-zero imaginary parts
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: y1,y2
  complex(kindr2) &   
     :: rslt ,oy1,oy2
!
   oy1 = 1-y1
   oy2 = 1-y2
   rslt = logc2( qonv(-y2)/qonv(-y1) )/y1 &
        + logc2( qonv(oy2)/qonv(oy1) )/oy1
   end function


   function plnr( y1,y2 ,p1,p2 ,aa,bb,cc ) result(rslt)
!*******************************************************************
!                   /   a    \          /   a    \
!            p1*log |--------| - p2*log |--------| 
!                   \ b*y1+c /          \ b*y2+c /
! 2*pi*imag* -------------------------------------
!                           y1 - y2
! 
! p1,p2 are logical, to be interpreted as 0,1 in the formula above 
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: y1,y2 ,aa,bb,cc
   logical         ,intent(in) :: p1,p2
  complex(kindr2) &   
     :: rslt ,x1,x2,xx
   type(qmplx_type) :: q1,q2
!
   if (p1) then
     x1 = bb*y1 + cc
     xx = aa/x1
     if (aimag(xx).eq.RZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop plnr: ' &
         ,'aa/x1 has zero imaginary part'
     endif
     q1 = qonv(xx)
   endif
   if (p2) then
     x2 = bb*y2 + cc
     xx = aa/x2
     if (aimag(xx).eq.RZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop plnr: ' &
         ,'aa/x2 has zero imaginary part'
     endif
     q2 = qonv(xx)
   endif
   if (p1) then
     if (p2) then
       rslt = logc2( q2/q1 ) * 2*IPI*bb/x2
     else
       rslt = logc( q1 ) * 2*IPI/(y1-y2)
     endif
   elseif (p2) then
     rslt = logc( q2 ) * 2*IPI/(y2-y1) ! minus sign
   else
     rslt = 0
   endif
   end function


end module


module avh_olo_dp
  use avh_olo_units
  use avh_olo_dp_print
  use avh_olo_dp_prec
!
  implicit none
  private
  public :: olo_unit ,olo_scale ,olo_onshell ,olo_setting
  public :: olo_precision
  public :: olo_a0 ,olo_b0 ,olo_db0 ,olo_b11 ,olo_c0 ,olo_d0
  public :: olo_an ,olo_bn
  public :: olo
  public :: olo_get_scale ,olo_get_onshell ,olo_get_precision
!
  integer ,public ,parameter :: olo_kind=kindr2    
!
  real(kindr2) &  
         ,save :: onshellthrs
  logical,save :: nonzerothrs = .false.
!
  real(kindr2) &  
         ,save :: muscale
!
  character(99) ,parameter :: warnonshell=&
       'it seems you forgot to put some input explicitly on shell. ' &
     //'You may  call olo_onshell  to cure this.'
!
  logical ,save :: initz=.true.
!
  interface olo_a0
    module procedure a0_r,a0rr,a0_c,a0cr
  end interface 
  interface olo_an
    module procedure an_r,anrr,an_c,ancr
  end interface 
  interface olo_b0
    module procedure b0rr,b0rrr,b0rc,b0rcr,b0cc,b0ccr
  end interface 
  interface olo_db0
    module procedure db0rr,db0rrr,db0rc,db0rcr,db0cc,db0ccr
  end interface 
  interface olo_b11
    module procedure b11rr,b11rrr,b11rc,b11rcr,b11cc,b11ccr
  end interface 
  interface olo_bn
    module procedure bnrr,bnrrr,bnrc,bnrcr,bncc,bnccr
  end interface 
  interface olo_c0
    module procedure c0rr,c0rrr,c0rc,c0rcr,c0cc,c0ccr
  end interface 
  interface olo_d0
    module procedure d0rr,d0rrr,d0rc,d0rcr,d0cc,d0ccr
  end interface 
!
  interface olo
    module procedure a0_r,a0rr,a0_c,a0cr
    module procedure an_r,anrr,an_c,ancr
    module procedure b0rr,b0rrr,b0rc,b0rcr,b0cc,b0ccr
    module procedure b11rr,b11rrr,b11rc,b11rcr,b11cc,b11ccr
    module procedure bnrr,bnrrr,bnrc,bnrcr,bncc,bnccr
    module procedure c0rr,c0rrr,c0rc,c0rcr,c0cc,c0ccr
    module procedure d0rr,d0rrr,d0rc,d0rcr,d0cc,d0ccr
  end interface 

contains

 
  subroutine init( ndec )
!*******************************************************************
!*******************************************************************
  use avh_olo_version
  integer,optional,intent(in) :: ndec
!
  call olo_version
!
  initz = .false.
!
  if (present(ndec)) then
    call olo_precision( ndec )
  else
    call olo_precision( 15 )
  endif
!
  onshellthrs = 0
  muscale = 1
  if (.not.nonzerothrs) onshellthrs = neglig(prcpar)
!
  end subroutine
 
 
  recursive subroutine olo_precision( ndec )
!*******************************************************************
!*******************************************************************
  use avh_olo_dp_olog  ,only: update_olog
  use avh_olo_dp_dilog ,only: update_dilog
  use avh_olo_dp_bnlog ,only: update_bnlog
  integer ,intent(in) :: ndec
  logical :: newprc
  if (initz) then
    call init( ndec )
  else
    call set_precision( newprc )       
    if (newprc) then
      call update_olog
      call update_dilog
      call update_bnlog
    endif
    if (.not.nonzerothrs) onshellthrs = neglig(prcpar)
  endif
  end subroutine

 
  subroutine olo_unit( val ,message )
!*******************************************************************
!*******************************************************************
  integer     ,intent(in) :: val
  character(*),intent(in),optional :: message
  if (initz) call init
  if (present(message)) then ;call set_unit( message ,val )
  else                       ;call set_unit( 'all'   ,val )
  endif
  end subroutine
 
 
  subroutine olo_scale( val )
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) ,intent(in) :: val
  if (initz) call init
  muscale = convert(val)
  end subroutine
 
 
  subroutine olo_onshell( thrs )
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) ,intent(in) :: thrs
  if (initz) call init
  nonzerothrs = .true.
  onshellthrs = convert(thrs)
  end subroutine


  function olo_get_precision() result(rslt)
!*******************************************************************
!*******************************************************************
  use avh_olo_dp_prec ,only: ndecim,prcpar
  integer :: rslt
  if (initz) call init
  rslt = ndecim(prcpar)
  end function

  function olo_get_scale() result(rslt)
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) :: rslt
  if (initz) call init
  rslt = adble(muscale)
  end function

  function olo_get_onshell() result(rslt)
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) :: rslt
  if (initz) call init
  rslt = adble(onshellthrs)
  end function


  subroutine olo_setting( iunit )
!*******************************************************************
!*******************************************************************
  integer,optional,intent(in) :: iunit
  integer :: nunit
  if (initz) call init
  nunit = munit
  if (present(iunit)) nunit = iunit
  if (nunit.le.0) return
!
  write(nunit,*) 'MESSAGE from OneLOop: real kind parameter =',trim(myprint(kindr2)) 
  write(nunit,*) 'MESSAGE from OneLOop: number of decimals  =',trim(myprint(ndecim(prcpar)))
!
  if (nonzerothrs) then
    write(nunit,*) 'MESSAGE from OneLOop: on-shell threshold =',trim(myprint(onshellthrs,12))
  else
    write(nunit,*) 'MESSAGE from OneLOop: on-shell threshold is not set'
  endif
!
  write(nunit,*) 'MESSAGE from OneLOop: default scale (mu, not mu^2) =',trim(myprint(muscale,12))
!
  end subroutine
 
 
!*******************************************************************
!
!           C   / d^(Dim)q
! rslt = ------ | -------- 
!        i*pi^2 / (q^2-mm)
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  mm = mass squared
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine a0_c( rslt ,mm )
!
  use avh_olo_dp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: mm
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine a0cr( rslt ,mm ,rmu )
!
  use avh_olo_dp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine a0_r( rslt ,mm  )
!
  use avh_olo_dp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: mm
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine a0rr( rslt ,mm ,rmu )
!
  use avh_olo_dp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine


  subroutine an_c( rslt ,rank ,mm )
!
  use avh_olo_dp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: mm
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine

  subroutine ancr( rslt ,rank ,mm ,rmu )
!
  use avh_olo_dp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine

  subroutine an_r( rslt ,rank ,mm  )
!
  use avh_olo_dp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: mm
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine

  subroutine anrr( rslt ,rank ,mm ,rmu )
!
  use avh_olo_dp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine


!*******************************************************************
!
!           C   /      d^(Dim)q
! rslt = ------ | --------------------
!        i*pi^2 / [q^2-m1][(q+k)^2-m2]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  pp = k^2, m1,m2 = mass squared
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine b0cc( rslt ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0ccr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rc( rslt ,pp ,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rcr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rr( rslt ,pp ,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rrr( rslt ,pp ,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

!*******************************************************************
! Derivative of B0
!*******************************************************************
  subroutine db0cc( rslt ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: dbub0
  use avh_olo_dp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0ccr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: dbub0
  use avh_olo_dp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rc( rslt ,pp ,m1,m2 )
!
  use avh_olo_dp_bub ,only: dbub0
  use avh_olo_dp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rcr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: dbub0
  use avh_olo_dp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rr( rslt ,pp ,m1,m2 )
!
  use avh_olo_dp_bub ,only: dbub0
  use avh_olo_dp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rrr( rslt ,pp ,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: dbub0
  use avh_olo_dp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine



!*******************************************************************
! Return the Papparino-Veltman functions b11,b00,b1,b0 , for
!
!      C   /      d^(Dim)q
!   ------ | -------------------- = b0
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
!      C   /    d^(Dim)q q^mu
!   ------ | -------------------- = p^mu b1
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
!      C   /  d^(Dim)q q^mu q^nu
!   ------ | -------------------- = g^{mu,nu} b00 + p^mu p^nu b11
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine b11cc( b11,b00,b1,b0 ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11ccr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rc( b11,b00,b1,b0 ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rcr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rr( b11,b00,b1,b0 ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rrr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine


  subroutine bncc( rslt ,rank ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnccr( rslt ,rank ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrc( rslt ,rank ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrcr( rslt ,rank ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrr( rslt ,rank ,pp,m1,m2 )
!
  use avh_olo_dp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrrr( rslt ,rank ,pp,m1,m2 ,rmu )
!
  use avh_olo_dp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine


!*******************************************************************
! calculates
!               C   /               d^(Dim)q
!            ------ | ---------------------------------------
!            i*pi^2 / [q^2-m1] [(q+k1)^2-m2] [(q+k1+k2)^2-m3]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps)
!             * GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
!
! input:  p1=k1^2, p2=k2^2, p3=(k1+k2)^2,  m1,m2,m3=squared masses
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine c0cc( rslt ,p1,p2,p3 ,m1,m2,m3 )
  use avh_olo_dp_tri
  use avh_olo_dp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
!
  complex(kindr2) &   
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0ccr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
  use avh_olo_dp_tri
  use avh_olo_dp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  complex(kindr2) &   
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rc( rslt ,p1,p2,p3 ,m1,m2,m3 )
  use avh_olo_dp_tri
  use avh_olo_dp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
!
  real(kindr2) &  
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rcr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
  use avh_olo_dp_tri
  use avh_olo_dp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rr( rslt ,p1,p2,p3 ,m1,m2,m3 )
  use avh_olo_dp_tri
  use avh_olo_dp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3
!
  real(kindr2) &  
    :: pp(3)
  real(kindr2) &  
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rrr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
  use avh_olo_dp_tri
  use avh_olo_dp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(3)
  real(kindr2) &  
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine


!*******************************************************************
! calculates
!
!    C   /                      d^(Dim)q
! ------ | --------------------------------------------------------
! i*pi^2 / [q^2-m1][(q+k1)^2-m2][(q+k1+k2)^2-m3][(q+k1+k2+k3)^2-m4]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps)
!             * GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
!
! input:  p1=k1^2, p2=k2^2, p3=k3^2, p4=(k1+k2+k3)^2, 
!         p12=(k1+k2)^2, p23=(k2+k3)^2, 
!         m1,m2,m3,m4=squared masses
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  avh_olo_dp_onshell  to find out how this
! routines decides when to return IR-divergent cases.
!*******************************************************************

  subroutine d0cc( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
  use avh_olo_dp_box
  use avh_olo_dp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
!
  complex(kindr2) &   
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) ,0 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0ccr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
  use avh_olo_dp_box
  use avh_olo_dp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  complex(kindr2) &   
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) ,0 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rc( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
  use avh_olo_dp_box
  use avh_olo_dp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
!
  real(kindr2) &  
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rcr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
  use avh_olo_dp_box
  use avh_olo_dp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
  use avh_olo_dp_box
  use avh_olo_dp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3,m4
!
  real(kindr2) &  
    :: pp(6)
  real(kindr2) &  
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rrr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
  use avh_olo_dp_box
  use avh_olo_dp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3,m4
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(6)
  real(kindr2) &  
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

end module


module avh_olo_qp_kinds
  integer ,parameter :: kindr2=16 
end module


module avh_olo_qp_arrays
  use avh_olo_units
  use avh_olo_qp_kinds 
  implicit none
  private
  public :: shift1,shift2,shift3,resize,enlarge

! Increase the size of the last dimension by one,
! and move  x(...,n:nsize)  to  x(...,n+1:nsize+1).
  interface shift1 ! for x(:)
    module procedure shift1_r,shift1_i
  end interface
  interface shift2 ! for x(:,:)
    module procedure shift2_r,shift2_i
  end interface
  interface shift3 ! for x(:,:,:)
    module procedure shift3_r,shift3_i
  end interface

! Resize x to the new bounds. Anything that doesn't fit anymore is lost.
  interface resize
    module procedure resize1_r,resize2_r
  end interface

! Resize x to the maximum of the bounds it has and then new bounds.
  interface enlarge
    module procedure enlarge1_r,enlarge2_r
  end interface

contains

  subroutine shift1_r( xx ,nn )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:)
  integer        ,intent(in   ) :: nn
  real(kindr2) &  
    ,allocatable :: tt(:)
  integer ,parameter :: dm=1
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift1_r'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(dm):ub(dm)))
  xx(lb(dm):nn-1) = tt(lb(dm):nn-1)
  xx(nn+1:ub(dm)) = tt(nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift1_i( xx ,nn )
  integer ,allocatable ,intent(inout) :: xx(:)
  integer              ,intent(in   ) :: nn
  integer ,allocatable :: tt(:)
  integer ,parameter :: dm=1
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift1_i'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(dm):ub(dm)))
  xx(lb(dm):nn-1) = tt(lb(dm):nn-1)
  xx(nn+1:ub(dm)) = tt(nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift2_r( xx ,nn )
  real(kindr2) &  
          ,allocatable ,intent(inout) :: xx(:,:)
  integer              ,intent(in   ) :: nn
  real(kindr2) &  
          ,allocatable :: tt(:,:)
  integer ,parameter :: dm=2
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift2_r'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(dm):ub(dm)))
  xx(:,lb(dm):nn-1) = tt(:,lb(dm):nn-1)
  xx(:,nn+1:ub(dm)) = tt(:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift2_i( xx ,nn )
  integer ,allocatable ,intent(inout) :: xx(:,:)
  integer              ,intent(in   ) :: nn
  integer ,allocatable :: tt(:,:)
  integer ,parameter :: dm=2
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift2_i'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(dm):ub(dm)))
  xx(:,lb(dm):nn-1) = tt(:,lb(dm):nn-1)
  xx(:,nn+1:ub(dm)) = tt(:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift3_r( xx ,nn )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:,:,:)
  integer        ,intent(in   ) :: nn
  real(kindr2) &  
    ,allocatable :: tt(:,:,:)
  integer ,parameter :: dm=3
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift3_r'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  xx(:,:,lb(dm):nn-1) = tt(:,:,lb(dm):nn-1)
  xx(:,:,nn+1:ub(dm)) = tt(:,:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

  subroutine shift3_i( xx ,nn )
  integer ,allocatable ,intent(inout) :: xx(:,:,:)
  integer              ,intent(in   ) :: nn
  integer ,allocatable :: tt(:,:,:)
  integer ,parameter :: dm=3
  integer :: lb(dm),ub(dm)
  if (.not.allocated(xx)) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop shift3_i'
    stop
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  tt = xx
  deallocate(xx)
  ub(dm) = ub(dm)+1
  allocate(xx(lb(1):ub(1),lb(2):ub(2),lb(dm):ub(dm)))
  xx(:,:,lb(dm):nn-1) = tt(:,:,lb(dm):nn-1)
  xx(:,:,nn+1:ub(dm)) = tt(:,:,nn:ub(dm)-1)
  deallocate(tt)
  end subroutine

 
  subroutine resize1_r( xx ,l1,u1 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:)
  integer        ,intent(in   ) :: l1,u1
  real(kindr2) &  
    ,allocatable :: tt(:)
  integer :: lb(1),ub(1)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1)))
  tt = xx
  deallocate(xx)
  allocate( xx(l1:u1) )
  lb(1)=max(l1,lb(1)) ;ub(1)=min(u1,ub(1))
  xx(lb(1):ub(1)) = tt(lb(1):ub(1))
  deallocate(tt)
  end subroutine 

  subroutine resize2_r( xx ,l1,u1 ,l2,u2 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:,:)
  integer        ,intent(in   ) :: l1,u1,l2,u2
  real(kindr2) &  
    ,allocatable :: tt(:,:)
  integer :: lb(2),ub(2)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1,l2:u2))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  allocate(tt(lb(1):ub(1),lb(2):ub(2)))
  tt = xx
  deallocate(xx)
  allocate( xx(l1:u1,l2:u2) )
  lb(1)=max(l1,lb(1)) ;ub(1)=min(u1,ub(1))
  lb(2)=max(l2,lb(2)) ;ub(2)=min(u2,ub(2))
  xx(lb(1):ub(1),lb(2):ub(2)) = &
  tt(lb(1):ub(1),lb(2):ub(2))
  deallocate(tt)
  end subroutine 


  subroutine enlarge1_r( xx ,l1,u1 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:)
  integer        ,intent(in   ) :: l1,u1
  real(kindr2) &  
    ,allocatable :: tt(:)
  integer :: lb(1),ub(1)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  if (lb(1).le.l1.and.u1.le.ub(1)) return
  if (lb(1).gt.ub(1)) then
    deallocate( xx )
    allocate( xx(min(l1,lb(1)):max(u1,ub(1))) )
    return
  endif
  allocate(tt(lb(1):ub(1)))
  tt = xx
  deallocate(xx)
  allocate( xx(min(l1,lb(1)):max(u1,ub(1))) )
  xx(lb(1):ub(1)) = tt(lb(1):ub(1))
  deallocate(tt)
  end subroutine 

  subroutine enlarge2_r( xx ,l1,u1 ,l2,u2 )
  real(kindr2) &  
    ,allocatable ,intent(inout) :: xx(:,:)
  integer        ,intent(in   ) :: l1,u1,l2,u2
  real(kindr2) &  
    ,allocatable :: tt(:,:)
  integer :: lb(2),ub(2)
  if (.not.allocated(xx)) then
    allocate(xx(l1:u1,l2:u2))
    return
  endif
  lb=lbound(xx) ;ub=ubound(xx)
  if (lb(1).le.l1.and.u1.le.ub(1).and. &
      lb(2).le.l2.and.u2.le.ub(2)      ) return
  if (lb(1).gt.ub(1).or.lb(2).gt.ub(2)) then
    deallocate( xx )
    allocate( xx(min(l1,lb(1)):max(u1,ub(1))  &
                ,min(l2,lb(2)):max(u2,ub(2))) )
    return
  endif
  allocate(tt(lb(1):ub(1),lb(2):ub(2)))
  tt = xx
  deallocate(xx)
  allocate( xx(min(l1,lb(1)):max(u1,ub(1))  &
              ,min(l2,lb(2)):max(u2,ub(2))) )
  xx(lb(1):ub(1),lb(2):ub(2)) = &
  tt(lb(1):ub(1),lb(2):ub(2))
  deallocate(tt)
  end subroutine 

end module


module avh_olo_qp_prec
  use avh_olo_qp_kinds

  implicit none
  public
  private :: IMAG,acmplx_r,acmplx_rr,acmplx_ir,acmplx_ri,acmplx_c

  integer ,save :: prcpar=0
  integer ,save :: ndecim(1)
  real(kindr2) &
          ,save :: epsilo(1),neglig(1)

  real(kindr2) &
    ,save :: RZRO ,RONE ,EPSN ,EPSN2 ,TWOPI ,ONEPI
  complex(kindr2) &
    ,save :: IEPS ,CZRO ,CONE ,IMAG ,PISQo24 ,IPI

  interface acmplx
    module procedure acmplx_r,acmplx_rr,acmplx_ir,acmplx_ri,acmplx_c
  end interface

contains


  subroutine set_precision( newprc )
!***********************************************************************
!***********************************************************************
  use avh_olo_units
  logical ,intent(out) :: newprc
  integer :: ndec                                  
  if (prcpar.eq.1) then                    
    newprc = .false.                             
    return                                       
  endif
  prcpar = 1                                   
  call set_epsn
  newprc = .true.                              
  RZRO=0
  RONE=1
  IMAG=cmplx(0,1,kind=kind(IMAG))
  CZRO=RZRO
  CONE=RONE
  ONEPI=4*atan(RONE)
  TWOPI=2*ONEPI
  PISQo24=CONE*ONEPI*ONEPI/24
  IPI=IMAG*ONEPI
  EPSN2= EPSN*EPSN
  IEPS= EPSN2*IMAG
!
  contains
!
  subroutine set_epsn
  EPSN = epsilon(EPSN)                         
  ndec = -log10(EPSN)                            
  ndecim(prcpar) = ndec                          
  epsilo(prcpar) = EPSN                        
  neglig(prcpar) = EPSN*10**(ndec/7)       
  end subroutine
!
  end subroutine


  function adble(xx) result(rslt)
!***********************************************************************
! Turn real(kindr2) into kind(1d0)
!***********************************************************************
  real(kindr2) ,intent(in) :: xx
  real(kind(1d0)) :: rslt
  rslt = real(xx,kind=kind(rslt))
  end function

  function convert(xx) result(rslt)
!***********************************************************************
! Turn kind(1d0) into real(kindr2)
!***********************************************************************
  real(kind(1d0)) ,intent(in) :: xx
  real(kindr2) :: rslt
  rslt = real(xx,kind=kind(rslt))
  end function

  function areal(zz) result(rslt)
!***********************************************************************
! Get real part of a complex
!***********************************************************************
  complex(kindr2) &
    ,intent(in) :: zz
  real(kindr2) &
    :: rslt
  rslt = zz
  end function

  function acmplx_r(xx) result(rslt)
!***********************************************************************
! Turn a real into a complex
!***********************************************************************
  real(kindr2) &
    ,intent(in) :: xx
  complex(kindr2) &
    :: rslt
  rslt = xx
  end function
  
  function acmplx_rr(xx,yy) result(rslt)
!***********************************************************************
! Turn two reals into one complex
!***********************************************************************
  real(kindr2) &
    ,intent(in) :: xx,yy
  complex(kindr2) &
    :: rslt
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
  function acmplx_ri(xx,yy) result(rslt)
!***********************************************************************
! Turn a real and an integer into one complex
!***********************************************************************
  real(kindr2) &
           ,intent(in) :: xx
  integer  ,intent(in) :: yy
  complex(kindr2) &
    :: rslt
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
  function acmplx_ir(xx,yy) result(rslt)
!***********************************************************************
! Turn an integer and a real into one complex
!***********************************************************************
  integer ,intent(in) :: xx
  real(kindr2) &
          ,intent(in) :: yy
  complex(kindr2) &
    :: rslt
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
  function acmplx_c(zz) result(rslt)
!***********************************************************************
! Replaces the real part of zz by its absolute value
!***********************************************************************
  complex(kindr2) &
    ,intent(in) :: zz
  complex(kindr2) &
    :: rslt
  real(kindr2) &
    :: xx,yy
  xx = zz
  xx = abs(xx)
  yy = aimag(zz)
  rslt = cmplx(xx,yy,kind=kind(rslt))
  end function
  
end module


module avh_olo_qp_print
  use avh_olo_qp_prec
  implicit none
  private
  public :: myprint

  integer ,parameter :: novh=10 !maximally 6 decimals for exponent
  integer ,parameter :: nxtr=4  !extra decimals

  interface myprint
    module procedure printr,printc,printi
  end interface

contains

  function printc( zz ,ndec ) result(rslt)
  complex(kindr2) &   
    ,intent(in) :: zz
  integer,optional,intent(in) :: ndec
  character((ndecim(prcpar)+nxtr+novh)*2+3) :: rslt
  if (present(ndec)) then
    rslt = '('//trim(printr(areal(zz),ndec)) &
         //','//trim(printr(aimag(zz),ndec)) &
         //')'
  else
    rslt = '('//trim(printr(areal(zz))) &
         //','//trim(printr(aimag(zz))) &
         //')'
  endif
  rslt = adjustl(rslt)
  end function

  function printr( xx_in ,ndec_in ) result(rslt)
  real(kindr2) &  
                  ,intent(in) :: xx_in
  integer,optional,intent(in) :: ndec_in
  character(ndecim(prcpar)+nxtr+novh  ) :: rslt
  character(ndecim(prcpar)+nxtr+novh+1) :: cc
  character(10) :: aa,bb
  integer :: ndec
  real(kindr2) :: xx     
  xx = xx_in
  if (present(ndec_in)) then ;ndec=ndec_in
                        else ;ndec=ndecim(prcpar)+nxtr
  endif
  write(aa,'(i10)') min(len(cc),ndec+novh+1) ;aa=adjustl(aa)
  write(bb,'(i10)') min(len(cc),ndec       ) ;bb=adjustl(bb)
  aa = '(e'//trim(aa)//'.'//trim(bb)//')'
  write(cc,aa) xx  ;cc=adjustl(cc)
  if (cc(1:2).eq.'-0') then ;rslt = '-'//cc(3:len(cc))
  else                      ;rslt = ' '//cc(2:len(cc))
  endif
  end function

  function printi( ii ) result(rslt)
  integer ,intent(in) :: ii
  character(ndecim(prcpar)) :: rslt
  character(ndecim(prcpar)) :: cc
  character(10) :: aa
  write(aa,'(i10)') ndecim(prcpar) ;aa=adjustl(aa)
  aa = '(i'//trim(aa)//')'
  write(cc,aa) ii ;cc=adjustl(cc)
  if (cc(1:1).ne.'-') then ;rslt=' '//cc
  else                     ;rslt=cc 
  endif
  end function

end module


module avh_olo_qp_auxfun
  use avh_olo_units
  use avh_olo_qp_prec

  implicit none
  private
  public :: mysqrt,eta5,eta3,eta2,sgnIm,sgnRe,kallen
  public :: solabc,rfun,rfun0,solabc_rcc

  interface mysqrt
    module procedure mysqrt_c,mysqrt_cr,mysqrt_ci
  end interface

  interface eta5
    module procedure eta5_0
  end interface
  interface eta3
    module procedure eta3_r,eta3_0
  end interface
  interface eta2
    module procedure eta2_r,eta2_0
  end interface

  interface sgnIm
    module procedure sgnIm_c,sgnIm_ci
  end interface
  interface sgnRe
    module procedure sgnRe_c,sgnRe_r,sgnRe_ri
  end interface

contains


  function mysqrt_c(xx) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! negative imaginary.
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt ,zz
  real(kindr2) &  
    :: xim,xre
  xim = aimag(xx)
  if (xim.eq.RZRO) then
    xre = areal(xx)
    if (xre.ge.RZRO) then
      zz = acmplx(sqrt(xre),0)
    else
      zz = acmplx(0,-sqrt(-xre))
    endif
  else
    zz = sqrt(xx)
  endif
  rslt = zz
  end function

  function mysqrt_cr(xx,sgn) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! imaginary and has the same sign as  sgn .
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  real(kindr2) &  
    ,intent(in) :: sgn
  complex(kindr2) &   
    :: rslt ,zz
  real(kindr2) &  
    :: xim,xre
  xim = aimag(xx)
  if (xim.eq.RZRO) then
    xre = areal(xx)
    if (xre.ge.RZRO) then
      zz = acmplx(sqrt(xre),0)
    else
      zz = acmplx(0,sign(sqrt(-xre),sgn))
    endif
  else
    zz = sqrt(xx)
  endif
  rslt = zz
  end function

  function mysqrt_ci(xx,sgn) result(rslt)
!*******************************************************************
! Returns the square-root of xx .
! If  Im(xx)  is equal zero and  Re(xx)  is negative, the result is
! imaginary and has the same sign as  sgn .
!*******************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: sgn
  complex(kindr2) &   
    :: rslt ,zz
  real(kindr2) &  
    :: xim,xre,hh
  xim = aimag(xx)
  if (xim.eq.RZRO) then
    xre = areal(xx)
    if (xre.ge.RZRO) then
      zz = acmplx(sqrt(xre),0)
    else
      hh = sgn
      zz = acmplx(0,sign(sqrt(-xre),hh))
    endif
  else
    zz = sqrt(xx)
  endif
  rslt = zz
  end function


  subroutine solabc( x1,x2 ,dd ,aa,bb,cc ,imode )
!*******************************************************************
! Returns the solutions  x1,x2  to the equation  aa*x^2+bb*x+cc=0
! Also returns  dd = aa*(x1-x2)
! If  imode=/=0  it uses  dd  as input as value of  sqrt(b^2-4*a*c)
!*******************************************************************
  complex(kindr2) &   
    ,intent(out)   :: x1,x2
  complex(kindr2) &   
    ,intent(inout) :: dd
  complex(kindr2) &   
    ,intent(in) :: aa,bb,cc
  integer         ,intent(in) :: imode
  complex(kindr2) &   
    :: qq,hh
  real(kindr2) &  
    :: r1,r2

  if (aa.eq.CZRO) then
    if (bb.eq.CZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop solabc: ' &
        ,'no solutions, returning 0'
      x1 = 0
      x2 = 0
      dd = 0
    else
      x1 = -cc/bb
      x2 = x1
      dd = bb
    endif
  elseif (cc.eq.CZRO) then
    dd = -bb
    x1 = dd/aa
    x2 = 0
  else
    if (imode.eq.0) dd = sqrt(bb*bb - 4*aa*cc)
    qq = -bb+dd
    hh = -bb-dd
    r1 = abs(qq)
    r2 = abs(hh)
    if (r1.ge.r2) then
      x1 = qq/(2*aa)
      x2 = (2*cc)/qq
    else
      qq = hh
      x2 = qq/(2*aa)
      x1 = (2*cc)/qq
    endif
  endif
  end subroutine


  subroutine solabc_rcc( x1,x2 ,aa,bb,cc )
!*******************************************************************
! Tested
!*******************************************************************
  intent(out) :: x1,x2
  intent(in ) :: aa,bb,cc
  complex(kindr2) &   
    :: x1,x2,bb,cc ,t1,t2
  real(kindr2) &  
    :: aa,xx,yy,pp,qq,uu,vv,pq1,pq2,uv1,uv2,dd,xd1,xd2,yd1,yd2 &
      ,gg,hh,rx1,rx2,ix1,ix2
  if (aa.eq.RZRO) then
    if (bb.eq.CZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop solabc: ' &
        ,'no solutions, returning 0'
      x1 = 0
      x2 = 0
    else
      x1 = -cc/bb
      x2 = x1
    endif
  elseif (cc.eq.CZRO) then
    x1 = -bb/aa
    x2 = 0
  else
    t1 = cc/aa          ;xx= areal(t1) ;yy= aimag(t1)
    t2 = bb/(aa*2)      ;pp=-areal(t2) ;uu=-aimag(t2)
    t2 = sqrt(t2*t2-t1) ;qq= areal(t2) ;vv= aimag(t2)
    pq1=pp+qq ;uv1=uu+vv
    pq2=pp-qq ;uv2=uu-vv
    dd=pq1*pq1+uv1*uv1 ;xd1=xx/dd ;yd1=yy/dd
    dd=pq2*pq2+uv2*uv2 ;xd2=xx/dd ;yd2=yy/dd
    if (abs(pq1).gt.abs(pq2)) then
      rx1 = pq1
      gg=xd1*pq1 ;hh=yd1*uv1
      rx2 = gg+hh
      if (abs(rx2).lt.neglig(prcpar)*max(abs(gg),abs(hh))) rx2 = 0
    elseif (abs(pq2).gt.abs(pq1)) then
      rx2 = pq2
      gg=xd2*pq2 ;hh=yd2*uv2
      rx1 = gg+hh
      if (abs(rx1).lt.neglig(prcpar)*max(abs(gg),abs(hh))) rx1 = 0
    else
      rx1 = pq1
      rx2 = pq2
    endif
    if (abs(uv1).gt.abs(uv2)) then
      ix1 = uv1
      gg=yd1*pq1 ;hh=xd1*uv1
      ix2 = gg-hh
      if (abs(ix2).lt.neglig(prcpar)*max(abs(gg),abs(hh))) ix2 = 0
    elseif (abs(uv2).gt.abs(uv1)) then
      ix2 = uv2
      gg=yd2*pq2 ;hh=xd2*uv2
      ix1 = gg-hh
      if (abs(ix1).lt.neglig(prcpar)*max(abs(gg),abs(hh))) ix1 = 0
    else
      ix1 = uv1
      ix2 = uv2
    endif
    x1 = acmplx(rx1,ix1)
    x2 = acmplx(rx2,ix2)
  endif
  end subroutine


  subroutine rfun(rr,dd ,qq)
!*******************************************************************
! Returns  rr  such that  qq = rr + 1/rr  and  Im(rr)  has the same
! sign as  Im(qq) .
! If  Im(qq)  is zero, then  Im(rr)  is negative or zero.
! If  Im(rr)  is zero, then  |rr| > 1/|rr| .
! Also returns  dd = rr - 1/rr .
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rr,dd
  complex(kindr2) &   
    ,intent(in)  :: qq
  complex(kindr2) &   
    :: r2
  real(kindr2) &  
    :: aa,bb
  integer :: ir,ik
  dd = sqrt(qq*qq-4)
  rr = qq+dd
  r2 = qq-dd
  aa = abs(rr)
  bb = abs(r2)
  if (bb.gt.aa) then
    rr = r2
    dd = -dd
  endif
  aa = aimag(qq)
  bb = aimag(rr)
  if (aa.eq.RZRO) then
    if (bb.le.RZRO) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  else
    ik = sgnRe(aa)
    ir = sgnRe(bb)
    if (ir.eq.ik) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  endif
  end subroutine

  subroutine rfun0(rr ,dd,qq)
!*******************************************************************
! Like rfun, but now  dd  is input, which may get a minus sign
!*******************************************************************
  complex(kindr2) &   
    ,intent(out)   :: rr
  complex(kindr2) &   
    ,intent(inout) :: dd
  complex(kindr2) &   
    ,intent(in)  :: qq
  complex(kindr2) &   
    :: r2
  real(kindr2) &  
    :: aa,bb
  integer :: ir,ik
  rr = qq+dd
  r2 = qq-dd
  aa = abs(rr)
  bb = abs(r2)
  if (bb.gt.aa) then
    rr = r2
    dd = -dd
  endif
  aa = aimag(qq)
  bb = aimag(rr)
  if (aa.eq.RZRO) then
    if (bb.le.RZRO) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  else
    ik = sgnRe(aa)
    ir = sgnRe(bb)
    if (ir.eq.ik) then
      rr = rr/2
    else
      rr = 2/rr
      dd = -dd
    endif
  endif
  end subroutine


  function eta3_r( aa,sa ,bb,sb ,cc,sc ) result(rslt)
!*******************************************************************
! 2*pi*imag times the result of
!     theta(-Im(a))*theta(-Im(b))*theta( Im(c))
!   - theta( Im(a))*theta( Im(b))*theta(-Im(c))
! where a,b,c are interpreted as a+i|eps|sa, b+i|eps|sb, c+i|eps|sc
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb,cc
  real(kindr2) &  
    ,intent(in) :: sa,sb,sc
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: ima,imb,imc
  ima = aimag(aa)
  imb = aimag(bb)
  imc = aimag(cc)
  if (ima.eq.RZRO) ima = sa
  if (imb.eq.RZRO) imb = sb
  if (imc.eq.RZRO) imc = sc
  ima = sgnRe(ima)
  imb = sgnRe(imb)
  imc = sgnRe(imc)
  if (ima.eq.imb.and.ima.ne.imc) then
    rslt = acmplx(0,imc*TWOPI)
  else
    rslt = 0
  endif
  end function

  function eta3_0( aa ,bb ,cc ) result(rslt)
!*******************************************************************
! 2*pi*imag times the result of
!     theta(-Im(a))*theta(-Im(b))*theta( Im(c))
!   - theta( Im(a))*theta( Im(b))*theta(-Im(c))
! where a,b,c are interpreted as a+i|eps|sa, b+i|eps|sb, c+i|eps|sc
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb,cc
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: ima,imb,imc
  ima = sgnIm(aa)
  imb = sgnIm(bb)
  imc = sgnIm(cc)
  if (ima.eq.imb.and.ima.ne.imc) then
    rslt = acmplx(0,imc*TWOPI)
  else
    rslt = 0
  endif
  end function

  function eta5_0( aa ,b1,c1 ,b2,c2 ) result(rslt)
!*******************************************************************
! eta3(aa,b1,c1) - eta3(aa,b2,c2)
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,b1,c1 ,b2,c2
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: imaa,imb1,imc1,imb2,imc2
  imaa = sgnIm(aa)
  imb1 = sgnIm(b1)
  imb2 = sgnIm(b2)
  imc1 = sgnIm(c1)
  imc2 = sgnIm(c2)
  if (imaa.eq.imb1) then
    if (imaa.eq.imb2) then
      if (imc1.eq.imc2) then
        rslt = 0
      elseif (imaa.ne.imc1) then
        rslt = acmplx(0, imc1*TWOPI)
      else
        rslt = acmplx(0,-imc2*TWOPI)
      endif
    elseif (imaa.ne.imc1) then
      rslt = acmplx(0, imc1*TWOPI)
    else
      rslt = 0
    endif
  elseif (imaa.eq.imb2.and.imaa.ne.imc2) then
    rslt = acmplx(0,-imc2*TWOPI)
  else
    rslt = 0
  endif
  end function

  function eta2_r( aa,sa ,bb,sb ) result(rslt)
!*******************************************************************
! The same as  eta3, but with  c=a*b, so that
!   eta(a,b) = log(a*b) - log(a) - log(b)
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb
  real(kindr2) &  
    ,intent(in) :: sa,sb
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rea,reb,ima,imb,imab
  rea = areal(aa)  ;ima = aimag(aa)
  reb = areal(bb)  ;imb = aimag(bb)
  imab = rea*imb + reb*ima
  if (ima .eq.RZRO) ima = sa
  if (imb .eq.RZRO) imb = sb
  if (imab.eq.RZRO) imab = sign(rea,sb) + sign(reb,sa)
  ima  = sgnRe(ima)
  imb  = sgnRe(imb)
  imab = sgnRe(imab)
  if (ima.eq.imb.and.ima.ne.imab) then
    rslt = acmplx(0,imab*TWOPI)
  else
    rslt = 0
  endif
  end function
 
  function eta2_0( aa ,bb ) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: aa,bb
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rea,reb,ima,imb,imab
  rea = areal(aa)  ;ima = aimag(aa)
  reb = areal(bb)  ;imb = aimag(bb)
  rea = rea*imb
  reb = reb*ima
  imab = rea+reb
  ima  = sgnRe(ima)
  imb  = sgnRe(imb)
  imab = sgnRe(imab)
  if (ima.eq.imb.and.ima.ne.imab) then
    rslt = acmplx(0,imab*TWOPI)
  else
    rslt = 0
  endif
  end function 


  function kallen( p1,p2,p3 ) result(rslt)
!*******************************************************************
!  p1^2 + p2^2 + p3^2 - 2*p1*p2 - 2*p2*p3 - 2*p3*p1
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3
  complex(kindr2) &   
    :: rslt ,y1,y2,y3
  real(kindr2) &  
    :: b1,b2,b3
  y1=p2*p3 ;b1=areal(y1)
  y2=p3*p1 ;b2=areal(y2)
  y3=p1*p2 ;b3=areal(y3)
      if (b1.le.RZRO) then  ;rslt = (p1-p2-p3)**2 - 4*y1
  elseif (b2.le.RZRO) then  ;rslt = (p2-p3-p1)**2 - 4*y2
  elseif (b3.le.RZRO) then  ;rslt = (p3-p1-p2)**2 - 4*y3
  elseif (b1.le.b2.and.b1.le.b3) then  ;rslt = (p1-p2-p3)**2 - 4*y1
  elseif (b2.le.b3.and.b2.le.b1) then  ;rslt = (p2-p3-p1)**2 - 4*y2
                                 else  ;rslt = (p3-p1-p2)**2 - 4*y3
  endif
  end function


  function sgnIm_c(zz) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: zz
  integer :: rslt
  real(kindr2) &  
    :: imz
  imz = aimag(zz)
  if (imz.ge.RZRO) then ;rslt= 1
                   else ;rslt=-1
  endif
  end function

  function sgnIm_ci(zz,ii) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
          ,intent(in) :: zz
  integer ,intent(in) :: ii
  integer :: rslt
  real(kindr2) &  
    :: imz
  imz = aimag(zz)
  if     (imz.gt.RZRO) then ;rslt= 1
  elseif (imz.lt.RZRO) then ;rslt=-1
                       else ;rslt= sign(1,ii)
  endif
  end function

  function sgnRe_c(zz) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: zz
  integer :: rslt
  real(kindr2) &  
    :: rez
  rez = zz
  if (rez.ge.RZRO) then ;rslt= 1
                   else ;rslt=-1
  endif
  end function

  function sgnRe_r(rez) result(rslt)
!*******************************************************************
!*******************************************************************
  real(kindr2) &  
    ,intent(in) :: rez
  integer :: rslt
  if (rez.ge.RZRO) then ;rslt= 1
                   else ;rslt=-1
  endif
  end function

  function sgnRe_ri(rez,ii) result(rslt)
!*******************************************************************
!*******************************************************************
  real(kindr2) &  
          ,intent(in) :: rez
  integer ,intent(in) :: ii
  integer :: rslt
  if     (rez.gt.RZRO) then ;rslt= 1
  elseif (rez.lt.RZRO) then ;rslt=-1
                       else ;rslt=sign(1,ii)
  endif
  end function

end module


module avh_olo_qp_olog
!***********************************************************************
! Provides the functions
!   olog(x,n) = log(x) + n*pi*imag  
!   olog1(x,n) = olog(x,n)/(x-1)
!   olog2(x,n) = ( olog1(x,n) - 1 )/(x-1)
!   olog3(x,n) = ( olog2(x,n) + 1/2 )/(x-1)
! In the vicinity of x=1,n=0, the logarithm of complex argument is
! evaluated with a series expansion.
!***********************************************************************
  use avh_olo_units
  use avh_olo_qp_prec
  use avh_olo_qp_print
  use avh_olo_qp_auxfun
  implicit none
  private
  public :: update_olog,olog,olog1,olog2,olog3

  real(kindr2) &  
         ,allocatable,save :: thrs(:,:)
  integer,allocatable,save :: ntrm(:,:)
  integer,parameter :: nStp=6

  interface olog
    module procedure log_c,log_r
  end interface
  interface olog1
    module procedure log1_c,log1_r
  end interface
  interface olog2
    module procedure log2_c,log2_r
  end interface
  interface olog3
    module procedure log3_c,log3_r
  end interface

contains

  subroutine update_olog
!***********************************************************************
!***********************************************************************
  use avh_olo_qp_arrays
  real(kindr2) &  
    :: tt
  integer :: nn,mm,ii,jj
!  real(kind(1d0)) :: xx(6) !DEBUG
  if (allocated(thrs)) then
    call shift2( thrs ,prcpar )
    call shift2( ntrm ,prcpar )
  else
    allocate(thrs(1:nStp,1:1))
    allocate(ntrm(1:nStp,1:1))
    if (prcpar.ne.1) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop update_olog'
      stop
    endif
  endif
  if (prcpar.gt.1) then ;nn=ntrm(nStp,prcpar-1)-1
                   else ;nn=1
  endif
  do
    nn = nn+1
    mm = 2*nn-1
    tt = 1
    tt = (EPSN*mm)**(tt/(mm-1))
    tt = 2*tt/(1-tt)
! expansion from x=1+d with |d|=1/1000
    if (1000*tt.gt.RONE) exit
  enddo
  ntrm(nStp,prcpar) = nn
  thrs(nStp,prcpar) = tt
  nn = max(1,nint(nn*1d0/nStp))
  do ii=nStp-1,1,-1
    ntrm(ii,prcpar) = ntrm(ii+1,prcpar)-nn
    if (ntrm(ii,prcpar).le.1) then
      do jj=1,ii
        ntrm(jj,prcpar) = ntrm(ii,prcpar)
        thrs(jj,prcpar) = 0 
      enddo
      exit
    endif
    mm = 2*ntrm(ii,prcpar)-1
    tt = 1
    tt = (EPSN*mm)**(tt/(mm-1))
    thrs(ii,prcpar) = 2*tt/(1-tt)
  enddo
!  do ii=lbound(thrs,2),ubound(thrs,2) !DEBUG
!    do jj=1,nStp                      !DEBUG
!      xx(jj) = thrs(jj,ii)            !DEBUG
!    enddo                             !DEBUG
!    write(*,'(99e10.3)') xx(:)        !DEBUG
!    write(*,'(99i10)'  ) ntrm(:,ii)   !DEBUG
!  enddo                               !DEBUG
  end subroutine


  function log_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,iyy
!
  rex = areal(xx)
  imx = aimag(xx)
  iyy = iph
!
  if (abs(imx).le.EPSN*abs(rex)) then
    if (rex.ge.RZRO) then
      rslt = log_r( rex, iyy )
    else
      rslt = log_r(-rex, iyy+sgnRe(imx) )
    endif
    return
  endif
!
  if (mod(iyy,2).eq.0) then
    yy = acmplx(rex,imx)
  else
    yy = acmplx(-rex,-imx)
    iyy = iyy+sgnRe(imx)
  endif
!
  if (iyy.ne.0) then
    rslt = log(yy) + IPI*iyy
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = log(yy)
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = 2 * ( z + z^3/3 + z^5/5 + z^7/7 + ... )  
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,1,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = zz*rslt
  end function


  function log_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rr
  integer :: jj
!
  if (xx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log_r: ' &
       ,'xx =',trim(myprint(xx)),', returning 0'
    rslt = 0
    return
  elseif (xx.gt.RZRO) then ;rr= xx ;jj= iph
                      else ;rr=-xx ;jj= iph+1 ! log(-1)=i*pi
  endif
!
  rslt = log(rr) + IPI*jj
  end function


  function log1_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,jj
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (abs(imx).le.EPSN*abs(rex)) then
    if (rex.ge.RZRO) then
      rslt = log1_r( rex, iph )
    else
      rslt = log1_r(-rex, iph+sgnRe(imx) )
    endif
    return
  endif
!
  if (mod(iph,2).eq.0) then ;yy= xx ;jj=iph
                       else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  if (jj.ne.0) then
    rslt = ( log(yy) + IPI*jj )/(yy-1)
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = log(yy)/zz
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = 2/(y+1) * ( 1 + z^2/3 + z^4/5 + z^6/7 + ... )  
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,1,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = rslt/(yy+1)
  end function


  function log1_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: rr,yy
  integer :: jj
!  include 'avh_olo_qp_real.h90'
!    :: aa,zz,z2
!  integer :: nn,ii
!
  if (xx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log1_r: ' &
       ,'xx =',trim(myprint(xx)),', returning 0'
    rslt = 0
    return
  elseif (xx.gt.RZRO) then ;rr= xx ;jj=iph
                      else ;rr=-xx ;jj=iph+1 ! log(-1)=i*pi
  endif
!
  yy=rr ;if (mod(jj,2).ne.0) yy=-rr
!
  if (abs(yy-1).le.10*EPSN) then
    if (jj.ne.0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log1_r: ' &
        ,'rr,jj =',trim(myprint(rr)),jj,', putting jj to 0'
    endif
    rslt = 1 - (yy-1)/2
    return
  endif
!
  rslt = ( log(rr) + IPI*jj )/(yy-1)
  end function


  function log2_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  rslt = log2_c(xx*CONE,iph)
  end function


  function log2_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,jj
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (rex.eq.RZRO.and.imx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log2_c: ' &
       ,'xx = 0, returning 0'
    rslt = 0
    return
  endif
!
  if (mod(iph,2).eq.0) then ;yy= xx ;jj=iph
                       else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  if (jj.ne.0) then
    rslt = ( olog1(yy,jj) - 1 )/(yy-1)
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = (log(yy)/zz-1)/zz
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = -1/(y+1) + 2/(y+1)^2 * ( z/3 + z^3/5 + z^5/7 + ... )  
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,3,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = ( -1 + zz*rslt/(yy+1) )/(yy+1)
  end function


  function log3_r(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  rslt = log3_c(xx*CONE,iph)
  end function


  function log3_c(xx,iph) result(rslt)
!***********************************************************************
!***********************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,zz,z2,HLF
  real(kindr2) &  
    :: aa,rex,imx
  integer :: nn,ii,jj
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (rex.eq.RZRO.and.imx.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop log2_c: ' &
       ,'xx = 0, returning 0'
    rslt = 0
    return
  endif
!
  HLF = CONE/2
!
  if (mod(iph,2).eq.0) then ;yy= xx ;jj=iph
                       else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  if (jj.ne.0) then
    rslt = ( olog2(xx,jj) + HLF )/(yy-1)
    return
  endif
!
  zz = yy-1
  aa = abs(zz)
  if     (aa.ge.thrs(6,prcpar)) then
    rslt = ((log(yy)/zz-1)/zz+HLF)/zz
    return
  elseif (aa.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (aa.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (aa.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (aa.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (aa.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! convergence acceleration using  z=(y-1)/(y+1)
! rslt = 1/(2*(y+1)) + 2/(y+1)^3 * ( 1/3 + z^2/5 + z^4/7 + ... )
  zz = zz/(yy+1)
  z2 = zz*zz
  aa = 2
  nn = 2*nn-1
  rslt = aa/nn
  do ii=nn-2,3,-2
    rslt = aa/ii + z2*rslt
  enddo
  rslt = ( HLF + rslt/(yy+1)**2 )/(yy+1)
  end function

end module




module avh_olo_qp_dilog
!***********************************************************************
!                     /1    ln(1-zz*t)
!   dilog(xx,iph) = - |  dt ---------- 
!                     /0        t
! with  zz = 1 - xx*exp(imag*pi*iph)  [pi, NOT 2*pi]
!
!   dilog(x1,i1,x2,i2) = ( dilog(x1,i1)-dilog(x2,i2) )/( x1-x2 )
!
! Arguments xx,x1,x2, may be all real or all complex,
! arguments iph,i1,i2 must be all integer.
!***********************************************************************
  use avh_olo_units
  use avh_olo_qp_prec
  use avh_olo_qp_print
  use avh_olo_qp_auxfun
  use avh_olo_qp_arrays
  implicit none
  private
  public :: update_dilog,dilog

  real(kindr2) &  
         ,allocatable,save :: coeff(:)
  real(kindr2) &  
         ,allocatable,save :: thrs(:,:)
  integer,allocatable,save :: ntrm(:,:)
  integer,parameter :: nStp=6

  real(kindr2) &  
         ,allocatable :: bern(:),fact(:)

  interface dilog
    module procedure dilog_c,dilog_r,dilog2_c,dilog2_r
  end interface

contains

  subroutine update_dilog
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
    :: tt
  integer :: nn,ii,jj
  logical :: highestSoFar
!  real(kind(1d0)) :: xx(6) !DEBUG
!
  if (allocated(thrs)) then
    call shift2( thrs ,prcpar )
    call shift2( ntrm ,prcpar )
  else
    allocate(thrs(1:nStp,1:1))
    allocate(ntrm(1:nStp,1:1))
    if (prcpar.ne.1) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop update_dilog'
      stop
    endif
  endif
!
  highestSoFar = prcpar.eq.ubound(ntrm,2)
  if (highestSoFar) then
    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(0:-1)) ! allocate at size=0
  endif
!
  if (prcpar.gt.1) then ;nn=ntrm(nStp,prcpar-1)-1
                   else ;nn=2
  endif
!
  do
    nn = nn+1
    if (nn.gt.ubound(coeff,1)) call update_coeff( 2*nn )
    tt = 1
    tt = (EPSN/abs(coeff(nn)))**(tt/(2*nn))
! expansion parameter is smaller than 1.05
    if (100*tt.gt.105*RONE) exit
  enddo
!
  if (highestSoFar) call resize( coeff ,0,nn )
!
  ntrm(nStp,prcpar) = nn
  thrs(nStp,prcpar) = tt
  nn = max(1,nint(nn*1d0/nStp))
  do ii=nStp-1,1,-1
    ntrm(ii,prcpar) = ntrm(ii+1,prcpar)-nn
    if (ntrm(ii,prcpar).le.2) then
      do jj=1,ii
        ntrm(jj,prcpar) = max(2,ntrm(ii,prcpar))
        thrs(jj,prcpar) = 0 
      enddo
      exit
    endif
    jj = ntrm(ii,prcpar)
    tt = 1
    tt = (EPSN/abs(coeff(jj)))**(tt/(2*jj))
    thrs(ii,prcpar) = tt
  enddo
!
  if (allocated(bern)) deallocate(bern)
  if (allocated(fact)) deallocate(fact)
!
!  do ii=lbound(thrs,2),ubound(thrs,2) !DEBUG
!    do jj=1,nStp                      !DEBUG
!      xx(jj) = thrs(jj,ii)            !DEBUG
!    enddo                             !DEBUG
!    write(*,'(99e10.3)') xx(:)        !DEBUG
!    write(*,'(99i10)'  ) ntrm(:,ii)   !DEBUG
!  enddo                               !DEBUG
  end subroutine


  subroutine update_coeff( ncf )
!*******************************************************************
!   coeff(0)=-1/4
!   coeff(n)=bern(2*n)/(2*n+1)
!    bern(n)=bernoulli(n)/n!
!    fact(n)=n!
! DO NOT SKIP THE ODD bern IN THE RECURSIVE LOOP
! DO NOT PUT THE ODD bern TO ZERO
!*******************************************************************
  integer ,intent(in) :: ncf
  integer :: ii,jj,nbern,nold
!
  if (allocated(bern)) then ;nold=ubound(bern,1)
                       else ;nold=0
  endif
!
  nbern = 2*ncf
!
  call enlarge( bern  ,1,nbern   )
  call enlarge( fact  ,0,nbern+1 )
  call enlarge( coeff ,0,ncf     )
!
  fact(0) = 1
  do ii=nold+1,nbern+1
    fact(ii) = fact(ii-1)*ii
  enddo
!
  do ii=nold+1,nbern
    bern(ii) = -1/fact(ii+1)
    do jj=1,ii-1
      bern(ii) = bern(ii) - bern(jj)/fact(ii+1-jj)
    enddo
  enddo
!
  coeff(0) = 1
  coeff(0) =-coeff(0)/4
  do ii=nold+2,nbern,2
    coeff(ii/2) = bern(ii)/(ii+1)
  enddo
!
  end subroutine


  function dilog_c(xx,iph) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt ,yy,lyy,loy,zz,z2
  real(kindr2) &  
    :: rex,imx,az
  integer :: ii,jj,ntwo,odd,nn
  logical :: r_gt_1 , y_lt_h
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (abs(imx).le.EPSN*abs(rex)) then
    if (rex.ge.RZRO) then
      rslt = dilog_r( rex, iph )
    else
      rslt = dilog_r(-rex, iph+sgnRe(imx) )
    endif
    return
  endif
!
  if (rex.gt.RZRO) then ;yy= xx ;jj=iph
                   else ;yy=-xx ;jj=iph+sgnRe(imx)
  endif
!
  odd = mod(jj,2)
  ntwo = jj-odd
! 
  r_gt_1 = (rex*rex+imx*imx.gt.RONE)
  lyy = log(yy)
  if (odd.ne.0) yy = -yy
!
  if (r_gt_1) then
    yy   = 1/yy
    lyy  =-lyy
    ntwo =-ntwo
    odd  =-odd
  endif
  loy = log(1-yy)
!
  y_lt_h = (2*areal(yy).lt.RONE)
  if (y_lt_h) then ;zz=-loy
              else ;zz=-lyy
  endif
!
  az = abs(zz)
! if (az.gt.thrs(6,prcpar)) ERROR az to big 
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
  z2 = zz*zz
  rslt = coeff(nn)
  do ii=nn,2,-1
    rslt = coeff(ii-1) + z2*rslt
  enddo
  rslt = zz*( 1 + zz*( coeff(0) + zz*rslt ) )
!
  if (y_lt_h) then
    rslt = 4*PISQo24 - rslt - loy*(lyy+IPI*(ntwo+odd))
  else
    rslt = rslt - loy*IPI*ntwo
  endif
!
  if (r_gt_1) rslt = -rslt - (lyy+IPI*(ntwo+odd))**2/2
  end function



  function dilog_r(xx,iph) result(rslt)
!*******************************************************************
!*******************************************************************
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: iph
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: yy,lyy,loy,zz,z2,liox,az
  integer :: jj,ii,ntwo,odd,nn
  logical :: r_gt_1 , y_lt_h
!
  if (xx.eq.RZRO) then
    rslt = 4*PISQo24
    return
  elseif (xx.gt.RZRO) then ;yy= xx ;jj=iph
                      else ;yy=-xx ;jj=iph+1 ! log(-1)=i*pi
  endif
!
  odd = mod(jj,2)
  ntwo = jj-odd
! 
  if (yy.eq.RONE.and.odd.eq.0) then
    if (ntwo.ne.0) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog_r: ' &
        ,'|x|,iph = ',trim(myprint(yy)),',',jj,', returning 0'
    endif
    rslt = 0
    return
  endif
!
  r_gt_1 = (yy.gt.RONE)
  lyy = log(yy)
  if (odd.ne.0) yy = -yy
!
  if (r_gt_1) then
    yy   = 1/yy
    lyy  =-lyy
    ntwo =-ntwo
    odd  =-odd
  endif
  loy = log(1-yy) ! log(1-yy) is always real
!
  y_lt_h = (2*yy.lt.RONE)
  if (y_lt_h) then
    zz = -loy ! log(1-yy) is real
  else
    zz = -lyy ! yy>0.5 => log(yy) is real
  endif
!
  az = abs(zz)
! if (az.gt.thrs(6,prcpar)) ERROR az to big 
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
  z2 = zz*zz
  liox = coeff(nn)
  do ii=nn,2,-1
    liox = coeff(ii-1) + z2*liox
  enddo
  liox = zz*( 1 + zz*( coeff(0) + zz*liox ) )
!
  rslt = acmplx(liox)
!
  if (y_lt_h) then
    rslt = 4*PISQo24 - rslt - acmplx(loy*lyy,loy*ONEPI*(ntwo+odd))
  else
    rslt = rslt + acmplx( 0 ,-loy*ONEPI*ntwo )
  endif
!
  if (r_gt_1) rslt = -rslt - acmplx(lyy,ONEPI*(ntwo+odd))**2/2
  end function


  function dilog2_c( x1,i1 ,x2,i2 ) result(rslt)
!*******************************************************************
!*******************************************************************
  use avh_olo_qp_olog
  complex(kindr2) &   
          ,intent(in) :: x1,x2
  integer ,intent(in) :: i1,i2
  complex(kindr2) &   
    :: rslt ,y1,y2 ,ff,gg,logr1,logr2,logo1,logo2,r1,r2,rr
  real(kindr2) &  
    :: eps ,re1,im1,re2,im2,a1,a2,aa,ao1,ao2
  integer :: j1,j2,ii,nn,oo
  integer,parameter :: pp(-1:1,-1:1)=&
                      reshape((/-2,-2,2 ,-2,0,2 ,-2,2,2/),(/3,3/))
!
  re1=areal(x1) ;re2=areal(x2)
  im1=aimag(x1) ;im2=aimag(x2)
!
  if (abs(im1).le.EPSN*abs(re1).and.abs(im2).le.EPSN*abs(re2)) then
    if (re1.ge.RZRO) then
      if (re2.ge.RZRO) then
        rslt = dilog2_r( re1,i1 , re2,i2 )
      else
        rslt = dilog2_r( re1,i1 ,-re2,i2+sgnRe(im2) )
      endif
    elseif (re2.ge.RZRO) then
      rslt = dilog2_r(-re1,i1+sgnRe(im1) , re2,i2 )
    else
      rslt = dilog2_r(-re1,i1+sgnRe(im1) ,-re2,i2+sgnRe(im2) )
    endif
    return
  endif
!
  if (re1.ge.RZRO) then ;r1= x1 ;j1=i1
                   else ;r1=-x1 ;j1=i1+sgnRe(im1,1)
  endif
  if (re2.ge.RZRO) then ;r2= x2 ;j2=i2
                   else ;r2=-x2 ;j2=i2+sgnRe(im2,1)
  endif
!
  a1=abs(r1) ;a2=abs(r2)
  if (a1.gt.a2) then
    aa=a1;a1=a2;a2=aa
    rr=r1;r1=r2;r2=rr
    ii=j1;j1=j2;j2=ii
  endif
!
  oo=mod(j1,2) ;nn=j1-oo ;y1=r1 ;if (oo.ne.0) y1=-y1
  oo=mod(j2,2) ;nn=j2-oo ;y2=r2 ;if (oo.ne.0) y2=-y2
!
  eps = 10*EPSN
!
  if (j1.ne.j2) then
    if (r1.eq.r2) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'j1,j2,r1-r2',j1,j2,',',trim(myprint(r1-r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_c j1=/=j2,r1=r2' !DEBUG
      return
    else
      rslt = ( dilog_c(r1,j1)-dilog_c(r2,j2) )/(y1-y2)
!      write(*,*) 'dilog2_c j1=/=j2' !DEBUG
      return
    endif
  endif
!
  if (a1.lt.eps) then
    if (a2.lt.eps) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r1,r2 =',trim(myprint(r1)),',',trim(myprint(r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_c r1<eps,r2<eps' !DEBUG
      return
    else
      rslt = (dilog_c(r2,j2)-4*PISQo24)/y2
!      write(*,*) 'dilog2_c r1<eps' !DEBUG
      return
    endif
  endif
!
  logr1=log(r1) ;logr2=log(r2)
!
  ao1=abs(1-y1) ;ao2=abs(1-y2)
  if (10*ao1.lt.RONE.or.10*ao2.lt.RONE) then
    aa = abs(r1/r2-1)
    if (10*aa.gt.RONE) then
      rslt = (dilog_c(r1,j1)-dilog_c(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_c ||1-y1|/|1-y2|-1|>0.1' !DEBUG
      return
    elseif (oo.eq.0.and.ao1.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r1,oo,nn =',trim(myprint(r1)),',',oo,nn,', putting nn=0'
      if (ao2.lt.eps) then
        rslt = -1
!        write(*,*) 'dilog2_c |1-y1|' !DEBUG
        return
      else
        y1=1-eps ;nn=0 ;logr1=0 ;r1=1-eps
      endif
    elseif (oo.eq.0.and.ao2.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r2,oo,nn =',trim(myprint(r2)),',',oo,nn,', putting nn=0'
      y2=1-eps ;nn=0 ;logr2=0 ;r2=1-eps
    endif
  else
    aa = abs((logr1+oo*IPI)/(logr2+oo*IPI)-1)
    if (10*aa.gt.RONE) then
      rslt = (dilog_c(r1,j1)-dilog_c(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_c |logr1/logr2-1|>0.1',logr1,logr2 !DEBUG
      return
    elseif (aa.lt.eps) then
      ii = 0
      if (a1.gt.RONE) ii = ii + (nn+pp(oo,sgnIm(y2)))
      if (a2.gt.RONE) ii = ii - (nn+pp(oo,sgnIm(y2)))
      ii = nn*ii
      if (ii.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_c: ' &
        ,'r1,r2,nn =',trim(myprint(r1)),',',trim(myprint(r2)),',',nn &
        ,', putting nn=0'
      rslt = -olog1(y2,0)
!      write(*,*) 'dilog2_c |logr1/lorg2|<eps' !DEBUG
      return
    endif
  endif
!
  if (a1.gt.RONE) then
    y1=1/y1 ;logr1=-logr1
    y2=1/y2 ;logr2=-logr2
    nn=-nn ;oo=-oo
  endif
!
  ff=y1/y2         ;ff=-olog1(ff,0)/y2
  gg=(1-y1)/(1-y2) ;gg=-olog1(gg,0)/(1-y2)
!
  if (2*areal(y1).ge.RONE) then
!    write(*,*) 'dilog2_c re>1/2' !DEBUG
    rslt = ff*sumterms_c(-logr1,-logr2) - nn*IPI*gg
  else
!    write(*,*) 'dilog2_c re<1/2' !DEBUG
    logo1 = log(1-y1)
    logo2 = log(1-y2)
    rslt = gg*( sumterms_c(-logo1,-logo2) - (nn+oo)*IPI - logr2 ) + ff*logo1
  endif
!
  if (a1.gt.RONE) then !implies also r2>1
!    write(*,*) 'dilog2_c r1>1,r2>1' !DEBUG
    rslt = y1*y2*( rslt - ff*((logr1+logr2)/2 + (nn+oo)*IPI) )
  elseif (a2.gt.RONE.and.nn.ne.0) then
!    write(*,*) 'dilog2_c r1<1,r2>1',oo,sgnIm(y2)!DEBUG
    rslt = rslt - 12*nn*( nn + pp(oo,sgnIm(y2)) )*PISQo24/(y1-y2)
  endif
!
  end function


  function dilog2_r( x1,i1 ,x2,i2 ) result(rslt)
!*******************************************************************
!*******************************************************************
  use avh_olo_qp_olog
  real(kindr2) &  
          ,intent(in) :: x1,x2
  integer ,intent(in) :: i1,i2
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: y1,y2 ,ff,gg,logr1,logr2,logo1,logo2
  real(kindr2) &  
    :: eps,r1,r2,rr,ro1,ro2
  integer :: j1,j2,ii,nn,oo
!
  if (x1.ge.RZRO) then ;r1= x1 ;j1=i1
                  else ;r1=-x1 ;j1=i1+1 ! log(-1)=i*pi
  endif
  if (x2.ge.RZRO) then ;r2= x2 ;j2=i2
                  else ;r2=-x2 ;j2=i2+1 ! log(-1)=i*pi
  endif
!
  if (r1.gt.r2) then
    rr=r1;r1=r2;r2=rr
    ii=j1;j1=j2;j2=ii
  endif
!
  oo=mod(j1,2) ;nn=j1-oo ;y1=r1 ;if (oo.ne.0) y1=-y1
  oo=mod(j2,2) ;nn=j2-oo ;y2=r2 ;if (oo.ne.0) y2=-y2
!
  eps = 10*EPSN
!
  if (j1.ne.j2) then
    if (r1.eq.r2) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'j1,j2,r1-r2',j1,j2,',',trim(myprint(r1-r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_r j1=/=j2,r1=r2' !DEBUG
      return
    else
      rslt = ( dilog_r(r1,j1)-dilog_r(r2,j2) )/(y1-y2)
!      write(*,*) 'dilog2_r j1=/=j2' !DEBUG
      return
    endif
  endif
!
  if (r1.lt.eps) then
    if (r2.lt.eps) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r1,r2 =',trim(myprint(r1)),',',trim(myprint(r2)),', returning 0'
      rslt = 0
!      write(*,*) 'dilog2_r r1<eps,r2<eps' !DEBUG
      return
    else
      rslt = (dilog_r(r2,j2)-4*PISQo24)/y2
!      write(*,*) 'dilog2_r r1<eps' !DEBUG
      return
    endif
  endif
!
  logr1=log(r1) ;logr2=log(r2)
!
  ro1=abs(1-y1) ;ro2=abs(1-y2)
  if (10*ro1.lt.RONE.or.10*ro2.lt.RONE) then
    rr = abs(r1/r2-1)
    if (10*rr.gt.RONE) then
      rslt = (dilog_r(r1,j1)-dilog_r(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_r ||1-y1|/|1-y2|-1|>0.1' !DEBUG
      return
    elseif (oo.eq.0.and.ro1.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r1,oo,nn =',trim(myprint(r1)),',',oo,nn,', putting nn=0'
      if (ro2.lt.eps) then
        rslt = -1
!        write(*,*) 'dilog2_r |1-y1|' !DEBUG
        return
      else
        y1=1-eps ;nn=0 ;logr1=0 ;r1=1-eps
      endif
    elseif (oo.eq.0.and.ro2.lt.eps) then
      if (nn.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r2,oo,nn =',trim(myprint(r2)),',',oo,nn,', putting nn=0'
      y2=1-eps ;nn=0 ;logr2=0 ;r2=1-eps
    endif
  else
    rr = abs((logr1+oo*IPI)/(logr2+oo*IPI)-1)
    if (10*rr.gt.RONE) then
      rslt = (dilog_r(r1,j1)-dilog_r(r2,j2))/(y1-y2)
!      write(*,*) 'dilog2_r |logr1/logr2-1|>0.1',logr1,logr2 !DEBUG
      return
    elseif (rr.lt.eps) then
      ii = 0
      if (r1.gt.RONE) ii = ii + (nn+2*oo)
      if (r2.gt.RONE) ii = ii - (nn+2*oo)
      ii = nn*ii
      if (ii.ne.0.and.eunit.gt.0) write(eunit,*) 'ERROR in OneLOop dilog2_r: ' &
        ,'r1,r2,nn =',trim(myprint(r1)),',',trim(myprint(r2)),',',nn &
        ,', putting nn=0'
      rslt = -olog1(y2,0)
!      write(*,*) 'dilog2_r |logr1/lorg2|<eps' !DEBUG
      return
    endif
  endif
!
  if (r1.gt.RONE) then
    y1=1/y1 ;logr1=-logr1
    y2=1/y2 ;logr2=-logr2
    nn=-nn ;oo=-oo
  endif
!
  ff=y1/y2         ;ff=-olog1(ff,0)/y2
  gg=(1-y1)/(1-y2) ;gg=-olog1(gg,0)/(1-y2)
!
  if (2*y1.ge.RONE) then
!    write(*,*) 'dilog2_r re>1/2' !DEBUG
    rslt = ff*sumterms_r(-logr1,-logr2) - nn*IPI*gg
  else
!    write(*,*) 'dilog2_r re<1/2' !DEBUG
    logo1 = log(1-y1)
    logo2 = log(1-y2)
    rslt = gg*( sumterms_r(-logo1,-logo2) - (nn+oo)*IPI - logr2 ) + ff*logo1
  endif
!
  if (r1.gt.RONE) then !implies also r2>1
!    write(*,*) 'dilog2_r r1>1,r2>1' !DEBUG
    rslt = y1*y2*( rslt - ff*((logr1+logr2)/2 + (nn+oo)*IPI) )
  elseif (r2.gt.RONE.and.nn.ne.0) then
!    write(*,*) 'dilog2_r r1<1,r2>1' !DEBUG
    rslt = rslt - 12*nn*PISQo24*(nn+2*oo)/(y1-y2)
  endif
!
  end function


  function sumterms_c( z1,z2 ) result(rslt)
!***********************************************************************
! ( f(z1)-f(z2) )/( z1-z2 ), where
! f(z)= z + c0*z^2 + c1*z^3 + c2*z^5 + c3*z^7 + ...
!***********************************************************************
  complex(kindr2) &   
    ,intent(in) :: z1,z2
  complex(kindr2) &   
    :: rslt,yy,zz
  real(kindr2) &  
    :: az
  integer :: ii,nn
  az = max(abs(z1),abs(z2))
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
! calculates all z(i)=(z1^i-z2^i)/(z1-z2) numerically stable
!  zz(1) = 1
!  yy    = 1
!  do ii=2,2*nn+1
!    yy = z2*yy
!    zz(ii) = z1*zz(ii-1) + yy
!  enddo
  zz = 1
  yy = 1
  rslt = zz
  yy = z2*yy
  zz = z1*zz+yy
  rslt = rslt + coeff(0)*zz
  do ii=1,nn
    yy = z2*yy
    zz = z1*zz+yy
    rslt = rslt + coeff(ii)*zz
    yy = z2*yy
    zz = z1*zz+yy
  enddo
  end function  


  function sumterms_r( z1,z2 ) result(rslt)
!***********************************************************************
! ( f(z1)-f(z2) )/( z1-z2 ), where
! f(z)= z + c0*z^2 + c1*z^3 + c2*z^5 + c3*z^7 + ...
!***********************************************************************
  real(kindr2) &  
    ,intent(in) :: z1,z2
  real(kindr2) &  
    :: rslt,yy,zz
  real(kindr2) &  
    :: az
  integer :: ii,nn
  az = max(abs(z1),abs(z2))
  if     (az.ge.thrs(5,prcpar)) then ;nn=ntrm(6,prcpar)
  elseif (az.ge.thrs(4,prcpar)) then ;nn=ntrm(5,prcpar)
  elseif (az.ge.thrs(3,prcpar)) then ;nn=ntrm(4,prcpar)
  elseif (az.ge.thrs(2,prcpar)) then ;nn=ntrm(3,prcpar)
  elseif (az.ge.thrs(1,prcpar)) then ;nn=ntrm(2,prcpar)
                                else ;nn=ntrm(1,prcpar)
  endif
  zz = 1
  yy = 1
  rslt = zz
  yy = z2*yy
  zz = z1*zz+yy
  rslt = rslt + coeff(0)*zz
  do ii=1,nn
    yy = z2*yy
    zz = z1*zz+yy
    rslt = rslt + coeff(ii)*zz
    yy = z2*yy
    zz = z1*zz+yy
  enddo
  end function  

end module


module avh_olo_qp_bnlog
!***********************************************************************
!                      /1    
!   bnlog(n,x) = (n+1) |  dt t^n ln(1-t/x) 
!                      /0 
!***********************************************************************
  use avh_olo_units
  use avh_olo_qp_prec
  use avh_olo_qp_auxfun
  use avh_olo_qp_arrays
  use avh_olo_qp_olog
  use avh_olo_qp_print
  implicit none
  private
  public :: update_bnlog,bnlog

  real(kindr2) &  
         ,allocatable,save :: coeff(:,:)
  real(kindr2) &  
         ,allocatable,save :: thrs(:,:,:)
  integer,allocatable,save :: ntrm(:,:,:)
  integer,parameter :: nStp=6
  integer,parameter :: rank=4
  integer,parameter :: aCoef(0:rank,0:rank)=reshape((/ &
                         1, 0, 0, 0, 0 & ! 1
                       , 1, 2, 0, 0, 0 & ! 1/2,1
                       , 2, 3, 6, 0, 0 & ! 1/3,1/2,1
                       , 3, 4, 6,12, 0 & ! 1/4,1/3,1/2,1
                       ,12,15,20,30,60 & ! 1/5,1/4,1/3,1/2,1
                       /),(/rank+1,rank+1/))

  interface bnlog
    module procedure bnlog_c,bnlog_r
  end interface

contains


  subroutine update_bnlog
!***********************************************************************
!***********************************************************************
  real(kindr2) &  
    :: tt
  integer :: nn,ii,jj,n1,nmax,irank
  logical :: highestSoFar
!  real(kind(1d0)) :: xx(6) !DEBUG
!
  if (allocated(thrs)) then
    call shift3( thrs ,prcpar )
    call shift3( ntrm ,prcpar )
  else
    allocate(thrs(1:nStp,0:rank,1:1))
    allocate(ntrm(1:nStp,0:rank,1:1))
    if (prcpar.ne.1) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop update_bnlog'
      stop
    endif
  endif
!
  highestSoFar = prcpar.eq.ubound(ntrm,3)
!
  if (highestSoFar) then
    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(0:-1,0:2)) ! allocate at size=0
  endif
!
  nmax = 0
!
  do irank=0,rank
!
    n1 = 2+irank
!
    if (prcpar.gt.1) then ;nn=ntrm(nStp,irank,prcpar-1)-1
                     else ;nn=n1
    endif
!  
    do
      nn = nn+1
      if (highestSoFar.and.nn.gt.ubound(coeff,1)) call update_coeff( 2*nn )
      tt = 1
      tt = (EPSN*abs(coeff(n1,irank)/coeff(nn,irank)))**(tt/(nn-n1))
      if (8*(irank+1)*tt.gt.RONE) exit
    enddo
!
    if (nn.gt.nmax) nmax=nn
!  
    ntrm(nStp,irank,prcpar) = nn
    thrs(nStp,irank,prcpar) = tt
    nn = max(1,nint(nn*1d0/nStp))
    do ii=nStp-1,1,-1
      ntrm(ii,irank,prcpar) = ntrm(ii+1,irank,prcpar)-nn
      if (ntrm(ii,irank,prcpar).le.n1) then
        do jj=1,ii
          ntrm(jj,irank,prcpar) = max(n1,ntrm(ii,irank,prcpar))
          thrs(jj,irank,prcpar) = 0 
        enddo
        exit
      endif
      jj = ntrm(ii,irank,prcpar)
      tt = 1
      tt = (EPSN*abs(coeff(n1,irank)/coeff(jj,irank)))**(tt/(jj-n1))
      thrs(ii,irank,prcpar) = tt
    enddo
!  
  enddo!irank=1,nrank
!  
  if (highestSoFar) call resize( coeff ,2,nmax ,0,rank )
!
!  do ii=lbound(thrs,3),ubound(thrs,3)        !DEBUG
!  do irank=0,rank                            !DEBUG
!    do jj=1,nStp                             !DEBUG
!      xx(jj) = thrs(jj,irank,ii)             !DEBUG
!    enddo                                    !DEBUG
!    write(*,'(i2,99e10.3)') irank,xx(:)      !DEBUG
!    write(*,'(2x,99i10)'  ) ntrm(:,irank,ii) !DEBUG
!  enddo                                      !DEBUG
!  enddo                                      !DEBUG
  end subroutine


  subroutine update_coeff( ncf )
!*******************************************************************
! Coefficients of the expansion of
!   f(n,x) = -int( t^n*log(1-t*x) ,t=0..1 )
! in terms of log(1-x)
!*******************************************************************
  integer ,intent(in) :: ncf
  integer :: ii,jj
  real(kindr2) &  
    :: fact,tt(rank)
!
  call enlarge( coeff ,2,ncf ,0,rank )
!
  do jj=0,rank
  do ii=2,1+jj
    coeff(ii,jj) = 0
  enddo
  enddo
  fact = 1
  do ii=1,rank ;tt(ii)=1 ;enddo
  do ii=2,ncf
    fact = fact*ii
    coeff(ii,0) = (ii-1)/fact
    if (ii.eq.2) cycle
    do jj=1,rank ;tt(jj)=tt(jj)*(jj+1) ;enddo
    coeff(ii,1) = coeff(ii,0)*(1-tt(1))
    if (ii.eq.3) cycle
    coeff(ii,2) = coeff(ii,0)*(1-2*tt(1)+tt(2))
    if (ii.eq.4) cycle
    coeff(ii,3) = coeff(ii,0)*(1-3*tt(1)+3*tt(2)-tt(3))
    if (ii.eq.5) cycle
    coeff(ii,4) = coeff(ii,0)*(1-4*tt(1)+6*tt(2)-4*tt(3)+tt(4))
!   if (ii.eq.n+1) cycle
!   coeff(ii,n) = coeff(ii,0)
!               * ( 1 - binom(n,1)*tt(1) + binom(n,2)*tt(2)...)
  enddo
!
  end subroutine


  function bnlog_c( irank ,xx ) result(rslt)
!*******************************************************************
!*******************************************************************
  integer ,intent(in) :: irank
  complex(kindr2) &   
    ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt,yy,omx
  real(kindr2) &  
    :: aa,rex,imx
  integer :: ii,nn
!
  rex = areal(xx)
  imx = aimag(xx)
!
  if (abs(imx).le.EPSN*abs(rex)) then
    rslt = bnlog_r( irank ,rex ,sgnRe(imx,1) )
    return
  endif
!
  if (abs(xx-1).le.EPSN*10) then
    aa = 1
    rslt = -1
    do ii=2,irank+1
      rslt = rslt - aa/ii
    enddo
    return
  endif
!
  yy = olog(1-1/xx,0)
  aa = abs(yy)
  if     (aa.ge.thrs(6,irank,prcpar)) then
     omx = 1
    rslt = aCoef(irank,irank)
    do ii=irank,1,-1
       omx = 1 + xx*omx
      rslt = aCoef(ii-1,irank) + xx*rslt
    enddo
     omx = (1-xx)*omx
    rslt = omx*yy - rslt/aCoef(irank,irank)
!    if     (irank.eq.0) then
!      rslt = (1-xx)*yy - 1
!    elseif (irank.eq.1) then
!      rslt = (1-xx)*(1+xx)*yy - (1+xx*2)/2
!    elseif (irank.eq.2) then
!      rslt = (1-xx)*(1+xx*(1+xx))*yy - (2+xx*(3+xx*6))/6
!    elseif (irank.eq.3) then
!      rslt = (1-xx)*(1+xx*(1+xx*(1+xx)))*yy &
!           - (3+xx*(4+xx*(6+xx*12)))/12
!    elseif (irank.eq.4) then
!      rslt = (1-xx)*(1+xx*(1+xx*(1+xx*(1+xx))))*yy &
!           - (12+xx*(15+xx*(20+xx*(30+xx*60))))/60
!    endif
    return
  elseif (aa.ge.thrs(5,irank,prcpar)) then ;nn=ntrm(6,irank,prcpar)
  elseif (aa.ge.thrs(4,irank,prcpar)) then ;nn=ntrm(5,irank,prcpar)
  elseif (aa.ge.thrs(3,irank,prcpar)) then ;nn=ntrm(4,irank,prcpar)
  elseif (aa.ge.thrs(2,irank,prcpar)) then ;nn=ntrm(3,irank,prcpar)
  elseif (aa.ge.thrs(1,irank,prcpar)) then ;nn=ntrm(2,irank,prcpar)
                                      else ;nn=ntrm(1,irank,prcpar)
  endif
!
  rslt = coeff(nn,irank)
  do ii=nn-1,2+irank,-1
    rslt = coeff(ii,irank) + yy*rslt
  enddo
  rslt = -(irank+1)*rslt*yy*(yy*xx)**(irank+1)
!
  aa = areal(rslt)
  if (abs(aimag(rslt)).le.EPSN*abs(aa)) rslt = acmplx(aa)
!
  end function


  function bnlog_r( irank ,xx ,sgn ) result(rslt)
!*******************************************************************
!*******************************************************************
  integer ,intent(in) :: irank
  real(kindr2) &  
          ,intent(in) :: xx
  integer ,intent(in) :: sgn
  complex(kindr2) &   
    :: rslt
  real(kindr2) &  
    :: yy,aa,omx
  integer :: ii,nn
  logical :: y_lt_0
!
  if (abs(xx).eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop bnlog_r: ' &
      ,'argument xx=',trim(myprint(xx,8)),', returning 0'
    rslt = 0
    return
  elseif (abs(xx-1).le.EPSN*10) then
    aa = 1
    rslt = -1
    do ii=2,irank+1
      rslt = rslt - aa/ii
    enddo
    return
  endif
!
  yy = 1-1/xx
  y_lt_0 = (yy.lt.RZRO)
  if (y_lt_0) then 
    yy = log(-yy)
    aa = sqrt(yy*yy+ONEPI*ONEPI)
  else
    yy = log( yy)
    aa = abs(yy)
  endif
!
  omx = 1
  do ii=irank,1,-1
    omx = 1+xx*omx
  enddo
  omx = (1-xx)*omx ! (1-x^{rank+1})
!
  if     (aa.ge.thrs(6,irank,prcpar)) then
    rslt = aCoef(irank,irank)
    do ii=irank,1,-1
      rslt = aCoef(ii-1,irank) + xx*rslt
    enddo
    rslt = omx*yy - rslt/aCoef(irank,irank)
!    if     (irank.eq.0) then
!      rslt = omx*yy - 1
!    elseif (irank.eq.1) then
!      rslt = omx*yy - (1+xx*2)/2
!    elseif (irank.eq.2) then
!      rslt = omx*yy - (2+xx*(3+xx*6))/6
!    elseif (irank.eq.3) then
!      rslt = omx*yy - (3+xx*(4+xx*(6+xx*12)))/12
!    elseif (irank.eq.4) then
!      rslt = omx*yy - (12+xx*(15+xx*(20+xx*(30+xx*60))))/60
!    endif
    if (y_lt_0) rslt = rslt + sgn*omx*IPI
    return
  elseif (aa.ge.thrs(5,irank,prcpar)) then ;nn=ntrm(6,irank,prcpar)
  elseif (aa.ge.thrs(4,irank,prcpar)) then ;nn=ntrm(5,irank,prcpar)
  elseif (aa.ge.thrs(3,irank,prcpar)) then ;nn=ntrm(4,irank,prcpar)
  elseif (aa.ge.thrs(2,irank,prcpar)) then ;nn=ntrm(3,irank,prcpar)
  elseif (aa.ge.thrs(1,irank,prcpar)) then ;nn=ntrm(2,irank,prcpar)
                                      else ;nn=ntrm(1,irank,prcpar)
  endif
!
  aa = coeff(nn,irank)
  do ii=nn-1,2+irank,-1
    aa = coeff(ii,irank) + yy*aa
  enddo
  rslt = -(irank+1)*aa*yy*(yy*xx)**(irank+1)
  if (y_lt_0) rslt = rslt + sgn*omx*IPI
!  
  end function

end module


module avh_olo_qp_qmplx
  use avh_olo_units
  use avh_olo_qp_prec
  use avh_olo_qp_auxfun
  use avh_olo_qp_olog
  use avh_olo_qp_dilog

  implicit none
  private
  public :: qmplx_type,qonv,directly,sheet,logc,logc2,logc3,li2c,li2c2
  public :: operator (*) ,operator (/)

  type :: qmplx_type
  complex(kindr2) &   
          :: c
  integer :: p
  end type

  interface qonv
    module procedure qonv_cr,qonv_ci,qonv_c,qonv_i
  end interface

  interface operator (*)
    module procedure prduct_qq,prduct_qr
  end interface
  interface operator (/)
    module procedure ratio_qq,ratio_qr
  end interface

contains


  function qonv_cr(xx,sgn) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz  becomes the
! sign of  sgn .
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  real(kindr2) &  
    ,intent(in) :: sgn
  type(qmplx_type) :: rslt
  real(kindr2) &  
    :: xre,xim
  xre = areal(xx)
  if (xre.ge.RZRO) then
    rslt%c = xx
    rslt%p = 0
  else
    xim = aimag(xx)
    if (xim.eq.RZRO) then
      rslt%c = -xre
      rslt%p = sgnRe(sgn)
    else
      rslt%c = -xx
      rslt%p = sgnRe(xim) ! xim = -Im(rslt%c)
    endif
  endif
  end function

  function qonv_ci(xx,sgn) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz  becomes the
! sign of  sgn .
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  integer         ,intent(in) :: sgn
  type(qmplx_type) :: rslt
  real(kindr2) &  
    :: xre,xim
  xre = areal(xx)
  if (xre.ge.RZRO) then
    rslt%c = xx
    rslt%p = 0
  else
    xim = aimag(xx)
    if (xim.eq.RZRO) then
      rslt%c = -xre
      rslt%p = sign(1,sgn)
    else
      rslt%c = -xx
      rslt%p = sgnRe(xim) ! xim = -Im(rslt%c)
    endif
  endif
  end function

  function qonv_c(xx) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz=1
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  type(qmplx_type) :: rslt
  real(kindr2) &  
    :: xre,xim
  xre = areal(xx)
  if (xre.ge.RZRO) then
    rslt%c = xx
    rslt%p = 0
  else
    xim = aimag(xx)
    if (xim.eq.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop qonv_c: ' &
        ,'negative input with undefined sign for the imaginary part, ' &
        ,'putting +ieps'
      rslt%c = -xre
      rslt%p = 1
    else
      rslt%c = -xx
      rslt%p = sgnRe(xim) ! xim = -Im(rslt%c)
    endif
  endif
  end function

  function qonv_i(xx) result(rslt)
!*******************************************************************
! zz=rslt%c ,iz=rslt%p
! Determine  zz,iz  such that  xx = zz*exp(iz*imag*pi)  and  Re(zz)
! is positive. If  Im(x)=0  and  Re(x)<0  then  iz=1
!*******************************************************************
  integer ,intent(in) :: xx
  type(qmplx_type) :: rslt
  if (xx.ge.0) then
    rslt%c = xx
    rslt%p = 0
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop qonv_i: ' &
      ,'negative input with undefined sign for the imaginary part, ' &
      ,'putting +ieps'
    rslt%c = -xx
    rslt%p = 1
  endif
  end function

  function directly(xx,ix) result(rslt)
!*******************************************************************
!*******************************************************************
  complex(kindr2) &   
    ,intent(in) :: xx
  integer         ,intent(in) :: ix
  type(qmplx_type) :: rslt
  rslt%c = xx
  rslt%p = ix
  end function


  function sheet(xx) result(ii)
!*******************************************************************
! Returns the number of the Riemann-sheet (times 2) for the complex
! number  xx*exp(ix*imag*pi) . The real part of xx is assumed to be
! positive or zero. Examples:
! xx=1+imag, ix=-1 -> ii= 0 
! xx=1+imag, ix= 1 -> ii= 2 
! xx=1-imag, ix=-1 -> ii=-2 
! xx=1-imag, ix= 1 -> ii= 0 
! xx=1     , ix= 1 -> ii= 0  convention that log(-1)=pi on
! xx=1     , ix=-1 -> ii=-2  the principal Riemann-sheet
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  integer :: ii,jj
  real(kindr2) &  
    :: xim
  jj = mod(xx%p,2)
  ii = xx%p-jj
  xim = aimag(xx%c)
  if (xim.le.RZRO) then ! also xim=0 <==> log(-1)=pi, not -pi
    if (jj.eq.-1) ii = ii-2
  else
    if (jj.eq. 1) ii = ii+2
  endif
  end function


  function prduct_qq(yy,xx) result(zz)
!*******************************************************************
! Return the product  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy,xx
  type(qmplx_type) :: zz
  zz%c = yy%c*xx%c
  zz%p = yy%p+xx%p
  if (areal(zz%c).lt.RZRO) then
    zz%p = zz%p + sgnIm(xx%c)
    zz%c = -zz%c
  endif
  end function

  function prduct_qr(yy,xx) result(zz)
!*******************************************************************
! Return the product  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy
  real(kindr2) &  
    ,intent(in) :: xx
  type(qmplx_type) :: zz
  zz%c = yy%c*abs(xx)
  zz%p = yy%p
  end function

  function ratio_qq(yy,xx) result(zz)
!*******************************************************************
! Return the ratio  zz  of  yy  and  xx  
! keeping track of (the multiple of pi of) the phase %p such that
! the real part of  zz%c  remains positive 
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy,xx
  type(qmplx_type) :: zz
  zz%c = yy%c/xx%c
  zz%p = yy%p-xx%p
  if (areal(zz%c).lt.RZRO) then
    zz%p = zz%p - sgnIm(xx%c)
    zz%c = -zz%c
  endif
  end function

  function ratio_qr(yy,xx) result(zz)
!*******************************************************************
!*******************************************************************
  type(qmplx_type) ,intent(in) :: yy
  real(kindr2) &  
    ,intent(in) :: xx
  type(qmplx_type) :: zz
  zz%c = yy%c/abs(xx)
  zz%p = yy%p
  end function


  function logc(xx) result(rslt)
!*******************************************************************
! log(xx)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = olog(acmplx(xx%c),xx%p)
  rslt = olog(xx%c,xx%p)
  end function

  function logc2(xx) result(rslt)
!*******************************************************************
! log(xx)/(1-xx)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = -olog1(acmplx(xx%c),xx%p)
  rslt = -olog1(xx%c,xx%p)
  end function

  function logc3(xx) result(rslt)
!*******************************************************************
!  ( log(xx)/(1-xx) + 1 )/(1-xx)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = olog2(acmplx(xx%c),xx%p)
  rslt = olog2(xx%c,xx%p)
  end function

  function li2c(xx) result(rslt)
!*******************************************************************
!    /1    ln(1-(1-xx)*t)
!  - |  dt -------------- 
!    /0        t
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx
  complex(kindr2) &   
    :: rslt
!  rslt = dilog(acmplx(xx%c),xx%p)
  rslt = dilog(xx%c,xx%p)
  end function

  function li2c2(xx,yy) result(rslt)
!*******************************************************************
! ( li2(xx) - li2(yy) )/(xx-yy)
!*******************************************************************
  type(qmplx_type) ,intent(in) :: xx,yy
  complex(kindr2) &   
    :: rslt
!  rslt = dilog( acmplx(xx%c),xx%p ,acmplx(yy%c),yy%p )
!  write(*,*) 'li2c2 x:',xx%c,xx%p !DEBUG
!  write(*,*) 'li2c2 y:',yy%c,yy%p !DEBUG
  rslt = dilog( xx%c,xx%p ,yy%c,yy%p )
!  write(*,*) 'li2c2 out:',rslt !DEBUG
  end function


end module


module avh_olo_qp_bub
  use avh_olo_units
  use avh_olo_qp_prec
  use avh_olo_qp_auxfun
  use avh_olo_qp_bnlog
  use avh_olo_qp_qmplx
  use avh_olo_qp_olog
  implicit none
  private
  public :: tadp ,tadpn ,bub0 ,dbub0 ,bub1 ,bub11 ,bub111 ,bub1111

contains

  subroutine tadp( rslt ,mm ,amm ,rmu2 )
!*******************************************************************
! The 1-loop scalar 1-point function.
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
    ,intent(in)  :: amm,rmu2
!
  rslt(2) = 0
  if (amm.eq.RZRO.or.mm.eq.CZRO) then
    rslt(1) = 0
    rslt(0) = 0
  else
    rslt(1) = mm
    rslt(0) = mm - mm*logc( qonv(mm/rmu2,-1) )
  endif
  end subroutine


  subroutine tadpn( rslt ,rank ,mm ,amm ,rmu2 )
!*******************************************************************
! The 1-loop tensor 1-point functions.
!   rslt(:,0) = A0
!   rslt(:,1) = A00
!   rslt(:,2) = A0000  etc.
! For input  rank  only  rslt(:,0:rank/2)  is filled.
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
    ,intent(in)  :: amm,rmu2
  integer ,intent(in) :: rank
  complex(kindr2) &   
    :: aa
  real(kindr2) &  
    :: bb
  integer :: ii
!
  do ii=0,rank
    rslt(2,ii) = 0
    rslt(1,ii) = 0
    rslt(0,ii) = 0
  enddo
  if (amm.eq.RZRO.or.mm.eq.CZRO) then
    return
  else
    rslt(1,0) = mm
    rslt(0,0) = mm - mm*logc( qonv(mm/rmu2,-1) )
    aa = 1
    bb = 0
    do ii=1,rank/2
      aa = aa*mm/(2*(ii+1))
      bb = bb + RONE/(ii+1)
      rslt(1,ii) = aa*( rslt(1,0) )
      rslt(0,ii) = aa*( rslt(0,0) + mm*bb )
    enddo
  endif
  end subroutine


!*******************************************************************
! Return the Passarino-Veltman functions
!
!      C   /      d^(Dim)q
!   ------ | -------------------- = b0
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!      C   /    d^(Dim)q q^mu
!   ------ | -------------------- = p^mu b1
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!      C   /  d^(Dim)q q^mu q^nu
!   ------ | -------------------- = g^{mu,nu} b00 + p^mu p^nu b11
!   i*pi^2 / [q^2-m0][(q+p)^2-m1]
!
!   etc.
!
! Based on the formulas from
! A. Denner, M. Dittmaier, Nucl.Phys. B734 (2006) 62-115
!*******************************************************************

  subroutine bub0( b0 &
                  ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
  endif
!
  b0(2) = 0
  b0(1) = CONE
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
  endif
!
  end subroutine

  subroutine bub1( b1,b0 &
                  ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  end subroutine

  subroutine bub11( b11,b00,b1,b0 &
                   ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      b00(0)=0 ;b00(1)=0 ;b00(2)=0 
      b11(0)=0 ;b11(1)=0 ;b11(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
  b11(2) = 0      
  b11(1) = CONE/3 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
      b11(0) = lna/3 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
      b11(0) = ( lna - bnlog(2,x1) )/3 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
      b11(0) = ( lna*3 + 2 )/9 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
      b11(0) = ( (lna-bnlog(2,x1))*3 + 1 )/9 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    b11(0) = ( lna - bnlog(2,x1) - bnlog(2,x2) )/3 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  rank = 0 
  call tadpn( a0 ,rank ,m1 ,am1 ,rmu2 )
  x1 = (m1-m0)-pp
  x2 = 2*m0
  b00(2) = 0
  b00(1) = ( a0(1,0) - x1*b1(1) + x2*b0(1) )/6
  b00(0) = ( a0(0,0) - x1*b1(0) + x2*b0(0) + 4*b00(1) )/6
  end subroutine

  subroutine bub111( b111,b001,b11,b00,b1,b0 &
                    ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b111(0:2),b001(0:2),b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      b00(0)=0 ;b00(1)=0 ;b00(2)=0 
      b11(0)=0 ;b11(1)=0 ;b11(2)=0 
      b001(0)=0 ;b001(1)=0 ;b001(2)=0 
      b111(0)=0 ;b111(1)=0 ;b111(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
  b11(2) = 0      
  b11(1) = CONE/3 
  b111(2) = 0      
  b111(1) =-CONE/4 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
      b11(0) = lna/3 
      b111(0) =-lna/4 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
      b11(0) = ( lna - bnlog(2,x1) )/3 
      b111(0) =-( lna - bnlog(3,x1) )/4 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
      b11(0) = ( lna*3 + 2 )/9 
      b111(0) =-( lna*4 + 2 )/16 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
      b11(0) = ( (lna-bnlog(2,x1))*3 + 1 )/9 
      b111(0) =-( (lna-bnlog(3,x1))*4 + 1 )/16 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    b11(0) = ( lna - bnlog(2,x1) - bnlog(2,x2) )/3 
    b111(0) =-( lna - bnlog(3,x1) - bnlog(3,x2) )/4 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  rank = 0 
  rank = 1 
  call tadpn( a0 ,rank ,m1 ,am1 ,rmu2 )
  x1 = (m1-m0)-pp
  x2 = 2*m0
  b00(2) = 0
  b00(1) = ( a0(1,0) - x1*b1(1) + x2*b0(1) )/6
  b00(0) = ( a0(0,0) - x1*b1(0) + x2*b0(0) + 4*b00(1) )/6
  b001(2) = 0
  b001(1) = (-a0(1,0) - x1*b11(1) + x2*b1(1) )/8
  b001(0) = (-a0(0,0) - x1*b11(0) + x2*b1(0) + 4*b001(1) )/8
  end subroutine

  subroutine bub1111( b1111,b0011,b0000,b111,b001,b11,b00,b1,b0 &
                    ,pp,m0i,m1i ,app,am0i,am1i ,rmu2 )
  complex(kindr2) &   
    ,intent(out) :: b1111(0:2),b0011(0:2),b0000(0:2) &
                   ,b111(0:2),b001(0:2),b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp,m0i,m1i
  real(kindr2) &  
    ,intent(in)  :: app,am0i,am1i,rmu2
  complex(kindr2) &   
    :: m0,m1,a0(0:2,0:1),lna,x1,x2,lambda
  real(kindr2) &  
    :: am0,am1,maxm
  logical :: switch 
  integer :: rank
!
  maxm = max(am0i,am1i)
  if (maxm.eq.RZRO) then
    if (app.eq.RZRO) then
      b0(0)=0 ;b0(1)=0 ;b0(2)=0
      b1(0)=0 ;b1(1)=0 ;b1(2)=0 
      b00(0)=0 ;b00(1)=0 ;b00(2)=0 
      b11(0)=0 ;b11(1)=0 ;b11(2)=0 
      b001(0)=0 ;b001(1)=0 ;b001(2)=0 
      b111(0)=0 ;b111(1)=0 ;b111(2)=0 
      b0000(0)=0 ;b0000(1)=0 ;b0000(2)=0 
      b0011(0)=0 ;b0011(1)=0 ;b0011(2)=0 
      b1111(0)=0 ;b1111(1)=0 ;b1111(2)=0 
      return
    endif
  endif
!
  if (am1i.ge.maxm) then
    m0=m0i ;am0=am0i
    m1=m1i ;am1=am1i
    switch = .false. 
  else
    m0=m1i ;am0=am1i
    m1=m0i ;am1=am0i
    switch = .true. 
  endif
!
  b0(2) = 0
  b0(1) = CONE
  b1(2) = 0      
  b1(1) =-CONE/2 
  b11(2) = 0      
  b11(1) = CONE/3 
  b111(2) = 0      
  b111(1) =-CONE/4 
  b1111(2) = 0      
  b1111(1) = CONE/5 
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = lna
      b1(0) =-lna/2 
      b11(0) = lna/3 
      b111(0) =-lna/4 
      b1111(0) = lna/5 
    else
      lna = -logc(qonv(m1/rmu2,-1))
      x1 = (m1-am1*IEPS)/(m1-m0)
      b0(0) =   lna - bnlog(0,x1)
      b1(0) =-( lna - bnlog(1,x1) )/2 
      b11(0) = ( lna - bnlog(2,x1) )/3 
      b111(0) =-( lna - bnlog(3,x1) )/4 
      b1111(0) = ( lna - bnlog(4,x1) )/5 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
    else
      b1111(0) =-b1111(0)-4*b111(0)-6*b11(0)-4*b1(0)-b0(0) 
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  elseif (am0.eq.RZRO) then
    if (abs(pp-m1).le.am1*EPSN*10) then
      lna = -logc(qonv(m1/rmu2,-1))
      b0(0) = ( lna   + 2 )
      b1(0) =-( lna*2 + 2 )/4 
      b11(0) = ( lna*3 + 2 )/9 
      b111(0) =-( lna*4 + 2 )/16 
      b1111(0) = ( lna*5 + 2 )/25 
    else
      lna = -logc(qonv((m1-pp)/rmu2,-1))
      x1  = (pp-m1+am1*IEPS)/pp
      b0(0) = ( lna-bnlog(0,x1) + 1 )
      b1(0) =-( (lna-bnlog(1,x1))*2 + 1 )/4 
      b11(0) = ( (lna-bnlog(2,x1))*3 + 1 )/9 
      b111(0) =-( (lna-bnlog(3,x1))*4 + 1 )/16 
      b1111(0) = ( (lna-bnlog(4,x1))*5 + 1 )/25 
    endif
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1111(0) =-b1111(0)-4*b111(0)-6*b11(0)-4*b1(0)-b0(0) 
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  else
    lna = -logc(qonv(m0/rmu2,-1))
    call solabc( x1,x2 ,lambda ,pp ,(m1-m0)-pp ,m0-am0*IEPS ,0 )
    b0(0) = ( lna - bnlog(0,x1) - bnlog(0,x2) ) 
    b1(0) =-( lna - bnlog(1,x1) - bnlog(1,x2) )/2 
    b11(0) = ( lna - bnlog(2,x1) - bnlog(2,x2) )/3 
    b111(0) =-( lna - bnlog(3,x1) - bnlog(3,x2) )/4 
    b1111(0) = ( lna - bnlog(4,x1) - bnlog(4,x2) )/5 
    if (switch) then
      x2=m0;m0=m1;m1=x2
      b1111(0) =-b1111(0)-4*b111(0)-6*b11(0)-4*b1(0)-b0(0) 
      b111(0) =-b111(0)-3*b11(0)-3*b1(0)-b0(0) 
      b11(0) = b11(0)+2*b1(0)+b0(0) 
      b1(0) =-b0(0)-b1(0)
    endif
  endif
!
  rank = 0 
  rank = 1 
  rank = 2 
  call tadpn( a0 ,rank ,m1 ,am1 ,rmu2 )
  x1 = (m1-m0)-pp
  x2 = 2*m0
  b00(2) = 0
  b00(1) = ( a0(1,0) - x1*b1(1) + x2*b0(1) )/6
  b00(0) = ( a0(0,0) - x1*b1(0) + x2*b0(0) + 4*b00(1) )/6
  b001(2) = 0
  b001(1) = (-a0(1,0) - x1*b11(1) + x2*b1(1) )/8
  b001(0) = (-a0(0,0) - x1*b11(0) + x2*b1(0) + 4*b001(1) )/8
  b0000(2) = 0
  b0000(1) = ( a0(1,1) - x1*b001(1) + x2*b00(1) )/10
  b0000(0) = ( a0(0,1) - x1*b001(0) + x2*b00(0) + 4*b0000(1) )/10
  b0011(2) = 0
  b0011(1) = ( a0(1,0) - x1*b111(1) + x2*b11(1) )/10
  b0011(0) = ( a0(0,0) - x1*b111(0) + x2*b11(0) + 4*b0011(1) )/10
  end subroutine


!*******************************************************************
! Derivative of B0
! expects  m0<m1
! only finite case, so input must not be  m0=0 & m1=pp
!*******************************************************************

  subroutine dbub0( rslt &
                   ,pp,m0,m1 ,app,am0,am1 )
  complex(kindr2) &   
    ,intent(out) :: rslt
  complex(kindr2) &   
    ,intent(in)  :: pp,m0,m1
  real(kindr2) &  
    ,intent(in)  :: app,am0,am1
  complex(kindr2) &   
    :: ch,x1,x2,lambda
  real(kindr2) &  
    :: ax1,ax2,ax1x2,maxa
  type(qmplx_type) :: q1,q2,q1o,q2o
  integer :: sgn
!
  if (am1.eq.RZRO) then
    if (app.eq.RZRO) then
      rslt = 0
      return
    endif
  endif
!
  if (app.eq.RZRO) then
    if (abs(m0-m1).le.am1*EPSN*10) then
      rslt = 1/(6*m1)
    else
      ch = m0/m1
      rslt = ( CONE/2 - ch*olog3(ch,0) )/m1 
    endif
  elseif (am1.eq.RZRO) then
    rslt =-1/pp
  else
    call solabc( x1,x2 ,lambda ,pp ,(m0-m1)-pp ,m1 ,0 )
    sgn =-sgnRe(pp)*sgnRe(x2-x1)
    q1  = qonv(x1  , sgn)
    q1o = qonv(x1-1, sgn)
    q2  = qonv(x2  ,-sgn)
    q2o = qonv(x2-1,-sgn)
    ax1 = abs(x1)
    ax2 = abs(x2)
    ax1x2 = abs(x1-x2)
    maxa = max(ax1,ax2)
    if (ax1x2.lt.maxa*EPSN*10) then
      rslt = ( (x1+x2-1)*logc(q2/q2o) - 2 )/pp
    elseif (ax1x2*2.lt.maxa) then
      if     (x1.eq.CZRO.or.x1.eq.CONE) then
        rslt = ( (x1+x2-1)*logc(q2/q2o) - 1 )/pp
      elseif (x2.eq.CZRO.or.x2.eq.CONE) then
        rslt = ( (x1+x2-1)*logc(q1/q1o) - 1 )/pp
      else
        rslt = x1*(x1-1)*( logc2(q1o/q2o)/(x2-1) - logc2(q1/q2)/x2 ) &
             + (x1+x2-1)*logc(q2/q2o) - 1
        rslt = rslt/pp
      endif
    else
      rslt = 0
      if (ax1.ne.RZRO) then
        if (ax1.lt.2*RONE) then
          rslt = rslt - x1
          if (x1.ne.CONE) rslt = rslt - x1*logc2(q1/q1o)
        else
          rslt = rslt + x1/(x1-1)*logc3(q1/q1o)
        endif
      endif
      if (ax2.ne.RZRO) then
        if (ax2.lt.2*RONE) then
          rslt = rslt + x2
          if (x2.ne.CONE) rslt = rslt + x2*logc2(q2/q2o)
        else
          rslt = rslt - x2/(x2-1)*logc3(q2/q2o)
        endif
      endif
      rslt = rslt/lambda
    endif
  endif
!
  end subroutine


end module


module avh_olo_qp_tri
  use avh_olo_units
  use avh_olo_qp_prec
  use avh_olo_qp_auxfun
  use avh_olo_qp_qmplx
  implicit none
  private
  public :: tria0,tria1,tria2,tria3,tria4,trif0,trif1,trif2,trif3 &
           ,trif3HV &
           ,permtable,casetable,base
  integer ,parameter :: permtable(3,0:7)=reshape((/ &
       1,2,3 &! 0, 0 masses non-zero, no permutation
      ,1,2,3 &! 1, 1 mass non-zero,   no permutation
      ,3,1,2 &! 2, 1 mass non-zero,   1 cyclic permutation
      ,1,2,3 &! 3, 2 masses non-zero, no permutation
      ,2,3,1 &! 4, 1 mass non-zero,   2 cyclic permutations
      ,2,3,1 &! 5, 2 masses non-zero, 2 cyclic permutations
      ,3,1,2 &! 6, 2 masses non-zero, 1 cyclic permutation
      ,1,2,3 &! 7, 3 masses non-zero, no permutation
      /) ,(/3,8/))                     ! 0,1,2,3,4,5,6,7
  integer ,parameter :: casetable(0:7)=(/0,1,1,2,1,2,2,3/)
  integer ,parameter :: base(3)=(/4,2,1/)

contains

   subroutine tria4( rslt ,cpp,cm2,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /             d^(Dim)q
!            ------ | ----------------------------------
!            i*pi^2 / q^2 [(q+k1)^2-m2] [(q+k1+k2)^2-m3]
!
! with  k1^2=m2, k2^2=pp, (k1+k2)^2=m3.
! m2,m3 should NOT be identically 0d0.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cm2,cm3,cpp
  real(kindr2) &  
     ,intent(in)  :: rmu2
   type(qmplx_type) :: q23,qm3,q32
  complex(kindr2) &   
     :: sm2,sm3,k23,r23,d23,cc
!
   sm2 = mysqrt(cm2)
   sm3 = mysqrt(cm3)
   k23 = (cm2+cm3-cpp)/(sm2*sm3)
   call rfun( r23,d23, k23 )
   if (r23.eq.-CONE) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop tria4: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   q23 = qonv(r23,-1)
   qm3 = qonv(cm3/rmu2,-1)
   q32 = qonv(sm3)/qonv(sm2)
!
   rslt(2) = 0
   cc = logc2(q23) * r23/(1+r23)/(sm2*sm3)
   rslt(1) = -cc
   rslt(0) = cc*( logc(qm3) - logc(q23) ) &
           - li2c2(q32*q23,q32/q23) / cm2 &
           + li2c2(q23*q23,qonv(1)) * r23/(sm2*sm3)
   end subroutine


   subroutine tria3( rslt ,cp2,cp3,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  p2=k2^2, p3=(k1+k2)^2.
! mm should NOT be identically 0d0,
! and p2 NOR p3 should be identical to mm. 
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp2,cp3,cm3
  real(kindr2) &  
     ,intent(in)  :: rmu2
   type(qmplx_type) :: q13,q23,qm3,x1,x2
  complex(kindr2) &   
     :: r13,r23
!
   r13 = cm3-cp3
   r23 = cm3-cp2
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   qm3 = qonv(cm3,-1)
   x1 = q23/qm3
   x2 = q13/qm3
   rslt(2) = 0
   rslt(1) = -logc2( q23/q13 )/r13
   rslt(0) = -li2c2( x1,x2 )/cm3 &
           - rslt(1)*( logc(x1*x2)+logc(qm3/rmu2) )
   end subroutine


   subroutine tria2( rslt ,cp3,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  k1^2 = 0 , k2^2 = m3  and  (k1+k2)^2 = p3.
! mm should NOT be identically 0d0,
! and pp should NOT be identical to mm. 
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cm3
  real(kindr2) &  
     ,intent(in)  :: rmu2
   type(qmplx_type) :: q13,qm3,qxx
  complex(kindr2) &   
     :: r13,logm,z2,z1,z0,cc
!
   r13 = cm3-cp3
   q13 = qonv(r13,-1)
   qm3 = qonv(cm3,-1)
   logm = logc( qm3/rmu2 )
   qxx = qm3/q13
   z2 = 1 
   z2 = z2/2
   z1 = logc(qxx)
   z0 = PISQo24 + z1*z1/2 - li2c(qxx)
   cc = -1/r13
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine tria1( rslt ,cm3 ,rmu2 )
!*******************************************************************
! calculates
!               C   /          d^(Dim)q
!            ------ | -----------------------------
!            i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3]
!
! with  k1^2 = (k1+k2)^2 = m3.
! mm should NOT be identically 0d0.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cm3
  real(kindr2) &  
     ,intent(in)  :: rmu2
  complex(kindr2) &   
     :: zm
!
   zm = 1/(2*cm3)
   rslt(2) = 0
   rslt(1) = -zm
   rslt(0) = zm*( 2 + logc(qonv(cm3/rmu2,-1)) )
   end subroutine


   subroutine tria0( rslt ,cp ,ap ,rmu2 )
!*******************************************************************
! calculates
!               C   /         d^(Dim)q
!            ------ | ------------------------
!            i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  p1 = k1^2,  p2 = k2^2,  p3 = k3^2
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! If any of these numbers is IDENTICALLY 0d0, the corresponding
! IR-singular case is returned.
!*******************************************************************
   use avh_olo_qp_olog
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp(3)
  real(kindr2) &  
     ,intent(in)  :: ap(3),rmu2
  real(kindr2) &  
     :: pp(3),rp1,rp2,rp3
  complex(kindr2) &   
     :: log2,log3
   integer :: icase,i1,i2,i3
!
   pp(1)=areal(cp(1))
   pp(2)=areal(cp(2))
   pp(3)=areal(cp(3))
!
   icase = 0
   if (ap(1).gt.RZRO) icase = icase + base(1)
   if (ap(2).gt.RZRO) icase = icase + base(2)
   if (ap(3).gt.RZRO) icase = icase + base(3)
   rp1 = pp(permtable(1,icase))
   rp2 = pp(permtable(2,icase))
   rp3 = pp(permtable(3,icase))
   icase  = casetable(  icase)
!
   i1=0 ;if (-rp1.lt.RZRO) i1=-1
   i2=0 ;if (-rp2.lt.RZRO) i2=-1
   i3=0 ;if (-rp3.lt.RZRO) i3=-1
!
   if     (icase.eq.0) then
! 0 masses non-zero
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop tria0: ' &
       ,'all external masses equal zero, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
   elseif (icase.eq.1) then
! 1 mass non-zero
    log3 = olog( abs(rp3/rmu2) ,i3 )
    rslt(2) = 1/rp3
    rslt(1) = -log3/rp3
    rslt(0) = ( log3**2/2 - 2*PISQo24 )/rp3
  elseif (icase.eq.2) then
! 2 masses non-zero
    log2 = olog( abs(rp2/rmu2) ,i2 )
    log3 = olog( abs(rp3/rmu2) ,i3 )
    rslt(2) = 0
    rslt(1) = -olog1( abs(rp3/rp2) ,i3-i2 )/rp2
    rslt(0) = -rslt(1)*(log3+log2)/2
  elseif (icase.eq.3) then
! 3 masses non-zero
    call trif0( rslt ,cp(1),cp(2),cp(3) )
  endif
  end subroutine


   subroutine trif0( rslt ,p1,p2,p3 )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! equal zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1,p2,p3
   type(qmplx_type) :: q23,q24,q34,qx1,qx2
  complex(kindr2) &   
     :: r23,r24,r34,aa,bb,cc,dd,x1,x2
  real(kindr2) &  
     :: hh
!
   r23 = -p1
   r24 = -p3
   r34 = -p2
!
   aa = r34*r24
   bb = r24 + r34 - r23
   cc = 1
   hh = areal(r23)
   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,1 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1, hh)
   qx2 = qonv(x2,-hh)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   rslt(0) = li2c2( qx1*q34 ,qx2*q34 )*r34 &
           + li2c2( qx1*q24 ,qx2*q24 )*r24 &
           - logc2( qx1/qx2 )*logc( qx1*qx2 )/(x2*2) &
           - logc2( qx1/qx2 )*logc( q23 )/x2
!
   rslt(0) = rslt(0)/aa
   end subroutine


   subroutine trif1( rslt ,p1i,p2i,p3i ,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with one internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1i,p2i,p3i ,m3i 
   type(qmplx_type) :: q23,q24,q34,qm4,qx1,qx2,qss
  complex(kindr2) &   
     :: p2,p3,p4,p12,p23,m4,sm2,sm3,sm4 &
                     ,aa,bb,cc,dd,x1,x2,r23,r24,r34
  real(kindr2) &  
     :: mhh
   logical :: r24Not0,r34Not0
!
!   p1 = nul
   p2 = p1i
   p3 = p2i
   p4 = p3i
   p12 = p1i
   p23 = p3i
!   m1 = infinite
!   m2 = m1i = 0
!   m3 = m2i = 0
   m4 = m3i
!
   sm4 = mysqrt(m4)
   mhh = abs(sm4)
   sm3 = mhh
   sm2 = sm3
!
   r23 = (   -p2 -p2 *IEPS )/(sm2*sm3)
   r24 = ( m4-p23-p23*IEPS )/(sm2*sm4)
   r34 = ( m4-p3 -p3 *IEPS )/(sm3*sm4)
!
   r24Not0 = (abs(areal(r24))+abs(aimag(r24)).ge.neglig(prcpar))
   r34Not0 = (abs(areal(r34))+abs(aimag(r34)).ge.neglig(prcpar))
!
   aa = r34*r24 - r23
!
   if (aa.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif1: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   bb = r24/sm3 + r34/sm2 - r23/sm4
   cc = 1/(sm2*sm3)
!   hh = areal(r23)
!   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,0 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1) ! x2 SHOULD HAVE im. part
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
   qm4 = qonv(sm4,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   rslt(0) = -logc2( qx1/qx2 )*logc( qx1*qx2/(qm4*qm4) )/(x2*2) &
             -li2c2( qx1*qm4 ,qx2*qm4 )*sm4
!
   if (r34Not0) then
     qss = q34*mhh
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34*sm3
   endif
!
   if (r24Not0) then
     qss = q24*mhh
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r24*sm2
   endif
!
   rslt(0) = rslt(0) - logc2( qx1/qx2 )*logc( q23*(mhh*mhh) )/x2
!
   rslt(0) = rslt(0)/(aa*sm2*sm3*sm4)
   end subroutine


   subroutine trif2( rslt ,p1i,p2i,p3i ,m2i,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with two internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1i,p2i,p3i ,m2i,m3i
   type(qmplx_type) :: q23,q34,q24,qm2,qm3,qm4,qx1,qx2,qss,qy1,qy2
  complex(kindr2) &   
     :: p2,p3,p23,m2,m4,sm2,sm3,sm4,aa,bb,cc,dd,x1,x2 &
                     ,r23,k24,r34,r24,d24
   logical :: r23Not0,r34Not0
!
!   p1 = nul
   p2 = p3i
   p3 = p1i
!   p4 = p2i
!   p12 = p3i
   p23 = p2i
!   m1 = infinite
   m2 = m3i
!   m3 = m1i = 0
   m4 = m2i
!
!   sm1 = infinite
   sm2 = mysqrt(m2)
   sm3 = abs(sm2) !mysqrt(m3)
   sm4 = mysqrt(m4)
!
   r23 = (    m2-p2 -p2 *IEPS )/(sm2*sm3) ! p2
   k24 = ( m2+m4-p23-p23*IEPS )/(sm2*sm4) ! p2+p3
   r34 = (    m4-p3 -p3 *IEPS )/(sm3*sm4) ! p3
!
   r23Not0 = (abs(areal(r23))+abs(aimag(r23)).ge.neglig(prcpar))
   r34Not0 = (abs(areal(r34))+abs(aimag(r34)).ge.neglig(prcpar))
!
   call rfun( r24,d24 ,k24 )
!
   aa = r34/r24 - r23
!
   if (aa.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif2: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   bb = -d24/sm3 + r34/sm2 - r23/sm4
   cc = (sm4/sm2 - r24)/(sm3*sm4)
!   hh = areal(r23 - r24*r34)
!   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc(x1,x2,dd ,aa,bb,cc ,0)
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1 ) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1 ) ! x2 SHOULD HAVE im. part
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
   qm2 = qonv(sm2,-1)
   qm3 = qonv(sm3,-1)
   qm4 = qonv(sm4,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   qy1 = qx1/q24
   qy2 = qx2/q24
!
   rslt(0) = li2c2( qy1*qm2 ,qy2*qm2 )/r24*sm2
!
   if (x2.ne.CZRO) then ! better to put a threshold on cc 
     rslt(0) = rslt(0) + ( logc2( qy1/qy2 )*logc( qy1*qy2/(qm2*qm2) ) &
                          -logc2( qx1/qx2 )*logc( qx1*qx2/(qm4*qm4) ) )/(x2*2)
   endif
!
   rslt(0) = rslt(0) - li2c2( qx1*qm4 ,qx2*qm4 )*sm4
!
   if (r23Not0) then
     qss = q23*qm3/q24
     rslt(0) = rslt(0) - li2c2( qx1*qss ,qx2*qss )*r23*sm3/r24
   endif
!
   if (r34Not0) then
     qss = q34*qm3
     rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34*sm3
   endif
!
   rslt(0) = rslt(0)/(aa*sm2*sm3*sm4)
   end subroutine


   subroutine trif3( rslt ,p1i,p2i,p3i ,m1i,m2i,m3i )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! non-zero. Obtained from the formulas for 4-point functions in
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
! by sending one internal mass to infinity.
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p1i,p2i,p3i,m1i,m2i,m3i
   type(qmplx_type) :: q12,q13,q23,qm1,qm2,qm3,qx1,qx2,qz1,qz2,qtt
  complex(kindr2) &   
     :: p1,p2,p3,m1,m2,m3,sm1,sm2,sm3,aa,bb,cc,dd,x1,x2 &
                     ,k12,k13,k23,r12,r13,r23,d12,d13,d23 
  real(kindr2) &  
     :: h1,h2,h3
!
   h1 = -aimag(m1i)
   h2 = -aimag(m2i)
   h3 = -aimag(m3i)
   if (h2.ge.h1.and.h2.ge.h3) then
     p1=p3i ;p2=p1i ;p3=p2i ;m1=m3i ;m2=m1i ;m3=m2i
   else
     p1=p1i ;p2=p2i ;p3=p3i ;m1=m1i ;m2=m2i ;m3=m3i
   endif
!
   sm1 = mysqrt(m1)
   sm2 = mysqrt(m2)
   sm3 = mysqrt(m3)
!
   k12 = 0
   k13 = 0
   k23 = 0
   if (m1+m2.ne.p1) k12 = ( m1+m2-p1-p1*IEPS )/(sm1*sm2) ! p1
   if (m1+m3.ne.p3) k13 = ( m1+m3-p3-p3*IEPS )/(sm1*sm3) ! p1+p2 => p12
   if (m2+m3.ne.p2) k23 = ( m2+m3-p2-p2*IEPS )/(sm2*sm3) ! p2
!
   call rfun( r12,d12 ,k12 )
   call rfun( r13,d13 ,k13 )
   call rfun( r23,d23 ,k23 )
!
   aa = sm2/sm3 - k23 + r13*(k12 - sm2/sm1)
!
   if (aa.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   bb = d13/sm2 + k12/sm3 - k23/sm1
   cc = ( sm1/sm3 - 1/r13 )/(sm1*sm2)
!   hh = areal( (r13-sm1/sm3)/(sm1*sm2) )
!   dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
   call solabc( x1,x2,dd ,aa,bb,cc ,0 )
   x1 = -x1
   x2 = -x2
!
   qx1 = qonv(x1 ,1) ! x1 SHOULD HAVE im. part
   qx2 = qonv(x2 ,1) ! x2 SHOULD HAVE im. part
   q12 = qonv(r12,-1)
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   qm1 = qonv(sm1,-1)
   qm2 = qonv(sm2,-1)
   qm3 = qonv(sm3,-1)
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   qz1 = qx1*qm2
   qz2 = qx2*qm2
   rslt(0) = rslt(0) + ( li2c2( qz1*q12 ,qz2*q12 )*r12 &
                        +li2c2( qz1/q12 ,qz2/q12 )/r12 )*sm2
   qtt = q13*qm2
   qz1 = qx1*qtt
   qz2 = qx2*qtt
   rslt(0) = rslt(0) - ( li2c2( qz1*q23 ,qz2*q23 )*r23 &
                        +li2c2( qz1/q23 ,qz2/q23 )/r23 )*r13*sm2
   qz1 = qx1*q13
   qz2 = qx2*q13
   rslt(0) = rslt(0) + li2c2( qz1*qm3 ,qz2*qm3 )*r13*sm3 &
                     - li2c2( qx1*qm1 ,qx2*qm1 )*sm1
   if (x2.ne.CZRO) then
     rslt(0) = rslt(0) + ( logc2( qz1/qz2 )*logc( qz1*qz2/(qm3*qm3) ) &
                          -logc2( qx1/qx2 )*logc( qx1*qx2/(qm1*qm1) ) )/(x2*2)
   endif
!
   rslt(0) = rslt(0)/(aa*sm1*sm2*sm3)
   end subroutine
   

   subroutine trif3HV( rslt ,pp,mm ,ap ,smax ,lam )
!*******************************************************************
! Finite 1-loop scalar 3-point function with all internal masses
! non-zero. Based on the fomula of 't Hooft & Veltman
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: pp(3),mm(3)
  real(kindr2) &  
     ,intent(in)  :: ap(3),smax
  complex(kindr2) &   
     ,optional ,intent(in) :: lam
  complex(kindr2) &   
     :: p1,p2,p3,m1,m2,m3,slam,yy
  complex(kindr2) &   
     :: sm1,sm2,sm3
   type(qmplx_type) :: qm1,qm2,qm3
  real(kindr2) &  
     :: a12,a23,a31,thrs,a1,a2,a3
!
! Order squared momenta, first one smallest
   if     (ap(1).le.ap(2).and.ap(1).le.ap(3)) then
     if (ap(2).le.ap(3)) then
       a1=ap(1) ;a2=ap(2) ;a3=ap(3)
       p1=pp(1) ;p2=pp(2) ;p3=pp(3)
       m1=mm(1) ;m2=mm(2) ;m3=mm(3)
     else
       a1=ap(1) ;a2=ap(3) ;a3=ap(2)
       p1=pp(1) ;p2=pp(3) ;p3=pp(2)
       m1=mm(2) ;m2=mm(1) ;m3=mm(3)
     endif
   elseif (ap(2).le.ap(3).and.ap(2).le.ap(1)) then
     if (ap(3).le.ap(1)) then
       a1=ap(2) ;a2=ap(3) ;a3=ap(1)
       p1=pp(2) ;p2=pp(3) ;p3=pp(1)
       m1=mm(2) ;m2=mm(3) ;m3=mm(1)
     else
       a1=ap(2) ;a2=ap(1) ;a3=ap(3)
       p1=pp(2) ;p2=pp(1) ;p3=pp(3)
       m1=mm(3) ;m2=mm(2) ;m3=mm(1)
     endif
   else
     if (ap(1).le.ap(2)) then
       a1=ap(3) ;a2=ap(1) ;a3=ap(2)
       p1=pp(3) ;p2=pp(1) ;p3=pp(2)
       m1=mm(3) ;m2=mm(1) ;m3=mm(2)
     else
       a1=ap(3) ;a2=ap(2) ;a3=ap(1)
       p1=pp(3) ;p2=pp(2) ;p3=pp(1)
       m1=mm(1) ;m2=mm(3) ;m3=mm(2)
     endif
   endif
!
! Need to cut out negligible squared momenta
   thrs = smax*neglig(prcpar)
!
! Add infinitesimal imaginary parts to masses
   m1 = m1 - abs(areal(m1))*IEPS
   m2 = m2 - abs(areal(m2))*IEPS
   m3 = m3 - abs(areal(m3))*IEPS
!       
   if (a1.gt.thrs) then ! 3 non-zero squared momenta
     if (present(lam)) then ;slam=lam
                       else ;slam=kallen(p1,p2,p3)
     endif
     if (slam.eq.CZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3HV: ' &
         ,'threshold singularity, returning 0'
       rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
       return
     endif
     slam = mysqrt( slam ,1 )
     sm1=mysqrt(m1,-1) ;sm2=mysqrt(m2,-1) ;sm3=mysqrt(m3,-1)
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     rslt(0) = s3fun( p1,sm1,sm2 , (m2-m3)+p2    ,p3-p1-p2 ,p2 ,slam ) &
             - s3fun( p3,sm1,sm3 ,-(m1-m2)+p3-p2 ,p2-p1-p3 ,p1 ,slam ) &
             + s3fun( p2,sm2,sm3 ,-(m1-m2)+p3-p2 ,p1+p2-p3 ,p1 ,slam )
     rslt(0) = -rslt(0)/slam
!
   elseif (a2.gt.thrs) then ! 2 non-zero squared momenta
     if (p2.eq.p3) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop trif3HV: ' &
         ,'threshold singularity, returning 0'
       rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
       return
     endif
     sm1=mysqrt(m1,-1) ;sm2=mysqrt(m2,-1) ;sm3=mysqrt(m3,-1)
     yy = ( (m1-m2)-p3+p2 )/( p2-p3 )
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     rslt(0) = s3fun( p3,sm1,sm3 ,yy ) - s3fun( p2,sm2,sm3 ,yy )
     rslt(0) = rslt(0)/(p2-p3)
!
   elseif (a3.gt.thrs) then ! 1 non-zero squared momentum
     sm1=mysqrt(m1,-1) ;sm3=mysqrt(m3,-1)
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     yy = -( (m1-m2)-p3 )/p3
     rslt(0) = s3fun( p3,sm1,sm3 ,yy ) - s2fun( m2-m3 ,m3 ,yy )
     rslt(0) = -rslt(0)/p3
!
   else ! all squared momenta zero
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     a12=abs(m1-m2) ;a23=abs(m2-m3) ;a31=abs(m3-m1)
     if     (a12.ge.a23.and.a12.ge.a31) then
       if (a12.eq.RZRO) then ;rslt(0)=-1/(2*m3) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm3/qm1) - logc2(qm3/qm2) )/(m1-m2)
       endif
     elseif (a23.ge.a12.and.a23.ge.a31) then
       if (a23.eq.RZRO) then ;rslt(0)=-1/(2*m1) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm1/qm2) - logc2(qm1/qm3) )/(m2-m3)
       endif
     else
       if (a31.eq.RZRO) then ;rslt(0)=-1/(2*m2) ;else
       qm1=qonv(m1) ;qm2=qonv(m2) ;qm3=qonv(m3)
       rslt(0) = ( logc2(qm2/qm3) - logc2(qm2/qm1) )/(m3-m1)
       endif
     endif
   endif
!
   contains
!
     function s3fun( aa,s1,s2 ,t1,t2,t3,t4 ) result(rslt)
!***************************************************************
! int( ( ln(a*y^2+b*y+c) - ln(a*y0^2+b*y0+c) )/(y-y0) ,y=0..1 )
! with  b=s1^2-s2^2-aa  and  c=s2^2
! and with  y0  in terms of t1,t2,t3,t4 defined at the "present"
! function below.
! t4  should be  sqrt(lambda(aa,t2,t3))
!***************************************************************
  complex(kindr2) &   
       ,intent(in) :: aa,s1,s2,t1
  complex(kindr2) &   
       ,optional,intent(in) :: t2,t3
  complex(kindr2) &   
       ,optional,intent(inout) :: t4
  complex(kindr2) &   
       :: rslt ,cc,bb,dd,y0,y1,y2,zz,hh,alpha
  real(kindr2) &  
       :: rez,arez,aimz
     type(qmplx_type) :: q1,q2
!
     bb = (s1+s2)*(s1-s2)-aa
     cc = s2*s2
     dd = (aa-(s1+s2)**2)*(aa-(s1-s2)**2)
     dd = sqrt( dd )!+ sign(abs(dd),areal(aa))*IEPS )
     call solabc( y1,y2 ,dd ,aa,bb,cc ,1 )
!
     if (present(t4)) then
       call solabc( alpha,hh ,t4 ,aa,t2,t3 ,1 )
       y0 = -(t1+bb*alpha)/t4
     else
       y0 = t1
     endif
!
     q1 = qonv(y0-y1)
     q2 = qonv(y0-y2)
     rslt = li2c(qonv(-y1)/q1) - li2c(qonv(1-y1)/q1) &
          + li2c(qonv(-y2)/q2) - li2c(qonv(1-y2)/q2)
! Take some care about the imaginary part of  a*y0^2+b*y0+c=a*(y0-y1)*(y0-y2)
     zz = y0*(aa*y0+bb)
     rez=areal(zz)  ;arez=abs(rez) ;aimz=abs(aimag(zz))
     if (arez*EPSN*EPSN.le.aimz*neglig(prcpar).and.aimz.le.arez*neglig(prcpar)) then
! Here, the value of Imz is just numerical noise due to cancellations.
! Realize that |Imz|~eps^2 indicates there were no such cancellations,
! so the lower limit is needed in in the if-statement!
       zz = (rez + cc)/aa
     else
       zz = (zz + cc)/aa
     endif
     hh = eta3(-y1,-y2,cc/aa) - eta3(y0-y1,y0-y2,zz)
     if (areal(aa).lt.RZRO.and.aimag(zz).lt.RZRO) hh = hh - 2*IPI
     if (hh.ne.CZRO) rslt = rslt + hh*logc(qonv((y0-1)/y0,1))
!
     end function
!
     function s2fun( aa,bb ,y0 ) result(rslt)
!**************************************************
! int( ( ln(a*y+b) - ln(a*y0+b) )/(y-y0) ,y=0..1 )
!**************************************************
  complex(kindr2) &   
       ,intent(in) :: aa,bb,y0
  complex(kindr2) &   
       :: rslt ,y1,hh
     type(qmplx_type) :: q1
     y1 = -bb/aa
     q1 = qonv(y0-y1)
     rslt = li2c(qonv(-y1,-1)/q1) - li2c(qonv(1-y1,-1)/q1)
! aa may have imaginary part, so  theta(-aa)*theta(-Im(y0-y1))  is not
! sufficient and need the following:
     hh = eta5( aa ,-y1,bb ,y0-y1,aa*(y0-y1) )
     if (hh.ne.CZRO) rslt = rslt + hh*logc(qonv((y0-1)/y0,1))
     end function
!
   end subroutine


end module


module avh_olo_qp_box
  use avh_olo_units
  use avh_olo_qp_prec
  use avh_olo_qp_auxfun
  use avh_olo_qp_qmplx
  implicit none
  private
  public :: box00,box03,box05,box06,box07,box08,box09,box10,box11,box12 &
           ,box13,box14,box15,box16,boxf1,boxf2,boxf3,boxf5,boxf4 &
           ,permtable,casetable,base
  integer ,parameter ::  permtable(6,0:15)=reshape((/ &
     1,2,3,4 ,5,6 &! 0, 0 masses non-zero,           no perm
    ,1,2,3,4 ,5,6 &! 1, 1 mass non-zero,             no perm
    ,4,1,2,3 ,6,5 &! 2, 1 mass non-zero,             1 cyclic perm
    ,1,2,3,4 ,5,6 &! 3, 2 neighbour masses non-zero, no perm
    ,3,4,1,2 ,5,6 &! 4, 1 mass   non-zero,           2 cyclic perm's
    ,1,2,3,4 ,5,6 &! 5, 2 opposite masses non-zero,  no perm
    ,4,1,2,3 ,6,5 &! 6, 2 neighbour masses non-zero, 1 cyclic perm
    ,1,2,3,4 ,5,6 &! 7, 3 masses non-zero,           no perm
    ,2,3,4,1 ,6,5 &! 8, 1 mass   non-zero,           3 cyclic perm's
    ,2,3,4,1 ,6,5 &! 9, 2 neighbour masses non-zero, 3 cyclic perm's
    ,4,1,2,3 ,6,5 &!10, 2 opposite masses non-zero,  1 cyclic perm
    ,2,3,4,1 ,6,5 &!11, 3 masses non-zero,           3 cyclic perm's
    ,3,4,1,2 ,5,6 &!12, 2 neighbour masses non-zero, 2 cyclic perm's
    ,3,4,1,2 ,5,6 &!13, 3 masses non-zero,           2 cyclic perm's
    ,4,1,2,3 ,6,5 &!14, 3 masses non-zero,           1 cyclic perm
    ,1,2,3,4 ,5,6 &!15, 4 masses non-zero,           no perm
    /),(/6,16/)) !          0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
  integer ,parameter :: casetable(0:15)= &
                          (/0,1,1,2,1,5,2,3,1,2, 5, 3, 2, 3, 3, 4/)
  integer ,parameter :: base(4)=(/8,4,2,1/)
contains

   subroutine box16( rslt ,p2,p3,p12,p23 ,m2,m3,m4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                     d^(Dim)q
! ------ | ------------------------------------------------------
! i*pi^2 / q^2 [(q+k1)^2-m2] [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=m2, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=m4
! m2,m4 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p12,p23 ,m2,m3,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp12,cp23,cm2,cm3,cm4,sm1,sm2,sm3,sm4 &
                     ,r13,r23,r24,r34,d23,d24,d34,log24,cc
   type(qmplx_type) :: q13,q23,q24,q34,qss,qy1,qy2,qz1,qz2
!
   if (abs(m2).gt.abs(m4)) then
     cm2=m2 ;cm4=m4 ;cp2=p2 ;cp3=p3
   else
     cm2=m4 ;cm4=m2 ;cp2=p3 ;cp3=p2
   endif
   cm3=m3 ;cp12=p12 ;cp23=p23
!
   if (cp12.eq.cm3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box16: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm1 = abs(rmu)
   sm2 = mysqrt(cm2)
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
!
   r13 = (cm3-cp12)/(sm1*sm3)
   call rfun( r23,d23 ,(cm2+cm3-cp2 )/(sm2*sm3) )
   call rfun( r24,d24 ,(cm2+cm4-cp23)/(sm2*sm4) )
   call rfun( r34,d34 ,(cm3+cm4-cp3 )/(sm3*sm4) )
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
!
   if (r24.eq.-CONE) then 
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box16: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   qss = q23*q34
   qy1 = qss*q24
   qy2 = qss/q24
!
   qss = q23/q34
   qz1 = qss*q24
   qz2 = qss/q24
!
   qss = q13*q23
   qss = (qss*qss)/q24
!
   cc = 1/( sm2*sm4*(cp12-cm3) )
   log24 = logc2(q24)*r24/(1+r24)
   rslt(2) = 0
   rslt(1) = -log24
   rslt(0) = log24*logc(qss) + li2c2(q24*q24,qonv(1))*r24 &
           - li2c2(qy1,qy2)*r23*r34 - li2c2(qz1,qz2)*r23/r34
   rslt(1) = cc*rslt(1)
   rslt(0) = cc*rslt(0)
   end subroutine


   subroutine box15( rslt ,p2,p3,p12,p23 ,m2,m4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 [(q+k1)^2-m2] (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=m2, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=m4
! m2,m4 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p12,p23 ,m2,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp12,cp23,cm2,cm4,sm1,sm2,sm3,sm4 &
                     ,r13,r23,r24,r34,d24,log24,cc
   type(qmplx_type) :: q13,q23,q24,q34,qss,qz1,qz2
!
   if (abs(m2-p2).gt.abs(m4-p3)) then
     cm2=m2 ;cm4=m4 ;cp2=p2 ;cp3=p3
   else
     cm2=m4 ;cm4=m2 ;cp2=p3 ;cp3=p2
   endif
   cp12=p12 ;cp23=p23
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box15: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm1 = abs(rmu)
   sm2 = mysqrt(cm2)
   sm4 = mysqrt(cm4)
   sm3 = abs(sm2)
   r13 = (       -cp12)/(sm1*sm3)
   r23 = (cm2    -cp2 )/(sm2*sm3)
   r34 = (    cm4-cp3 )/(sm3*sm4)
   call rfun( r24,d24 ,(cm2+cm4-cp23)/(sm2*sm4) )
!
   if (r24.eq.-CONE) then 
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box15: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   q13 = qonv(r13,-1)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1)
!
   qss = q13/q23
   qss = (qss*qss)/q24
!
   cc = r24/(sm2*sm4*cp12)
   log24 = logc2(q24)/(1+r24)
   rslt(2) = 0
   rslt(1) = -log24
   rslt(0) = log24 * logc(qss) + li2c2(q24*q24,qonv(1))
   if (r34.ne.CZRO) then
     qss = q34/q23
     qz1 = qss*q24
     qz2 = qss/q24
     rslt(0) = rslt(0) - li2c2(qz1,qz2)*r34/(r23*r24)
   endif
   rslt(1) = cc*rslt(1)
   rslt(0) = cc*rslt(0)
   end subroutine


   subroutine box14( rslt ,cp12,cp23 ,cm2,cm4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 [(q+k1)^2-m2] (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=m2, k2^2=m2, k3^2=m4, (k1+k2+k3)^2=m4
! m2,m4 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp12,cp23,cm2,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: sm2,sm4,r24,d24,cc
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box14: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm2 = mysqrt(cm2)
   sm4 = mysqrt(cm4)
   call rfun( r24,d24 ,(cm2+cm4-cp23)/(sm2*sm4) )
!
   if (r24.eq.-CONE) then 
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box14: ' &
       ,'threshold singularity, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   cc = -2*logc2(qonv(r24,-1))*r24/(1+r24)/(sm2*sm4*cp12)
!
   rslt(2) = 0
   rslt(1) = cc
   rslt(0) = -cc*logc(qonv(-cp12/(rmu*rmu),-1))
   end subroutine


   subroutine box13( rslt ,p2,p3,p4,p12,p23 ,m3,m4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=p4
! m3,m4 should NOT be identically 0d0
! p4 should NOT be identical to m4
! p2 should NOT be identical to m3
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p4,p12,p23,m3,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp4,cp12,cp23,cm3,cm4,sm3,sm4,sm1,sm2 &
             ,r13,r14,r23,r24,r34,d34,cc,logd,li2d,loge,li2f,li2b,li2e
   type(qmplx_type) :: q13,q14,q23,q24,q34,qy1,qy2
  real(kindr2) &  
     :: h1,h2
!
   if (p12.eq.m3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box13: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (p23.eq.m4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box13: ' &
       ,'p23=m4, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   h1 = abs((m3-p12)*(m4-p23))
   h2 = abs((m3-p2 )*(m4-p4 ))
   if (h1.ge.h2) then
     cp2=p2  ;cp3=p3 ;cp4=p4  ;cp12=p12 ;cp23=p23 ;cm3=m3 ;cm4=m4
   else
     cp2=p12 ;cp3=p3 ;cp4=p23 ;cp12=p2  ;cp23=p4  ;cm3=m3 ;cm4=m4
   endif
!
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
   sm1 = abs(rmu)
   sm2 = sm1
!
   r13 = (cm3-cp12)/(sm1*sm3)
   r14 = (cm4-cp4 )/(sm1*sm4)
   r23 = (cm3-cp2 )/(sm2*sm3)
   r24 = (cm4-cp23)/(sm2*sm4)
   call rfun( r34,d34 ,(cm3+cm4-cp3)/(sm3*sm4) )
!
   q13 = qonv(r13,-1)
   q14 = qonv(r14,-1)
   q23 = qonv(r23,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1) 
!
   qy1 = q14*q23/q13/q24
   logd = logc2(qy1     )/(r13*r24)
   li2d = li2c2(qy1,qonv(1))/(r13*r24)
   loge = logc(q13)
!
   qy1 = q23/q24
   qy2 = q13/q14
   li2f = li2c2( qy1*q34,qy2*q34 )*r34/(r14*r24)
   li2b = li2c2( qy1/q34,qy2/q34 )/(r34*r14*r24)
   li2e = li2c2( q14/q24,q13/q23 )/(r23*r24)
!
   rslt(2) = 0
   rslt(1) = logd
   rslt(0) = li2f + li2b + 2*li2e - 2*li2d - 2*logd*loge
   cc = sm1*sm2*sm3*sm4
   rslt(1) = rslt(1)/cc
   rslt(0) = rslt(0)/cc
   end subroutine


   subroutine box12( rslt ,cp3,cp4,cp12,cp23 ,cm3,cm4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=m3, k3^2=p3, (k1+k2+k3)^2=p4
! m3,m4 should NOT be indentiallcy 0d0
! p4 should NOT be identical to m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cp4,cp12,cp23,cm3,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: sm3,sm4,sm1,sm2,r13,r14,r24,r34,d34,cc &
                     ,log13,log14,log24,log34,li2f,li2b,li2d
   type(qmplx_type) :: q13,q14,q24,q34,qyy
!
   if (cp12.eq.cm3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box12: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box12: ' &
       ,'p23=m4, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
   sm1 = abs(rmu)
   sm2 = sm1
!
   r13 = (cm3-cp12)/(sm1*sm3)
   r14 = (cm4-cp4 )/(sm1*sm4)
   r24 = (cm4-cp23)/(sm2*sm4)
   call rfun( r34,d34 ,(cm3+cm4-cp3)/(sm3*sm4) )
!
   q13 = qonv(r13,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   q34 = qonv(r34,-1) 
!
   log13 = logc(q13) 
   log14 = logc(q14) 
   log24 = logc(q24) 
   log34 = logc(q34) 
!
   qyy = q14/q13
   li2f = li2c(qyy*q34)
   li2b = li2c(qyy/q34)
   li2d = li2c(q14/q24)
!
   rslt(2) = 1
   rslt(2) = rslt(2)/2
   rslt(1) = log14 - log24 - log13
   rslt(0) = 2*log13*log24 - log14*log14 - log34*log34 &
           - 2*li2d - li2f - li2b - 3*PISQo24
   cc = (cm3-cp12)*(cm4-cp23) ! = sm1*sm2*sm3*sm4*r13*r24
   rslt(2) = rslt(2)/cc
   rslt(1) = rslt(1)/cc
   rslt(0) = rslt(0)/cc
   end subroutine


   subroutine box11( rslt ,cp3,cp12,cp23 ,cm3,cm4 ,rmu )
!*******************************************************************
! calculates
!
!    C   /                  d^(Dim)q
! ------ | -------------------------------------------------
! i*pi^2 / q^2 (q+k1)^2 [(q+k1+k2)^2-m3] [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=m3, k3^2=p3, (k1+k2+k3)^2=m4
! m3,m4 should NOT be indentiallcy 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cp12,cp23,cm3,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: sm3,sm4,sm1,sm2,r13,r24,r34,d34 &
                     ,cc,log13,log24,log34
!
   if (cp12.eq.cm3) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box11: ' &
       ,'p12=m3, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box11: ' &
       ,'p23=m4, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   sm3 = mysqrt(cm3)
   sm4 = mysqrt(cm4)
   sm1 = abs(rmu)
   sm2 = sm1
!
   r13 = (cm3-cp12)/(sm1*sm3)
   r24 = (cm4-cp23)/(sm2*sm4)
   call rfun( r34,d34 ,(cm3+cm4-cp3 )/(sm3*sm4) )
!
   log13 = logc(qonv(r13,-1)) 
   log24 = logc(qonv(r24,-1)) 
   log34 = logc(qonv(r34,-1)) 
!
   rslt(2) = 1
   rslt(1) = -log13-log24
   rslt(0) = 2*log13*log24 - log34*log34 - 14*PISQo24
   cc = (cm3-cp12)*(cm4-cp23) ! = sm1*sm2*sm3*sm4*r13*r24
   rslt(2) = rslt(2)/cc
   rslt(1) = rslt(1)/cc
   rslt(0) = rslt(0)/cc
   end subroutine


   subroutine box10( rslt ,p2,p3,p4,p12,p23 ,m4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=p4
! m4 should NOT be identically 0d0
! p2 should NOT be identically 0d0
! p4 should NOT be identical to m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p4,p12,p23,m4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: cp2,cp3,cp4,cp12,cp23,cm4,r13,r14,r23,r24,r34,z1,z0
   type(qmplx_type) :: q13,q14,q23,q24,q34,qm4,qxx,qx1,qx2
  real(kindr2) &  
     :: h1,h2
!
   if (p12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box10: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (p23.eq.m4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box10: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   h1 = abs(p12*(m4-p23))
   h2 = abs( p2*(m4-p4 ))
   if (h1.ge.h2) then
     cp2=p2  ;cp3=p3 ;cp4=p4  ;cp12=p12 ;cp23=p23 ;cm4=m4
   else
     cp2=p12 ;cp3=p3 ;cp4=p23 ;cp12=p2  ;cp23=p4  ;cm4=m4
   endif
!
   r23 =    -cp2
   r13 =    -cp12
   r34 = cm4-cp3
   r14 = cm4-cp4
   r24 = cm4-cp23
   q23 = qonv(r23,-1)
   q13 = qonv(r13,-1)
   q34 = qonv(r34,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   if (r34.ne.CZRO) then
     qx1 = q34/qm4
     qx2 = qx1*q14/q13
     qx1 = qx1*q24/q23
     z0 = -li2c2(qx1,qx2)*r34/(2*cm4*r23)
   else
     z0 = 0
   endif
!
   qx1 = q23/q13
   qx2 = q24/q14
   qxx = qx1/qx2
   z1 = -logc2(qxx)/r24
   z0 = z0 - li2c2(qx1,qx2)/r14
   z0 = z0 + li2c2(qxx,qonv(1))/r24
   z0 = z0 + z1*( logc(qm4/q24) - logc(qm4/(rmu*rmu))/2 )
!
   rslt(2) = 0
   rslt(1) = -z1/r13
   rslt(0) = -2*z0/r13
   end subroutine


   subroutine box09( rslt ,cp2,cp3,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=0, k2^2=p2, k3^2=p3, (k1+k2+k3)^2=m4
! m4 should NOT be identically 0d0
! p2 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp2,cp3,cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
  complex(kindr2) &   
     :: logm,log12,log23,li12,li23,z2,z1,z0,cc &
                     ,r13,r23,r24,r34
   type(qmplx_type) :: q13,q23,q24,q34,qm4,qxx
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box09: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box09: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   r23 =    -cp2
   r13 =    -cp12
   r34 = cm4-cp3
   r24 = cm4-cp23
   q23 = qonv(r23,-1)
   q13 = qonv(r13,-1)
   q34 = qonv(r34,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   logm  = logc(qm4/(rmu*rmu))
   qxx = q13/q23
   log12 = logc(qxx)
   li12  = li2c(qxx)
!
   qxx = q24/qm4
   log23 = logc(qxx)
   li23  = li2c(qxx*q34/q23)
!
   z2 = 1
   z2 = z2/2
   z1 = -log12 - log23
   z0 = li23 + 2*li12 + z1*z1 + PISQo24
   cc = 1/(r13*r24)
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine box08( rslt ,cp3,cp4,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=k2^2=0, k3^2=p3, (k1+k2+k3)^2=p4
! mm should NOT be identically 0d0
! p3 NOR p4 should be identically m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp3,cp4,cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q13,q14,q24,q34,qm4,qxx,qx1,qx2,qx3
  complex(kindr2) &   
     :: r13,r14,r24,r34,z1,z0,cc
  real(kindr2) &  
     :: rmu2
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box08: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box08: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   rmu2 = rmu*rmu
   r13 =    -cp12
   r34 = cm4-cp3
   r14 = cm4-cp4
   r24 = cm4-cp23
   q13 = qonv(r13,-1)
   q34 = qonv(r34,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   qx1 = q34/q24
   qx2 = q14/q24
   qx3 = q13/rmu2
   z1 = logc(qx1*qx2/qx3)
   z0 = 2*( logc(q24/rmu2)*logc(qx3) - (li2c(qx1)+li2c(qx2)) )
!
   qx1 = q34/rmu2
   qx2 = q14/rmu2
   qxx = qx1*qx2/qx3
   z0 = z0 - logc(qx1)**2 - logc(qx2)**2 &
           + logc(qxx)**2/2 + li2c(qm4/qxx/rmu2)
!
   cc = 1/(r13*r24)
   rslt(2) = cc
   rslt(1) = cc*z1
   rslt(0) = cc*( z0 - 6*PISQo24 )
   end subroutine


   subroutine box07( rslt ,cp4,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=k2^2=0, k3^2=m4, (k1+k2+k3)^2=p4
! m3 should NOT be identically 0d0
! p4 should NOT be identically m4
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp4,cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q13,q14,q24,qm4
  complex(kindr2) &   
     :: r13,r14,r24,logm,log12,log23,log4,li423 &
                     ,z2,z1,z0,cc
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box07: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box07: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   r13 =    -cp12
   r14 = cm4-cp4
   r24 = cm4-cp23
   q13 = qonv(r13,-1)
   q14 = qonv(r14,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   logm  = logc(qm4/(rmu*rmu))
   log12 = logc(q13/qm4)
   log23 = logc(q24/qm4)
   log4  = logc(q14/qm4)
   li423 = li2c(q14/q24)
!
   z2 = 3
   z2 = z2/2
   z1 = -2*log23 - log12 + log4
   z0 = 2*(log12*log23 - li423) - log4*log4 - 13*PISQo24
   cc = 1/(r13*r24)
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine box06( rslt ,cp12,cp23 ,cm4 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | --------------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 [(q+k1+k2+k3)^2-m4]
!
! with  k1^2=k2^2=0, k3^2=(k1+k2+k3)^2=m4
! m3 should NOT be identically 0d0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp12,cp23,cm4
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q13,q24,qm4
  complex(kindr2) &   
     :: r13,r24,logm,log1,log2,z2,z1,z0,cc
!
   if (cp12.eq.CZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box06: ' &
       ,'p12=0, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
   if (cp23.eq.cm4) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop box06: ' &
       ,'p23=mm, returning 0'
     rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
     return
   endif
!
   r13 =    -cp12
   r24 = cm4-cp23
   q13 = qonv(r13,-1)
   q24 = qonv(r24,-1)
   qm4 = qonv(cm4,-1)
!
   logm = logc(qm4/(rmu*rmu))
   log1 = logc(q13/qm4)
   log2 = logc(q24/qm4)
!
   z2 = 2
   z1 = -2*log2 - log1
   z0 = 2*(log2*log1 - 8*PISQo24)
   cc = 1/(r13*r24)
   rslt(2) = cc*z2
   rslt(1) = cc*(z1 - z2*logm)
   rslt(0) = cc*(z0 + (z2*logm/2-z1)*logm)
   end subroutine


   subroutine box03( rslt ,p2,p4,p5,p6 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | ---------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 (q+k1+k2+k3)^2
!
! with  k1^2=k3^2=0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p4,p5,p6 
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) :: q2,q4,q5,q6,q26,q54,qy
  complex(kindr2) &   
     :: logy
  real(kindr2) &  
     :: rmu2
!
   rmu2 = rmu*rmu
   q2 = qonv(-p2,-1)
   q4 = qonv(-p4,-1)
   q5 = qonv(-p5,-1)
   q6 = qonv(-p6,-1)
   q26 = q2/q6
   q54 = q5/q4
   qy = q26/q54
   logy = logc2(qy)/(p5*p6)
   rslt(1) = logy
   rslt(0) = li2c2(q6/q4,q2/q5)/(p4*p5) &
           + li2c2(q54,q26)/(p4*p6)     &
           - li2c2(qonv(1),qy)/(p5*p6) &
           - logy*logc(q54*q2*q6/(rmu2*rmu2))/2
   rslt(2) = 0
   rslt(1) = 2*rslt(1)
   rslt(0) = 2*rslt(0)
   end subroutine


   subroutine box05( rslt ,p2,p3,p4,p5,p6 ,rmu )
!*******************************************************************
! calculates
!
!     C   /               d^(Dim)q
!  ------ | ---------------------------------------
!  i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 (q+k1+k2+k3)^2
!
! with  k1^2=0
!*******************************************************************
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: p2,p3,p4,p5,p6
  real(kindr2) &  
     ,intent(in)  :: rmu
   type(qmplx_type) ::q2,q3,q4,q5,q6 ,q25,q64,qy,qz
  complex(kindr2) &   
     :: logy
  real(kindr2) &  
     :: rmu2
!
   rmu2 = rmu*rmu
   q2 = qonv(-p2,-1)
   q3 = qonv(-p3,-1)
   q4 = qonv(-p4,-1)
   q5 = qonv(-p5,-1)
   q6 = qonv(-p6,-1)
   q25 = q2/q5
   q64 = q6/q4
   qy = q25/q64
   qz = q64*q2*q5*q6*q6/q3/q3/(rmu2*rmu2)
!
   logy = logc2(qy)/(p5*p6)
   rslt(2) = 0
   rslt(1) = logy
   rslt(0) = li2c2(q64,q25)/(p4*p5) &
           - li2c2(qonv(1),qy)/(p5*p6) &
           - logy*logc(qz)/4
   rslt(0) = 2*rslt(0)
   end subroutine


   subroutine box00( rslt ,cp ,api ,rmu )
!*******************************************************************
! calculates
!               C   /              d^(Dim)q
!            ------ | ---------------------------------------
!            i*pi^2 / q^2 (q+k1)^2 (q+k1+k2)^2 (q+k1+k2+k3)^2
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  p1 = k1^2,  p2 = k2^2,  p3 = k3^2,  p4 = (k1+k2+k3)^2,
!         p12 = (k1+k2)^2,  p23 = (k2+k3)^2
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! If any of these numbers is IDENTICALLY 0d0, the corresponding
! IR-singular case is returned.
!*******************************************************************
   use avh_olo_qp_olog
   use avh_olo_qp_dilog
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: cp(6)
  real(kindr2) &  
     ,intent(in)  :: api(6),rmu
  complex(kindr2) &   
     :: log3,log4,log5,log6,li24,li25,li26 &
                     ,li254,li263
  real(kindr2) &  
     :: rp1,rp2,rp3,rp4,rp5,rp6,pp(6),ap(6),gg,ff,hh,arg,rmu2
   integer :: icase,sf,sgn,i3,i4,i5,i6
   integer ,parameter :: base(4)=(/8,4,2,1/)
!
   rmu2 = rmu*rmu
   ff = api(5)*api(6)
   gg = api(2)*api(4)
   hh = api(1)*api(3)
   if     (ff.ge.gg.and.ff.ge.hh) then
     pp(1)=areal(cp(1)) ;ap(1)=api(1)
     pp(2)=areal(cp(2)) ;ap(2)=api(2)
     pp(3)=areal(cp(3)) ;ap(3)=api(3)
     pp(4)=areal(cp(4)) ;ap(4)=api(4)
     pp(5)=areal(cp(5)) ;ap(5)=api(5)
     pp(6)=areal(cp(6)) ;ap(6)=api(6)
   elseif (gg.ge.ff.and.gg.ge.hh) then
     pp(1)=areal(cp(1)) ;ap(1)=api(1)
     pp(2)=areal(cp(6)) ;ap(2)=api(6)
     pp(3)=areal(cp(3)) ;ap(3)=api(3)
     pp(4)=areal(cp(5)) ;ap(4)=api(5)
     pp(5)=areal(cp(4)) ;ap(5)=api(4)
     pp(6)=areal(cp(2)) ;ap(6)=api(2)
   else
     pp(1)=areal(cp(5)) ;ap(1)=api(5)
     pp(2)=areal(cp(2)) ;ap(2)=api(2)
     pp(3)=areal(cp(6)) ;ap(3)=api(6)
     pp(4)=areal(cp(4)) ;ap(4)=api(4)
     pp(5)=areal(cp(1)) ;ap(5)=api(1)
     pp(6)=areal(cp(3)) ;ap(6)=api(3)
   endif
!
   icase = 0
   if (ap(1).gt.RZRO) icase = icase + base(1)
   if (ap(2).gt.RZRO) icase = icase + base(2)
   if (ap(3).gt.RZRO) icase = icase + base(3)
   if (ap(4).gt.RZRO) icase = icase + base(4)
   rp1 = pp(permtable(1,icase))
   rp2 = pp(permtable(2,icase))
   rp3 = pp(permtable(3,icase))
   rp4 = pp(permtable(4,icase))
   rp5 = pp(permtable(5,icase))
   rp6 = pp(permtable(6,icase))
   icase = casetable(   icase)
!
   i3=0 ;if (-rp3.lt.RZRO) i3=-1
   i4=0 ;if (-rp4.lt.RZRO) i4=-1
   i5=0 ;if (-rp5.lt.RZRO) i5=-1
   i6=0 ;if (-rp6.lt.RZRO) i6=-1
!
   if     (icase.eq.0) then
! 0 masses non-zero
     gg = 1/( rp5 * rp6 )
     log5 = olog(abs(rp5/rmu2),i5)
     log6 = olog(abs(rp6/rmu2),i6)
     rslt(2) = gg*( 4 )
     rslt(1) = gg*(-2*(log5 + log6) )
     rslt(0) = gg*( log5**2 + log6**2 - olog(abs(rp5/rp6),i5-i6)**2 - 32*PISQo24 )
   elseif (icase.eq.1) then
! 1 mass non-zero
     gg = 1/( rp5 * rp6 )
     ff =  gg*( rp5 + rp6 - rp4 )
     log4 = olog(abs(rp4/rmu2),i4)
     log5 = olog(abs(rp5/rmu2),i5)
     log6 = olog(abs(rp6/rmu2),i6)
     sf = sgnRe(ff)
     sgn = 0
       arg = rp4*ff 
       if (arg.lt.RZRO) sgn = sf
       li24 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp5*ff 
       if (arg.lt.RZRO) sgn = sf
       li25 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp6*ff 
       if (arg.lt.RZRO) sgn = sf
       li26 = dilog(abs(arg),sgn)
     rslt(2) = gg*( 2 )
     rslt(1) = gg*( 2*(log4-log5-log6) )
     rslt(0) = gg*( log5**2 + log6**2 - log4**2 - 12*PISQo24 &
                   + 2*(li25 + li26 - li24) )
   elseif (icase.eq.2) then
! 2 neighbour masses non-zero
     gg = 1/( rp5 * rp6 )
     ff =  gg*( rp5 + rp6 - rp4 )
     log3 = olog(abs(rp3/rmu2),i3)
     log4 = olog(abs(rp4/rmu2),i4)
     log5 = olog(abs(rp5/rmu2),i5)
     log6 = olog(abs(rp6/rmu2),i6)
     li254 = dilog( abs(rp4/rp5) ,i4-i5 )
     li263 = dilog( abs(rp3/rp6) ,i3-i6 )
     sf = sgnRe(ff)
     sgn = 0
       arg = rp4*ff 
       if (arg.lt.RZRO) sgn = sf
       li24 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp5*ff 
       if (arg.lt.RZRO) sgn = sf
       li25 = dilog(abs(arg),sgn)
     sgn = 0
       arg = rp6*ff 
       if (arg.lt.RZRO) sgn = sf
       li26 = dilog(abs(arg),sgn)
     rslt(2) = gg
     rslt(1) = gg*( log4 + log3 - log5 - 2*log6 )
     rslt(0) = gg*( log5**2 + log6**2 - log3**2 - log4**2 &
                   + (log3 + log4 - log5)**2/2 &
                   - 2*PISQo24 + 2*(li254 - li263 + li25 + li26 - li24) )
   elseif (icase.eq.5) then
! 2 opposite masses non-zero
     call box03( rslt ,acmplx(rp2),acmplx(rp4) &
                      ,acmplx(rp5),acmplx(rp6) ,rmu )
   elseif (icase.eq.3) then
! 3 masses non-zero
     call box05( rslt ,acmplx(rp2),acmplx(rp3) &
                      ,acmplx(rp4),acmplx(rp5) &
                      ,acmplx(rp6) ,rmu )
   elseif (icase.eq.4) then
! 4 masses non-zero
     call boxf0( rslt ,acmplx(rp1),acmplx(rp2) &
                      ,acmplx(rp3),acmplx(rp4) &
                      ,acmplx(rp5),acmplx(rp6) )
   endif
   end subroutine

  
  subroutine boxf0( rslt ,p1,p2,p3,p4,p12,p23 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with all internal masses
! equal zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23
  type(qmplx_type) :: q12,q13,q14,q23,q24,q34,qx1,qx2,qss
  complex(kindr2) &   
    :: aa,bb,cc,dd,x1,x2,ss,r12,r13,r14,r23,r24,r34
  real(kindr2) &  
    :: hh
!
  r12 = -p1  !  p1
  r13 = -p12 !  p1+p2
  r14 = -p4  !  p1+p2+p3
  r23 = -p2  !  p2
  r24 = -p23 !  p2+p3
  r34 = -p3  !  p3      
!
  aa = r34*r24
!
  if (r13.eq.CZRO.or.aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf0: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = r13*r24 + r12*r34 - r14*r23
  cc = r12*r13
  hh = areal(r23)
  dd = mysqrt( bb*bb - 4*aa*cc , -areal(aa)*hh )
  call solabc(x1,x2,dd ,aa,bb,cc ,1)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 , hh)
  qx2 = qonv(x2 ,-hh)
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  qss = q34/q13
  rslt(0) = rslt(0) + li2c2(qx1*qss,qx2*qss) * r34/r13
!
  qss = q24/q12
  rslt(0) = rslt(0) + li2c2(qx1*qss,qx2*qss) * r24/r12
!
  ss = -logc2(qx1/qx2) / x2
  rslt(0) = rslt(0) + ss*( logc(qx1*qx2)/2 - logc(q12*q13/q14/q23) )
!
  rslt(0) = -rslt(0) / aa
  end subroutine


  subroutine boxf1( rslt ,p1,p2,p3,p4,p12,p23 ,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with one internal mass
! non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23 ,m4
  type(qmplx_type) :: qx1,qx2,qss,q12,q13,q14,q23,q24,q34
  complex(kindr2) &   
    :: smm,sm4,aa,bb,cc,dd,x1,x2,r12,r13,r14,r23,r24,r34
  logical :: r12zero,r13zero,r14zero
!
  sm4 = mysqrt(m4)
  smm = abs(sm4) 
!
  r12 = ( m4-p4 -p4 *IEPS )/(smm*sm4)
  r13 = ( m4-p23-p23*IEPS )/(smm*sm4)
  r14 = ( m4-p3 -p3 *IEPS )/(smm*sm4)
  r23 = (   -p1 -p1 *IEPS )/(smm*smm)
  r24 = (   -p12-p12*IEPS )/(smm*smm)
  r34 = (   -p2 -p2 *IEPS )/(smm*smm)
!
  r12zero=(abs(areal(r12))+abs(aimag(r12)).lt.neglig(prcpar))
  r13zero=(abs(areal(r13))+abs(aimag(r13)).lt.neglig(prcpar))
  r14zero=(abs(areal(r14))+abs(aimag(r14)).lt.neglig(prcpar))
!
  aa = r34*r24
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf1: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = r13*r24 + r12*r34 - r14*r23
  cc = r12*r13 - r23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 ,1 )
  qx2 = qonv(x2 ,1 )
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  if (r12zero.and.r13zero) then
    qss = qx1*qx2*q34*q24/q23
    qss = qss*qss
    rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
  else
    if (r13zero) then
      qss = q34*q12/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    else
      qss = q34/q13
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34/r13
    endif
    if (r12zero) then
      qss = q24*q13/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    else
      qss = q24/q12
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r24/r12
    endif
    if (.not.r12zero.and..not.r13zero) then
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( q12*q13/q23 )/x2
    endif
  endif
!
  if (.not.r14zero) then
    rslt(0) = rslt(0) - li2c2( qx1*q14 ,qx2*q14 )*r14
  endif
!
  rslt(0) = -rslt(0)/(aa*smm*smm*smm*sm4)
  end subroutine


  subroutine boxf5( rslt ,p1,p2,p3,p4,p12,p23, m2,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with two opposite internal
! masses non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m2,m4
  call boxf2( rslt ,p12,p2,p23,p4,p1,p3 ,m2,m4 )
  end subroutine


  subroutine boxf2( rslt ,p1,p2,p3,p4,p12,p23 ,m3,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with two adjacent internal
! masses non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m3,m4
  type(qmplx_type) :: qx1,qx2,qss,q12,q13,q14,q23,q24,q34
  complex(kindr2) &   
    :: smm,sm3,sm4,aa,bb,cc,dd,x1,x2 &
                    ,r12,r13,r14,r23,r24,r34,d14,k14
  logical :: r12zero,r13zero,r24zero,r34zero
!
  sm3 = mysqrt(m3)
  sm4 = mysqrt(m4)
!
  smm = abs(sm3)
!
  r12 = (    m4-p4 -p4 *IEPS )/(smm*sm4)
  r13 = (    m4-p23-p23*IEPS )/(smm*sm4)
  k14 = ( m3+m4-p3 -p3 *IEPS )/(sm3*sm4)
  r23 = (      -p1 -p1 *IEPS )/(smm*smm)
  r24 = (    m3-p12-p12*IEPS )/(smm*sm3)
  r34 = (    m3-p2 -p2 *IEPS )/(smm*sm3)
!
  r12zero = (abs(areal(r12))+abs(aimag(r12)).lt.neglig(prcpar))
  r13zero = (abs(areal(r13))+abs(aimag(r13)).lt.neglig(prcpar))
  r24zero = (abs(areal(r24))+abs(aimag(r24)).lt.neglig(prcpar))
  r34zero = (abs(areal(r34))+abs(aimag(r34)).lt.neglig(prcpar))
!
  if (r12zero.and.r24zero) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf2: ' &
       ,'m4=p4 and m3=p12, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
  if (r13zero.and.r34zero) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf2: ' &
       ,'m4=p23 and m3=p2, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  call rfun( r14,d14 ,k14 )
!
  aa = r34*r24 - r23
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf2: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = r13*r24 + r12*r34 - k14*r23
  cc = r12*r13 - r23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 ,1 )
  qx2 = qonv(x2 ,1 )
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  rslt(0) = rslt(0) - li2c2( qx1*q14 ,qx2*q14 )*r14
  rslt(0) = rslt(0) - li2c2( qx1/q14 ,qx2/q14 )/r14
!
  if (r12zero.and.r13zero) then
    qss = qx1*qx2*q34*q24/q23
    qss = qss*qss
    rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
  else
    if (r13zero) then
      qss = q34*q12/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    elseif (.not.r34zero) then
      qss = q34/q13
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34/r13
    endif
    if (r12zero) then
      qss = q24*q13/q23
      qss = qx1*qx2*qss*qss
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( qss )/(x2*2)
    elseif (.not.r24zero) then
      qss = q24/q12
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r24/r12
    endif
    if (.not.r12zero.and..not.r13zero) then
      rslt(0) = rslt(0) + logc2( qx1/qx2 )*logc( q12*q13/q23 )/x2 
    endif
  endif
!
  rslt(0) = -rslt(0)/(aa*smm*smm*sm3*sm4)
  end subroutine


  subroutine boxf3( rslt ,pp ,mm )
!*******************************************************************
! Finite 1-loop scalar 4-point function with three internal masses
! non-zero.
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: pp(6),mm(4)
  integer :: j
  integer ,parameter :: ip(6)=(/4,5,2,6,3,1/)
  integer ,parameter :: im(4)=(/4,1,3,2/)
  integer ,parameter :: ic(4,6)=reshape((/1,2,3,4 ,2,3,4,1 ,3,4,1,2 &
                                  ,4,1,2,3 ,5,6,5,6 ,6,5,6,5/),(/4,6/))
!
  if     (mm(1).eq.CZRO) then ;j=3
  elseif (mm(2).eq.CZRO) then ;j=4
  elseif (mm(3).eq.CZRO) then ;j=1
  else                        ;j=2
  endif
  call boxf33( rslt ,pp(ic(j,ip(1))) ,pp(ic(j,ip(2))) ,pp(ic(j,ip(3))) &
                    ,pp(ic(j,ip(4))) ,pp(ic(j,ip(5))) ,pp(ic(j,ip(6))) &
                    ,mm(ic(j,im(1))) ,mm(ic(j,im(2))) ,mm(ic(j,im(4))) )
  end subroutine

  subroutine boxf33( rslt ,p1,p2,p3,p4,p12,p23, m1,m2,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with three internal masses
! non-zero, and m3=0. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m1,m2,m4
  type(qmplx_type) :: qx1,qx2,qss,q12,q13,q14,q23,q24,q34,qy1,qy2
  complex(kindr2) &   
    :: sm1,sm2,sm3,sm4 ,aa,bb,cc,dd,x1,x2 &
                    ,r12,r13,r14,r23,r24,r34,d12,d14,d24,k12,k14,k24
  logical ::r13zero,r23zero,r34zero
!
  sm1 = mysqrt(m1)
  sm2 = mysqrt(m2)
  sm4 = mysqrt(m4)
  sm3 = abs(sm2)
!
  k12 = ( m1+m2-p1 -p1 *IEPS )/(sm1*sm2) ! p1
  r13 = ( m1   -p12-p12*IEPS )/(sm1*sm3) ! p1+p2
  k14 = ( m1+m4-p4 -p4 *IEPS )/(sm1*sm4) ! p1+p2+p3
  r23 = ( m2   -p2 -p2 *IEPS )/(sm2*sm3) ! p2
  k24 = ( m2+m4-p23-p23*IEPS )/(sm2*sm4) ! p2+p3
  r34 = (    m4-p3 -p3 *IEPS )/(sm3*sm4) ! p3
!
  r13zero = (abs(areal(r13))+abs(aimag(r13)).lt.neglig(prcpar))
  r23zero = (abs(areal(r23))+abs(aimag(r23)).lt.neglig(prcpar))
  r34zero = (abs(areal(r34))+abs(aimag(r34)).lt.neglig(prcpar))
!
  if (r13zero) then
    if     (r23zero) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf33: ' &
       ,'m4=p4 and m3=p12, returning 0'
      rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
      return
    elseif (r34zero) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf33: ' &
       ,'m2=p1 and m3=p12, returning 0'
      rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
      return
    endif
  endif
!
  call rfun( r12,d12 ,k12 )
  call rfun( r14,d14 ,k14 )
  call rfun( r24,d24 ,k24 )
!
  aa = r34/r24 - r23
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf33: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = -r13*d24 + k12*r34 - k14*r23
  cc = k12*r13 + r24*r34 - k14*r24*r13 - r23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
  x1 = -x1
  x2 = -x2
!
  qx1 = qonv(x1 ,1 ) ! x1 SHOULD HAVE im. part
  qx2 = qonv(x2 ,1 ) ! x2 SHOULD HAVE im. part
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  qy1 = qx1/q24
  qy2 = qx2/q24
  rslt(0) = rslt(0) + li2c2( qy1*q12 ,qy2*q12 )/r24*r12
  rslt(0) = rslt(0) + li2c2( qy1/q12 ,qy2/q12 )/r24/r12
  rslt(0) = rslt(0) - li2c2( qx1*q14 ,qx2*q14 )*r14
  rslt(0) = rslt(0) - li2c2( qx1/q14 ,qx2/q14 )/r14
!
  if (.not.r13zero) then
    if (.not.r23zero) then
      qss = q23/q13/q24
      rslt(0) = rslt(0) - li2c2( qx1*qss ,qx2*qss )*r23/(r13*r24)
    endif
    if (.not.r34zero) then
      qss = q34/q13
      rslt(0) = rslt(0) + li2c2( qx1*qss ,qx2*qss )*r34/r13
    endif
  else
    rslt(0) = rslt(0) - logc2( qx1/qx2 )*logc( q23/q24/q34 )/x2 
  endif
!
  rslt(0) = -rslt(0)/(aa*sm1*sm2*sm3*sm4)
  end subroutine


  subroutine boxf4( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
!*******************************************************************
! Finite 1-loop scalar 4-point function with all internal masses
! non-zero. Based on the formulas from
! A. Denner, U. Nierste, R. Scharf, Nucl.Phys.B367(1991)637-656
!*******************************************************************
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2) 
  complex(kindr2) &   
    ,intent(in) :: p1,p2,p3,p4,p12,p23,m1,m2,m3,m4
  type(qmplx_type) :: q12,q13,q14,q23,q24,q34,qx1,qx2,qy1,qy2,qtt
  complex(kindr2) &   
    :: sm1,sm2,sm3,sm4 ,aa,bb,cc,dd,x1,x2,tt &
                    ,k12,k13,k14,k23,k24,k34 &
                    ,r12,r13,r14,r23,r24,r34 &
                    ,d12,d13,d14,d23,d24,d34
  real(kindr2) &  
    :: h1,h2
!
  sm1 = mysqrt(m1)
  sm2 = mysqrt(m2)
  sm3 = mysqrt(m3)
  sm4 = mysqrt(m4)
!
  k12 = ( m1+m2-p1 -p1 *IEPS)/(sm1*sm2) ! p1
  k13 = ( m1+m3-p12-p12*IEPS)/(sm1*sm3) ! p1+p2
  k14 = ( m1+m4-p4 -p4 *IEPS)/(sm1*sm4) ! p1+p2+p3
  k23 = ( m2+m3-p2 -p2 *IEPS)/(sm2*sm3) ! p2
  k24 = ( m2+m4-p23-p23*IEPS)/(sm2*sm4) ! p2+p3
  k34 = ( m3+m4-p3 -p3 *IEPS)/(sm3*sm4) ! p3
!
  call rfun( r12,d12 ,k12 )
  call rfun( r13,d13 ,k13 )
  call rfun( r14,d14 ,k14 )
  call rfun( r23,d23 ,k23 )
  call rfun( r24,d24 ,k24 )
  call rfun( r34,d34 ,k34 )
!
  aa = k34/r24 + r13*k12 - k14*r13/r24 - k23
!
  if (aa.eq.CZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxf4: ' &
       ,'threshold singularity, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  bb = d13*d24 + k12*k34 - k14*k23
  cc = k12/r13 + r24*k34 - k14*r24/r13 - k23
  call solabc(x1,x2,dd ,aa,bb,cc ,0)
!
  h1 = areal(k23 - r13*k12 - r24*k34 + r13*r24*k14)
  h2 = h1*areal(aa)*areal(x1)
  h1 = h1*areal(aa)*areal(x2)
!
  qx1 = qonv(-x1,-h1) ! x1 should have im. part
  qx2 = qonv(-x2,-h2) ! x2 should have im. part
  q12 = qonv(r12,-1)
  q13 = qonv(r13,-1)
  q14 = qonv(r14,-1)
  q23 = qonv(r23,-1)
  q24 = qonv(r24,-1)
  q34 = qonv(r34,-1)
!
  rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
  qy1 = qx1/q24
  qy2 = qx2/q24
  rslt(0) = rslt(0) + ( li2c2( qy1*q12 ,qy2*q12 )*r12 &
                      + li2c2( qy1/q12 ,qy2/q12 )/r12 )/r24
  tt = r13/r24
  qtt = qonv(tt,-areal(r24) )
  qy1 = qx1*qtt
  qy2 = qx2*qtt
  rslt(0) = rslt(0) - ( li2c2( qy1*q23 ,qy2*q23 )*r23 &
                      + li2c2( qy1/q23 ,qy2/q23 )/r23 )*tt
  qy1 = qx1*q13
  qy2 = qx2*q13
  rslt(0) = rslt(0) + ( li2c2( qy1*q34 ,qy2*q34 )*r34 &
                      + li2c2( qy1/q34 ,qy2/q34 )/r34 )*r13
!
  rslt(0) = rslt(0) - ( li2c2( qx1*q14 ,qx2*q14 )*r14 &
                      + li2c2( qx1/q14 ,qx2/q14 )/r14 )
!
  rslt(0) = -rslt(0)/(aa*sm1*sm2*sm3*sm4)
  end subroutine

end module


module avh_olo_qp_boxc
   use avh_olo_units
   use avh_olo_qp_prec
   use avh_olo_qp_auxfun
   use avh_olo_qp_qmplx
   implicit none
   private
   public :: boxc

contains

   subroutine boxc( rslt ,pp_in ,mm_in ,ap_in ,smax )
!*******************************************************************
! Finite 1-loop scalar 4-point function for complex internal masses
! Based on the formulas from
!   Dao Thi Nhung and Le Duc Ninh, arXiv:0902.0325 [hep-ph]
!   G. 't Hooft and M.J.G. Veltman, Nucl.Phys.B153:365-401,1979 
!*******************************************************************
   use avh_olo_qp_box ,only: base,casetable,ll=>permtable
  complex(kindr2) &   
     ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
     ,intent(in)  :: pp_in(6),mm_in(4)
  real(kindr2) &  
     ,intent(in)  :: ap_in(6),smax
  complex(kindr2) &   
     :: pp(6),mm(4)
  real(kindr2) &  
     :: ap(6),aptmp(6),rem,imm,hh
  complex(kindr2) &   
     :: a,b,c,d,e,f,g,h,j,k,dpe,epk,x1,x2,sdnt,o1,j1,e1 &
       ,dek,dpf,def,dpk,abc,bgj,jph,cph
   integer :: icase,jcase,ii
!
   rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
!
   hh = neglig(prcpar)*smax
   do ii=1,6
     if (ap_in(ii).ge.hh) then ;ap(ii)=ap_in(ii)
                          else ;ap(ii)=0
     endif
   enddo
!
   do ii=1,4
     if (ap(ii).eq.RZRO) then ;pp(ii)=0
                         else ;pp(ii)=pp_in(ii)
     endif
   enddo
   if (ap(5).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
       ,' |s| too small, putting it by hand'
     ap(5) = hh
     pp(5) = acmplx(sign(hh,areal(pp_in(5))))
   else
     pp(5) = pp_in(5)
   endif
   if (ap(6).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
       ,' |t| too small, putting it by hand'
     ap(6) = hh
     pp(6) = acmplx(sign(hh,areal(pp_in(6))))
   else
     pp(6) = pp_in(6)
   endif
!
   do ii=1,4
     rem = areal(mm_in(ii))
     imm = aimag(mm_in(ii))
     hh = EPSN*abs(rem)
     if (abs(imm).lt.hh) imm = -hh
     mm(ii) = acmplx(rem,imm)
   enddo
!
   icase = 0
   do ii=1,4
     if (ap(ii).gt.RZRO) icase = icase + base(ii)
   enddo
!
   if (icase.lt.15) then
! at least one exernal mass equal zero
     jcase = casetable(icase)
     if (jcase.eq.0.or.jcase.eq.1.or.jcase.eq.5) then
! two opposite masses equal zero
       a = pp(ll(5,icase)) - pp(ll(1,icase))
       c = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(3,icase))
       g = pp(ll(2,icase))
       h = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
       d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
       e = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
       f = mm(ll(4,icase))
       k = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
       dpe = (mm(ll(1,icase)) - mm(ll(4,icase))) - pp(ll(4,icase))
       dpk = (mm(ll(2,icase)) - mm(ll(4,icase))) - pp(ll(6,icase))
       dpf = mm(ll(3,icase)) - pp(ll(3,icase))
       rslt(0) = t13fun( a,c,g,h ,d,e,f,k ,dpe,dpk,dpf )
     else
       a = pp(ll(3,icase))
       b = pp(ll(2,icase))
       c = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
       h = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(6,icase)) + pp(ll(2,icase))
       j = pp(ll(5,icase)) - pp(ll(1,icase)) - pp(ll(2,icase))
       d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
       e = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
       k = (mm(ll(1,icase)) - mm(ll(2,icase))) + pp(ll(6,icase)) - pp(ll(4,icase))
       f = mm(ll(4,icase))
       cph = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(3,icase))
       dpe = (mm(ll(2,icase)) - mm(ll(4,icase))) - pp(ll(6,icase))
       epk = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
       dek = (mm(ll(1,icase)) - mm(ll(4,icase))) - pp(ll(4,icase))
       dpf = mm(ll(3,icase)) - pp(ll(3,icase))
       rslt(0) = tfun( a,b  ,c  ,h,j ,d,e  ,f ,k ,dpe,dpf ) &
               - tfun( a,b+j,cph,h,j ,d,epk,f ,k ,dek,dpf )
     endif
   else
! no extenal mass equal zero
     if    (areal((pp(5)-pp(1)-pp(2))**2-4*pp(1)*pp(2)).gt.RZRO)then ;icase=0 !12, no permutation
     elseif(areal((pp(6)-pp(2)-pp(3))**2-4*pp(2)*pp(3)).gt.RZRO)then ;icase=8 !23, 1 cyclic permutation
     elseif(areal((pp(4)-pp(5)-pp(3))**2-4*pp(5)*pp(3)).gt.RZRO)then ;icase=4 !34, 2 cyclic permutations
     elseif(areal((pp(4)-pp(1)-pp(6))**2-4*pp(1)*pp(6)).gt.RZRO)then ;icase=2 !41, 3 cyclic permutations
     else
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
         ,'no positive lambda, returning 0'
       return
     endif
     a = pp(ll(3,icase))
     b = pp(ll(2,icase))
     g = pp(ll(1,icase))
     c = pp(ll(6,icase)) - pp(ll(2,icase)) - pp(ll(3,icase))
     h = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(6,icase)) + pp(ll(2,icase))
     j = pp(ll(5,icase)) - pp(ll(1,icase)) - pp(ll(2,icase))
     d = (mm(ll(3,icase)) - mm(ll(4,icase))) - pp(ll(3,icase))
     e = (mm(ll(2,icase)) - mm(ll(3,icase))) - pp(ll(6,icase)) + pp(ll(3,icase))
     k = (mm(ll(1,icase)) - mm(ll(2,icase))) + pp(ll(6,icase)) - pp(ll(4,icase))
     f = mm(ll(4,icase))
     abc = pp(ll(6,icase))
     bgj = pp(ll(5,icase))
     jph = pp(ll(4,icase)) - pp(ll(1,icase)) - pp(ll(6,icase))
     cph = pp(ll(4,icase)) - pp(ll(5,icase)) - pp(ll(3,icase))
     dpe = (mm(ll(2,icase)) - mm(ll(4,icase))) - pp(ll(6,icase))
     epk = (mm(ll(1,icase)) - mm(ll(3,icase))) + pp(ll(3,icase)) - pp(ll(4,icase))
     dek = (mm(ll(1,icase)) - mm(ll(4,icase))) - pp(ll(4,icase))
     dpf = mm(ll(3,icase)) - pp(ll(3,icase))
     def = mm(ll(2,icase)) - pp(ll(6,icase))
     call solabc( x1,x2 ,sdnt ,g,j,b ,0 )
     if (aimag(sdnt).ne.RZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop boxc: ' &
         ,'no real solution for alpha, returning 0'
       return
     endif
!BAD        if (abs(areal(x1)).gt.abs(areal(x2))) then
     if (abs(areal(x1)).lt.abs(areal(x2))) then !BETTER
       sdnt = x1
       x1 = x2
       x2 = sdnt
     endif
     o1 = 1-x1
     j1 = j+2*g*x1
     e1 = e+k*x1
     rslt(0) =   -tfun( abc,g  ,jph,c+2*b+(h+j)*x1, j1   ,dpe,k  ,f,e1 ,dek,def ) &
             + o1*tfun( a  ,bgj,cph,c+h*x1        , o1*j1,d  ,epk,f,e1 ,dek,dpf ) &
             + x1*tfun( a  ,b  ,c  ,c+h*x1        ,-j1*x1,d  ,e  ,f,e1 ,dpe,dpf )
   endif
   end subroutine


   function t13fun( aa,cc,gg,hh ,dd,ee,ff,jj ,dpe,dpj,dpf ) result(rslt)
!*******************************************************************
! /1   /x                             y
! | dx |  dy -----------------------------------------------------
! /0   /0    (gy^2 + hxy + dx + jy + f)*(ax^2 + cxy + dx + ey + f)
!
! jj should have negative imaginary part
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: aa,cc,gg,hh ,dd,ee,ff,jj ,dpe,dpj,dpf
  complex(kindr2) &   
     :: rslt ,kk,ll,nn,y1,y2,sdnt
!
!
   kk = hh*aa - cc*gg
   ll = aa*dd + hh*ee - dd*gg - cc*jj
   nn = dd*(ee - jj) + (hh - cc)*(ff-IEPS*abs(areal(ff)))
   call solabc( y1,y2 ,sdnt ,kk,ll,nn ,0 )
!
   rslt = - s3fun( y1,y2 ,CZRO,CONE ,aa   ,ee+cc,dpf ) &
          + s3fun( y1,y2 ,CZRO,CONE ,gg   ,jj+hh,dpf ) &
          - s3fun( y1,y2 ,CZRO,CONE ,gg+hh,dpj  ,ff  ) &
          + s3fun( y1,y2 ,CZRO,CONE ,aa+cc,dpe  ,ff  )
!
   rslt = rslt/kk
   end function


   function t1fun( aa,cc,gg,hh ,dd,ee,ff,jj ,dpe ) result(rslt)
!*******************************************************************
! /1   /x                         1
! | dx |  dy ----------------------------------------------
! /0   /0    (g*x + h*x + j)*(a*x^2 + c*xy + d*x + e*y + f)
!
! jj should have negative imaginary part
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: aa,cc,gg,hh ,dd,ee,ff,jj,dpe
  complex(kindr2) &   
     ::rslt ,kk,ll,nn,y1,y2,sdnt
!
!
   kk = hh*aa - cc*gg
   ll = hh*dd - cc*jj - ee*gg
   nn = hh*(ff-IEPS*abs(areal(ff))) - ee*jj
   call solabc( y1,y2 ,sdnt ,kk,ll,nn ,0 )
!
   rslt = - s3fun( y1,y2 ,CZRO,CONE ,aa+cc,dpe  ,ff ) &
          + s3fun( y1,y2 ,CZRO,CONE ,CZRO ,gg+hh,jj ) &
          - s3fun( y1,y2 ,CZRO,CONE ,CZRO ,gg   ,jj ) &
          + s3fun( y1,y2 ,CZRO,CONE ,aa   ,dd   ,ff )
!
   rslt = rslt/kk
   end function


   function tfun( aa,bb,cc ,gin,hin ,dd,ee,ff ,jin ,dpe ,dpf ) result(rslt)
!*******************************************************************
! /1   /x                             1
! | dx |  dy ------------------------------------------------------
! /0   /0    (g*x + h*x + j)*(a*x^2 + b*y^2 + c*xy + d*x + e*y + f)
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: aa,bb,cc ,gin,hin ,dd,ee,ff ,jin ,dpe ,dpf
  complex(kindr2) &   
     :: rslt ,gg,hh,jj,zz(2),beta,tmpa(2),tmpb(2) &
       ,tmpc(2),kiz(2),ll,nn,kk,y1,y2,yy(2,2),sdnt
  real(kindr2) &  
     :: ab1,ab2,ac1,ac2,abab,acac,abac,det,ap1,ap2 &
                  ,apab,apac,x1(2,2),x2(2,2),xmin
   integer :: iz,iy,izmin,sj
   logical :: pp(2,2),p1,p2
!
!
   sj = sgnIm(jin,-1)
   gg = -sj*gin
   hh = -sj*hin
   jj = -sj*jin
!
   if     (bb.eq.CZRO) then
     rslt = -sj*t1fun( aa,cc,gg,hh ,dd,ee,ff,jj ,dpe )
     return
   elseif (aa.eq.CZRO) then
     rslt = -sj*t1fun( bb+cc,-cc,-gg-hh,gg, -dpe-2*(bb+cc),dd+cc &
                      ,dpe+bb+cc+ff,gg+hh+jj ,-ee-2*bb-cc )
     return
   endif
!
   call solabc( zz(1),zz(2) ,sdnt ,bb,cc,aa ,0 )
   if (abs(zz(1)).gt.abs(zz(2))) then
     beta = zz(1)
     zz(1) = zz(2)
     zz(2) = beta
   endif
!
   do iz=1,2
     beta = zz(iz)
     tmpa(iz) = gg + beta*hh
     tmpb(iz) = cc + 2*beta*bb
     tmpc(iz) = dd + beta*ee
     kiz(iz) =        bb*tmpa(iz)               - hh*tmpb(iz)
     ll      =        ee*tmpa(iz) - hh*tmpc(iz) - jj*tmpb(iz)
     nn      = (ff-IEPS*abs(areal(ff)))*tmpa(iz) - jj*tmpc(iz)
     call solabc( yy(iz,1),yy(iz,2) ,sdnt ,kiz(iz),ll,nn ,0 )
     if (abs(aimag(beta)).ne.RZRO) then
       ab1 = areal(-beta)
       ab2 = aimag(-beta)
       ac1 = ab1+1 !areal(1-beta)
       ac2 = ab2   !aimag(1-beta)
       abab = ab1*ab1 + ab2*ab2
       acac = ac1*ac1 + ac2*ac2
       abac = ab1*ac1 + ab2*ac2
       det = abab*acac - abac*abac
       do iy=1,2
         ap1 = areal(yy(iz,iy))
         ap2 = aimag(yy(iz,iy))
         apab = ap1*ab1 + ap2*ab2
         apac = ap1*ac1 + ap2*ac2
         x1(iz,iy) = ( acac*apab - abac*apac )/det
         x2(iz,iy) = (-abac*apab + abab*apac )/det
       enddo
     else
       do iy=1,2
         x1(iz,iy) = -1
         x2(iz,iy) = -1
       enddo
     endif
   enddo
   xmin = 1
   izmin = 2
   do iz=1,2
   do iy=1,2
     if ( x1(iz,iy).ge.RZRO.and.x2(iz,iy).ge.RZRO &
                 .and.x1(iz,iy)+x2(iz,iy).le.RONE ) then
       pp(iz,iy) = .true.
       if (x1(iz,iy).lt.xmin) then
         xmin = x1(iz,iy)
         izmin = iz
       endif
       if (x2(iz,iy).lt.xmin) then
         xmin = x2(iz,iy)
         izmin = iz
       endif
     else
       pp(iz,iy) = .false.
     endif
   enddo
   enddo
   iz = izmin+1
   if (iz.eq.3) iz = 1
!
   beta = zz(iz)
   kk = kiz(iz)
   y1 = yy(iz,1)
   y2 = yy(iz,2)
   p1 = pp(iz,1)
   p2 = pp(iz,2)
!
   rslt =   s3fun( y1,y2 ,beta ,CONE      ,CZRO    ,hh   ,gg+jj  ) &
          - s3fun( y1,y2 ,CZRO ,CONE-beta ,CZRO    ,gg+hh,   jj  ) &
          + s3fun( y1,y2 ,CZRO ,    -beta ,CZRO    ,gg   ,   jj  ) &
          - s3fun( y1,y2 ,beta ,CONE      ,bb      ,cc+ee,aa+dpf ) &
          + s3fun( y1,y2 ,CZRO ,CONE-beta ,aa+bb+cc,dpe  ,ff     ) &
          - s3fun( y1,y2 ,CZRO ,    -beta ,aa      ,dd   ,ff     )
!
   sdnt = plnr( y1,y2 ,p1,p2, tmpa(iz),tmpb(iz),tmpc(iz) )
   if (aimag(beta).le.RZRO) then ;rslt = rslt + sdnt
                            else ;rslt = rslt - sdnt
   endif
!
   rslt = -sj*rslt/kk
   end function


   function s3fun( y1i,y2i ,dd,ee ,aa,bb,cin ) result(rslt)
!*******************************************************************
! Calculate
!            ( S3(y1i) - S3(y2i) )/( y1i - y2i )
! where
!               /1    ee * ln( aa*x^2 + bb*x + cc )
!       S3(y) = |  dx -----------------------------
!               /0           ee*x - y - dd
!
! y1i,y2i should have a non-zero imaginary part
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) ::  y1i,y2i ,dd,ee ,aa,bb,cin
  complex(kindr2) &   
     :: rslt ,y1,y2,fy1y2,z1,z2,tmp,cc
  real(kindr2) &  
     ::rea,reb,rez1,rez2,imz1,imz2,simc,hh
!
!
   if (ee.eq.CZRO) then
     rslt = 0
     return
   endif
!
   cc = cin
   rea = abs(aa)
   reb = abs(bb)
   simc = abs(cc)
   if (simc.lt.10*neglig(prcpar)*min(rea,reb)) cc = 0
!
   simc = aimag(cc)
   if (simc.eq.RZRO) then
     simc = aimag(bb)
     if (simc.eq.RZRO) simc = -1
   endif
   simc = sgnRe(simc)
!
   y1 = (dd+y1i)/ee
   y2 = (dd+y2i)/ee
   if (aimag(y1).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'y1 has zero imaginary part'
   endif
   if (aimag(y2).eq.RZRO) then
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'y2 has zero imaginary part'
   endif
   fy1y2 = r0fun( y1,y2 )
!
   if     (aa.ne.CZRO) then
!
!     call solabc( z1,z2 ,tmp ,aa,bb,cc ,0 )
     call solabc_rcc( z1,z2 ,areal(aa),bb,cc )
     rea  = sgnRe(aa)
     rez1 = areal(z1)
     rez2 = areal(z2) 
     imz1 = aimag(z1) ! sign(Im(a*z1*z2)) = simc
     imz2 = aimag(z2)
     hh = abs(EPSN2*rez1)
!     if (abs(imz1).lt.EPSN*hh) imz1 = simc*rea*sgnRe(rez2)*hh
     if (imz1.eq.RZRO) imz1 = simc*rea*sgnRe(rez2)*hh
     hh = abs(EPSN2*rez2)
!     if (abs(imz2).lt.EPSN*hh) imz2 = simc*rea*sgnRe(rez1)*hh
     if (imz2.eq.RZRO) imz2 = simc*rea*sgnRe(rez1)*hh
     z1 = acmplx( rez1,imz1)
     z2 = acmplx( rez2,imz2)
     rslt = fy1y2 * ( logc(qonv(aa,simc)) &
                    + eta3( -z1,-imz1,-z2,-imz2,CZRO,simc*rea ) ) &
          + r1fun( z1,y1,y2,fy1y2 ) &
          + r1fun( z2,y1,y2,fy1y2 )
!
   elseif (bb.ne.CZRO) then
!
     z1 = -cc/bb ! - i|eps|Re(b)
     reb  = areal(bb)
     rez1 = areal(z1)
     imz1 = aimag(z1)
     if (abs(imz1).eq.RZRO) then
       imz1 = -simc*reb*abs(EPSN2*rez1/reb)
       z1 = acmplx( rez1,imz1)
     endif
     rslt = fy1y2 * ( logc(qonv(bb,simc)) &
                    + eta3(bb,simc ,-z1,-imz1 ,cc,simc) ) &
          + r1fun( z1,y1,y2,fy1y2 )
!
   elseif (cc.ne.CZRO) then
!
     rslt = logc( qonv(cc,simc) )*fy1y2
!
   else!if (aa=bb=cc=0)
!
     if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop s3fun: ' &
       ,'cc equal zero, returning 0'
     rslt = 0
!
   endif
!
   rslt = rslt/ee
   end function


   function r1fun( zz,y1,y2,fy1y2 ) result(rslt)
!*******************************************************************
! calculates  ( R1(y1,z) - R1(y2,z) )/( y1 - y2 )
! where
!                          /     / 1-y \       / 1-z \ \
!      R1(y,z) = ln(y-z) * | log |-----| - log |-----| |
!                          \     \ -y  /       \ -z  / / 
!
!                      /    y-z \       /    y-z \
!                - Li2 |1 - ----| + Li2 |1 - ----|
!                      \    -z  /       \    1-z /
!
!                                     / 1-y1 \       / 1-y2 \
!                                 log |------| - log |------| 
! input fy1y2 should be equal to      \  -y1 /       \  -y2 /
!                                 ---------------------------
!                                           y1 - y2
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: y1,y2,zz,fy1y2
  complex(kindr2) &   
     :: rslt ,oz
   type(qmplx_type) :: q1z,q2z,qq
  real(kindr2) &  
     :: h12,hz1,hz2,hzz,hoz
   logical :: zzsmall,ozsmall
!
!
   oz = 1-zz
   h12 = abs(y1-y2)
   hz1 = abs(y1-zz)
   hz2 = abs(y2-zz)
   hzz = abs(zz)
   hoz = abs(oz)
   q1z = qonv(y1-zz)
   q2z = qonv(y2-zz)
!
   zzsmall = .false.
   ozsmall = .false.
   if     (hzz.lt.hz1.and.hzz.lt.hz2.and.hzz.lt.hoz) then ! |z| < |y1-z|,|y2-z|
     zzsmall = .true.
     rslt = fy1y2*logc( q1z ) &
          - ( logc(q1z*q2z)/2 + logc(qonv((y2-1)/y2)) &
                                     - logc(qonv(oz)) )*logc2(q1z/q2z)/(y2-zz)
   elseif (hoz.lt.hz1.and.hoz.lt.hz2) then ! |1-z| < |y1-z|,|y2-z|
     ozsmall = .true.
     rslt = fy1y2*logc( q1z ) &
          - (-logc(q1z*q2z)/2 + logc(qonv((y2-1)/y2)) &
                                    + logc(qonv(-zz)) )*logc2(q1z/q2z)/(y2-zz)
   elseif (h12.le.hz2.and.hz2.le.hz1) then ! |y1-y2| < |y2-z| < |y1-z|
     rslt = fy1y2*logc( q1z ) - r0fun( y2,zz )*logc2( q1z/q2z )        
   elseif (h12.le.hz1.and.hz1.le.hz2) then ! |y1-y2| < |y2-z| < |y1-z|
     rslt = fy1y2*logc( q2z ) - r0fun( y1,zz )*logc2( q2z/q1z )        
   else!if(hz1.lt.h12.or.hz2.lt.h12) then ! |y2-z|,|y1-z| < |y1-y2|
     rslt = 0
     if (hz1.ne.RZRO) rslt = rslt + (y1-zz)*logc( q1z )*r0fun( y1,zz )
     if (hz2.ne.RZRO) rslt = rslt - (y2-zz)*logc( q2z )*r0fun( y2,zz )
     rslt = rslt/(y1-y2)
   endif
!
   if (zzsmall) then ! |z| < |y1-z|,|y2-z|
     qq  = qonv(-zz)
     rslt = rslt + ( li2c( qq/q1z ) - li2c( qq/q2z ) )/(y1-y2)
   else
     qq  = qonv(-zz)
     rslt = rslt + li2c2( q1z/qq ,q2z/qq )/zz
   endif
!
   if (ozsmall) then ! |1-z| < |y1-z|,|y2-z|
     qq  = qonv(oz)
     rslt = rslt - ( li2c( qq/q1z ) - li2c( qq/q2z ) )/(y1-y2)
   else
     qq = qonv(oz)
     rslt = rslt + li2c2( q1z/qq ,q2z/qq )/oz
   endif
   end function


   function r0fun( y1,y2 ) result(rslt)
!*******************************************************************
!      / 1-y1 \       / 1-y2 \
!  log |------| - log |------| 
!      \  -y1 /       \  -y2 /
!  ---------------------------
!            y1 - y2
!
! y1,y2 should have non-zero imaginary parts
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: y1,y2
  complex(kindr2) &   
     :: rslt ,oy1,oy2
!
   oy1 = 1-y1
   oy2 = 1-y2
   rslt = logc2( qonv(-y2)/qonv(-y1) )/y1 &
        + logc2( qonv(oy2)/qonv(oy1) )/oy1
   end function


   function plnr( y1,y2 ,p1,p2 ,aa,bb,cc ) result(rslt)
!*******************************************************************
!                   /   a    \          /   a    \
!            p1*log |--------| - p2*log |--------| 
!                   \ b*y1+c /          \ b*y2+c /
! 2*pi*imag* -------------------------------------
!                           y1 - y2
! 
! p1,p2 are logical, to be interpreted as 0,1 in the formula above 
!*******************************************************************
  complex(kindr2) &   
     ,intent(in) :: y1,y2 ,aa,bb,cc
   logical         ,intent(in) :: p1,p2
  complex(kindr2) &   
     :: rslt ,x1,x2,xx
   type(qmplx_type) :: q1,q2
!
   if (p1) then
     x1 = bb*y1 + cc
     xx = aa/x1
     if (aimag(xx).eq.RZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop plnr: ' &
         ,'aa/x1 has zero imaginary part'
     endif
     q1 = qonv(xx)
   endif
   if (p2) then
     x2 = bb*y2 + cc
     xx = aa/x2
     if (aimag(xx).eq.RZRO) then
       if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop plnr: ' &
         ,'aa/x2 has zero imaginary part'
     endif
     q2 = qonv(xx)
   endif
   if (p1) then
     if (p2) then
       rslt = logc2( q2/q1 ) * 2*IPI*bb/x2
     else
       rslt = logc( q1 ) * 2*IPI/(y1-y2)
     endif
   elseif (p2) then
     rslt = logc( q2 ) * 2*IPI/(y2-y1) ! minus sign
   else
     rslt = 0
   endif
   end function


end module


module avh_olo_qp
  use avh_olo_units
  use avh_olo_qp_print
  use avh_olo_qp_prec
!
  implicit none
  private
  public :: olo_unit ,olo_scale ,olo_onshell ,olo_setting
  public :: olo_precision
  public :: olo_a0 ,olo_b0 ,olo_db0 ,olo_b11 ,olo_c0 ,olo_d0
  public :: olo_an ,olo_bn
  public :: olo
  public :: olo_get_scale ,olo_get_onshell ,olo_get_precision
!
  integer ,public ,parameter :: olo_kind=kindr2    
!
  real(kindr2) &  
         ,save :: onshellthrs
  logical,save :: nonzerothrs = .false.
!
  real(kindr2) &  
         ,save :: muscale
!
  character(99) ,parameter :: warnonshell=&
       'it seems you forgot to put some input explicitly on shell. ' &
     //'You may  call olo_onshell  to cure this.'
!
  logical ,save :: initz=.true.
!
  interface olo_a0
    module procedure a0_r,a0rr,a0_c,a0cr
  end interface 
  interface olo_an
    module procedure an_r,anrr,an_c,ancr
  end interface 
  interface olo_b0
    module procedure b0rr,b0rrr,b0rc,b0rcr,b0cc,b0ccr
  end interface 
  interface olo_db0
    module procedure db0rr,db0rrr,db0rc,db0rcr,db0cc,db0ccr
  end interface 
  interface olo_b11
    module procedure b11rr,b11rrr,b11rc,b11rcr,b11cc,b11ccr
  end interface 
  interface olo_bn
    module procedure bnrr,bnrrr,bnrc,bnrcr,bncc,bnccr
  end interface 
  interface olo_c0
    module procedure c0rr,c0rrr,c0rc,c0rcr,c0cc,c0ccr
  end interface 
  interface olo_d0
    module procedure d0rr,d0rrr,d0rc,d0rcr,d0cc,d0ccr
  end interface 
!
  interface olo
    module procedure a0_r,a0rr,a0_c,a0cr
    module procedure an_r,anrr,an_c,ancr
    module procedure b0rr,b0rrr,b0rc,b0rcr,b0cc,b0ccr
    module procedure b11rr,b11rrr,b11rc,b11rcr,b11cc,b11ccr
    module procedure bnrr,bnrrr,bnrc,bnrcr,bncc,bnccr
    module procedure c0rr,c0rrr,c0rc,c0rcr,c0cc,c0ccr
    module procedure d0rr,d0rrr,d0rc,d0rcr,d0cc,d0ccr
  end interface 

contains

 
  subroutine init( ndec )
!*******************************************************************
!*******************************************************************
  use avh_olo_version
  integer,optional,intent(in) :: ndec
!
  call olo_version
!
  initz = .false.
!
  if (present(ndec)) then
    call olo_precision( ndec )
  else
    call olo_precision( 15 )
  endif
!
  onshellthrs = 0
  muscale = 1
  if (.not.nonzerothrs) onshellthrs = neglig(prcpar)
!
  end subroutine
 
 
  recursive subroutine olo_precision( ndec )
!*******************************************************************
!*******************************************************************
  use avh_olo_qp_olog  ,only: update_olog
  use avh_olo_qp_dilog ,only: update_dilog
  use avh_olo_qp_bnlog ,only: update_bnlog
  integer ,intent(in) :: ndec
  logical :: newprc
  if (initz) then
    call init( ndec )
  else
    call set_precision( newprc )       
    if (newprc) then
      call update_olog
      call update_dilog
      call update_bnlog
    endif
    if (.not.nonzerothrs) onshellthrs = neglig(prcpar)
  endif
  end subroutine

 
  subroutine olo_unit( val ,message )
!*******************************************************************
!*******************************************************************
  integer     ,intent(in) :: val
  character(*),intent(in),optional :: message
  if (initz) call init
  if (present(message)) then ;call set_unit( message ,val )
  else                       ;call set_unit( 'all'   ,val )
  endif
  end subroutine
 
 
  subroutine olo_scale( val )
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) ,intent(in) :: val
  if (initz) call init
  muscale = convert(val)
  end subroutine
 
 
  subroutine olo_onshell( thrs )
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) ,intent(in) :: thrs
  if (initz) call init
  nonzerothrs = .true.
  onshellthrs = convert(thrs)
  end subroutine


  function olo_get_precision() result(rslt)
!*******************************************************************
!*******************************************************************
  use avh_olo_qp_prec ,only: ndecim,prcpar
  integer :: rslt
  if (initz) call init
  rslt = ndecim(prcpar)
  end function

  function olo_get_scale() result(rslt)
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) :: rslt
  if (initz) call init
  rslt = adble(muscale)
  end function

  function olo_get_onshell() result(rslt)
!*******************************************************************
!*******************************************************************
  real(kind(1d0)) :: rslt
  if (initz) call init
  rslt = adble(onshellthrs)
  end function


  subroutine olo_setting( iunit )
!*******************************************************************
!*******************************************************************
  integer,optional,intent(in) :: iunit
  integer :: nunit
  if (initz) call init
  nunit = munit
  if (present(iunit)) nunit = iunit
  if (nunit.le.0) return
!
  write(nunit,*) 'MESSAGE from OneLOop: real kind parameter =',trim(myprint(kindr2)) 
  write(nunit,*) 'MESSAGE from OneLOop: number of decimals  =',trim(myprint(ndecim(prcpar)))
!
  if (nonzerothrs) then
    write(nunit,*) 'MESSAGE from OneLOop: on-shell threshold =',trim(myprint(onshellthrs,12))
  else
    write(nunit,*) 'MESSAGE from OneLOop: on-shell threshold is not set'
  endif
!
  write(nunit,*) 'MESSAGE from OneLOop: default scale (mu, not mu^2) =',trim(myprint(muscale,12))
!
  end subroutine
 
 
!*******************************************************************
!
!           C   / d^(Dim)q
! rslt = ------ | -------- 
!        i*pi^2 / (q^2-mm)
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  mm = mass squared
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine a0_c( rslt ,mm )
!
  use avh_olo_qp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: mm
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine a0cr( rslt ,mm ,rmu )
!
  use avh_olo_qp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine a0_r( rslt ,mm  )
!
  use avh_olo_qp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: mm
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine a0rr( rslt ,mm ,rmu )
!
  use avh_olo_qp_bub ,only: tadp
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop a0: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadp( rslt ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    write(punit,*) 'a0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'a0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'a0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine


  subroutine an_c( rslt ,rank ,mm )
!
  use avh_olo_qp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: mm
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine

  subroutine ancr( rslt ,rank ,mm ,rmu )
!
  use avh_olo_qp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine

  subroutine an_r( rslt ,rank ,mm  )
!
  use avh_olo_qp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: mm
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = muscale 
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine

  subroutine anrr( rslt ,rank ,mm ,rmu )
!
  use avh_olo_qp_bub ,only: tadpn
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: mm
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss
  real(kindr2) &  
    :: am,hh,mulocal,mulocal2
  integer :: ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop An: '//warnonshell
  if (initz) call init
!
  mulocal = rmu     
!
  am = abs(mm)
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (am.lt.hh) am = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(am,mulocal2)
    if (RZRO.lt.am.and.am.lt.hh) write(wunit,*) warning
  endif
!
  ss = mm
  call tadpn( rslt ,rank ,ss ,am ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' mm:',trim(myprint(mm))
    do ii=0,rank/2
    write(punit,*) 'A(2,',trim(myprint(ii)),'):',trim(myprint(rslt(2,ii)))
    write(punit,*) 'A(1,',trim(myprint(ii)),'):',trim(myprint(rslt(1,ii)))
    write(punit,*) 'A(0,',trim(myprint(ii)),'):',trim(myprint(rslt(0,ii)))
    enddo
  endif
  end subroutine


!*******************************************************************
!
!           C   /      d^(Dim)q
! rslt = ------ | --------------------
!        i*pi^2 / [q^2-m1][(q+k)^2-m2]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps) * exp(gamma_Euler*eps)
!
! input:  pp = k^2, m1,m2 = mass squared
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine b0cc( rslt ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0ccr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rc( rslt ,pp ,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rcr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rr( rslt ,pp ,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine b0rrr( rslt ,pp ,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub0
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop b0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub0( rslt ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

!*******************************************************************
! Derivative of B0
!*******************************************************************
  subroutine db0cc( rslt ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: dbub0
  use avh_olo_qp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0ccr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: dbub0
  use avh_olo_qp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rc( rslt ,pp ,m1,m2 )
!
  use avh_olo_qp_bub ,only: dbub0
  use avh_olo_qp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rcr( rslt ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: dbub0
  use avh_olo_qp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop db0: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b0: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rr( rslt ,pp ,m1,m2 )
!
  use avh_olo_qp_bub ,only: dbub0
  use avh_olo_qp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine db0rrr( rslt ,pp ,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: dbub0
  use avh_olo_qp_olog ,only: olog
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2,ch
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2,ssr2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop db0: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (am1.gt.am2) then
    ch=r1 ; r1=r2 ; r2=ch
    hh=am1;am1=am2;am2=hh
  endif
  ssr2 = abs(ss-r2)
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
    if (ssr2.lt.hh) ssr2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,am2))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.ssr2.and.ssr2.lt.hh) write(wunit,*) warning
  endif
!
  rslt(0)=0;rslt(1)=0;rslt(2)=0
!
  if (am1.eq.RZRO.and.ssr2.eq.RZRO) then
    rslt(1) =-1/(2*ss)
    rslt(0) =-( 1 + olog(mulocal2/ss,0)/2 )/ss
  else
    call dbub0( rslt(0) ,ss,r1,r2 ,app,am1,am2 )
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'db0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'db0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'db0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine



!*******************************************************************
! Return the Papparino-Veltman functions b11,b00,b1,b0 , for
!
!      C   /      d^(Dim)q
!   ------ | -------------------- = b0
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
!      C   /    d^(Dim)q q^mu
!   ------ | -------------------- = p^mu b1
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
!      C   /  d^(Dim)q q^mu q^nu
!   ------ | -------------------- = g^{mu,nu} b00 + p^mu p^nu b11
!   i*pi^2 / [q^2-m1][(q+p)^2-m2]
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine b11cc( b11,b00,b1,b0 ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11ccr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rc( b11,b00,b1,b0 ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rcr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop b11: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rr( b11,b00,b1,b0 ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine

  subroutine b11rrr( b11,b00,b1,b0 ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub11
!
  complex(kindr2) &   
    ,intent(out) :: b11(0:2),b00(0:2),b1(0:2),b0(0:2)
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop b11: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  call bub11( b11,b00,b1,b0 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' pp:',trim(myprint(pp))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) 'b11(2):',trim(myprint(b11(2)))
    write(punit,*) 'b11(1):',trim(myprint(b11(1)))
    write(punit,*) 'b11(0):',trim(myprint(b11(0)))
    write(punit,*) 'b00(2):',trim(myprint(b00(2)))
    write(punit,*) 'b00(1):',trim(myprint(b00(1)))
    write(punit,*) 'b00(0):',trim(myprint(b00(0)))
    write(punit,*) ' b1(2):',trim(myprint(b1(2) ))
    write(punit,*) ' b1(1):',trim(myprint(b1(1) ))
    write(punit,*) ' b1(0):',trim(myprint(b1(0) ))
    write(punit,*) ' b0(2):',trim(myprint(b0(2) ))
    write(punit,*) ' b0(1):',trim(myprint(b0(1) ))
    write(punit,*) ' b0(0):',trim(myprint(b0(0) ))
  endif
  end subroutine


  subroutine bncc( rslt ,rank ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnccr( rslt ,rank ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  complex(kindr2) &   
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = areal(ss)
  if (aimag(ss).ne.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'ss has non-zero imaginary part, putting it to zero.'
    ss = acmplx( app )
  endif
  app = abs(app)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrc( rslt ,rank ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrcr( rslt ,rank ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  complex(kindr2) &   
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = areal(r1)
  hh  = aimag(r1)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r1 has positive imaginary part, switching its sign.'
    r1 = acmplx( am1 ,-hh )
  endif
  am1 = abs(am1) + abs(hh)
!
  am2 = areal(r2)
  hh  = aimag(r2)
  if (hh.gt.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'r2 has positive imaginary part, switching its sign.'
    r2 = acmplx( am2 ,-hh )
  endif
  am2 = abs(am2) + abs(hh)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrr( rslt ,rank ,pp,m1,m2 )
!
  use avh_olo_qp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine

  subroutine bnrrr( rslt ,rank ,pp,m1,m2 ,rmu )
!
  use avh_olo_qp_bub ,only: bub0,bub1,bub11,bub111,bub1111
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:,0:)   
  real(kindr2) &  
    ,intent(in)  :: pp
  real(kindr2) &  
    ,intent(in)  :: m1,m2
  real(kindr2) &  
   ,intent(in)  :: rmu       
  integer,intent(in) :: rank
!
  complex(kindr2) &   
    :: ss,r1,r2
  real(kindr2) &  
    :: app,am1,am2,hh,mulocal,mulocal2
  character(26+99) ,parameter :: warning=&
                     'WARNING from OneLOop bn: '//warnonshell
  if (initz) call init
  ss = pp
  r1 = m1
  r2 = m2
!
  app = abs(pp)
!
  am1 = abs(m1)
  am2 = abs(m2)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (nonzerothrs) then
    hh = onshellthrs
    if (app.lt.hh) app = 0
    if (am1.lt.hh) am1 = 0
    if (am2.lt.hh) am2 = 0
  elseif (wunit.gt.0) then
    hh = onshellthrs*max(app,max(am1,max(am2,mulocal2)))
    if (RZRO.lt.app.and.app.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am1.and.am1.lt.hh) write(wunit,*) warning
    if (RZRO.lt.am2.and.am2.lt.hh) write(wunit,*) warning
  endif
!
  if     (rank.eq.0) then
    call bub0( rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.1) then
    call bub1( rslt(:,1),rslt(:,0) &
              ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.2) then
    call bub11( rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
               ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.3) then
    call bub111( rslt(:,5),rslt(:,4) &
                ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  elseif (rank.eq.4) then
    call bub1111( rslt(:,8),rslt(:,7),rslt(:,6) &
                 ,rslt(:,5),rslt(:,4) &
                 ,rslt(:,3),rslt(:,2),rslt(:,1),rslt(:,0) &
                 ,ss,r1,r2 ,app,am1,am2 ,mulocal2 )
  else
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop Bn: ' &
      ,'rank=',rank,' not implemented'
  endif
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) 'pp:',trim(myprint(pp))
    write(punit,*) 'm1:',trim(myprint(m1))
    write(punit,*) 'm2:',trim(myprint(m2))
    if (rank.ge.0) then
    write(punit,*) 'b0(2):',trim(myprint(rslt(2,0) ))
    write(punit,*) 'b0(1):',trim(myprint(rslt(1,0) ))
    write(punit,*) 'b0(0):',trim(myprint(rslt(0,0) ))
    if (rank.ge.1) then
    write(punit,*) 'b1(2):',trim(myprint(rslt(2,1) ))
    write(punit,*) 'b1(1):',trim(myprint(rslt(1,1) ))
    write(punit,*) 'b1(0):',trim(myprint(rslt(0,1) ))
    if (rank.ge.2) then
    write(punit,*) 'b00(2):',trim(myprint(rslt(2,2)))
    write(punit,*) 'b00(1):',trim(myprint(rslt(1,2)))
    write(punit,*) 'b00(0):',trim(myprint(rslt(0,2)))
    write(punit,*) 'b11(2):',trim(myprint(rslt(2,3)))
    write(punit,*) 'b11(1):',trim(myprint(rslt(1,3)))
    write(punit,*) 'b11(0):',trim(myprint(rslt(0,3)))
    if (rank.ge.3) then
    write(punit,*) 'b001(2):',trim(myprint(rslt(2,4)))
    write(punit,*) 'b001(1):',trim(myprint(rslt(1,4)))
    write(punit,*) 'b001(0):',trim(myprint(rslt(0,4)))
    write(punit,*) 'b111(2):',trim(myprint(rslt(2,5)))
    write(punit,*) 'b111(1):',trim(myprint(rslt(1,5)))
    write(punit,*) 'b111(0):',trim(myprint(rslt(0,5)))
    if (rank.ge.4) then
    write(punit,*) 'b0000(2):',trim(myprint(rslt(2,6)))
    write(punit,*) 'b0000(1):',trim(myprint(rslt(1,6)))
    write(punit,*) 'b0000(0):',trim(myprint(rslt(0,6)))
    write(punit,*) 'b0011(2):',trim(myprint(rslt(2,7)))
    write(punit,*) 'b0011(1):',trim(myprint(rslt(1,7)))
    write(punit,*) 'b0011(0):',trim(myprint(rslt(0,7)))
    write(punit,*) 'b1111(2):',trim(myprint(rslt(2,8)))
    write(punit,*) 'b1111(1):',trim(myprint(rslt(1,8)))
    write(punit,*) 'b1111(0):',trim(myprint(rslt(0,8)))
    endif;endif;endif;endif;endif
  endif
  end subroutine


!*******************************************************************
! calculates
!               C   /               d^(Dim)q
!            ------ | ---------------------------------------
!            i*pi^2 / [q^2-m1] [(q+k1)^2-m2] [(q+k1+k2)^2-m3]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps)
!             * GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
!
! input:  p1=k1^2, p2=k2^2, p3=(k1+k2)^2,  m1,m2,m3=squared masses
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  subroutine olo_onshell  to find out how
! this routine decides when to return IR-divergent cases.
!*******************************************************************

  subroutine c0cc( rslt ,p1,p2,p3 ,m1,m2,m3 )
  use avh_olo_qp_tri
  use avh_olo_qp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
!
  complex(kindr2) &   
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0ccr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
  use avh_olo_qp_tri
  use avh_olo_qp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  complex(kindr2) &   
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rc( rslt ,p1,p2,p3 ,m1,m2,m3 )
  use avh_olo_qp_tri
  use avh_olo_qp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
!
  real(kindr2) &  
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rcr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
  use avh_olo_qp_tri
  use avh_olo_qp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(3)
  complex(kindr2) &   
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = areal(mm(ii))
    hh     = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rr( rslt ,p1,p2,p3 ,m1,m2,m3 )
  use avh_olo_qp_tri
  use avh_olo_qp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3
!
  real(kindr2) &  
    :: pp(3)
  real(kindr2) &  
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine c0rrr( rslt ,p1,p2,p3 ,m1,m2,m3 ,rmu )
  use avh_olo_qp_tri
  use avh_olo_qp_auxfun ,only: kallen
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(3)
  real(kindr2) &  
    :: mm(3)
  complex(kindr2) &   
    :: ss(3),rr(3),lambda
  real(kindr2) &  
    :: smax,ap(3),am(3),as(3),ar(3),hh,s1r2,s2r3,s3r3
  real(kindr2) &  
    :: mulocal,mulocal2
  integer :: icase,ii
  character(25+99) ,parameter :: warning=&
                     'WARNING from OneLOop c0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  smax = 0
!
  do ii=1,3
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,3
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop c0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,3
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,3
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  icase = 0
  do ii=1,3
    if (am(ii).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(permtable(1,icase)) ;as(1)=ap(permtable(1,icase))
  ss(2)=pp(permtable(2,icase)) ;as(2)=ap(permtable(2,icase))
  ss(3)=pp(permtable(3,icase)) ;as(3)=ap(permtable(3,icase))
  rr(1)=mm(permtable(1,icase)) ;ar(1)=am(permtable(1,icase))
  rr(2)=mm(permtable(2,icase)) ;ar(2)=am(permtable(2,icase))
  rr(3)=mm(permtable(3,icase)) ;ar(3)=am(permtable(3,icase))
  icase = casetable(icase)
!
  s1r2 = abs(ss(1)-rr(2))
  s2r3 = abs(ss(2)-rr(3))
  s3r3 = abs(ss(3)-rr(3))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r3.lt.hh) s3r3 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r3.and.s3r3.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.3) then
! 3 non-zero internal masses
    lambda = kallen( ss(1),ss(2),ss(3) )
    if (areal(lambda).lt.RZRO) then
      call trif3HV( rslt ,ss,rr ,as ,smax ,lambda )
    else
      call trif3( rslt ,ss(1),ss(2),ss(3) ,rr(1),rr(2),rr(3) )
    endif
  elseif (icase.eq.2) then
! 2 non-zero internal masses
    if (s1r2.ne.RZRO.or.s3r3.ne.RZRO) then
      call trif2( rslt ,ss(1),ss(2),ss(3) ,rr(2),rr(3) )
    else
      call tria4( rslt ,ss(2) ,rr(2),rr(3) ,mulocal2 )
    endif
  elseif (icase.eq.1) then
! 1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      call trif1( rslt ,ss(1),ss(2),ss(3), rr(3) )
    elseif (s2r3.ne.RZRO) then
      if   (s3r3.ne.RZRO) then
        call tria3( rslt ,ss(2),ss(3) ,rr(3) ,mulocal2 )
      else
        call tria2( rslt ,ss(2) ,rr(3) ,mulocal2 )
      endif
    elseif (s3r3.ne.RZRO) then
      call tria2( rslt ,ss(3) ,rr(3) ,mulocal2 )
    else
      call tria1( rslt ,rr(3) ,mulocal2 )
    endif
  else
! 0 non-zero internal masses
    call tria0( rslt ,ss ,as ,mulocal2 )
  endif
! exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) 'c0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'c0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'c0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine


!*******************************************************************
! calculates
!
!    C   /                      d^(Dim)q
! ------ | --------------------------------------------------------
! i*pi^2 / [q^2-m1][(q+k1)^2-m2][(q+k1+k2)^2-m3][(q+k1+k2+k3)^2-m4]
!
! with  Dim = 4-2*eps
!         C = pi^eps * mu^(2*eps)
!             * GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
!
! input:  p1=k1^2, p2=k2^2, p3=k3^2, p4=(k1+k2+k3)^2, 
!         p12=(k1+k2)^2, p23=(k2+k3)^2, 
!         m1,m2,m3,m4=squared masses
! output: rslt(0) = eps^0   -coefficient
!         rslt(1) = eps^(-1)-coefficient
!         rslt(2) = eps^(-2)-coefficient
!
! Check the comments in  avh_olo_qp_onshell  to find out how this
! routines decides when to return IR-divergent cases.
!*******************************************************************

  subroutine d0cc( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
  use avh_olo_qp_box
  use avh_olo_qp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
!
  complex(kindr2) &   
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) ,0 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0ccr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
  use avh_olo_qp_box
  use avh_olo_qp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  complex(kindr2) &   
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  complex(kindr2) &   
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = areal(pp(ii))
    if (aimag(pp(ii)).ne.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'momentum with non-zero imaginary part, putting it to zero.'
      pp(ii) = acmplx( ap(ii) ,0 )
    endif
    ap(ii) = abs(ap(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rc( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
  use avh_olo_qp_box
  use avh_olo_qp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
!
  real(kindr2) &  
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rcr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
  use avh_olo_qp_box
  use avh_olo_qp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  complex(kindr2) &   
    ,intent(in)  :: m1,m2,m3,m4
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(6)
  complex(kindr2) &   
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = areal(mm(ii))
    hh = aimag(mm(ii))
    if (hh.gt.RZRO) then
      if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
        ,'mass-squared has positive imaginary part, switching its sign.'
      mm(ii) = acmplx( am(ii) ,-hh )
    endif
    am(ii) = abs(am(ii)) + abs(hh)
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 )
  use avh_olo_qp_box
  use avh_olo_qp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3,m4
!
  real(kindr2) &  
    :: pp(6)
  real(kindr2) &  
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = muscale 
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

  subroutine d0rrr( rslt ,p1,p2,p3,p4,p12,p23 ,m1,m2,m3,m4 ,rmu )
  use avh_olo_qp_box
  use avh_olo_qp_boxc
!
  complex(kindr2) &   
    ,intent(out) :: rslt(0:2)
  real(kindr2) &  
    ,intent(in)  :: p1,p2,p3,p4,p12,p23
  real(kindr2) &  
    ,intent(in)  :: m1,m2,m3,m4
  real(kindr2) &  
    ,intent(in)  :: rmu      
!
  real(kindr2) &  
    :: pp(6)
  real(kindr2) &  
    :: mm(4)
  complex(kindr2) &   
    :: ss(6),rr(4)
  real(kindr2) &  
    :: smax,ap(6),am(4),as(6),ar(4),s1r2,s2r2,s2r3,s3r4,s4r4
  real(kindr2) &  
    :: mulocal,mulocal2,small,hh,min13,min24,min56
  integer :: icase,ii,jj
  logical :: useboxc
  integer ,parameter :: lp(6,3)=&
           reshape((/1,2,3,4,5,6 ,5,2,6,4,1,3 ,1,6,3,5,4,2/),(/6,3/))
  integer ,parameter :: lm(4,3)=&
           reshape((/1,2,3,4     ,1,3,2,4     ,1,2,4,3    /),(/4,3/))
  character(25+99) ,parameter :: warning=&
                 'WARNING from OneLOop d0: '//warnonshell
  if (initz) call init
  pp(1) = p1
  pp(2) = p2
  pp(3) = p3
  pp(4) = p4
  pp(5) = p12
  pp(6) = p23
  mm(1) = m1
  mm(2) = m2
  mm(3) = m3
  mm(4) = m4
  smax = 0
!
  do ii=1,6
    ap(ii) = abs(pp(ii))
    if (ap(ii).gt.smax) smax = ap(ii)
  enddo
!
  do ii=1,4
    am(ii) = abs(mm(ii))
    if (am(ii).gt.smax) smax = am(ii)
  enddo
!
  small = 0
  do ii=1,6
    hh = abs(ap(ii))
    if (hh.gt.small) small=hh
  enddo
  small = small*neglig(prcpar)
!
  mulocal = rmu     
!
  mulocal2 = mulocal*mulocal
!
  if (smax.eq.RZRO) then
    if (eunit.gt.0) write(eunit,*) 'ERROR in OneLOop d0: ' &
      ,'all input equal zero, returning 0'
    rslt(0)=0 ;rslt(1)=0 ;rslt(2)=0
    return
  endif
!
  if (mulocal2.gt.smax) smax = mulocal2
!
  if (nonzerothrs) then
    hh = onshellthrs
    do ii=1,4
      if (ap(ii).lt.hh) ap(ii) = 0
      if (am(ii).lt.hh) am(ii) = 0
    enddo
  else
    hh = onshellthrs*smax
    if (wunit.gt.0) then
    do ii=1,4
      if (RZRO.lt.ap(ii).and.ap(ii).lt.hh) write(wunit,*) warning
      if (RZRO.lt.am(ii).and.am(ii).lt.hh) write(wunit,*) warning
    enddo
    endif
  endif
!
  jj = 1
  min56 = min(ap(5),ap(6))
  if (min56.lt.hh) then
    if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
      ,'input does not seem to represent hard kinematics, '&
      ,'trying to permutate'
    min13=min(ap(1),ap(3))
    min24=min(ap(2),ap(4))
    if     (min13.gt.min24.and.min13.gt.min56) then ;jj=2
    elseif (min24.gt.min13.and.min24.gt.min56) then ;jj=3
    else
      if (wunit.gt.0) write(wunit,*) 'WARNING from OneLOop d0: ' &
        ,'no permutation helps, errors might follow'
    endif
  endif
!
  icase = 0
  do ii=1,4
    if (am(lm(ii,jj)).gt.RZRO) icase = icase + base(ii)
  enddo
  ss(1)=pp(lp(permtable(1,icase),jj)) ;as(1)=ap(lp(permtable(1,icase),jj))
  ss(2)=pp(lp(permtable(2,icase),jj)) ;as(2)=ap(lp(permtable(2,icase),jj))
  ss(3)=pp(lp(permtable(3,icase),jj)) ;as(3)=ap(lp(permtable(3,icase),jj))
  ss(4)=pp(lp(permtable(4,icase),jj)) ;as(4)=ap(lp(permtable(4,icase),jj))
  ss(5)=pp(lp(permtable(5,icase),jj)) ;as(5)=ap(lp(permtable(5,icase),jj))
  ss(6)=pp(lp(permtable(6,icase),jj)) ;as(6)=ap(lp(permtable(6,icase),jj))
  rr(1)=mm(lm(permtable(1,icase),jj)) ;ar(1)=am(lm(permtable(1,icase),jj))
  rr(2)=mm(lm(permtable(2,icase),jj)) ;ar(2)=am(lm(permtable(2,icase),jj))
  rr(3)=mm(lm(permtable(3,icase),jj)) ;ar(3)=am(lm(permtable(3,icase),jj))
  rr(4)=mm(lm(permtable(4,icase),jj)) ;ar(4)=am(lm(permtable(4,icase),jj))
  icase = casetable(icase)
!
  s1r2 = abs(areal(ss(1)-rr(2))) + abs(aimag(ss(1)-rr(2)))
  s2r2 = abs(areal(ss(2)-rr(2))) + abs(aimag(ss(2)-rr(2)))
  s2r3 = abs(areal(ss(2)-rr(3))) + abs(aimag(ss(2)-rr(3)))
  s3r4 = abs(areal(ss(3)-rr(4))) + abs(aimag(ss(3)-rr(4)))
  s4r4 = abs(areal(ss(4)-rr(4))) + abs(aimag(ss(4)-rr(4)))
  if (nonzerothrs) then
    if (s1r2.lt.hh) s1r2 = 0
    if (s2r2.lt.hh) s2r2 = 0
    if (s2r3.lt.hh) s2r3 = 0
    if (s3r4.lt.hh) s3r4 = 0
    if (s4r4.lt.hh) s4r4 = 0
  elseif (wunit.gt.0) then
    if (RZRO.lt.s1r2.and.s1r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r2.and.s2r2.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s2r3.and.s2r3.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s3r4.and.s3r4.lt.hh) write(wunit,*) warning
    if (RZRO.lt.s4r4.and.s4r4.lt.hh) write(wunit,*) warning
  endif
!
  if     (icase.eq.4) then
!4 non-zero internal masses
    useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
               .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
               .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
               .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
               .or.(     areal(ss(1)).ge.-small  &
                    .and.areal(ss(2)).ge.-small  &
                    .and.areal(ss(3)).ge.-small  &
                    .and.areal(ss(4)).ge.-small) &
               .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
    if (useboxc) then
      call boxc( rslt ,ss,rr ,as ,smax )
    else
      call boxf4( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(1),rr(2),rr(3),rr(4) )
    endif
  elseif (icase.eq.3) then
!3 non-zero internal masses
    if (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      useboxc = (    (ar(1).ne.RZRO.and.aimag(rr(1)).ne.RZRO) &
                 .or.(ar(2).ne.RZRO.and.aimag(rr(2)).ne.RZRO) &
                 .or.(ar(3).ne.RZRO.and.aimag(rr(3)).ne.RZRO) &
                 .or.(ar(4).ne.RZRO.and.aimag(rr(4)).ne.RZRO) &
                 .or.(     areal(ss(1)).ge.-small  &
                      .and.areal(ss(2)).ge.-small  &
                      .and.areal(ss(3)).ge.-small  &
                      .and.areal(ss(4)).ge.-small) &
                 .or.(areal(ss(5)).ge.-small.and.areal(ss(6)).ge.-small))
      if (useboxc) then
        call boxc( rslt ,ss,rr ,as ,smax )
      else
        call boxf3( rslt, ss,rr )
      endif
    else
      call box16( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.5) then
!2 non-zero internal masses, opposite case
    if     (s1r2.ne.RZRO.or.s4r4.ne.RZRO) then
      if     (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
        call boxf5( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(2),rr(4) )
      else
        call box15( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
      endif
    elseif (s2r2.ne.RZRO.or.s3r4.ne.RZRO) then
      call box15( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    else
      call box14( rslt ,ss(5),ss(6) ,rr(2),rr(4) ,mulocal )
    endif
  elseif (icase.eq.2) then
!2 non-zero internal masses, adjacent case
    if     (as(1).ne.RZRO) then
      call boxf2( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) )
    elseif (s2r3.ne.RZRO) then
      if     (s4r4.ne.RZRO) then
        call box13( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
      else
        call box12( rslt ,ss(3),ss(2),ss(6),ss(5) ,rr(4),rr(3) ,mulocal )
      endif
    elseif (s4r4.ne.RZRO) then
      call box12( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    else
      call box11( rslt ,ss(3),ss(5),ss(6) ,rr(3),rr(4) ,mulocal )
    endif
  elseif (icase.eq.1) then
!1 non-zero internal mass
    if     (as(1).ne.RZRO) then
      if      (as(2).ne.RZRO) then
        call boxf1( rslt ,ss(1),ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) )
      else
        if     (s3r4.ne.RZRO) then
          call box10( rslt ,ss(1),ss(4),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box09( rslt ,ss(1),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      endif
    elseif (as(2).ne.RZRO) then
      if      (s4r4.ne.RZRO) then
        call box10( rslt ,ss(2),ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box09( rslt ,ss(2),ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    else
      if     (s3r4.ne.RZRO) then
        if     (s4r4.ne.RZRO) then
          call box08( rslt ,ss(3),ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
        else
          call box07( rslt ,ss(3),ss(5),ss(6) ,rr(4) ,mulocal )
        endif
      elseif (s4r4.ne.RZRO) then
        call box07( rslt ,ss(4),ss(5),ss(6) ,rr(4) ,mulocal )
      else
        call box06( rslt ,ss(5),ss(6) ,rr(4) ,mulocal )
      endif
    endif
  else
!0 non-zero internal mass
    call box00( rslt ,ss ,as ,mulocal )
  endif
!exp(eps*gamma_EULER) -> GAMMA(1-2*eps)/GAMMA(1-eps)^2/GAMMA(1+eps)
  rslt(0) = rslt(0) + 2*PISQo24*rslt(2)
!
  if (punit.gt.0) then
    if (nonzerothrs) write(punit,*) 'onshell:',trim(myprint(onshellthrs))
    write(punit,*) 'muscale:',trim(myprint(mulocal))
    write(punit,*) ' p1:',trim(myprint(p1))
    write(punit,*) ' p2:',trim(myprint(p2))
    write(punit,*) ' p3:',trim(myprint(p3))
    write(punit,*) ' p4:',trim(myprint(p4))
    write(punit,*) 'p12:',trim(myprint(p12))
    write(punit,*) 'p23:',trim(myprint(p23))
    write(punit,*) ' m1:',trim(myprint(m1))
    write(punit,*) ' m2:',trim(myprint(m2))
    write(punit,*) ' m3:',trim(myprint(m3))
    write(punit,*) ' m4:',trim(myprint(m4))
    write(punit,*) 'd0(2):',trim(myprint(rslt(2)))
    write(punit,*) 'd0(1):',trim(myprint(rslt(1)))
    write(punit,*) 'd0(0):',trim(myprint(rslt(0)))
  endif
  end subroutine

end module


module avh_olo

  use avh_olo_dp ,only: &
     olo_dp_kind=>olo_kind &
    ,olo_dp_scale=>olo_get_scale &
    ,olo_dp_onshell=>olo_get_onshell &
    ,olo_dp_precision=>olo_get_precision &
    ,olo,olo_a0,olo_an,olo_b0,olo_db0,olo_b11,olo_bn,olo_c0,olo_d0
  use avh_olo_qp ,only: &
     olo_qp_kind=>olo_kind &
    ,olo_qp_scale=>olo_get_scale &
    ,olo_qp_onshell=>olo_get_onshell &
    ,olo_qp_precision=>olo_get_precision &
    ,olo,olo_a0,olo_an,olo_b0,olo_db0,olo_b11,olo_bn,olo_c0,olo_d0

  implicit none

contains

  subroutine olo_unit( val ,message )
  use avh_olo_version
  use avh_olo_units ,only: set_unit
  integer     ,intent(in) :: val
  character(*),intent(in),optional :: message
  call olo_version
  if (present(message)) then ;call set_unit( message ,val )
  else                       ;call set_unit( 'all'   ,val )
  endif
  end subroutine

  subroutine olo_precision( ndec )
  use avh_olo_dp ,only: dp_sub=>olo_precision 
  use avh_olo_qp ,only: qp_sub=>olo_precision 
  integer ,intent(in) :: ndec
  call dp_sub( ndec ) 
  call qp_sub( ndec ) 
  end subroutine

  subroutine olo_scale( val )
  use avh_olo_dp ,only: dp_sub=>olo_scale 
  use avh_olo_qp ,only: qp_sub=>olo_scale 
  real(kind(1d0)) ,intent(in) :: val
  call dp_sub( val ) 
  call qp_sub( val ) 
  end subroutine

  subroutine olo_onshell( val )
  use avh_olo_dp ,only: dp_sub=>olo_onshell 
  use avh_olo_qp ,only: qp_sub=>olo_onshell 
  real(kind(1d0)) ,intent(in) :: val
  call dp_sub( val ) 
  call qp_sub( val ) 
  end subroutine

  subroutine olo_setting( iunit )
  use avh_olo_units
  use avh_olo_version
  integer,optional,intent(in) :: iunit
  integer :: nunit
  call olo_version
  nunit = munit
  if (present(iunit)) nunit = iunit
  if (nunit.le.0) return
  write(nunit,*) 'ERROR in OneLOop: subroutine olo_setting  is not available,'
  write(nunit,*) 'ERROR in OneLOop: use  function olo_get_scale  etc. instead.'
  end subroutine

end module

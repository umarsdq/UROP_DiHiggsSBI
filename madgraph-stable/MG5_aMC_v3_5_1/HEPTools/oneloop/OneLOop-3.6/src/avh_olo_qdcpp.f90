!!
!! Copyright (C) 2015 Andreas van Hameren. 
!!
!! This file is part of OneLOop-3.6.
!!
!! OneLOop-3.6 is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! OneLOop-3.6 is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with OneLOop-3.6.  If not, see <http://www.gnu.org/licenses/>.
!!


module avh_olo_prec
  use qdmodule

  implicit none
  public
  private :: IMAG,acmplx_r,acmplx_rr,acmplx_ir,acmplx_ri,acmplx_c
  private :: qcFROMi
  private :: prduct_qr_i,prduct_i_qr
  private :: prduct_qc_i,prduct_i_qc
  private :: ratio_qr_i,ratio_i_qr
  private :: ratio_qc_i,ratio_i_qc
  private :: plus_qr_i,plus_i_qr
  private :: plus_qc_i,plus_i_qc
  private :: minus_qr_i,minus_i_qr
  private :: minus_qc_i,minus_i_qc

  integer ,save :: prcpar=0
  integer ,save :: ndecim(1)
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
          ,save :: epsilo(1),neglig(1)

  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
    ,save :: RZRO ,RONE ,EPSN ,EPSN2 ,TWOPI ,ONEPI
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    ,save :: IEPS ,CZRO ,CONE ,IMAG ,PISQo24 ,IPI

  protected :: prcpar,ndecim,epsilo,neglig      !]PROTECTED
  protected :: RZRO,RONE,EPSN,EPSN2,TWOPI,ONEPI !]PROTECTED
  protected :: IEPS,CZRO,CONE,PISQo24,IPI       !]PROTECTED

  interface acmplx
    module procedure acmplx_r,acmplx_rr,acmplx_ir,acmplx_ri,acmplx_c
  end interface

  interface sqrt
    module procedure myCsqrt
  end interface

  interface operator (/)
    module procedure ratio_qc_i,ratio_i_qc
  end interface
  interface operator (+)
    module procedure plus_qc_i,plus_i_qc
  end interface
  interface operator (-)
    module procedure minus_qr_i,minus_i_qr
    module procedure minus_qc_i,minus_i_qc
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
  IMAG=cmplx(0d0,1d0,kind=kind(1d0))
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
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
    :: ten
  ten = 10                                       
  ndec = 32 !|RCTYPE=ddtype
!#  ndec = 64 !|RCTYPE=qdtype
  EPSN = ten**(-ndec)                            
  ndecim(prcpar) = ndec                         
  epsilo(prcpar) = EPSN                         
  neglig(prcpar) = EPSN*ten**(ndec/7)            
  end subroutine
!
  end subroutine


  function adble(xx) result(rslt)
!***********************************************************************
! Turn qd_real into kind(1d0)
!***********************************************************************
  type(dd_real) ,intent(in) :: xx !|RCTYPE=ddtype
!#  type(qd_real) ,intent(in) :: xx !|RCTYPE=qdtype
  real(kind(1d0)) :: rslt
  rslt = xx
  end function

  function convert(xx) result(rslt)
!***********************************************************************
! Turn kind(1d0) into qd_real
!***********************************************************************
  real(kind(1d0)) ,intent(in) :: xx
  type(dd_real) :: rslt !|RCTYPE=ddtype
!#  type(qd_real) :: rslt !|RCTYPE=qdtype
  rslt = ddreal(xx) !|RCTYPE=ddtype
!#  rslt = qdreal(xx) !|RCTYPE=qdtype
  end function

  function areal(zz) result(rslt)
!***********************************************************************
! Get real part of a complex
!***********************************************************************
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    ,intent(in) :: zz
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
    :: rslt
  rslt = ddreal(zz) !|RCTYPE=ddtype
!#  rslt = qdreal(zz) !|RCTYPE=qdtype
  end function

  function acmplx_r(xx) result(rslt)
!***********************************************************************
! Turn a real into a complex
!***********************************************************************
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
    ,intent(in) :: xx
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    :: rslt
  rslt = xx
  end function
  
  function acmplx_rr(xx,yy) result(rslt)
!***********************************************************************
! Turn two reals into one complex
!***********************************************************************
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
    ,intent(in) :: xx,yy
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    :: rslt
  rslt = xx + yy*IMAG
  end function
  
  function acmplx_ri(xx,yy) result(rslt)
!***********************************************************************
! Turn a real and an integer into one complex
!***********************************************************************
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
          ,intent(in) :: xx
  integer ,intent(in) :: yy
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    :: rslt
  rslt = xx + ddreal(yy)*IMAG !|RCTYPE=ddtype
!#  rslt = xx + qdreal(yy)*IMAG !|RCTYPE=qdtype
  end function
  
  function acmplx_ir(xx,yy) result(rslt)
!***********************************************************************
! Turn an integer and a real into one complex
!***********************************************************************
  integer ,intent(in) :: xx
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
          ,intent(in) :: yy
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    :: rslt
  rslt = ddreal(xx) + yy*IMAG !|RCTYPE=ddtype
!#  rslt = qdreal(xx) + yy*IMAG !|RCTYPE=qdtype
  end function
  
  function acmplx_c(zz) result(rslt)
!***********************************************************************
! Replaces the real part of zz by its absolute value
!***********************************************************************
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    ,intent(in) :: zz
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    :: rslt
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
    :: xx,yy
  xx = ddreal(zz) !|RCTYPE=ddtype
!#  xx = qdreal(zz) !|RCTYPE=qdtype
  xx = abs(xx)
  yy = aimag(zz)
  rslt = xx + yy*IMAG
  end function


  function myCsqrt(zz) result(rslt)
!***********************************************************************
!* The complex square-root
!***********************************************************************
  intent(in) :: zz
  type(dd_complex) & !|RCTYPE=ddtype
!#  type(qd_complex) & !|RCTYPE=qdtype
    :: zz,rslt
  type(dd_real) & !|RCTYPE=ddtype
!#  type(qd_real) & !|RCTYPE=qdtype
    :: aa,bb,cc,xx,yy,y2
  xx = areal(zz)
  yy = aimag(zz)
  if (yy.eq.RZRO) then
    if     (xx.lt.RZRO) then ;rslt=sqrt(-xx)*IMAG
    elseif (xx.gt.RZRO) then ;rslt=sqrt( xx)
                        else ;rslt=0
    endif
    return
  endif
  y2 = yy*yy
  aa = sqrt(xx*xx+y2)
  if (xx.ge.RZRO) then ;bb=aa+xx
                  else ;bb=y2/(aa-xx)
  endif
  cc = sqrt(bb/2)
  rslt = cc + IMAG*yy/(cc*2)
  end function
 
  
  include 'avh_olo_intrf.h90'  

end module

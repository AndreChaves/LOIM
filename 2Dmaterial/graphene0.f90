module graphene_susc
implicit none

integer, PARAMETER :: LONG_INT = SELECTED_INT_KIND(9)
integer, PARAMETER :: LATT_INT = SELECTED_INT_KIND(3) 
integer, PARAMETER :: SMAL_INT = SELECTED_INT_KIND(1) 
!   integer, parameter :: dp = selected_real_kind(14) 

integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
real(fkind),parameter ::  fermi_velocity=1.*10**6, &
                          pi=acos(-1.0),planck_constant=6.62*10.**(-34), &
                          electric_charge=1.6*10.**(-19),hbar=6.58211928*10.**(-16) ,& !hbar(eVs) 
                          magnetic_permeability=4.*pi*10.**(-7)
 complex(fkind),parameter :: I=(0.,1.) 

contains


function func1(u) result(z)
  integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
  real(fkind),parameter :: pi=acos(-1.0)
  complex(fkind),parameter :: I=(0.,1.) 
  complex(fkind) :: u,z
  
  z=u*sqrt(u**2-1.)-acosh(u)
  if(aimag(u) .ne. 0) then
     z=z-I*pi/2.
  endif
end function func1

function func2(u) result(z)
  integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
  complex(fkind),parameter :: I=(0.,1.) 
  complex(fkind) :: u,z
  z=u*sqrt(1.-u**2)-acos(u)
end function func2

function polarizibility(y,x) result(chi)
  integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
  complex(fkind),parameter :: I=(0.,1.) 
  complex(fkind) ::  y,chi
  real(fkind) :: y2,x
  y2=real(y)
  if(x>y2) then
     if(x+y2<2.) then
         chi =polarizibility1A(y,x)
     else
         chi = polarizibility2A(y,x)
      endif
  else 
       if(x+y2<2.) then
         chi =polarizibility1B(y,x)
       else
         if(y2>x+2.) then
             chi = polarizibility3B(y,x)
         else
             chi = polarizibility2B(y,x)
         endif
       endif
  endif
end function

function polarizibility3B(y,x) result(chi)
 integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
 complex(fkind),parameter :: I=(0.,1.) 
 complex(fkind) ::y,z1,z2,chi
 real(fkind) :: x

 z1=(2.+y)/x
 z2=(y-2.)/x
 chi=-2./pi+1./(4.*pi)*x**2/sqrt(y**2-x**2)*(func1(z1)-func1(z2))-I/4./sqrt(y**2-x**2)*x**2
!  chi=-0.25*I*x**2/sqrt(y**2-x**2)
end function 

function polarizibility3A(y,x) result(chi)
 integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
 complex(fkind),parameter :: I=(0.,1.) 
 real(fkind),parameter :: pi=acos(-1.0)
 complex(fkind) ::y,z1,z2,chi
 real(fkind) ::x

 z1=(2.+y)/x
 z2=(y-2.)/x
 chi=-2./pi+1./(4.*pi)*x**2/sqrt(x**2-y**2)*(func2(z1)-func2(z2))
end function 

function polarizibility1A(y,x) result(chi)
  integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
  complex(fkind),parameter :: I=(0.,1.)
  real(fkind),parameter :: pi=acos(-1.0)  
  complex(fkind) :: y,z1,z2,chi
  real(fkind) :: x
  
  z1=(2.-y)/x
  z2=(2.+y)/x
  chi=-2./pi+I/(4.*pi)*x**2/sqrt(x**2-y**2)*(func1(z1)-func1(z2))
end function

function polarizibility2A(y,x) result(chi)
  integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
  real(fkind),parameter :: pi=acos(-1.0)
  complex(fkind),parameter :: I=(0.,1.) 
  complex(fkind) :: y,z1,z2,chi
  real(fkind) :: x
  
  z1=(2.-y)/x
  z2=(2.+y)/x
  chi=-2./pi+1./(4.*pi)*x**2/sqrt(x**2-y**2)*(func2(z1)-I*func1(z2))
end function

function polarizibility1B(y,x) result(chi)
  integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
  complex(fkind),parameter :: I=(0.,1.) 
  real(fkind),parameter :: pi=acos(-1.0)
  complex(fkind) :: y,z1,z2,chi
  real(fkind) :: x
  
  z1=(2.-y)/x
  z2=(2.+y)/x
  chi=-2./pi+1./(4.*pi)*x**2/sqrt(y**2-x**2)*(func1(z2)-func1(z1))
end function

function polarizibility2B(y,x) result(chi)
  integer, PARAMETER :: fkind = SELECTED_REAL_KIND(15,307)
  complex(fkind),parameter :: I=(0.,1.) 
  real(fkind),parameter :: pi=acos(-1.0)
  complex(fkind) :: y,z1,z2,chi
  real(fkind) :: x
  
  z1=(2.-y)/x
  z2=(2.+y)/x
  chi=-2./pi+1./(4.*pi)*x**2/sqrt(y**2-x**2)*(func1(z2)+I*func2(z1))
end function

end module

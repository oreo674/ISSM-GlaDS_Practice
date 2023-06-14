!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Olivier Gaglairdini
! *  Email:   
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *   2016/09/02 
! *****************************************************************************
!> USF_ValleySHMIOP.F90
!>
!> Compute Zs and Zb for the valley geometry of SHMIP tests

!> Needs the following parameters in the Constants section :
!> Valley Parameter : The controling parameter for the overdeepening 
         
FUNCTION GetZb ( Model, nodenumber, x) RESULT(Zb)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(ValueList_t), POINTER :: Constants
   INTEGER :: nodenumber
   REAL(KIND=dp) :: x, y, bed, Para, Zb         
   LOGICAL :: GotIt        
   
   x = Model % Nodes % x (nodenumber)
   y = Model % Nodes % y (nodenumber)
   
   ! Get the constants
   Constants => GetConstants()
   Para = GetConstReal( Constants, 'Valley Parameter', GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(a)')'Keyword >Valley Parameter< not found in Constant section'
      CALL FATAL('GetZb', Message)
   END IF
   
   Zb = bed(x,y,Para)

END FUNCTION GetZb

FUNCTION GetZs ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(ValueList_t), POINTER :: Constants
   INTEGER :: nodenumber
   REAL(KIND=dp) :: x, y, surface, bed, para, Zs         
   LOGICAL :: GotIt        
   
   x = Model % Nodes % x (nodenumber)
   y = Model % Nodes % y (nodenumber)

   ! Get the constants
   Constants => GetConstants()
   Para = GetConstReal( Constants, 'Valley Parameter', GotIt )
   IF (.NOT.GotIt) THEN
      WRITE(Message,'(a)')'Keyword >Valley Parameter< not found in Constant section'
      CALL FATAL('GetZb', Message)
   END IF
   
   Zs = bed(x,y,Para) + MAX(1.0,surface(x,y)-bed(x,y,Para))

END FUNCTION GetZs

FUNCTION surface(x,y)
   USE types
   IMPLICIT NONE
   REAL(KIND=dp) :: x, y, surface
   surface = 100.0*(x+200.0)**0.25 + x/60.0 - 2.0e10**0.25 + 1.0
END FUNCTION surface   

FUNCTION f(x,para)
   USE types
   IMPLICIT NONE
   REAL(KIND=dp) :: x, para, f, surface 
   f = (surface(6.0e3_dp,0.0_dp) - para*6.0e3)/6.0e3**2.0 * x**2.0 + para*x
END FUNCTION f         

FUNCTION g(y)
   USE types
   IMPLICIT NONE
   REAL(KIND=dp) :: y, g
   g = 0.5e-6 * abs(y)**3.0 
END FUNCTION g         

FUNCTION h(x, para)
   USE types
   IMPLICIT NONE
   REAL(KIND=dp) :: x, para, h, surface, f
   REAL(KIND=dp), PARAMETER :: para_bench = 0.05
   h = (-4.5*x/6.0e3 + 5.0) * (surface(x,0.0_dp) - f(x,para)) &
    / (surface(x,0.0_dp) - f(x, para_bench) + 1.0e-10 ) 
END FUNCTION h         

FUNCTION bed(x, y, para)
   USE types
   IMPLICIT NONE
   REAL(KIND=dp) :: x, y, para, bed, f, g, h
   bed = f(x, para) + g(y) * h(x, para)
END FUNCTION bed       


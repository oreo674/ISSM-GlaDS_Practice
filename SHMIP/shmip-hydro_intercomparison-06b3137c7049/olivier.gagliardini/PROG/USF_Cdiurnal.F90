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
! *   2016/10/24 
! *****************************************************************************
!> USF_Cdiurnal.F90
!>
!> Compute runoff input in moulins for the SHMIP C tests

!> Needs the following parameters in the Constants section :
!> ra parameter : The controling parameter for the amplitude of runoff 
         
FUNCTION GetMeltInput ( Model, nodenumber, t) RESULT(melt)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(ValueList_t), POINTER :: Constants
   TYPE(Variable_t), POINTER :: TimeVar
   INTEGER :: nodenumber
   REAL(KIND=dp) :: t, melt, ra         
   REAL(KIND=dp), Parameter :: dayinsec = 24*60*60
   REAL(KIND=dp), Parameter :: yearinsec = 365*dayinsec
   REAL(KIND=dp), Parameter :: moulinIn = 0.9*yearinsec
   
   ! Get the constants
   Constants => GetConstants()
   ra = ListGetConstReal( Constants, 'ra Parameter', UnFoundFatal=.TRUE. )
   
   TimeVar => VariableGet( Model % Mesh % Variables, 'Time' )
   t = TimeVar % Values(1) * yearinsec 

   melt = MAX(0.0_dp, moulinIn*(1.0-ra*SIN(2.0*Pi*t/dayinsec))) 
END FUNCTION GetMeltInput



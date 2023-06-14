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
!> USF_Dseasonal.F90
!>
!> Compute basal melt for the SHMIP D tests

!> Needs the following parameters in the Constants section :
!> DT parameter : The delta temperature parameter for the amplitude of runoff 
         
FUNCTION GetMeltInput ( Model, nodenumber, input) RESULT(melt)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(ValueList_t), POINTER :: Constants
   INTEGER :: nodenumber
   REAL(KIND=dp) :: input(2) 
   REAL(KIND=dp) :: t, melt, DT, zs, temp     
   REAL(KIND=dp), Parameter :: dayinsec = 24*60*60
   REAL(KIND=dp), Parameter :: yearinsec = 365*dayinsec
   REAL(KIND=dp), Parameter :: lr = -0.0075_dp 
   REAL(KIND=dp), Parameter :: DDF = 0.01*365 !in m/K/year
   REAL(KIND=dp), Parameter :: basal = 7.93e-11*yearinsec 
   
   ! Get the constants
   Constants => GetConstants()
   DT = ListGetConstReal( Constants, 'DT Parameter', UnFoundFatal=.TRUE. )
   
   zs = input(1) 
   t = input(2)
   
   temp = -16.0*COS(2.0*Pi*t) - 5.0 + DT  
   melt = MAX(0.0_dp, (zs*lr+temp)*DDF) + basal 
END FUNCTION GetMeltInput


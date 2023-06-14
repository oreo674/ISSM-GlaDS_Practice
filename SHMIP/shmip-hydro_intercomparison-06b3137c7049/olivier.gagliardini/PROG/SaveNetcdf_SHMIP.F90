!------------------------------------------------------------------------------
! A solver to save in NETCDF the SHMIP experiments
! From an example by Fabien, OG - August 2016
!------------------------------------------------------------------------------
      SUBROUTINE SAVENETCDF( Model,Solver,dt,TransientSimulation )
      USE DefUtils
      USE MainUtils
      USE Netcdf
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t),TARGET :: Solver  !< Linear & nonlinear equation solver options
      TYPE(Model_t) :: Model    !< All model information (mesh, materials, BCs, etc...)
      REAL(KIND=dp) :: dt       !< Timestep size for time dependent simulations
      LOGICAL :: TransientSimulation      !< Steady state or transient simulation
!------------------------------------------------------------------------------
      
      TYPE(ValueList_t), POINTER :: Params, Constants
      TYPE(Element_t),POINTER :: Element, Edge
      TYPE(Variable_t), POINTER :: TimeVar

      ! Output Elmer Variables
      CHARACTER(LEN=MAX_NAME_LEN) :: SheetThicknessName, ChannelAreaName, ZbName
      TYPE(Variable_t), POINTER :: ThickSol, AreaSol, NSol, ZbSol,  &
          qSol, ZsSol, QcSol, hstoreSol
      REAL(KIND=dp), POINTER :: ThickSolution(:), ZsSolution(:), &
           NSolution(:), AreaSolution(:), ZbSolution(:), qSolution(:), &
           QcSolution(:), hstoreSolution(:)
      INTEGER, POINTER :: ZbPerm(:), ZsPerm(:), QcPerm(:), &
             ThickPerm(:), NPerm(:), AreaPerm(:), qPerm(:), &
             hstorePerm(:)
   
      LOGICAL :: GotIt, FirstTime = .TRUE., Found

      ! Define the required Dimensions 
      INTEGER :: D_nodenr                ! number of nodes
      INTEGER :: D_cellnr                ! number of cells (elements)
      INTEGER :: D_edgenr             ! number of channels
      INTEGER, PARAMETER :: D_dim = 2               ! dimension 
      INTEGER, PARAMETER :: D_nodespercell = 3      ! number of nodes per cell
      INTEGER, PARAMETER :: D_nodesperedge = 2   ! number of nodes to define a channel
      INTEGER :: D_time                             ! number of output time               

      ! Define the required Variables
      INTEGER, ALLOCATABLE :: V_nodenr(:)        
      INTEGER, ALLOCATABLE :: V_cellnr(:) 
      INTEGER, ALLOCATABLE :: V_edgenr(:) 
      INTEGER, DIMENSION(D_dim) :: V_dim 
      INTEGER, DIMENSION(D_nodespercell) :: V_nodespercell
      INTEGER, DIMENSION(D_nodesperedge) :: V_nodesperedge
      REAL(KIND=DP) :: V_time
      REAL(KIND=DP), ALLOCATABLE :: V_xy(:,:)        ! Nodes coordinates 
      REAL(KIND=DP), ALLOCATABLE :: V_index2(:,:)    ! Cell coordinates
      REAL(KIND=DP), ALLOCATABLE :: V_index_ch(:,:)  ! Midpoint channel coodinates
      INTEGER, ALLOCATABLE :: V_cellconnect(:,:)     ! Cell connectivity  
      INTEGER, ALLOCATABLE :: V_channelconnect(:,:)  ! Channel connectivity  
      REAL(KIND=DP), ALLOCATABLE :: V_H(:)           ! Ice Thickness  
      REAL(KIND=DP), ALLOCATABLE :: V_B(:)           ! Bed elevation  
      REAL(KIND=DP), ALLOCATABLE :: V_N(:)           ! Effective Pressure
      REAL(KIND=DP), ALLOCATABLE :: V_hs(:)          ! layer thickness   
      REAL(KIND=DP), ALLOCATABLE :: V_he(:)      ! Storage layer thickness   
      REAL(KIND=DP), ALLOCATABLE :: V_qs(:)          ! discharge distributed system
      REAL(KIND=DP), ALLOCATABLE :: V_S(:)           ! Channel cross-sectional area
      REAL(KIND=DP), ALLOCATABLE :: V_Q(:)           ! Discharge channel           

      REAL(KIND=DP) :: qx, qy, StartTime
      INTEGER :: t,i,n, Ntime
      REAL(KIND=DP), PARAMETER :: Yearinsec = 365.0_dp*24.0_dp*3600.0_dp

      INTEGER :: ncid     !netcdf file ID
      INTEGER :: dimid(6) ! Number of Dimensions in the netcdf
      INTEGER :: varid(13) ! Number of Variables in the netcdf

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='SAVE_NETCDF'
      CHARACTER(LEN=MAX_NAME_LEN) :: fname ! Name of the netcdf file

      ! Save what is required                                  
      SAVE D_nodenr, D_cellnr, D_edgenr
      SAVE V_nodenr, V_cellnr, V_edgenr, V_time, V_xy, &
         V_index2, V_index_ch, &
         V_cellconnect, V_channelconnect, V_H, V_B, V_N, V_hs, &
         V_he, V_qs, V_S, V_Q
      SAVE FirstTime, Ntime, D_time, ncid, StartTime
      SAVE SheetThicknessName, ChannelAreaName, ZbName, fname
      SAVE varid, dimid
     
      ! Will be done only the first time 
      IF (FirstTime) THEN
         D_nodenr = Solver % Mesh % NumberOfNodes
         D_cellnr = Solver % NumberOfActiveElements
         D_edgenr = Solver % Mesh % NumberOfEdges 

         ALLOCATE( V_nodenr(D_nodenr), & 
                   V_cellnr(D_cellnr), &
                   V_edgenr(D_edgenr), &
                   V_xy(D_nodenr,D_dim), &
                   V_index2(D_cellnr,D_dim), &
                   V_index_ch(D_edgenr,D_dim), &
                   V_cellconnect(D_cellnr,D_nodespercell), &
                   V_channelconnect(D_edgenr,D_nodesperedge), &
                   V_H(D_nodenr), V_B(D_nodenr), V_N(D_nodenr), &
                   V_hs(D_nodenr), V_he(D_nodenr), V_qs(D_cellnr), &
                   V_S(D_edgenr), V_Q(D_edgenr) )

         ! Read parameters in the sif file 
         Params => Solver % Values
         fname=ListGetString(Params, 'OutputFile Name', GotIt)
         IF (.Not.GotIt) &
            CALL FATAL(SolverName,'<OutputFile Name> Not Found')
         D_time = GetInteger(Params, 'Number of Outputs', Gotit)
         IF (.NOT.GotIt) & 
            CALL FATAL(SolverName,'<Number of Outputs> Not Found')
         StartTime = GetConstReal(Params, 'First Time Saved', Gotit)
         If (.Not.GotIt) StartTime = 0.0_dp

         ! Read some more constant  
         Constants => GetConstants()
         ChannelAreaName = GetString( Constants, &
                   'Channel Area Variable Name', Found )
         IF(.NOT.Found) THEN        
            CALL WARN(SolverName,'Keyword >Channel Area Variable Name<  &
                        & not found in section Constants')
            CALL WARN(SolverName,'Taking default value >Channel Area<')
            WRITE(ChannelAreaName,'(A)') 'Channel Area'
         END IF

         SheetThicknessName = GetString( Constants, &
               'Sheet Thickness Variable Name', Found )
         IF(.NOT.Found) THEN        
            CALL WARN(SolverName,'Keyword >Sheet Thickness Variable Name< &
                       & not found in section Constants')
            CALL WARN(SolverName,'Taking default value >Sheet Thickness<')
            WRITE(SheetThicknessName,'(A)') 'Sheet Thickness'
         END IF

         ZbName = GetString( Constants, 'Bedrock Variable Name', Found )
         IF(.NOT.Found) THEN        
            CALL WARN(SolverName,'Keyword >Bedrock Variable Name< & 
                       & not found in section Constants')
            CALL WARN(SolverName,'Taking default value >zb<')
            WRITE(ZbName,'(A)') 'Zb'
         END IF

         ZsSol => VariableGet( Solver % Mesh % Variables, 'Zs', & 
                                                UnfoundFatal = .TRUE. )
         ZsPerm     => ZsSol % Perm
         ZsSolution => ZsSol % Values

         ZbSol => VariableGet( Solver % Mesh % Variables, ZbName, &
                                                UnfoundFatal = .TRUE. )
         ZbPerm     => ZbSol % Perm
         ZbSolution => ZbSol % Values

         ! Create the netcf 
         CALL CreateNetcdfFile

         ! Fill in what is done only at first time
         ! Node number and coordinates
         DO i = 1, D_nodenr
            V_nodenr(i)=i-1
            V_xy(i,1) = Solver % Mesh % Nodes % x(i)
            V_xy(i,2) = Solver % Mesh % Nodes % y(i)
         END DO

         ! Element connectivity
         DO t = 1, D_cellnr
            V_cellnr(t)=t-1
            Element =>  GetActiveElement(t)
            n = GetElementNOFNodes()
            IF (n/=3) CALL FATAL(SolverName,&
                           'Sorry solver only support 3 nodes elements')
            V_cellconnect(t,1:n) = Element % NodeIndexes(1:n)-1
            V_index2(t,1) = SUM(Solver%Mesh%Nodes%x(Element%NodeIndexes(1:n)))/3.0
            V_index2(t,2) = SUM(Solver%Mesh%Nodes%y(Element%NodeIndexes(1:n)))/3.0
         END DO

         ! Edge connectivity
         DO t = 1, D_edgenr
            V_edgenr(t)=t-1
            Edge => Solver % Mesh % Edges(t)
            n = Edge % TYPE % NumberOfNodes
            IF (n/=2) CALL FATAL(SolverName, & 
                 'Work only for edge element of type 202')
            V_channelconnect(t,1:n) = Edge % NodeIndexes(1:n)-1     
            V_index_ch(t,1) = SUM(Solver%Mesh%Nodes%x(Edge%NodeIndexes(1:n)))/2.0
            V_index_ch(t,2) = SUM(Solver%Mesh%Nodes%y(Edge%NodeIndexes(1:n)))/2.0
         END DO

         ! Thickness elevation in [m]
         V_H(1:D_nodenr) = ZsSolution(ZsPerm(1:D_nodenr))-ZbSolution(ZbPerm(1:D_nodenr))

         ! Bedrock elevation in [m]
         V_B(1:D_nodenr) = ZbSolution(ZbPerm(1:D_nodenr))

         V_dim = (/0,1/)
         V_nodespercell = (/0,1,2/)
         V_nodesperedge = (/0,1/)

         ! Save time independant variable in netcdf file 
         CALL handle_err(nf90_put_var(ncid,varid(2),V_xy))
         CALL handle_err(nf90_put_var(ncid,varid(3),V_index2))
         CALL handle_err(nf90_put_var(ncid,varid(4),V_index_ch))
         CALL handle_err(nf90_put_var(ncid,varid(5),V_channelconnect))
         CALL handle_err(nf90_put_var(ncid,varid(6),V_H ))
         CALL handle_err(nf90_put_var(ncid,varid(7),V_B ))

         ! Close  Netcdf file
         CALL handle_err(NF90_CLOSE(ncid))

         Ntime = 1
         FirstTime = .FALSE. 
      END IF


      ! Save the current time [s] 
      TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )

      ! Start saving when StartTime is reached
      IF ((TimeVar % Values(1) - StartTime)/StartTime < -1.0e-6 ) RETURN

      V_time = TimeVar % Values(1) * Yearinsec  
      write(*,*)'V_time', V_time
      write(*,*)'Ntime', Ntime

      ! Will now save the time-dependant variable for that time step
      NSol => VariableGet( Solver % Mesh % Variables, 'Effective Pressure', & 
                                                UnfoundFatal = .TRUE. )

      NPerm     => NSol % Perm
      NSolution => NSol % Values

      ThickSol => VariableGet( Solver % Mesh % Variables,SheetThicknessName, & 
                                                UnfoundFatal = .TRUE. )
      ThickPerm     => ThickSol % Perm
      ThickSolution => ThickSol % Values

      qSol => VariableGet( Solver % Mesh % Variables, 'Sheet Discharge', &
                                                UnfoundFatal = .TRUE. )
      qPerm     => qSol % Perm
      qSolution => qSol % Values

      hstoreSol => VariableGet( Solver % Mesh % Variables, & 
          'Sheet Storage', UnfoundFatal = .TRUE. )
      hstorePerm     => hstoreSol % Perm
      hstoreSolution => hstoreSol % Values

      AreaSol => VariableGet( Solver % Mesh % Variables, ChannelAreaName, & 
                                                UnfoundFatal = .TRUE. )
      AreaPerm     => AreaSol % Perm
      AreaSolution => AreaSol % Values

      QcSol => VariableGet( Solver % Mesh % Variables, 'Channel Flux', & 
                                                UnfoundFatal = .TRUE. )
      QcPerm     => QcSol % Perm
      QcSolution => QcSol % Values

      ! Effective presssure in [Pa]
      V_N(1:D_nodenr) = NSolution(NPerm(1:D_nodenr))*1.0e6
      ! Sheet thickness in [m]
      V_hs(1:D_nodenr) = ThickSolution(ThickPerm(1:D_nodenr))
      ! sheet storgae in [m]
      V_he(1:D_nodenr) = hstoreSolution(hstorePerm(1:D_nodenr))
      ! Channel flux in [m3/s]
      V_Q(1:D_edgenr) = QcSolution(QcPerm(D_nodenr+1:D_nodenr+D_edgenr)) / Yearinsec
      ! Channel cross-sectional area in [m]
      V_S(1:D_edgenr) = AreaSolution(AreaPerm(D_nodenr+1:D_nodenr+D_edgenr))
      ! sheet discharge in [m2/s]
      ! calculate the norm + the value in the cell (=mean nodal values)
      DO t = 1, D_cellnr
         Element =>  GetActiveElement(t)
         n = GetElementNOFNodes()
         qx = 0.0_dp
         qy = 0.0_dp
         DO i = 1, n
            qx = qx + qSolution(D_dim*qPerm(Element % NodeIndexes(i))-1)
            qy = qy + qSolution(D_dim*qPerm(Element % NodeIndexes(i)))
         END DO
         qx = qx / n
         qy = qy / n
         V_qs(t) = SQRT(qx*qx+qy*qy)/ Yearinsec
      END DO
      ! If a nodal value
      !V_qs(1:D_nodenr) = SQRT(qSolution(D_dim*qPerm(1:D_nodenr)-1)**2.0 + &
      !       qSolution(D_dim*qPerm(1:D_nodenr))**2.0 ) / Yearinsec 

      ! Open  Netcdf file
      CALL handle_err(nf90_open(TRIM(fname),NF90_WRITE,ncid)) 

      ! Save all time dependent variables for that time
      CALL handle_err(nf90_put_var(ncid,varid(1),V_time, start = (/Ntime/)))
      CALL handle_err(nf90_put_var(ncid,varid(8),V_N, start = (/1,Ntime/)))
      CALL handle_err(nf90_put_var(ncid,varid(9),V_hs, start = (/1,Ntime/)))
      CALL handle_err(nf90_put_var(ncid,varid(10),V_he, start = (/1,Ntime/)))
      CALL handle_err(nf90_put_var(ncid,varid(11),V_qs, start = (/1,Ntime/)))
      CALL handle_err(nf90_put_var(ncid,varid(12),V_S, start = (/1,Ntime/)))
      CALL handle_err(nf90_put_var(ncid,varid(13),V_Q, start = (/1,Ntime/)))

      ! Close  Netcdf file
      CALL handle_err(NF90_CLOSE(ncid))

      ! Update Ntime
      Ntime = Ntime + 1

      CONTAINS

      SUBROUTINE CreateNetcdfFile
      ! Create file 
         CALL handle_err(nf90_create(TRIM(fname),NF90_CLOBBER,ncid)) 
      ! Define Dimension
         CALL handle_err(nf90_def_dim(ncid,'time',D_time,dimid(1)))
         CALL handle_err(nf90_def_dim(ncid,'dim',D_dim,dimid(2)))
         CALL handle_err(nf90_def_dim(ncid,'n_nodes_ch',D_nodesperedge,dimid(3)))
         CALL handle_err(nf90_def_dim(ncid,'index1',D_nodenr,dimid(4)))
         CALL handle_err(nf90_def_dim(ncid,'index2',D_cellnr,dimid(5)))
         CALL handle_err(nf90_def_dim(ncid,'index_ch',D_edgenr,dimid(6)))

      ! Define Variables        
         CALL handle_err(nf90_def_var(ncid,'time',NF90_DOUBLE,dimid(1),varid(1)))
         CALL handle_err(nf90_put_att(ncid,varid(1),'units','s'))
         CALL handle_err(nf90_put_att(ncid,varid(1),'long_name','time'))

         CALL handle_err(nf90_def_var(ncid,'coords1',NF90_DOUBLE,(/dimid(4),dimid(2)/),varid(2)))
         CALL handle_err(nf90_put_att(ncid,varid(2),'units','m'))
         CALL handle_err(nf90_put_att(ncid,varid(2),'long_name','node coordinates'))

         CALL handle_err(nf90_def_var(ncid,'coords2',NF90_DOUBLE,(/dimid(5),dimid(2)/),varid(3)))
         CALL handle_err(nf90_put_att(ncid,varid(3),'units','m'))
         CALL handle_err(nf90_put_att(ncid,varid(3),'long_name','cell midpoint coordinates'))

         CALL handle_err(nf90_def_var(ncid,'coords_ch',NF90_DOUBLE,(/dimid(6),dimid(2)/),varid(4)))
         CALL handle_err(nf90_put_att(ncid,varid(4),'units','m'))
         CALL handle_err(nf90_put_att(ncid,varid(4),'long_name', &
                 &   'channel midpoint coordinates'))

         CALL handle_err(nf90_def_var(ncid,'connect_ch',NF90_INT,(/dimid(6),dimid(3)/),varid(5)))
         CALL handle_err(nf90_put_att(ncid,varid(5),'units',''))
         CALL handle_err(nf90_put_att(ncid,varid(5),'long_name', &
                                &   'channel connectivity'))

         CALL handle_err(nf90_def_var(ncid,'H',NF90_DOUBLE,dimid(4),varid(6)))
         CALL handle_err(nf90_put_att(ncid,varid(6),'units','m'))
         CALL handle_err(nf90_put_att(ncid,varid(6), &
                                 'long_name','ice thickness'))

         CALL handle_err(nf90_def_var(ncid,'B',NF90_DOUBLE,dimid(4),varid(7)))
         CALL handle_err(nf90_put_att(ncid,varid(7),'units','m'))
         CALL handle_err(nf90_put_att(ncid,varid(7), &
                                 'long_name','bed elevation'))

         CALL handle_err(nf90_def_var(ncid,'N',NF90_DOUBLE,(/dimid(4),dimid(1)/),varid(8)))
         CALL handle_err(nf90_put_att(ncid,varid(8),'units','Pa'))
         CALL handle_err(nf90_put_att(ncid,varid(8), &
                                 'long_name','effective pressure'))

         CALL handle_err(nf90_def_var(ncid,'h',NF90_DOUBLE,(/dimid(4),dimid(1)/),varid(9)))
         CALL handle_err(nf90_put_att(ncid,varid(9),'units','m'))
         CALL handle_err(nf90_put_att(ncid,varid(9), &
                                 'long_name','water sheet thickness'))

         CALL handle_err(nf90_def_var(ncid,'hstore',NF90_DOUBLE,(/dimid(4),dimid(1)/),varid(10)))
         CALL handle_err(nf90_put_att(ncid,varid(10),'units','m'))
         CALL handle_err(nf90_put_att(ncid,varid(10), &
                  'long_name','stored water effective layer thickness'))

         CALL handle_err(nf90_def_var(ncid,'q',NF90_DOUBLE,(/dimid(5),dimid(1)/),varid(11)))
         CALL handle_err(nf90_put_att(ncid,varid(11),'units','m^2/s'))
         CALL handle_err(nf90_put_att(ncid,varid(11), &
                          'long_name','water sheet discharge'))

         CALL handle_err(nf90_def_var(ncid,'S',NF90_DOUBLE,(/dimid(6),dimid(1)/),varid(12)))
         CALL handle_err(nf90_put_att(ncid,varid(12),'units','m^2'))
         CALL handle_err(nf90_put_att(ncid,varid(12), &
                       'long_name','channel cross-sectional area'))
                                  
         CALL handle_err(nf90_def_var(ncid,'Q',NF90_DOUBLE,(/dimid(6),dimid(1)/),varid(13)))
         CALL handle_err(nf90_put_att(ncid,varid(13),'units','m^3/s'))
         CALL handle_err(nf90_put_att(ncid,varid(13), &
                                    'long_name','channel discharge'))

         ! Global attribute
         CALL handle_err(nf90_put_att(ncid,NF90_GLOBAL,'title',fname))
         CALL handle_err(nf90_put_att(ncid,NF90_GLOBAL,'meshtype','unstructured'))
         CALL handle_err(nf90_put_att(ncid,NF90_GLOBAL,'dimension','2D'))
         CALL handle_err(nf90_put_att(ncid,NF90_GLOBAL,'channels_on_edges','yes'))
         CALL handle_err(nf90_put_att(ncid,NF90_GLOBAL,'institution', &
             'Olivier Gagliardini, IGE UGA CNRS, Grenoble, France'))
         CALL handle_err(nf90_put_att(ncid,NF90_GLOBAL,'source', &
              'Elmer/Ice Version ' // GetVersion() &  
              // ', Rev: ' // GetRevision() & 
              // ', Compiled: ' // GetCompilationDate() ))
         CALL handle_err(nf90_put_att(ncid,NF90_GLOBAL,'references','http://shmip.bitbucket.io/'))
         

     ! End of definitionss
         CALL handle_err(NF90_ENDDEF(ncid))

     END SUBROUTINE CreateNetcdfFile

     SUBROUTINE handle_err(status)
        INTEGER, INTENT (in) :: status
        IF(status /= nf90_noerr) THEN
           PRINT *, trim(nf90_strerror(status))
           CALL FATAL(SolverName,'NETCDF ERROR')
        END IF
     END SUBROUTINE handle_err
       END SUBROUTINE SAVENETCDF

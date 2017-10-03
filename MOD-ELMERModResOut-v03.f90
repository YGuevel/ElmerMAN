!> \ingroup IMNS

! # YG@2016
! # yann.guevel@univ-ubs.fr
! # Part of ANM continuation/detection/branch-switching in ELMER FEM
! # Jean-Marc Cadou / Gregory girault / Yann Guevel
! # Institut de Recherche Dupuy de Lome - Universite Bretagne Sud



MODULE ELMERModResOut




CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE ResultOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils

  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

  LOGICAL :: SaveGid, SaveVTK, SaveOpenDx, SaveGmsh, &
      SaveVTU, SaveEP, SaveAny, ListSet = .FALSE., ActiveMesh, &
      SomeMeshSaved
  INTEGER :: i,nInterval=1, nstep=0, OutputCount = 0, MeshDim,MeshLevel,nlen
  INTEGER, POINTER :: OutputIntervals(:), TimeSteps(:)

  TYPE(Mesh_t), POINTER :: Mesh, iMesh, ListMesh
  CHARACTER(10) :: OutputFormat
  CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix, MeshName, iMeshName
  LOGICAL :: SubroutineVisited=.FALSE.,Found
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: ModelVariables
    
  SAVE SubroutineVisited, OutputCount, ListSet, MeshDim, ListMesh

  INTERFACE
    RECURSIVE SUBROUTINE ElmerPostOutputSolver( Model, Solver,dt,TransientSimulation,ONOEfound )
      USE DefUtils
      IMPLICIT NONE
!------------------------------------------------------------------------------
      TYPE(Solver_t) :: Solver
      TYPE(Model_t) :: Model
      REAL(KIND=dp) :: dt
      LOGICAL, OPTIONAL :: ONOEfound
      LOGICAL :: TransientSimulation
    END SUBROUTINE ElmerPostOutputSolver
  END INTERFACE

  CALL Info( 'ResultOutputSolver', '-------------------------------------')

  Params => Solver % Values
  SaveGid = GetLogical(Params,'Gid Format',Found)
  SaveGmsh = GetLogical(Params,'Gmsh Format',Found)
  SaveVTK = GetLogical(Params,'VTK Format',Found)
  SaveVTU = GetLogical(Params,'VTU Format',Found)
  SaveOpenDx = GetLogical(Params,'Dx Format',Found)
  SaveEP = GetLogical(Params,'Elmerpost Format',Found)

  OutputFormat = GetString( Params, 'Output Format', Found )
  IF(Found) THEN
    IF( OutputFormat == "gid" )THEN
      SaveGid = .TRUE.
    ELSE IF( OutputFormat == "vtk" )THEN
      SaveVTK = .TRUE.
    ELSE IF( OutputFormat == "vtu" )THEN
      SaveVTU = .TRUE.
    ELSE IF( OutputFormat == "dx" )THEN
      SaveOpenDx = .TRUE.
    ELSE IF( OutputFormat == "gmsh" )THEN
      SaveGmsh = .TRUE.
    ELSE IF( OutputFormat == "elmerpost" )THEN
      SaveEP = .TRUE.
    ELSE
      CALL Warn( 'ResultOutputSolver', &
                 'Unknown output format "' // TRIM(OutputFormat) // '"' )
      CALL Warn( 'ResultOutputSolver', &
                 'Available formats are "GiD", "VTK", "VTU", "DX", "gmsh" and "elmerpost"' )
      RETURN
    END IF
  END IF

  IF( .NOT. SubroutineVisited ) THEN
    IF ( GetLogical(Params,'Show Variables',Found) ) THEN
      CALL CreateListForSaving( Model, Params,.TRUE. )    
    END IF
  END IF

  SaveAny = SaveGid .OR. SaveVTK .OR. SaveVTU .OR. SaveOpenDX .OR. SaveGmsh .OR. SaveEp
  IF(.NOT. SaveAny ) THEN
    CALL Warn('ResultOutputSolver','No output format given, assuming VTU')
    SaveVTU = .TRUE.
  END IF

  FilePrefix = GetString( Params,'Output File Name',Found)
!   write(*,*) "RESOUT - FilePrefix=",FilePrefix
!   write(*,*) "RESOUT - Found=",Found  
  IF(.NOT. Found) THEN
    FilePrefix = 'Case'
    CALL ListAddString( Params,'Output File Name',FilePrefix)
  END IF
  
  IF( .NOT. SubroutineVisited ) THEN 
    CALL Info('ResultOutputSolve','Saving with prefix: '//TRIM(FilePrefix))
  END IF        


  ! The idea of this is that the independent subroutines may be called 
  ! with different data sets and still maintaining the standard output calling convention
  OutputCount = OutputCount + 1
  CALL ListAddInteger( Params,'Output Count',OutputCount)


  ! Finally go for it and write desired data
  ! Some formats requite that the list of variables is explicitely given
  !-----------------------------------------

  MeshLevel = GetInteger( Params,'Output Mesh Level',Found)
  SomeMeshSaved = .FALSE.

  iMesh => Model % Meshes
  DO WHILE( ASSOCIATED(iMesh) )
    
    IF ( .NOT. iMesh % OutputActive ) THEN
      CALL Info('ResultOutputSolver','Skipping mesh: '//TRIM(iMesh % Name), Level=10 )
      iMesh => iMesh % next
      CYCLE 
    END IF    

    CALL Info('ResultOutputSolver','Working on mesh: '//TRIM(iMesh % Name), Level=7 )
    WRITE(Message,'(A,I0)') 'Dimension of mesh: ',iMesh % MeshDim 
    CALL Info('ResultOutputSolver',Message) 
    IF( iMesh % MeshDim < 2 ) THEN
      CALL Info('ResultOutputSolver','Skipping meshes with too low dimension')
      iMesh => iMesh % next
      CYCLE
    END IF


    ! Optionally skip the writing of given meshes 
    !---------------------------------------------------------------    
    MeshName = GetString( Params,'Mesh Name',Found )
    IF(Found) THEN
      nlen = StringToLowerCase( iMeshName, iMesh % Name ) 
      nlen = LEN_TRIM(MeshName)
      IF( MeshName(1:nlen) /= iMeshName(1:nlen) ) THEN
        iMesh => iMesh % next
        CYCLE
      END IF
    END IF

    CALL SetCurrentMesh( Model, iMesh )
    ModelVariables => Model % Variables
    Model % Variables => iMesh % variables 

    IF( .NOT. ListSet ) THEN
      CALL Info('ResultOutputSolver','Creating list for saving - if not present')
      CALL CreateListForSaving( Model, Params,.TRUE. )    
      ListSet = .TRUE.
    ELSE IF( MeshDim /= Model % Mesh % MeshDim .OR. .NOT. ASSOCIATED(iMesh, ListMesh)) THEN
      CALL Info('ResultOutputSolver','Recreating list for saving')
      CALL CreateListForSaving( Model, Params,.TRUE.,.TRUE.)
    END IF

    MeshDim = Model % Mesh % MeshDim
    ListMesh => iMesh

    ! In case there are multiple mesh levels one may also save coarser ones
    !----------------------------------------------------------------------
    Mesh => iMesh
    DO i=1,MeshLevel
      Mesh => Mesh % Parent
      IF (.NOT.ASSOCIATED(Mesh)) EXIT
    END DO
    IF ( ASSOCIATED(Mesh)) THEN

    CALL SetCurrentMesh( Model, Mesh )
    Model % Variables => Mesh % variables 
    SomeMeshSaved = .TRUE.

    IF( SaveGid ) THEN
      CALL Info( 'ResultOutputSolver','Saving in GiD format' )    
      CALL GiDOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
    IF( SaveGmsh ) THEN
      CALL Info( 'ResultOutputSolver','Saving in new gmsh format' )        
      CALL GmshOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
!     IF( SaveVTK ) THEN
!       CALL Info( 'ResultOutputSolver','Saving in legacy VTK format' )            
!       CALL VtkOutputSolver( Model,Solver,dt,TransientSimulation )
!     END IF
    IF( SaveVTU ) THEN
      CALL Info( 'ResultOutputSolver','Saving in unstructured VTK XML (.vtu) format' )               
      CALL VtuOutputSolver( Model,Solver,dt,TransientSimulation )
    END IF
!     IF( SaveOpenDx ) THEN
!       CALL Info( 'ResultOutputSolver','Saving in OpenDX format' )                   
!       CALL DXOutputSolver( Model,Solver,dt,TransientSimulation )
!     END IF
!     IF( SaveEP ) THEN
!       CALL Info( 'ResultOutputSolver','Saving in ElmerPost format' )                   
!       CALL ElmerPostOutputSolver( Model,Solver,dt,TransientSimulation )
!     END IF

    CALL Info( 'ResultOutputSolver', '-------------------------------------')
    END IF

    iMesh => iMesh % Next
  END DO

  IF( .NOT. SomeMeshSaved ) THEN
    OutputCount = OutputCount - 1
  END IF
  Model % Variables => ModelVariables

  SubroutineVisited = .TRUE.

END SUBROUTINE ResultOutputSolver


!------------------------------------------------------------------------------
!> A fork of the ElmerPost output utility which also includes use of masks.
!> \ingroup Solvers
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE ElmerPostOutputSolver( Model, Solver,dt,TransientSimulation,ONOEfound )
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL, OPTIONAL :: ONOEfound
    LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
    INTEGER :: TimeCount
    CHARACTER(LEN=512) :: PostFile, FilePrefix
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(Variable_t), POINTER :: Var,Var1,Displacement=>NULL(),MeshUpdate=>NULL(),MaskVar
    
    CHARACTER(LEN=512) :: Row
    CHARACTER(MAX_NAME_LEN) :: Str, DateStr, VarName, Txt
    
    LOGICAL :: gotIt, FreeSurfaceFlag, MoveBoundary, MeshMoved, MaskExists, Found
    INTEGER :: ScalarFields, VectorFields, ONOECount
    REAL(KIND=dp) :: Time,MeshScale, NodeCoords(3)
    REAL(KIND=dp), POINTER :: Values(:), Values2(:), Values3(:)
    INTEGER :: jj,n1,n2,ii,i,j,k,l,n,m,q,Node,DOFs,TimeStep, On_nodes_on_elements, &
         NumberOfNodes, NumberOfElements, ind, Vari, MeshDim, ExtCount
    INTEGER, POINTER :: MaskPerm(:), MaskOrder(:), TimeSteps(:)
    LOGICAL :: MaskAllocated 
    TYPE(Mesh_t), POINTER :: Mesh

    INTEGER, SAVE :: ParallelNodes, SavedCount = 0
    LOGICAL, SAVE :: Visited = .FALSE.
!------------------------------------------------------------------------------
    
    MeshDim = Model % Mesh % MeshDim
    Mesh => Model % Mesh

    IF ( .NOT.PRESENT(ONOEfound) ) THEN
      ExtCount = GetInteger( Solver % Values,'Output Count',GotIt)
    
      IF( GotIt ) THEN
        SavedCount = ExtCount 
      ELSE
        SavedCount = SavedCount + 1
      END IF
      Visited = ( SavedCount > 1 )
    END IF

    Timesteps => ListGetIntegerArray( CurrentModel % Simulation, &
        'Timestep Intervals', GotIt )
    IF( GotIt ) THEN
      TimeCount = SUM(TimeSteps)
    ELSE 
      TimeCount = GetInteger( Model % Simulation,'Steady State Max Iterations')
    END IF
    
    FilePrefix = GetString( Solver % Values,'Output File Name',GotIt )
    IF ( .NOT.GotIt ) FilePrefix = "Output"
    
    IF(INDEX(FilePrefix,'.') == 0) THEN
      WRITE( Postfile,'(A,A)') TRIM(FilePrefix),".ep"
    ELSE
      PostFile = FilePrefix
    END IF

    IF ( PRESENT(ONOEfound) ) PostFile = TRIM(Postfile) // ".el"
    
    IF ( INDEX( PostFile, ':') == 0 .AND. PostFile(1:1) /= '/' .AND. &
        PostFile(1:1) /= Backslash ) THEN
      
      IF ( LEN_TRIM(OutputPath) > 0 ) THEN
        IF ( Visited )  THEN
          OPEN( PostFileUnit,File=TRIM(OutputPath) // '/' // &
              TRIM(PostFile), POSITION='APPEND' )
        ELSE
          OPEN( PostFileUnit,File=TRIM(OutputPath) // '/' // &
              TRIM(PostFile),STATUS='UNKNOWN' )
        END IF
      ELSE
        IF ( Visited ) THEN
          OPEN( PostFileUnit,File=TRIM(PostFile),POSITION='APPEND' )
        ELSE
          OPEN( PostFileUnit,File=TRIM(PostFile),STATUS='UNKNOWN' )
        ENDIF
      END IF
    ELSE
      IF ( Visited  ) THEN
        OPEN( PostFileUnit,File=TRIM(PostFile),POSITION='APPEND' )
      ELSE
        OPEN( PostFileUnit,File=TRIM(PostFile),STATUS='UNKNOWN' )
      END IF
    END IF


    FreeSurfaceFlag = .FALSE.
    MoveBoundary    = .FALSE.
    DO i=1,CurrentModel % NumberOfBCs
      FreeSurfaceFlag = FreeSurfaceFlag .OR. GetLogical( &
          CurrentModel % BCs(i) % Values,'Free Surface', GotIt )
      IF ( FreeSurfaceFlag ) THEN
        MoveBoundary =  GetLogical( &
            CurrentModel % BCs(i) % Values,'Internal Move Boundary', GotIt )
        
        IF ( .NOT. GotIt ) MoveBoundary = .TRUE.
        
        FreeSurfaceFlag = FreeSurfaceFlag .AND. MoveBoundary
      END IF
      
      IF ( FreeSurfaceFlag ) EXIT
    END DO
    
!------------------------------------------------------------------------------
! Initialize stuff for masked saving
!------------------------------------------------------------------------------
    Str = GetString( Solver % Values,'Mask Variable',MaskExists)
    MaskAllocated = .FALSE.
    IF( MaskExists ) THEN
      MaskVar => VariableGet(Solver % Mesh % Variables,TRIM(Str))
      IF( ASSOCIATED(MaskVar)) THEN
        MaskPerm => MaskVar % Perm
        MaskExists = ASSOCIATED(MaskPerm)
      END IF
    ELSE
      IF( MeshDim == 2 ) THEN
        Str = GetString( Solver % Values,'2D Mask Name',GotIt)    
      ELSE IF( MeshDim == 3 ) THEN  
        Str = GetString( Solver % Values,'3D Mask Name',GotIt)    
      END IF
      IF(.NOT. GotIt) Str = GetString( Solver % Values,'Mask Name',GotIt) 
      IF( GotIt ) THEN
        ALLOCATE( MaskPerm( Model % NumberOfNodes ) ) 
        CALL MakePermUsingMask( Model,Solver,Mesh,Str, &
            .FALSE., MaskPerm, NumberOfNodes )
        ParallelNodes = NINT( ParallelReduction( 1.0_dp * NumberOfNodes ) )
        IF( ParallelNodes > 0 ) THEN
          MaskExists = .TRUE.
          MaskAllocated = .TRUE.
        ELSE
          CALL Warn('ElmerPostOutputSolver','Given mask not active: '//TRIM(Str) )
          DEALLOCATE( MaskPerm )
        END IF
      END IF
    END IF
    

    IF(MaskExists) THEN
      CALL Info('ElmerPostOutputSolver','Using > '// TRIM(Str) // ' < as mask variable')
      NumberOfNodes = MAXVAL(MaskPerm)
      ALLOCATE(MaskOrder(NumberOfNodes))
      DO i=1,SIZE(MaskPerm)
        j = MaskPerm(i)
        IF(j > 0) MaskOrder(j) = i
      END DO
      NumberOfElements = 0
      IF ( PRESENT(ONOEfound) ) NumberOfNodes=0
      DO i=1,Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
        Element => Model % Elements(i)
        IF( ALL(MaskPerm(Element % NodeIndexes)>0)) THEN
          NumberOfElements = NumberOfElements + 1
        END IF
        IF ( PRESENT(ONOEfound) ) &
          NumberOfNodes=NumberOfNodes+Element % TYPE % NumberOfNodes
      END DO
    ELSE
      NumberOfNodes = Model % NumberOfNodes
      NumberOfElements =  Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
      IF ( PRESENT(ONOEfound) ) THEN
        NumberOfNodes=0
        DO i=1,Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
          Element => Model % Elements(i)
          NumberOfNodes=NumberOfNodes+Element % TYPE % NumberOfNodes
        END DO
      END IF
    END IF
 
!------------------------------------------------------------------------------
!   Count degrees of freedom to be saved
!------------------------------------------------------------------------------
    
    Dofs = 0
    ScalarFields = 0
    On_nodes_on_elements = 0
    DO Vari = 1, 999
      WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
      VarName=ListGetString( Solver % Values, Txt, Found )
      IF ( .NOT. Found ) EXIT
      Var => VariableGet( Solver % Mesh % Variables, Varname )
      IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
        On_Nodes_on_elements=On_Nodes_on_elements+1
        IF ( PRESENT(ONOEfound) ) Dofs=Dofs+1
      ELSE
        IF ( .NOT.PRESENT(ONOEfound) ) Dofs=Dofs+1
      END IF
      ScalarFields = ScalarFields + 1
    END DO
    
    VectorFields = 0
    DO Vari = 1, 999
      WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
      VarName=ListGetString( Solver % Values, Txt, Found )
      IF ( .NOT. Found ) EXIT
      Var => VariableGet( Solver % Mesh % Variables, Varname )
      IF(.NOT.ASSOCIATED(Var)) & 
        Var => VariableGet( Solver % Mesh % Variables, TRIM(Varname)//' 1' )
      IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
        On_Nodes_on_elements=On_Nodes_on_elements+1
        IF ( PRESENT(ONOEfound) ) Dofs=Dofs+3
      ELSE
        IF ( .NOT.PRESENT(ONOEfound) ) Dofs=Dofs+3
      END IF
      VectorFields = VectorFields+1
    END DO
    
    IF( ScalarFields + VectorFields == 0 ) GOTO 10
    
!------------------------------------------------------------------------------
! Write header to output
!------------------------------------------------------------------------------
    IF ( .NOT. Visited ) THEN
      WRITE(PostFileUnit,'(i10,i10,i7,i7)',ADVANCE='NO' ) NumberOfNodes, &
          NumberOfElements, DOFs, TimeCount
      
      DO Vari = 1, ScalarFields
        WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
        VarName = ListGetString( Solver % Values, TRIM(Txt) )
        Var => VariableGet( Solver % Mesh % Variables,VarName )
        IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
          IF ( .NOT.PRESENT(ONOEfound) ) CYCLE
        ELSE
          IF ( PRESENT(ONOEfound) ) CYCLE
        END IF

        VarName(1:1) = CHAR(ICHAR(VarName(1:1))-ICHAR('a')+ICHAR('A'))
        k = LEN_TRIM(VarName)
        DO j=1,k
          IF ( VarName(j:j) == ' ' ) VarName(j:j) = '.'
        END DO
       
        WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' scalar: '//TRIM(VarName)
      END DO
      
      DO Vari = 1, VectorFields
        WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
        Var => VariableGet( Solver % Mesh % Variables,VarName )
        IF(.NOT.ASSOCIATED(Var)) & 
          Var => VariableGet( Solver % Mesh % Variables, TRIM(Varname)//' 1' )
        IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
          IF ( .NOT.PRESENT(ONOEfound) ) CYCLE
        ELSE
          IF ( PRESENT(ONOEfound) ) CYCLE
        END IF


        VarName(1:1) = CHAR(ICHAR(VarName(1:1))-ICHAR('a')+ICHAR('A'))
        k = LEN_TRIM(VarName)
        DO j=1,k
          IF ( VarName(j:j) == ' ' ) VarName(j:j) = '.'
        END DO

        WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) ' vector: '//TRIM(VarName)
      END DO
      
      WRITE(PostFileUnit,'()')
      DateStr = FormatDate()
      WRITE( PostFileUnit, '("#File started at: ",A)' ) TRIM(DateStr)
!------------------------------------------------------------------------------
!   Coordinates
!------------------------------------------------------------------------------
!
      MeshScale = 1.0_dp
      DO i=1,Model % NumberOfSolvers
        IF ( Model % Solvers(i) % Variable % NameLen <= 0 ) CYCLE
        
        IF (Model % Solvers(i) % Variable % Name &
            (1:Model % Solvers(i) % Variable % NameLen) == 'displacement') THEN
          MeshMoved = ListGetLogical( Model % Solvers(i) % Values, &
              'Displace Mesh', Gotit )

          IF ( .NOT. GotIt ) &
            MeshMoved=.NOT.EigenOrHarmonicAnalysis(Model % Solvers(i))

          IF ( .NOT.MeshMoved ) MeshScale = 0.0_dp
        END IF
      END DO
      
      IF ( PRESENT(ONOEfound) ) THEN
        DO i=1,Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
          Element => Model % Elements(i)
          IF(MaskExists) THEN
            IF( ANY(MaskPerm(Element % NodeIndexes)<=0)) CYCLE
          END IF

          DO j=1,Element % TYPE % NumberOfNodes
            NodeCoords(1) = Model % Nodes % x(Element % NodeIndexes(j))
            NodeCoords(2) = Model % Nodes % y(Element % NodeIndexes(j))
            NodeCoords(3) = Model % Nodes % z(Element % NodeIndexes(j))

            IF( MeshDim <= 2 ) THEN           
              WRITE(PostFileUnit,'(2ES16.7E3,A)') NodeCoords(1:2), ' 0.0'
            ELSE
              WRITE(PostFileUnit,'(3ES16.7E3)') NodeCoords(1:3)
            END IF
          END DO
        END DO
      ELSE
        DO ii=1,NumberOfNodes
          i = ii
          IF(MaskExists) i = MaskOrder(i)

          NodeCoords(1) = Model % Nodes % x(i)
          NodeCoords(2) = Model % Nodes % y(i)
          NodeCoords(3) = Model % Nodes % z(i)
        
          IF ( ASSOCIATED(Displacement) ) THEN
            k = Displacement % Perm(i)
          
            IF ( k > 0 ) THEN
              k = Displacement % DOFs * (k-1)
              NodeCoords(1) = NodeCoords(1) - MeshScale*Displacement % Values(k+1)
              IF( Displacement % DOFs >= 2) THEN
                NodeCoords(2) = NodeCoords(2) - MeshScale*Displacement % Values(k+2)
                END IF
              IF( Displacement % DOFs == 3) THEN
                NodeCoords(3) = NodeCoords(3) - MeshScale*Displacement % Values(k+3)
              END IF
            ELSE
              IF ( ASSOCIATED( MeshUpdate ) ) k = MeshUpdate % Perm(i)
            
              k = MeshUpdate % DOFs * (k-1)
              NodeCoords(1) = NodeCoords(1) - MeshUpdate % Values(k+1)
              IF( MeshUpdate % DOFs >= 2) THEN
                NodeCoords(2) = NodeCoords(2) - MeshUpdate % Values(k+2)
              END IF
              IF( MeshUpdate % DOFs == 3) THEN
                NodeCoords(3) = NodeCoords(3) - MeshUpdate % Values(k+3)
              END IF
            END IF
          END IF

          IF( MeshDim <= 2 ) THEN           
            WRITE(PostFileUnit,'(2ES16.7E3,A)') NodeCoords(1:2), ' 0.0'
          ELSE
            WRITE(PostFileUnit,'(3ES16.7E3)') NodeCoords(1:3)
          END IF
        END DO
      END IF

!------------------------------------------------------------------------------
! Elements
!------------------------------------------------------------------------------
      WRITE(PostFileUnit,'(a)') '#group all'
      ONOECount=0
      DO i=1,Model % NumberOfBulkElements
        Element => Model % Elements(i)
        
        IF(MaskExists) THEN
          IF( ANY(MaskPerm(Element % NodeIndexes)<=0)) CYCLE
        END IF
        
        k = Element % BodyId
        gotIt = .FALSE.
        IF ( k >= 1 .AND. k <= Model % NumberOfBodies ) THEN
          Str = ListGetString( Model % Bodies(k) % Values,'Name',gotIt )
        END IF
        
        IF ( gotIt ) THEN
          k = LEN_TRIM(Str)
          DO j=1,k
            IF ( Str(j:j) == ' ' ) Str(j:j) = '.'
          END DO
          
          WRITE( PostFileUnit,'(a)',ADVANCE='NO' )  Str(1:k)
        ELSE
          WRITE(PostFileUnit,'(a)',ADVANCE='NO' ) 'body'//TRIM(I2S(k))//' '
        END IF
        
        WRITE(PostFileUnit,'(i5)', ADVANCE='NO') Element % TYPE % ElementCode
        n = 0
        DO j=1,Element % TYPE % NumberOfNodes,4
          DO k=1,MIN(4,Element % TYPE % NumberOfNodes-n)
            n = n + 1
            ind = Element % NodeIndexes(n)
            IF(PRESENT(ONOEfound)) THEN
              ONOECount = ONOECount+1
              ind = ONOECount
            ELSE IF(MaskExists) THEN
              ind = MaskPerm(ind)
            END IF
            WRITE(PostFileUnit, '(i8)', ADVANCE='NO')  ind - 1
          END DO
          WRITE( PostFileUnit,'(a)' ) ''
        END DO
      END DO

      DO i=Model % NumberOfBulkElements + 1,Model % NumberOfBulkElements + &
          Model % NumberOfBoundaryElements
        
        Element => Model % Elements(i)
        
        IF(MaskExists) THEN
          IF( ANY(MaskPerm(Element % NodeIndexes) <= 0)) CYCLE
        END IF
        
        k = Element % BoundaryInfo % Constraint
        
        gotIt = .FALSE.
        IF ( k >= 1 .AND. k <= Model % NumberOfBCs ) THEN
          Str = ListGetString( Model % BCs(k) % Values,'Name',gotIt )
        END IF
        
        IF ( gotIt ) THEN
          k = LEN_TRIM(Str)
          DO j=1,k
            IF ( Str(j:j) == ' ' ) Str(j:j) = '.'
          END DO
          
          WRITE( PostFileUnit,'(a)',ADVANCE='NO' )  Str(1:k)
        ELSE
          WRITE( PostFileUnit,'(a)',ADVANCE='NO' ) 'Constraint'//TRIM(I2S(k))//' '
        END IF
        
        WRITE(PostFileUnit,'(i5)', ADVANCE='NO') Element % TYPE % ElementCode
        DO k=1,Element % TYPE % NumberOfNodes
          ind = Element % NodeIndexes(k)
          IF(PRESENT(ONOEfound)) THEN
            ONOECount=ONOECount+1
            ind = ONOECount
          ELSE IF(MaskExists) THEN
            ind = MaskPerm(ind)
          END IF
          WRITE( PostFileUnit, '(i8)', ADVANCE='NO' )  ind-1
        END DO
        WRITE( PostFileUnit,'(a)' ) ''
      END DO
      WRITE(PostFileUnit,'(a)') '#endgroup all'
    END IF
   
!------------------------------------------------------------------------------
!  Save resulst on a timestep (or steady state iteration step)
!------------------------------------------------------------------------------
 
    TimeStep   = SavedCount
    Var => VariableGet( Model % Variables, 'Time' )        
    Time = 1.0d0
    IF ( ASSOCIATED(Var) ) Time = Var % Values(1)

    
    WRITE( PostFileUnit,'(a,i7,i7,ES16.7E3)' ) '#time ',SavedCount,Timestep,Time

    IF ( PRESENT(ONOEfound) ) THEN
       n1=Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
    ELSE
       n1=1
       n2=NumberOfNodes
    END IF

    DO jj=1,n1

      IF( PRESENT(ONOEfound) ) THEN
        Element => Model % Elements(jj)
        IF(MaskExists) THEN
          IF ( ANY(MaskPerm(Element % NodeIndexes)<=0) ) CYCLE
        END IF
        n2=Element % TYPE % NumberOfNodes

        Parent => NULL()
        IF (ASSOCIATED(Element % BoundaryInfo) ) THEN
          Parent => Element % BoundaryInfo % Left
          n = 0
          IF(ASSOCIATED(Parent)) THEN
            DO l=1,GetElementNOFNodes(Parent)
              DO m=1,n2
                IF (Element % NodeIndexes(m)==Parent % NodeIndexes(l)) n=n+1
              END DO
            END DO
          END IF
          IF (n/=n2) Parent => Element % BoundaryInfo % Right
        END IF
      END IF

      DO ii=1,n2
      
        IF ( PRESENT(ONOEfound) ) THEN
          IF ( ASSOCIATED(Parent) ) THEN
            DO l=1,GetElementNOFNOdes(Parent)
              IF (Element % NodeIndexes(ii)==Parent % NodeIndexes(l)) EXIT
            END DO
            i = Parent % DGIndexes(l)
          ELSE
            i = Element % DGIndexes(ii)
          END IF
        ELSE
          i = ii
          IF(MaskExists) i = MaskOrder(i)
        END IF
      
        DO Vari = 1, ScalarFields
          WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
          Var => VariableGet( Model % Variables, TRIM(VarName) ) 
          IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
            IF (.NOT.PRESENT(ONOEfound)) CYCLE
          ELSE IF ( PRESENT(ONOEfound) ) THEN
            CYCLE
          END IF
          Values => Var % Values

          k = i
          IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
          IF ( k > 0 ) THEN
            WRITE(PostFileUnit,'(ES16.7E3)',ADVANCE='NO') Values(k)
          ELSE
            WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
          END IF
        END DO

        DO Vari = 1, VectorFields

          WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
          Var => VariableGet( Model % Variables, VarName ) 

          IF( ASSOCIATED(Var) ) THEN
            IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
              IF (.NOT.PRESENT(ONOEfound) ) CYCLE
            ELSE IF ( PRESENT(ONOEfound) ) THEN
              CYCLE
            END IF
            k = i
            Dofs = Var% Dofs
            IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)
          
            IF( k == 0 ) THEN
              ! Check for the presence of secondary variable (in practice 'mesh update')
              WRITE(Txt,'(A,I0,A)') 'Vector Field ',Vari,' Complement'
              VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
              IF( GotIt ) THEN
                Var => VariableGet( Model % Variables, VarName ) 
                k = Var % Perm(i)
              END IF
            END IF
          
            IF( k > 0 ) THEN
              DO j=1,dofs
                WRITE(PostFileUnit,'(ES16.7E3)',ADVANCE='NO') Var % Values(dofs*(k-1)+j)
              END DO
              IF( dofs < 3 ) THEN
                WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0'
              END IF
            ELSE 
              WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'
            END IF
          
          ELSE
            ! Check for the presence of component vectors given by its components i=1,2,3
            Var => VariableGet( Model % Variables, TRIM(VarName)//' 1' ) 
            IF( ASSOCIATED(Var)) THEN
              IF ( Var % TYPE == Variable_on_nodes_on_elements ) THEN
                IF (.NOT.PRESENT(ONOEfound) ) CYCLE
              ELSE IF ( PRESENT(ONOEfound) ) THEN
                CYCLE
              END IF

              k = i
              IF ( ASSOCIATED(Var % Perm) ) k = Var % Perm(k)

              IF( k == 0 ) THEN
                WRITE(Txt,'(A,I0,A)') 'Vector Field ',Vari,' Complement'
                VarName = ListGetString( Solver % Values, TRIM(Txt),GotIt)
                IF( GotIt ) THEN
                  Var => VariableGet( Model % Variables, TRIM(VarName)//' 1' ) 
                  k = Var % Perm(i)
                END IF
              END IF
              IF( k > 0 ) THEN
                Values => Var % Values
                Var => VariableGet( Model % Variables, TRIM(VarName)//' 2' ) 
                Values2 => Var % Values
                IF( MeshDim == 2 ) THEN
                  WRITE(PostFileUnit,'(2ES16.7E3,A)',ADVANCE='NO') Values(k),&
                      Values2(k),' 0.0'               
                ELSE
                  Var => VariableGet( Model % Variables, TRIM(VarName)//' 3' ) 
                  Values3 => Var % Values
                  WRITE(PostFileUnit,'(3ES16.7E3)',ADVANCE='NO') Values(k),&
                      Values2(k), Values3(k)
                END IF
              ELSE
                WRITE(PostFileUnit,'(A)',ADVANCE='NO') ' 0.0 0.0 0.0'              
              END IF
            END IF
          END IF
        END DO
        WRITE(PostFileUnit,'()')
      END DO
    END DO

!------------------------------------------------------------------------------
!   We are done here close the files and deallocate
!------------------------------------------------------------------------------
10  CONTINUE

    CLOSE(PostFileUnit)
    
    IF(MaskExists) DEALLOCATE(MaskOrder)
    IF(MaskAllocated) DEALLOCATE( MaskPerm )

    IF ( .NOT. PRESENT(ONOEfound) .AND. On_nodes_on_elements>0 ) &
      CALL ElmerPostOutputSolver(Model,Solver,dt,TransientSimulation,.TRUE.)

  END SUBROUTINE ElmerPostOutputSolver
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Saves results in a format understood by the pre-/postprocessing software GiD.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE GiDOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: Solution
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp), POINTER :: Values(:)
  COMPLEX(KIND=dp), POINTER :: CValues(:)
  TYPE(Variable_t), POINTER :: TimeVariable
  TYPE(ValueList_t), POINTER :: SolverParams

  LOGICAL :: Found, CoordinatesWritten = .FALSE.
  LOGICAL :: EigenAnalysis = .FALSE., FirstTimeStep

  INTEGER :: i,j,k,m,n,dim, Code, body_id, ElementCounter, Nloop, Loop, ExtCount
  INTEGER :: tensorComponents, SteadyStep = 0
  INTEGER, PARAMETER :: MaxElemCode = 1000
  INTEGER :: ListElemTypes(MaxElemCode)


  CHARACTER(LEN=1024) :: OutputFile, ResFile, MshFile, Txt, Family, &
       ScalarFieldName, VectorFieldName, TensorFieldName, CompName
  CHARACTER(LEN=1024) :: Txt2, Txt3

  INTEGER :: PyramidMap613(14,4), PyramidMap605(2,4)
  INTEGER :: WedgeMap706(3,4), WedgeMap715(21,4)

!------------------------------------------------------------------------------

  PyramidMap605(1,:) = (/ 3, 5, 4, 1 /)
  PyramidMap605(2,:) = (/ 3, 5, 2, 1 /)

  PyramidMap613(1,:)  = (/ 7, 8, 3, 12 /)
  PyramidMap613(2,:)  = (/ 10, 11, 12, 5 /)
  PyramidMap613(3,:)  = (/ 10, 13, 12, 5 /)
  PyramidMap613(4,:)  = (/ 9, 10, 13, 11 /)
  PyramidMap613(5,:)  = (/ 9, 13, 11, 12 /)
  PyramidMap613(6,:)  = (/ 9, 10, 11, 1 /)
  PyramidMap613(7,:)  = (/ 9, 6, 11, 1 /)
  PyramidMap613(8,:)  = (/ 9, 6, 11, 12 /)
  PyramidMap613(9,:)  = (/ 9, 8, 12, 4 /)
  PyramidMap613(10,:) = (/ 9, 13, 12, 4 /)
  PyramidMap613(11,:) = (/ 7, 9, 8, 12 /)
  PyramidMap613(12,:) = (/ 7, 9, 6, 12 /)
  PyramidMap613(13,:) = (/ 7, 6, 11, 12 /)
  PyramidMap613(14,:) = (/ 7, 6, 11, 2 /)

  WedgeMap706(1,:) = (/ 5, 4, 3, 1 /)
  WedgeMap706(2,:) = (/ 5, 3, 2, 1 /)
  WedgeMap706(3,:) = (/ 5, 6, 4, 3 /)

  WedgeMap715(1,:) = (/ 10, 11, 5, 2 /)
  WedgeMap715(2,:) = (/ 12, 11, 6, 3 /)
  WedgeMap715(3,:) = (/ 12, 10, 4, 1 /)
  WedgeMap715(4,:) = (/ 7, 8, 11, 2 /)
  WedgeMap715(5,:) = (/ 7, 10, 11, 2 /)
  WedgeMap715(6,:) = (/ 13, 14, 11, 5 /)
  WedgeMap715(7,:) = (/ 13, 10, 11, 5 /)
  WedgeMap715(8,:) = (/ 9, 10, 8, 11 /)
  WedgeMap715(9,:) = (/ 9, 7, 10, 8 /)
  WedgeMap715(10,:) = (/ 9, 12, 10, 11 /)
  WedgeMap715(11,:) = (/ 9, 12, 10, 1 /)
  WedgeMap715(12,:) = (/ 9, 7, 10, 1 /)
  WedgeMap715(13,:) = (/ 9, 8, 11, 3 /)
  WedgeMap715(14,:) = (/ 9, 12, 11, 3 /)
  WedgeMap715(15,:) = (/ 15, 12, 10, 11 /)
  WedgeMap715(16,:) = (/ 15, 12, 10, 4 /)
  WedgeMap715(17,:) = (/ 15, 10, 14, 11 /)
  WedgeMap715(18,:) = (/ 15, 13, 10, 14 /)
  WedgeMap715(19,:) = (/ 15, 13, 10, 4 /)
  WedgeMap715(20,:) = (/ 15, 14, 11, 6 /)
  WedgeMap715(21,:) = (/ 15, 12, 11, 6 /)

  SolverParams => GetSolverParams()
  EigenAnalysis = GetLogical( SolverParams, 'Eigen Analysis', Found )

  
  ExtCount = ListGetInteger( Solver % Values,'Output Count',Found)
  IF( Found ) THEN
    SteadyStep = ExtCount 
  ELSE
    SteadyStep = SteadyStep + 1
  END IF
  FirstTimeStep = ( SteadyStep == 1 )


  OutputFile = GetString( Solver % Values, 'Output File Name', Found )
  IF(.NOT. Found) OutputFile = 'Output'

  WRITE(ResFile,'(A,A)') TRIM(OutputFile),'.flavia.res'
  WRITE(MshFile,'(A,A)') TRIM(OutputFile),'.flavia.msh'


  CALL Info('GidOutputSolver','Writing result for GiD postprocessing')
  CALL Info('GidOutputSolver','res-file = :'//TRIM(ResFile) )
  CALL Info('GidOutputSolver','msh-file = :'//TRIM(MshFile) )


  ! Write the GiD msh-file:
  !------------------------
  dim = CoordinateSystemDimension()
  IF( CoordinatesWritten ) GOTO 10

  OPEN( UNIT=10, FILE=MshFile )

  ! First check how many element types are involved in the analysis:
  !-----------------------------------------------------------------
  ListElemTypes = 0
  ElementCounter = 0
!  DO i = 1, Solver % NumberOfActiveElements
!     Element => GetActiveElement(i)
  DO i = 1, Model % NumberOfBulkElements !+ Model % NumberOfBoundaryElements
     Element => Model % Mesh % Elements(i)

     Code = Element % TYPE % ElementCode
     ListElemTypes(Code) = ListElemTypes(Code)+1
  END DO
  PRINT *,'Total number of elements =',SUM(ListElemTypes)

  ! Write the different element types in different blocks:
  !-------------------------------------------------------
  DO i = 1,MaxElemCode
     IF(ListElemTypes(i) == 0) CYCLE
     PRINT *,ListElemTypes(i),'elements of type',i
     n = MOD(i,100)
     IF( INT(i/100) == 1 ) Family = 'Point'
     IF( INT(i/100) == 2 ) Family = 'Linear'
     IF( INT(i/100) == 3 ) Family = 'Triangle'
     IF( INT(i/100) == 4 ) Family = 'Quadrilateral'
     IF( INT(i/100) == 5 ) Family = 'Tetrahedra'
     IF( INT(i/100) == 6 ) THEN
        Family = 'Tetrahedra' ! PYRAMIDS WILL BE SPLITTED
        n = 4                 ! INTO LINEAR TETRAHEDRA
     END IF
     IF( INT(i/100) == 7 ) THEN
        Family = 'Tetrahedra' ! WEDGES WILL BE SPLITTED
        n = 4                 ! INTO LINEAR TETRAHEDRA
     END IF
     IF( INT(i/100) == 8 ) Family = 'Hexahedra'

     WRITE(Txt,'(A,I1,A,A,A,I0)') 'MESH "Elmer Mesh" dimension ',&
             dim,' ElemType ', TRIM(Family),' Nnode ', n

     WRITE(10,'(A)') TRIM(Txt)

     ! Write all node coordinates in the first block:
     !-----------------------------------------------
     IF( .NOT.CoordinatesWritten ) THEN
        WRITE(10,'(A)') 'Coordinates'
        DO j = 1, Model % Mesh % NumberOfNodes
           WRITE(10,'(I6,3ES16.7E3)') j, &
                Model % Mesh % Nodes % x(j), &
                Model % Mesh % Nodes % y(j), &
                Model % Mesh % Nodes % z(j)
        END DO
        WRITE(10,'(A)') 'end coordinates'
        WRITE(10,'(A)') ' '
        CoordinatesWritten = .TRUE.
     END IF

     ! Write the element connectivity tables:
     !---------------------------------------
     WRITE(10,'(A)') 'Elements'
!     DO j = 1, Solver % NumberOfActiveElements
!        Element => GetActiveElement( j )
     DO j = 1, Model % NumberOfBulkElements !+ Model % NumberOfBoundaryElements
        Element => Model % Mesh % Elements(j)

        Code = Element % TYPE % ElementCode
        IF( Code /= i ) CYCLE
        body_id = Element % BodyId

        IF( Code == 613 ) THEN
           ! 13 noded pyramids will be splitted into 14 linear tetraheda
           !------------------------------------------------------------
           DO m = 1,14
              ElementCounter = ElementCounter + 1
              WRITE(10,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(PyramidMap613(m,:)), body_id
           END DO
        ELSEIF( Code == 605 ) THEN
           ! 5 noded pyramids will be splitted into 2 linear tetraheda
           !----------------------------------------------------------
           DO m = 1,2
              ElementCounter = ElementCounter + 1
              WRITE(10,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(PyramidMap605(m,:)), body_id
           END DO           
        ELSEIF( Code == 706 ) THEN
           ! 6 noded wedges will be splitted into 3 linear tetraheda
           !---------------------------------------------------------
           DO m = 1,3
              ElementCounter = ElementCounter + 1
              WRITE(10,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(WedgeMap706(m,:)), body_id
           END DO           
        ELSEIF( Code == 715 ) THEN
           ! 15 noded wedges will be splitted into 21 linear tetraheda
           !----------------------------------------------------------
           DO m = 1,21
              ElementCounter = ElementCounter + 1
              WRITE(10,'(100I10)') ElementCounter, &
                   Element % NodeIndexes(WedgeMap715(m,:)), body_id
           END DO           
        ELSE
           ! Standard elements are understood by GiD as such
           !------------------------------------------------
           ElementCounter = ElementCounter + 1 
           WRITE(10,'(100I10)') ElementCounter, Element % NodeIndexes, body_id
        END IF

     END DO
     WRITE(10,'(A)') 'end elements'
     WRITE(10,'(A)') ' '
  END DO

10 CONTINUE
  
  ! Write the GiD res-file:
  !------------------------
  IF( FirstTimeStep ) THEN
    OPEN(UNIT=10, FILE=ResFile )
    WRITE(10,'(A)') 'GiD Post Result File 1.0'
  ELSE
    OPEN(UNIT=10, FILE=ResFile, POSITION='APPEND' )
  END IF

  Nloop = 1
  IF( EigenAnalysis ) THEN
     Nloop = GetInteger( Solver % Values, 'Eigen System Values', Found )
     IF( .NOT. Found ) Nloop = GetInteger( Solver % Values, 'Number Of EigenModes', Found )
     IF( .NOT.Found ) Nloop = 1
   END IF

   DO Loop = 1, Nloop
     PRINT *,'------------'

     ! First scalar fields:
     !----------------------
     DO i = 1, 999
       WRITE(Txt,'(A,I0)') 'Scalar Field ',i
       
       ScalarFieldName = GetString( Solver % Values, TRIM(Txt), Found )
       IF(.NOT. Found) EXIT
       
       Solution => VariableGet( Model % Mesh % Variables, ScalarFieldName )
       IF( .NOT.ASSOCIATED( Solution ) ) THEN
         PRINT *,'Scalar field "',TRIM(ScalarFieldName),'" not found'
       ELSE
         PRINT *,'Scarar field',i,'= "',TRIM(ScalarFieldName),'"'
         Perm => Solution % Perm
         
         IF( .NOT.EigenAnalysis ) THEN
           Values => Solution % Values
         ELSE
           Cvalues => Solution % EigenVectors(Loop,:)
         END IF
         
         IF( TransientSimulation ) THEN
           TimeVariable => VariableGet( Model % Mesh % Variables, 'Time' )
           PRINT *,'Current time=',TimeVariable % Values(1)
           WRITE(10,'(A,A,A,ES16.7E3,A)') 'Result "',&
               TRIM(ScalarFieldName),'" "Transient analysis" ', &
               TimeVariable % Values(1) ,' Scalar OnNodes'
         ELSE
           IF( .NOT.EigenAnalysis ) THEN
             WRITE(10,'(A,A,A,I2,A)') 'Result "',&
                 TRIM(ScalarFieldName),'" "Steady analysis"',SteadyStep,' Scalar OnNodes'
           ELSE
             WRITE(10,'(A,A,A,I2,A)') 'Result "',&
                 TRIM(ScalarFieldName),'" "Eigen analysis"',Loop,' Scalar OnNodes'
           END IF
         END IF
         
         WRITE(10,'(A,A,A)') 'ComponentNames "',TRIM(ScalarFieldName),'"'
         WRITE(10,'(A)') 'Values'
         DO j = 1, Model % Mesh % NumberOfNodes
           k = Perm(j)
           IF( .NOT.EigenAnalysis ) THEN
             WRITE(10,'(I6,ES16.7E3)') j, Values(k)
           ELSE
             WRITE(10,'(I6,ES16.7E3)') j, REAL(CValues(k))
           END IF
         END DO
         WRITE(10,'(A)') 'end values'
         WRITE(10,'(A)') ' '
       END IF
     END DO

     ! Then vector fields:
     !--------------------
     DO i = 1, 999
       WRITE(Txt,'(A,I0)') 'Vector Field ',i
       
       VectorFieldName = GetString( Solver % Values, TRIM(Txt), Found )
       IF(.NOT. Found) EXIT
       PRINT *,'Vector field',i,'= "',TRIM(VectorFieldName),'"'
       
       IF( TransientSimulation ) THEN
         TimeVariable => VariableGet( Model % Mesh % Variables, 'Time' )
         PRINT *,'Current time=',TimeVariable % Values(1)
         WRITE(10,'(A,A,A,ES16.7E3,A)') 'Result "',&
             TRIM(VectorFieldName),'" "Transient analysis" ', &
             TimeVariable % Values(1) ,' Vector OnNodes'
       ELSE
         IF( .NOT.EigenAnalysis ) THEN
           WRITE(10,'(A,A,A,I2,A)') 'Result "',&
               TRIM(VectorFieldName),'" "Steady analysis"',SteadyStep,' Vector OnNodes'
         ELSE
           WRITE(10,'(A,A,A,I2,A)') 'Result "',&
               TRIM(VectorFieldName),'" "Eigen analysis"',Loop,' Vector OnNodes'
         END IF
       END IF
       
       WRITE(Txt,'(A)') 'ComponentNames '
       DO j = 1, dim
         IF(j<Dim) THEN
           WRITE(Txt,'(A,A,A,I2,A)' ) &
               TRIM(Txt), ' "', TRIM(VectorFieldName),j,'",'
         ELSE
           WRITE(Txt,'(A,A,A,I2,A)' ) &
               TRIM(Txt), ' "', TRIM(VectorFieldName),j,'"'
         END IF
       END DO
       WRITE(10,'(A)') TRIM(Txt)
       WRITE(10,'(A)') 'Values'
       
       DO j = 1, Model % Mesh % NumberOfNodes
         WRITE(Txt2,'(I10)') j
         DO k = 1,dim
           
           ! Check if vector field components have been defined explicitely:
           !----------------------------------------------------------------
           WRITE(Txt3,'(A,I1,A,I1)') 'Vector Field ',i,' component ',k
           CompName = GetString( Solver % Values, TRIM(Txt3), Found )
           IF( Found ) THEN
             WRITE(Txt,'(A)') TRIM(CompName)
           ELSE
             WRITE(Txt,'(A,I2)') TRIM(VectorFieldName), k
           END IF
           
           IF( j==1 ) PRINT *, TRIM(Txt3),' = "', TRIM(Txt),'"'
           
           Solution => VariableGet( Model % Mesh % Variables, TRIM(Txt) )
           IF( .NOT.ASSOCIATED( Solution ) ) THEN
             PRINT *,'Vector field component',k,' not found'
           ELSE
             Perm => Solution % Perm
             
             IF( .NOT.EigenAnalysis ) THEN
               Values => Solution % Values
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), Values( Perm(j) )
             ELSE
               CValues => Solution % Eigenvectors(Loop,:)
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), REAL(CValues( Perm(j) ) )
             END IF
             
           END IF
         END DO
         WRITE(10,'(A)') TRIM(Txt2)
         
       END DO
       WRITE(10,'(A)') 'end values'
     END DO
     

     ! Finally tensor fields:
     !-----------------------
     DO i = 1, 999
       WRITE(Txt,'(A,I0)') 'Tensor Field ',i
       TensorFieldName = GetString( Solver % Values, TRIM(Txt), Found )
       IF(.NOT. Found) EXIT
       
       PRINT *,'Tensor field',i,'= "',TRIM(TensorFieldName),'"'
       
       IF( TransientSimulation ) THEN
         TimeVariable => VariableGet( Model % Mesh % Variables, 'Time' )
         PRINT *,'Current time=',TimeVariable % Values(1)
         WRITE(10,'(A,A,A,ES16.7E3,A)') 'Result "',&
             TRIM(TensorFieldName),'" "Transient analysis" ', &
             TimeVariable % Values(1) ,' Matrix OnNodes'
       ELSE
         IF( .NOT.EigenAnalysis ) THEN
           WRITE(10,'(A,A,A,I2,A)') 'Result "',&
               TRIM(TensorFieldName),'" "Steady analysis"',SteadyStep,' Matrix OnNodes'
         ELSE
           WRITE(10,'(A,A,A,I2,A)') 'Result "',&
               TRIM(TensorFieldName),'" "Eigen analysis"',Loop,' Matrix OnNodes'
         END IF
         
       END IF
       
       WRITE(Txt,'(A)') 'ComponentNames '
       IF( dim == 2 ) THEN
         TensorComponents = 3
       ELSE
         TensorComponents = 6
       END IF
       
       DO j = 1, TensorComponents
         WRITE(Txt3,'(A,I1,A,I1)') 'Tensor Field ',i,' component ',j
         CompName = GetString( Solver % Values, TRIM(Txt3), Found )

         IF( Found ) THEN
           WRITE(Txt2,'(A)') TRIM(CompName)
         ELSE
           WRITE(Txt2,'(A,A,I1)') TRIM(TensorFieldName), ' ', j
         END IF

         IF(j<TensorComponents) THEN
           WRITE(Txt,'(A,A,A,A)' ) &
               TRIM(Txt), ' "', TRIM(Txt2), '",'
         ELSE
           WRITE(Txt,'(A,A,A,A)' ) &
               TRIM(Txt), ' "', TRIM(txt2), '"'
         END IF
       END DO
       WRITE(10,'(A)') TRIM(Txt)
       WRITE(10,'(A)') 'Values'
       
       DO j = 1, Model % Mesh % NumberOfNodes
         WRITE(Txt2,'(I10)') j
         DO k = 1,TensorComponents
           
           ! Check if tensor field components have been defined explicitely:
           !----------------------------------------------------------------
           WRITE(Txt3,'(A,I1,A,I1)') 'Tensor Field ',i,' component ',k
           CompName = GetString( Solver % Values, TRIM(Txt3), Found )
           IF( Found ) THEN
             WRITE(Txt,'(A)') TRIM(CompName)
           ELSE
             WRITE(Txt,'(A,A,I1)') TRIM(TensorFieldName), ' ', k
           END IF
           
           IF( j==1 ) PRINT *, TRIM(Txt3),' = "', TRIM(Txt),'"'
           
           Solution => VariableGet( Model % Mesh % Variables, TRIM(Txt) )
           IF( .NOT.ASSOCIATED( Solution ) ) THEN
             PRINT *,'Tensor field component',k,' not found'
           ELSE
             Perm => Solution % Perm
             !Values => Solution % Values
             !WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), Values( Perm(j) )
             
             IF( .NOT.EigenAnalysis ) THEN
               Values => Solution % Values
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), Values( Perm(j) )
             ELSE
               CValues => Solution % Eigenvectors(Loop,:)
               WRITE(Txt2,'(A,ES16.7E3)') TRIM(Txt2), REAL(CValues( Perm(j) ) )
             END IF

           END IF
         END DO
         WRITE(10,'(A)') TRIM(Txt2)
         
       END DO
       WRITE(10,'(A)') 'end values'
     END DO
     
   END DO ! Nloop


   CLOSE(10)
  
   PRINT *,'Output complete.'

!------------------------------------------------------------------------------
 END SUBROUTINE GiDOutputSolver
!------------------------------------------------------------------------------
  



!------------------------------------------------------------------------------
!> Saves results in ascii format understood by the pre-/postprocessing software Gmsh.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE GmshOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: Perm(:)
  REAL(KIND=dp), POINTER :: Values(:),Values2(:),Values3(:)
  REAL(KIND=dp) :: Vector(3), Time
  COMPLEX(KIND=dp), POINTER :: CValues(:)
  TYPE(Variable_t), POINTER :: Solution, TimeVariable
  TYPE(ValueList_t), POINTER :: SolverParams
  
  LOGICAL :: Found, GotField, FileAppend, AlterTopology, MaskExists
  LOGICAL :: EigenAnalysis = .FALSE., EigenActive, ComponentVector
  
  INTEGER :: VisitedTimes = 0, ExtCount
  INTEGER :: i,j,k,l,m,n,nsize,dim,dofs,ElmerCode, GmshCode,body_id, Vari, Rank, truedim
  INTEGER :: Tag, NumberOfAllElements, BCOffSet
  INTEGER, PARAMETER :: MaxElemCode = 827
  INTEGER :: ElmerToGmshType(MaxElemCode), GmshToElmerType(21), GmshIndexes(27) 
  INTEGER, POINTER :: NodeIndexes(:), ElmerIndexes(:), MaskPerm(:)

  INTEGER, PARAMETER :: LENGTH = 1024
  CHARACTER(LEN=LENGTH) :: OutputFile, Txt, FieldName, CompName
    
  SAVE VisitedTimes
  
!------------------------------------------------------------------------------

  CALL Info('GmshOutputSolver','Saving results in Gmsh format')

  ExtCount = ListGetInteger( Solver % Values,'Output Count',Found)
  IF( Found ) THEN
    VisitedTimes = ExtCount
  ELSE
    VisitedTimes = VisitedTimes + 1
  END IF

  GmshToElmerType = (/ 202, 303, 404, 504, 808, 706, 605, 203, 306, 409, &
      510, 827, 0, 0, 101, 408, 820, 715, 613, 0, 310 /)
  ElmerToGmshType = 0

  DO i=1,SIZE(GmshToElmerType)
    j = GmshToElmerType(i)
    IF( j > 0 ) ElmerToGmshType(j) = i
  END DO

  SolverParams => GetSolverParams()
  EigenAnalysis = GetLogical( SolverParams, 'Eigen Analysis', Found )
  FileAppend = GetLogical( SolverParams,'File Append',Found)
  IF(.NOT. Found) FileAppend = .TRUE.
  AlterTopology = GetLogical( SolverParams,'Alter Topology',Found)
  
  Txt = ListGetString( SolverParams,'Mask Variable',MaskExists)
  IF( MaskExists ) THEN
    Solution => VariableGet(Model % Variables,TRIM(Txt))
    IF( ASSOCIATED(Solution)) MaskPerm => Solution % Perm
    MaskExists = ASSOCIATED(MaskPerm)
  END IF

  OutputFile = GetString( Solver % Values, 'Output File Name', Found )
  IF( .NOT.Found ) OutputFile = 'Output.msh'

  IF(INDEX(OutputFile,'.') == 0) WRITE( OutputFile,'(A,A)') TRIM(OutputFile),".msh"

  dim = CoordinateSystemDimension()
  IF( VisitedTimes > 1 ) THEN
    IF( AlterTopology ) THEN
      OutputFile = NextFreeFilename( OutputFile )
      CALL Info('GmshOutputSolver','Writing mesh and data to a new file: '//TRIM(OutputFile))
    ELSE IF( FileAppend ) THEN      
      CALL Info('GmshOutputSolver','Appending data to the same file: '//TRIM(OutputFile))
      OPEN(UNIT=10, FILE=OutputFile, POSITION='APPEND' )      
      GOTO 10
    ELSE
      OutputFile = NextFreeFilename( OutputFile )          
      CALL Info('GmshOutputSolver','Writing data to a new file: '//TRIM(OutputFile))
      OPEN(UNIT=10, FILE=OutputFile )
      WRITE(10,'(A)') '$MeshFormat'
      WRITE(10,'(A)') '2.0 0 8'
      WRITE(10,'(A)') '$EndMeshFormat'          
      GOTO 10    
    END IF
  END IF


  ! Save the header
  !-------------------------------------------------
  CALL Info('GsmhOutputSolver','Saving results to file: '//TRIM(OutputFile))
  OPEN(UNIT=10, FILE=OutputFile )
  
  WRITE(10,'(A)') '$MeshFormat'
  WRITE(10,'(A)') '2.0 0 8'
  WRITE(10,'(A)') '$EndMeshFormat'    
  

  ! Save the mesh nodes
  !-------------------------------------------------
  CALL Info('GmshOutputSolver','Writing the mesh nodes')
  IF( MaskExists ) THEN
    nsize = MAXVAL( MaskPerm ) 
  ELSE
    nsize = Model % NumberOfNodes
  END IF

  WRITE(10,'(A)') '$Nodes'
  WRITE(10,'(I8)') nsize
  IF( dim == 3 ) THEN
    DO i = 1, Model % NumberOfNodes
      IF( MaskExists ) THEN
        IF( MaskPerm(i) == 0 ) CYCLE
      END IF      
      WRITE(10,'(I8,3ES16.7E3)') i,Model % Nodes % x(i),Model % Nodes % y(i), Model % Nodes % z(i)
    END DO
  ELSE 
    DO i = 1, Model % NumberOfNodes
      IF( MaskExists ) THEN
        IF( MaskPerm(i) == 0 ) CYCLE
      END IF            
      WRITE(10,'(I8,2ES16.7E3,A)') i,Model % Nodes % x(i),Model % Nodes % y(i),' 0.0' 
    END DO
  END IF
  WRITE(10,'(A)') '$EndNodes'

  ! Save the mesh elements
  !-------------------------------------------------
  CALL Info('GmshOutputSolver','Writing the mesh elements')
  NumberOfAllElements = Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
  
  IF( MaskExists ) THEN
    nsize = 0
    DO i=1,NumberOfAllElements
      Element => Model % Mesh % Elements(i)
      ElmerIndexes => Element % NodeIndexes
      IF( ANY(MaskPerm(ElmerIndexes) == 0) ) CYCLE
      nsize = nsize  + 1
    END DO
  ELSE
    nsize = NumberOfAllElements
  END IF

  BCOffSet = 100
  DO WHILE( BCOffset <= Model % NumberOfBodies ) 
    BCOffset = 10 * BCOffset
  END DO

  WRITE(10,'(A)') '$Elements'
  WRITE(10,'(I8)') nsize
  DO i = 1, NumberOfAllElements
    Element => Model % Mesh % Elements(i)
    ElmerCode = Element % TYPE % ElementCode
    ElmerIndexes => Element % NodeIndexes
    
    IF( MaskExists ) THEN
      IF( ANY(MaskPerm(ElmerIndexes) == 0) ) CYCLE     
    END IF

    GmshCode = ElmerToGmshType(ElmerCode)
    IF( GmshCode == 0 ) THEN
      CALL Warn('GmshOutputSolver','Gmsh element index not found!')
      CYCLE
    END IF

    IF( i <= Model % NumberOfBulkElements ) THEN
      Tag = Element % BodyId
    ELSE
      Tag = GetBCId( Element ) + BCOffset
    END IF

    WRITE(10,'(I8,I3,I3,I5,I5)',ADVANCE='NO') i,GmshCode,2,Tag,Tag
    k = MOD(ElmerCode,100)

    CALL ElmerToGmshIndex(ElmerCode,ElmerIndexes,GmshIndexes)

    DO j=1,k-1
      WRITE(10,'(I8)',ADVANCE='NO') GmshIndexes(j)
    END DO
    WRITE(10,'(I8)') GmshIndexes(k)
  END DO
  WRITE(10,'(A)') '$EndElements'

  ! With a mask the list of physical entities should be checked
  !-------------------------------------------------------------
  IF(.NOT. MaskExists ) THEN
    nsize = Model % NumberOfBodies + Model % NumberOfBCs
    WRITE(10,'(A)') '$PhysicalNames'
    WRITE(10,'(I8)') nsize
    DO i=1,Model % NumberOfBodies 
      Txt = ListGetString( Model % Bodies(i) % Values,'Name',Found)
      IF( Found ) THEN
        WRITE(10,'(I8,A)') i,'"'//TRIM(Txt)//'"'
      ELSE
        WRITE(10,'(I8,A,I0,A)') i,'"Body ',i,'"'       
      END IF
    END DO
    DO i=1,Model % NumberOfBCs
      Txt = ListGetString( Model % BCs(i) % Values,'Name',Found)
      IF( Found ) THEN
        WRITE(10,'(I8,A)') i+BCOffset,'"'//TRIM(Txt)//'"'
      ELSE
        WRITE(10,'(I8,A,I0,A)') i+BCOffset,'"Boundary Condition ',i,'"'               
      END IF
    END DO
    WRITE(10,'(A)') '$EndPhysicalNames'
  END IF


10 CONTINUE


  ! Time is needed
  !-------------------------------------------------
  TimeVariable => VariableGet( Model % Variables, 'Time' )        
  Time = TimeVariable % Values(1)
  
  ! Loop over different type of variables
  !-------------------------------------------------
  CALL Info('GmshOutputSolver','Writing the nodal data')
  DO Rank = 0,2
    DO Vari = 1, 999
      IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
      IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
      IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari

      FieldName = GetString( Solver % Values, TRIM(Txt), Found )
      IF(.NOT. Found) EXIT 
      IF( Rank == 2) THEN
        CALL Warn('GmshOutputSolver','Not implemented yet for tensors!')
        CYCLE
      END IF

      ComponentVector = .FALSE.
      Solution => VariableGet( Model % Mesh % Variables, FieldName )
      IF(ASSOCIATED(Solution)) THEN
        Values => Solution % Values
        Perm => Solution % Perm
        dofs = Solution % DOFs
      ELSE
        IF( Rank == 1 ) THEN
          Solution => VariableGet( Model % Mesh % Variables, FieldName//' 1' )
          IF( ASSOCIATED( Solution ) ) THEN
            ComponentVector = .TRUE.
            Values => Solution % Values
            Perm => Solution % Perm
            dofs = 1
            Solution => VariableGet( Model % Mesh % Variables, FieldName//' 2' )
            IF( ASSOCIATED(Solution)) THEN
              Values2 => Solution % Values
              dofs = 2
            END IF            
            Solution => VariableGet( Model % Mesh % Variables, FieldName//' 3' )
            IF( ASSOCIATED(Solution)) THEN
              Values3 => Solution % Values
              dofs = 3
            END IF
          END IF
        END IF
        IF( .NOT. ASSOCIATED(Solution)) THEN
          CALL Warn('GsmhOutputSolver','Variable not present: '//TRIM(FieldName))
          CYCLE
        END IF
      END IF
      IF( ASSOCIATED(Solution % EigenVectors) ) THEN
        CALL Warn('GmshOutputSolver','Eigenvectors related to field: '//TRIM(FieldName))
        CALL Warn('GmshOutputSolver','Eigenvectors saving yet not supported')
      END IF

      truedim = MIN(dofs, dim)
      IF( MaskExists ) THEN
        nsize = 0
        DO i=1,SIZE(Perm)
          IF( Perm(i)==0 .OR. MaskPerm(i) == 0 ) CYCLE
          nsize = nsize + 1
        END DO
        IF( nsize == 0 ) THEN
          CALL Warn('GmshOutputSolver','No dofs with the current mask for saving: '//TRIM(FieldName))         
        END IF
      ELSE
        nsize = SIZE(Values) / Dofs
      END IF

      
      WRITE(10,'(A)') '$NodeData'
      WRITE(10,'(A)') '1'
      WRITE(10,'(A)') '"'//TRIM(FieldName)//'"'
      WRITE(10,'(A)') '1'

      ! Gmsh starts steady state indexes from zero, hence deductions by one
      IF( TransientSimulation ) THEN
        WRITE(10,'(ES16.7E3)') Time
      ELSE
        WRITE(10,'(ES16.7E3)') Time - 1.0_dp
      END IF
      WRITE(10,'(A)') '3'
      WRITE(10,'(I8)') VisitedTimes-1
      IF(Rank == 0) THEN
        WRITE(10,'(A)') '1'
      ELSE IF(Rank == 1) THEN
        WRITE(10,'(A)') '3'
      ELSE 
        WRITE(10,'(A)') '9'
      END IF     
      WRITE(10,'(I8)') nsize
     
      DO i=1,SIZE(Perm) 
        j = Perm(i)
        IF( j == 0) CYCLE
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0 ) CYCLE
        END IF
        
        IF( Rank == 0 ) THEN
          WRITE(10,'(I8,ES16.7E3)') i,Values(j)
        ELSE IF(Rank == 1) THEN
          IF( ComponentVector ) THEN
            IF( truedim == 2 ) THEN
              WRITE(10,'(I8,2ES16.7E3,A)') i,&
                  Values(j),Values2(j),' 0.0'
            ELSE
              WRITE(10,'(I8,3ES16.7E3)') i,&
                  Values(j),Values2(j),Values3(j)
            END IF
          ELSE
            IF( truedim == 2 ) THEN
              WRITE(10,'(I8,2ES16.7E3,A)') i,&
                  Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),' 0.0'
            ELSE
              WRITE(10,'(I8,3ES16.7E3)') i,&
                  Values(dofs*(j-1)+1),Values(dofs*(j-1)+2),Values(dofs*(j-1)+3)
            END IF           
          END IF
        END IF
      END DO
      WRITE(10,'(A)') '$EndNodeData'

    END DO
  END DO
  
      
  IF(.FALSE.) THEN
    WRITE(10,'(A)') '$ElementData'
    WRITE(10,'(A)') '$EndElementData'
  END IF
  
  IF(.FALSE.) THEN
    WRITE(10,'(A)') '$ElementNodeData'
    WRITE(10,'(A)') '$EndElementNodeData'
  END IF
  
  CLOSE(10)
  
  CALL Info('GmshOutputSolver','Gmsh output complete')

CONTAINS



  SUBROUTINE ElmerToGmshIndex(Code,ElmerIndexes,GmshIndexes)

    INTEGER :: Code
    INTEGER :: ElmerIndexes(:),GmshIndexes(:)
    INTEGER :: i,n
    LOGICAL :: reorder, Visited = .FALSE.

    INTEGER, TARGET :: order510(10),order613(13),order715(15),order820(20)
    INTEGER, POINTER :: order(:)

    SAVE Visited

    IF(.NOT. Visited ) THEN
      order510(:) = (/ 0,1,2,3,4,5,6,7,9,8 /)
      order613(:) = (/ 0,1,2,3,4,5,8,10,6,7,9,11,12 /)
      order715(:) = (/ 0,1,2,3,4,5,6,9,7,8,10,11,12,14,13 /)
      order820(:) = (/ 0,1,2,3,4,5,6,7,8,11,12,9,10,12,14,15,16,18,19,17 /)
      Visited = .TRUE.
    END IF

    reorder = .FALSE.

    SELECT CASE( Code )
      
    CASE (510)
      reorder = .TRUE.
      order => order510
      
    CASE (613)
      reorder = .TRUE.
      order => order613
      
    CASE (715)
      reorder = .TRUE.
      order => order715
      
    CASE (820)
      reorder = .TRUE.
      order => order820
     
    CASE DEFAULT
      
    END SELECT

    n = MOD(Code,100) 
    IF( reorder ) THEN
      DO i=1,n 
        GmshIndexes(order(i)+1) = ElmerIndexes(i)
      END DO
    ELSE
      GmshIndexes(1:n) = ElmerIndexes(1:n)      
    END IF


  END SUBROUTINE ElmerToGmshIndex

!------------------------------------------------------------------------------
END SUBROUTINE GmshOutputSolver
!------------------------------------------------------------------------------
  

! ! ! ! ! ! 
! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !> Subroutine for legacy VTK output. 
! ! ! ! ! ! !> Note that this has been replaced by the more concurrent XML VTK cased formats.
! ! ! ! ! ! !> \ingroup Solvers
! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! SUBROUTINE VtkOutputSolver( Model,Solver,dt,TransientSimulation )
! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! !   
! ! ! ! ! !   USE DefUtils 
! ! ! ! ! !   USE VtkLegacyFile
! ! ! ! ! !   
! ! ! ! ! !   IMPLICIT NONE
! ! ! ! ! !   TYPE(Solver_t) :: Solver
! ! ! ! ! !   TYPE(Model_t) :: Model
! ! ! ! ! !   REAL(dp) :: dt
! ! ! ! ! !   LOGICAL :: TransientSimulation
! ! ! ! ! !   
! ! ! ! ! !   INTEGER, SAVE :: nTime = 0
! ! ! ! ! !   LOGICAL :: GotIt
! ! ! ! ! !   CHARACTER(MAX_NAME_LEN), SAVE :: FilePrefix
! ! ! ! ! !   
! ! ! ! ! !   ! Avoid compiler warings about unused variables
! ! ! ! ! !   IF ( TransientSimulation ) THEN; ENDIF
! ! ! ! ! !     IF ( dt > 0.0 ) THEN; ENDIF
! ! ! ! ! !       
! ! ! ! ! !       IF ( nTime == 0 ) THEN
! ! ! ! ! !         FilePrefix = GetString( Solver % Values,'Output File Name',GotIt )
! ! ! ! ! !         IF ( .NOT.GotIt ) FilePrefix = "Output"
! ! ! ! ! !       END IF
! ! ! ! ! !       nTime = nTime + 1
! ! ! ! ! !       
! ! ! ! ! !       CALL WriteData( TRIM(FilePrefix), Model, nTime )
! ! ! ! ! !       
! ! ! ! ! ! 
! ! ! ! ! !     CONTAINS
! ! ! ! ! !       
! ! ! ! ! !       SUBROUTINE WriteData( Prefix, Model, nTime )
! ! ! ! ! !         CHARACTER(*), INTENT(IN) :: Prefix
! ! ! ! ! !         TYPE(Model_t) :: Model
! ! ! ! ! !         INTEGER, INTENT(IN) :: nTime
! ! ! ! ! !         CHARACTER(MAX_NAME_LEN) :: VtkFile
! ! ! ! ! !         TYPE(Mesh_t), POINTER :: Mesh
! ! ! ! ! !         TYPE(Variable_t), POINTER :: Var
! ! ! ! ! !         INTEGER :: i, j, k
! ! ! ! ! !         LOGICAL :: EigAnal
! ! ! ! ! !         REAL(dp), POINTER :: OrigValues(:)
! ! ! ! ! !         INTEGER :: OrigDOFs
! ! ! ! ! !         CHARACTER(MAX_NAME_LEN) :: Dir
! ! ! ! ! !         
! ! ! ! ! !         Mesh => Model % Mesh
! ! ! ! ! !           
! ! ! ! ! !         IF (LEN_TRIM(Mesh % Name) > 0 ) THEN
! ! ! ! ! !           Dir = TRIM(Mesh % Name) // "/"
! ! ! ! ! !         ELSE
! ! ! ! ! !           Dir = "./"
! ! ! ! ! !         END IF
! ! ! ! ! !           
! ! ! ! ! !         EigAnal = .FALSE.
! ! ! ! ! !         
! ! ! ! ! !         Solvers: DO i = 1, Model % NumberOfSolvers
! ! ! ! ! !           EigAnal = ListGetLogical( Model % Solvers(i) % Values, &
! ! ! ! ! !               "Eigen Analysis", GotIt )
! ! ! ! ! !           Var => Model % Solvers(i) % Variable
! ! ! ! ! !           IF ( EigAnal .AND. ASSOCIATED(Var % EigenValues) ) THEN
! ! ! ! ! !             DO j = 1, Model % Solvers(i) % NOfEigenValues
! ! ! ! ! !               OrigValues => Var % Values
! ! ! ! ! !               OrigDOFs = Var % DOFs
! ! ! ! ! !               
! ! ! ! ! !               IF ( Model % Solvers(i) % Matrix % COMPLEX ) THEN
! ! ! ! ! !                 Var % DOFs = Var % DOFs*2
! ! ! ! ! !                 ALLOCATE( Var % Values(2*SIZE(Var%EigenVectors,2)) )
! ! ! ! ! !                 FORALL ( k = 1:SIZE(Var % Values)/2 )
! ! ! ! ! !                   Var%Values(2*k-1) = REAL(Var%EigenVectors(j,k))
! ! ! ! ! !                   Var%Values(2*k) = AIMAG(Var%EigenVectors(j,k))
! ! ! ! ! !                 END FORALL
! ! ! ! ! !               ELSE
! ! ! ! ! !                 ALLOCATE( Var % Values(SIZE(Var % EigenVectors,2)) )
! ! ! ! ! !                 Var % Values = Var % EigenVectors(j,:)
! ! ! ! ! !               END IF
! ! ! ! ! !               
! ! ! ! ! !               WRITE( VtkFile, '(A,A,I4.4,"_",I3.3,".vtk")' ) &
! ! ! ! ! !                   TRIM(Dir), Prefix, nTime, j
! ! ! ! ! !               CALL WriteVtkLegacyFile( VtkFile, Model, .FALSE. )
! ! ! ! ! !               
! ! ! ! ! !               DEALLOCATE( Var % Values )
! ! ! ! ! !               Var % Values => OrigValues
! ! ! ! ! !               Var % DOFs = OrigDOFs
! ! ! ! ! !             END DO
! ! ! ! ! !             EXIT Solvers
! ! ! ! ! !           END IF
! ! ! ! ! !         END DO Solvers
! ! ! ! ! !         
! ! ! ! ! !         IF ( .NOT.EigAnal ) THEN
! ! ! ! ! !           WRITE( VtkFile,'(A,A,I4.4,".vtk")' ) TRIM(Dir),Prefix,nTime
! ! ! ! ! !           CALL WriteVtkLegacyFile( VtkFile, Model, .TRUE. )
! ! ! ! ! !         END IF
! ! ! ! ! !         
! ! ! ! ! !       END SUBROUTINE WriteData
! ! ! ! ! !       
! ! ! ! ! !     END SUBROUTINE VtkOutputSolver







!------------------------------------------------------------------------------
!> Subroutine for saving the results in XML based VTK format (VTU). Both ascii and binary
!> output is available, in single or double precision. The format is understood by 
!> visualization softwares Paraview and ViSit, for example.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE VtuOutputSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

  USE DefUtils 
  USE MeshUtils
  USE ElementDescription
  USE AscBinOutputUtils
  
  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(dp) :: dt
  LOGICAL :: TransientSimulation
  
  INTEGER, SAVE :: nTime = 0
  LOGICAL :: GotIt, Hit, Parallel, FixedMesh
  CHARACTER(MAX_NAME_LEN), SAVE :: FilePrefix
  CHARACTER(MAX_NAME_LEN) :: BaseFile, VtuFile, PvtuFile, PvdFile, DataSetFile
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Variable_t), POINTER :: Var
  INTEGER :: i, j, k, l, n, Partitions, Part, ExtCount, FileindexOffSet, MeshDim, PrecBits, &
             PrecSize, IntSize, FileIndex
  CHARACTER(MAX_NAME_LEN) :: Dir
  LOGICAL :: Visited = .FALSE.
  REAL(KIND=dp) :: DoubleWrk
  REAL :: SingleWrk

  LOGICAL :: MaskExists, BinaryOutput, AsciiOutput, SinglePrec, NoFileindex
  CHARACTER(MAX_NAME_LEN) :: Str, MaskName
  TYPE(Variable_t), POINTER :: MaskVar
  INTEGER, POINTER :: MaskPerm(:), NodePerm(:),NodeOrder(:)
  INTEGER :: NumberOfNodes, NumberOfElements, ParallelNodes, ParallelElements
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(ValueList_t),POINTER :: Params
  INTEGER :: MaxModes, BCOffset, ElemFirst, ElemLast
  INTEGER, POINTER :: ActiveModes(:), Indexes(:)
  LOGICAL :: GotActiveModes, EigenAnalysis, WriteIds, SaveBoundariesOnly, SaveBulkOnly, &
      GotMaskName, AllNodesUsed, SaveElemental, SaveNodal, GotMaskCond
  LOGICAL, ALLOCATABLE :: ActiveElem(:)
  REAL(KIND=dp), ALLOCATABLE :: MaskCond(:)

! Parameters for buffered binary output
  INTEGER :: BufferSize

  Params => GetSolverParams()

  ExtCount = GetInteger( Params,'Output Count',GotIt)
  IF( GotIt ) THEN
    nTime = ExtCount
  ELSE
    nTime = nTime + 1
  END IF

  FileIndexOffset = GetInteger( Params,'Fileindex offset',GotIt)
  FileIndex = nTime + FileIndexOffset

  BinaryOutput = GetLogical( Params,'Binary Output',GotIt)
  IF( GotIt ) THEN
    AsciiOutput = .NOT. BinaryOutput
  ELSE
    AsciiOutput = GetLogical( Params,'Ascii Output',GotIt)
    BinaryOutput = .NOT. AsciiOutput
  END IF

  SaveElemental = GetLogical( Params,'Save Elemental Fields',GotIt)
  IF(.NOT. GotIt) SaveElemental = .TRUE.

  SaveNodal = GetLogical( Params,'Save Nodal Fields',GotIt) 
  IF(.NOT. GotIt) SaveNodal = .TRUE.

  SinglePrec = GetLogical( Params,'Single Precision',GotIt) 
  IF( SinglePrec ) THEN
    CALL Info('VtuOutputSolver','Using single precision arithmetics in output!',Level=7)
  END IF

  IF( SinglePrec ) THEN
    PrecBits = 32
    PrecSize = KIND( SingleWrk ) 
  ELSE
    PrecBits = 64
    PrecSize = KIND( DoubleWrk ) 
  END IF
  IntSize = KIND(i)

  Partitions = ParEnv % PEs
  Part = ParEnv % MyPE
  Parallel = (Partitions > 1) .OR. GetLogical(Params,'Enforce Parallel format',GotIt)

  NoFileindex = GetLogical( Params,'No Fileindex',GotIt)

  IF ( nTime == 1 ) THEN
    FilePrefix = GetString( Params,'Output File Name',GotIt )
    IF ( .NOT.GotIt ) FilePrefix = "Output"
    CALL Info('VtuOutputSolver','Saving results in VTK XML format with prefix: '//TRIM(FilePrefix))
    WRITE( Message,'(A,I0)') 'Saving number of partitions: ',Partitions
    CALL Info('VtuOutputSolver', Message )
  END IF

  FixedMesh = ListGetLogical(Params,'Fixed Mesh',GotIt)
  Mesh => Model % Mesh
  MeshDim = Mesh % MeshDim


  !------------------------------------------------------------------------------
  ! Initialize stuff for masked saving
  !------------------------------------------------------------------------------
  GotMaskName = .FALSE.
  Str = GetString( Params,'Mask Variable',MaskExists)
  IF( MaskExists ) THEN
    MaskVar => VariableGet(Model % Variables,TRIM(Str))
    IF( ASSOCIATED(MaskVar)) MaskPerm => MaskVar % Perm
    MaskExists = ASSOCIATED(MaskPerm)
    IF( MaskExists ) THEN
      CALL Info('VtuOutputSolver','Using > '// TRIM(Str) // ' < as mask variable')
    END IF
  ELSE
    ! Check if there is an additional mask name given
    IF( MeshDim == 2 ) THEN
      MaskName = GetString( Params,'2D Mask Name',GotIt)    
    ELSE IF( MeshDim == 3 ) THEN  
      MaskName = GetString( Params,'3D Mask Name',GotIt)    
    END IF
    IF(.NOT. GotIt) MaskName = GetString( Params,'Mask Name',GotIt) 
    GotMaskName = GotIt
  END IF

  GotMaskCond = .FALSE.
  IF( .NOT. GotMaskName ) THEN
    MaskName = GetString( Params,'Mask Condition',GotMaskCond)
    IF( GotMaskCond ) THEN
      n = Mesh % MaxElementNodes
      ALLOCATE( MaskCond(n) )
    END IF
  END IF

  SaveBoundariesOnly = GetLogical( Params,'Save Boundaries Only',GotIt ) 
  SaveBulkOnly = GetLogical( Params,'Save Bulk Only',GotIt ) 
  
  NumberOfNodes = Mesh % NumberOfNodes
  IF( MaskExists ) THEN
    NumberOfNodes = COUNT( MaskPerm(1:NumberOfNodes) > 0 ) 
  END IF
  NumberOfElements = 0

  IF( NumberOfNodes > 0) THEN
    ElemFirst = HUGE( ElemFirst )
    ElemLast = 0 
    ALLOCATE(NodePerm(Mesh % NumberOfNodes))
    NodePerm = 0

    ALLOCATE(ActiveElem(Mesh % NumberOfBulkElements + & 
        Mesh % NumberOfBoundaryElements))
    ActiveElem = .FALSE.
    
    ! Count the true number of elements and mark the 1st and last element
    !-----------------------------------------------------------------------
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      
      IF( i > Mesh % NumberOfBulkElements ) THEN
        IF( SaveBulkOnly ) CYCLE
      ELSE
        IF( SaveBoundariesOnly ) CYCLE
      END IF
      
      CurrentElement => Mesh % Elements(i)
      Model % CurrentElement => CurrentElement

      IF( GetElementFamily( CurrentElement ) == 1 ) CYCLE          
      
      IF( MaskExists ) THEN
        IF( ANY(MaskPerm(CurrentElement % NodeIndexes) <= 0) ) CYCLE
      END IF

      IF( GotMaskName ) THEN
        Hit = .FALSE.
        IF( i <= Mesh % NumberOfBulkElements ) THEN
          l = CurrentElement % BodyId
          k = ListGetInteger( Model % Bodies(l) % Values,'Body Force',GotIt)
          IF( GotIt ) THEN
            Hit = ListGetLogical( Model % BodyForces(k) % Values, TRIM(MaskName), GotIt)
          END  IF
          IF( .NOT. Hit ) THEN
            k = ListGetInteger( Model % Bodies(l) % Values,'Equation',GotIt)
            IF( GotIt ) THEN
              Hit = ListGetLogical( Model % Equations(k) % Values, TRIM(MaskName), GotIt)
            END IF
          END IF
        ELSE
          DO l=1, Model % NumberOfBCs
            IF ( Model % BCs(l) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
            Hit = ListGetLogical(Model % BCs(l) % Values, MaskName, GotIt ) 
            EXIT
          END DO
        END IF
        IF(.NOT. Hit ) CYCLE
      END IF

      IF( GotMaskCond ) THEN
        n = CurrentElement % TYPE % NumberOfNodes
        Indexes => CurrentElement % NodeIndexes

        IF( i <= Mesh % NumberOfBulkElements ) THEN
          l = CurrentElement % BodyId
          k = ListGetInteger( Model % Bodies(l) % Values,'Body Force',GotIt)
          IF( GotIt ) THEN
            MaskCond(1:n) = ListGetReal( Model % BodyForces(k) % Values, TRIM(MaskName), &
                n, Indexes, GotIt)
          END  IF

          IF( .NOT. Hit ) THEN
            k = ListGetInteger( Model % Bodies(l) % Values,'Equation',GotIt)
            IF( GotIt ) THEN
              MaskCond(1:n) = ListGetReal( Model % Equations(k) % Values, TRIM(MaskName), &
                  n, Indexes, GotIt)
            END IF
          END IF
        ELSE
          GotIt = .FALSE.
          IF( ASSOCIATED( CurrentElement % BoundaryInfo ) ) THEN
            DO l=1, Model % NumberOfBCs
              IF ( Model % BCs(l) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
              MaskCond(1:n) = ListGetReal(Model % BCs(l) % Values, MaskName, &
                  n, Indexes, GotIt ) 
              EXIT
            END DO
          END IF
        END IF
        IF( .NOT. GotIt ) CYCLE
        IF( .NOT. ALL(MaskCond(1:n) > 0.0_dp ) ) CYCLE
      END IF

      NodePerm( CurrentElement % NodeIndexes ) = 1
      ActiveElem(i) = .TRUE.

      NumberOfElements = NumberOfElements + 1
      ElemFirst = MIN( ElemFirst, i )
      ElemLast = MAX( ElemLast, i )
    END DO
    

    ! Make the numbered list of nodes used. 
    NumberOfNodes = COUNT( NodePerm > 0 ) 
    AllNodesUsed = ( NumberOfNodes == Mesh % NumberOfNodes )

    IF( AllNodesUsed ) THEN
      DEALLOCATE( NodePerm ) 
    ELSE
      ALLOCATE( NodeOrder( NumberOfNodes ) ) 
      NodeOrder = 0
      j = 0
      DO i=1,Mesh % NumberOfNodes
        IF( NodePerm(i) > 0 ) THEN
          j = j + 1       
          NodePerm(i) = j
          NodeOrder(j) = i
        END IF
      END DO
    END IF
  END IF

  ! The partition is active for saving if there are any nodes 
  ! to write. There can be no elements nor dofs without nodes.
  CALL ParallelActive( NumberOfNodes > 0 )

  IF( nTime == 1 ) THEN
    ParallelNodes = NINT( ParallelReduction( 1.0_dp * NumberOfNodes ) )
    ParallelElements = NINT( ParallelReduction( 1.0_dp * NumberOfElements ) )
    WRITE( Message,'(A,I8)') 'Number of nodes to save:',ParallelNodes
    CALL Info('VtuOutputSolver',Message)
    WRITE( Message,'(A,I8)') 'Number of elements to save:',ParallelElements
    CALL Info('VtuOutputSolver',Message)
  END IF

  IF( BinaryOutput ) THEN
    BufferSize = GetInteger( Params,'Binary Output Buffer Size',GotIt)
    IF( .NOT. GotIt ) BufferSize = MAX( NumberOfNodes, NumberOfElements )
  END IF

  BaseFile = FilePrefix
  IF ( .NOT. FileNameQualified(FilePrefix) ) THEN
    Dir = GetString( Params,'Output Directory',GotIt) 
    IF( GotIt ) THEN
      IF( LEN_TRIM(Dir) > 0 ) THEN
        BaseFile = TRIM(Dir)// '/' //TRIM(FilePrefix)
        CALL MakeDirectory( TRIM(Dir) // CHAR(0) )
      END IF
    ELSE 
      BaseFile = TRIM(OutputPath) // '/' // TRIM(Mesh % Name) // '/' //TRIM(FilePrefix)
    END IF
  END IF
  CALL Info('VtuOutputSolver','Full filename base is: '//TRIM(Basefile), Level=10 )


  ActiveModes => ListGetIntegerArray( Params,'Active EigenModes',GotActiveModes ) 
  IF( GotActiveModes ) THEN
    MaxModes = SIZE( ActiveModes )
  ELSE
    MaxModes = GetInteger( Params,'Number of EigenModes',GotIt)
    IF(.NOT. GotIt) MaxModes = GetInteger( Params,'Eigen System Values',GotIt)
    IF(.NOT. GotIt) THEN
      DO i=1,Model % NumberOfSolvers
        MaxModes = MAX( MaxModes, &
            GetInteger( Model % Solvers(i) % Values,'Eigen System Values', GotIt ) )
        MaxModes = MAX( MaxModes, &
            GetInteger( Model % Solvers(i) % Values,'Harmonic System Values', GotIt ) )       
      END DO
    END IF     
  END IF

  EigenAnalysis = ListGetLogical( Params,'Eigen Analysis',GotIt)
  IF( EigenAnalysis ) THEN
    CALL Info('VtuOutputSolver','Saving each eigenmode to different file')
    FileIndex = 1
  END IF

  WriteIds = GetLogical( Params,'Save Geometry Ids',GotIt)  
  IF( ElemFirst <= Mesh % NumberOfBulkElements ) THEN
    BCOffset = 100
    DO WHILE( BCOffset <= Model % NumberOfBodies ) 
      BCOffset = 10 * BCOffset
    END DO
  ELSE
    BCOffset = 0
  END IF


 100   CONTINUE

  IF(Parallel) THEN
    IF( NoFileindex ) THEN
      WRITE( PvtuFile,'(A,".pvtu")' ) TRIM(BaseFile)
    ELSE IF( FileIndex < 10000 ) THEN
      WRITE( PvtuFile,'(A,I4.4,".pvtu")' ) TRIM(BaseFile),FileIndex
    ELSE   
      WRITE( PvtuFile,'(A,I0,".pvtu")' ) TRIM(BaseFile),FileIndex
    END IF
    CALL Info('VtuOutputSolver','Writing the pvtu file: '//TRIM(PvtuFile), Level=10)
    CALL WritePvtuFile( PvtuFile, Model )
    CALL Info('VtuOutputSolver','Finished writing pvtu file',Level=12)
  END IF


  ! Write the Vtu file with all the data
  !--------------------------------------------------------------------------
  IF( NumberOfNodes > 0 ) THEN
    IF ( Parallel ) THEN
      IF( NoFileindex ) THEN
        WRITE( VtuFile,'(A,I4.4,A,".vtu")' ) TRIM(BaseFile),Part+1,"par"
      ELSE IF( FileIndex < 10000 ) THEN
        WRITE( VtuFile,'(A,I4.4,A,I4.4,".vtu")' ) TRIM(BaseFile),Part+1,"par",&
            FileIndex
      ELSE
        WRITE( VtuFile,'(A,I4.4,A,I0,".vtu")' ) TRIM(BaseFile),Part+1,"par",&
            FileIndex
      END IF
    ELSE
      IF( NoFileindex ) THEN
        WRITE( VtuFile,'(A,".vtu")' ) TRIM(BaseFile)
      ELSE IF( FileIndex < 10000 ) THEN
        WRITE( VtuFile,'(A,I4.4,".vtu")' ) TRIM(BaseFile),FileIndex
      ELSE
        WRITE( VtuFile,'(A,I0,".vtu")' ) TRIM(BaseFile),FileIndex
      END IF
    END IF

    CALL Info('VtuOutputSolver','Writing the vtu file: '//TRIM(VtuFile),Level=7)
    CALL WriteVtuFile( VtuFile, Model, FixedMesh )
    CALL Info('VtuOutputSolver','Finished writing vtu file',Level=12)
  END IF

  ! For transient simulation write a holder for the timesteps
  !-----------------------------------------------------------
  IF( GetLogical( Params,'Vtu Time Collection', GotIt ) ) THEN
    IF( TransientSimulation .AND. .NOT. NoFileIndex ) THEN
      WRITE( PvdFile,'(A,".pvd")' ) TRIM(BaseFile)
      IF( Parallel ) THEN
        DataSetFile = PvtuFile
      ELSE
        DataSetFile = VtuFile
      END IF
      CALL Info('VtuOutputSolver','Writing the pvd file: '//TRIM(DataSetFile),Level=10)
      CALL WritePvdFile( PvdFile, DataSetFile, FileIndex, Model )
      CALL Info('VtuOutputSolver','Finished writing pvd file',Level=12)     
    END IF
  END IF


  IF( EigenAnalysis ) THEN
    FileIndex = FileIndex + 1
    IF( FileIndex <= MaxModes ) GOTO 100
  END IF

  IF( NumberOfNodes > 0 ) THEN
    IF( .NOT. AllNodesUsed ) THEN
      DEALLOCATE( NodeOrder, NodePerm ) 
    END IF
    DEALLOCATE( ActiveElem ) 
  END IF

  CALL Info('VtuOutputSolver','All done for now',Level=10)     


CONTAINS




  SUBROUTINE WriteVtuFile( VtuFile, Model, RemoveDisp )
!!!!  
!!!!
!!!! HCK TO GET 827 H27 VTU CODE 29  
!!!!
!!!!
    CHARACTER(LEN=*), INTENT(IN) :: VtuFile
    TYPE(Model_t) :: Model 
    LOGICAL, INTENT(IN) :: RemoveDisp
    INTEGER, PARAMETER :: VtuUnit = 58
    TYPE(Variable_t), POINTER :: Var,Var1
    CHARACTER(LEN=512) :: str
    INTEGER :: i,ii,j,k,dofs,Rank,cumn,n,dim,vari,sdofs,dispdofs, disp2dofs, Offset, &
        NoFields, IndField, iField, NoModes
    CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, TensorFieldName, &
        FieldName, FieldName2, OutStr
    CHARACTER :: lf
    LOGICAL :: ScalarsExist, VectorsExist, Found,&
              ComponentVector, ComplementExists, Use2
    LOGICAL :: WriteData, WriteXML, L, Buffered
    TYPE(Variable_t), POINTER :: Solution
    INTEGER, POINTER :: Perm(:), Perm2(:), DispPerm(:), Disp2Perm(:)
    REAL(KIND=dp), POINTER :: Values(:), DispValues(:), Disp2Values(:), Values2(:), Values3(:)
    REAL(KIND=dp) :: x,y,z, val
    INTEGER(KIND=dp), ALLOCATABLE :: ElemInd(:)
    REAL(KIND=dp), ALLOCATABLE :: ElemVectVal(:)
    COMPLEX(KIND=dp), POINTER :: EigenVectors(:,:)
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Element_t), POINTER :: CurrentElement, Parent
    TYPE(ValueList_t), POINTER :: Params
    !--------------------------------------------------------------
    !--------------------------------------------------------------    
    INTEGER, ALLOCATABLE :: BigIndex(:),Newind(:)
    !--------------------------------------------------------------
    !--------------------------------------------------------------


    ! Initialize the auxiliary module for buffered writing
    !--------------------------------------------------------------
    CALL AscBinWriteInit( AsciiOutput, SinglePrec, VtuUnit, BufferSize )
    
    ! Linefeed character
    !-----------------------------------
    lf = CHAR(10)

    Offset = 0
    WriteXML = .TRUE.
    WriteData = AsciiOutput
    Params => GetSolverParams()
    Buffered = .TRUE.

    
    ! This is a hack to ensure that the streamed saving will cover the whole file
    !----------------------------------------------------------------------------
    IF(.TRUE.) THEN
      OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'formatted', STATUS='unknown' )
      WRITE( VtuUnit,'(A)') ' '
      CLOSE( VtuUnit ) 
    END IF


    ! This format works both for ascii and binary output
    !-------------------------------------------------------------------------
    OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'unformatted', ACCESS = 'stream', STATUS='unknown' )
    
    Solver => Model % Solver

    ! VTU seemingly only works with 3D cases, so enforce it
    dim = 3

    
    WRITE( OutStr,'(A)') '<?xml version="1.0"?>'//lf
    CALL AscBinStrWrite( OutStr ) 

    IF ( LittleEndian() ) THEN
      OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
    ELSE
      OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'//lf
    END IF
    CALL AscBinStrWrite( OutStr )
    WRITE( OutStr,'(A)') '  <UnstructuredGrid>'//lf
    CALL AscBinStrWrite( OutStr )
    WRITE( OutStr,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NumberOfNodes,&
        '" NumberOfCells="',NumberOfElements,'">'//lf
    CALL AscBinStrWrite( OutStr )
    
    ! nodewise information
    !-------------------------------------
    ScalarFieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
    VectorFieldName = GetString( Params,'Vector Field 1',VectorsExist)
    IF( .NOT. ( ScalarsExist .OR. VectorsExist) ) THEN
      CALL Warn('WriteVtuFile','Are there really no scalars or vectors?')
    END IF
    
    WRITE( OutStr,'(A)') '      <PointData>'//lf
    CALL AscBinStrWrite( OutStr )

    DispDofs = 0
    Disp2Dofs = 0
    IF(RemoveDisp) THEN
      Solution => VariableGet( Model % Mesh % Variables, 'Displacement')
      IF( ASSOCIATED( Solution ) ) THEN
        Solver => Solution % Solver
        L = GetLogical( GetSolverParams(Solver),'Displace Mesh',Found)
        IF(.NOT.Found) L=.NOT.EigenOrHarmonicAnalysis(Solver)
        IF (L) THEN
          DispPerm => Solution % Perm
          DispValues => Solution % Values
          DispDofs = Solution % Dofs
        END IF
      END IF
      
      Solution => VariableGet( Model % Mesh % Variables, 'Mesh Update')
      IF( ASSOCIATED( Solution ) ) THEN
        Disp2Perm => Solution % Perm
        Disp2Values => Solution % Values
        Disp2Dofs = Solution % Dofs
      END IF
    END IF

    
! ! ! ! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*         
! ! ! ! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*         
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !     
! ! !      ALLOCATE( BigIndex(NumberOfNodes) )
! ! ! ! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*         
! ! ! ! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*         
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                 
! ! !      ! LOOP over elements in order to create the correct permutation due to 820 or 829
! ! !      DO i = ElemFirst, ElemLast
! ! !         IF(.NOT. ActiveElem(i) ) CYCLE
! ! ! 
! ! !         CurrentElement => Model % Elements(i)
! ! !         n = GetElementNOFNodes(CurrentElement)
! ! !        
! ! !         ALLOCATE(Newind(n)) 
! ! ! !    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -
! ! !         IF (Elmer2VTKElement(CurrentElement % TYPE % ElementCode) == 29) THEN
! ! !         ! Elmer 827 -> VTU 29        
! ! !             Newind(1:n) = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16,24,22,21,23,25,26,27/)
! ! !         ElseIf (Elmer2VTKElement(CurrentElement % TYPE % ElementCode) == 25) THEN
! ! !         ! 820 -> 25
! ! !             Newind(1:n) = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16/)
! ! !         ELSE
! ! !             Newind(1:n) = (/ (i, i = 1, n) /)   ! nothing change
! ! !         ENDIF
! ! ! !    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -        
! ! !         DO j=1,n
! ! !            BigIndex( CurrentElement % NodeIndexes(j) ) =  CurrentElement % NodeIndexes(Newind(j))
! ! !         END DO
! ! ! !    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -        
! ! !         DEALLOCATE(Newind)            
! ! !       END DO
! ! !      
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !          
! ! ! ! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*         
! ! ! ! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 

    ! When the data is 'appended' two loops will be taken and the data will be written
    ! on the second loop. Offset is the position in the appended data after the '_' mark.
    !------------------------------------------------------------------------------------
100 Offset = 0
    
  IF( SaveNodal ) THEN
    DO Rank = 0,2
      DO Vari = 1, 999
        IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
        IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari
        
        FieldName = GetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        IF(Rank == 2) THEN
          CALL Fatal('VtuOutputSolver','Do the tensors')
        END IF
        
        !---------------------------------------------------------------------
        ! Find the variable with the given name in the normal manner 
        !---------------------------------------------------------------------
        Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
        ComponentVector = .FALSE.
        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
          IF( ASSOCIATED(Solution)) THEN 
            ComponentVector = .TRUE.
          ELSE
            WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
            CALL Warn('WriteVtuXMLFile', Txt)
            CYCLE
          END IF
        END IF
        IF ( Solution % TYPE == Variable_on_nodes_on_elements ) CYCLE
        
        IF( ASSOCIATED(Solution % EigenVectors) ) THEN
          NoModes = SIZE( Solution % EigenValues )
          IF( ComponentVector ) THEN
            CALL Warn('WriteVtuXMLFile','Eigenmodes cannot be given componentwise!')
            CYCLE
          ELSE IF( EigenAnalysis ) THEN
            IF( GotActiveModes ) THEN
              IndField = ActiveModes( FileIndex ) 
            ELSE
              IndField = FileIndex
            END IF
            IF( IndField > NoModes ) THEN
              WRITE( Message,'(A,I0,A,I0,A)') 'Too few eigenmodes (',&
                     IndField,',',NoModes,') in '//TRIM(FieldName)       
              CALL Warn('WriteVtuXMLFile',Message)
              CYCLE
            END IF
            NoModes = 1
            NoFields = 1
          ELSE    
            IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
            NoFields = NoModes
          END IF
          EigenVectors => Solution % EigenVectors
        ELSE
          NoModes = 0 
          NoFields = 1
        END IF
        
        Perm => Solution % Perm
        dofs = Solution % DOFs
        Values => Solution % Values

        !---------------------------------------------------------------------
        ! Some vectors are defined by a set of components (either 2 or 3)
        !---------------------------------------------------------------------
        IF( ComponentVector ) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
          IF( ASSOCIATED(Solution)) THEN
            Values2 => Solution % Values
            dofs = 2
          END IF
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
          IF( ASSOCIATED(Solution)) THEN
            Values3 => Solution % Values
            dofs = 3
          END IF
        END IF

        !---------------------------------------------------------------------
        ! There may be special complementary variables such as 
        ! displacement & mesh update 
        !---------------------------------------------------------------------
        ComplementExists = .FALSE.
        IF( NoModes > 0 ) THEN
          IF(Rank==0) WRITE(Txt,'(A,I0,A)') 'Scalar Field ',Vari,' Complement'
          IF(Rank==1) WRITE(Txt,'(A,I0,A)') 'Vector Field ',Vari,' Complement'
          IF(Rank==2) WRITE(Txt,'(A,I0,A)') 'Tensor Field ',Vari,' Complement'
          
          FieldName2 = GetString( Params, TRIM(Txt), Found )
          IF( Found ) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName2))
            IF( ASSOCIATED(Solution)) THEN 
              Values2 => Solution % Values
              Perm2 => Solution % Perm 
              ComplementExists = .TRUE.
            ELSE
              CALL Warn('WriteVTUFile','Complement does not exist:'//TRIM(FieldName2))
            END IF
          END IF
        END IF

        IF( dofs > 1 ) THEN
          sdofs = MAX(dofs,dim)
        ELSE
          sdofs = 1
        END IF
        

        !---------------------------------------------------------------------
        ! Finally save the field values 
        !---------------------------------------------------------------------
        DO iField = 1, NoFields          

          IF( Nomodes > 0 .AND. .NOT. EigenAnalysis ) THEN
            IF( GotActiveModes ) THEN
              IndField = ActiveModes( iField ) 
            ELSE
              IndField = iField
            END IF
          END IF

          IF( WriteXML ) THEN
            IF( NoModes == 0 .OR. EigenAnalysis ) THEN
              WRITE( OutStr,'(A,I0,A)') '        <DataArray type="Float',PrecBits,'" Name="'//TRIM(FieldName)
            ELSE 
              WRITE( OutStr,'(A,I0,A,I0)') '        <DataArray type="Float',PrecBits,'" Name="'//&
                  TRIM(FieldName)//' mode',IndField
            END IF
            CALL AscBinStrWrite( OutStr )
            
            WRITE( OutStr,'(A,I0,A)') '" NumberOfComponents="',sdofs,'"'          
            CALL AscBinStrWrite( OutStr ) 
            
            IF( AsciiOutput ) THEN
              WRITE( OutStr,'(A)') ' format="ascii">'//lf
              CALL AscBinStrWrite( OutStr ) 
            ELSE
              WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
              CALL AscBinStrWrite( OutStr ) 
            END IF
          END IF
          
          IF( BinaryOutput ) THEN
            k = NumberOfNodes * PrecSize * sdofs
            Offset = Offset + IntSize + k
          END IF
         
   
          !---------------------------------------------------------------------
          ! Data may also be appended and then its saved on the second sweep
          !---------------------------------------------------------------------
          IF( WriteData ) THEN
            
            IF( BinaryOutput ) WRITE( VtuUnit ) k
            
            DO ii = 1, NumberOfNodes
              IF( AllNodesUsed ) THEN
                i = ii 
              ELSE
                i = NodeOrder(ii)
              END IF

              IF( ASSOCIATED( Perm ) ) THEN
                j = Perm(i)
              ELSE
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*                      
                j = i
!                 j = BigIndex(i)
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*        
              END IF
              
              Use2 = .FALSE.
              IF( ComplementExists ) THEN
                IF( j == 0 ) THEN
                  Use2 = .TRUE. 
                  j = Perm2(i)
                END IF
              END IF
              
              DO k=1,sdofs              
                IF(j==0 .OR. k > dofs) THEN
                  val = 0.0_dp
                ELSE IF( ComponentVector ) THEN
                  IF( k == 1 ) val = Values(j)
                  IF( k == 2 ) val = Values2(j)
                  IF( k == 3 ) val = Values3(j)
                ELSE IF( Use2 ) THEN
                  val = Values2(dofs*(j-1)+k)              
                ELSE IF( NoModes > 0 ) THEN
                  val = EigenVectors(IndField,dofs*(j-1)+k)                              
                ELSE
                  val = Values(dofs*(j-1)+k)              
                END IF
                
                CALL AscBinRealWrite( val )

              END DO
            END DO

            CALL AscBinRealWrite( 0.0_dp, .TRUE. )

          END IF
        
          IF( AsciiOutput ) THEN
            WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
            CALL AscBinStrWrite( OutStr ) 
          END IF
        END DO
      END DO
    END DO
  END IF ! IF( SaveNodal )
  
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      </PointData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    
    ! Elementwise information
    !-------------------------------------
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      <CellData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    
  IF( SaveElemental ) THEN
    DO Rank = 0,1
      DO Vari = 1, 999

        IF( Rank == 0 ) THEN
          WRITE(Txt,'(A,I0)') 'Scalar Field Elemental ',Vari
        ELSE
          WRITE(Txt,'(A,I0)') 'Vector Field Elemental ',Vari
        END IF
        FieldName = GetString( Params, TRIM(Txt), Found )
        L = Found

        IF(.NOT. Found) THEN
          IF( Rank == 0 ) THEN
            WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
          ELSE
            WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
          END IF
          FieldName = GetString( Params, TRIM(Txt), Found )
        END IF

        IF(.NOT. Found) EXIT

        !---------------------------------------------------------------------
        ! Find the variable with the given name in the normal manner 
        !---------------------------------------------------------------------
        Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
        ComponentVector = .FALSE.

        ! If we are looking for a vector just one dofs wont do!
        ! This circumvents a problem somewhere else in the code. 
        IF( ASSOCIATED( Solution ) ) THEN
          IF( Rank > 0 .AND. Solution % Dofs <= 1 ) NULLIFY( Solution ) 
        END IF

        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
          IF( ASSOCIATED(Solution)) THEN 
            ComponentVector = .TRUE.
          ELSE 
            IF( L ) THEN
              WRITE(Txt, '(A,A)') 'Nonexistent elemental variable: ',TRIM(FieldName)
              CALL Warn('WriteVtuXMLFile', Txt)
            END IF
            CYCLE
          END IF
        END IF

        IF (Solution % TYPE /= Variable_on_nodes_on_elements ) CYCLE

        Perm => Solution % Perm
        Dofs = Solution % DOFs
        Values => Solution % Values

       !---------------------------------------------------------------------
        ! Some vectors are defined by a set of components (either 2 or 3)
        !---------------------------------------------------------------------
        IF( ComponentVector ) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
          IF( ASSOCIATED(Solution)) THEN
            Values2 => Solution % Values
            dofs = 2
          END IF
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
          IF( ASSOCIATED(Solution)) THEN
            Values3 => Solution % Values
            dofs = 3
          END IF
        END IF

        IF( dofs > 1 ) THEN
          sdofs = MAX(dofs,dim)
        ELSE
          sdofs = 1
        END IF

        !---------------------------------------------------------------------
        ! Finally save the field values 
        !---------------------------------------------------------------------
        IF( WriteXML ) THEN
          CALL Info('WriteVtuFile','Writing variable: '//TRIM(FieldName) )
          WRITE( OutStr,'(A,I0,A)') '        <DataArray type="Float',PrecBits,'" Name="'//TRIM(FieldName)
          CALL AscBinStrWrite( OutStr )
          
          WRITE( OutStr,'(A,I0,A)') '" NumberOfComponents="',sdofs,'"'          
          CALL AscBinStrWrite( OutStr ) 
          
          IF( AsciiOutput ) THEN
            WRITE( OutStr,'(A)') ' format="ascii">'//lf
            CALL AscBinStrWrite( OutStr ) 
          ELSE
            WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
            CALL AscBinStrWrite( OutStr ) 
          END IF
        END IF
        
        IF( BinaryOutput ) THEN
          k = PrecSize * sdofs * NumberOfElements
          Offset = Offset + IntSize + k
        END IF

        
        IF( WriteData ) THEN
          IF( BinaryOutput ) WRITE( VtuUnit ) k

          ALLOCATE(ElemVectVal(sdofs), ElemInd(Model % Mesh % MaxElementDOFS))
          DO i = ElemFirst, ElemLast
            IF( .NOT. ActiveElem(i) ) CYCLE
            CurrentElement => Model % Elements(i)

            ElemVectVal = 0._dp; ElemInd = 0
            n = GetElementNOFNodes(CurrentElement)

            IF ( ASSOCIATED(CurrentElement % BoundaryInfo) .AND. .NOT. &
                 ASSOCIATED(CurrentElement % DGIndexes) ) THEN

              Parent => CurrentElement % BoundaryInfo % Left
              IF (.NOT.ASSOCIATED(Parent) ) &
                Parent => CurrentElement % BoundaryInfo % Right

              IF ( ASSOCIATED(Parent) ) THEN
                IF (ASSOCIATED(Parent % DGIndexes) ) THEN
                  DO j=1,n
                    DO k=1,Parent % TYPE % NumberOfNodes
                      IF(Currentelement % NodeIndexes(j) == Parent % NodeIndexes(k)) &
                      ElemInd(j) = Perm( Parent % DGIndexes(k) )
                    END DO
                  END DO
                END IF
              END IF
            ELSE
              ElemInd(1:n) = Perm( CurrentElement % DGIndexes(1:n) )
            END IF

            IF ( ALL(ElemInd(1:n)>0)) THEN
              DO k=1,sdofs
                IF( k > dofs ) THEN
                  ElemVectVal(k) = 0.0_dp
                ELSE IF(ComponentVector) THEN
                  IF (k==1) ElemVectVal(k) = SUM(Values(ElemInd(1:n)))/n
                  IF (k==2) ElemVectVal(k) = SUM(Values2(ElemInd(1:n)))/n
                  IF (k==3) ElemVectVal(k) = SUM(Values3(ElemInd(1:n)))/n
                ELSE
                  ElemVectVal(k) = SUM(Values(dofs*(ElemInd(1:n)-1)+k))/n
                END IF
              END DO
            END IF

            DO k=1,sdofs
              CALL AscBinRealWrite( ElemVectVal(k) )
            END DO
          END DO

          CALL AscBinRealWrite( 0.0_dp, .TRUE. )
          DEALLOCATE(ElemVectVal,ElemInd)
        END IF
        
        IF( AsciiOutput ) THEN
          WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END DO
    END DO
  END IF  ! IF( SaveElemental )

    !---------------------------------------------------------------------
    ! If requested write the body and bc indexes
    !---------------------------------------------------------------------
    IF( WriteIds ) THEN
      !---------------------------------------------------------------------
      ! Finally save the field values 
      !---------------------------------------------------------------------
      IF( WriteXML ) THEN
        CALL Info('WriteVtuFile','Writing body and BC indexes')

        WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="GeometryIds"'
        CALL AscBinStrWrite( OutStr )
          
        IF( AsciiOutput ) THEN
          WRITE( OutStr,'(A)') ' format="ascii">'//lf
          CALL AscBinStrWrite( OutStr ) 
        ELSE
          WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF
        
      IF( BinaryOutput ) THEN
        k = IntSize * NumberOfElements
        Offset = Offset + IntSize + k
      END IF
              
      IF( WriteData ) THEN        
        IF( BinaryOutput ) WRITE( VtuUnit ) k
        
        DO i = ElemFirst, ElemLast
          IF(.NOT. ActiveElem(i)) CYCLE

          CurrentElement => Model % Elements(i)
          
          IF( i <= Mesh % NumberOfBulkElements ) THEN
            j = CurrentElement % BodyId
          ELSE
            j = GetBCId( CurrentElement ) + BCOffset
          END IF
          
          CALL AscBinIntegerWrite( j )
        END DO        
        CALL AscBinIntegerWrite( 0, .TRUE. )
      END IF
        
      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      </CellData>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF        

    ! Coordinates of each point
    !-------------------------------------
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      <Points>'//lf
      CALL AscBinStrWrite( OutStr ) 

      WRITE( OutStr,'(A,I0,A,I0,A)') '        <DataArray type="Float',PrecBits,'" NumberOfComponents="',dim,'"'
      CALL AscBinStrWrite( OutStr )       

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' format="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF


    IF( BinaryOutput ) THEN
      k = dim * NumberOfNodes * PrecSize
      Offset = Offset + IntSize + k
    END IF

    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k 

      DO ii = 1, NumberOfNodes
        IF( AllNodesUsed ) THEN
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*             
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*             
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*             
          i = ii 
!             i = BigIndex(ii)
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*     
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*     
! ! ! ! !-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*     
        ELSE
          i = NodeOrder(ii)
        END IF

        x = Model % Mesh % Nodes % x(i)
        y = Model % Mesh % Nodes % y(i)
        z = Model % Mesh % Nodes % z(i)
        
        ! If displacement field is active remove the displacement from the coordinates
        IF( dispdofs > 0 .OR. disp2dofs > 0) THEN
          j = 0
          IF(dispdofs > 0) THEN
            j = DispPerm(i)
            IF( j > 0 ) THEN
              x = x - DispValues(dispdofs*(j-1)+1)
              y = y - DispValues(dispdofs*(j-1)+2)
              IF(dispdofs == 3) z = z - DispValues(dispdofs*(j-1)+3)
            END IF
          END IF
          IF(disp2dofs > 0 .AND. j==0) THEN
            j = Disp2Perm(i)
            IF( j > 0 ) THEN
              x = x - Disp2Values(disp2dofs*(j-1)+1)
              y = y - Disp2Values(disp2dofs*(j-1)+2)
              IF(disp2dofs == 3) z = z - Disp2Values(disp2dofs*(j-1)+3)
            END IF
          END IF
        END IF
        
        CALL AscBinRealWrite( x )
        CALL AscBinRealWrite( y )
        CALL AscBinRealWrite( z )
      END DO

      CALL AscBinRealWrite( 0.0_dp, .TRUE.)
    END IF

    IF( AsciiOutput ) THEN   
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      </Points>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    ! Write out the mesh
    !-------------------------------------
! ! !             <DataArray type="Int32" Name="connectivity" format="ascii">
! 0 8 44 17 32 53 98 80 9 45 46 19 33 54 99 81 55 100 101 84 56 102 103 85 47 104 105 32 53 98 80 4 20 89 29 55 100 101 84 34 59
!         </DataArray>

    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '      <Cells>'//lf
      CALL AscBinStrWrite( OutStr ) 
      
      WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="connectivity"'
      CALL AscBinStrWrite( OutStr ) 

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' format="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( BinaryOutput ) THEN
      ! The offset needs to be summed over all elements, this is just the size factor
      k = 0
      DO i = ElemFirst, ElemLast
        IF( .NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
        n = GetElementNOFNodes(CurrentElement)
    
        k = k + n * IntSize
      END DO
      Offset = Offset + k + IntSize
    END IF


    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k

      DO i = ElemFirst, ElemLast
        IF(.NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
        n = GetElementNOFNodes(CurrentElement)
       
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                 
        ALLOCATE(Newind(n))            
        IF (Elmer2VTKElement(CurrentElement % TYPE % ElementCode) == 29) THEN
        ! Elmer 827 -> VTU 29        
            Newind(1:n) = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16,24,22,21,23,25,26,27/)
   
        ElseIf (Elmer2VTKElement(CurrentElement % TYPE % ElementCode) == 25) THEN
        ! 820 -> 25
            Newind(1:n) = (/1,2,3,4,5,6,7,8,9,10,11,12,17,18,19,20,13,14,15,16/)
        ELSE
            Newind(1:n) = (/ (i, i = 1, n) /)   ! nothing change
        ENDIF            

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                        
        DO j=1,n
          IF( AllNodesUsed ) THEN
! !             CALL AscBinIntegerWrite( CurrentElement % NodeIndexes(j) - 1 )
            CALL AscBinIntegerWrite( CurrentElement % NodeIndexes(Newind(j)) - 1 )
          ELSE
            CALL AscBinIntegerWrite( NodePerm( CurrentElement % NodeIndexes(j)) - 1 )
          END IF
        END DO
        DEALLOCATE(Newind)            

      END DO
      CALL AscBinIntegerWrite( 0, .TRUE. ) 
    END IF

    IF( AsciiOutput ) THEN
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF


    ! Offsets for element indexes 
    !-------------------------------------------------------------------

    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="offsets"'
      CALL AscBinStrWrite( OutStr ) 
      
      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' format="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( BinaryOutput ) THEN
      k = NumberOfElements * IntSize
      Offset = Offset + IntSize + k
    END IF


    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k 

      cumn = 0
      DO i = ElemFirst, ElemLast
        IF( .NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
        n = CurrentElement % TYPE % NumberOfNodes
        cumn = cumn + n

        CALL AscBinIntegerWrite( cumn )
      END DO
        
      CALL AscBinIntegerWrite( 0, .TRUE.)
      
    END IF


    IF( AsciiOutput ) THEN   
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    IF( WriteXML ) THEN
      WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="types"'
      CALL AscBinStrWrite( OutStr ) 

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') ' FORMAT="ascii">'//lf
        CALL AscBinStrWrite( OutStr ) 
      ELSE
        WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    IF( BinaryOutput ) THEN
      k = NumberOfElements * IntSize
      Offset = Offset + IntSize + k
    END IF


    IF( WriteData ) THEN
      IF( BinaryOutput ) WRITE( VtuUnit ) k

      DO i = ElemFirst, ElemLast
        IF( .NOT. ActiveElem(i) ) CYCLE

        CurrentElement => Model % Elements(i)
! !         WRITE (*,*) "CurrentElement % TYPE % ElementCode=",CurrentElement % TYPE % ElementCode
        n = Elmer2VTKElement(CurrentElement % TYPE % ElementCode)

        CALL AscBinIntegerWrite( n )
      END DO

      CALL AscBinIntegerWrite( 0, .TRUE. )     
    END IF

    IF( AsciiOutput ) THEN
      WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF
    IF( WriteXml ) THEN
      WRITE( OutStr,'(A)') '      </Cells>'//lf
      CALL AscBinStrWrite( OutStr ) 
      WRITE( OutStr,'(A)') '    </Piece>'//lf
      CALL AscBinStrWrite( OutStr ) 
      WRITE( OutStr,'(A)') '  </UnstructuredGrid>'//lf
      CALL AscBinStrWrite( OutStr ) 
    END IF

    IF( BinaryOutput ) THEN
      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '<AppendedData encoding="raw">'//lf                    
        CALL AscBinStrWrite( OutStr )           
        WRITE( VtuUnit ) '_'
        
        WriteXML = .FALSE.
        WriteData = .TRUE.
        GOTO 100
      ELSE
        WRITE( OutStr,'(A)') lf//'</AppendedData>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
    END IF

    WRITE( OutStr,'(A)') '</VTKFile>'//lf
    CALL AscBinStrWrite( OutStr ) 

    WRITE( OutStr,'(A)') ' '
    CALL AscBinStrWrite( OutStr ) 

    CLOSE( VtuUnit )

    CALL AscBinWriteFree()

    
! ! ! ! ! !     DEALLOCATE(BigIndex)
  
  END SUBROUTINE WriteVtuFile
    


  SUBROUTINE WritePvdFile( PvdFile, DataSetFile, nTime, Model )
    CHARACTER(LEN=*), INTENT(IN) :: PvdFile, DataSetFile
    INTEGER :: nTime, RecLen = 0
    TYPE(Model_t) :: Model     
    INTEGER, PARAMETER :: VtuUnit = 58
    INTEGER :: n
    REAL(KIND=dp) :: time
    CHARACTER :: lf
    CHARACTER(LEN=MAX_NAME_LEN) :: Str

    SAVE RecLen

    lf = CHAR(10)

    IF( ParEnv % PEs > 1 ) THEN
      IF( ParEnv % MyPE > 0 ) RETURN
    END IF
    time = GetTime()
    IF( GetLogical( Params,'Vtu time previous') ) THEN
      time = time - GetTimestepSize()
    END IF


    IF( nTime == 1 .OR. Reclen == 0 ) THEN
      ! Find the maximum record length (modulo four)
      WRITE( Str,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>'
      n = LEN_TRIM( Str ) 

      WRITE( Str,'(A,ES16.7,A)') '<DataSet timestep="',time,&
        '" group="" part="0" file="'//TRIM(DataSetFile)//'"/>'
      n = MAX( LEN_TRIM( Str ), n ) 
      
      RecLen = ((n/4)+1)*4
    END IF

    IF( nTime == 1 ) THEN
      OPEN( UNIT=VtuUnit, FILE=PvdFile, form = 'formatted', STATUS='REPLACE', &
          ACCESS='DIRECT', ACTION='WRITE', RECL=RecLen)

      IF ( LittleEndian() ) THEN
        WRITE( VtuUnit,'(A)',REC=1) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>'//lf
      ELSE
        WRITE( VtuUnit,'(A)',REC=1) '<VTKFile type="Collection" version="0.1" byte_order="BigEndian"><Collection>'//lf
      END IF     
    ELSE
      OPEN( UNIT=VtuUnit, FILE=PvdFile, form = 'formatted', STATUS='OLD', &
          ACCESS='DIRECT', ACTION='READWRITE', RECL=RecLen)     
    END IF

    WRITE( VtuUnit,'(A,ES12.3,A)',REC=nTime+1) '<DataSet timestep="',time,&
        '" group="" part="0" file="'//TRIM(DataSetFile)//'"/>'//lf
    WRITE( VtuUnit,'(A)',REC=nTime+2) '</Collection></VTKFile>'//lf

    CLOSE( VtuUnit )

  END SUBROUTINE WritePvdFile


  SUBROUTINE WritePvtuFile( VtuFile, Model )
    CHARACTER(LEN=*), INTENT(IN) :: VtuFile
    TYPE(Model_t) :: Model 
    INTEGER, PARAMETER :: VtuUnit = 58
    TYPE(Variable_t), POINTER :: Var,Var1
    CHARACTER(LEN=512) :: str
    INTEGER :: i,j,k,dofs,Rank,cumn,n,dim,vari,sdofs
    CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, TensorFieldName, &
        FieldName, FullName, ShortName
    LOGICAL :: ScalarsExist, VectorsExist, Found, VeloFlag, ComponentVector, &
               AllActive, ThisActive, L
    LOGICAL, POINTER :: ActivePartition(:)
    TYPE(Variable_t), POINTER :: Solution
    INTEGER, POINTER :: Perm(:)
    INTEGER :: Active, NoActive, ierr, NoFields, NoModes, IndField, iField
    REAL(KIND=dp), POINTER :: Values(:)
    COMPLEX(KIND=dp), POINTER :: EigenVectors(:,:)
    TYPE(Element_t), POINTER :: CurrentElement

    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    

    ThisActive = ( NumberOfNodes > 0 ) 
   
    Active = 0
    IF( ThisActive ) Active = 1         

    NoActive = NINT( ParallelReduction( 1.0_dp * Active ) )
    IF( NoActive == 0 ) THEN
      CALL Warn('WritePvtuFile','No active partitions for saving')
    END IF
    AllActive = ( NoActive == Partitions )
    
    WRITE( Message,'(A,I0,A,I0,A)') 'Number of active partitions is ',&
        NoActive,' (out of ',Partitions,')'
    CALL Info('WritePvtuFile',Message,Level=10)

    IF(.NOT. AllActive ) THEN
      IF( Part == 0 ) THEN
        ALLOCATE( ActivePartition(Partitions))
        ActivePartition = .FALSE.
        ActivePartition(1) = ThisActive 
      END IF

      DO i=2,Partitions
        IF( i == Part + 1 ) THEN
          CALL MPI_SEND(ThisActive, 1, MPI_LOGICAL, &
               0, 1000, MPI_COMM_WORLD, ierr )
        END IF
        IF( Part == 0 ) THEN
          CALL MPI_RECV( ActivePartition(i), 1, MPI_LOGICAL, &
               i-1, 1000, MPI_COMM_WORLD, status, ierr )
        END IF
      END DO
    END IF

    IF( Part > 0 ) RETURN
    CALL Info('WritePvtuFile','List of active partitions was composed',Level=12)


    OPEN( UNIT=VtuUnit, FILE=VtuFile, form = 'formatted', STATUS='UNKNOWN' )
    dim = 3

    IF ( LittleEndian() ) THEN
      WRITE( VtuUnit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    ELSE
      WRITE( VtuUnit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="BigEndian">'
    END IF
    WRITE( VtuUnit,'(A)') '  <PUnstructuredGrid>'
    
    ! nodewise information
    !-------------------------------------
    ScalarFieldName = GetString( Params,'Scalar Field 1',ScalarsExist)
    VectorFieldName = GetString( Params,'Vector Field 1',VectorsExist)

  IF( SaveNodal ) THEN
    IF( ScalarsExist .AND. VectorsExist) THEN
      WRITE( VtuUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)&
          //'" Vectors="'//TRIM(VectorFieldName)//'">'
    ELSE IF( ScalarsExist ) THEN
      WRITE( VtuUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)//'">'
    ELSE IF( VectorsExist ) THEN
      WRITE( VtuUnit,'(A)') '    <PPointData Vectors="'//TRIM(VectorFieldName)//'">'
    ELSE
      CALL Warn('WritePvtuFile','Are there really no scalars or vectors?')
    END IF

    
    DO Rank = 0,2
      DO Vari = 1, 999
        IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
        IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        IF(Rank==2) WRITE(Txt,'(A,I0)') 'Tensor Field ',Vari
        
        FieldName = GetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        IF(Rank == 2) THEN
          CALL Warn('WritePvtuFile','Do the tensors')
          EXIT
        END IF
        
        Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
        ComponentVector = .FALSE.

        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
          IF( ASSOCIATED(Solution)) THEN 
            ComponentVector = .TRUE.
          ELSE
            WRITE(Txt, '(A,A)') 'Nonexistent variable 2: ',TRIM(FieldName)
            CALL Warn('WriteVtuXMLFile', Txt)
            CYCLE
          END IF
        END IF

        IF( Solution % TYPE == Variable_on_nodes_on_elements ) CYCLE

        IF( ASSOCIATED(Solution % EigenVectors)) THEN
           NoModes = SIZE( Solution % EigenValues )
           IF( ComponentVector ) THEN
             CALL Warn('WritePvtuXMLFile','Eigenmodes cannot be given componentwise!')
             CYCLE
           ELSE IF( EigenAnalysis ) THEN
             IF( GotActiveModes ) THEN
               IndField = ActiveModes( FileIndex ) 
             ELSE
               IndField = FileIndex
             END IF
             IF( IndField > NoModes ) THEN
               CALL Warn('WriteVtuXMLFile','Too few eigenmodes!')
               CYCLE
             END IF
             NoModes = 1
             NoFields = 1
           ELSE   
             IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
             NoFields = NoModes
           END IF
        ELSE
          NoModes = 0 
          NoFields = 1
        END IF
  
        dofs = Solution % DOFs
        IF( ComponentVector ) THEN
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
          IF( ASSOCIATED(Solution)) dofs = 2
          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
          IF( ASSOCIATED(Solution)) dofs = 3
        END IF
        
        IF( dofs > 1 ) THEN
          sdofs = MAX(dofs,3)
        ELSE
          sdofs = 1
        END IF

        DO iField = 1, NoFields

          IF( NoModes == 0 .OR. EigenAnalysis ) THEN
            FullName = TRIM( FieldName ) 
          ELSE          
            IF( GotActiveModes ) THEN
              IndField = ActiveModes( iField ) 
            ELSE
              IndField = iField
            END IF
            WRITE( FullName,'(A,I0)') TRIM( FieldName )//' mode',IndField
          END IF

          IF( AsciiOutput ) THEN
            WRITE( VtuUnit,'(A,A,I1,A)') '      <PDataArray type="Float64" Name="'//TRIM(FullName), &
                '" NumberOfComponents="',sdofs,'" format="ascii"/>'  
          ELSE 
            WRITE( VtuUnit,'(A,I0,A,A,I0,A)') '      <PDataArray type="Float',PrecBits,'" Name="'//TRIM(FullName), &
                '" NumberOfComponents="',sdofs,'" format="appended"/>'  
          END IF
        END DO

      END DO
    END DO
    WRITE( VtuUnit,'(A)') '    </PPointData>'
  END IF


    ! Elementwise information
    !-------------------------------------
    WRITE( VtuUnit,'(A)') '    <PCellData>'

  IF( SaveElemental ) THEN
    IF( ScalarsExist .OR. VectorsExist ) THEN
      DO Rank = 0,2
        DO Vari = 1, 999

          IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field Elemental ',Vari
          IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field Elemental ',Vari
          FieldName = GetString( Params, TRIM(Txt), Found )
          L = Found 
          
          IF(.NOT. Found) THEN          
            IF(Rank==0) WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
            IF(Rank==1) WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
            FieldName = GetString( Params, TRIM(Txt), Found )
          END IF
          IF(.NOT. Found ) EXIT

          Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName))
          ComponentVector = .FALSE.

          IF(.NOT. ASSOCIATED(Solution)) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 1')
            IF( ASSOCIATED(Solution)) THEN 
              ComponentVector = .TRUE.
            ELSE 
              IF( L ) THEN
                WRITE(Txt, '(A,A)') 'Nonexistent elemental variable 2: ',TRIM(FieldName)
                CALL Warn('WriteVtuXMLFile', Txt)
              END IF
              CYCLE
            END IF
          END IF
          
          IF( Solution % TYPE /= Variable_on_nodes_on_elements ) CYCLE

          IF( ASSOCIATED(Solution % EigenVectors)) THEN
            NoModes = SIZE( Solution % EigenValues )
            IF( ComponentVector ) THEN
              CALL Warn('WritePvtuXMLFile','Eigenmodes cannot be given componentwise!')
              CYCLE
            ELSE IF( EigenAnalysis ) THEN
              IF( GotActiveModes ) THEN
                IndField = ActiveModes( FileIndex ) 
              ELSE
                IndField = FileIndex
              END IF
              IF( IndField > NoModes ) THEN
                CALL Warn('WriteVtuXMLFile','Too few eigenmodes!')
                CYCLE
              END IF
              NoModes = 1
              NoFields = 1
            ELSE          
              IF( MaxModes > 0 ) NoModes = MIN( MaxModes, NoModes )
              NoFields = NoModes
            END IF
          ELSE
            NoModes = 0 
            NoFields = 1
          END IF
          
          dofs = Solution % DOFs
          IF( ComponentVector ) THEN
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 2')
            IF( ASSOCIATED(Solution)) dofs = 2
            Solution => VariableGet( Model % Mesh % Variables, TRIM(FieldName)//' 3')
            IF( ASSOCIATED(Solution)) dofs = 3
          END IF
          
          IF( dofs > 1 ) THEN
            sdofs = MAX(dofs,3)
          ELSE
            sdofs = 1
          END IF


          DO iField = 1, NoFields
            
            IF( NoModes == 0 .OR. EigenAnalysis ) THEN
              FullName = TRIM( FieldName ) 
            ELSE          
              IF( GotActiveModes ) THEN
                IndField = ActiveModes( iField ) 
              ELSE
                IndField = iField
              END IF
              WRITE( FullName,'(A,I0)') TRIM( FieldName )//' mode',IndField
            END IF
            
            IF( AsciiOutput ) THEN
              WRITE( VtuUnit,'(A,A,I1,A)') '      <PDataArray type="Float64" Name="'//TRIM(FullName), &
                  '" NumberOfComponents="',sdofs,'" format="ascii"/>'  
            ELSE 
              WRITE( VtuUnit,'(A,I0,A,A,I0,A)') '      <PDataArray type="Float',PrecBits,'" Name="'//TRIM(FullName), &
                  '" NumberOfComponents="',sdofs,'" format="appended"/>'  
            END IF
          END DO

        END DO
      END DO
    END IF
  END IF

    ! Body and BC indexes
    IF( WriteIds ) THEN
      IF( AsciiOutput ) THEN
        WRITE( VtuUnit,'(A)') '      <PDataArray type="Int32" Name="GeometryIds" format="ascii"/>'
      ELSE
        WRITE( VtuUnit,'(A)') '      <PDataArray type="Int32" Name="GeometryIds" format="appended"/>'        
      END IF
    END IF

    WRITE( VtuUnit,'(A)') '    </PCellData>'

    ! Coordinates of each point
    !-------------------------------------
    WRITE( VtuUnit,'(A)') '    <PPoints>'
    IF( AsciiOutput ) THEN
      WRITE( VtuUnit,'(A,I0,A)') '      <DataArray type="Float64" NumberOfComponents="',dim,'" format="ascii"/>'    
    ELSE 
      WRITE( VtuUnit,'(A,I0,A,I0,A)') '      <DataArray type="Float',PrecBits,'" NumberOfComponents="',dim,'" format="appended"/>'    
    END IF
    WRITE( VtuUnit,'(A)') '    </PPoints>' 

 
    ! Write the pieces to the file 
    !-------------------------------------
    j = INDEX( FilePrefix,'/') 
    IF( j == 0 ) THEN
      ShortName = FilePrefix
    ELSE
      ShortName = FilePrefix(j+1:)
    END IF

    DO i=1,Partitions
      IF(.NOT. AllActive ) THEN
        IF( .NOT. ActivePartition(i)) CYCLE
      END IF

      IF( NoFileindex ) THEN
        WRITE( VtuUnit,'(A,I4.4,A,A)' ) '    <Piece Source="'//&
            TRIM(ShortName),i,"par",'.vtu"/>'
      ELSE IF( FileIndex < 10000 ) THEN
        WRITE( VtuUnit,'(A,I4.4,A,I4.4,A)' ) '    <Piece Source="'//&
            TRIM(ShortName),i,"par",FileIndex,'.vtu"/>'        
      ELSE
        WRITE( VtuUnit,'(A,I4.4,A,I0,A)' ) '    <Piece Source="'//&
            TRIM(ShortName),i,"par",FileIndex,'.vtu"/>'        
      END IF
    END DO

    WRITE( VtuUnit,'(A)') '  </PUnstructuredGrid>'
    WRITE( VtuUnit,'(A)') '</VTKFile>'

    CLOSE( VtuUnit )

    IF(.NOT. AllActive ) DEALLOCATE( ActivePartition ) 
  
  END SUBROUTINE WritePvtuFile

! http://www.vtk.org/pipermail/vtk-developers/2007-June/004606.html
! These cell types are supported:
!     VTK_POLY_VERTEX                       2
!     VTK_QUAD                              9
!     VTK_HEXAHEDRON                        12
!     VTK_QUADRATIC_TRIANGLE                22
!     VTK_QUADRATIC_QUAD                    23
!     VTK_QUADRATIC_HEXAHEDRON              25
!     VTK_QUADRATIC_WEDGE                   26
!     VTK_BIQUADRATIC_QUAD                  28
!     VTK_TRIQUADRATIC_HEXAHEDRON           29
!     VTK_QUADRATIC_LINEAR_QUAD             30
!     VTK_QUADRATIC_LINEAR_WEDGE            31
!     VTK_BIQUADRATIC_QUADRATIC_WEDGE       32
!     VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON  33
!     VTK_CUBIC_EDGE                        34
!     VTK_BICUBIC_QUAD                      36
!     VTK_TRICUBIC_HEXAHEDRON               37
  FUNCTION Elmer2VtkElement( ElmerCode ) RESULT ( VTKCode )
    INTEGER :: ElmerCode, VTKCode
!     write (*,*) "ElmerCode=",ElmerCode
    SELECT CASE (ElmerCode)
    CASE( 101 )
      VTKCode = 1
    CASE( 202 )
      VTKCode = 3
    CASE( 203 )
      VTKCode = 21
    CASE( 303 )
      VTKCode = 5
    CASE( 306 )
      VTKCode = 22
    CASE( 404 )
      VTKCode = 9
    CASE( 408 )
      VTKCode = 23   ! VTK_QUADRATIC_QUAD
    CASE( 409 )
      VTKCode = 28   ! VTK_BIQUADRATIC_QUAD
    CASE( 504 )
      VTKCode = 10
    CASE( 510 )
      VTKCode = 24
    CASE( 605 )
      VTKCode = 14
    CASE( 613 )
      VTKCode = 27
    CASE( 706 )
      VTKCode = 13
    CASE( 715 ) 
      VTKCode = 26
    CASE( 808 )
      VTKCode = 12  ! VTK_HEXAHEDRON
    CASE( 820 )
      VTKCode = 25  ! VTK_QUADRATIC_HEXAHEDRON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    CASE( 827 )
      VTKCode = 29  !      VTK_TRIQUADRATIC_HEXAHEDRON           29
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    CASE DEFAULT
      WRITE(Message,'(A,I0)') 'Not implemented for elementtype: ',ElmerCode
      CALL Fatal('Elmer2VtkElement','Not Implemented for elementtype')
      
    END SELECT
  END FUNCTION Elmer2VtkElement


!------------------------------------------------------------------------------
END SUBROUTINE VtuOutputSolver
!------------------------------------------------------------------------------
END MODULE ELMERModResOut



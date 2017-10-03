!> \ingroup IMNS
!> \}
! # YG@2016
! # yann.guevel@univ-ubs.fr
! # Part of ANM continuation/detection/branch-switching in ELMER FEM
! # Jean-Marc Cadou / Gregory girault / Yann Guevel
! # Institut de Recherche Dupuy de Lome - Universite Bretagne Sud


!------------------------------------------------------------------------------
!>  Utility routines for ANM CFD computations
!------------------------------------------------------------------------------
MODULE ANMToolBox

!------------------------------------------------------------------------------

  USE SolverUtils
  USE DefUtils

  USE FEMtoolBox    

  
!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
    

CONTAINS

!  ----  SERIES  ----
! + TANGENT and RHS OPERATORS        : MANTangentOperatorCompose
! + FQk                              : MANFQLocalOPTI
! + AMAX                             : ANMseriesRoV
! + U* = U/L                         : ANMExtractULambda          
! + U(amax) and L(amax)              : ANMSolPol

! TODO:
!  -> PADE V2 routine de racine plus stable!!!!
!  -> Multiple solution over one step : help visualisation of branches


! DONE 
! +                                  : ANMPadeBifDetect
! + Gram Schmidt Modif ORTHOG        : ANMOrthogGSm
! + COMPUTE DENOMINATOR FOR PADE     : DENPADE
! + GET ROOTS OF PADE, SMALLST POLE  : ANMPadeRACINES
! + COMPUTE ROOTS OF POLYNOMIAL      : BAI
! + COMPUTE ROOTS OF POLYNOMIAL BAI  : RESOUT1
! + COMPUTE ROOTS OF POLYNOMIAL BAI  : RESOUT2
! + COMPUTE ROOTS OF POLYNOMIAL BAI  : RAPHSON
! + U(AmaxPad,PADE), L(AmaxPad,PADE) : ANMSolPad
! + Range of Validity RATIONAL REP   : ANMPadeRoV


!  ----  BIFURCATION  ----
!  -> INDICATOR : 
!     -- f alea
!     -- Fnl
!     -- DeltaU_o init avec mu_o=1
!     -- K X_k = F_k  => mu_k => Delta U_k
!     -- Rayon {Delta U}  VS Rayon {U}

!  -> INDICATOR : 
! + Detection CM 2012          : BifurcationDetection
! 
!  ----  BIFURCATION  ----
!  + BIF1SolveAugSystem

!------------------------------------------------------------------------------
! SERIES COMPUTATION





   
   
   
   
   
   
   
 
 
 
 
 
! U and Lambda as  K U = lambda F , but K U* = F is solved  -> U* = U/lambda




! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! ! ! CalCondLimCoeff Reads boundary condition coutour values
! ! ! ! ! ! ! ! ! ! ! ! ! !   - Constraints : like a 0 on velocity Dof
! ! ! ! ! ! ! ! ! ! ! ! ! !   - Load : given velocity profile on contour
! ! ! ! ! ! ! ! ! ! ! ! ! ! Argument COND
! ! ! ! ! ! ! ! ! ! ! ! ! ! 0 : put a given vector for both Constraints and Load
! ! ! ! ! ! ! ! ! ! ! ! ! ! 1 : multiply Both contraints and load by a scalar
! ! ! ! ! ! ! ! ! ! ! ! ! ! 2 : Only Constraints (CALIMP EVE)
! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! !> ANM : For ANM orders > 1 set to Zero all the RHS DoF in Boundary 
! ! ! ! ! ! ! ! ! ! ! ! ! !>       RHS Velocity components <- 0 
! ! ! ! ! ! ! ! ! ! ! ! ! !>       EVE CALIMP and CALCONDILM
! ! ! ! ! ! ! ! ! ! ! ! ! !> Sets the Dirichlet conditions related to the variables of the active solver.
! ! ! ! ! ! ! ! ! ! ! ! ! !>
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! e1 bndry p1 p2 type n1 ... nn
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! e2 bndry p1 p2 type n1 ... nn
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ...
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! en bndry p1 p2 type n1 ... nn
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ++ The first integer is again the identification number of the boundary element. 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ++ Next the identification number of the part of the boundary where this element is located is given.
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ++ Whether the boundary element can be represented as the side of a parent element defined in the file 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! mesh.elements is indicated using the two parent element numbers p1 and p2. 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! -> If the boundary element  is located on an outer boundary of the body, it has only one parent 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! element and either of these  two integers may be set to be zero. 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! It is also possible that both parent element numbers are zeros. 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ++ Finally the element type code and element nodes are listed.
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! !   SUBROUTINE CalCondLimCoeff( Vector, Coeff, Cond, VectBCInit, USolver,Ux,UOffset,OffDiagonalMatrix,PiolaCurlTransform )
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! !             USE Adaptive
! ! ! ! ! ! ! ! ! ! ! ! !             USE SolverUtils
! ! ! ! ! ! ! ! ! ! ! ! !             USE DefUtils
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      IMPLICIT NONE
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !      INTEGER, OPTIONAL :: UOffset
! ! ! ! ! ! ! ! ! ! ! ! !      REAL(KIND=dp) :: Coeff ! reel dp 0 ou "lambda"
! ! ! ! ! ! ! ! ! ! ! ! !      INTEGER :: Cond ! 0 si on replaque les CIs et BCs , 1 sinon
! ! ! ! ! ! ! ! ! ! ! ! !      REAL(KIND=dp), OPTIONAL, DIMENSION(:) :: VectBCInit
! ! ! ! ! ! ! ! ! ! ! ! !      LOGICAL, OPTIONAL :: OffDiagonalMatrix
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Variable_t), OPTIONAL, TARGET :: Ux
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Solver_t), OPTIONAL, TARGET :: USolver
! ! ! ! ! ! ! ! ! ! ! ! !      LOGICAL, OPTIONAL :: PiolaCurlTransform  ! An additional argument for indicating that
! ! ! ! ! ! ! ! ! ! ! ! !                                               ! the solution is expanded in terms of H(curl)-
! ! ! ! ! ! ! ! ! ! ! ! !                                               ! conforming basis functions defined via the 
! ! ! ! ! ! ! ! ! ! ! ! !                                               ! Piola transform.  
! ! ! ! ! ! ! ! ! ! ! ! ! !--------------------------------------------------------------------------------------------     
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Matrix_t), POINTER   :: A
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Variable_t), POINTER :: x
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Solver_t), POINTER :: Solver
! ! ! ! ! ! ! ! ! ! ! ! !      REAL(KIND=dp), POINTER    :: b(:)
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      REAL(KIND=dp) :: xx, temp
! ! ! ! ! ! ! ! ! ! ! ! !      REAL(KIND=dp), POINTER :: DiagScaling(:)
! ! ! ! ! ! ! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE :: Work(:), STIFF(:,:)
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      INTEGER, ALLOCATABLE :: lInd(:), gInd(:)
! ! ! ! ! ! ! ! ! ! ! ! !      INTEGER :: i,j, k, kk, l, m, n,nd, nb, mb, nn, ni, nj, &
! ! ! ! ! ! ! ! ! ! ! ! !           DOF, local, numEdgeDofs,istat, n_start, Offset, iter, &     
! ! ! ! ! ! ! ! ! ! ! ! !           n1, n2, n3, n4
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      LOGICAL :: Flag,Found, ConstantValue, ScaleSystem
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(ValueList_t), POINTER :: BC, ptr, Params
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Element_t), POINTER :: Element, Parent, Edge, Face, SaveElement
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      CHARACTER(LEN=MAX_NAME_LEN) :: name
! ! ! ! ! ! ! ! ! ! ! ! !      LOGICAL :: BUpd
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      INTEGER::iii=0
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !     TYPE(Nodes_t) :: Nodes
! ! ! ! ! ! ! ! ! ! ! ! !     REAL(KIND=dp) :: xip,yip,zip,s,DetJ,Load
! ! ! ! ! ! ! ! ! ! ! ! !         REAL(KIND=dp) :: Basis(50)
! ! ! ! ! ! ! ! ! ! ! ! !     TYPE(GaussIntegrationPoints_t) :: IP
! ! ! ! ! ! ! ! ! ! ! ! !     LOGICAL :: stat
! ! ! ! ! ! ! ! ! ! ! ! !     INTEGER :: p,q,t,ndd,bcid
! ! ! ! ! ! ! ! ! ! ! ! !     
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! !    Right Hand Side for ANM high order 
! ! ! ! ! ! ! ! ! ! ! ! !      REAL(KIND=dp) :: Vector(:)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !      SAVE gInd, lInd, STIFF, Work
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      IF ( PRESENT( USolver ) ) THEN
! ! ! ! ! ! ! ! ! ! ! ! !         Solver => USolver
! ! ! ! ! ! ! ! ! ! ! ! !      ELSE
! ! ! ! ! ! ! ! ! ! ! ! !         Solver => CurrentModel % Solver
! ! ! ! ! ! ! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !      Params => GetSolverParams(Solver)
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! !---------------------------------------------------     
! ! ! ! ! ! ! ! ! ! ! ! !      A => Solver % Matrix
! ! ! ! ! ! ! ! ! ! ! ! !      b => A % RHS
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !      IF ( PRESENT(Ux) ) THEN
! ! ! ! ! ! ! ! ! ! ! ! !        x => Ux
! ! ! ! ! ! ! ! ! ! ! ! !      ELSE
! ! ! ! ! ! ! ! ! ! ! ! !        x => Solver % Variable
! ! ! ! ! ! ! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! ! ! ! ! ! !      IF(.NOT.ALLOCATED(A % ConstrainedDOF)) &
! ! ! ! ! ! ! ! ! ! ! ! !        ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
! ! ! ! ! ! ! ! ! ! ! ! !      A % ConstrainedDOF = .FALSE.
! ! ! ! ! ! ! ! ! ! ! ! ! !---------------------------------------------------     
! ! ! ! ! ! ! ! ! ! ! ! !      Offset = 0
! ! ! ! ! ! ! ! ! ! ! ! !      IF(PRESENT(UOffset)) Offset=UOffset
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      n = Solver % Mesh % MaxElementDOFs
! ! ! ! ! ! ! ! ! ! ! ! !      IF ( .NOT. ALLOCATED( gInd ) ) THEN
! ! ! ! ! ! ! ! ! ! ! ! !        ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
! ! ! ! ! ! ! ! ! ! ! ! !        IF ( istat /= 0 ) &
! ! ! ! ! ! ! ! ! ! ! ! !           CALL Fatal('ANM::CalCondLimCoeff','Memory allocation failed.' )
! ! ! ! ! ! ! ! ! ! ! ! !      ELSE IF ( SIZE(gInd) < n ) THEN
! ! ! ! ! ! ! ! ! ! ! ! !        DEALLOCATE( gInd, lInd, STIFF, Work )
! ! ! ! ! ! ! ! ! ! ! ! !        ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
! ! ! ! ! ! ! ! ! ! ! ! !        IF ( istat /= 0 ) &
! ! ! ! ! ! ! ! ! ! ! ! !           CALL Fatal('ANM::CalCondLimCoeff','Memory allocation failed.' )
! ! ! ! ! ! ! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! ! ! ! ! ! !      IF ( x % DOFs > 1 ) THEN
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !         CALL SetDirichletBoundaries( CurrentModel,A, b, GetVarName(x),-1,x % DOFs,x % Perm )
! ! ! ! ! ! ! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! ! ! ! ! ! !      CALL Info('ANM::CalCondLimCoeff', &
! ! ! ! ! ! ! ! ! ! ! ! !             'Setting Dirichlet boundary conditions for the ANM RHS', Level=5)
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! !      WRITE(*,*) "----------------------------------------------------------"
! ! ! ! ! ! ! ! ! ! ! ! ! !      WRITE(*,*) "RECAP TEST 3D"
! ! ! ! ! ! ! ! ! ! ! ! ! !      WRITE(*,*) "x % DOFs = ",x % DOFs
! ! ! ! ! ! ! ! ! ! ! ! ! !      WRITE(*,*) "Solver % Mesh % NumberOfBoundaryElements = ",Solver % Mesh % NumberOfBoundaryElements
! ! ! ! ! ! ! ! ! ! ! ! ! !      WRITE(*,*) "----------------------------------------------------------"            
! ! ! ! ! ! ! ! ! ! ! ! !      ! Set Dirichlet dofs for edges and faces
! ! ! ! ! ! ! ! ! ! ! ! !      ConstantValue = .FALSE.
! ! ! ! ! ! ! ! ! ! ! ! !      DO DOF=1,x % DOFs
! ! ! ! ! ! ! ! ! ! ! ! !         name = x % name
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !         IF ( x % DOFs > 1 ) name = ComponentName(name,DOF)
! ! ! ! ! ! ! ! ! ! ! ! !         ! clear bc face & edge dofs
! ! ! ! ! ! ! ! ! ! ! ! !         SaveElement => CurrentModel % CurrentElement
! ! ! ! ! ! ! ! ! ! ! ! !         DO i=1,Solver % Mesh % NumberOfBoundaryElements
! ! ! ! ! ! ! ! ! ! ! ! !            Element => GetBoundaryElement(i)
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ActiveBoundaryElement() ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !            ! Get parent element:
! ! ! ! ! ! ! ! ! ! ! ! !            ! -------------------
! ! ! ! ! ! ! ! ! ! ! ! !            Parent => Element % BoundaryInfo % Left
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ASSOCIATED( Parent ) ) &
! ! ! ! ! ! ! ! ! ! ! ! !              Parent => Element % BoundaryInfo % Right
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ASSOCIATED(Parent) )   CYCLE
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !            BC => GetBC()
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ASSOCIATED(BC) ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! !            ptr => ListFind(BC, Name,Found )
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ASSOCIATED(ptr) ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! !            ConstantValue =  ptr % PROCEDURE == 0 .AND. &
! ! ! ! ! ! ! ! ! ! ! ! !              ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! !            ! Get indexes for boundary and values for dofs associated to them
! ! ! ! ! ! ! ! ! ! ! ! !            n = GetElementNOFNodes()
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( isActivePElement(Parent)) THEN
! ! ! ! ! ! ! ! ! ! ! ! !              CALL getBoundaryIndexes( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )
! ! ! ! ! ! ! ! ! ! ! ! !            ELSE
! ! ! ! ! ! ! ! ! ! ! ! !              CYCLE 
! ! ! ! ! ! ! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            ! Contribute this boundary to global system
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            ! (i.e solve global boundary problem)
! ! ! ! ! ! ! ! ! ! ! ! !            DO k=n+1,numEdgeDofs
! ! ! ! ! ! ! ! ! ! ! ! ! !              write(*,*) 'gInd',gind(k)
! ! ! ! ! ! ! ! ! ! ! ! !              nb = x % Perm( gInd(k) )
! ! ! ! ! ! ! ! ! ! ! ! !              IF ( nb <= 0 ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! !              nb = Offset + x % DOFs * (nb-1) + DOF
! ! ! ! ! ! ! ! ! ! ! ! !            END DO
! ! ! ! ! ! ! ! ! ! ! ! !         END DO
! ! ! ! ! ! ! ! ! ! ! ! !         CurrentModel % CurrentElement => SaveElement
! ! ! ! ! ! ! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! ! ! ! ! ! !  
! ! ! ! ! ! ! ! ! ! ! ! !      ! Set Dirichlet dofs for edges and faces
! ! ! ! ! ! ! ! ! ! ! ! !     
! ! ! ! ! ! ! ! ! ! ! ! !      DO DOF=1,x % DOFs
! ! ! ! ! ! ! ! ! ! ! ! !         name = x % name
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !         IF (x % DOFs>1) name=ComponentName(name,DOF)
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) 'DOFs=',DOF,'Name=', name
! ! ! ! ! ! ! ! ! ! ! ! !         CALL SetNodalLoads( CurrentModel,A, b, &
! ! ! ! ! ! ! ! ! ! ! ! !             Name,DOF,x % DOFs,x % Perm ) ! , Offset ) not yet ?
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !         CALL SetDirichletBoundaries( CurrentModel, A, b, &
! ! ! ! ! ! ! ! ! ! ! ! !              Name, DOF, x % DOFs, x % Perm, Offset, OffDiagonalMatrix )
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! !       Dirichlet BCs for face & edge DOFs:
! ! ! ! ! ! ! ! ! ! ! ! ! !       -----------------------------------
! ! ! ! ! ! ! ! ! ! ! ! !         SaveElement => CurrentModel % CurrentElement
! ! ! ! ! ! ! ! ! ! ! ! !         
! ! ! ! ! ! ! ! ! ! ! ! !         DO i=1,Solver % Mesh % NumberOfBoundaryElements
! ! ! ! ! ! ! ! ! ! ! ! ! !            write(*,*) "i=",i
! ! ! ! ! ! ! ! ! ! ! ! !            Element => GetBoundaryElement(i)
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ActiveBoundaryElement() ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! ! !            write(*,*) "GetBC()"
! ! ! ! ! ! ! ! ! ! ! ! !            BC => GetBC()
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ASSOCIATED(BC) ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ListCheckPresent(BC, Name) .AND. &
! ! ! ! ! ! ! ! ! ! ! ! !                 .NOT. ListCheckPresent(BC, TRIM(Name)//' {e}') .AND. &
! ! ! ! ! ! ! ! ! ! ! ! !                 .NOT. ListCheckPresent(BC, TRIM(Name)//' {f}') ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! !            ! Get parent element:
! ! ! ! ! ! ! ! ! ! ! ! !            ! -------------------
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! !            Parent => Element % BoundaryInfo % Left
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ASSOCIATED( Parent ) ) THEN
! ! ! ! ! ! ! ! ! ! ! ! !                Parent => Element % BoundaryInfo % Right
! ! ! ! ! ! ! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( .NOT. ASSOCIATED( Parent ) )   CYCLE
! ! ! ! ! ! ! ! ! ! ! ! ! !            write(*,*) " ASSOCIATED( Parent ) OK"
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! ! !---------------------------------------------------     
! ! ! ! ! ! ! ! ! ! ! ! ! !!! TEST COMMENTE avant pour 2D valider
! ! ! ! ! ! ! ! ! ! ! ! ! !!!! TEST ISACTIVEPELEMENT...?
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) "isPelement(Parent)=",isPelement(Parent)
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) Parent % Type % ElementCode / 100
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) CurrentModel % Solver % Def_Dofs(:,:,:)
! ! ! ! ! ! ! ! ! ! ! ! ! !         WRITE(*,*) "=====Element % NodeIndexes(1:n)=",Element % NodeIndexes(1:n)
! ! ! ! ! ! ! ! ! ! ! ! !         
! ! ! ! ! ! ! ! ! ! ! ! ! !!!!!!
! ! ! ! ! ! ! ! ! ! ! ! ! !!! TEST DECOMMENTE pour 3D
! ! ! ! ! ! ! ! ! ! ! ! ! !           WRITE(*,*) "CYCLE if .NOT.isActivePElement(Parent))=",.NOT.isActivePElement(Parent)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            IF (.NOT.isActivePElement(Parent)) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! !            ptr => ListFind(BC, Name,Found )
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! !            WRITE (*,*) "BC, NAME=",Name
! ! ! ! ! ! ! ! ! ! ! ! ! ! !            bcid=GetBCId( Element )
! ! ! ! ! ! ! ! ! ! ! ! ! ! !            WRITE (*,*) "BC, id=",bcid
! ! ! ! ! ! ! ! ! ! ! ! ! ! !            WRITE (*,*) "ptr % FValues",ptr % FValues 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! !      ! Get nodes of boundary elements parent and gauss points for boundary
! ! ! ! ! ! ! ! ! ! ! ! ! ! !     CALL GetElementNodes( Nodes, Element )
! ! ! ! ! ! ! ! ! ! ! ! ! ! !     IP = GaussPoints( Element )
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! !     DO t=1,IP % n
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        stat = ElementInfo( Element, Nodes, IP % u(t), &
! ! ! ! ! ! ! ! ! ! ! ! ! ! !           IP % v(t), IP % w(t), DetJ, Basis )
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        s = IP % s(t) * DetJ
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        ! Get value of boundary condition
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        ndd = Element % TYPE % NumberOfNodes
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        xip = SUM( Basis(1:ndd) * Nodes % x(1:ndd) )
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        yip = SUM( Basis(1:ndd) * Nodes % y(1:ndd) )
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        zip = SUM( Basis(1:ndd) * Nodes % z(1:ndd) )
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        Load = ListGetConstReal( BC, Name, x=xip,y=yip,z=zip )
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        write(*,*) "Name=",TRIM(Name)," - BcID=",GetBCId( Element )," - Load(",xip,",",yip,",",zip,")=",Load
! ! ! ! ! ! ! ! ! ! ! ! ! ! !     END DO
! ! ! ! ! ! ! ! ! ! ! ! ! ! !        
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************           
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) "COUCOU4"           
! ! ! ! ! ! ! ! ! ! ! ! !               ! Number of nodes for this element
! ! ! ! ! ! ! ! ! ! ! ! !               n = Element % TYPE % NumberOfNodes
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) "n = Element % TYPE % NumberOfNodes =",n          
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! Get indexes for boundary and values for dofs associated to them
! ! ! ! ! ! ! ! ! ! ! ! ! !> Calculate global indexes of boundary dofs for given element and its boundary.  
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !               CALL getBoundaryIndexesMAN( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )
! ! ! ! ! ! ! ! ! ! ! ! !         
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) "=====gInd=",gInd
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               CALL LocalBcBDOFs( BC, Element, numEdgeDofs, Name, STIFF, Work )
! ! ! ! ! ! ! ! ! ! ! ! ! !         write(*,*) "COUCOU6"
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !        
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! !                 DO l=1,n
! ! ! ! ! ! ! ! ! ! ! ! !                   nb = x % Perm( gInd(l) )
! ! ! ! ! ! ! ! ! ! ! ! ! !                   write(*,*) "nb=",nb
! ! ! ! ! ! ! ! ! ! ! ! ! ! 20140916 - MODIF GIND PAR NODEINDEX
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                   nb = x % Perm(Element % NodeIndexes(l))
! ! ! ! ! ! ! ! ! ! ! ! !                   
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !                   nb = Offset + x % DOFs * (nb-1) + DOF
! ! ! ! ! ! ! ! ! ! ! ! ! ! !                   write(*,*) "nb=",nb
! ! ! ! ! ! ! ! ! ! ! ! !                   IF ( Cond == 0 ) THEN
! ! ! ! ! ! ! ! ! ! ! ! ! !                      COUNTOUR <- COUTOUR GIVEN
! ! ! ! ! ! ! ! ! ! ! ! !                     Vector(nb) = VectBCInit(nb)
! ! ! ! ! ! ! ! ! ! ! ! !                   ELSE IF ( Cond == 1 ) THEN
! ! ! ! ! ! ! ! ! ! ! ! ! !                      USED AS CALCONDLIM + CALIMP IF Coef==0                     
! ! ! ! ! ! ! ! ! ! ! ! ! !                       vector <- coeff X COUNTOUR
! ! ! ! ! ! ! ! ! ! ! ! ! !                       vector <- 0
! ! ! ! ! ! ! ! ! ! ! ! ! !                       vector <- Lambda X Ud on COUNTOUR
! ! ! ! ! ! ! ! ! ! ! ! !                     IF(VectBCInit(nb).GT.1e-15)THEN
! ! ! ! ! ! ! ! ! ! ! ! !                        Vector(nb) = Coeff * VectBCInit(nb)
! ! ! ! ! ! ! ! ! ! ! ! !                     ELSE
! ! ! ! ! ! ! ! ! ! ! ! !                       Vector(nb) = 0.0_dp
! ! ! ! ! ! ! ! ! ! ! ! !                     ENDIF
! ! ! ! ! ! ! ! ! ! ! ! ! !                      write(*,*) TRIM(Name),l," Coeff=",Coeff,"x % Values(nb)=",x % Values(nb)
! ! ! ! ! ! ! ! ! ! ! ! ! !                      Vector(nb) = Coeff * x % Values(nb)
! ! ! ! ! ! ! ! ! ! ! ! ! !                      write(*,*) TRIM(Name),l,"Vector(nb)=",Vector(nb)
! ! ! ! ! ! ! ! ! ! ! ! !                     
! ! ! ! ! ! ! ! ! ! ! ! !                   ELSE IF ( Cond == 2 ) THEN
! ! ! ! ! ! ! ! ! ! ! ! ! !                      USED AS CALIMP without user "load"
! ! ! ! ! ! ! ! ! ! ! ! ! !                      IF Value is 0 on given BC then put 0 on current vector
! ! ! ! ! ! ! ! ! ! ! ! ! !                      IF Value is not 0 for the BC then don't touch the given vector
! ! ! ! ! ! ! ! ! ! ! ! !                     IF (VectBCInit(nb).LT.1e-15) THEN
! ! ! ! ! ! ! ! ! ! ! ! !                       Vector(nb) = 0.0_dp
! ! ! ! ! ! ! ! ! ! ! ! !                     ENDIF
! ! ! ! ! ! ! ! ! ! ! ! !                   END IF
! ! ! ! ! ! ! ! ! ! ! ! ! !                   write(*,*) "Vector(nb)",Vector(nb)
! ! ! ! ! ! ! ! ! ! ! ! !                 END DO   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ******************************************************
! ! ! ! ! ! ! ! ! ! ! ! !           
! ! ! ! ! ! ! ! ! ! ! ! !            ConstantValue =  ptr % PROCEDURE == 0 .AND. &
! ! ! ! ! ! ! ! ! ! ! ! !              ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR
! ! ! ! ! ! ! ! ! ! ! ! !            IF ( ConstantValue ) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !         END DO
! ! ! ! ! ! ! ! ! ! ! ! !         CurrentModel % CurrentElement => SaveElement
! ! ! ! ! ! ! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      CALL Info('ANM::CalCondLimCoeff','Dirichlet boundary conditions are set for the ANMD RHS', Level=5)
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! CONTAINS
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! ! !> Calculate global indexes of boundary dofs for given element and its boundary.
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! !    SUBROUTINE getBoundaryIndexesMAN( Mesh, Element, Parent, Indexes, indSize )
! ! ! ! ! ! ! ! ! ! ! ! ! !!!!!!!
! ! ! ! ! ! ! ! ! ! ! ! !             USE Adaptive
! ! ! ! ! ! ! ! ! ! ! ! !             USE SolverUtils
! ! ! ! ! ! ! ! ! ! ! ! !             USE DefUtils
! ! ! ! ! ! ! ! ! ! ! ! ! !!!!!!!
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! !    Type(Mesh_t) :: Mesh
! ! ! ! ! ! ! ! ! ! ! ! ! !      INPUT: Finite element mesh containing edges and faces of elements
! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! !    Type(Element_t) :: Element
! ! ! ! ! ! ! ! ! ! ! ! ! !      INPUT: Boundary element to get indexes for
! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! !    Type(Element_t) :: Parent
! ! ! ! ! ! ! ! ! ! ! ! ! !      INPUT: Parent of boundary element to get indexes for
! ! ! ! ! ! ! ! ! ! ! ! ! !
! ! ! ! ! ! ! ! ! ! ! ! ! !    INTEGER :: Indexes(:)
! ! ! ! ! ! ! ! ! ! ! ! ! !      OUTPUT: Calculated indexes of boundary element in global system
! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! !    INTEGER :: indSize
! ! ! ! ! ! ! ! ! ! ! ! ! !      OUTPUT: Size of created index vector, i.e. how many indexes were created
! ! ! ! ! ! ! ! ! ! ! ! ! !        starting from index 1
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! !      IMPLICIT NONE
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      ! Parameters
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Mesh_t) :: Mesh
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Element_t) :: Parent
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Element_t), POINTER :: Element
! ! ! ! ! ! ! ! ! ! ! ! !      INTEGER :: indSize, Indexes(:)
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !      ! Variables
! ! ! ! ! ! ! ! ! ! ! ! !      TYPE(Element_t), POINTER :: Edge, Face
! ! ! ! ! ! ! ! ! ! ! ! !      INTEGER :: i,j,n
! ! ! ! ! ! ! ! ! ! ! ! ! !      write(*,*) "GBI INIT"
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      ! Clear indexes
! ! ! ! ! ! ! ! ! ! ! ! !      Indexes = 0
! ! ! ! ! ! ! ! ! ! ! ! !      n = Element % TYPE % NumberOfNodes
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !      ! Nodal indexes
! ! ! ! ! ! ! ! ! ! ! ! !      Indexes(1:n) = Element % NodeIndexes(1:n)
! ! ! ! ! ! ! ! ! ! ! ! ! !      write(*,*) "Indexes(1:n)=",Indexes(1:n)
! ! ! ! ! ! ! ! ! ! ! ! ! !      write(*,*) "Parent % TYPE % DIMENSION=",Parent % TYPE % DIMENSION
! ! ! ! ! ! ! ! ! ! ! ! !      ! Assign rest of indexes if neccessary
! ! ! ! ! ! ! ! ! ! ! ! !      SELECT CASE(Parent % TYPE % DIMENSION)
! ! ! ! ! ! ! ! ! ! ! ! !      CASE (1)
! ! ! ! ! ! ! ! ! ! ! ! !        indSize = n 
! ! ! ! ! ! ! ! ! ! ! ! !      CASE (2)
! ! ! ! ! ! ! ! ! ! ! ! ! !      write(*,*) "GBI2D - Element % BDOFs=",Element % BDOFs
! ! ! ! ! ! ! ! ! ! ! ! !         ! Add index for each bubble dof in edge
! ! ! ! ! ! ! ! ! ! ! ! !         DO i=1,Element % BDOFs
! ! ! ! ! ! ! ! ! ! ! ! !            n = n+1
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! !            IF (SIZE(Indexes) < n) THEN
! ! ! ! ! ! ! ! ! ! ! ! !               CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
! ! ! ! ! ! ! ! ! ! ! ! !               RETURN
! ! ! ! ! ! ! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! ! ! ! ! ! ! !         WRITE(*,*) "Element % PDefs % localNumber=",Element % PDefs % localNumber
! ! ! ! ! ! ! ! ! ! ! ! ! !         WRITE(*,*) "Parent % EdgeIndexes(Element % PDefs % localNumber)=",Parent % EdgeIndexes(Element % PDefs % localNumber)
! ! ! ! ! ! ! ! ! ! ! ! !  
! ! ! ! ! ! ! ! ! ! ! ! !            Indexes(n) = Mesh % NumberOfNodes + &
! ! ! ! ! ! ! ! ! ! ! ! !                 (Parent % EdgeIndexes(Element % PDefs % localNumber)-1) * Mesh % MaxEdgeDOFs + i
! ! ! ! ! ! ! ! ! ! ! ! !         END DO
! ! ! ! ! ! ! ! ! ! ! ! !      
! ! ! ! ! ! ! ! ! ! ! ! !         indSize = n 
! ! ! ! ! ! ! ! ! ! ! ! !      CASE (3)
! ! ! ! ! ! ! ! ! ! ! ! ! !      write(*,*) "GBI3D"
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         WRITE(*,*) "Element % PDefs % localNumber=",Element % PDefs % localNumber
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         WRITE(*,*) "Parent % FaceIndexes(previous)=",Parent % FaceIndexes(Element % PDefs % localNumber)
! ! ! ! ! ! ! ! ! ! ! ! !         ! Get boundary face
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         Face => Mesh % Faces( Parent % FaceIndexes(Element % PDefs % localNumber) )
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         WRITE(*,*) "GBI3D -  Face % TYPE % NumberOfEdges", Face % TYPE % NumberOfEdges
! ! ! ! ! ! ! ! ! ! ! ! !         ! Add indexes of faces edges 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         DO i=1, Face % TYPE % NumberOfEdges
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            Edge => Mesh % Edges( Face % EdgeIndexes(i) )
! ! ! ! ! ! ! ! ! ! ! ! !            
! ! ! ! ! ! ! ! ! ! ! ! !            ! If edge has no dofs jump to next edge
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            IF (Edge % BDOFs <= 0) CYCLE
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            DO j=1,Edge % BDOFs
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               n = n + 1
! ! ! ! ! ! ! ! ! ! ! ! !               
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               IF (SIZE(Indexes) < n) THEN
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                  CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                  RETURN
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               END IF
! ! ! ! ! ! ! ! ! ! ! ! !               
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               Indexes(n) = Mesh % NumberOfNodes +&
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                   ( Face % EdgeIndexes(i)-1)*Mesh % MaxEdgeDOFs + j
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            END DO
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         END DO
! ! ! ! ! ! ! ! ! ! ! ! !                
! ! ! ! ! ! ! ! ! ! ! ! !         ! Add indexes of faces bubbles
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         DO i=1,Face % BDOFs
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            n = n + 1
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            IF (SIZE(Indexes) < n) THEN
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               RETURN
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !            Indexes(n) = Mesh % NumberOfNodes + &
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                 Mesh % NumberOfEdges * Mesh % MaxEdgeDOFs + &
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                 (Parent % FaceIndexes( Element % PDefs % localNumber )-1) * Mesh % MaxFaceDOFs + i
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !         END DO        
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !         indSize = n
! ! ! ! ! ! ! ! ! ! ! ! !      CASE DEFAULT
! ! ! ! ! ! ! ! ! ! ! ! !         CALL Fatal('DefUtils::getBoundaryIndexes','Unsupported dimension')
! ! ! ! ! ! ! ! ! ! ! ! !      END SELECT
! ! ! ! ! ! ! ! ! ! ! ! ! !      WRITE(*,*) "Indexes after GBI SIZE(Indexes)=",SIZE(Indexes)
! ! ! ! ! ! ! ! ! ! ! ! ! !      DO i=1,SIZE(Indexes)
! ! ! ! ! ! ! ! ! ! ! ! ! !        WRITE(*,*) Indexes(i)
! ! ! ! ! ! ! ! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! !    END SUBROUTINE getBoundaryIndexesMAN
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! !   END SUBROUTINE CalCondLimCoeff
! ! ! ! ! ! ! ! ! ! ! ! ! !------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Get U and Lambda with U*
!------------------------------------------------------------------------------
 SUBROUTINE ANMExtractULambda( Lambda, UMan, OrdreMAN, NoPressure, VMan ) 
! UMan = (VMan Pressure)
      USE DirectSolve
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
      REAL(KIND=dp), DIMENSION(:,:) :: UMan,VMan
      REAL(KIND=dp), DIMENSION(:) :: Lambda, NoPressure
      INTEGER :: OrdreMAN
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: NU1,  projvect1 
!------------------------------------------------------------------------------       
      
      ! On extrait la vitesse
      CALL PressureFree( UMan(:,1), NoPressure, VMan(:,1) )
      ! Calcul de la norme sans les termes de pression
      NU1 = DSQRT( DOT_PRODUCT( VMan(:,1),VMan(:,1) ))      
      
      IF (OrdreMAN == 1) THEN
      ! K X1 = F
      ! X1 = U1 / L1
      ! <U1,U1> + L1*L1 = 1 <=> <X1,X1> + 1 = 1/(L1*L1)
      ! => L1 = 1/ SQRT(<X1,X1> + 1)
        Lambda( 1 ) = 1._dp/( NU1*NU1  + 1._dp )
        Lambda( 1 ) = DSQRT( Lambda(1) )
!       U1=lambda1*U1*
        UMan(:,1) = UMan(:,1) * Lambda(1)
        
        CALL PressureFree( UMan(:,1), NoPressure, VMan(:,1) )
        
      ELSE
        ! K Xp = -FQp
        ! K Up = Lp F - FQp
        ! WITH F = K U1/L1
        ! K Up = Lp K U1/L1 - FQp
        ! => K ( Up - Lp U1/L1) = -FQp
        ! On extrait la vitesse
        CALL PressureFree( UMan(:,OrdreMAN), NoPressure, VMan(:,OrdreMAN) )
        ! Produit scalaire de Uordrep et Uordre1
        projvect1 = DOT_PRODUCT( VMan(:,OrdreMAN) , VMan(:,1) )
        ! Calcul de lambdap
        Lambda( OrdreMAN ) = -( Lambda( 1 ) * projvect1 ) / ( Lambda(1)*Lambda(1) + NU1*NU1 )
        
        UMan(:,OrdreMAN)   = UMan(:,OrdreMAN) + UMan(:,1)*( Lambda(OrdreMAN)/Lambda(1) )
        
        CALL PressureFree( UMan(:,OrdreMAN), NoPressure, VMan(:,OrdreMAN) )

      END IF
      
        
 END SUBROUTINE ANMExtractULambda


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! RANGE OF VALIDITY OF THE SERIES
  FUNCTION ANMseriesRoV( VMan, ANMorder, Tolerance ) RESULT( amax )
  
     REAL(KIND=dp), DIMENSION(:,:) :: VMan
     REAL(KIND=dp) :: Tolerance, amax
     INTEGER :: ANMorder

     REAL(KIND=dp), ALLOCATABLE :: UTemp(:)
     REAL(KIND=dp) :: poww, norm1, normp
     
     norm1 = DSQRT( DOT_PRODUCT( VMan(:,1),VMan(:,1) ))
     WRITE(*,*) 'ANMseriesRoV :  ||U1||=',norm1
     normp = DSQRT( DOT_PRODUCT( VMan(:,ANMorder),VMan(:,ANMorder) ) )
     WRITE(*,*) 'ANMseriesRoV :  ||UN||=',normp
     poww = 1._dp / ( ANMorder-1) 
     amax = ( Tolerance * ( norm1 / normp ) ) ** poww
  
  END FUNCTION ANMseriesRoV
  
 
! SOLUTION POLY
!------------------------------------------------------------------------------
  SUBROUTINE ANMSolPol( U0, lambda0, UMan, lambda, ANMorder, amax )
!>  NEXT SOLUTION IS FORMED : U0 L0 with the actual computed series and amax  

     USE DefUtils

     IMPLICIT NONE
      
     REAL(KIND=dp), DIMENSION(:,:) :: UMan
     REAL(KIND=dp), DIMENSION(:) :: lambda, U0
     REAL(KIND=dp) :: amax, lambda0
     INTEGER :: ANMorder

     INTEGER :: i
     REAL(KIND=dp) :: pamax, lambdatemp
     REAL(KIND=dp), DIMENSION(SIZE(U0)) :: Utemp
        
     pamax = 1._dp
     Utemp = U0
     lambdatemp = lambda0
        
     DO i = 1 , ANMorder
       pamax = pamax * amax
       Utemp = Utemp + pamax * UMan( : , i )
       lambdatemp = lambdatemp + pamax * lambda( i )
     END DO
        
     U0 = Utemp
     lambda0 = lambdatemp
                
  END SUBROUTINE ANMSolPol
  
! SOLUTION POLY
!------------------------------------------------------------------------------
  SUBROUTINE ANMMultipleSolPol( U0, lambda0, UMan, lambda, ANMorder, amax,       &
                                NBSolByStep, NODEOUT, Step, FlowPerm, NSDOFs, dim)
!>  NEXT SOLUTION IS FORMED : U0 L0 with the actual computed series and amax  

     USE DefUtils

     IMPLICIT NONE
      
     REAL(KIND=dp), DIMENSION(:,:) :: UMan
     REAL(KIND=dp), DIMENSION(:) :: lambda, U0
     REAL(KIND=dp) :: amax, lambda0
     INTEGER :: ANMorder,NBSolByStep,NODEOUT,Step,NSDOFs,dim,FlowPerm(:)

     INTEGER :: i,k
     REAL(KIND=dp) :: pamax, lambdatemp,da,amaxtmp,CDBx,CDBy,CDBz
     REAL(KIND=dp), DIMENSION(SIZE(U0)) :: Utemp
        

     if (NBSolByStep < 1) NBSolByStep = 1
     da = amax/NBSolByStep
     
!      pamax = 1._dp       
     amaxtmp = 0.0_dp
     DO K = 1,NBSolByStep
       Utemp = U0
       lambdatemp = lambda0     
       pamax = 1._dp        
       amaxtmp = K * da
       DO i = 1 , ANMorder
         pamax = pamax * amaxtmp
         Utemp      = Utemp      + pamax * UMan( : , i )
         lambdatemp = lambdatemp + pamax * lambda( i )
       END DO

       CDBx = Utemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 1) )
       CDBy = Utemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 2) )
       CDBz = 0.0_dp
       if (dim == 3)  CDBz = Utemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 3) )       
              
       WRITE(144,903)  Step,lambdatemp,CDBx,CDBy,CDBz
       CALL FLUSH(144)
       
     ENDDO
     U0 = Utemp
     lambda0 = lambdatemp           
  903 FORMAT( 1X, I4, 2X, G15.8, 2X, G15.8, 2X, G15.8, 2X, G15.8)     
  END SUBROUTINE ANMMultipleSolPol
 
 
 
!------------------------------------------------------------------------------       
! ooooooooo.         .o.       oooooooooo.   oooooooooooo 
! `888   `Y88.      .888.      `888'   `Y8b  `888'     `8 
!  888   .d88'     .8"888.      888      888  888         
!  888ooo88P'     .8' `888.     888      888  888oooo8    
!  888           .88ooo8888.    888      888  888    "    
!  888          .8'     `888.   888     d88'  888       o 
! o888o        o88o     o8888o o888bood8P'   o888ooooood8 

! 1 - ORTHOGONALISER U -> AlphaPad
! 2 - COEFF DENOMINATEURS di
! 3 - RACINE PADE
! 4 - Coeff Padé
! 5 - Rayon padé

! TOOLBOX PADE: Gramm Schmidt, pole, ...
!
  SUBROUTINE ANMPadeBifDetect(amaxPad,U0,lambda0,UMan,VMan,lambda,VTMPOrtho,VSOLTMP,OMan, &
                              RAC_MIN,LambdaCRIT,U0Pad,L0Pad)
! EXTERNE
      REAL(KIND=dp) :: amaxPad, U0(:), lambda0,lambda(:), VSOLTMP(:) , LambdaCRIT,Apoly, L0Pad
      REAL(KIND=dp), DIMENSION(:,:) :: UMan(:,:),VMan(:,:),VTMPOrtho(:,:),U0Pad(:)
      INTEGER :: OMan      
! INTERNE      
      REAL(KIND=dp) :: AlphaPad(OMan,OMan), DN(OMan),DNM1(OMan-1),DD(OMan)
      REAL(KIND=dp) :: RAC(OMan),RACZRHQR(OMan),RAC_MIN_BIS
      REAL(KIND=dp) :: EPSILONPAD,RANGEMAX,RAC_MIN,RAC_MIN_C,DELTA,Ltmp
      REAL(KIND=dp) :: CRIT , NB_RAC_ZRHQR
      INTEGER :: I,J,K,N,NDL,NB_RAC,NB_RAC_POS,NB_RAC_C , NB_RAC_POS_BIS
      LOGICAL :: DIVERGE
      
      Apoly = amaxPad
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                          
!       WRITE(*,*) ' ANM PADE : ORTHOGONALISATION DE GRAM-SCHMIT '
      CALL ANMOrthogGSm( VTMPOrtho, AlphaPad, OMan )
!       WRITE(*,*) ' ANM PADE : COEFF di pour PADE '
      CALL DENPADE(AlphaPad,DN,OMan,OMan-1)
      CALL DENPADE(AlphaPad,DNM1,OMan-1,OMan-2)

      WRITE(*,*) ' ANM PADE : RECHERCHE DES POLE DES PADE '       
      DD = DN
      EPSILONPAD=1.D-10
      RANGEMAX=100.0_dp * Apoly
! !       CALL ANMPadeRACINESzrhqr(DD, OMan-1, EPSILONPAD, RAC_MIN, RAC_MIN_C, DIVERGE, & 
! !                     RANGEMAX , NB_RAC, NB_RAC_POS, NB_RAC_C,RAC,            &
! !                                 RACZRHQR,RAC_MIN_BIS,NB_RAC_ZRHQR)
! ! ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                        
      CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,LambdaCRIT,OMan-1,RAC_MIN,DN,DELTA)
             write (*,*) "*-*-*- ANM PADE : RACZR_MIN=",RAC_MIN," LracMin=",LambdaCRIT
      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                 
                                
                                
      CALL ANMPadeRACINES(DD, OMan-1, EPSILONPAD, RAC_MIN, RAC_MIN_C, DIVERGE, & 
                    RANGEMAX , NB_RAC, NB_RAC_POS, NB_RAC_C,RAC) !,               &
!                     RACZRHQR,RAC_MIN_BIS,NB_RAC_POS_BIS)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                        
      CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,LambdaCRIT,OMan-1,RAC_MIN,DD,DELTA)
             write (*,*) "*-*-*- ANM PADE : RACBAI_MIN=",RAC_MIN," LracMin=",LambdaCRIT
      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                    
      IF(NB_RAC_POS.GT.0)  THEN
        RANGEMAX= RAC_MIN
      ELSE 
        RANGEMAX= 100.0_dp * Apoly
      ENDIF
      WRITE(*,*) ' ANM PADE : RECHERCHE DU Rayon de Validité dans ',Apoly,RANGEMAX
      CRIT = 1.D-8
!     RAYON PADE DANS [Apoly , RAC_MIN ] 
      CALL ANMPadeRoV(CRIT,VMan,DN,DNM1,lambda0,U0,lambda,Apoly,OMan,amaxPad,RANGEMAX, &
                      U0Pad,L0Pad) 
      WRITE(*,*) ' ANM PADE : Rayon de Validité : amaxPad = ',amaxPad
!                     
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                    
!     TEST RACINE ET LAMBDA CRITIQUE POUR DEBUG et pour voir
      if (NB_RAC_POS.GT.0._dp) then
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         write(*,*) "*-*-*- ANM PADE     PADE's POLE "      
         ! VERIFIER POINT LIMITE
         !test lambda critique
         write(*,*) "*-*-*- ANM PADE     lambda0 ",lambda0
         
!          write(*,*) 'ANM PADE : lambda0 = ',lambda0
         DO I=1,NB_RAC_POS
           IF (RAC(I).GT.0._dp ) THEN
             CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,Ltmp,OMan-1,RAC(I),DD,DELTA)
             write (*,*) '*-*-*- ANM PADE : LambdaPADEavant (',RAC(I),') = ',Ltmp
           ENDIF
         ENDDO  
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         
      endif      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
  END SUBROUTINE ANMPadeBifDetect


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *          
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *   
! IN  : serie, solution
! OUT : range of validity of Padé approx
! side product : test of every real pole (to see....)
  SUBROUTINE ANMPadeBifDetectZRHQR(amaxPad, U0,lambda0, UMan,VMan,lambda,      &
                                   VTMPOrtho,VSOLTMP,OMan, RAC_MIN,LambdaCRIT, &
                                   U0Pad,L0Pad,Step,Apoly)
! EXTERNE
      REAL(KIND=dp) :: amaxPad, U0(:), lambda0,lambda(:), VSOLTMP(:) , LambdaCRIT,Apoly, L0Pad
      REAL(KIND=dp) :: RAC_MIN
      REAL(KIND=dp), DIMENSION(:,:) :: UMan(:,:),VMan(:,:),VTMPOrtho(:,:),U0Pad(:)
      INTEGER :: OMan,Step
! INTERNE      
      REAL(KIND=dp) :: AlphaPad(OMan,OMan), DN(OMan),DNM1(OMan-1),DD(OMan)
      REAL(KIND=dp) :: RAC(OMan)
      REAL(KIND=dp) :: EPSILONPAD, RANGEMAX, RAC_BAI_MIN, RAC_MIN_C, DELTA, Ltmp
      REAL(KIND=dp) :: CRIT
      !
      REAL(KIND=dp) :: RAC_R_POS_MIN
      REAL(KIND=dp) :: TAB_RAC_R_POS(OMan), TAB_RAC_R_NEG(OMan)
      INTEGER       :: NB_RAC_R_POS, NB_RAC_R_NEG
      !
      INTEGER       :: I,J,K,N,NDL,NB_RAC,NB_RAC_POS,NB_RAC_C
      LOGICAL       :: DIVERGE
      
      RAC_MIN = 0._dp
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                          
!       WRITE(*,*) ' ANM PADE : ORTHOGONALISATION DE GRAM-SCHMIT '
      CALL ANMOrthogGSm( VTMPOrtho, AlphaPad, OMan )
!       WRITE(*,*) ' ANM PADE : COEFF di pour PADE '
      CALL DENPADE(AlphaPad,DN,OMan,OMan-1)
      CALL DENPADE(AlphaPad,DNM1,OMan-1,OMan-2)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                          

      WRITE(*,*) ' ANM PADE ZRHQR : RECHERCHE DES POLE DES PADE '       
      DD = DN
      EPSILONPAD=1.D-10
      RANGEMAX=100.0_dp * Apoly
      
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -      
! On souhaite avoir ici:
!   + Plus petite racine réelle positive
!   + pour test : racine réelle pos et neg pour vérif efficacité double defaut dans série

! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! POLE PADE ZRHQR
      CALL ANMPadeRACINESzrhqr(DD, OMan-1, RAC_R_POS_MIN,       &
                               NB_RAC_R_POS , TAB_RAC_R_POS  ,  &
                               NB_RAC_R_NEG , TAB_RAC_R_NEG     )
!----------------------------
                               
      LambdaCRIT = 0.0_dp
!       CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,LambdaCRIT,OMan-1,RAC_R_POS_MIN,DN,DELTA)
      CALL ANMSolPad(U0,lambda0,UMan,lambda,VSOLTMP,LambdaCRIT,OMan-1,RAC_R_POS_MIN,DN,DELTA)
      write (*,*) "*-*-*- ANM PADE ZRHQR ANMSolPad : RAC_R_POS_MIN=",RAC_R_POS_MIN," LracMin=",LambdaCRIT      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -      


!----------------------------
! POLE PADE BAIRSTOW
      CALL ANMPadeRACINES(DD, OMan-1, EPSILONPAD, RAC_BAI_MIN, RAC_MIN_C, DIVERGE, & 
                    RANGEMAX , NB_RAC, NB_RAC_POS, NB_RAC_C,RAC) !,               &
!       CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,LambdaCRIT,OMan-1,RAC_MIN,DD,DELTA)
!              write (*,*) "*-*-*- ANM PADE : RACBAI_MIN=",RAC_MIN," LracMin=",LambdaCRIT
                    
!----------------------------
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -      

!
!
! -> dL/da (RACMIN) if RACMIN is Acritical /!\:
!                  0 => LIMIT POINT  
!              NOT 0 => BIFURCATION
! -> AMXPAD  in [0,RACMIN]
!
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
     IF ( (NB_RAC_POS.GT.0).AND.(NB_RAC_R_POS.GT.0) ) THEN
       RAC_MIN = MIN(RAC_R_POS_MIN, RAC_BAI_MIN)
     ELSEIF ( (NB_RAC_POS.EQ.0).AND.(NB_RAC_R_POS.GT.0) ) THEN
       RAC_MIN = RAC_R_POS_MIN
     ELSEIF ( (NB_RAC_POS.GT.0).AND.(NB_RAC_R_POS.EQ.0) ) THEN
       RAC_MIN = RAC_BAI_MIN
     ENDIF
     
     IF(RAC_MIN.GT.0)  THEN
        RANGEMAX = RAC_MIN        
      ELSE
        RANGEMAX= 100.0_dp * Apoly
      ENDIF
      WRITE(*,*) 'ANM PADE : AmaxPad in [',Apoly,',',RANGEMAX,']'
      CRIT = 1.D-8
!     RAYON PADE DANS [Apoly , RAC_MIN ] 
!       CALL ANMPadeRoV(CRIT,VMan,DN,DNM1,lambda0,U0,lambda,Apoly,OMan,amaxPad,RANGEMAX, &
!                       U0Pad,L0Pad) 
      CALL ANMPadeRoV(CRIT,UMan,DN,DNM1,lambda0,U0,lambda,Apoly,OMan,amaxPad,RANGEMAX, &
                      U0Pad,L0Pad)                       
      WRITE(*,*) 'ANM PADE : AmaxPol = ',Apoly
      WRITE(*,*) 'ANM PADE : AmaxPad = ',amaxPad
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

!                     
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   
! TABLEAU VERIF POUR TEST
!     TEST RACINE ET LAMBDA CRITIQUE POUR DEBUG et pour voir
      if (NB_RAC_POS.GT.0._dp) then
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         write(*,*) "*-*-*- ANM PADE BAI    PADE's POLE "      
         ! VERIFIER POINT LIMITE
         !test lambda critique
         write(*,*) "*-*-*- ANM PADE BAI    lambda0 ",lambda0
         write(147,*) "# - - - - - - - - - BAI - - - - - - - - - - "
!          write(*,*) 'ANM PADE : lambda0 = ',lambda0
         DO I=1,NB_RAC_POS
           IF (RAC(I).GT.0._dp ) THEN
!              CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,Ltmp,OMan-1,RAC(I),DD,DELTA)
             CALL ANMSolPad(U0,lambda0,UMan,lambda,VSOLTMP,Ltmp,OMan-1,RAC(I),DD,DELTA)
             write (*,*) '*-*-*- ANM PADE BAI: LambdaPADEavant (',RAC(I),') = ',Ltmp
             write(147,905) Step,lambda0,Apoly,amaxPad,RAC(I),Ltmp
             CALL FLUSH(147)       
             
           ENDIF
         ENDDO  
      endif      
      
      if (NB_RAC_R_POS.GT.0._dp) then
         write(147,*) "# - - - - - - - - - ZRHQR POS - - - - - - - - - - "               
         write(*,*) "*- ZRHQR *- ZRHQR  *- ZRHQR  *- ZRHQR  *- ZRHQR  *- "
         write(*,*) "*-*-*- ANM PADE ZRHQR    PADE's POLE "      
         ! VERIFIER POINT LIMITE
         !test lambda critique
         write(*,*) "*-*-*- ANM PADE ZRHQR    lambda0 ",lambda0
!          write(*,*) 'ANM PADE : lambda0 = ',lambda0
         DO I=1,NB_RAC_POS
!              CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,Ltmp,OMan-1,TAB_RAC_R_POS(I),DD,DELTA)
             CALL ANMSolPad(U0,lambda0,UMan,lambda,VSOLTMP,Ltmp,OMan-1,TAB_RAC_R_POS(I),DD,DELTA)
             write (*,*) '*-*-*- ANM PADE ZRHQR POS : LambdaPADEavant (',TAB_RAC_R_POS(I),') = ',Ltmp
             write(147,905) Step,lambda0,Apoly,amaxPad,TAB_RAC_R_POS(I),Ltmp
             CALL FLUSH(147)       
             
         ENDDO
      endif       
         
      if (NB_RAC_R_NEG.GT.0._dp) then         
         write(147,*) "# - - - - - - - - - ZRHQR NEG - - - - - - - - - - "         
         DO I=1,NB_RAC_R_NEG
!              CALL ANMSolPad(U0,lambda0,VMan,lambda,VSOLTMP,Ltmp,OMan-1,TAB_RAC_R_NEG(I),DD,DELTA)
             CALL ANMSolPad(U0,lambda0,UMan,lambda,VSOLTMP,Ltmp,OMan-1,TAB_RAC_R_NEG(I),DD,DELTA)
             write (*,*) '*-*-*- ANM PADE ZRHQR NEG : LambdaPADEavant (',TAB_RAC_R_NEG(I),') = ',Ltmp
             write(147,905) Step,lambda0,Apoly,amaxPad,TAB_RAC_R_POS(I),Ltmp
             CALL FLUSH(147)                
         ENDDO           
      endif       
         
      write(*,*) "*- ZRHQR *- ZRHQR  *- ZRHQR  *- ZRHQR  *- ZRHQR  *- "
      write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"   
      
  905 FORMAT( 1X,I4,2X,G15.8,2X,G15.8,2X,G15.8,2X,'|',2X,G15.8,2X,G15.8) 
         
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

  END SUBROUTINE ANMPadeBifDetectZRHQR
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *      
  
  
  

! ORTHOG
  SUBROUTINE ANMOrthogGSm( V, AlphaPad, N )
!       SUBROUTINE ORTHOGS_as5(V,ALPHA,NDL,N)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%                                                                 %
!%  ORTHOGONALISATION DE GRAM-SCHMIDT (produit scalaire ordinaire) %
!%                                                                 %
!%  entree V      :  N vecteurs de dimension NDL                   %
!%                                                                 %
!%  sortie V      :  N vecteurs de dimension NDL orthonormaux      %
!%                                                                 %
!%     ALPHA      :  Coefficients d'orthogonalisation              %
!%                                                                 %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       DIMENSION V(NDL,1),AlphaPad(N,N)      
      REAL(KIND=dp), DIMENSION(:,:) :: V,AlphaPad
      INTEGER :: I,J,K,N
      
!.................//////INITIALISATION////..........................    
      AlphaPad = 0.0_dp
!      
      AlphaPad(1,1) = DSQRT(DOT_PRODUCT(V(:,1),V(:,1)))
!       WRITE(6,*) 'AlphaPad(1,1)=',AlphaPad(1,1) 

      V(:,1)=V(:,1)/AlphaPad(1,1)

!......//////pour chacun des vecteurs a partir du deuxieme////.....
      DO  I=2,N 
        DO  J=1,I-1   
          AlphaPad(I,J)=DOT_PRODUCT(V(:,I),V(:,J))
          AlphaPad(J,I) = AlphaPad(I,J)
!           WRITE(6,*) 'AlphaPad(',I,J,')=',AlphaPad(I,J) 
            V(:,I)=V(:,I)-AlphaPad(I,J)*V(:,J)
        ENDDO! 
        AlphaPad(I,I)=DSQRT(DOT_PRODUCT(V(:,I),V(:,I)))
!         WRITE(6,*) 'AlphaPad(',I,I,')=',AlphaPad(I,I) 
        V(:,I)=V(:,I)/AlphaPad(I,I)      
      ENDDO
!......//////impression des AlphaPad ////.....     
!     WRITE(6,50)
!     DO 15 I=1,N
!     DO 15 J=1,I
!15   WRITE(6,60) I,J,AlphaPad(I,J)
!50   FORMAT( ///,' ORTHOGONALISATION DE GRAM-SCHMIT ',/,
!    *            ' -------------------------------- ',//)
!60   FORMAT( ' AlphaPad(',I2,',',I2,')=',G15.8)     
      END SUBROUTINE ANMOrthogGSm
  
  
      SUBROUTINE DENPADE(AlphaPad,DEN,NORDRE,NPAD)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%                                                                 %
!%   COEEF Di POUR COEEF PADE 1 + a d1 + a^2 d2 +......            %
!%                                                                 %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER :: I,J,K,NORDRE,NPAD
      REAL(KIND=dp) :: AlphaPad(NORDRE,NORDRE),DEN(0:NPAD)
        
!     
      DEN = 0.0_dp
      DEN(0)=1.0_dp
      DO I=1,NPAD
        K=NPAD+1-I
        DO J=0,I-1
          DEN(I)=DEN(I)+ DEN(J)*AlphaPad(NPAD+1-J,K)
        ENDDO
        DEN(I)=DEN(I)/(-AlphaPad(K,K))
      ENDDO
!......//////impression des DEN ////.....     
!       WRITE(6,50)
!       DO I=0,NPAD-1
!         WRITE(6,60) I,DEN(I)
!       ENDDO
!  50   FORMAT( ///,' DENOMINATEUR DES APPROX. DE PADE ',/, &
!                  ' -------------------------------- ',//)
!  60   FORMAT( ' D(',I2,')=',G15.8) 
      END SUBROUTINE DENPADE
      
!---------------------------------------------------------------

!
      subroutine ANMPadeRACINES(P,NDEGRE,EPSILON2,RAC_MIN,RAC_MIN_C,DIVERGE, & 
                                AM,NB_RAC,NB_RAC_POS,NB_RAC_C,RAC )!,           &
! ! !                                 RACZRHQR,RAC_MIN_BIS,NB_RAC_ZRHQR)
! EXTERNE     
        INTEGER :: NDEGRE,NB_RAC,NB_RAC_POS,NB_RAC_C
        REAL(KIND=dp) :: P(0:NDEGRE), E1(NDEGRE)
        REAL(KIND=dp) :: EPSILON2,RAC_MIN,RAC_MIN_C,AM
        REAL(KIND=dp) :: RAC(NDEGRE),RACZRHQR(NDEGRE)
        LOGICAL :: DIVERGE
! INTERNE        
        REAL(KIND=dp) :: RAC_C(NDEGRE),PREEL(NDEGRE),PCOMP(NDEGRE),RACPOS(NDEGRE)
! ! !         REAL(KIND=dp) :: RAC_C_BIS(NDEGRE),rtr(NDEGRE),rti(NDEGRE),RACPOS_BIS(NDEGRE)
! ! !         REAL(KIND=dp) :: NB_RAC_POS_BIS, NB_RAC_C_BIS,NB_RAC_ZRHQR, RAC_MIN_BIS
        INTEGER :: I,K
        
        NB_RAC=0
        NB_RAC_POS=0
        NB_RAC_C=0
        RAC=0.0_dp
        
! ! !         RACZRHQR=0.0_dp
! ! !         NB_RAC_ZRHQR=0
        

        ! RECHERCHE DES RACINES
        call BAI(P,NDEGRE,EPSILON2,RAC,RAC_C,NB_RAC,NB_RAC_POS,NB_RAC_C,DIVERGE, &
                 AM,PREEL,PCOMP)

        ! TRAITEMENTS DES RACINES                 
        WRITE(*,*) 'ANM PADE RACINE : NB_RAC_POS = ',NB_RAC_POS
        WRITE(*,*) 'ANM PADE RACINE : NB_RAC_C   = ',NB_RAC_C        
        ! PLUS PETITE RACINE POSITIVE REELLE
        IF (NB_RAC_POS.GT.0) THEN
          K=1
          RACPOS=1D+20
          NB_RAC = NB_RAC_POS
          NB_RAC_POS = 0
          DO I=1,NB_RAC
            IF (RAC(I).GT.0._dp) THEN
              RACPOS(K)=RAC(I)
              NB_RAC_POS = NB_RAC_POS + 1
              K=K+1
            ENDIF
          ENDDO
          RAC_MIN = minval(RACPOS(1:K))
        WRITE(*,*) 'ANM PADE RACINE : RAC_MIN    = ',RAC_MIN
        ENDIF
! ! ! !         CALL TRI(PREEL,PCOMP,RAC,NB_RAC,NB_RAC_C,
! ! ! !      &  NB_RAC_POS,NDEGRE)
! ! ! !         if (nb_rac_pos.gt.0) then
! ! ! !         call  RACINE_MIN(RAC,NB_RAC_POS,RAC_MIN)
! ! ! !         else
! ! ! !         rac_min=10.d0*AM
! ! ! !         endif 
! ! ! !         
! ! ! !         if (nb_rac_C.gt.0) then
! ! ! !         call  RACINE_MIN(RAC_C,NB_RAC_C,RAC_MIN_C)
! ! ! !         else
! ! ! !         rac_min_C=0.d0
! ! ! !         endif 
! ! ! ! C------- je multiplie par 2 car a+-ib
! ! ! !         NB_RAC_C=2*NB_RAC_C
        IF (nb_rac_pos.gt.0) then 
! ! ! !            WRITE(6,*) 'RACINES REELLES : ',nb_rac_pos
           DO I=1,NB_RAC
             WRITE(6,*) 'REELLE(',I,') = ',RAC(I),RAC(I),0.0
!          WRITE(11,*) RAC(I),0.0
             PREEL(NB_RAC_C+I)=RAC(I)
             PCOMP(NB_RAC_C+I)=0.D0
             RAC_C(NB_RAC_C+I)=RAC(I)
           ENDDO 
! ! ! ! C        ELSE
! ! ! ! C           WRITE(6,*) 'QUE DES RACINES COMPLEXES'
        ENDIF
        

! ! ! !   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
! ! ! ! RACINE AVEC ZRHQR
! ! ! !   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 
! ! ! ! ! !   Preparation des donnees avant  appel
! ! !         DO I=1,NDEGRE
! ! !           E1(I) = P(I-1)
! ! !         ENDDO
! ! ! !       APPEL DE ZRHQR 
! ! !         CALL zrhqr(E1,NDEGRE-1,rtr,rti)
! ! ! !       TRI DES INFOS
! ! ! ! TABLEAU DES RACINES POUR TEST DOUBLE VUE : AVANT ET ARRIERE
! ! ! !    ReyCritique avant----------Step------------RecCritique apres
! ! !         DO I=1,NDEGRE
! ! !           IF ( ( ABS( rtr(I) ).GT.0._dp ).AND.( rti(I).LT.1D-15) ) THEN
! ! !             NB_RAC_ZRHQR=NB_RAC_ZRHQR + 1 
! ! !             RACZRHQR(NB_RAC_ZRHQR)
! ! !           ENDIF
! ! !         ENDDO
! ! ! !   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -        
! ! !         K=1
! ! !         RACPOS_BIS=1D+20
! ! !         NB_RAC_POS_BIS = 0
! ! !         DO I=1,NB_RAC_ZRHQR
! ! !             IF (RACZRHQR(I).GT.0._dp) THEN
! ! !               RACPOS_BIS(K)=RACZRHQR(I)
! ! !               NB_RAC_ZRHQR = NB_RAC_ZRHQR + 1
! ! !               K=K+1
! ! !             ENDIF
! ! !           ENDDO
! ! !           RAC_MIN_BIS = minval(RACPOS_BIS(1:K))        
! ! !         ENDDO
! ! ! !   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -                    
! ! !         
        END subroutine ANMPadeRACINES

! - - - - - - - - - - - - - - - - - - - - 
! Min Real POS Pole  - Check Almost Real zi<<1 but !=0
! Tab of sorted Real POS Poles
! Tab of sorted Real NEG Poles
!
      subroutine ANMPadeRACINESzrhqr(P, NDEGRE, RAC_R_POS_MIN,      &
                                     NB_RAC_R_POS , TAB_RAC_R_POS , &
                                     NB_RAC_R_NEG , TAB_RAC_R_NEG   )
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -                                       
! EXTERNE     
        INTEGER       :: NDEGRE
        REAL(KIND=dp) :: P(0:NDEGRE)
        REAL(KIND=dp) :: RAC_R_POS_MIN, RAC_R_NEG_MIN
        REAL(KIND=dp) :: TAB_RAC_R_POS(NDEGRE),TAB_RAC_R_NEG(NDEGRE)
        INTEGER       :: NB_RAC_R_POS, NB_RAC_R_NEG
! INTERNE        
        REAL(KIND=dp) :: E1(NDEGRE)
        REAL(KIND=dp) :: rtr(NDEGRE),rti(NDEGRE)
        REAL(KIND=dp) :: ZERO, SIGNTEST, TMP
        INTEGER       :: I , INDMIN, INDMAX
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -                 

        RAC_R_POS_MIN=0._dp
        RAC_R_NEG_MIN=0._dp
        
        ZERO = 0._dp
               
        NB_RAC_R_POS=0            !
        TAB_RAC_R_POS=0._dp       ! Sorted array of real positive pole : RAC_R_POS_MIN=TAB(1)
        NB_RAC_R_NEG=0            !
        TAB_RAC_R_NEG=0._dp       ! Sorted array of real negative pole : RAC_R_NEG_MIN=TAB(1)
        
! RACINE AVEC ZRHQR
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 
! ! !   Preparation des donnees avant  appel
        DO I=1,NDEGRE
          E1(I) = P(I-1)
        ENDDO
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -         
!       APPEL DE ZRHQR 
!       zrhqr  out  : 2 tab rtr rti        
!       Z(K) = rtr(K) + i rti(K)
        CALL zrhqr(E1,NDEGRE-1,rtr,rti)
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 
!       TRI DES INFOS
        DO I=1,NDEGRE-1
          WRITE(6,*) 'z(',I,') = ',rtr(I),' + i ',rti(I)
          IF (.NOT.(DABS(rti(i)).GT.ZERO)) THEN
            WRITE(6,*) '-- VERIF : zreal(',i,') = ',rtr(i),' + i ',rti(i)          
            ! -- TEST POS OR NEG --
            SIGNTEST = DABS(rtr(I)) * rtr(I)
            IF ( SIGNTEST.GT.ZERO ) THEN          
              NB_RAC_R_POS = NB_RAC_R_POS + 1
              TAB_RAC_R_POS(NB_RAC_R_POS) = rtr(I)
            ELSEIF ( SIGNTEST.LT.ZERO ) THEN
              NB_RAC_R_NEG = NB_RAC_R_NEG + 1
              TAB_RAC_R_NEG(NB_RAC_R_NEG) = rtr(I)            
            ENDIF ! END TEST SIGN
          ENDIF   ! END TEST REAL
        ENDDO
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -         
        IF (NB_RAC_R_POS==0) THEN
          WRITE(*,*) 'AUCUNE RACINE RELLE POSITIVE'
          RAC_R_POS_MIN=0.0_dp
        ELSE 
!    -    -    -    -    -    SORT REAL POS POLES  1 < 2  -    -    -    -    -    -    -        
! Algo : tri selection : echange (actuel , plus petit du tableau restant)
         DO I=1,NB_RAC_R_POS
           TMP = TAB_RAC_R_POS(I)                          ! TMP <- Tab(I)
           INDMIN = I+MINLOC(TAB_RAC_R_POS(I:NB_RAC_R_POS),1)  ! find indice of min in last part of array
           TAB_RAC_R_POS(I) = TAB_RAC_R_POS(INDMIN)        ! Tab(I) <- Min
           TAB_RAC_R_POS(INDMIN) = TMP                     ! Tab Min <- Tab(I) stored
         ENDDO
         RAC_R_POS_MIN = TAB_RAC_R_POS(1) 
        ENDIF
         
!    -    -    -    -    -    SORT REAL NEG POLES  |-1| < |-2|  -    -    -    -    -    -    -                
        IF (NB_RAC_R_NEG.GT.0) THEN
          DO I=1,NB_RAC_R_NEG
            TMP = TAB_RAC_R_NEG(I)                          ! TMP <- Tab(I)
            INDMAX = I+MAXLOC(TAB_RAC_R_NEG(I:NB_RAC_R_NEG),1)  ! find indice of min in last part of array
            TAB_RAC_R_NEG(I) = TAB_RAC_R_NEG(INDMAX)        ! Tab(I) <- Min
            TAB_RAC_R_NEG(INDMAX) = TMP                     ! Tab Min <- Tab(I) stored         
          ENDDO   
          RAC_R_NEG_MIN = TAB_RAC_R_NEG(1)        
         ENDIF
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - 

      END subroutine ANMPadeRACINESzrhqr
        
        
        

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! RECHERCHE DES RACINES D'UN POLYNOME D'ORDRE ELEVE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 ! ! !    rarines_zrhqr.f
! ! ! 
! ! ! COMM  Preparation des donnees avant  appel
! ! ! 
! ! !       DO I=1,NORDRE-1
! ! !           E1(I) = DEN(I-1)
! ! !        ENDDO
! ! ! COMM  APPEL DE ZRHQR       
! ! !        CALL zrhqr(E1,NORDRE-2,rtr,rti) 

! ! ! COMM  Subroutines necessaires au calcul des racines

      SUBROUTINE zrhqr(a,m,rtr,rti)
! C    USES balanc,hqr
! C    Find all the roots of a polynomial with real coefficients,    m+1
! C    i=1 a(i)xi-1, given the degree
! C    m and the coefficients a(1:m+1). The method is to construct an upper Hessenberg matrix
! c    whose eigenvalues are the desired roots, and then use the routines balanc and hqr. The
! c    real and imaginary parts of the roots are returned in rtr(1:m) and rti(1:m), respectively.      
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)       
!
      INTEGER :: m,MAXM
      REAL(KIND=dp) :: a(m+1),rtr(m),rti(m)
!
      REAL(KIND=dp) :: xr,xi
      INTEGER :: j,k            
      REAL(KIND=dp), ALLOCATABLE :: hess(:,:)

      MAXM=80      
      ALLOCATE( hess(MAXM,MAXM) ) 
      hess = 0.0_dp
      
! c      write(6,*) 'MAXM =',MAXM
! c      write(6,*) 'a(',m+1,')=',a(m+1)
      if (m.gt.MAXM.or.a(m+1).eq.0.) then
      write(6,*) 'dans zrhqr m=',m
      write(6,*) 'a(',m+1,') =',a(m+1)
      endif
      
      if (m.gt.MAXM.or.a(m+1).eq.0.) THEN
       WRITE(*,*) 'bad args in zrhqr'
      ELSE
         do 12 k=1,m    
            hess(1,k)=-a(m+1-k)/a(m+1)
            do 11 j=2,m
               hess(j,k)=0.d0
  11        CONTINUE
           if (k.ne.m) hess(k+1,k)=1.
  12        CONTINUE
         call balanc(hess,m,MAXM)   
         call hqr(hess,m,MAXM,rtr,rti)
         do 14 j=2,m   
            xr=rtr(j)
            xi=rti(j)
            do 13 k=j-1,1,-1
               if(rtr(k).le.xr)goto 1
                 rtr(k+1)=rtr(k)
                 rti(k+1)=rti(k)
   13       CONTINUE            
            k=0
   1    rtr(k+1)=xr
        rti(k+1)=xi
   14   CONTINUE
      ENDIF
        DEALLOCATE( hess ) 
        return
        END


      SUBROUTINE balanc(a,n,np)
      INTEGER :: n,np
      REAL(KIND=dp) :: a(np,np)      
      
      INTEGER :: i,j,last
      REAL(KIND=dp) :: c,f,g,r,s
      REAL(KIND=dp) :: RADIX , SQRDX
      RADIX=2.
      SQRDX=RADIX**2      
! Given an n by n matrix a stored in an array of physical dimensions np by np, this routine
! replaces it by a balanced matrix with identical eigenvalues. A symmetric matrix is already
! balanced and is unaffected by this procedure. The parameter RADIX should be the machine's
! floating-point radix.
!       INTEGER i,j,last
! REAL c,f,g,r,s
 1    continue
       last=1
       do 14 i=1,n  
       c=0.D0
       r=0.D0
       do 11 j=1,n
          if(j.ne.i)then
            c= c + dabs(a(j,i))
            r=r+abs(a(i,j))
           endif
  11   CONTINUE       
       if(c.ne.0..and.r.ne.0.) then    
          g=r/RADIX
          f=1.d0
          s=c+r
 2    if(c.lt.g)then    
            f=f*RADIX    
            c=c*SQRDX
      goto 2
      endif
      g=r*RADIX
  3    if(c.gt.g)then
         f=f/RADIX
         c=c/SQRDX
      goto 3
      endif
! ! c      if((c+r)/f.lt.0.95*s)then
      if((c+r)/f.lt.(95d-2)*s)then
      last=0
      g=1.d0/f
      do 12 j=1,n    
          a(i,j)=a(i,j)*g
 12   CONTINUE
      do 13 j=1,n
          a(j,i)=a(j,i)*f
 13   CONTINUE
      endif
      endif
 14   CONTINUE    
      if(last.eq.0)  goto 1
      return
      END

      SUBROUTINE hqr(a,n,np,wr,wi)
      INTEGER :: n,np      
      REAL(KIND=dp) :: a(np,np),wi(np),wr(np)
!  Finds all eigenvalues of an n by n upper Hessen
!  array. On input a can be exactly as output fro
!  The real and imaginary parts of the eigenvalues
      REAL(KIND=dp) :: anorm,p,q,r,s,t,u,v,w,x,y,z
      INTEGER :: i,j,k,l,m,last,its,nn
      
       anorm=0.d0
       do 12 i=1,n
          do 11 j=max(i-1,1),n
              anorm=anorm+dabs(a(i,j))
 11        CONTINUE 
 12    CONTINUE
       nn=n
       t=0.d0
 1    if(nn.ge.1) then
       its=0
2     do 13 l=nn,2,-1
        s=dabs(a(l-1,l-1))+dabs(a(l,l))
        if(s.eq.0.) s=anorm
        if(dabs(a(l,l-1))+s.eq.s) then
           a(l,l-1)=0.d0
        goto 3
        endif
 13   CONTINUE 
      l=1
3     x=a(nn,nn)
      if(l.eq.nn)then
         wr(nn)=x+t
         wi(nn)=0.d0
         nn=nn-1

      else
         y=a(nn-1,nn-1)
         w=a(nn,nn-1)*a(nn-1,nn)
         if(l.eq.nn-1)then    
            p=5.d-1*(y-x)
            q=p**2+w
            z=dsqrt(dabs(q))
            x=x+t
            if(q.ge.0.)then
               z=p+sign(z,p)
               wr(nn)=x+z
               wr(nn-1)=wr(nn)
               if(z.ne.0.)wr(nn)=x-w/z
                  wi(nn)=0.d0
                  wi(nn-1)=0.d0
               else    
                  wr(nn)=x+p
                  wr(nn-1)=wr(nn)
                  wi(nn)=z
                  wi(nn-1)=-z
               endif
               nn=nn-2
             else    
! ! c             if(its.eq.40)pause 'too many iterations in hqr'
             if(its.eq.40) RETURN
                         
             if(its.eq.10.or.its.eq.20)then    
             t=t+x
             do 14 i=1,nn
                   a(i,i)=a(i,i)-x
  14         CONTINUE 
             s=dabs(a(nn,nn-1))+dabs(a(nn-1,nn-2))
             x=75.d-2*s
             y=x
             w=-0.4375d0*s**2
             endif
           its=its+1
           do 15 m=nn-2,l,-1    
              z=a(m,m)    
              r=x-z
              s=y-z
              p=(r*s-w)/a(m+1,m)+a(m,m+1)    
              q=a(m+1,m+1)-z-r-s
              r=a(m+2,m+1)
              s=dabs(p)+dabs(q)+dabs(r)
              p=p/s
              q=q/s
              r=r/s
              if(m.eq.l)goto 4
              u=dabs(a(m,m-1))*(dabs(q)+abs(r))
              v=dabs(p)*(dabs(a(m-1,m-1))+dabs(z)+dabs(a(m+1,m+1)))
              if(u+v.eq.v)goto 4    
 15          CONTINUE 
 4           do 16 i=m+2,nn
               a(i,i-2)=0.d0
               if (i.ne.m+2) a(i,i-3)=0.d0
 16           CONTINUE 
            do 19 k=m,nn-1     
               if(k.ne.m)then    
                  p=a(k,k-1)    
                  q=a(k+1,k-1)
                  r=0.d0
                if(k.ne.nn-1)r=a(k+2,k-1)
                x=dabs(p)+dabs(q)+dabs(r)
                if(x.ne.0.)then
                   p=p/x    
                   q=q/x
                   r=r/x
                endif
               endif
               s=dsign(dsqrt(p**2+q**2+r**2),p)
               if(s.ne.0.)then
               if(k.eq.m)then
                  if(l.ne.m)a(k,k-1)=-a(k,k-1)
               else
                  a(k,k-1)=-s*x
               endif
               p=p+s    
               x=p/s
               y=q/s
               z=r/s
               q=q/p
               r=r/p
               do 17 j=k,nn    
                  p=a(k,j)+q*a(k+1,j)
                  if(k.ne.nn-1)then
                     p=p+r*a(k+2,j)
                     a(k+2,j)=a(k+2,j)-p*z
                  endif
                  a(k+1,j)=a(k+1,j)-p*y
                  a(k,j)=a(k,j)-p*x
 17            CONTINUE 
               do 18 i=l,min(nn,k+3)   
                  p=x*a(i,k)+y*a(i,k+1)
                  if(k.ne.nn-1)then
                      p=p+z*a(i,k+2)
                      a(i,k+2)=a(i,k+2)-p*r
                  endif
                  a(i,k+1)=a(i,k+1)-p*q
                  a(i,k)=a(i,k)-p
  18            CONTINUE 
              endif
  19       CONTINUE
           goto 2    
           endif
           endif
           goto 1    
           endif
           return
         END
           
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ! ! ! ! ! ! ! ! ! ! RACINEs PADE AVANT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE BAI(A,NDEGRE,EPSILON2,RAC,RAC_C,NB_RAC,NB_RAC_POS,NB_RAC_C,DIVERGE, & 
                       AM,PREEL,PCOMP)
! EXTERNE                            
        INTEGER :: NDEGRE,NB_RAC,NB_RAC_POS,NB_RAC_C
        REAL(KIND=dp) :: EPSILON2,AM
        REAL(KIND=dp) :: A(0:NDEGRE),B(0:NDEGRE),RAC(NDEGRE),RAC_C(NDEGRE),PREEL(1),PCOMP(1)
        LOGICAL :: DIVERGE
! INTERNE        
        LOGICAL :: REP
        INTEGER :: N,iter,ii,K
        REAL(KIND=dp) :: S,p
!
        N=NDEGRE
!-------test si poly > 2
        DIVERGE=.FALSE.
        if (N.gt.2) then
          S=AM+AM
          p=AM*AM
          DO WHILE(N.GT.2.AND.(.NOT.DIVERGE))
            S=AM+AM
            p=AM*AM
            REP=.FALSE.
            iter=0
            ii=0
            DO WHILE((.NOT.rep).and.ii.lt.500)
              iter=0
              S=(AM-dble(ii)*AM/500.d0)*2.D0
              p=(AM-dble(ii)*AM/500.D0)**2.d0-(AM-dble(ii)*AM/500.d0)
              DO WHILE((.NOT.REP).and.iter.le.100)
                CALL RAPHSON(A,B,N,S,P,EPSILON2,REP)
                iter=iter+1
              ENDDO
              ii=ii+1
            ENDDO
!       write (*,*) 'nb_iterations', iter
!-------test si n-r a converge: si oui on continue 
!-------sinon on sort 
            if (REP) then
              CALL RESOUT2(1.D0,-S,P,RAC,RAC_C,NB_RAC,NB_RAC_POS, NB_RAC_C,PREEL,PCOMP,NDEGRE)
              DO K=0,N-2
                A(K)=B(K+2)
              enddo
              N=N-2
            else
!       WRITE(*,*) 'NR n a pas trouve toutes les racines'
              DIVERGE=.TRUE.
            endif
!-------fin test n-r convergence        
          enddo   
        
          IF (N.eq.2.AND.(.NOT.DIVERGE)) THEN
            CALL RESOUT2(A(2),A(1),A(0),RAC,RAC_C,NB_RAC,NB_RAC_POS,NB_RAC_C,PREEL,PCOMP,NDEGRE)
          else
            IF (N.eq.1.AND.(.NOT.DIVERGE)) THEN
              CALL RESOUT1(A,NDEGRE,RAC,NB_RAC,NB_RAC_POS)
            ENDIF
          endif
!-------fin test si poly > 2
        ELSE
!-------si le polynome de depart < 3
          IF (N.eq.2) THEN
            CALL RESOUT2(A(2),A(1),A(0),RAC,RAC_C,NB_RAC,NB_RAC_POS,NB_RAC_C,PREEL,PCOMP,NDEGRE)
          else
            IF (N.eq.1) THEN
              CALL RESOUT1(A,NDEGRE,RAC,NB_RAC,NB_RAC_POS)
            ENDIF
          endif
!-------fin si le polynome de depart < 3
!-------fin si poly >2
        ENDIF   
        END SUBROUTINE BAI
        
        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SUBROUTINE RESOUT1(A,NDEGRE,RAC,NB_RAC,NB_RAC_POS)
!EXTERNE
        INTEGER :: NDEGRE,NB_RAC,NB_RAC_POS      
        REAL(KIND=dp) :: A(0:NDEGRE),RAC(NDEGRE)
!INTERNE
        REAL(KIND=dp) :: RAPP,ZERO

        RAPP=-A(0)/A(1)
        ZERO = 1d-15
        IF (ABS(A(1)).GT.ZERO) THEN
!         IF (RAPP.GE.0.D0) THEN 
          RAC(NB_RAC_POS+1)=RAPP
          NB_RAC_POS=NB_RAC_POS+1
!         ENDIF
          NB_RAC=NB_RAC+1
        ENDIF
        END SUBROUTINE RESOUT1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SUBROUTINE RESOUT2(A,B,C,RAC,RAC_C,NB_RAC,NB_RAC_POS,NB_RAC_C,PREEL,PCOMP,NDEGRE)
!EXTERNE
        INTEGER :: NDEGRE,NB_RAC,NB_RAC_POS,NB_RAC_C        
        REAL(KIND=dp) :: A,B,C,RAC(NDEGRE),RAC_C(NDEGRE),PREEL(NDEGRE),PCOMP(NDEGRE)
!INTERNE
        REAL(KIND=dp) :: DELTA,ZERO,RAC1,RAC2,xx
        
        DELTA=B*B-4.D0*A*C
        ZERO = 1d-15
        IF (DELTA.GE.0.D0) THEN
          IF (A.LE.ZERO) THEN
            RAC(NB_RAC_POS+1)= -1.D0 * C / B
            NB_RAC=NB_RAC+1
          ELSE
            RAC1= (-B+DSQRT(DELTA))/(2.D0*A)
            RAC2=(-B-DSQRT(DELTA))/(2.D0*A)
            RAC(NB_RAC_POS+1)=RAC1
            NB_RAC_POS=NB_RAC_POS+1
        
            RAC(NB_RAC_POS+1)=RAC2
            NB_RAC_POS=NB_RAC_POS+1
        
            NB_RAC=NB_RAC+2
          ENDIF
        ELSE 
          c=dabs(DSQRT(-DELTA)/(2.D0*A))
          xx=-B/(2.D0*A)
!       on tient pas compte des racines complexes
          RAC_C(NB_RAC_C+1)=dsqrt(xx*xx+c*c)
          PREEL(NB_RAC_C+1)= xx
          PCOMP(NB_RAC_C+1)= c
          PREEL(NB_RAC_C+2)= xx
          PCOMP(NB_RAC_C+2)= -1.D0 * c
          NB_RAC_C=NB_RAC_C+1
          NB_RAC=NB_RAC+2
        ENDIF       
        END SUBROUTINE RESOUT2
        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SUBROUTINE RAPHSON(A,B,N,S,P,EPSILON2,REP)
        
        INTEGER :: N
        REAL(KIND=dp) :: A(0:N),B(0:N),S,P,EPSILON2
        LOGICAL REP
        
        INTEGER :: K        
        REAL(KIND=dp) ::  DBDS(0:N),DBDP(0:N),AS,AP,val,DET

        AS=S
        AP=P
        B(N)=A(N)
        B(N-1)=A(N-1)+S*B(N)
        DO K=N-2,1,-1
          B(K)=A(K)+S*B(K+1)-P*B(K+2)
        ENDDO
        B(0)=A(0)-P*B(2)
        DBDS(N)=0.D0
        DBDS(N-1)=B(N)+S*DBDS(N)
        DO K=N-2,1,-1
          DBDS(K)=B(K+1)+S*DBDS(K+1)-P*DBDS(K+2)
        enddo
        DBDS(0)=-P*DBDS(2)
        
        DBDP(N)=0.D0
        DBDP(N-1)=S*DBDP(N)
        DO K=N-2,1,-1
          DBDP(K)=S*DBDP(K+1)-B(K+2)-P*DBDP(K+2)
        enddo
        DBDP(0)=-B(2)-P*DBDP(2)
        DET=DBDS(1)*DBDP(0)-DBDS(0)*DBDP(1)
        S=AS+(-DBDP(0)*B(1)+DBDP(1)*B(0))/DET
        P=AP+(DBDS(0)*B(1)-DBDS(1)*B(0))/DET
        val=((S-AS)*(S-AS)+(P-AP)*(P-AP))
        if (val.lt.EPSILON2) then
          REP=.true.
        endif
        END SUBROUTINE RAPHSON      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



! SOLUTION PADE
        SUBROUTINE ANMSolPad(U0,lambda0,UMan,lambda,Utmp,Ltmp,NPAD,A,DEN,DELTA)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%                                                                 %
!%   ON FORME LA SOLUTION  Lambda,U                                %
!%                                                                 %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        REAL(KIND=dp), DIMENSION(:,:) :: UMan
        REAL(KIND=dp), DIMENSION(:) :: lambda, U0 ,Utmp
        REAL(KIND=dp) :: A, lambda0, Ltmp, DEN(0:NPAD), DELTA
        INTEGER :: NPAD

        INTEGER :: I,J
        REAL(KIND=dp) :: AN(NPAD),P(0:NPAD),COEF
        
        Ltmp = lambda0 + A*lambda(1)
        Utmp   = U0 + A*UMan(:,1)

        AN(1) = A
        DO I=2,NPAD
          AN(I) = AN(I-1)*A
        ENDDO
        P(0) = DEN(0)
        DO I=1,NPAD-1
          P(I) = P(I-1) + AN(I)*DEN(I)
        ENDDO
        DELTA = P(NPAD-1)

        DO J=2,NPAD
          COEF = P(NPAD-J)/DELTA 
!       WRITE(6,50) A,J,COEF        
          COEF = COEF*AN(J)
          Ltmp = Ltmp + lambda(J)*COEF
          Utmp = Utmp + UMan(:,J)*COEF
        ENDDO
        END SUBROUTINE ANMSolPad

        
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
! DU/DA en A et DL/DA en A
!
        SUBROUTINE ANMDerivatePolySer(Uder,Lder,UMan,lambda,NORDRE,A)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! IN
        REAL(KIND=dp), DIMENSION(:,:) :: UMan  ! SERIE U
        REAL(KIND=dp), DIMENSION(:) :: lambda  ! SERIE Lambda
        INTEGER :: NORDRE                      ! SIZE OF SERIE
        REAL(KIND=dp) :: A                     ! RoV for Der evaluation
! OUT        
        REAL(KIND=dp), DIMENSION(:) :: Uder        
        REAL(KIND=dp) :: Lder
! LOCAL
        INTEGER :: I,J
        REAL(KIND=dp) :: PAM1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
        Lder = 0.0_dp
        Uder = 0.0_dp
        PAM1 = 1.0_dp                     !  A**(I-1)
        DO I=1,NORDRE
          Lder = Lder + I * PAM1 * lambda(I)        
          Uder = Uder + I * PAM1 * UMan(:,I)
          PAM1 = PAM1 * A            !  A**(I-1)
        ENDDO

        END SUBROUTINE ANMDerivatePolySer        
        
        
        
! ! ! ! ! !         
! ! ! ! ! !       SUBROUTINE VERIF_POINT_LIMITE(GAMA,Ua,NORDRE,NDL,A,DUa,
! ! ! ! ! !      &      DGAMA,VUTI,TROUVE,CALCULS)
! ! ! ! ! ! COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ! ! ! ! ! COMM %   Subroutine qui calcule la dérivée de LAMBDA, on forme 
! ! ! ! ! ! COMM %                   le Padé de cette dérivée, %
! ! ! ! ! ! COMM %             %
! ! ! ! ! ! COMM %     On se sert de Ua (et donc de sa derivee) juste pour 
! ! ! ! ! ! COMM %            calculer les coefficients des Pade                   %
! ! ! ! ! ! COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ! ! ! ! !       IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
! ! ! ! ! !       DIMENSION GAMA(1),Ua(NDL,*),DUa(NDL,1),DGAMA(1)
! ! ! ! ! !       DIMENSION ALPHA(NORDRE-1,NORDRE-1),VUTI(NDL,1),DU(NDL)
! ! ! ! ! !       DIMENSION DUa0(NDL),DEN(0:NORDRE-1)
! ! ! ! ! !       DIMENSION ALPHA2(NORDRE-2,NORDRE-2)
! ! ! ! ! !       LOGICAL TROUVE,CALCULS
! ! ! ! ! ! comm      1ere chose a faire : calcul des derivees        
! ! ! ! ! ! commm    ordre ZERO de la dérivée
! ! ! ! ! !       DGAMA0 = GAMA(1)
! ! ! ! ! !       DO I=1,NDL
! ! ! ! ! !          DUa0(I) = Ua(I,1)
! ! ! ! ! !       ENDDO
! ! ! ! ! ! COMM    Ordre superieur   
! ! ! ! ! !       DO I=2,NORDRE
! ! ! ! ! !          DGAMA(I-1) = I*1.d0*GAMA(I)
! ! ! ! ! !              DO J=1,NDL
! ! ! ! ! !                  DUa(J,I-1)= I*1.d0*Ua(J,I)
! ! ! ! ! !          ENDDO
! ! ! ! ! !       ENDDO
! ! ! ! ! !       CALL SHIFTD(DUa,VUTI,NDL*NORDRE-1)
! ! ! ! ! !       CALL ORTHOGS_as5(VUTI,ALPHA,NDL,NORDRE-1)
! ! ! ! ! !       CALL DENPADE_VERIF(ALPHA,DEN,NORDRE-1,NORDRE-2)                  
! ! ! ! ! ! 
! ! ! ! ! ! 
! ! ! ! ! ! comm      On verifie la valeur de  la derivee de LAMBDA
! ! ! ! ! ! 
! ! ! ! ! !       CALL SOLPAD2(DGAMA0,DUa0,DGAMA,DUa,DG,DU,
! ! ! ! ! !      *                  NDL,NORDRE-2,A,DEN,DELTA)         
! ! ! ! ! !       WRITE(6,*) 'Valeur de la derivee LAMBDA au pt singulier :',DG
! ! ! ! ! !                   
! ! ! ! ! !            DGAMA0 = DGAMA(1)
! ! ! ! ! !            DO I=1,NDL
! ! ! ! ! !            DUa0(I) = DUa(I,1)
! ! ! ! ! !            ENDDO
! ! ! ! ! ! COMM    Ordre superieur   
! ! ! ! ! !            DO I=2,NORDRE-1
! ! ! ! ! !              DGAMA(I-1) = I*1.d0*DGAMA(I)
! ! ! ! ! !                  DO J=1,NDL
! ! ! ! ! !                         DUa(J,I-1)= I*1.d0*DUa(J,I)
! ! ! ! ! !              ENDDO
! ! ! ! ! !            ENDDO
! ! ! ! ! !            CALL SHIFTD(DUa,VUTI,NDL*NORDRE-2)
! ! ! ! ! !            CALL ORTHOGS_as5(VUTI,ALPHA2,NDL,NORDRE-2)
! ! ! ! ! !            CALL DENPADE_VERIF(ALPHA2,DEN,NORDRE-2,NORDRE-3)  
! ! ! ! ! ! 
! ! ! ! ! !               CALL SOLPAD2(DGAMA0,DUa0,DGAMA,DUa,DDG,DU,
! ! ! ! ! !      *                  NDL,NORDRE-3,A,DEN,DELTA)         
! ! ! ! ! !       WRITE(6,*) 'Valeur de la derivee seconde :',DDG
! ! ! ! ! !       ZERO = 1.D-5
! ! ! ! ! !       IF (DABS(DG).LT.ZERO) THEN
! ! ! ! ! !              WRITE(6,*) 'Le POINT SINGULIER EST UN POINT LIMITE'
! ! ! ! ! !          TROUVE = .FALSE.
! ! ! ! ! !          CALCULS = .TRUE.
! ! ! ! ! !       ELSE
! ! ! ! ! ! !          WRITE(6,*) 'Le POINT SINGULIER EST UN POINT de BIFURCATION'            
! ! ! ! ! !              WRITE(6,*) "Le POINT SINGULIER N'EST PAS UN POINT LIMITE"
! ! ! ! ! !       ENDIF
! ! ! ! ! !       END               
! ! ! ! ! !         
! ! ! !       


!> RANGE OF VALIDITY OF RATIONNAL REPRESENTATION
      SUBROUTINE ANMPadeRoV(CRIT,VMan,DN,DNM1,lambda0,U0,lambda,AmaxPoly,OMan,amaxPad,RAC_MIN,  &
                            U0Pad, L0Pad) 

      ! EXTERNE
      INTEGER :: OMan          
      REAL(KIND=dp) :: CRIT, amaxPad, U0(:), lambda0, lambda(:), AmaxPoly, RAC_MIN
      REAL(KIND=dp) :: VMan(:,:)
      REAL(KIND=dp) :: DN(OMan),DNM1(OMan-1)
      REAL(KIND=dp) :: U0Pad(:), L0Pad
! INTERNE
      INTEGER :: NDL
      REAL(KIND=dp) :: Ltmp1,Ltmp2, XNORP, XNOR1,RAP, diff1, diff2 , borninf, bornsup , ZERO,DELTA
      REAL(KIND=dp), ALLOCATABLE :: UPN(:),UTRA(:),UTRA2(:)
      LOGICAL :: fin_dicho

      NDL=SIZE(VMan(:,1))
      ALLOCATE( UPN(NDL) , UTRA(NDL) )
      
      ZERO = 1.d-16
      borninf=AmaxPoly
      bornsup=RAC_MIN
      

      fin_dicho=.false.
      do while(.NOT.fin_dicho)
        amaxPad=(borninf+bornsup)/2.d0
        write(*,*) 'ANM PADE RoV : Atest = ',amaxPad
! - - -
!       PN(Atest)
        CALL ANMSolPad(U0,lambda0,VMan,lambda,UPN,Ltmp1,OMan-1,amaxPad,DN,DELTA)
        UTRA2=UPN-U0
        XNORP=DSQRT(DOT_PRODUCT(UTRA2,UTRA2))    
        write(*,*) 'ANM PADE RoV : XNORP=',XNORP
! - - -        
!       PN-1(Atest)        
        CALL ANMSolPad(U0,lambda0,VMan,lambda,UTRA,Ltmp2,OMan-2,amaxPad,DNM1,DELTA)                
!       || PN(Atest) - PN-1(Atest) ||
        UTRA = UPN - UTRA        
        XNOR1 = DSQRT(DOT_PRODUCT(UTRA,UTRA))
        write(*,*) 'ANM PADE RoV : XNOR1=',XNOR1        
! - - -
!        || PN(Atest) - PN-1(Atest) || / || PN(Atest) ||
        RAP=XNOR1/XNORP
        write(*,*) 'ANM PADE RoV : RAP = XNOR1/XNORP=',RAP        
! - - -        
        IF (XNORP.LT.ZERO) THEN
          fin_dicho=.true.
        ENDIF

!         if (RAP.lt.(CRIT/10.d0)) then
!           borninf=amaxPad
!         endif
!         if (RAP.gt.CRIT) then
!           bornsup=amaxPad
!         endif
        if (RAP.lt.CRIT) then
          write(*,*) 'ANM PADE RoV : DICHO :  - - - - - - - - - - - - - A DROITE'
          borninf=amaxPad
        else 
          write(*,*) 'ANM PADE RoV : DICHO :  - - - - - - - - - - - - - A GAUCHE'        
          bornsup=amaxPad
        endif
        
        diff1=dabs(amaxPad-AmaxPoly)
        diff2=dabs(amaxPad-RAC_MIN)
        if ((RAP.le.CRIT.and.RAP.ge.(CRIT/10.d0)).or.(diff1.lt.10.d-6).or.(diff2.lt.10.d-6) ) then
          fin_dicho=.true.
        endif
      enddo
      
      U0Pad = UPN
      L0Pad = Ltmp1
      
      DEALLOCATE( UPN, UTRA )      
      END SUBROUTINE ANMPadeRoV

      
      
      
      
      
      
      

!------------------------------------------------------------------------------
! BIFURCATION DETECTION
!    ___           _                     _       _ ____   ___  _ ____  
!   / __\___   ___| |__   /\/\   ___  __| | __ _| |___ \ / _ \/ |___ \ 
!  / /  / _ \ / __| '_ \ /    \ / _ \/ _` |/ _` | | __) | | | | | __) |
! / /__| (_) | (__| | | / /\/\ \  __/ (_| | (_| | |/ __/| |_| | |/ __/ 
! \____/\___/ \___|_| |_\/    \/\___|\__,_|\__,_|_|_____|\___/|_|_____| 
! COCHELIN MEDALE 2012
!------------------------------------------------------------------------------
 SUBROUTINE BIFdetectCritCochMed( VMan, ordreman, NoPressure, calcalphabif, sumcrit, &
                                                 epsilon1, epsilon2, alphaBIF,       &
                                                 Lambda )
 !------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
    IMPLICIT NONE
     REAL(KIND=dp), DIMENSION(:,:) :: VMan
     REAL(KIND=dp), DIMENSION(:) :: calcalphabif, NoPressure, sumcrit, Lambda
     REAL(KIND=dp) :: epsilon1, epsilon2, alphaBIF,alphaBIFbis
     INTEGER :: ordreman
!------------------------------------------------------------------------------
     REAL(KIND=dp), ALLOCATABLE :: VpOrtho(:,:)
     REAL(KIND=dp) :: power, projvect, normVn, sum1, sum2 , rapp
     REAL(KIND=dp) :: normVpOrtho, normVp
     REAL(KIND=dp) :: normVnBis,Lbissq,projvectB,calcalphabifBis(size(calcalphabif)),Lbis
     INTEGER :: ordretemp, n, i, j
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------ 
      n=SIZE(VMan(:,1))
      ALLOCATE( VpOrtho(n,3) )
      alphaBIF = 0.0_dp
      
            
!     CHOIX TECHNIQUE : + AlphaBif calculé avec U = (V,0)  de (V, Pression) 
!                         avec Lambda, sans pression
!                         -> on a vérifié que Lambda subit le défaut aussi.
!                       + Astuce : On pourrait ajouter Lambda à la place de la derniere pression
!                         afin de calculer alphaBif avec (V,Lambda)
!
!     Norme L2 de Un 
      normVn = DSQRT( DOT_PRODUCT( VMan(:,ordreman) , VMan(:,ordreman) ) )
! TEST AVEC LAMBDA      
      Lbis = Lambda(ordreman)*Lambda(ordreman)
      normVnBis = DSQRT( DOT_PRODUCT( VMan(:,ordreman) , VMan(:,ordreman) ) + Lbis)
      
      calcalphabif = 0.0_dp
!     Calcul de alpha_p et Up_ortho
      DO i=ordreman-3, ordreman-1
        j = i - ( ordreman-3 ) + 1
        projvect = DOT_PRODUCT( VMan(:,i), VMan(:,ordreman) )        
        calcalphabif(j) = projvect / (normVn*normVn)
        VpOrtho(:,j) = VMan(:,i) - calcalphabif(j)*VMan(:,ordreman)
! TEST AVEC LAMBDA        
        Lbis = Lambda(i)*Lambda(ordreman)      
        projvectB = DOT_PRODUCT( VMan(:,i), VMan(:,ordreman) ) + Lbis
        calcalphabifBis(j) = projvectB / (normVnBis*normVnBis)        
      END DO

      do i=1,3
        print*,'alpha(',i,')= ',calcalphabif(i),calcalphabifBis(i)
      end do
      
!     Vérification du critère de détection (Cochelin Medale 2013 )

!     ERREUR RAISON
      sum1 = 0.0_dp
      DO i=ordreman-3, ordreman-2
        j = i - ( ordreman-3 ) + 1
        power = 1._dp/( ordreman - i )
        rapp =  DABS( calcalphabif( j ) )**power /  DABS( calcalphabif( 3 ) ) 
        sum1 = sum1 + ( rapp - 1.0_dp )**2
       END DO

!     COLINEAIRE       
      sum2 = 0.0_dp      
      DO i=ordreman-3,ordreman-1
        j = i - ( ordreman-3 ) + 1
        normVp = DSQRT( DOT_PRODUCT( VMan(:,i),VMan(:,i) ))
        normVpOrtho = DSQRT( DOT_PRODUCT(   VpOrtho(:,j), VpOrtho(:,j) ))
        sum2 = sum2 + normVpOrtho / normVp
      END DO
        
       print*,'BIFdetectCritCochMed : sum1 Alpha = ',sum1
       print*,'BIFdetectCritCochMed : sum2 Colin = ',sum2
!      alphaBIF = 0.0_dp
!      print*,'alpha n-1 = ', calcalphabif(3)
      IF ( sum1 < epsilon1 .AND. sum2 < epsilon2 ) THEN
!       Lambda ( a = alpha3 ) = Lambda critique       
        alphaBIF = calcalphabif( 3 )
      ELSE
        alphaBIF = 0.0_dp
      END IF
      
      sumcrit = 0.0_dp
      sumcrit( 1 ) = sum1
      sumcrit( 2 ) = sum2
      DEALLOCATE( VpOrtho )
 END SUBROUTINE BIFdetectCritCochMed

 
 
 
!------------------------------------------------------------------------------
  SUBROUTINE BIFCMserieprop( VMan, Lambda, ordreman, alphaBIF , VManProp , LambdaProp,Vt1orth )
! EQUATION 22 bis  : U chapeau :  Cochelin Médale - JCP2013
     IMPLICIT NONE
     REAL(KIND=dp), DIMENSION(:,:) :: VMan , VManProp
     REAL(KIND=dp) :: alphaBIF ,LambdaProp(:),Lambda(:),Vt1orth(:)
     INTEGER :: ordreman
!
     REAL(KIND=dp) :: poww,xnorm
     INTEGER :: I,J,K
!
!       WRITE(*,*) "chapeau"
      VManProp = 0.0_dp
      LambdaProp = 0.0_dp
      DO I=1, ordreman-1
        poww = ordreman - I
        VManProp(:,I) = VMan(:,I) - VMan(:,ordreman) * alphaBIF**poww
        LambdaProp(I) = Lambda(I) - Lambda(ordreman) * alphaBIF**poww
      END DO
!       write(*,*) "BIFCMserieprop - "
      xnorm = DSQRT(DOT_PRODUCT(VMan(:,ordreman),VMan(:,ordreman)))
      Vt1orth = VMan(:,ordreman) / xnorm      
!       write(*,*) "BIFCMserieprop - xnorm=",xnorm
!       DO i=1,size(VMan(:,ordreman))
!         write(*,*) "BIFCMserieprop - VMan(",i,",",ordreman,")=",VMan(I,ordreman)
!         Ut1orth(I) = VMan(I,ordreman) / xnorm
!       ENDDO

 END SUBROUTINE BIFCMserieprop
!------------------------------------------------------------------------------

 
 
 
  
!          _    _  _____ __  __ ______ _   _ _______ ______ _____  
!     /\  | |  | |/ ____|  \/  |  ____| \ | |__   __|  ____|  __ \ 
!    /  \ | |  | | |  __| \  / | |__  |  \| |  | |  | |__  | |  | |
!   / /\ \| |  | | | |_ | |\/| |  __| | . ` |  | |  |  __| | |  | |
!  / ____ \ |__| | |__| | |  | | |____| |\  |  | |  | |____| |__| |
! /_/    \_\____/ \_____|_|  |_|______|_| \_|  |_|  |______|_____/  
!  
! Search for W and LeftMode in case of Simple BIFURCATION  
!
! | Ltc V |  W    F
! |       |     = 
! | VT  0 |  k    0
!
! | LtcT V |  PSI   0
! |        |      = 
! | VT   0 |  k     1
!
! CollectionMatrix X  CollectionSolution  = CollectionSolution
! 
!------------------------------------------------------------------------------

!>  This subroutine will solve the system with some linear restriction.
  SUBROUTINE BIF1SolveAugSystem( StiffMatrix, ForceVector, Solution,     &
                                 BifMode, CondVal, LEFT,                 &
                                 Norm, DOFs, Solver,MUMPSFICH,           &
                                 CollectionMatrix) 
!---------------------------------------------------------------------------------------------
    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
! - - - - - - - - - - - - - - -     
    USE HomeMadeSolvers
! - - - - - - - - - - - - - - -     
!------------------------------------------------------------------------------  
  IMPLICIT NONE
  TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 
                                         !< The restriction matrix is assumed to be in the EMatrix-field
  REAL(KIND=dp),TARGET :: ForceVector(:)        !< The right hand side of the linear equation
  REAL(KIND=dp),TARGET :: Solution(:)           !< Previous solution as input, new solution as output.
  REAL(KIND=dp) :: Norm                  !< The L2 norm of the solution.
  INTEGER :: DOFs                        !< Number of degrees of freedom of the equation.
  TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
!------------------------------------------------------------------------------
  REAL(KIND=dp),TARGET  :: BifMode(:)    ! Mode Bif
  REAL(KIND=dp)         :: CondVal        ! Condition : <Vec,W> = 0 ou <LeftMode,Vec> = 1
  LOGICAL               :: LEFT          ! Transpose : W(Flase) or LeftMode(True)


!------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: CollectionMatrix, RestMatrix, AddMatrix, &
       RestMatrixTranspose
  REAL(KIND=dp), POINTER ,CONTIGUOUS :: CollectionVector(:), RestVector(:),&
                 MultiplierValues(:),AddVector(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: CollectionSolution(:)
  INTEGER, ALLOCATABLE :: TmpRow(:)
  INTEGER :: NumberOfRows, NumberOfValues, MultiplierDOFs, istat
  INTEGER :: i, j, k, l
  TYPE(Variable_t), POINTER :: MultVar
  REAL(KIND=dp) :: scl,restv
  LOGICAL :: Found, ExportMultiplier, NotExplicit, Refactorize, EnforceDirichlet
  SAVE MultiplierValues
  CHARACTER(LEN=MAX_NAME_LEN) :: MultiplierName
  INTEGER :: MUMPSFICH
  REAL(KIND=dp) :: BMJ,TIC
!------------------------------------------------------------------------------
  
  
!!!!!! Commented for SysAug keep for branch  
!   CollectionMatrix => StiffMatrix % CollectionMatrix
!   IF(.NOT.ASSOCIATED(CollectionMatrix)) THEN
!     CollectionMatrix => AllocateMatrix()
!     CollectionMatrix % FORMAT = MATRIX_LIST
!   ELSE
!     DEALLOCATE(CollectionMatrix % RHS)
!     CollectionMatrix % Values = 0.0_dp
!   END IF
  
  
! !   Refactorize = .TRUE.
  
  NumberOfRows = StiffMatrix % NumberOfRows
  NumberOfRows = NumberOfRows + 1               ! CAS DE LA BIF SIMPLE : AUGMENTE d'UN VECT

  ALLOCATE( CollectionMatrix % RHS( NumberOfRows ), &
            CollectionSolution( NumberOfRows ),     &
            STAT = istat )  
  
!------------------------------------------------------------------------------  
  CollectionVector => CollectionMatrix % RHS
  CollectionVector = 0.0_dp
  CollectionSolution = 0.0_dp

!---------------------------------------------------------------------------------------------
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Put the RestMatrix to lower part of CollectionMatrix
! Enforce Exact Dirichlet BCs? StiffMatrix % ConstrainedDOF
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       RestMatrix => NULL()
! !       RestMatrix => StiffMatrix % ConstraintMatrix
!       RestMatrix => BifMode ! NE PEUX PAS FONCTIONNER CAR BIFMODe!=MATRIX
! ! ! ! AddToMatrixElement ( MAT,i,j,value )  (Matrix_list or CRS )
!       DO i=RestMatrix % NumberOfRows,1,-1         ! I = 1 sur un mode!

!!! IMPROVE : DO ONLY ONCE IF TOO LONG BUT HOW TO DO IT IF ALLREDAY FACTORIZED???
!! AddToMatrixElement TOO LONG to be used like that??
!! AddToMatrixElement for BifMode too long??
!!  
!!
!! CRS_AddToMatrixElement : k = CRS_Search for Values(k) = Values(k) + VALUE
      write(6,*) ' FORMAT = MATRIX_CRS //  MATRIX_LIST',CollectionMatrix % FORMAT
!                               CollectionMatrix ==> MATRIX_LIST = 4
!          CALL List_AddToMatrixElement( A % ListMatrix, i, j, VALUE )
!      Entry => List_GetMatrixIndex(List,k1,k2)

!        Entry % Value = Entry % Value + Value
      WRITE(6,*) '- BIF1SolveAugSystem - CollectionMatrix assembly : BifMode PARTS'
      CALL FLUSH(6)
      TIC = CPUTime()
      
        k=StiffMatrix % NumberOfRows              ! K = NDL
!         k=k+i                                     ! K = NDL + 1
        CALL AddToMatrixElement( CollectionMatrix, k+1 , k+1 , 0.0_dp )
!         DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIF1SolveAugSystem - CollectionMatrix assembly : BifMode  k+1 , k+1 TIME =',TIC
      CALL FLUSH(6)
      
      TIC = CPUTime()
        DO j=1,k
!         Mat(j,k) = BifMode(j)
!           CALL AddToMatrixElement( CollectionMatrix, RestMatrix % Cols(j) , k ,      &
!                                    RestMatrix % Values(j) )  
          BMJ=BifMode(j)
          CALL AddToMatrixElement( CollectionMatrix,   j , k+1 ,BMJ  )  
        END DO
      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIF1SolveAugSystem - CollectionMatrix assembly : BifMode j ,k+1,BMJ   TIME =',TIC
      CALL FLUSH(6)
      
      TIC = CPUTime()
!         DO j=1,k
        DO j=k,1,-1   !reverse because of allocation copy in background? List_GetMatrixIndex > List_EnlargeMatrix
!         Mat(k,j) = BifMode(j)   
          BMJ=BifMode(j)
          CALL AddToMatrixElement( CollectionMatrix, k+1 ,   j , BMJ )
        END DO
!       END DO
      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIF1SolveAugSystem - CollectionMatrix assembly : BifMode  k+1,j, BMJ    TIME =',TIC
      CALL FLUSH(6)
          
      CollectionVector(k+1) = CondVal
    

! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Put the StiffMatrix to upper part of CollectionMatrix
! Lt : W
! Lt Transpose LeftBifMode : POINT => CRS_Transpose( StiffMatrix ) ou i,j -> j,i
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,*) '- BIF1SolveAugSystem - CollectionMatrix assembly : Tangent Op PARTS'
      CALL FLUSH(6)
      TIC = CPUTime()
      
      IF (LEFT.EQV..FALSE.) THEN   
        DO i=StiffMatrix % NumberOfRows,1,-1      
          DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1
            CALL AddToMatrixElement( CollectionMatrix, i, StiffMatrix % Cols(j), StiffMatrix % Values(j) )
            CollectionVector(i) =  ForceVector(i)
          END DO
        END DO        
      ELSE
      ! TRANSPOSE BY HAND : i,j -> j,i
        DO i=StiffMatrix % NumberOfRows,1,-1
          DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1          
            CALL AddToMatrixElement( CollectionMatrix, StiffMatrix % Cols(j), i, StiffMatrix % Values(j) ) 
          END DO
        END DO
      ENDIF            
      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIF1SolveAugSystem - CollectionMatrix assembly : StiffMatrix   TIME =',TIC
      CALL FLUSH(6)      
      WRITE(6,*) '- BIF1SolveAugSystem - List_toCRSMatrix'
      CALL FLUSH(6)
      CALL List_toCRSMatrix(CollectionMatrix)

! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Assign values to CollectionVector
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
      j = StiffMatrix % NumberOfRows  
      CollectionSolution(1:j) = Solution(1:j)
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve the Collection-system 
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! TODO  : SOLVE WITH MUMPS
!       CALL SolveLinearSystem( CollectionMatrix, CollectionVector, &
!                               CollectionSolution, Norm, DOFs, Solver, StiffMatrix )
      WRITE(6,*) '- BIF1SolveAugSystem - HMumps_SolveLinearSystem'
      CALL FLUSH(6)
      
      
      ! LEFT      : CollectionMatrix % mumpsIDL % ICNTL(9) = 2
      ! Classic W : CollectionMatrix % mumpsIDL % ICNTL(9) = 1
! ! ! ! ICNTL(9) computes the solution using A or AT  augsys
! ! ! ! Phase: accessed by the host during the solve phase.
! ! ! ! Possible values :
! ! ! ! 1 : AX = B is solved.
! ! ! ! = 1 : AT X = B is solved.
! ! ! ! Default value: 1
! ! ! ! Related parameters: ICNTL(10), ICNTL(11), ICNTL(21), ICNTL(32)
! ! ! ! Remarks: when a forward elimination is performed during the factorization (see ICNTL(32))
! ! ! ! only ICNTL(9)=1 is allowed.
    
          

      CALL HMumps_SolveLinearSystem( CollectionMatrix, CollectionSolution, CollectionVector, &
                                     Solver , MUMPSFICH )                              
                              
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Separate the solution from CollectionSolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Solution = 0.0_dp
      i = 1
      j = StiffMatrix % NumberOfRows
      Solution(i:j) = CollectionSolution(i:j)
      restv=CollectionSolution(j+1)
      write(*,*) "BIF1SolveAugSystem - RestVect=",restv
!       WRITE(*,*) 'ASSOCIATED(CollectionMatrix) in= ',ASSOCIATED(CollectionMatrix)
      DEALLOCATE(CollectionSolution)
      CALL Info( 'BIF1SolveAugSystem', 'All done', Level=5 )

  END SUBROUTINE BIF1SolveAugSystem
    
!---------------------------------------------------------------------------------------------




!---------------------------------------------------------------------------------------------
!>  This subroutine will solve the system with some linear restriction.
!> Use Mumps for Ax=b and Atx=b 
!> 1 Assemble Aug Sys
!> 2 Solver either At or A
      ! LEFT      : CollectionMatrix % mumpsIDL % ICNTL(9) = 2
      ! Classic W : CollectionMatrix % mumpsIDL % ICNTL(9) = 1
! ! ! ! ICNTL(9) computes the solution using A or AT  augsys
! ! ! ! Phase: accessed by the host during the solve phase.
! ! ! ! Possible values :
! ! ! ! 1 : AX = B is solved.
! ! ! ! = 1 : AT X = B is solved.
! ! ! ! Default value: 1
! ! ! ! Related parameters: ICNTL(10), ICNTL(11), ICNTL(21), ICNTL(32)
! ! ! ! Remarks: when a forward elimination is performed during the factorization (see ICNTL(32))
! ! ! ! only ICNTL(9)=1 is allowed.
  SUBROUTINE BIFStaDimKer1ComputePsiW( StiffMatrix, ForceVector, BifMode,  &
                                 Norm, DOFs, Solver,MUMPSFICH,           &
                                 CollectionMatrix, abePsiL, abeW ) 
!---------------------------------------------------------------------------------------------
    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
! - - - - - - - - - - - - - - -     
    USE HomeMadeSolvers
! - - - - - - - - - - - - - - -     
!------------------------------------------------------------------------------  
  IMPLICIT NONE
  TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 
                                         !< The restriction matrix is assumed to be in the EMatrix-field
  REAL(KIND=dp),TARGET :: ForceVector(:)        !< The right hand side of the linear equation
  REAL(KIND=dp) :: Norm                  !< The L2 norm of the solution.
  INTEGER :: DOFs                        !< Number of degrees of freedom of the equation.
  TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
!------------------------------------------------------------------------------
  REAL(KIND=dp),TARGET  :: BifMode(:)    ! Mode Bif
!------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: CollectionMatrix, RestMatrix, AddMatrix, &
       RestMatrixTranspose
  REAL(KIND=dp), POINTER ,CONTIGUOUS :: CollectionVector(:), RestVector(:),&
                 MultiplierValues(:),AddVector(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: CollectionSolution(:)
  INTEGER, ALLOCATABLE :: TmpRow(:)
  INTEGER :: NumberOfRows, NumberOfValues, MultiplierDOFs, istat
  INTEGER :: i, j, k, l
  TYPE(Variable_t), POINTER :: MultVar
  REAL(KIND=dp) :: scl,restv
  LOGICAL :: Found, ExportMultiplier, NotExplicit, Refactorize, EnforceDirichlet
  SAVE MultiplierValues
  CHARACTER(LEN=MAX_NAME_LEN) :: MultiplierName
  INTEGER :: MUMPSFICH
  REAL(KIND=dp) :: BMJ,TIC
!------------------------------------------------------------------------------
  REAL(KIND=dp),TARGET  :: abePsiL(:),  abeW(:)   ! OUTPUT VECTORS
!------------------------------------------------------------------------------

  
  
  NumberOfRows = StiffMatrix % NumberOfRows
  NumberOfRows = NumberOfRows + 1               ! CAS DE LA BIF SIMPLE : AUGMENTE d'UN VECT

  ALLOCATE( CollectionMatrix % RHS( NumberOfRows ), &
            CollectionSolution( NumberOfRows ),     &
            STAT = istat )  
  
!------------------------------------------------------------------------------  
  CollectionVector => CollectionMatrix % RHS
  CollectionVector = 0.0_dp
  CollectionSolution = 0.0_dp

!---------------------------------------------------------------------------------------------
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - CollectionMatrix assembly : BifMode PARTS'
      CALL FLUSH(6)
      TIC = CPUTime()
      
      k=StiffMatrix % NumberOfRows              ! K = NDL
!         k=k+i                                     ! K = NDL + 1
      CALL AddToMatrixElement( CollectionMatrix, k+1 , k+1 , 0.0_dp )
!         DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - CollectionMatrix assembly : BifMode  k+1 , k+1 TIME =',TIC
      CALL FLUSH(6)
      
      TIC = CPUTime()
      DO j=1,k
!         Mat(j,k) = BifMode(j)
!           CALL AddToMatrixElement( CollectionMatrix, RestMatrix % Cols(j) , k ,      &
!                                    RestMatrix % Values(j) )  
        BMJ=BifMode(j)
        CALL AddToMatrixElement( CollectionMatrix,   j , k+1 ,BMJ  )  
      END DO
      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - CollectionMatrix assembly : BifMode j ,k+1,BMJ   TIME =',TIC
      CALL FLUSH(6)
      
      TIC = CPUTime()
!         DO j=1,k
      DO j=k,1,-1   !reverse because of allocation copy in background? List_GetMatrixIndex > List_EnlargeMatrix
!         Mat(k,j) = BifMode(j)   
        BMJ=BifMode(j)
        CALL AddToMatrixElement( CollectionMatrix, k+1 ,   j , BMJ )
      END DO
!       END DO
      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - CollectionMatrix assembly : BifMode  k+1,j, BMJ    TIME =',TIC
      CALL FLUSH(6)
          
    

! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Put the StiffMatrix to upper part of CollectionMatrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - CollectionMatrix assembly : Tangent Op PARTS'
      CALL FLUSH(6)
      TIC = CPUTime()
     
      DO i=StiffMatrix % NumberOfRows,1,-1      
        DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1
          CALL AddToMatrixElement( CollectionMatrix, i, StiffMatrix % Cols(j), StiffMatrix % Values(j) )
        END DO
      END DO

      TIC = CPUTime() - TIC
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - CollectionMatrix assembly : StiffMatrix   TIME =',TIC
      CALL FLUSH(6)      
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - List_toCRSMatrix'
      CALL FLUSH(6)
      CALL List_toCRSMatrix(CollectionMatrix)


! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Solve the Collection-system 
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! TODO  : SOLVE WITH MUMPS
!       CALL SolveLinearSystem( CollectionMatrix, CollectionVector, &
!                               CollectionSolution, Norm, DOFs, Solver, StiffMatrix )
      WRITE(6,*) '- BIFStaDimKer1ComputePsiW - HMumps_SolveLinearSystem'
      CALL FLUSH(6)
      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       First Psi
!
! | LtcT V |  PSI   0
! |        |      = 
! | VT   0 |  k     1
!
      write(6,*) "BIFStaDimKer1ComputePsiW - - - - - - - - - - - - - - - - - - - - - - -"
      write(6,*) "BIFStaDimKer1ComputePsiW - >>       Compute Left Mode PSI          <<"
      write(6,*) "BIFStaDimKer1ComputePsiW - >> MUMPS V5 SOLVE TRANSPOSE ICNTL(9)!=1 <<"
      CALL FLUSH(6)

      CollectionVector = 0.0_dp
      CollectionVector(k+1) = 1.0_dp
      
!       CollectionMatrix % mumpsIDL % ICNTL(9) = 2
      CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .TRUE. ) 
      CALL HMumps_SolveLinearSystemTransp( CollectionMatrix, CollectionSolution, CollectionVector, &
                                     Solver , MUMPSFICH )
      CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .FALSE. )                                       
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Separate the solution from CollectionSolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
      abePsiL = 0.0_dp
      abePsiL(1:k) = CollectionSolution(1:k)
      restv=CollectionSolution(k+1)
      write(6,*) "BIFStaDimKer1ComputePsiW - PSI RestVect=",restv
      CALL FLUSH(6)

                                     
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       Second W ( Such that we solve Ax=b for the branch using the same matrix)
!
! | Ltc V |  W    F
! |       |     = 
! | VT  0 |  k    0
!
      write(6,*) "BIFStaDimKer1ComputePsiW - - - - - - - - - - - - - - - - - - - - - - -"
      write(6,*) "BIFStaDimKer1ComputePsiW - >> Compute vector W <<"
      CALL FLUSH(6)
      CollectionVector(1:k) =  ForceVector(1:k)
      CollectionVector(k+1) = 0.0_dp

!       CollectionMatrix % mumpsIDL % ICNTL(9) = 1
      CALL HMumps_SolveLinearSystem( CollectionMatrix, CollectionSolution, CollectionVector, &
                                     Solver , MUMPSFICH )   
      
                  
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Separate the solution from CollectionSolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - -
      abeW = 0.0_dp
      i = 1
      j = StiffMatrix % NumberOfRows
      abeW(i:j) = CollectionSolution(i:j)
      restv=CollectionSolution(j+1)
      write(6,*) "BIFStaDimKer1ComputePsiW - W RestVect=",restv
      CALL FLUSH(6)
!       WRITE(*,*) 'ASSOCIATED(CollectionMatrix) in= ',ASSOCIATED(CollectionMatrix)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      DEALLOCATE(CollectionSolution)
      CALL Info( 'BIFStaDimKer1ComputePsiW', 'Vectors computed', Level=5 )
      write(6,*) "BIFStaDimKer1ComputePsiW - - - - - - - - - - - - - - - - - - - - - - -"
      CALL FLUSH(6)

      
  END SUBROUTINE BIFStaDimKer1ComputePsiW
    
!---------------------------------------------------------------------------------------------








 
 ! ! ! TODO
!------------------------------------------------------------------------------
! VERIF POINT LIMITE
!  D L / Da en acritique !=0 => Bifurcation
! SERIES
! PADE
! !   FUNCTION IS_IT_A_LIMIT_POINT( LAMBDA, Acrit ,  NORDRE ) RESULT ( LIMITPOINT )
! !       IMPLICIT NONE
! !       INTEGER       :: NORDRE      
! !       REAL(KIND=dp) :: LAMBDA(NORDRE), Acrit
! !       LOGICAL       :: LIMITPOINT
! ! !      
! !       REAL(KIND=dp) :: lambdaDER1(NORDRE-1)
! !       INTEGER       :: i
! ! !      
! !       LIMITPOINT = .FALSE.
! !       
! ! ! FIRST DERIVATIVE      
! !       DO i=1, size(NoPressure)
! !         UTemp(i) = Solution(i)*NoPressure(i)
! !       END DO
! !       
! ! ! EVALUATION AT CRITICAL PARAMTER        
! !       pamax = 1._dp
! !       lambdatemp = 0.0_dp
! !       DO i = 1 , NORDRE-1
! !        pamax = pamax * Acrit
! !        lambdatemp = lambdatemp + pamax * lambdaDER1( i )
! !       END DO
! !         
! !  END FUNCTION IS_IT_A_LIMIT_POINT




! ! ! TODO
!------------------------------------------------------------------------------

! !    SUBROUTINE BIFs_ComputeTangents(UBif,LBif,BifMode,                   &
! !                                    UTa1,LTa1,UTb1,LTb1)
! ! 
! !  
! !    END SUBROUTINE 

 
 

!---------------------------------------------------------------------------------------
! ______                      _               
! | ___ \                    | |              
! | |_/ /_ __ __ _ _ __   ___| |__   ___  ___ 
! | ___ \ '__/ _` | '_ \ / __| '_ \ / _ \/ __|
! | |_/ / | | (_| | | | | (__| | | |  __/\__ \
! \____/|_|  \__,_|_| |_|\___|_| |_|\___||___/
!
! TODO : Calcul des Aij une seul fois car Q(A,B) coute du temps avec 10Mddl !
  
   SUBROUTINE BIFs_Compute_serie_one_branch( SerUbranch, SerLbranch , VectUT, LT,       &
                                            BifMode,abeW,abePsiL, NDL,NSDOFs,OMan,      &
                                            FQMan, FQManTemp, USAV, GradSAV,            &
                                            Uelex, Ueley, Uelez, Velex, Veley, Velez,   &
                                            DensityTMP, Material,                       &
                                            Solver,StiffMatrix, FlowSolution_Init,      &
                                            MUMPSFICH,Uchaptmp,AUGStiffMatrix,          &
                                            Vectmpa,Vectmpb,FlowPerm)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
! - - - - - - - - - - - - - - -     
    USE HomeMadeSolvers
    USE DiscreteOperators    
! - - - - - - - - - - - - - - - 

!----------------------
  IMPLICIT NONE
! EXTER
      TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
      TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 

      REAL(KIND=dp) :: LT
      INTEGER       :: NDL,NSDOFs,MUMPSFICH
      REAL(KIND=dp) :: VectUT(:) , BifMode(:) , abeW(:), abePsiL(:)
      REAL(KIND=dp) :: FQMan(:),FQManTemp(:), USAV(:,:,:,:), GradSAV(:,:,:,:,:)
      REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:),Velex(:), Veley(:), Velez(:)
      REAL(KIND=dp) :: DensityTMP(:),FlowSolution_Init(:),Uchaptmp(:)
      ! OUT EXTER
      REAL(KIND=dp) :: SerUbranch(:,:),SerLbranch(:)
      TYPE(ValueList_t),POINTER :: Material
      REAL(KIND=dp)          :: Vectmpa(:), Vectmpb(:)
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
! LOCALE
      TYPE(Matrix_t),POINTER :: AUGStiffMatrix
      REAL(KIND=dp), POINTER,CONTIGUOUS :: AUGRHS(:)
      INTEGER, POINTER       :: FlowPerm(:)   
      INTEGER       :: p,OMan,j      
      REAL(KIND=dp) :: A11,A12,A21,A22,B1,B2,XDENTMP,XNUMTP,Eta(OMan)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -

      WRITE(*,*) "BrSw - BIFs_Compute_serie_one_branch"
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
!       WRITE(*,*) "FlowPerm       => Solver % Variable % Perm"
!       FlowPerm       => Solver % Variable % Perm
!        WRITE(*,*) "Associated(StiffMatrix % CollectionMatrix)=",Associated(StiffMatrix % CollectionMatrix)
!       WRITE(*,*) "AUGStiffMatrix => StiffMatrix % CollectionMatrix"   
!       AUGStiffMatrix => StiffMatrix % CollectionMatrix
!       WRITE(*,*) "AUGRHS         => AUGStiffMatrix % RHS"      
      AUGRHS         => AUGStiffMatrix % RHS
!       WRITE(*,*) "AUGRHS = 0"      
      AUGRHS   = 0.0_dp    
      USAV     = 0.0_dp
      GradSAV  = 0.0_dp        
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
! Ut1 , Lt1 connus, on cherche la branche comme série ANM
      SerUbranch( : , 1 ) = VectUT
      SerLbranch( 1 )     = LT
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
! Ubranch_p = L_p W + E_p BifMode + Uchap_p   
      DO p=2, OMan   
      WRITE(*,*) "BrSw Branch serie Order=",p
      
!    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -
! 1 ) Uchap p  : Ltc  Uchaptmp p = - FQ  (sys aug like W)
! UCHAP
        Uchaptmp = 0.0_dp
        FQMan    = 0.0_dp
        CALL ANMPFrhsFQ(FQMan, SerUbranch, p, NSDOFs, FlowPerm , &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez , &
                         DensityTMP,Material,FlowSolution_Init)
!          WRITE(*,*) "BrSw - NORME FQMan = ",DSQRT(DOT_PRODUCT(FQMan,FQMan))
         ! Resolution systeme augmenté identique à Ltc W = F Sauf Second membre

         AUGRHS(1:NDL) = FQMan(1:NDL)   ! avec FQman = -Sum Qcrois
         ! SOLVE
!          WRITE(*,*) "BrSw - AUGRHS = - sum Q(Ubran)"      
         
         CALL HMumps_SolveLinearSystem( AUGStiffMatrix , Uchaptmp, AUGRHS, Solver , &
                                        MUMPSFICH)
!          WRITE(*,*) "BrSw - HMumps_SolveLinearSystem DONE"      
!          CALL CalCondLimCoeff( Uchaptmp,  0.0_dp, 0, FlowSolution_Init  )      
!Modif 20150415
         CALL CalCondLimCoeff( Uchaptmp,  0.0_dp, 1, FlowSolution_Init  )      ! 0.0_dp, 1 => put 0 : CALCOND and CALIMP
! 0 : put a given vector for both Constraints and Load
! 1 : multiply Both contraints and load by a scalar
! 2 : Only Constraints (CALIMP EVE)                                              
! UCHAP Branch OK       
!    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -
! 2 ) Lp et Ep : path parameter projection a = <u-u0,u1>
!
!    A11 * Lp + A12 * Ep = B1
!    A21 * Lp + A22 * Ep = B2
!
! => Lp = ( B1 - B2 * A12/A22 ) / ( A11 - A22 * A12/A21 )
! => Ep = ( B2 -  A21 * LP) / A22
!
! a   A11 = < abePsiL, Q(Ut1,W) + Q(W,Ut1)>
!          WRITE(*,*) "BrSw - A11 = < abePsiL, Q(Ut1,W) + Q(W,Ut1)>"
          CALL ANMQAB( Vectmpa , abeW, VectUT, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)
          CALL ANMQAB( Vectmpb , VectUT,abeW, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)   
          A11 = DOT_PRODUCT(abePsiL,Vectmpa)
          A11 = A11 + DOT_PRODUCT(abePsiL,Vectmpb)                         
          
! b   A12 = < abePsiL, Q(Ut1 ,BM) + Q(BM,Ut1)>
!           WRITE(*,*) "BrSw - A12 = < abePsiL, Q(Ut1 ,BM) + Q(BM,Ut1)>"
          CALL ANMQAB( Vectmpa , VectUT,BifMode, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)
          CALL ANMQAB( Vectmpb , BifMode,VectUT, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)   
          A12 = DOT_PRODUCT(abePsiL,Vectmpa)
          A12 = A12 + DOT_PRODUCT(abePsiL,Vectmpb)  
          
! d   A21 = <  W , Ut1  >
!           WRITE(*,*) "BrSw - A21 = <  W , Ut1  >"
!           A21 = DOT_PRODUCT(abeW,VectUT)
!  Pilota U lambda : ajout + L1
          A21 = DOT_PRODUCT(abeW,VectUT) + SerLbranch( 1 )

! e   A22 = < BM , Ut1 >
!           WRITE(*,*) "BrSw -  A22 = < BM , Ut1 >"
          A22 = DOT_PRODUCT(BifMode,VectUT)

! c     B1 = - < abePsiL, Q( Ut1 , Uchap p ) + Q( Uchap p , Ut1 ) > 
!           WRITE(*,*) "BrSw - -g1 = - < abePsiL, Q( Ut1 , Uchap p ) - Q( Uchap p , Ut1 ) > "
          CALL ANMQAB( Vectmpa , VectUT, Uchaptmp, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)
          CALL ANMQAB( Vectmpb , VectUT, Uchaptmp, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)                         
         B1 = DOT_PRODUCT(abePsiL,Vectmpa)
         B1 = B1 + DOT_PRODUCT(abePsiL,Vectmpb)
!        -  Sum_{j=2}^{p-1}  < abePsiL, Q(Uj , U_{P+1-j}) ! RHS FQ MODIFIED
!           WRITE(*,*) "BrSw  - < abePsiL, ANMBranchrhsFQ2> "
         IF (p.GT.2) THEN
           CALL ANMBranchrhsFQ( Vectmpa, SerUbranch, p, NSDOFs, FlowPerm,    &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez, &
                         DensityTMP, Material, FlowSolution_Init)              
           B1 = B1 +  DOT_PRODUCT(abePsiL,Vectmpa)
         ENDIF
         B1 = -B1
        
!        B2 = - <Uchap p , Ut1>
         B2 = DOT_PRODUCT(Uchaptmp,VectUT)
         B2= -B2
!    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -

         WRITE(*,*) "_______________________________________"
         WRITE(*,*) "A11=",A11," A12=",A12
         WRITE(*,*) "A21=",A21," A22=",A22
         WRITE(*,*) "B1=",B1
         WRITE(*,*) "B2=",B2
         WRITE(*,*) "_______________________________________"  
!          CALL SolveLinSys2x2(A,x,b)
! => Lp = ( B1 - B2 * A12/A22 ) / ( A11 - A21 * A12/A22 )
         XNUMTP  = A12 / A22
         XNUMTP  = B2 * XNUMTP
         XNUMTP  = B1 - XNUMTP
         XDENTMP = A12 / A22
         XDENTMP = A21 * XDENTMP
         XDENTMP = A11 - XDENTMP
         SerLbranch(p) = XNUMTP /  XDENTMP
         WRITE(*,*) 'L(',p,',) =',SerLbranch(p)  
         
! => Ep = ( B2 -  A21 * LP) / A22         
         XNUMTP =  A21 * SerLbranch(p)
         XNUMTP =  B2 - XNUMTP
         Eta(p) =  XNUMTP / A22
         WRITE(*,*) 'Eta(',p,')=',Eta(p)
         WRITE(*,*) "_______________________________________"         
         
!====> U(p) = L(p) W + E(p) BM + Uchap(p)
         SerUbranch(:,p) = SerLbranch(p)*abeW + Eta(p)*BifMode + Uchaptmp
!          CALL CalCondLimCoeff( SerUbranch(:,p),  0.0_dp, 0, FlowSolution_Init  )      
         CALL CalCondLimCoeff( SerUbranch(:,p),  0.0_dp, 1, FlowSolution_Init  ) ! 0.0_dp, 1 => put 0 : CALCOND and CALIMP     
         
       END DO ! fin boucle sur les ordres
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                            

!    { UBranch,LBranch } series

   END SUBROUTINE 
 
                           
                           
                           
                           
                           
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
! WITH PITCHFORK : W is sym tangent and Phi(bifmode) is asym tangent
! Resolution is simpler(?)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -  
   SUBROUTINE BIFs_Compute_serie_one_branch_SYM( SerUbranch, SerLbranch , VectUT, LT,       &
                                            BifMode,abeW,abePsiL, NDL,NSDOFs,OMan,      &
                                            FQMan, FQManTemp, USAV, GradSAV,            &
                                            Uelex, Ueley, Uelez, Velex, Veley, Velez,   &
                                            DensityTMP, Material,                       &
                                            Solver,StiffMatrix, FlowSolution_Init,      &
                                            MUMPSFICH,Uchaptmp,AUGStiffMatrix,          &
                                            Vectmpa,Vectmpb,FlowPerm,SYMBRANCH,abeCoB)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
! - - - - - - - - - - - - - - -     
    USE HomeMadeSolvers
    USE DiscreteOperators    
! - - - - - - - - - - - - - - - 
!----------------------
  IMPLICIT NONE
! EXTER
      TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
      TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 

      REAL(KIND=dp) :: LT
      INTEGER       :: NDL,NSDOFs,MUMPSFICH
      REAL(KIND=dp) :: VectUT(:) , BifMode(:) , abeW(:), abePsiL(:)
      REAL(KIND=dp) :: FQMan(:),FQManTemp(:), USAV(:,:,:,:), GradSAV(:,:,:,:,:)
      REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:),Velex(:), Veley(:), Velez(:)
      REAL(KIND=dp) :: DensityTMP(:),FlowSolution_Init(:),Uchaptmp(:)
      ! OUT EXTER
      REAL(KIND=dp) :: SerUbranch(:,:),SerLbranch(:)
      TYPE(ValueList_t),POINTER :: Material
      REAL(KIND=dp)          :: Vectmpa(:), Vectmpb(:)
      LOGICAL       :: SYMBRANCH
      REAL(KIND=dp) :: abeCoB
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
! LOCAL
      TYPE(Matrix_t),POINTER :: AUGStiffMatrix
      REAL(KIND=dp), POINTER,CONTIGUOUS :: AUGRHS(:)
      INTEGER, POINTER       :: FlowPerm(:)   
      INTEGER       :: p,OMan,j      
      REAL(KIND=dp) :: A11,A12,A21,A22,B1,B2,XDENTMP,XNUMTP,Eta(OMan), NPSIL
      REAL(KIND=dp) :: NTMP
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -

      WRITE(*,*) "BrSw - BIFs_Compute_serie_one_branch_SYM"
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
!       WRITE(*,*) "FlowPerm       => Solver % Variable % Perm"
!       FlowPerm       => Solver % Variable % Perm
!        WRITE(*,*) "Associated(StiffMatrix % CollectionMatrix)=",Associated(StiffMatrix % CollectionMatrix)
!       WRITE(*,*) "AUGStiffMatrix => StiffMatrix % CollectionMatrix"   
!       AUGStiffMatrix => StiffMatrix % CollectionMatrix
!       WRITE(*,*) "AUGRHS         => AUGStiffMatrix % RHS"      
      AUGRHS         => AUGStiffMatrix % RHS
!       WRITE(*,*) "AUGRHS = 0"      
      AUGRHS   = 0.0_dp    
      USAV     = 0.0_dp
      GradSAV  = 0.0_dp        
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -
! Ut1 , Lt1 connus, on cherche la branche comme série ANM
      SerUbranch( : , 1 ) = VectUT
      SerLbranch( 1 )     = LT
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -


! Ubranch_p = L_p W + E_p BifMode + Uchap_p   
      DO p=2, OMan   
      WRITE(*,*) "BrSw Branch serie Order=",p
      
!    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -
! 1 ) Uchap p  : Ltc  Uchaptmp p = - FQ  (sys aug like W)
! UCHAP
        Uchaptmp = 0.0_dp
        FQMan    = 0.0_dp
        AUGRHS   = 0.0_dp    
         
        CALL ANMPFrhsFQ(FQMan, SerUbranch, p, NSDOFs, FlowPerm , &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez , &
                         DensityTMP,Material,FlowSolution_Init)
!          WRITE(*,*) "BrSw - NORME FQMan = ",DSQRT(DOT_PRODUCT(FQMan,FQMan))
         ! Resolution systeme augmenté identique à Ltc W = F Sauf Second membre

         AUGRHS(1:NDL) = FQMan(1:NDL)   ! avec FQman = -Sum Qcrois
         ! SOLVE
!          WRITE(*,*) "BrSw - AUGRHS = - sum Q(Ubran)"      
         
         CALL HMumps_SolveLinearSystem( AUGStiffMatrix , Uchaptmp, AUGRHS, Solver , &
                                        MUMPSFICH)
!          WRITE(*,*) "BrSw - HMumps_SolveLinearSystem DONE"      
!          CALL CalCondLimCoeff( Uchaptmp,  0.0_dp, 0, FlowSolution_Init  )      
!Modif 20150415
         CALL CalCondLimCoeff( Uchaptmp,  0.0_dp, 1, FlowSolution_Init  )      ! 0.0_dp, 1 => put 0 : CALCOND and CALIMP
! 0 : put a given vector for both Constraints and Load
! 1 : multiply Both contraints and load by a scalar
! 2 : Only Constraints (CALIMP EVE)                                              
! UCHAP Branch OK       
!    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -
! 2 ) Lp et Ep : path parameter projection a = <u-u0,u1>, ou complet avec Lambda?
!
!    A11 * Lp + A12 * Ep = B1
!    A21 * Lp + A22 * Ep = B2
!
! => Lp = ( B1 - B2 * A12/A22 ) / ( A11 - A22 * A12/A21 )
! => Ep = ( B2 -  A21 * LP) / A22
!

! c     B1 = - < abePsiL, Q( Ut1 , Uchap p ) + Q( Uchap p , Ut1 ) > 
!           WRITE(*,*) "BrSw - -g1 = - < abePsiL, Q( Ut1 , Uchap p ) - Q( Uchap p , Ut1 ) > "
          CALL ANMQAB( Vectmpa , VectUT, Uchaptmp, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)
          CALL ANMQAB( Vectmpb , Uchaptmp, VectUT, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)                         
         B1 = DOT_PRODUCT(abePsiL,Vectmpa)
         B1 = B1 + DOT_PRODUCT(abePsiL,Vectmpb)
!        -  Sum_{j=2}^{p-1}  < abePsiL, Q(Uj , U_{P+1-j}) ! RHS FQ MODIFIED
!           WRITE(*,*) "BrSw  - < abePsiL, ANMBranchrhsFQ2> "
         IF (p.GE.3) THEN
           CALL ANMBranchrhsFQ( Vectmpa, SerUbranch, p, NSDOFs, FlowPerm,    &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez, &
                         DensityTMP, Material, FlowSolution_Init)              
           B1 = B1 +  DOT_PRODUCT(abePsiL,Vectmpa)
         ENDIF
         B1 = -B1
        
!        B2 = - <Uchap p , Ut1>
         B2 = DOT_PRODUCT(Uchaptmp,VectUT)
         B2= -B2
!    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -    -
! On peut normer LeftMode pour assurer les projections? en test si ok faire aussi avant pour coeff ABE
!          NPSIL = DSQRT(DOT_PRODUCT(abePsiL,abePsiL))
!          abeCoB=abeCoB/NPSIL         
         abeCoB=abeCoB       
         IF (SYMBRANCH) THEN
           WRITE(*,*) "________________SYM BRANCH____(L1W,L1)______________________"
           SerLbranch(p) = B2 / (DOT_PRODUCT(VectUT,VectUT) + SerLbranch(1))
           WRITE(*,*) 'L(',p,',) =',SerLbranch(p)  
           Eta(p) =  B1/abeCoB
           WRITE(*,*) 'Eta(',p,')=',Eta(p)
         ELSE
           WRITE(*,*) "_______________ASYM BRANCH____(BifMode,0)___________________"                           
           SerLbranch(p) = B1/abeCoB
           WRITE(*,*) 'L(',p,',) =',SerLbranch(p)  

           Eta(p) =  B2
           WRITE(*,*) 'Eta(',p,')=',Eta(p)
         ENDIF
         WRITE(*,*) "_______________________________________"         
         
!====> U(p) = L(p) W + E(p) BM + Uchap(p)
         SerUbranch(:,p) = SerLbranch(p)*abeW + Eta(p)*BifMode + Uchaptmp
!          CALL CalCondLimCoeff( SerUbranch(:,p),  0.0_dp, 0, FlowSolution_Init  )      
         CALL CalCondLimCoeff( SerUbranch(:,p),  0.0_dp, 1, FlowSolution_Init  ) ! 0.0_dp, 1 => put 0 : CALCOND and CALIMP     
         NTMP = DOT_PRODUCT(SerUbranch(:,p),SerUbranch( : , 1 )) + SerLbranch(p) * SerLbranch(1)
         WRITE(*,*) 'VERIF COND PROJECTION <up,u1>+lpl1=',NTMP
          
       END DO ! fin boucle sur les ordres
! - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                            

!    { UBranch,LBranch } series

   END SUBROUTINE BIFs_Compute_serie_one_branch_SYM
 
 
 
 
 
 
 
!------------------------------------------------------------------------------
! > Multiply the boundary conditions values by Lamda0 : Lambda x Ud
! See ListGetReal to understand :
!     LIST_TYPE_VARIABLE_SCALAR      :  F(i) = ptr % Coeff * ExecRealFunction (MATC?) 
!      F(1) = ptr % Coeff * F(1)
!     LIST_TYPE_CONSTANT_SCALAR_PROC :  F(i) = Ptr % Coeff * ExecConstRealFunction(MATC?  x,y,z
! ------------------------------------------------------------------------------
     SUBROUTINE BoundaryConditionMultCOEFF(Model, Lambda0, Reset)
     
        TYPE(Model_t)                :: Model          !< The current model structure
        REAL(KIND=dp)                :: Lambda0        !< Coeff as Lambda X Ud
! INTERNAL       
        TYPE(Solver_t), POINTER      :: Solver
        INTEGER                      :: DOF, i
        LOGICAL                      :: gotIt, Reset
        CHARACTER(LEN=MAX_NAME_LEN)  :: Name           !< Name of the dof to be set
        TYPE(ValueList_t), POINTER   :: ValueList, ptr
        
        Solver => Model % Solver
! RESET TODO : GET MAX ON CONTOUR AND NORMALIZE => THEN YOU HAVE LAMBDA0        
! MAX FOR A COMPONENT
        IF (Reset) Lambda0=1._dp 
!       DO DOF        
        DO DOF=1,Solver % Variable % DOFs
!         GET NAME       
          Name = Solver % Variable % Name
          IF (Solver % Variable % DOFs>1) Name=ComponentName(Name,DOF)
!         DO BC
          DO i=1,Model % NumberOfBCs
!           write(*,*) "DOF=",DOF," Name=",TRIM(Name),"BCs(",i,")"        
            ValueList => Model % BCs(i) % Values
            ptr => ListFind(ValueList,Name,gotIt)
            IF ( .NOT.ASSOCIATED(ptr) ) cycle
!           write(*,*) "ptr % Coeff=",ptr % Coeff
            ptr % Coeff = Lambda0 * ptr % Coeff
!           write(*,*) "ptr % Coeff=",ptr % Coeff
          END DO
        END DO
        
     END SUBROUTINE BoundaryConditionMultCOEFF
! ------------------------------------------------------------------------------

 
! ------------------------------------------------------------------------------
!      SUBROUTINE SENSCO(V,U,IP,APOSITIF,RAYON,NPOI,DA,NDL)
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!       LOGICAL APOSITIF 
!       DIMENSION V(NDL,1),U(1)
!       IF(IP.EQ.1) THEN
!        APOSITIF=.TRUE.
!        DA=RAYON/DBLE(NPOI)          
!       ELSE      
!         PROJ=PROSCA(V(1,1),U,NDL)
!         IF(PROJ.GT.0) THEN
!          APOSITIF=.TRUE.
!          WRITE(6,*) 'Dans Sensco positif'
!         DA=RAYON/DBLE(NPOI)       
!         ELSE
!          APOSITIF=.FALSE.
!          WRITE(6,*) 'Dans Sensco negatif'                
!          DA=-RAYON/DBLE(NPOI)
!         ENDIF
!       ENDIF
!       END
! ------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! UTILS
!       NORM = DSQRT( DOT_PRODUCT( VECT, VECT ) )

!>  On retire le terme de pression dans le vecteur solution pour le calcul des 
!       normes, lambda et P.S
! PRESSURE FREE of the vector (U,V,P,U,V,P,...) => (U,V,0,U,V,0....) to compute range of VALIDITY
!------------------------------------------------------------------------------
  SUBROUTINE PressureFree( Solution, NoPressure, UTemp )
      IMPLICIT NONE
      REAL(KIND=dp), DIMENSION(:) :: Solution, UTemp
      REAL(KIND=dp), DIMENSION(:) :: NoPressure
      INTEGER :: n,i
       
      DO i=1, size(NoPressure)
        UTemp(i) = Solution(i)*NoPressure(i)
      END DO
     
 END SUBROUTINE PressureFree

!------------------------------------------------------------------------------
 character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
!------------------------------------------------------------------------------
 
 
 
 
 
 
 
END MODULE ANMToolBox

!> \} IMNS
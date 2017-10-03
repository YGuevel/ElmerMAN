!> \ingroup IMNS
!> \}

!> \ingroup IMNS

! # YG@2016
! # yann.guevel@univ-ubs.fr
! # Part of ANM continuation/detection/branch-switching in ELMER FEM
! # Jean-Marc Cadou / Gregory girault / Yann Guevel
! # Institut de Recherche Dupuy de Lome - Universite Bretagne Sud


!------------------------------------------------------------------------------
!>  Utility routines for FEM computations
!------------------------------------------------------------------------------
MODULE FEMtoolBox


! + DIRICHLET MAISON                 : CalCondLimCoeff  (contains getBoundaryIndexesMAN)


!------------------------------------------------------------------------------

  USE SolverUtils
  USE DefUtils

!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
    

CONTAINS


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!------------------------------------------------------------------------------
! CalCondLimCoeff Reads boundary condition coutour values
!   - Constraints : like a 0 on velocity Dof
!   - Load : given velocity profile on contour
! Argument COND
! 0 : put a given vector for both Constraints and Load
! 1 : multiply Both contraints and load by a scalar
! 2 : Only Constraints (CALIMP EVE)
!
!> ANM : For ANM orders > 1 set to Zero all the RHS DoF in Boundary 
!>       RHS Velocity components <- 0 
!>       EVE CALIMP and CALCONDILM
!> Sets the Dirichlet conditions related to the variables of the active solver.
!>
! ! ! e1 bndry p1 p2 type n1 ... nn
! ! ! e2 bndry p1 p2 type n1 ... nn
! ! ! ...
! ! ! en bndry p1 p2 type n1 ... nn
! ! ! ++ The first integer is again the identification number of the boundary element. 
! ! ! ++ Next the identification number of the part of the boundary where this element is located is given.
! ! ! ++ Whether the boundary element can be represented as the side of a parent element defined in the file 
! ! ! mesh.elements is indicated using the two parent element numbers p1 and p2. 
! ! ! -> If the boundary element  is located on an outer boundary of the body, it has only one parent 
! ! ! element and either of these  two integers may be set to be zero. 
! ! ! It is also possible that both parent element numbers are zeros. 
! ! ! ++ Finally the element type code and element nodes are listed.
!------------------------------------------------------------------------------------------
  SUBROUTINE CalCondLimCoeff( Vector, Coeff, Cond, VectBCInit, USolver,Ux,UOffset,OffDiagonalMatrix,PiolaCurlTransform )
!------------------------------------------------------------------------------------------
            USE Adaptive
            USE SolverUtils
            USE DefUtils

     IMPLICIT NONE
     
     INTEGER, OPTIONAL :: UOffset
     REAL(KIND=dp) :: Coeff ! reel dp 0 ou "lambda"
     INTEGER :: Cond ! 0 si on replaque les CIs et BCs , 1 sinon
     REAL(KIND=dp), OPTIONAL, DIMENSION(:) :: VectBCInit
     LOGICAL, OPTIONAL :: OffDiagonalMatrix
     TYPE(Variable_t), OPTIONAL, TARGET :: Ux
     TYPE(Solver_t), OPTIONAL, TARGET :: USolver
     LOGICAL, OPTIONAL :: PiolaCurlTransform  ! An additional argument for indicating that
                                              ! the solution is expanded in terms of H(curl)-
                                              ! conforming basis functions defined via the 
                                              ! Piola transform.  
!--------------------------------------------------------------------------------------------     
     TYPE(Matrix_t), POINTER   :: A
     TYPE(Variable_t), POINTER :: x
     TYPE(Solver_t), POINTER :: Solver
     REAL(KIND=dp), POINTER    :: b(:)
     
     

     REAL(KIND=dp) :: xx, temp
     REAL(KIND=dp), POINTER :: DiagScaling(:)
     REAL(KIND=dp), ALLOCATABLE :: Work(:), STIFF(:,:)

     INTEGER, ALLOCATABLE :: lInd(:), gInd(:)
     INTEGER :: i,j, k, kk, l, m, n,nd, nb, mb, nn, ni, nj, &
          DOF, local, numEdgeDofs,istat, n_start, Offset, iter, &     
          n1, n2, n3, n4

     LOGICAL :: Flag,Found, ConstantValue, ScaleSystem
     TYPE(ValueList_t), POINTER :: BC, ptr, Params
     TYPE(Element_t), POINTER :: Element, Parent, Edge, Face, SaveElement

     CHARACTER(LEN=MAX_NAME_LEN) :: name
     LOGICAL :: BUpd

     INTEGER::iii=0
! ******************************************************
! ******************************************************
! ******************************************************
     
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: xip,yip,zip,s,DetJ,Load
        REAL(KIND=dp) :: Basis(50)
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL :: stat
    INTEGER :: p,q,t,ndd,bcid
    
! ******************************************************
! ******************************************************
! ******************************************************


! ******************************************************
!    Right Hand Side for ANM high order 
     REAL(KIND=dp) :: Vector(:)
! ******************************************************
     
     SAVE gInd, lInd, STIFF, Work

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF
     
     Params => GetSolverParams(Solver)

!---------------------------------------------------     
     A => Solver % Matrix
     b => A % RHS
     
     
     IF ( PRESENT(Ux) ) THEN
       x => Ux
     ELSE
       x => Solver % Variable
     END IF
     IF(.NOT.ALLOCATED(A % ConstrainedDOF)) &
       ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
     A % ConstrainedDOF = .FALSE.
!---------------------------------------------------     
     Offset = 0
     IF(PRESENT(UOffset)) Offset=UOffset

     n = Solver % Mesh % MaxElementDOFs
     IF ( .NOT. ALLOCATED( gInd ) ) THEN
       ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
       IF ( istat /= 0 ) &
          CALL Fatal('ANM::CalCondLimCoeff','Memory allocation failed.' )
     ELSE IF ( SIZE(gInd) < n ) THEN
       DEALLOCATE( gInd, lInd, STIFF, Work )
       ALLOCATE( gInd(n), lInd(n), STIFF(n,n), Work(n), stat=istat )
       IF ( istat /= 0 ) &
          CALL Fatal('ANM::CalCondLimCoeff','Memory allocation failed.' )
     END IF
     IF ( x % DOFs > 1 ) THEN

        CALL SetDirichletBoundaries( CurrentModel,A, b, GetVarName(x),-1,x % DOFs,x % Perm )
     END IF
     CALL Info('ANM::CalCondLimCoeff', &
            'Setting Dirichlet boundary conditions for the ANM RHS', Level=5)

!      WRITE(*,*) "----------------------------------------------------------"
!      WRITE(*,*) "RECAP TEST 3D"
!      WRITE(*,*) "x % DOFs = ",x % DOFs
!      WRITE(*,*) "Solver % Mesh % NumberOfBoundaryElements = ",Solver % Mesh % NumberOfBoundaryElements
!      WRITE(*,*) "----------------------------------------------------------"            
     ! Set Dirichlet dofs for edges and faces
     ConstantValue = .FALSE.
     DO DOF=1,x % DOFs
        name = x % name

        IF ( x % DOFs > 1 ) name = ComponentName(name,DOF)
        ! clear bc face & edge dofs
        SaveElement => CurrentModel % CurrentElement
        DO i=1,Solver % Mesh % NumberOfBoundaryElements
           Element => GetBoundaryElement(i)

           IF ( .NOT. ActiveBoundaryElement() ) CYCLE

           ! Get parent element:
           ! -------------------
           Parent => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( Parent ) ) &
             Parent => Element % BoundaryInfo % Right

           IF ( .NOT. ASSOCIATED(Parent) )   CYCLE

           BC => GetBC()
           IF ( .NOT. ASSOCIATED(BC) ) CYCLE
           
           ptr => ListFind(BC, Name,Found )

           
           IF ( .NOT. ASSOCIATED(ptr) ) CYCLE
           ConstantValue =  ptr % PROCEDURE == 0 .AND. &
             ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR
           
           ! Get indexes for boundary and values for dofs associated to them
           n = GetElementNOFNodes()
           IF ( isActivePElement(Parent)) THEN
             CALL getBoundaryIndexes( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )
           ELSE
             CYCLE 
           END IF
           
! ! !            ! Contribute this boundary to global system
! ! !            ! (i.e solve global boundary problem)
           DO k=n+1,numEdgeDofs
!              write(*,*) 'gInd',gind(k)
             nb = x % Perm( gInd(k) )
             IF ( nb <= 0 ) CYCLE
             nb = Offset + x % DOFs * (nb-1) + DOF
           END DO
        END DO
        CurrentModel % CurrentElement => SaveElement
     END DO
 
     ! Set Dirichlet dofs for edges and faces
    
     DO DOF=1,x % DOFs
        name = x % name

        IF (x % DOFs>1) name=ComponentName(name,DOF)
!         write(*,*) 'DOFs=',DOF,'Name=', name
        CALL SetNodalLoads( CurrentModel,A, b, &
            Name,DOF,x % DOFs,x % Perm ) ! , Offset ) not yet ?

        CALL SetDirichletBoundaries( CurrentModel, A, b, &
             Name, DOF, x % DOFs, x % Perm, Offset, OffDiagonalMatrix )

!       Dirichlet BCs for face & edge DOFs:
!       -----------------------------------
        SaveElement => CurrentModel % CurrentElement
        
        DO i=1,Solver % Mesh % NumberOfBoundaryElements
!            write(*,*) "i=",i
           Element => GetBoundaryElement(i)
           IF ( .NOT. ActiveBoundaryElement() ) CYCLE
!            write(*,*) "GetBC()"
           BC => GetBC()
! ******************************************************
! ******************************************************
           
           IF ( .NOT. ASSOCIATED(BC) ) CYCLE
           IF ( .NOT. ListCheckPresent(BC, Name) .AND. &
                .NOT. ListCheckPresent(BC, TRIM(Name)//' {e}') .AND. &
                .NOT. ListCheckPresent(BC, TRIM(Name)//' {f}') ) CYCLE
           ! Get parent element:
           ! -------------------
           
           Parent => Element % BoundaryInfo % Left
           IF ( .NOT. ASSOCIATED( Parent ) ) THEN
               Parent => Element % BoundaryInfo % Right
           END IF

           IF ( .NOT. ASSOCIATED( Parent ) )   CYCLE
!            write(*,*) " ASSOCIATED( Parent ) OK"
           
!---------------------------------------------------     
!!! TEST COMMENTE avant pour 2D valider
!!!! TEST ISACTIVEPELEMENT...?
!         write(*,*) "isPelement(Parent)=",isPelement(Parent)
!         write(*,*) Parent % Type % ElementCode / 100
!         write(*,*) CurrentModel % Solver % Def_Dofs(:,:,:)
!         WRITE(*,*) "=====Element % NodeIndexes(1:n)=",Element % NodeIndexes(1:n)
        
!!!!!!
!!! TEST DECOMMENTE pour 3D
!           WRITE(*,*) "CYCLE if .NOT.isActivePElement(Parent))=",.NOT.isActivePElement(Parent)
! ! !            IF (.NOT.isActivePElement(Parent)) CYCLE
           ptr => ListFind(BC, Name,Found )

! ! ! ******************************************************
! ! ! ******************************************************
! ! ! ******************************************************
! ! ! ******************************************************
! !            WRITE (*,*) "BC, NAME=",Name
! !            bcid=GetBCId( Element )
! !            WRITE (*,*) "BC, id=",bcid
! !            WRITE (*,*) "ptr % FValues",ptr % FValues 
! ! 
! !      ! Get nodes of boundary elements parent and gauss points for boundary
! !     CALL GetElementNodes( Nodes, Element )
! !     IP = GaussPoints( Element )
! ! 
! !     DO t=1,IP % n
! !        stat = ElementInfo( Element, Nodes, IP % u(t), &
! !           IP % v(t), IP % w(t), DetJ, Basis )
! ! 
! !        s = IP % s(t) * DetJ
! ! 
! !        ! Get value of boundary condition
! !        ndd = Element % TYPE % NumberOfNodes
! !        xip = SUM( Basis(1:ndd) * Nodes % x(1:ndd) )
! !        yip = SUM( Basis(1:ndd) * Nodes % y(1:ndd) )
! !        zip = SUM( Basis(1:ndd) * Nodes % z(1:ndd) )
! !        Load = ListGetConstReal( BC, Name, x=xip,y=yip,z=zip )
! !        write(*,*) "Name=",TRIM(Name)," - BcID=",GetBCId( Element )," - Load(",xip,",",yip,",",zip,")=",Load
! !     END DO
! !        
! ! ! ******************************************************
! ! ! ******************************************************
! ! ! ******************************************************
! ! ! ******************************************************           
! ! ! ******************************************************
!         write(*,*) "COUCOU4"           
              ! Number of nodes for this element
              n = Element % TYPE % NumberOfNodes
!         write(*,*) "n = Element % TYPE % NumberOfNodes =",n          

! Get indexes for boundary and values for dofs associated to them
!> Calculate global indexes of boundary dofs for given element and its boundary.  

              CALL getBoundaryIndexesMAN( Solver % Mesh, Element, Parent, gInd, numEdgeDofs )
        
!         write(*,*) "=====gInd=",gInd
! ! ! ! ! ! ! ! !               CALL LocalBcBDOFs( BC, Element, numEdgeDofs, Name, STIFF, Work )
!         write(*,*) "COUCOU6"

       
! ******************************************************
                DO l=1,n
                  nb = x % Perm( gInd(l) )
!                   write(*,*) "nb=",nb
! 20140916 - MODIF GIND PAR NODEINDEX
! ! !                   nb = x % Perm(Element % NodeIndexes(l))
                  

                  nb = Offset + x % DOFs * (nb-1) + DOF
! !                   write(*,*) "nb=",nb
                  IF ( Cond == 0 ) THEN
!                      COUNTOUR <- COUTOUR GIVEN
                    Vector(nb) = VectBCInit(nb)
                  ELSE IF ( Cond == 1 ) THEN
!                      USED AS CALCONDLIM + CALIMP IF Coef==0                     
!                       vector <- coeff X COUNTOUR
!                       vector <- 0
!                       vector <- Lambda X Ud on COUNTOUR
                    IF(VectBCInit(nb).GT.1e-15)THEN
                       Vector(nb) = Coeff * VectBCInit(nb)
                    ELSE
                      Vector(nb) = 0.0_dp
                    ENDIF
!                      write(*,*) TRIM(Name),l," Coeff=",Coeff,"x % Values(nb)=",x % Values(nb)
!                      Vector(nb) = Coeff * x % Values(nb)
!                      write(*,*) TRIM(Name),l,"Vector(nb)=",Vector(nb)
                    
                  ELSE IF ( Cond == 2 ) THEN
!                      USED AS CALIMP without user "load"
!                      IF Value is 0 on given BC then put 0 on current vector
!                      IF Value is not 0 for the BC then don't touch the given vector
                    IF (VectBCInit(nb).LT.1e-15) THEN
                      Vector(nb) = 0.0_dp
                    ENDIF
                  END IF
!                   write(*,*) "Vector(nb)",Vector(nb)
                END DO   
! ******************************************************
          
           ConstantValue =  ptr % PROCEDURE == 0 .AND. &
             ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR
           IF ( ConstantValue ) CYCLE


        END DO
        CurrentModel % CurrentElement => SaveElement
     END DO



     CALL Info('ANM::CalCondLimCoeff','Dirichlet boundary conditions are set for the ANMD RHS', Level=5)
!------------------------------------------------------------------------------

CONTAINS
!------------------------------------------------------------------------------
!> Calculate global indexes of boundary dofs for given element and its boundary.
!------------------------------------------------------------------------------
   SUBROUTINE getBoundaryIndexesMAN( Mesh, Element, Parent, Indexes, indSize )
!!!!!!!
            USE Adaptive
            USE SolverUtils
            USE DefUtils
!!!!!!!
!------------------------------------------------------------------------------
!
!    Type(Mesh_t) :: Mesh
!      INPUT: Finite element mesh containing edges and faces of elements
!
!    Type(Element_t) :: Element
!      INPUT: Boundary element to get indexes for
!
!    Type(Element_t) :: Parent
!      INPUT: Parent of boundary element to get indexes for
!
!    INTEGER :: Indexes(:)
!      OUTPUT: Calculated indexes of boundary element in global system
! 
!    INTEGER :: indSize
!      OUTPUT: Size of created index vector, i.e. how many indexes were created
!        starting from index 1
!------------------------------------------------------------------------------
     IMPLICIT NONE

     ! Parameters
     TYPE(Mesh_t) :: Mesh
     TYPE(Element_t) :: Parent
     TYPE(Element_t), POINTER :: Element
     INTEGER :: indSize, Indexes(:)
     
     ! Variables
     TYPE(Element_t), POINTER :: Edge, Face
     INTEGER :: i,j,n
!      write(*,*) "GBI INIT"

     ! Clear indexes
     Indexes = 0
     n = Element % TYPE % NumberOfNodes

     ! Nodal indexes
     Indexes(1:n) = Element % NodeIndexes(1:n)
!      write(*,*) "Indexes(1:n)=",Indexes(1:n)
!      write(*,*) "Parent % TYPE % DIMENSION=",Parent % TYPE % DIMENSION
     ! Assign rest of indexes if neccessary
     SELECT CASE(Parent % TYPE % DIMENSION)
     CASE (1)
       indSize = n 
     CASE (2)
!      write(*,*) "GBI2D - Element % BDOFs=",Element % BDOFs
        ! Add index for each bubble dof in edge
        DO i=1,Element % BDOFs
           n = n+1
           
           IF (SIZE(Indexes) < n) THEN
              CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
              RETURN
           END IF
!         WRITE(*,*) "Element % PDefs % localNumber=",Element % PDefs % localNumber
!         WRITE(*,*) "Parent % EdgeIndexes(Element % PDefs % localNumber)=",Parent % EdgeIndexes(Element % PDefs % localNumber)
 
           Indexes(n) = Mesh % NumberOfNodes + &
                (Parent % EdgeIndexes(Element % PDefs % localNumber)-1) * Mesh % MaxEdgeDOFs + i
        END DO
     
        indSize = n 
     CASE (3)
!      write(*,*) "GBI3D"
! ! ! ! !         WRITE(*,*) "Element % PDefs % localNumber=",Element % PDefs % localNumber
! ! ! ! !         WRITE(*,*) "Parent % FaceIndexes(previous)=",Parent % FaceIndexes(Element % PDefs % localNumber)
        ! Get boundary face
! ! ! ! !         Face => Mesh % Faces( Parent % FaceIndexes(Element % PDefs % localNumber) )
! ! ! ! !         WRITE(*,*) "GBI3D -  Face % TYPE % NumberOfEdges", Face % TYPE % NumberOfEdges
        ! Add indexes of faces edges 
! ! ! ! !         DO i=1, Face % TYPE % NumberOfEdges
! ! ! ! !            Edge => Mesh % Edges( Face % EdgeIndexes(i) )
           
           ! If edge has no dofs jump to next edge
! ! ! ! !            IF (Edge % BDOFs <= 0) CYCLE

! ! ! ! !            DO j=1,Edge % BDOFs
! ! ! ! !               n = n + 1
              
! ! ! ! !               IF (SIZE(Indexes) < n) THEN
! ! ! ! !                  CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
! ! ! ! !                  RETURN
! ! ! ! !               END IF
              
! ! ! ! !               Indexes(n) = Mesh % NumberOfNodes +&
! ! ! ! !                   ( Face % EdgeIndexes(i)-1)*Mesh % MaxEdgeDOFs + j
! ! ! ! !            END DO
! ! ! ! !         END DO
               
        ! Add indexes of faces bubbles
! ! ! ! !         DO i=1,Face % BDOFs
! ! ! ! !            n = n + 1

! ! ! ! !            IF (SIZE(Indexes) < n) THEN
! ! ! ! !               CALL Warn('DefUtils::getBoundaryIndexes','Not enough space reserved for indexes')
! ! ! ! !               RETURN
! ! ! ! !            END IF

! ! ! ! !            Indexes(n) = Mesh % NumberOfNodes + &
! ! ! ! !                 Mesh % NumberOfEdges * Mesh % MaxEdgeDOFs + &
! ! ! ! !                 (Parent % FaceIndexes( Element % PDefs % localNumber )-1) * Mesh % MaxFaceDOFs + i
! ! ! ! !         END DO        

        indSize = n
     CASE DEFAULT
        CALL Fatal('DefUtils::getBoundaryIndexes','Unsupported dimension')
     END SELECT
!      WRITE(*,*) "Indexes after GBI SIZE(Indexes)=",SIZE(Indexes)
!      DO i=1,SIZE(Indexes)
!        WRITE(*,*) Indexes(i)
!      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE getBoundaryIndexesMAN
!------------------------------------------------------------------------------
  END SUBROUTINE CalCondLimCoeff
!------------------------------------------------------------------------------




 
 
END MODULE FEMtoolBox

!> \} IMNS
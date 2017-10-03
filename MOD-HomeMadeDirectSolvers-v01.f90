!> \ingroup IMNS
!> \}


! # YG@2016
! # yann.guevel@univ-ubs.fr
! # Part of ANM continuation/detection/branch-switching in ELMER FEM
! # Jean-Marc Cadou / Gregory girault / Yann Guevel
! # Institut de Recherche Dupuy de Lome - Universite Bretagne Sud

!------------------------------------------------------------------------------
!>   
!>  DIRECT SOLVER : MUMPS SERIAL
!>  
!>  
!>
!> TODO : 
!>  
!------------------------------------------------------------------------------
MODULE HomeMadeSolvers

!------------------------------------------------------------------------------

  USE SolverUtils
  USE DefUtils

!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
    
CONTAINS




!------------------------------------------------------------------------------
!> Solves local linear system using MUMPS direct solver. 
!>  -> EXTRACTED AND ADAPTED FROM DIRECTSOLVE.SRC
!------------------------------------------------------------------------------
  SUBROUTINE HMumps_SolveLinearSystem( A, x, b, Solver , MUMPSFICH )
! - If first time : Free objects, create and factorize Matrix
! - Fill the RHS
! - Solve
! - Fill x with the Result
! - Every MUMPS information is flushed in UNIT (MUMPSFICH)
!------------------------------------------------------------------------------
     IMPLICIT NONE
! ! ! ! ! ! ! ! ! !     INCLUDE './mpiflocal.h'

     TYPE(Matrix_t) :: A
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp), TARGET :: x(:), b(:)

     INTEGER :: i,MUMPSFICH
     LOGICAL :: Factorize, FreeFactorize, stat

     ! Refactorize local matrix if needed
     Factorize = ListGetLogical( Solver % Values, &
                                'Linear System Refactorize', stat )
     IF (.NOT. stat) Factorize = .TRUE.

     IF (Factorize .OR. .NOT. ASSOCIATED(A % mumpsIDL)) THEN
       CALL HMumpsCreateFactorize(Solver, A , MUMPSFICH)
     END IF

     ! Set RHS
     A % mumpsIDL % NRHS = 1
     A % mumpsIDL % LRHS = A % mumpsIDL % n
     DO i=1,A % NumberOfRows
       A % mumpsIDL % RHS(i) = b(i)
     END DO
     WRITE(*,*) " /\/\/ HMUMPS - ||b|| = ",DSQRT(DOT_PRODUCT(b,b))     
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, b, 1, A % mumpsIDL % RHS, 1)
    WRITE(*,*) " /\/\/ HMUMPS - Solve A x = b"

! ! ! ! ICNTL(9) computes the solution using A or AT  augsys
! ! ! ! Phase: accessed by the host during the solve phase.
! ! ! ! Possible values :
! ! ! ! 1 : AX = B is solved.
! ! ! ! = 1 : AT X = B is solved.
! ! ! ! Default value: 1
! ! ! ! Related parameters: ICNTL(10), ICNTL(11), ICNTL(21), ICNTL(32)
! ! ! ! Remarks: when a forward elimination is performed during the factorization (see ICNTL(32))
! ! ! ! only ICNTL(9)=1 is allowed.
    
    
      ! LEFT      : CollectionMatrix % mumpsIDL % ICNTL(9) = 2
      ! Classic W : CollectionMatrix % mumpsIDL % ICNTL(9) = 1
     A  % mumpsIDL % ICNTL(9) = 1 
     ! SOLUTION PHASE
     A % mumpsIDL % job = 3
     CALL DMumps(A % mumpsIDL)

     ! TODO: If solution is not local, redistribute the solution vector here

     ! Set local solution
     DO i=1,A % NumberOfRows
       x(i)=A % mumpsIDL % RHS(i)
     END DO
     WRITE(*,*) " /\/\/ HMUMPS - ||x|| = ",DSQRT(DOT_PRODUCT(x,x))
     
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, A % mumpsIDL % RHS, 1, x, 1)
!------------------------------------------------------------------------------
  END SUBROUTINE HMumps_SolveLinearSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solves local linear system using MUMPS direct solver. 
!>  -> EXTRACTED AND ADAPTED FROM DIRECTSOLVE.SRC
!------------------------------------------------------------------------------
  SUBROUTINE HMumps_SolveLinearSystemTransp( A, x, b, Solver , MUMPSFICH )
! - If first time : Free objects, create and factorize Matrix
! - Fill the RHS
! - Solve
! - Fill x with the Result
! - Every MUMPS information is flushed in UNIT (MUMPSFICH)
!------------------------------------------------------------------------------
     IMPLICIT NONE
! ! ! ! ! ! ! ! ! !     INCLUDE './mpiflocal.h'

     TYPE(Matrix_t) :: A
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp), TARGET :: x(:), b(:)

     INTEGER :: i,MUMPSFICH
     LOGICAL :: Factorize, FreeFactorize, stat

     ! Refactorize local matrix if needed
     Factorize = ListGetLogical( Solver % Values, &
                                'Linear System Refactorize', stat )
     IF (.NOT. stat) Factorize = .TRUE.

     IF (Factorize .OR. .NOT. ASSOCIATED(A % mumpsIDL)) THEN
       CALL HMumpsCreateFactorize(Solver, A , MUMPSFICH)
     END IF

     ! Set RHS
     A % mumpsIDL % NRHS = 1
     A % mumpsIDL % LRHS = A % mumpsIDL % n
     DO i=1,A % NumberOfRows
       A % mumpsIDL % RHS(i) = b(i)
     END DO
     WRITE(*,*) " /\/\/ HMUMPS - ||b|| = ",DSQRT(DOT_PRODUCT(b,b))     
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, b, 1, A % mumpsIDL % RHS, 1)
    WRITE(*,*) " /\/\/ HMUMPS - Solve A x = b"

! ! ! ! ICNTL(9) computes the solution using A or AT  augsys
! ! ! ! Phase: accessed by the host during the solve phase.
! ! ! ! Possible values :
! ! ! ! 1 : AX = B is solved.
! ! ! ! = 1 : AT X = B is solved.
! ! ! ! Default value: 1
! ! ! ! Related parameters: ICNTL(10), ICNTL(11), ICNTL(21), ICNTL(32)
! ! ! ! Remarks: when a forward elimination is performed during the factorization (see ICNTL(32))
! ! ! ! only ICNTL(9)=1 is allowed.
    
    
      ! LEFT      : CollectionMatrix % mumpsIDL % ICNTL(9) = 2
      ! Classic W : CollectionMatrix % mumpsIDL % ICNTL(9) = 1
     A  % mumpsIDL % ICNTL(9) = 2 
     ! SOLUTION PHASE
     A % mumpsIDL % job = 3
     CALL DMumps(A % mumpsIDL)

     ! TODO: If solution is not local, redistribute the solution vector here

     ! Set local solution
     DO i=1,A % NumberOfRows
       x(i)=A % mumpsIDL % RHS(i)
     END DO
     WRITE(*,*) " /\/\/ HMUMPS - ||x|| = ",DSQRT(DOT_PRODUCT(x,x))
     
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, A % mumpsIDL % RHS, 1, x, 1)
!------------------------------------------------------------------------------
  END SUBROUTINE HMumps_SolveLinearSystemTransp
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Free local Mumps variables and solver internal allocations
!>  -> EXTRACTED AND ADAPTED FROM DIRECTSOLVE.SRC
!------------------------------------------------------------------------------
  SUBROUTINE HMumps_Free_mumpsIDL(A)
!------------------------------------------------------------------------------
       IMPLICIT NONE

       TYPE(Matrix_t) :: A
       INTEGER :: IERR
!     INCLUDE './mpiflocal.h'

       IF (ASSOCIATED(A % mumpsIDL)) THEN
          ! Deallocate Mumps structures
          DEALLOCATE( A % mumpsIDL % irn, A % mumpsIDL % jcn, &
             A % mumpsIDL % a, A % mumpsIDL % rhs)

          ! Free Mumps internal allocations
          A % mumpsIDL % job = -2
          CALL DMumps(A % mumpsIDL)
          DEALLOCATE(A % mumpsIDL)

          A % mumpsIDL => Null()
       END IF
!------------------------------------------------------------------------------
  END SUBROUTINE HMumps_Free_mumpsIDL
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Create and Factorize Matrix with Mumps
! MumpsLocal_Factorize : modified
!>  -> EXTRACTED AND ADAPTED FROM DIRECTSOLVE.SRC
!------------------------------------------------------------------------------
  SUBROUTINE HMumpsCreateFactorize(Solver, A , MUMPSFICH)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A

!     INCLUDE './mpiflocal.h'

    INTEGER :: i, j, n, nz, allocstat, icntlft, ptype, nzloc,IERR,MUMPSFICH
    LOGICAL :: matpd, matsym, nullpiv, stat
    
    IF ( ASSOCIATED(A % mumpsIDL) ) THEN
         CALL HMumps_Free_mumpsIDL(A)
    END IF

    ALLOCATE(A % mumpsIDL)
    ! INITIALIZATION PHASE

    ! Initialize local instance of Mumps
!     ./mpif_stub.h:       PARAMETER (MPI_COMM_SELF=1)
!     A % mumpsIDL % COMM = MPI_COMM_SELF
    A % mumpsIDL % COMM = 1 ! Dummy MPI
    A % mumpsIDL % PAR  = 1 ! Host (=self) takes part in factorization

!   CFD : Matrix is unsymmetric
    A % mumpsIDL % SYM = 0

! _ _ _ _ _ MUMPS Initialize _ _ _ _ _ _ _
    A % mumpsIDL % JOB  = -1 ! Initialize
    CALL DMumps(A % mumpsIDL)
    
    
    ! FACTORIZE PHASE
    ! Set stdio parameters
    A % mumpsIDL % ICNTL(1)  = MUMPSFICH ! 6  ! Error messages to stdout
    A % mumpsIDL % ICNTL(2)  = MUMPSFICH ! 6 ! -1 No diagnostic and warning messages
    A % mumpsIDL % ICNTL(3)  = MUMPSFICH !6 ! -1 No statistic messages
    A % mumpsIDL % ICNTL(4)  = 3  ! Print only error messages
    
    ! Set matrix format
    A % mumpsIDL % ICNTL(5)  = 0 ! Assembled matrix format
    A % mumpsIDL % ICNTL(18) = 0 ! Centralized matrix
    A % mumpsIDL % ICNTL(21) = 0 ! Centralized dense solution phase

    ! Set permutation strategy for Mumps
! ! !  7 : Based on the structural symmetry of the input matrix and on the availability of the numerical values, the value of ICNTL(6) is automatically chosen by the software.
!     ptype = ListGetInteger(Solver % Values, &
!                                 'Mumps Permutation Type', stat)
!     IF (stat) THEN
!       A % mumpsIDL % ICNTL(6) = ptype
!     END IF

! _ _ _ _ _ FILL THE MATRIX _ _ _ _ _ _ _
    ! TODO: Change this if system is larger than local.
    ! For larger than local systems define global->local numbering
    n = A % NumberofRows
    nz = A % Rows(A % NumberOfRows+1)-1
    A % mumpsIDL % N  = n
    A % mumpsIDL % NZ = nz
    
    ! Allocate rows and columns for MUMPS
    ALLOCATE( A % mumpsIDL % IRN(nz), &
          A % mumpsIDL % JCN(nz), &
          A % mumpsIDL % A(nz), STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('HMumpsCreateFactorize', &
            'Memory allocation for MUMPS row and column indices failed.')
    END IF

    ! Set matrix for Mumps (unsymmetric case)
      DO i=1,A % NumberOfRows
         DO j=A % Rows(i),A % Rows(i+1)-1
            A % mumpsIDL % IRN(j) = i
         END DO
       END DO
       ! Set columns and values
       DO i=1,A % mumpsIDL % nz
         A % mumpsIDL % JCN(i) = A % Cols(i)
       END DO
       DO i=1,A % mumpsIDL % nz
         A % mumpsIDL % A(i) = A % Values(i)
       END DO

    icntlft = ListGetInteger(Solver % Values, 'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % mumpsIDL % ICNTL(14) = icntlft
    END IF
    
! _ _ _ _ _ MUMPS Analysis _ _ _ _ _ _ _
    WRITE(*,*) " - - - - HMUMPS - Analysis - - - - - "

    A % mumpsIDL % JOB = 1 ! Perform analysis
    CALL DMumps(A % mumpsIDL)
    CALL Flush(MUMPSFICH)   
    
    ! Check return status
    IF (A % mumpsIDL % INFO(1)<0) THEN
      CALL Fatal('HMumpsCreateFactorize','Mumps analysis phase failed')
    END IF

! _ _ _ _ _ MUMPS factorization _ _ _ _ _ _ _    
    WRITE(*,*) " - - - - HMUMPS - factorization - - - - - "
    A % mumpsIDL % JOB = 2 ! Perform factorization
    CALL DMumps(A % mumpsIDL)
    CALL Flush(MUMPSFICH)
    
    ! Check return status
    IF (A % mumpsIDL % INFO(1)<0) THEN
      CALL Fatal('HMumpsCreateFactorize','Mumps factorize phase failed')
    END IF

    ! Allocate RHS
    ALLOCATE(A % mumpsIDL % RHS(A % mumpsIDL % N), STAT=allocstat)
    IF (allocstat /= 0) THEN
         CALL Fatal('HMumpsCreateFactorize', &
                 'Memory allocation for MUMPS solution vector and RHS failed.' )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE HMumpsCreateFactorize
!-----------------------------------------------------------------



!! MUMPS COMPUTE NULL SPACE BASIS

!* Null pivot (ICNTL(24)) and null space detection ICNTL(25)) improved !for unsymmetric matrices

! - - - - - - - - - RANK deficiency - - - - - - - - - - - - - - - -
!  : Null pivot (ICNTL(24)
! MUMPS gives the possibility to detect null pivots of a matrix during factorization. This option is controlled
! by ICNTL(24). The number of null pivots provides an estimate of the rank deficiency.
! At the solution phase, one of the possible solution of the deficient system AX = B can be computed
! using the control parameter ICNTL(25) (Subsection 3.5).
! It is also possible to compute all or a part of the null vectors (that is the vectors solving AX = 0)
! associated to these null pivots using the same control parameter ICNTL(25).
!
! ICNTL(24 = 1: Null pivot row/column detection.
! CNTL(3) is used to compute the threshold to decide if a pivot row is “null”.
! > 0.0: thres = CNTL(3) ×kApre k
! = 0.0: thres = eps × 10−5 × kApre k  
! < 0.0: thres = |CN T L(3)|

! where Apre is the preprocessed matrix to be factorized (see Equation (5)), eps is the machine precision and k.k is the infinite norm.


! - - - - - - - - - null space basis - - - - - - - - - - - - - - - -
! ICNTL(25) allows the computation of a solution of a deficient matrix and also of a null space basis.
! -1: The complete null space basis is computed.
!
! Remarks: Null space basis computation can be used when a zero-pivot detection option was
! requested (ICNTL(24) 6= 0) during the factorization and when the matrix was found to be deficient
! (INFOG(28) > 0).














END MODULE HomeMadeSolvers
!> \} IMNS
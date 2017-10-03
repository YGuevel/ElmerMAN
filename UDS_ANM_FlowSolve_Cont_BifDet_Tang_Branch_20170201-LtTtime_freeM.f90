!> \ingroup Solvers
!> \{

! # YG@2016
! # yann.guevel@univ-ubs.fr
! # Part of ANM continuation/detection/branch-switching in ELMER FEM
! # Jean-Marc Cadou / Gregory girault / Yann Guevel
! # Institut de Recherche Dupuy de Lome - Universite Bretagne Sud


!------------------------------------------------------------------------------
!   ____  _____  _____     _     ____ __   __  ____  _____   _   _____  _____ 
!  / ___||_   _|| ____|   / \   |  _ \\ \ / / / ___||_   _| / \ |_   _|| ____|
!  \___ \  | |  |  _|    / _ \  | | | |\ V /  \___ \  | |  / _ \  | |  |  _|  
!   ___) | | |  | |___  / ___ \ | |_| | | |    ___) | | | / ___ \ | |  | |___ 
!  |____/  |_|  |_____|/_/   \_\|____/  |_|   |____/  |_|/_/   \_\|_|  |_____|
!   ____  ___   _   _  _____  ___  _   _  _   _    _   _____  ___  ___   _   _ 
!  / ___|/ _ \ | \ | ||_   _||_ _|| \ | || | | |  / \ |_   _||_ _|/ _ \ | \ | |
! | |   | | | ||  \| |  | |   | | |  \| || | | | / _ \  | |   | || | | ||  \| |
! | |___| |_| || |\  |  | |   | | | |\  || |_| |/ ___ \ | |   | || |_| || |\  |
!  \____|\___/ |_| \_|  |_|  |___||_| \_| \___//_/   \_\|_|  |___|\___/ |_| \_|
!------------------------------------------------------------------------------
! ___.   .__  _____                            __  .__               
! \_ |__ |__|/ ____\_ _________   ____ _____ _/  |_|__| ____   ____  
!  | __ \|  \   __\  |  \_  __ \_/ ___\\__  \\   __\  |/  _ \ /    \ 
!  | \_\ \  ||  | |  |  /|  | \/\  \___ / __ \|  | |  (  <_> )   |  \
!  |___  /__||__| |____/ |__|    \___  >____  /__| |__|\____/|___|  /
!      \/                            \/     \/                    \/        
!------------------------------------------------------------------------------
! .__ .___.___..___ __ .___.._..__..  .
! |  \[__   |  [__ /  `  |   | |  ||\ |
! |__/[___  |  [___\__.  |  _|_|__|| \|
!------------------------------------------------------------------------------
! .__ .__ .__..  . __ .  .   __..  .._..___. __ .  .._..  ..__ 
! [__)[__)[__]|\ |/  `|__|  (__ |  | |   |  /  `|__| | |\ |[ __
! [__)|  \|  || \|\__.|  |  .__)|/\|_|_  |  \__.|  |_|_| \|[_./
!------------------------------------------------------------------------------
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
!
! TODO LIST:
!201504
! Pade V2                                       -> EnTest
! BifOrLimitPoint : pade geom                   -> EnCours
! Mumps V5 : bif multiple ; null basis ....     -> Trick -fPic
! Residual correction : ElmerLoop 1step 1step.. : no possiblo : equation NS + NS forbiden
! BifMulitple : NB tangents, NB mode , NB leftMode, .... ABE CF Zhen Mei D4

! CHANGE LOG:
! - - - - - - 
! 2016 - 11 : LtC FreeMatrix if Newton correction ;)
! 2016 - 10 : AlphaBif : bifurcation behind : new Amax with new series
! 2016 - 10 : Newton Correction Hard coded -> to be moved in Modules?
! 2015 - 05 : Memory ok, RIKS correction ok, Assembly AUGSYS faster
! 2015 - 04 : Bifurcated Branches : 
! 2015 - 02 : MUMPS V4 ok
! 2015 - 02 : ABE Tangents : W PSI AugSys
! 2015 - jan : steady state bifurcation indicator cadou 2006 guevel 2011 ( Tri 1996 fluid, vanucci98 boutyour93 )
! 2014 - oct modularity   : fonction utilitaires dans un module
! elmerf90 -c ELMERModResOut.f90 ANMModToolBoxV00.f90
! elmerf90 ANMModToolBoxV00.o ELMERModResOut.o FLowSolveANMPathFollowingWithToolBox20141024.f90 -o FlowSolveANMPathFollowing_V03_3D_MOD.so
! cp FlowSolveANMPathFollowing_V03_3D_MOD.so *.mod CALCS/

! 2014 - sept FQopt OK    : optimisation FQ en conservant les Ur et Grad Up-r déjà calculé.
! 2014 - sept 3D ok       :  ajout routine getBoundaryIndexesMAN en commentant la partie 3D.....
! 2014 - juin Tual Allain : stage 2D
!
!


!------------------------------------------------------------------------------
!  _____     _______ _    _   ______ ____  _      _      ______          _______ _   _  _____ 
! |  __ \ /\|__   __| |  | | |  ____/ __ \| |    | |    / __ \ \        / /_   _| \ | |/ ____|
! | |__) /  \  | |  | |__| | | |__ | |  | | |    | |   | |  | \ \  /\  / /  | | |  \| | |  __ 
! |  ___/ /\ \ | |  |  __  | |  __|| |  | | |    | |   | |  | |\ \/  \/ /   | | | . ` | | |_ |
! | |  / ____ \| |  | |  | | | |   | |__| | |____| |___| |__| | \  /\  /   _| |_| |\  | |__| |
! |_| /_/    \_\_|  |_|  |_| |_|    \____/|______|______\____/   \/  \/   |_____|_| \_|\_____|                                                                                            
!------------------------------------------------------------------------------
   SUBROUTINE ANMPathFollowing( Model,Solver,dt,TransientSimulation)
!------------------------------------------------------------------------------
    USE NavierStokes
    USE NavierStokesGeneral
    USE NavierStokesCylindrical
    USE Adaptive
    USE DefUtils
    USE FreeSurface
    
!------------------------------------------------------------------------------    
! MODULE FAIT MAISON
    USE ANMToolBox
    USE ELMERModResOut
    USE DiscreteOperators
    USE HomeMadeSolvers
    USE FlowSolveCorrection
    USE ANMSingularities
    
!------------------------------------------------------------------------------
    IMPLICIT NONE

     TYPE(Model_t) :: Model
     TYPE(Solver_t), TARGET :: Solver

     REAL(KIND=dp) :: dt
     LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Matrix_t),POINTER :: StiffMatrix

     INTEGER :: i,j,k,n,nb,nd,t,iter,LocalNodes,istat,ii,jj,kk,IO,NDL

     TYPE(ValueList_t),POINTER :: Material, BC, BodyForce, Equation
     TYPE(Element_t),POINTER :: Element

     REAL(KIND=dp) :: RelativeChange,UNorm,Gravity(3),AngularVelocity(3), &
       Tdiff,s,Relaxation,NewtonTol,NonlinearTol, &
       ReferencePressure=0.0, SpecificHeatRatio, &
       PseudoCompressibilityScale=1.0, NonlinearRelax,FreeSTol

     INTEGER :: NSDOFs,NewtonIter,NonlinearIter,FreeSIter

     TYPE(Variable_t), POINTER :: DensitySol, TimeVar
     TYPE(Variable_t), POINTER :: FlowSol, TempSol, MeshSol

     INTEGER, POINTER :: FlowPerm(:),TempPerm(:), MeshPerm(:)
     REAL(KIND=dp), POINTER :: FlowSolution(:), Temperature(:), &
        gWork(:,:), ForceVector(:), LayerThickness(:), &
           SurfaceRoughness(:),MeshVelocity(:)

     REAL(KIND=dp), POINTER :: TempPrev(:)

     LOGICAL :: Stabilize,NewtonLinearization = .FALSE., GotForceBC, GotIt, &
                  MBFlag, Convect  = .TRUE., NormalTangential, RelaxBefore, &
                  divDiscretization, GradPDiscretization, ComputeFree=.FALSE., &
                  Transient, Rotating

     CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, StabilizeFlag, VarName
     CHARACTER(LEN=MAX_NAME_LEN) :: LocalCoords, FlowModel
     INTEGER :: CompressibilityModel, ModelCoords, ModelDim
     INTEGER :: body_id,bf_id,eq_id,DIM
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)


     INTEGER, SAVE :: Timestep, SaveTimestep=-1
     REAL(KIND=dp), ALLOCATABLE, SAVE :: pDensity0(:), pDensity1(:)
     REAL(KIND=dp), ALLOCATABLE, SAVE :: PseudoPressure(:),PSolution(:),DensityTMP(:),NodalMuTMP(:)

     LOGICAL :: AllocationsDone = .FALSE., FreeSurfaceFlag,             &
         PseudoPressureExists, PseudoCompressible, Bubbles,             &
         Porous =.FALSE., PotentialForce=.FALSE., Hydrostatic=.FALSE.,  &
         MagneticForce =.FALSE., UseLocalCoords, PseudoPressureUpdate

! BEFORE :          
!       REAL(KIND=dp) :: at,at0,at1,totat,st,totst,CPUTime,RealTime  
! ! ! ! FLowSolveANMPathFollowingWithToolBox_LtQResOK_TANGENTs_20150227-MUMPS_clssc.f90:135.56:
! ! ! !       REAL(KIND=dp) :: at,at0,at1,totat,st,totst,CPUTime,RealTime
! ! ! ! FLowSolveANMPathFollowingWithToolBox_LtQResOK_TANGENTs_20150227-MUMPS_clssc.f90:66.8:
! ! ! !     USE NavierStokes
! ! ! ! Error: Symbol 'cputime' at (1) conflicts with symbol from module 'loadmod', use-associated at (2)
! AFTER  : 
      REAL(KIND=dp) :: at,at0,at1,totat,st,totst

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Déclaration test, lecture des données dans le .SIF
!------------------------------------------------------------------------------

        REAL(KIND=dp) :: col , row , val
        REAL(KIND=dp), ALLOCATABLE :: Uelex(:), Ueley(:), Uelez(:), PressureMan(:) 
        REAL(KIND=dp), ALLOCATABLE :: Velex(:), Veley(:), Velez(:)
                                      

        LOGICAL :: arg ! on en a besoin pour les fonctions GET"type_de_données"
        
        ! OMan: ordre MAN, r: ordre intermédiaire pour le calcul du 2nd membre MAN
        INTEGER :: OMan, ordreman, Step, ANMStep, r, p
        ! DIAG BIF
        REAL(KIND=dp) :: CritrDiagBif,CDBx,CDBy,CDBz
        INTEGER :: NODEOUT,NdOutCOMP,NBSolByStep

        REAL(KIND=dp) ,ALLOCATABLE, TARGET :: UMan(:,:),VMan(:,:),VTMPOrtho(:,:),VSOLTMP(:),VChap(:,:)
        REAL(KIND=dp) ,ALLOCATABLE :: FQMan(:), UTemp(:), lambda(:),      &
                                      NormUMan(:), NormFQ(:), U0(:), NoPressure(:),     &
                                      FQManTemp(:),  FlowSolutionTemp(:), &
                                      FlowSolution_Init(:),                             &
                                      USAV(:,:,:,:), GradSAV(:,:,:,:,:), &
                                      U0Pad(:)
        REAL(KIND=dp) :: Ltemp,DirectionCont
        REAL(KIND=dp) ,ALLOCATABLE :: VORIENT(:)
        INTEGER :: SENSPARCOURS
        LOGICAL :: SYMRESTART


                                      
!------------------------------------------------------------------------------                                      
! PADE
        REAL(KIND=dp) ,ALLOCATABLE :: AlphaPad(:,:), DN(:) ! PADE
        REAL(KIND=dp) :: amaxPad,PPole,L0Pad,L0Padcrit

!------------------------------------------------------------------------------        
! Steady Bifurcation indicator
        REAL(KIND=dp) ,ALLOCATABLE :: INDDU0init(:), INDDU0(:), INDDUserie(:,:), INDMuserie(:), &
                                      INDFNL(:), INDRHS(:), INDFalea(:),INDYf(:),INDXk(:),      &
                                      INDFNLtemp(:) ,                                           &
                                      INDDUSAV(:,:,:,:), GradINDDUSAV( :,:,:,:,: )  ,           &
                                      INDUelex(:), INDUeley(:), INDUelez(:)     
        REAL(KIND=dp) :: INDMu0init, INDMu0 , INDDENUM , INDMU0NUM , INDFaleaMagnitude
        REAL(KIND=dp) :: amaxpolyUL, amaxpolyVL, amaxpolyINDDU , amaxpolyBr
        
        REAL(KIND=dp), POINTER  :: SaveValues(:)
        LOGICAL :: Found,ResOUTornot

!------------------------------------------------------------------------------
! COCHELIN MEDALE 2012
        REAL(KIND=dp) ,ALLOCATABLE :: calcalphabif(:) ! pour calculer la distance à la bifurcation::: Cochelin Medale
        REAL(KIND=dp) :: EPSILON1, EPSILON2, sumcrit(2) , LambdaCritPad
        ! V de U = V P <=> U = V 0 
        REAL(KIND=dp) ,ALLOCATABLE :: VManProp(:,:), LambdaVProp(:), VBifMode(:) , VBif(:)
        REAL(KIND=dp) :: LVBif, alphaVBIF
        ! U = V P
        REAL(KIND=dp) ,ALLOCATABLE :: UManProp(:,:), LambdaUProp(:), UBifMode(:), UBif(:)
        REAL(KIND=dp) :: LUBif, alphaUBIF
        ! - - 
        REAL(KIND=dp) ,ALLOCATABLE :: BifMode(:)
        REAL(KIND=dp) :: alphaBIF
!
        INTEGER :: m, N1, N2, N3
        INTEGER :: NBEL,NBBE

        REAL(KIND=dp) :: TolMan, amax, lambda0, L2UMan, L2FQMan, normMAN,&
                         norm1, norm2, Coeff
        LOGICAL :: Ex_Logical
        CHARACTER(LEN=MAX_NAME_LEN) :: Filename, Filename2, Name

!        SAVE lambda0
        SAVE FlowSolution_Init
!  
        CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix
        
        REAL(KIND=dp) :: alphaBIFSENSPARCOURS  ! for limit point case with inversion of direction...
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! BRANCH SWITCHING
! TANGENT Lt_c ModeBif W ModeGauche Q 
! Ut = a * W + b * ModeBif
! PsiL * Q(Ut,Ut) = 0  => ABE
! MatrixVectorMultiply( A,u,v ) V = Au
! TransposeMatrixVectorMultiply( A,u,v ) V = AT u
        REAL(KIND=dp) ,ALLOCATABLE :: abeW(:), abePsiL(:)
        REAL(KIND=dp) ,ALLOCATABLE :: abeQAA(:), abeQBB(:), abeQAB(:), abeQBA(:)
        REAL(KIND=dp) ,ALLOCATABLE :: Vtmp(:)
        REAL(KIND=dp) ,ALLOCATABLE :: UBranch(:,:), LBranch(:)
        REAL(KIND=dp) ,ALLOCATABLE :: abeUT(:,:), abeLT(:)        
        REAL(KIND=dp)              :: abeCoA,abeCoB,abeCoC,Discr,DeuxA
        REAL(KIND=dp)              :: abeL1,abeL2, abeE1,abeE2,abeL1temp,abeL2temp
        REAL(KIND=dp)              :: XNW,XNBifM ,XPWW,XPDD,XRES,XTMP,XTEST
        INTEGER                    :: KKI,IT,NBTANGENTS , BFN
        TYPE(Matrix_t), POINTER    :: StiffAugW      
        CHARACTER(LEN=MAX_NAME_LEN) :: StabBis
!------------------------------------------------------------------------------
        INTEGER :: MUMPSFICH, ABETFICH, RESIDUALFICH, BIFSTADIAGFICH 
        INTEGER :: CRITGEOMFICH, PADEFICH , BIFSTAINDICFICH, TIMEFILEFICH,DBSKPFICH
!------------------------------------------------------------------------------
        REAL(KIND=dp) ,ALLOCATABLE,TARGET :: ResVec(:), ResVecTMP(:)
        REAL(KIND=dp)              :: FNORMG, ResNorm,NORMUO
!------------------------------------------------------------------------------
        INTEGER :: SaveCount
        LOGICAL :: Binary,SaveAll, BIFSTAGEOMDETECT
!------------------------------------------------------------------------------
        TYPE(Variable_t), POINTER :: Var 
        REAL(KIND=dp), POINTER    :: ResidualPointer(:)
        REAL(KIND=dp) :: ManualResidualU0j, HMResidualU0jp1, HkResidualNorm,    &
                         ManualResidualU0jstab, HMResidualU0jp1NoC
!------------------------------------------------------------------------------
        REAL(KIND=dp) :: INDCADMu0, INDCADNum, INDCADDenum 
        REAL(KIND=dp) :: INDVANMu0, INDVANNum, INDVANDenum  
        REAL(KIND=dp) :: INNORMMu0, INNORMNum, INNORMDenum 
        REAL(KINd=dp), ALLOCATABLE :: INDDU0initA(:)        
        REAL(KIND=dp) :: INDCADMu0A, INDCADNumA, &
                         INDVANMu0A, INDVANNumA, &
                         INNORMMu0A, INNORMNumA, INDMu0initA
!------------------------------------------------------------------------------
        LOGICAL :: FACTOTODO
!------------------------------------------------------------------------------
! DERIVATE OF SERIES
        REAL(KINd=dp), ALLOCATABLE :: Uder(:)        
        REAL(KIND=dp) :: Lder
!------------------------------------------------------------------------------
!NORME TEST
        REAL(KIND=dp) :: NLeftMode,PLB,NFaug,NTMP1,NTMP2,NW
! PROSCA TEST        
        REAL(KIND=dp) :: LtB,LMW,LtW,BMW, ZPF
!------------------------------------------------------------------------------
        LOGICAL :: SYMBRANCH,PITCHFORKbifDETECTED
!------------------------------------------------------------------------------        
        REAL(KIND=dp) :: RiksDL,RiksNTMP, PrecRIKS
        LOGICAL       :: RIKSCORRECTION, NEWTCORRECTION,CallNewtonCorr
        INTEGER       :: RIC,RICMAX, NEWTMAX,NEWTITER
!------------------------------------------------------------------------------
        REAL(KIND=dp) ,ALLOCATABLE :: URestart(:)
        REAL(KIND=dp)              :: LRestart
!------------------------------------------------------------------------------        
        REAL(KIND=dp)  :: NORMtest1,NORMtest2
        INTEGER        :: idbs
!------------------------------------------------------------------------------        
        REAL(KIND=dp) :: TIC, TAC,TOC,  NORMTMP
!------------------------------------------------------------------------------        


        
        
        
        
        
        TOC = CPUTime()

        

!------------------------------------------------------------------------------
!    Get variables needed for solving the system
!------------------------------------------------------------------------------

       IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN


!______________________________________________________________________________
!------------------------------------------------------------------------------
!    Lecture des paramètres MAN dans le .SIF
       WRITE(*,*) '----------------------------------------------------------------------'
       WRITE(*,*) ''
       WRITE(*,*) '  -  -  -  -  -  -  ANM PARAMETERS FROM SIF FILE  -  -  -  -  -  -  -'
       WRITE(*,*) ''
       
       ANMStep           =   GetInteger( Solver%Values , 'ANM Step', arg )         ! NBSTEPs       
       OMan              =   GetInteger( Solver%Values , 'ANM Order', arg )        ! ANM ORDER
       TolMan            = GetConstReal( Solver%Values , 'ANM Tolerance', arg )    ! PREC FOR Amax
!    
       NODEOUT           =   GetInteger( Solver%Values , 'ANM NodeOut', arg )      ! 4 DIAG BIF
       NdOutCOMP         =   GetInteger( Solver%Values , 'ANM NodeOutComp', arg )  ! useless now
       ResOUTornot       =   GetLogical( Solver%Values , 'ANM ResOUTornot', arg )  ! Parav or not
!        
       NBSolByStep       =   GetInteger( Solver%Values , 'ANM NBSolByStep', arg )  ! NBsol on begin Bif branches
       BIFSTAGEOMDETECT  = .FALSE.
       BIFSTAGEOMDETECT  = GetLogical( Solver%Values   , 'ANM Steady Bifurcation GeomDetection', arg ) ! Activation of the geom criterion to detect bifsta
       EPSILON1          = GetConstReal( Solver%Values , 'ANM Steady Bifurcation GeomDetection EPSA', arg )
       EPSILON2          = GetConstReal( Solver%Values , 'ANM Steady Bifurcation GeomDetection EPSV', arg )
       SYMRESTART        =   GetLogical( Solver%Values , 'ANM Bifurcation SYM BRANCH RESTART', arg )  ! FB FB FB
       INDFaleaMagnitude = GetConstReal( Solver%Values , 'ANM Steady Bifurcation Indicator Perturbation', arg )
       PrecRIKS          = GetConstReal( Solver%Values , 'ANM Steady RIKS Correction Tolerance', arg )
       CallNewtonCorr    =   GetLogical( Solver%Values , 'ANM Steady USE NEWTON Correction', arg )  ! FB FB FB
       NEWTMAX           =   GetInteger( Solver%Values , 'ANM Steady NEWT Correction MAXITER', arg )  ! FB FB FB
       
       
       
       
       WRITE(*,*) 'ANM Step       = ', ANMStep
       WRITE(*,*) 'ANM Order      = ', OMan
       WRITE(*,*) 'ANM Tolerance  = ', TolMan
       WRITE(*,*) 'ANM NodeOut    = ', NODEOUT
!        WRITE(*,*) 'initial lambda        = ', lambda0
       WRITE(*,*) ''
       WRITE(*,*) '----------------------------------------------------------------------'       
!      Fin lecture du .SIF POUR LES MOT CLES ANM
!______________________________________________________________________________
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!    Check for local coordinate system
     LocalCoords = GetString( Solver % Values, 'Solver Coordinate System', &
          UseLocalCoords )

     IF ( UseLocalCoords ) THEN
        ModelCoords = Coordinates
        ModelDim = Model % DIMENSION
        SELECT CASE ( LocalCoords )
           CASE( 'cartesian 2d' )
              Coordinates = 1
              Model % DIMENSION = 2
              CALL Info( 'ManFlowSolve', 'Solver Coordinate System is cartesian 2d', LEVEL=31 )
           CASE( 'cartesian 3d' )
              Coordinates = 1
              Model % DIMENSION = 3
              CALL Info( 'ManFlowSolve', 'Solver Coordinate System is cartesian 3d', LEVEL=31 )
           CASE( 'axi symmetric' )
              Coordinates = 4
              Model % DIMENSION = 2
              CALL Info( 'ManFlowSolve', 'Solver Coordinate System is axi symmetric', LEVEL=31 )
           CASE( 'cylindric symmetric' )
              Coordinates = 3
              Model % DIMENSION = 3
              CALL Info( 'ManFlowSolve', 'Solver Coordinate System is cylindric symmetric', LEVEL=31 )
           CASE DEFAULT
              CALL Warn( 'ManFlowSolve', 'Solver coordinate system not recognised, using original' )
        END SELECT
     END IF

     ! check for Flow model, one of 'full', 'no convection', 'stokes':
     ! ---------------------------------------------------------------
     Transient = TransientSimulation
     Convect = .TRUE.
     FlowModel = GetString( GetSolverParams(), 'Flow Model', Gotit )
     
     SELECT CASE(FlowModel)
     CASE('no convection')
       Convect = .FALSE.
     CASE('stokes')
       Convect = .FALSE.
       Transient = .FALSE.
     CASE DEFAULT
       FlowModel = 'full'
     END SELECT

!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -     
     DIM = CoordinateSystemDimension()

     FlowSol => Solver % Variable
     NSDOFs         =  FlowSol % DOFs
     FlowPerm       => FlowSol % Perm
     FlowSolution   => FlowSol % Values
     
     VarName = GetVarName(FlowSol)
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

     LocalNodes = COUNT( FlowPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

     TempSol => VariableGet( Solver % Mesh % Variables, 'Temperature' )
     IF ( ASSOCIATED( TempSol ) ) THEN
       TempPerm     => TempSol % Perm
       Temperature  => TempSol % Values
       IF( Transient ) THEN
         IF ( ASSOCIATED(TempSol % PrevValues) ) TempPrev => TempSol % PrevValues(:,1)
       END IF
     END IF

     MeshSol => VariableGet( Solver % Mesh % Variables, 'Mesh Velocity')
     NULLIFY( MeshVelocity )
     IF ( ASSOCIATED( MeshSol ) ) THEN
       MeshPerm     => MeshSol % Perm
       MeshVelocity => MeshSol % Values
     END IF

!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -     
     DensitySol => VariableGet( Solver % Mesh % Variables, 'Density' )
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
     StiffMatrix => Solver % Matrix
     ForceVector => StiffMatrix % RHS
     UNorm = Solver % Variable % Norm
!  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN

       N = Solver % Mesh % MaxElementDOFs
!------------------------------------------------------------------------------       
       IF( AllocationsDone ) THEN
          DEALLOCATE(                                &
               UMan, VMan, VTMPOrtho,VSOLTMP,Vchap,  &
               FQMan, UTemp, Lambda,   U0Pad,        &
               U0, NormUMan, NormFQ,         &
               NoPressure, FQManTemp,    &
               FlowSolutionTemp, FlowSolution_Init,  &
               calcalphabif, VORIENT,                &
               AlphaPad,DN,                          &
               VManProp,LambdaVProp, VBif, VBifMode, &
               UManProp,LambdaUProp, UBif, UBifMode, &
               UBranch, LBranch ,                    &
               abeUT, abeLT ,                        &
               INDDU0init,                           &
               INDDU0initA,                           &
               INDDU0,                               &
               INDDUserie,                           &
               INDMuserie,                           &
               INDFNL,                               &
               INDRHS,                               &
               INDFalea,                             &
               INDYf,                                &
               INDXk,                                &
               INDFNLtemp,                           &   
               PSolution,                            &
               DensityTMP,  NodalMuTMP,              &
               Indexes,                              &               
               abeW, abePsiL, abeQAA, abeQBB,        &               
               abeQAB, abeQBA,                       & 
               BifMode, Vtmp,                        &
               ResVec,ResVecTMP,                     &
               Uder,URestart,                        &
               STAT=istat )
       END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
       IF( AllocationsDone ) THEN
         DEALLOCATE(                                   &
               Uelex , Ueley , Uelez,                  &
               Velex , Veley , Velez,                  &
               INDUelex , INDUeley , INDUelez,         &
               PressureMan,                            &
               USAV,  GradSAV,                         &
               INDDUSAV, GradINDDUSAV,                 &  
               STAT=istat )
       END IF
       
! Nombre d'élement en boundary condition
       NBBE=GetNOFBoundaryElements()
! nombre d'élement actif en solver Navier stokes       
       NBEL=GetNOFActive()
       ALLOCATE( Uelex(N), Ueley(N), Uelez(N),              &
                 Velex(N), Veley(N), Velez(N),              &       
                 PressureMan(N),                            &
                 USAV   ( NBEL , N , OMan , 3  ) ,          &
                 GradSAV( NBEL , N , OMan , 3 , 3 ),        &
                 INDUelex(N) , INDUeley(N) , INDUelez(N),   &               
                 INDDUSAV   ( NBEL , N , OMan+1 , 3  ) ,    &
                 GradINDDUSAV( NBEL , N , OMan+1 , 3 , 3 ), &
                 DensityTMP(N),                             &
                 NodalMuTMP(N),                             &
                 Indexes( N ),                              &                                  
                 STAT=istat )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------



       NDL = StiffMatrix % NumberOfRows
       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - "
       WRITE(*,*) " NDL = ",NDL
       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - "
       ALLOCATE( U0( NDL ), FlowSolution_Init( NDL ), NoPressure( NDL ),                   &
                 VORIENT( NDL ),                                                           &
                 UTemp( NDL ), FlowSolutionTemp( NDL ), VSOLTMP( NDL ), Vtmp(NDL),         &
                 FQMan( NDL ),FQManTemp( 2*NSDOFs*N ),                                     &
                 NormUMan( OMan ), NormFQ( OMan-1 ),                                       &
                 U0Pad( NDL ), AlphaPad(OMan,OMan),   DN(OMan),                            & 
                 ResVec(NDL), ResVecTMP(NDL),                                              &
                 UBif( NDL ),VBif( NDL ), Uder(NDL),                                       &
                 BifMode(NDL), UBifMode(NDL), VBifMode(NDL), abeW(NDL), abePsiL(NDL),      &
                 UMan( NDL , OMan ),  VMan( NDL , OMan ),  Lambda( OMan ),                 &
                 calcalphabif( 3 ),                                                        &
                 UManProp( NDL , OMan ), VManProp( NDL , OMan ) , VTMPOrtho( NDL , OMan ), &
                 LambdaUProp( NDL ), LambdaVProp( NDL ),                                   &
                 abeQAA(NDL), abeQBB(NDL), abeQAB(NDL), abeQBA(NDL),                       &
                 UBranch( NDL , OMan ),LBranch( NDL ),                                     &
                 URestart( NDL ),                                                          &
                 INDDU0(NDL),INDDU0init(NDL),INDDU0initA(NDL),                             &
                 INDDUserie(NDL , OMan+1),INDMuserie( OMan+1),                             &
                 INDFNL(NDL), INDRHS(NDL), INDFalea(NDL), INDYf(NDL), INDXk(NDL),          &
                 INDFNLtemp(NDL),                                                          &
                 PSolution( SIZE( FlowSolution ) ),                                        &  
                 STAT=istat )
!------------------------------------------------------------------------------

!______________________________________________________________________________
!------------------------------------------------------------------------------
       
!! VECTOR noPressure in ordre to extract velocity from FlowSolution vector by scalar product projection
!     NoPressure = (1, 1, 0, ... ) to avoid the pressure in order to calculate the velocity norm
       NoPressure = 1._dp
!        NoPressure(DIM+1:NDL:NSDOFs)=0._dp
       DO i=NSDOFs,NDL,NSDOFs
         NoPressure(i) = 0.0_dp
       END DO
       
       
!------------------------------------------------------------------------------       
       AlphaBIF = 0._dp       
       BFN = 0
!------------------------------------------------------------------------------
       
       
! !______________________________________________________________________________  
! .-----.-----.----.|  |_.--.--.----.|  |--.---.-.|  |_|__|.-----.-----.
! |  _  |  -__|   _||   _|  |  |   _||  _  |  _  ||   _|  ||  _  |     |
! |   __|_____|__|  |____|_____|__|  |_____|___._||____|__||_____|__|__|
! |__|                                                                  
!   ___                        
! .'  _|.-----.----.----.-----.
! |   _||  _  |   _|  __|  -__|
! |__|  |_____|__| |____|_____|
!                     
!      FORCE ALEATOIRE POUR LE CALCUL DE L INDICATEUR DE BIFURCATION       
       INDFalea = INDFaleaMagnitude

!        INDFalea = 0.001_dp
!      FORCE ALEATOIRE POUR INDICATEUR
!      ON UTILISE NoPressure comme vecteur sans pression....
!        INDFalea = NoPressure * 10.
!      CONDITON FORCE Nulle sur les contours : chargement et blocage
! TEST SANS PRESSION
       CALL PressureFree( INDFalea, NoPressure, INDFalea )
       CALL CalCondLimCoeff( INDFalea, 0.0_dp, 1, FlowSolution )
       WRITE(*,*) "IND - NORME INDFalea = ",DSQRT(DOT_PRODUCT(INDFalea,INDFalea))

! ! ! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$                     
! ! ! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
! ! !        SaveValues => FlowSolution
! ! !        
! ! !        FlowSolution = INDFalea
! ! !        FlowSol % Values => FlowSolution
! ! !        
! ! !        FilePrefix = 'INDICAT-INDFalea-'
! ! !        CALL ListAddString( Solver % Values,'Output File Name',FilePrefix)
! ! !        
! ! !        CALL ResultOutputSolver( Model,Solver,dt,Transient )
! ! !        
! ! !        FlowSolution = SaveValues
! ! !        FlowSol % Values => FlowSolution
! ! !        
! ! !        write(*,*) "-ListGetString -F ",ListGetString( Solver % Values, 'Output File Name', Found )
! ! !        CALL ListRemove( Solver % Values,'Output File Name')
! ! ! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$       
! ! ! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$                     
       
!______________________________________________________________________________       
       
!______________________________________________________________________________
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! ! ! !        Drag = 0.0d0

       PseudoPressureExists = .FALSE.
       DO k=1,Model % NumberOfMaterials
         Material => Model % Materials(k) % Values
         CompressibilityFlag = ListGetString( Material, &
             'Compressibility Model', GotIt)
         IF (gotIt .AND. CompressibilityFlag == 'artificial compressible') THEN
            PseudoPressureExists = .TRUE.
         END IF
       END DO

       IF ( PseudoPressureExists ) THEN
          IF ( AllocationsDone ) THEN
             DEALLOCATE( PseudoPressure )
          END IF
          n = SIZE( FlowSolution ) / NSDOFs
          ALLOCATE( PseudoPressure(n),STAT=istat ) 
       END IF

       IF ( istat /= 0 ) THEN
         CALL Fatal( 'ManFlowSolve','Memory allocation error, Aborting.' )
       END IF

!------------------------------------------------------------------------------

       AllocationsDone = .TRUE.
     END IF
!------------------------------------------------------------------------------

     TimeVar => VariableGet( Solver % Mesh % Variables, 'Timestep')
     Timestep = NINT(Timevar % Values(1))
     IF ( SaveTimestep /= Timestep ) THEN
       IF ( ALLOCATED(pDensity0) ) pDensity0 = pDensity1
       SaveTimestep=Timestep
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

     gWork => ListGetConstRealArray( Model % Constants,'Gravity',GotIt)
     IF ( GotIt ) THEN
       Gravity = gWork(1:3,1)*gWork(4,1)
     ELSE
       Gravity    =  0.00D0
       Gravity(2) = -9.81D0
     END IF
!------------------------------------------------------------------------------

     Bubbles   = ListGetLogical( Solver % Values,'Bubbles',GotIt )
     Stabilize = ListGetLogical( Solver % Values,'Stabilize',GotIt )

     StabilizeFlag = ListGetString( Solver % Values, &
           'Stabilization Method', GotIt )
     IF ( .NOT. GotIt ) THEN
       IF ( Stabilize ) THEN
          StabilizeFlag = 'stabilized'
       ELSE IF ( Bubbles  ) THEN
          StabilizeFlag = 'bubbles'
       ELSE
          StabilizeFlag = 'stabilized'
       END IF
     END IF

     IF ( StabilizeFlag == 'bubbles' ) Bubbles = .TRUE.

     DivDiscretization = ListGetLogical( Solver % Values, &
              'Div Discretization', GotIt )

     GradPDiscretization = ListGetLogical( Solver % Values, &
              'Gradp Discretization', GotIt )

     NonlinearTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance',minv=0.0d0 )

     NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance', minv=0.0d0 )

     NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations', minv=0 )
     IF ( NewtonIter == 0 ) NewtonLinearization = .TRUE.


     IF (GetLogical( GetSolverParams(), &
         'Nonlinear System Reset Newton',  GotIt)) NewtonLinearization=.FALSE.

     NonlinearIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Max Iterations', minv=0 )

     IF ( .NOT. ListCheckPresent( Solver % Values, &
        'Nonlinear System Norm Dofs' ) ) THEN
       CALL ListAddInteger( Solver % Values, 'Nonlinear System Norm DOFs', NSDOFs-1 )
     END IF

     FreeSTol = ListGetConstReal( Solver % Values, &
        'Free Surface After Tolerance', GotIt, minv=0.0d0 )
     IF ( .NOT. GotIt ) FreeSTol = HUGE(1.0d0)

     FreeSIter = ListGetInteger( Solver % Values, &
        'Free Surface After Iterations', GotIt, minv=0 )
     IF ( .NOT. GotIt ) FreeSIter = 0
!------------------------------------------------------------------------------
!    We do our own relaxation...
!------------------------------------------------------------------------------
     NonlinearRelax = GetCReal( Solver % Values, &
        'Nonlinear System Relaxation Factor', GotIt )
     IF ( .NOT. GotIt ) NonlinearRelax = 1.0d0

     CALL ListAddConstReal( Solver % Values, &
            'Nonlinear System Relaxation Factor', 1.0d0 )

     IF ( NonlinearRelax /= 1._dp ) &
       CALL ListAddLogical( Solver % Values, 'Skip Compute Nonlinear Change', .TRUE. )
!------------------------------------------------------------------------------
!    Check if free surfaces present
!------------------------------------------------------------------------------
     FreeSurfaceFlag = .FALSE.
     DO i=1,Model % NumberOfBCs
       FreeSurfaceFlag = FreeSurfaceFlag .OR. ListGetLogical( &
          Model % BCs(i) % Values,'Free Surface', GotIt )
       IF ( FreeSurfaceFlag ) EXIT
     END DO

     CALL CheckCircleBoundary()
!------------------------------------------------------------------------------

     totat = 0.0d0
     totst = 0.0d0

     ! Initialize the pressure to be used in artificial compressibility 
     IF(PseudoPressureExists) THEN
       PseudoPressure = FlowSolution(NSDOFs:SIZE(FlowSolution):NSDOFs)

       WRITE(Message,'(A,T25,E15.4)') 'PseudoPressure mean: ',&
           SUM(PseudoPressure)/SIZE(PseudoPressure)
       CALL Info('ManFlowSolve',Message,Level=5)

       PseudoCompressibilityScale = ListGetConstReal( Model % Simulation, &
           'Artificial Compressibility Scaling',GotIt)      

       IF(.NOT.GotIt) PseudoCompressibilityScale = 1.0
       IF(Transient) THEN
         PseudoCompressibilityScale = PseudoCompressibilityScale / dt
       END IF
       PseudoPressureUpdate = ListGetLogical( Model % Simulation, &
           'Pseudo Pressure Update',GotIt)
       IF (.NOT.GotIt) PseudoPressureUpdate = .FALSE.
     END IF

!------------------------------------------------------------------------------
! END CONFIG SOLVER

!______________________________________________________________________________
!******************************************************************************
!******************************************************************************
!******************************************************************************
     OPEN(unit=142,file='SauvegardeRayonSeries.dat',status='unknown',form='formatted')
     
     CRITGEOMFICH = 143
     BIFSTADIAGFICH=144
     MUMPSFICH = 145   
     ABETFICH = 146
     PADEFICH = 147
     RESIDUALFICH = 148
     BIFSTAINDICFICH = 149
     TIMEFILEFICH = 150
     DBSKPFICH    = 151
     
     OPEN(unit=CRITGEOMFICH,file='CritereBifurcartion.dat',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     WRITE(CRITGEOMFICH,900)           
     
     OPEN(unit=BIFSTADIAGFICH,file='BifSta_Mu0_Diag.dat',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     WRITE(BIFSTADIAGFICH,902)      
     
     OPEN(unit=MUMPSFICH,file='MUMPS_.log',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     
     OPEN(unit=ABETFICH,file='ABE_TANGENTS.log',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     write(ABETFICH,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
     
     OPEN(unit=PADEFICH,file='PADE.log',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     WRITE(PADEFICH,904)
     
     OPEN(unit=RESIDUALFICH,file='RESIDUAL.log',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     WRITE(RESIDUALFICH,906)     
     
     OPEN(unit=BIFSTAINDICFICH,file='BifstaIndicators.dat',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     WRITE(BIFSTAINDICFICH,908)          

     OPEN(unit=TIMEFILEFICH,file='Time.dat',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     WRITE(TIMEFILEFICH,*) '# EXCUTION TIMES FILE'
     
     
     OPEN(unit=DBSKPFICH,file='DBSKPlot.dat',status='unknown',form='formatted',ACCESS='SEQUENTIAL')     
     WRITE(DBSKPFICH,*) '#  K       Dk'
     
     
! ! !      CONFIG RESULT OUTPUT VTU FILES
!        FilePrefix = 'EXP3D-'
!        CALL ListAddString( Solver % Values,'Output File Name',FilePrefix)
!        CALL ListAddLogical( Solver % Values,'Ascii Output',.TRUE.)       
        CALL ListAddLogical( Solver % Values,'BinaryOutput',.TRUE.)              
!   SaveBulkOnly = GetLogical( Params,'Save Bulk Only',GotIt ) 
       CALL ListAddLogical( Solver % Values,'Save Bulk Only',.TRUE.)       
     

        
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               ! ADD STUFF IN VTU FILES ! TEST
! ! ! ! !            ResVec = 0._dp
! ! ! ! !            ResidualPointer => ResVec( 1: SIZE(ResVec) : NSDOFs )
! ! ! ! !            CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! !                   'Residual 1', 1,  ResidualPointer, FlowPerm) 
! ! ! ! ! 
! ! ! ! !            ResidualPointer => ResVec( 2: SIZE(ResVec) : NSDOFs )
! ! ! ! !            CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! !                   'Residual 2', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! ! 
! ! ! ! !            IF ( NSDOFs == 3 ) THEN
! ! ! ! !              ResidualPointer => ResVec( 3: SIZE(ResVec) : NSDOFs )
! ! ! ! !              CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! !                   'Residual Div', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! !            ELSE
! ! ! ! !              ResidualPointer => ResVec( 3: SIZE(ResVec) : NSDOFs )
! ! ! ! !              CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh,  Solver, &
! ! ! ! !                   'Residual 3', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! !              
! ! ! ! !              ResidualPointer => ResVec( 4: SIZE(ResVec) : NSDOFs )
! ! ! ! !              CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! !                   'Residual Div', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! !            ENDIF
! ! ! ! !            ResidualPointer => ResVec
! ! ! ! !            CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! !                 'Residual', NSDOFs, ResidualPointer ,FlowPerm )           
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!        FlowSolution_Init = FlowSolution
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = U0
!        FlowSol % Values => FlowSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!        FlowSolution = -1.d0
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = FlowSolution_Init
!        FlowSol % Values => FlowSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
     
      TOC = CPUTime() - TOC
      WRITE(TIMEFILEFICH,*) 'Step INIT - INIT',TOC
      CALL FLUSH(TIMEFILEFICH)

!           _   _ __  __    _____ _               _                       
!     /\   | \ | |  \/  |  / ____| |             | |                      
!    /  \  |  \| | \  / | | (___ | |_ ___ _ __   | |     ___   ___  _ __  
!   / /\ \ | . ` | |\/| |  \___ \| __/ _ \ '_ \  | |    / _ \ / _ \| '_ \ 
!  / ____ \| |\  | |  | |  ____) | ||  __/ |_) | | |___| (_) | (_) | |_) |
! /_/    \_\_| \_|_|  |_| |_____/ \__\___| .__/  |______\___/ \___/| .__/ 
!                                        | |                       | |    
!                                        |_|                       |_|    
     DO Step=1, ANMStep
       st = CPUTIme()
     
       CALL Info( 'ManFlowSolve', '-------------------------------------', Level=4 )
       WRITE( Message, * ) 'ANM STEP = ', Step!
       CALL Info( 'ManFlowSolve', ' ', Level=4 )
!------------------------------------------------------------------------------                  
      
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!        Pour conserver le profil générique et les CLs      
       IF( Step == 1 ) THEN
! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART       
! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART       
! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART       
!       RESTART :  U0  - Read from file
!                  L0  - Read from NODE4Lambda? or Check MAx on coutour ? or given by user?
!                  Ud  - FlowSolution_Init = Normalized FlowSolution on COUNTOUR by Lambda0 
!            => Correction?


!       INIT THEN FlowSolution_Init = FlowSolution( From Dirichlet of BC )
         FlowSolution_Init = FlowSolution
         VORIENT = 1._dp
         SENSPARCOURS = 1

         TimeVar => VariableGet( Solver % Mesh % Variables, 'Time')
         
         WRITE(6,*) "TimeVar % Values = ",TimeVar % Values
         lambda0 = TimeVar % Values(1) - 1.0_dp                     !!!!!!!! time has +1  DANGER ??? !!!!!!!!
         WRITE(6,*) "lambda0 = ",lambda0
         CALL FLUSH(6)      
! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART! RESTART RESTART         
       ENDIF
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!          U0 = FlowSolution
!          FlowSolution = FlowSolution_Init
!          FlowSol % Values => FlowSolution
!          CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!          FlowSolution = U0
!          FlowSol % Values => FlowSolution
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       !Remise des conditions limites générique (on espère ...)
! !        IF( Step /= 1 ) CALL CalCondLimCoeff( FlowSolution , 0.0_dp, 0, FlowSolution_Init )

       
       
       
!------------------------------------------------------------------------------         
       Element => GetActiveElement(1)
       NodeIndexes => Element % NodeIndexes

       n = GetElementNOFNodes()
       nb = GetElementNOFBDOFs()
       nd = GetElementDOFs( Indexes )
       
!------------------------------------------------------------------------------    
       IF (PseudoPressureExists .AND. PseudoPressureUpdate) &
         PseudoPressure = FlowSolution(NSDOFs:SIZE(FlowSolution):NSDOFs)

       at  = CPUTime()
       at0 = RealTime()
       at1 = RealTime()
       
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! CF MODULE DISCRETE OPERATORS
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Lt + Boundary condition
!------------------------------------------------------------------------------
       CALL Info( 'ManFlowSolve',Message, Level=4 )
       CALL Info( 'ManFlowSolve','-------------------------------------', Level=4 )
       CALL Info( 'ManFlowSolve', ' ', Level=4 )
       CALL Info( 'ManFlowSolve','Starting Assembly...', Level=4 )    
!------------------------------------------------------------------------------
!        WRITE(*,*) 'CALL OperatorsLtF'
!------------------------------------------------------------------------------  
       WRITE(6,*)  ' - - - - - - - - - - - - - - - - - - - - - - - -'
       NORMUO=DSQRT( DOT_PRODUCT( FlowSolution,FlowSolution ))    
       WRITE(6,*)  '|   ||U0||       = ',NORMUO
! TEST PICARD NEWTON DIFFERENCE FOR PRESSURE       
       CALL CalCondLimCoeff( FlowSolution ,   0._dp, 0, lambda0 * FlowSolution_Init )
!        CALL CalCondLimCoeff( FlowSolution , lambda0, 1,           FlowSolution_Init )

       TIC = CPUTime() 
       CALL OperatorsLtF( StiffMatrix, ForceVector, FlowSolution, NSDOFs,       &
                          MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
                          MeshSol, DensitySol, LocalNodes, dt, Transient,       &
                          PseudoCompressibilityScale, TempSol ,TempPerm,        &
                          Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &           
                          CompressibilityModel,DivDiscretization,               &
                          GradPDiscretization,NewtonLinearization,              &
                          Gravity                                               &
                        )
       FACTOTODO=.TRUE.
        
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- Lt_assembly',TIC
       CALL FLUSH(TIMEFILEFICH)

        
!TEST 20150317 - Remove Dirichlet from here, and put it outside to check residual without Diric on Stiff   
       TIC = CPUTime() 
       NORMTMP = DSQRT( DOT_PRODUCT( ForceVector,ForceVector ))
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- DSQRT(DOT_PRODUCT(',TIC       
       WRITE(6,*)  '|              ||Fav|| =', NORMTMP


!        WRITE(6,*)  'avant DSQRT( DOT_PRODUCT( ForceVector,ForceVector ))'
!        WRITE(6,*)  '|              ||Fap|| =',DSQRT( DOT_PRODUCT( ForceVector,ForceVector ))                          

!        WRITE(*,*)  '|      ||Lambda0 X F|| =',DSQRT( DOT_PRODUCT( lambda0*ForceVector,lambda0*ForceVector ))                                 
         
!        WRITE(*,*) 'OK OperatorsLtF'
                    
!------------------------------------------------------------------------------
!      R(U0,L0) =  L(U0) + Q(U0,U0) - L0 * F  ! F = 0
!------------------------------------------------------------------------------
       TIC = CPUTime() 
       CALL MatrixVectorMultiply(StiffMatrix,FlowSolution,ResVec)
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- MatrixVectorMultiply',TIC             
       NORMTMP = DSQRT( DOT_PRODUCT( ResVec,ResVec ))  
       WRITE(*,*)  '|RRB         ||Lt X U0)|| =',NORMTMP
       

       IF (Stabilize) THEN
              CALL  ANMQABSTAB( abeQAB , FlowSolution, FlowSolution, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                           
                         DensityTMP,Material,NodalMuTMP)
              NORMTMP = DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))  
              WRITE(*,*)  '|RRB  test||Qstab(U0,U0)|| =', NORMTMP
              ResVec = ResVec - abeQAB
       ELSE
              CALL  ANMQAB( abeQAA , FlowSolution, FlowSolution, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                           
                         DensityTMP,Material)
              NORMTMP = DSQRT( DOT_PRODUCT( abeQAA,abeQAA ))
              WRITE(*,*)  '|RRB         ||Q(U0,U0)|| =',NORMTMP
              ResVec = ResVec - abeQAA
       ENDIF
! ! ! ! !      R = [ L(U0) + Q(U0,U0) + Q(U0,U0) ] - Q(U0,U0) - L0 * F
! ! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
! ! ! !        ManualResidualU0j = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
! ! ! !        WRITE(*,*)  '|RRB ||L(U0) + Q(U0,U0)|| =',ManualResidualU0j
! ! ! ! !      

       TIC = CPUTime()
       CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )      
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- CalCondLimCoeff',TIC         
       
       ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
       WRITE(*,*)  '|RESIDUAL as ||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab
       WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,ManualResidualU0jstab
       
       RIKSCORRECTION=.FALSE.
       IF(ManualResidualU0jstab.GT.PrecRIKS) THEN
         IF (CallNewtonCorr.EQV..TRUE.) Then 
           NEWTCORRECTION=.TRUE.
         ELSE         
           RIKSCORRECTION=.TRUE.
         ENDIF
       ENDIF  
       
!      
       WRITE(RESIDUALFICH,*) ''
! !        ResVec = ResVec - lambda0*ForceVector
! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )
!        WRITE(*,*)  '|RR           ||ResVec|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))    
!        
! ! ! ! ! !        
! ! ! ! ! !        WRITE(*,*)  '| - - - - - COMPUTERESIDUALVECTOR BEFORE STEP - - - - - - -'
! ! ! ! ! ! 
! ! ! ! ! !        ResVec=0.0_dp
! ! ! ! ! !        CALL HMFlowResidual( U0, lambda0, ResVec, ResVecTMP, Model, Solver % Mesh , &
! ! ! ! ! !                             FlowPerm, NSDOFs, FNORMG, ResNorm )
! ! ! ! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )
! ! ! ! ! !        WRITE(*,*)  '|RRB      ||NODALFORCES|| =',FNORMG              
! ! ! ! ! ! ! DSQRT (Sum_elem [ s * (Element % hK**2 * SUM(Residual(1:dim)**2) + Residual(dim+1)**2 ) ]
! ! ! ! ! ! ! Element % hK : characterstic length of element
! ! ! ! ! ! ! Residual(dim+1) => continuity equation residual
! ! ! ! ! !        WRITE(*,*)  '|RRB     COMPUTED ResNorm =',ResNorm
! ! ! ! ! !        WRITE(*,*)  '|RRB   ||HMFlowResidual|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))
! ! ! ! ! !        CALL FLUSH(6)       
       
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   
!# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
! ! ! !        CALL BoundaryConditionMultCOEFF(Model, 1._dp , .TRUE.)
! ! ! ! !   FOR CORRECTION....
! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>

!  ___ ___ ___ ___ ___ _  _ _    ___ _____ 
! |   \_ _| _ \_ _/ __| || | |  | __|_   _|
! | |) | ||   /| | (__| __ | |__| _|  | |  
! |___/___|_|_\___\___|_||_|____|___| |_|  
!                                          
       CALL CalCondLimCoeff( FlowSolution , 0._dp, 0, FlowSolution_Init )
       TIC = CPUTime()       
       CALL DefaultDirichletBCs()       
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- DefaultDirichletBCs',TIC             
!        CHECK 
!        http://milamin.sourceforge.net/technical-notes/boundary-conditions-new
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\       
       WRITE(6,*)  'apres DefaultDirichletBCs -> OK'   
       WRITE(6,*)  'apres DSQRT( DOT_PRODUCT( ForceVector,ForceVector ))'
       CALL FLUSH(6)       
       
       
       TIC = CPUTime()

! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>       
! ______ _____ _   __ _____ 
! | ___ \_   _| | / //  ___|
! | |_/ / | | | |/ / \ `--. 
! |    /  | | |    \  `--. \
! | |\ \ _| |_| |\  \/\__/ /
! \_| \_|\___/\_| \_/\____/                                                     
!  _____ _________________ _____ _____ _____ _____ _____ _   _ 
! /  __ \  _  | ___ \ ___ \  ___/  __ \_   _|_   _|  _  | \ | |
! | /  \/ | | | |_/ / |_/ / |__ | /  \/ | |   | | | | | |  \| |
! | |   | | | |    /|    /|  __|| |     | |   | | | | | | . ` |
! | \__/\ \_/ / |\ \| |\ \| |___| \__/\ | |  _| |_\ \_/ / |\  |
!  \____/\___/\_| \_\_| \_\____/ \____/ \_/  \___/ \___/\_| \_/
        
       IF (RIKSCORRECTION) THEN
        RIC=0
        RICMAX=5
        DO WHILE((ManualResidualU0jstab.GE.PrecRIKS).AND.(RIC.LT.RICMAX))
         RIC=RIC+1
         WRITE(6,*) "----------------------------------------------------------------"
         WRITE(6,*) "-                RIKS CORRECTION ",RIC
         WRITE(6,*) "- PrecRIKS     = ",PrecRIKS
         WRITE(6,*) "- ResidualSTAB = ",ManualResidualU0jstab
         WRITE(6,*) "----------------------------------------------------------------"
!        PREDICTOR : KU1 = L1F
!        COMPUTE ORDER 1 MAN
         CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO )
         CALL HMumps_SolveLinearSystem( StiffMatrix, UMan(:,1), ForceVector, Solver , MUMPSFICH)
         CALL ANMExtractULambda( Lambda, UMan, 1, NoPressure,VMan)          
         FACTOTODO=.FALSE.
         CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO )

!        CORRECTOR K Dubar = -R0   Dubar <=> ResVecTMP pour ne pas declarer un tab en plus
         ResVecTMP=0.0_dp
         CALL CalCondLimCoeff( ResVec , 0._dp, 1, FlowSolution_Init )
         CALL HMumps_SolveLinearSystem( StiffMatrix, ResVecTMP, -ResVec, Solver , MUMPSFICH)
         CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )
         
         WRITE(6,*) ' - RIKS : ||DuBar||=',DSQRT(dot_product(ResVecTMP,ResVecTMP))
         
!        Dl = <Dubar,U1>/<U1,U1>
! VANNUCI98
!          RiksDL = -1.0_dp * DOT_PRODUCT(ResVecTMP,UMan(:,1))/DOT_PRODUCT(UMan(:,1),UMan(:,1))
         RiksDL = -Lambda(1) * DOT_PRODUCT(ResVecTMP,UMan(:,1))
         WRITE(6,*) ' - RIKS : RiksDL=',RiksDL
         WRITE(6,*) ' - RIKS : lambda0av=',lambda0         
         lambda0      = lambda0 + RiksDL
         WRITE(6,*) ' - RIKS : lambda0ap=',lambda0
         
!        DU = DL U1 + Dubar
! VANNUCI98
!          ResVecTMP = ResVecTMP + RiksDL * UMan(:,1)   
         ResVecTMP = ResVecTMP + RiksDL * UMan(:,1) /Lambda(1)
         WRITE(6,*) ' - RIKS : ||Du||=',DSQRT(dot_product(ResVecTMP,ResVecTMP))
         CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )   
         FlowSolution = FlowSolution + ResVecTMP
         CALL CalCondLimCoeff( FlowSolution ,   0._dp, 0, lambda0 * FlowSolution_Init )
         

! FlowSolution = FlowSolution + CORRECTION
! lambda0      = lambda0      + CORRECTION
! NEW RESIDUAL
!  -> new assembly
!        CALL CalCondLimCoeff( FlowSolution , lambda0, 1,           FlowSolution_Init )

         CALL OperatorsLtF( StiffMatrix, ForceVector, FlowSolution, NSDOFs,       &
                          MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
                          MeshSol, DensitySol, LocalNodes, dt, Transient,       &
                          PseudoCompressibilityScale, TempSol ,TempPerm,        &
                          Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &           
                          CompressibilityModel,DivDiscretization,               &
                          GradPDiscretization,NewtonLinearization,              &
                          Gravity                                               &
                        )
         FACTOTODO=.TRUE.
         ResVec = 0.0_dp
         CALL MatrixVectorMultiply(StiffMatrix,FlowSolution,ResVec)
         WRITE(6,*)  '|RIKS         ||Lt X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
       
         CALL  ANMQABSTAB( abeQAB , FlowSolution, FlowSolution, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                           
                         DensityTMP,Material,NodalMuTMP)
         WRITE(6,*)  '|RIKS   ||Qstab(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))     
         ResVec = ResVec - abeQAB
         CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
         ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
         WRITE(6,*)  '|RIKS RESSTAB||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab      
         WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,ManualResidualU0jstab
         
         CALL FLUSH(6)       
         CALL CalCondLimCoeff( FlowSolution , 0._dp, 0, FlowSolution_Init )
         CALL DefaultDirichletBCs()   
       END DO
       WRITE(*,*) "---------------------------------------------------------------------"
       ENDIF      
! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
       
       
       IF (NEWTCORRECTION) THEN
        NEWTITER=0
        DO WHILE((ManualResidualU0jstab.GE.PrecRIKS).AND.(NEWTITER.LT.NEWTMAX))
         NEWTITER=NEWTITER+1
         WRITE(6,*) "----------------------------------------------------------------"
         WRITE(6,*) "-                NEWTON CORRECTION ",NEWTITER
         WRITE(6,*) "- PrecRIKS     = ",PrecRIKS
         WRITE(6,*) "- ResidualSTAB = ",ManualResidualU0jstab
         WRITE(6,*) "----------------------------------------------------------------"

!        CORRECTOR K Dubar = -R0   Dubar <=> ResVecTMP pour ne pas declarer un tab en plus
         ResVecTMP=0.0_dp
         CALL CalCondLimCoeff( ResVec , 0._dp, 1, FlowSolution_Init )
         CALL HMumps_SolveLinearSystem( StiffMatrix, ResVecTMP, -ResVec, Solver , MUMPSFICH)
         CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )
         
         WRITE(6,*) ' - NEWTON : ||Du||=',DSQRT(dot_product(ResVecTMP,ResVecTMP))

         CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )   
         FlowSolution = FlowSolution + ResVecTMP
         CALL CalCondLimCoeff( FlowSolution ,   0._dp, 0, lambda0 * FlowSolution_Init )
         

! FlowSolution = FlowSolution + CORRECTION
! NEW RESIDUAL
!  -> new assembly
!        CALL CalCondLimCoeff( FlowSolution , lambda0, 1,           FlowSolution_Init )

         CALL OperatorsLtF( StiffMatrix, ForceVector, FlowSolution, NSDOFs,       &
                          MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
                          MeshSol, DensitySol, LocalNodes, dt, Transient,       &
                          PseudoCompressibilityScale, TempSol ,TempPerm,        &
                          Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &           
                          CompressibilityModel,DivDiscretization,               &
                          GradPDiscretization,NewtonLinearization,              &
                          Gravity                                               &
                        )
         FACTOTODO=.TRUE.
         ResVec = 0.0_dp
         CALL MatrixVectorMultiply(StiffMatrix,FlowSolution,ResVec)
         WRITE(6,*)  '|NEWTON         ||Lt X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
       
         CALL  ANMQABSTAB( abeQAB , FlowSolution, FlowSolution, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                           
                         DensityTMP,Material,NodalMuTMP)
         WRITE(6,*)  '|NEWTON   ||Qstab(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))     
         ResVec = ResVec - abeQAB
         CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
         ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
         WRITE(6,*)  '|NEWTON RESSTAB||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab      
         WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,ManualResidualU0jstab
         
         CALL FLUSH(6)       
         CALL CalCondLimCoeff( FlowSolution , 0._dp, 0, FlowSolution_Init )
         CALL DefaultDirichletBCs()   
       END DO
       WRITE(*,*) "---------------------------------------------------------------------"
       ENDIF      
! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- NLCorrection',TIC    
       
       TIC = CPUTime()
!# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #        
! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! .__            .___.__               __                
! ! ! ! ! ! |__| ____    __| _/|__| ____ _____ _/  |_  ___________ 
! ! ! ! ! ! |  |/    \  / __ | |  |/ ___\\__  \\   __\/  _ \_  __ \
! ! ! ! ! ! |  |   |  \/ /_/ | |  \  \___ / __ \|  | (  <_> )  | \/
! ! ! ! ! ! |__|___|  /\____ | |__|\___  >____  /__|  \____/|__|   
! ! ! ! ! !         \/      \/         \/     \/        
! ! ! ! ! !  ----  STEADY BIFURCATION INDICATOR 
       WRITE(6,*) '--------------------STEADY BIFURCATION INDICATOR---------------------------'
       CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO ) !!!
       WRITE(6,*) '--------------------Factorisation en cours---------------------------'
       call flush(6)
       CALL HMumps_SolveLinearSystem( StiffMatrix, INDYf, INDFalea , Solver , MUMPSFICH)
       FACTOTODO=.FALSE.       
       CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO ) !!!
       
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- First_Facto_Solve',TIC           
       TIC = CPUTime()
       
! TEST SANS PRESSION
       CALL PressureFree( INDYf, NoPressure, INDYf )
!      CONDITON "0" sur les contours   sauf inlet    
!        CALL CalCondLimCoeff( INDYf, 0.0_dp, 2, FlowSolution_Init )
! 0,1 CALCOND + CALIMP
       CALL CalCondLimCoeff( INDYf, 0.0_dp, 1, FlowSolution_Init )       
       
       WRITE(*,*) "IND - NORME INDYf    = ",DSQRT(DOT_PRODUCT(INDYf,INDYf))
!        CALL PressureFree( INDYf, NoPressure, INDYf )            
       IF (Step==1) THEN
!         COMPUTE Delta U_0 init
!         K  ( DUOI / mu_o_init)  = Falea  et  mu_o_init = 1
          INDMu0init = 1._dp
          INDDU0init = INDYf

          INDCADNum = DOT_PRODUCT( INDDU0init, INDDU0init)
          INNORMNum = INDCADNum
! Vannucci HDR
          INDVANNum = DOT_PRODUCT( INDFalea, INDDU0init)
          
          INDCADMu0 = INDMu0init
          INDVANMu0 = INDMu0init
          INNORMMu0 = INDMu0init

          INDMu0 = INDMu0init          
          INDDU0 = INDDU0init
       ELSE 
!           WRITE(*,*) "IND - DOT_PRODUCT( INDDU0init, INDDU0init)  C= ",DOT_PRODUCT( INDDU0init, INDDU0init)
!           WRITE(*,*) "IND - DOT_PRODUCT(      INDYf, INDDU0init)  C= ",DOT_PRODUCT( INDYf, INDDU0init)               
!           WRITE(*,*) "IND - DOT_PRODUCT( INDFalea, INDDU0init)    V= ",DOT_PRODUCT( INDFalea, INDDU0init)     
!           WRITE(*,*) "IND - DOT_PRODUCT(   INDFalea   , INDYf)    V= ",DOT_PRODUCT( INDFalea, INDYf)           

          INDCADDenum = DOT_PRODUCT(INDYf, INDDU0init)
          INDVANDenum = DOT_PRODUCT( INDFalea, INDYf)            
          INNORMDenum = DOT_PRODUCT(INDYf, INDYf)  ! Norm condition as in Hopf or Steady 2S 4S

          INDCADMu0 = INDCADNum / INDCADDenum
          INDVANMu0 = INDVANNum / INDVANDenum
          INNORMMu0 = DSQRT(INNORMNum / INNORMDenum)       
       ENDIF
       WRITE(BIFSTAINDICFICH,909) Step,Lambda0,INDCADMu0,INDVANMu0,INNORMMu0
       CALL FLUSH(BIFSTAINDICFICH)
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- INDICATOR',TIC           
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------                            


      IF ( NonlinearRelax /= 1.0d0 ) PSolution = FlowSolution

      
       TAC = CPUTime()

!------------------------ Implémentation de la MAN  ---------------------------
!           _   _ __  __    _____ ______ _____  _____ ______ 
!     /\   | \ | |  \/  |  / ____|  ____|  __ \|_   _|  ____|
!    /  \  |  \| | \  / | | (___ | |__  | |__) | | | | |__   
!   / /\ \ | . ` | |\/| |  \___ \|  __| |  _  /  | | |  __|  
!  / ____ \| |\  | |  | |  ____) | |____| | \ \ _| |_| |____ 
! /_/    \_\_| \_|_|  |_| |_____/|______|_|  \_\_____|______|
!                                                            
!------------------------------------------------------------------------------
!    _   _  _ __  __              _           _ 
!   /_\ | \| |  \/  |  ___ _ _ __| |___ _ _  / |
!  / _ \| .` | |\/| | / _ \ '_/ _` / -_) '_| | |
! /_/ \_\_|\_|_|  |_| \___/_| \__,_\___|_|   |_|
       WRITE(6,*) 'Actual ANM Step ', Step
       L2UMan = DSQRT( DOT_PRODUCT( FlowSolution,FlowSolution ))
       WRITE(6,*) 'Norme Point Solution ===', L2UMan
! Résolution à l'odre 1
       UMan     = 0.0_dp
       VMan     = 0.0_dp
       NormUMan = 0.0_dp
       NormFQ   = 0.0_dp
       UTemp    = 0.0_dp
       USAV     = 0.0_dp
       GradSAV  = 0.0_dp

!  ___  ___  _ __   _____ 
! / __|/ _ \| |\ \ / / __|
! \__ \ (_) | |_\ V /| _| 
! |___/\___/|____\_/ |___|
!      Kt(U0) U1 = L1 F
!      Kt(U0)  X =    F
       WRITE(6,*) "PF - NORME ForceVector = ",DSQRT(DOT_PRODUCT(ForceVector,ForceVector))
       CALL FLUSH(6)       

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        
!! SOLVE DIRECT METHOD ELMER        
!        CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .TRUE. )
!        CALL SolveLinearSystem( StiffMatrix, ForceVector, UMan(:,1), normMAN, NSDOFs, Solver )
!        CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', .FALSE. )
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        
! DIRECT MUMPS UNSYM SERIAL VERSION (threaded lib?)


! DANGER SI INDICATEUR ALORS DEJA FACTO? SINON IL FAUT FACTO!!!
!
!
! DANGER REMETTRE DES QUE POSSIBLE        
       TIC = CPUTime() 
       CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO )
!      A x = b, Solver       
       CALL HMumps_SolveLinearSystem( StiffMatrix, UMan(:,1), ForceVector, Solver , MUMPSFICH)
       FACTOTODO=.FALSE.
       CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO )       
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANM_SOLVE_U1',TIC
       CALL FLUSH(TIMEFILEFICH)
       
       NORMTMP = DSQRT(DOT_PRODUCT(UMan(:,1),UMan(:,1)))
       WRITE(*,*) "PF - NORME UMan(:,1) = ",NORMTMP
       
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        
!      X = U1 / L1 
! ! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$       
!        FlowSolution = UMan(:,1)
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = U0
!        FlowSol % Values => FlowSolution
! ! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
       CALL ANMExtractULambda( Lambda, UMan, 1, NoPressure,VMan)
       
     
!------------------------------------------------------------------------------
! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$       
!        FlowSolution = UMan(:,1)
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = U0
!        FlowSol % Values => FlowSolution
! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
!------------------------------------------------------------------------------
!    _   _  _ __  __              _           _               
!   /_\ | \| |  \/  |  ___ _ _ __| |___ _ _  | |___  ___ _ __ 
!  / _ \| .` | |\/| | / _ \ '_/ _` / -_) '_| | / _ \/ _ \ '_ \
! /_/ \_\_|\_|_|  |_| \___/_| \__,_\___|_|   |_\___/\___/ .__/
!                                                       |_|   
!      Résolution des ordres suivants => 2 à OMan
!      On prépare les termes croisés pour le 2nd membre MAN
!      Boucle sur les ordres > 2
       DO IO=2, OMan   
         TIC = CPUTime()
!                       __      __ 
!         /\ |\ ||\/|  |__)|__|(_  
!        /--\| \||  |  | \ |  |__)        
         FQMan = 0.0_dp
         CALL ANMPFrhsFQ(FQMan, UMan, IO, NSDOFs, FlowPerm, &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez , &
                         DensityTMP,Material,FlowSolution_Init)
         TIC = CPUTime() - TIC
         WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANM_ANMPFrhsFQ',TIC
         NORMTMP = DSQRT(DOT_PRODUCT(FQMan,FQMan))
         WRITE(*,*) "PF - Order ",IO," NORME FQMan = ",NORMTMP
!-----------------------------------------------------------------------------
!         ___  ___  _ __   _____ 
!        / __|/ _ \| |\ \ / / __|
!        \__ \ (_) | |_\ V /| _| 
!        |___/\___/|____\_/ |___|                 
!------------------------------------------------------------------------------
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !        
!! SOLVE DIRECT METHOD ELMER    
!          normMAN = 0.0_dp
!          CALL SolveLinearSystem( StiffMatrix, FQMan, UMan( : , IO ), normMAN, NSDOFs, Solver )
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
         TIC = CPUTime()
         CALL HMumps_SolveLinearSystem( StiffMatrix, UMan(:,IO), FQMan , Solver , MUMPSFICH)
         TIC = CPUTime() - TIC
         WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANM_SOLVE_UK',TIC
         CALL FLUSH(TIMEFILEFICH)
         
!------------------------------------------------------------------------------
!        Calcul de lambda à l'ordre IO et correction de Up
         CALL ANMExtractULambda(Lambda, UMan, IO, NoPressure, VMan )
         
!------------------------------------------------------------------------------
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$       
!        FlowSolution = UMan(:,IO)
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = U0
!        FlowSol % Values => FlowSolution
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

       END DO ! fin boucle sur les ordres
!------------------------------------------------------------------------------       
       write(*,*) 'fin boucle sur les ordres'
!  ___ _  _ ___      _   _  _ __  __    ___         _           _               
! | __| \| |   \    /_\ | \| |  \/  |  / _ \ _ _ __| |___ _ _  | |___  ___ _ __ 
! | _|| .` | |) |  / _ \| .` | |\/| | | (_) | '_/ _` / -_) '_| | / _ \/ _ \ '_ \
! |___|_|\_|___/  /_/ \_\_|\_|_|  |_|  \___/|_| \__,_\___|_|   |_\___/\___/ .__/
!                                                                         |_|   
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

         TAC = CPUTime() - TAC
         WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANM_END_LOOP',TAC


       close(80)
! ! ! ! ! ! ! ! ! !        WRITE(*,*) '------------------------------------------------------------------------------'
! ! ! ! ! ! ! ! ! !        WRITE(*,*) '-------------- TEST AFFICHAGE SERIE LAMBDA + VORIENT'
! ! ! ! ! ! ! ! ! !        IO=1
! ! ! ! ! ! ! ! ! !        write(*,*) IO,Lambda(IO), DOT_PRODUCT(VORIENT,UMan(:,IO))       
! ! ! ! ! ! ! ! ! !        DO IO=2, OMan   
! ! ! ! ! ! ! ! ! !         write(*,*) IO,Lambda(IO), DOT_PRODUCT(VORIENT,UMan(:,IO)),Lambda(IO-1)/Lambda(IO)
! ! ! ! ! ! ! ! ! !        END DO 
! ! ! ! ! ! ! ! ! !        WRITE(*,*) '------------------------------------------------------------------------------'
! ! ! ! ! ! ! ! ! !        
! ! ! ! ! ! ! ! ! ! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ! ! ! ! ! ! ! ! ! ! TEST MIZU
! ! ! ! ! ! ! ! ! !        FlowSolution = -30.0_dp
! ! ! ! ! ! ! ! ! !        FlowSol % Values => FlowSolution
! ! ! ! ! ! ! ! ! !        CALL ResultOutputSolver( Model,Solver,dt,Transient ) 
! ! ! ! ! ! ! ! ! !        
! ! ! ! ! ! ! ! ! !        FlowSolution = UMan(:,1)
! ! ! ! ! ! ! ! ! !        FlowSol % Values => FlowSolution
! ! ! ! ! ! ! ! ! !        CALL ResultOutputSolver( Model,Solver,dt,Transient )
! ! ! ! ! ! ! ! ! !       
! ! ! ! ! ! ! ! ! !        DO IO=2, OMan   
! ! ! ! ! ! ! ! ! !         FlowSolution = UMan(:,IO) / DSQRT(DOT_PRODUCT(UMan(:,IO) ,UMan(:,IO) ))
! ! ! ! ! ! ! ! ! !        
! ! ! ! ! ! ! ! ! !         FlowSol % Values => FlowSolution
! ! ! ! ! ! ! ! ! !         CALL ResultOutputSolver( Model,Solver,dt,Transient )        
! ! ! ! ! ! ! ! ! !        END DO        
! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! !     
! ! ! ! ! ! ! ! ! ! ! TEST MIZU
! ! ! ! ! ! ! ! ! !        FlowSolution = -20.0_dp
! ! ! ! ! ! ! ! ! !        FlowSol % Values => FlowSolution
! ! ! ! ! ! ! ! ! !        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
! ! ! ! ! ! ! ! ! !        
! ! ! ! ! ! ! ! ! !        FlowSolution = U0
! ! ! ! ! ! ! ! ! !        FlowSol % Values => FlowSolution       
! ! ! ! ! ! ! ! ! ! ! TEST MIZU       



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TEST POUR DOMB SYKES
       WRITE(*,*) " ----------------TEST POUR DOMB SYKES L=",lambda0,"----------------------"
      do idbs=2,OMan
        NORMtest2 = DOT_PRODUCT( UMan(:,idbs)   , UMan(:,idbs-1) )
        NORMtest1 = DOT_PRODUCT( UMan(:,idbs-1) , UMan(:,idbs-1))
        NORMtest2 = NORMtest2 / NORMtest1
        WRITE(*,*)  ' Ck/Ck-1(',idbs,')= ',NORMtest2
      end do
      CALL ANMDombSykesPlot( Lambda0, Lambda, UMan, OMan,DBSKPFICH )       

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       
!------------------------------------------------------------------------------       
! QUOI FAIRE UNE FOIS LA SERIE CALCULEE
!   + RAYON DE SERIE POLY
!   + CRITERE CM2013 : Alpha / ReyC / Prochaines branches
!   + PADE : amélioration rayon et critère bif sta sur pole (ReyC / Prochaines branches)
!   + Indicateur simplifié : premier ordre pour savoir si bif ou pas

!------------------------------------------------------------------------------       
!------------------------------------------------------------------------------       
!------------------------------------------------------------------------------     
!------------------------------------------------------------------------------       
!     // | |     /|    //| |     // | |  \\ / / 
!    //__| |    //|   // | |    //__| |   \  /  
!   / ___  |   // |  //  | |   / ___  |   / /   
!  //    | |  //  | //   | |  //    | |  / /\\  
! //     | | //   |//    | | //     | | / /  \\        
!       Calcul du rayon amax
       WRITE(*,*) '---------------------------------------------'
       TIC = CPUTime()
       amaxpolyUL = ANMseriesRoV( UMan, OMan, TolMan )
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANMseriesRoV',TIC       
       WRITE(*,*) 'AMAX = f(V P) = ',amaxpolyUL       
       amaxpolyVL = ANMseriesRoV( VMan, OMan, TolMan )     
       WRITE(*,*) 'AMAX = f(V 0) = ',amaxpolyVL
!------------------------------------------------------------------------------       
!------------------------------------------------------------------------------       
!------------------------------------------------------------------------------       
!   ____       _            _    _____            _   
!  / __ \     (_)          | |  / ____|          | |  
! | |  | |_ __ _  ___ _ __ | |_| |     ___  _ __ | |_ 
! | |  | | '__| |/ _ \ '_ \| __| |    / _ \| '_ \| __|
! | |__| | |  | |  __/ | | | |_| |___| (_) | | | | |_ 
!  \____/|_|  |_|\___|_| |_|\__|\_____\___/|_| |_|\__|
!        DirectionCont = DOT_PRODUCT(UMan(:,1),VORIENT)
!        DirectionCont = DirectionCont / DABS(DirectionCont)
! MODIF 20161130 - pilotage au premier pas, le prosca avec VORIENT=(1 1 1 ) peut être négatif.....
       IF( Step > 1 ) THEN
         DirectionCont = DOT_PRODUCT(UMan(:,1),VORIENT)
         DirectionCont = DirectionCont / DABS(DirectionCont)       
         SENSPARCOURS = SENSPARCOURS * DirectionCont
         amaxpolyUL = SENSPARCOURS*amaxpolyUL
       ENDIF
       

       
       IF (SENSPARCOURS>0) THEN
         WRITE(*,*) 'FORWARD  DIRECTION FOR THE CONTINUATION'
       ELSE 
         WRITE(*,*) 'BACKWARD DIRECTION FOR THE CONTINUATION'     
       ENDIF
       TIC = CPUTime() 
       CALL ANMDerivatePolySer(VORIENT,Lder,UMan,Lambda,OMan,amaxpolyUL)
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANMDerivatePolySer',TIC  
! -----------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------       
!------------------------------------------------------------------------------       
       WRITE(*,*) '---------------------------------------------'
       WRITE(*,*) '---- STEADY STATE BIFURCATION INDICATORs ----'
       WRITE(*,*) '---------------------------------------------'
!------------------------------------------------------------------------------
! .__            .___.__               __                
! |__| ____    __| _/|__| ____ _____ _/  |_  ___________ 
! |  |/    \  / __ | |  |/ ___\\__  \\   __\/  _ \_  __ \
! |  |   |  \/ /_/ | |  \  \___ / __ \|  | (  <_> )  | \/
! |__|___|  /\____ | |__|\___  >____  /__|  \____/|__|   
!         \/      \/         \/     \/        
!  ----  STEADY BIFURCATION INDICATOR  ---- AFTER DIRICHLET! !!!!<!><!><!><!><!><!><!><!><!><!><!><!>
! TODO:
!  -> INDICATOR : 
!     -- f alea
!     -- Fnl
!     -- DeltaU_o init avec mu_o=1
!     -- K X_k = F_k  => mu_k => Delta U_k
!     -- Rayon {Delta U}  VS Rayon {U}
! Question : magnitude force?, dirichlet Lt?, calimp calcond?
! Vannucci 98  : Mu0 =       <      f , DU0ini > / <     f , K-1 f  >
! Cadou 2006   : Mu0 =       < DU0ini , DU0ini > / < K-1 f , DU0ini >
! Cadou 2012   : Mu0 = SQRT( < DU0ini , DU0ini > / < K-1 f , K-1 f > )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! INDICATOR : init
!      SOLVE 
!      Yf = K^-1 f
! ! ! ! ! !        CALL SolveLinearSystem( StiffMatrix, INDFalea, INDYf, normMAN, NSDOFs, Solver )
! ! ! ! ! ! ! ! ! ! ! ! !        CALL HMumps_SolveLinearSystem( StiffMatrix, INDYf, INDFalea , Solver , MUMPSFICH)
! ! ! ! ! ! ! ! ! ! ! ! ! ! TEST SANS PRESSION
! ! ! ! ! ! ! ! ! ! ! ! !        CALL PressureFree( INDYf, NoPressure, INDYf )
! ! ! ! ! ! ! ! ! ! ! ! ! !      CONDITON "0" sur les contours   sauf inlet    
! ! ! ! ! ! ! ! ! ! ! ! ! !        CALL CalCondLimCoeff( INDYf, 0.0_dp, 2, FlowSolution_Init )
! ! ! ! ! ! ! ! ! ! ! ! ! ! 0,1 CALCOND + CALIMP
! ! ! ! ! ! ! ! ! ! ! ! !        CALL CalCondLimCoeff( INDYf, 0.0_dp, 1, FlowSolution_Init )       
! ! ! ! ! ! ! ! ! ! ! ! !        
! ! ! ! ! ! ! ! ! ! ! ! !        WRITE(*,*) "IND - NORME INDYf    = ",DSQRT(DOT_PRODUCT(INDYf,INDYf))
! ! ! ! ! ! ! ! ! ! ! ! ! !        CALL PressureFree( INDYf, NoPressure, INDYf )            
! ! ! ! ! ! ! ! ! ! ! ! !        IF (Step==1) THEN
! ! ! ! ! ! ! ! ! ! ! ! ! !         COMPUTE Delta U_0 init
! ! ! ! ! ! ! ! ! ! ! ! ! !         K  ( DUOI / mu_o_init)  = Falea  et  mu_o_init = 1
! ! ! ! ! ! ! ! ! ! ! ! !           INDMu0init = 1._dp
! ! ! ! ! ! ! ! ! ! ! ! !           INDDU0init = INDYf
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !           INDCADNum = DOT_PRODUCT( INDDU0init, INDDU0init)
! ! ! ! ! ! ! ! ! ! ! ! !           INNORMNum = INDCADNum
! ! ! ! ! ! ! ! ! ! ! ! ! ! Vannucci HDR
! ! ! ! ! ! ! ! ! ! ! ! !           INDVANNum = DOT_PRODUCT( INDFalea, INDDU0init)
! ! ! ! ! ! ! ! ! ! ! ! !           
! ! ! ! ! ! ! ! ! ! ! ! !           INDCADMu0 = INDMu0init
! ! ! ! ! ! ! ! ! ! ! ! !           INDVANMu0 = INDMu0init
! ! ! ! ! ! ! ! ! ! ! ! !           INNORMMu0 = INDMu0init
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !           INDMu0 = INDMu0init          
! ! ! ! ! ! ! ! ! ! ! ! !           INDDU0 = INDDU0init
! ! ! ! ! ! ! ! ! ! ! ! !        ELSE 
! ! ! ! ! ! ! ! ! ! ! ! ! !           WRITE(*,*) "IND - DOT_PRODUCT( INDDU0init, INDDU0init)  C= ",DOT_PRODUCT( INDDU0init, INDDU0init)
! ! ! ! ! ! ! ! ! ! ! ! ! !           WRITE(*,*) "IND - DOT_PRODUCT(      INDYf, INDDU0init)  C= ",DOT_PRODUCT( INDYf, INDDU0init)               
! ! ! ! ! ! ! ! ! ! ! ! ! !           WRITE(*,*) "IND - DOT_PRODUCT( INDFalea, INDDU0init)    V= ",DOT_PRODUCT( INDFalea, INDDU0init)     
! ! ! ! ! ! ! ! ! ! ! ! ! !           WRITE(*,*) "IND - DOT_PRODUCT(   INDFalea   , INDYf)    V= ",DOT_PRODUCT( INDFalea, INDYf)           
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !           INDCADDenum = DOT_PRODUCT(INDYf, INDDU0init)
! ! ! ! ! ! ! ! ! ! ! ! !           INDVANDenum = DOT_PRODUCT( INDFalea, INDYf)            
! ! ! ! ! ! ! ! ! ! ! ! !           INNORMDenum = DOT_PRODUCT(INDYf, INDYf)  ! Norm condition as in Hopf or Steady 2S 4S
! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! ! ! !           INDCADMu0 = INDCADNum / INDCADDenum
! ! ! ! ! ! ! ! ! ! ! ! !           INDVANMu0 = INDVANNum / INDVANDenum
! ! ! ! ! ! ! ! ! ! ! ! !           INNORMMu0 = DSQRT(INNORMNum / INNORMDenum)       
! ! ! ! ! ! ! ! ! ! ! ! !        ENDIF
! ! ! ! ! ! ! ! ! ! ! ! !        WRITE(BIFSTAINDICFICH,909) Step,Lambda0,INDCADMu0,INDVANMu0,INNORMMu0
! ! ! ! ! ! ! ! ! ! ! ! !        CALL FLUSH(BIFSTAINDICFICH)
       
!        WRITE(*,*) "IND - NORME INDDU0   = ",DSQRT(DOT_PRODUCT(INDDU0,INDDU0))
! ! ! ! ! ! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! ! ! ! ! ! ! ! INDICATOR : compute series Du_k mu_k
! ! ! ! ! ! !        DO IO=1, OMan    
! ! ! ! ! ! !          ! COMPOSE RHS
! ! ! ! ! ! !          INDFNL = 0.0_dp
! ! ! ! ! ! ! !        Fnl_1 = - [ Q( INDU_0 , U1) + Q( U1 , INDU_0 ) ]         
! ! ! ! ! ! ! !        Fnl_k = - SUM_{i=1}^k [ Q( INDU_k-i , Ui) + Q( Ui , INDU_k-i ) ]
! ! ! ! ! ! ! !          DO KK=1,NDL
! ! ! ! ! ! ! !           if (INDDU0(KK).GT.1e-15) THEN
! ! ! ! ! ! ! !              WRITE(*,*) "INDDU0(",KK,")=",INDDU0(KK)
! ! ! ! ! ! ! !           ENDIF
! ! ! ! ! ! ! !          ENDDO
! ! ! ! ! ! ! 
! ! ! ! ! ! !          CALL ANMSBIrhsFnl(INDFNL, INDDU0, INDDUserie, IO, NSDOFs, FlowPerm, &
! ! ! ! ! ! !                          USAV, GradSAV,                                      &
! ! ! ! ! ! !                          INDDUSAV, GradINDDUSAV,                             &
! ! ! ! ! ! !                          INDFNLtemp, INDUelex , INDUeley , INDUelez,         &
! ! ! ! ! ! !                          Density,Material,NDL,FlowSolution_Init)       
! ! ! ! ! ! !          WRITE(*,*) "IND - NORME Fnl_k    = ",DSQRT(DOT_PRODUCT(INDFNL,INDFNL))
! ! ! ! ! ! ! 
! ! ! ! ! ! !          CALL SolveLinearSystem( StiffMatrix, INDFNL, INDXk, normMAN, NSDOFs, Solver )
! ! ! ! ! ! !          WRITE(*,*) "IND - NORME INDXk    = ",DSQRT(DOT_PRODUCT(INDXk,INDXk))
! ! ! ! ! ! !          CALL PressureFree( INDXk, NoPressure, INDXk )         
! ! ! ! ! ! !          WRITE(*,*) "IND - NORME INDXk    = ",DSQRT(DOT_PRODUCT(INDXk,INDXk))
! ! ! ! ! ! !          
! ! ! ! ! ! ! !          STOP ' / / / / / / / / STOP DEBUG INDICATOR / / / / / / / / /'
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! !        Calcul de Mu_k et Delta Uk
! ! ! ! ! ! !          INDDENUM = DOT_PRODUCT( INDYf, INDDU0init )
! ! ! ! ! ! !          INDMuserie(IO)   = -1._dp * DOT_PRODUCT( INDXk, INDDU0init ) / INDDENUM
! ! ! ! ! ! !          INDDUserie(:,IO) = INDMuserie(IO) * INDYf + INDXk
! ! ! ! ! ! ! 
! ! ! ! ! ! !        END DO ! fin boucle sur les ordres
! ! ! ! ! ! ! ! ! ! ! !        amaxpolyINDDU = ANMseriesRoV( INDUserie, OMan, TolMan )     
! ! ! ! ! ! ! 
! ! ! ! ! ! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! ! ! ! ! ! ! ! INDICATOR : range of validity
! ! ! ! ! ! ! ! .__    .  .
! ! ! ! ! ! ! ! [__) _ \  /
! ! ! ! ! ! ! ! |  \(_) \/ 
! ! ! ! ! ! !        amaxpolyINDDU = ANMseriesRoV( INDDUserie, OMan, TolMan )     
! ! ! ! ! ! !        WRITE(*,*) 'amaxpoly INDICATEUR = ',amaxpolyINDDU
! ! ! ! ! ! !        WRITE(*,*) '---------------------------------------------'
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------


       
!------------------------------------------------------------------------------       
!------------------------------------------------------------------------------       
! ooooooooo.         .o.       oooooooooo.   oooooooooooo 
! `888   `Y88.      .888.      `888'   `Y8b  `888'     `8 
!  888   .d88'     .8"888.      888      888  888         
!  888ooo88P'     .8' `888.     888      888  888oooo8    
!  888           .88ooo8888.    888      888  888    "    
!  888          .8'     `888.   888     d88'  888       o 
! o888o        o88o     o8888o o888bood8P'   o888ooooood8 
!       
!  CONTINUATION PADE : COEFF + POLE + RoV + New solution
!  BIF DETECT PADE   : COEFF + POLE + VERIF BIF  Point limite, Indicateur 
!  TODO PADE V2 (EN TEST)
       write(*,*) "PADE PADE PADE PADE PADE PADE PADE PADE PADE PADE PADE PADE "

       VTMPOrtho = VMan
       amaxPad = amaxpolyUL
       PPole=0.0_dp
       
       TIC = CPUTime()
!        CALL ANMPadeBifDetect(amaxPad,U0,lambda0,UMan,VMan,Lambda,VTMPOrtho,VSOLTMP,OMan,  &
!                              PPole,LambdaCritPad,U0Pad,L0Pad)
       CALL ANMPadeBifDetectZRHQR(amaxPad,U0,lambda0,UMan,VMan,Lambda,VTMPOrtho,VSOLTMP,OMan,  &
                             PPole,LambdaCritPad,U0Pad,L0Pad,Step,amaxpolyUL)
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANMPadeBifDetectZRHQR',TIC                               
! PPole : plus petit pole
! ! ! ! IS_IT_A_LIMIT_POINT( LAMBDA, Acrit ,  NORDRE ) RESULT ( LIMITPOINT )                             
       write(*,*) "PADE done PADE done PADE done PADE done PADE done PADE done "
!------------------------------------------------------------------------------       
!------------------------------------------------------------------------------       
       
       
       
       
!------------------------------------------------------------------------------
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         
!------------------------------------------------------------------------------       
       alphaBIF=0.0_dp
       IF (BIFSTAGEOMDETECT) THEN
!    ___           _                     _       _ ____   ___  _ ____  
!   / __\___   ___| |__   /\/\   ___  __| | __ _| |___ \ / _ \/ |___ \ 
!  / /  / _ \ / __| '_ \ /    \ / _ \/ _` |/ _` | | __) | | | | | __) |
! / /__| (_) | (__| | | / /\/\ \  __/ (_| | (_| | |/ __/| |_| | |/ __/ 
! \____/\___/ \___|_| |_\/    \/\___|\__,_|\__,_|_|_____|\___/|_|_____| 
!        EPSILON1=0.001_dp
!        EPSILON2=0.000001_dp
!        EPSILON1=0.000001_dp
!        EPSILON2=0.001_dp
         TIC = CPUTime()
         CALL BIFdetectCritCochMed( VMan, OMan, NoPressure, calcalphabif, sumcrit ,  &
                                           EPSILON1, EPSILON2, alphaVBIF, Lambda )
         TIC = CPUTime() - TIC
         WRITE(TIMEFILEFICH,*) 'Step ', Step, '- BIFdetectCritCochMed',TIC                                               
         alphaBIF = alphaVBIF   
         write(*,*) "*-*-*-   alphaBIF  = ",alphaBIF 
!********************************************************************************
! BIFURCATION DETECTION CM       
       ENDIF
!------------------------------------------------------------------------------
! SORTIE critère bifurcation DANS FICHIER
! AMAX 
! ALPHA
       write(*,*) "*-*-*-   SENSPARCOURS*alphaBIF  = ",SENSPARCOURS*alphaBIF 
       alphaBIFSENSPARCOURS = SENSPARCOURS*alphaBIF 
       
       
       IF ( alphaBIFSENSPARCOURS.EQ.0.0_dp ) THEN
        write(*,*) "NO CRITICAL POINT DETECTED"
        ! CONT CLASSIQUE PF 94
         amax = amaxpolyUL
         write(*,*) ' PF 94 : ANMSolPol amax=',amax
         IF (amax<0) THEN
           WRITE(*,*) " PF 94 : MOVING BACKWARD"
         ENDIF
         
         U0 = FlowSolution
         TIC = CPUTime()
         CALL ANMSolPol( U0, lambda0, UMan, Lambda, OMan, amax )               
         TIC = CPUTime() - TIC
         WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANMSolPol',TIC   
         
         
         
         
       ELSE IF ( alphaBIFSENSPARCOURS < 0 ) THEN
         ! CRITICAL POINT DETECTED BEHIND
         WRITE(*,*) ""
         WRITE(*,*) ""       
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         write(146,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         write(*,*) "*-*-*-      CRITICAL POINT DETECTED BEHIND"
         write(146,*) "*-*-*-      CRITICAL POINT DETECTED BEHIND"
         write(*,*) "*-*-*-   Geometric progression detected (Cochelin Medale 2013)"
         write(*,*) "*-*-*-   alphaBIFSENSPARCOURS  = ",alphaBIFSENSPARCOURS 
         write(146,*) "*-*-*-   alphaBIFSENSPARCOURS  = ",alphaBIFSENSPARCOURS 
         CALL FLUSH(146)  
       
!------------------------------------------------------------------------------       

! ___________.____       _____  __      __.____     ___________ _________ _________
! \_   _____/|    |     /  _  \/  \    /  \    |    \_   _____//   _____//   _____/
!  |    __)  |    |    /  /_\  \   \/\/   /    |     |    __)_ \_____  \ \_____  \ 
!  |     \   |    |___/    |    \        /|    |___  |        \/        \/        \
!  \___  /   |_______ \____|__  /\__/\  / |_______ \/_______  /_______  /_______  /
!      \/            \/       \/      \/          \/        \/        \/        \/ 
! - - - - - - : SERIE PROPRE ET MODE DE BIFURCATION       
! test it with pressure UMan UManProp  , BifUPMode
         TIC = CPUTime() 
         CALL BIFCMserieprop( UMan, Lambda, OMan, alphaBIF , UManProp , LambdaUProp, UBifMode )
         TIC = CPUTime() - TIC
         WRITE(TIMEFILEFICH,*) 'Step ', Step, '- BIFCMserieprop',TIC          
         CALL CalCondLimCoeff( UBifMode, 0.0_dp , 1 , FlowSolution_Init)  
         ! 0 : put a given vector for both Constraints and Load
         ! 1 : multiply Both contraints and load by a scalar
         ! 2 : Only Constraints (CALIMP EVE)  
         
         
         
         
         !!!!! 20170215 : Should check if another ProgGeom is seen..... or not?
         !!!!! SE3D_E3_Ai7 padé see at 308 Rec 314 at 311 Rec 314 .... is this a bif?
       
!     // | |     /|    //| |     // | |  \\ / / 
!    //__| |    //|   // | |    //__| |   \  /  
!   / ___  |   // |  //  | |   / ___  |   / /   
!  //    | |  //  | //   | |  //    | |  / /\\  
! //     | | //   |//    | | //     | | / /  \\        
!       Calcul du rayon amax
       WRITE(*,*) '---------------------------------------------'
       WRITE(*,*) 'AMAX AVANT = ',amaxpolyUL   
       amaxpolyUL = ANMseriesRoV( UManProp, OMan-1, TolMan )     
       WRITE(*,*) 'AMAX UManProp = f(V P) = ',amaxpolyUL   
         amax = amaxpolyUL

         
!   _________      .__          __  .__               
!  /   _____/ ____ |  |  __ ___/  |_|__| ____   ____  
!  \_____  \ /  _ \|  | |  |  \   __\  |/  _ \ /    \ 
!  /        (  <_> )  |_|  |  /|  | |  (  <_> )   |  \
! /_______  /\____/|____/____/ |__| |__|\____/|___|  /
!         \/                                       \/ 

         U0 = FlowSolution
         CALL ANMSolPol( U0, lambda0, UMan, Lambda, OMan, amax )
         
!------------------------------------------------------------------------------       
         
         
         
         
       

       ELSE IF ( alphaBIFSENSPARCOURS > 0 ) THEN       
         ! CRITICAL POINT DETECTED AHEAD
         ! LP  : LIMIT POINT       - dimKer=1 - Phi - Uder=?Phi  and Lder=0?
         ! SBP : Simple BifPoint   - dimKer=1 - Phi - ABE W - PitchFork?
         ! MBP : multiple BifPoint - dimKer>1 MumpsV5
         WRITE(*,*) ""
         WRITE(*,*) ""       
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         write(146,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
         write(*,*) "*-*-*-      CRITICAL POINT DETECTED AHEAD"
         write(146,*) "*-*-*-      CRITICAL POINT DETECTED AHEAD"
         write(*,*) "*-*-*-   Geometric progression detected (Cochelin Medale 2013)"
         write(*,*) "*-*-*-   alphaBIFSENSPARCOURS  = ",alphaBIFSENSPARCOURS 
         write(146,*) "*-*-*-   alphaBIFSENSPARCOURS  = ",alphaBIFSENSPARCOURS 
         CALL FLUSH(146)       
         BFN = BFN + 1
! ! ! ! ! ! ! ! 
! ___________.____       _____  __      __.____     ___________ _________ _________
! \_   _____/|    |     /  _  \/  \    /  \    |    \_   _____//   _____//   _____/
!  |    __)  |    |    /  /_\  \   \/\/   /    |     |    __)_ \_____  \ \_____  \ 
!  |     \   |    |___/    |    \        /|    |___  |        \/        \/        \
!  \___  /   |_______ \____|__  /\__/\  / |_______ \/_______  /_______  /_______  /
!      \/            \/       \/      \/          \/        \/        \/        \/ 
! - - - - - - : SERIE PROPRE ET MODE DE BIFURCATION
         CALL BIFCMserieprop( VMan, Lambda, OMan, alphaBIF , VManProp , LambdaVProp, VBifMode )
         CALL CalCondLimCoeff( VBifMode, 0.0_dp , 1 , FlowSolution_Init)       
         
! test it with pressure UMan UManProp  , BifUPMode
         TIC = CPUTime()
         CALL BIFCMserieprop( UMan, Lambda, OMan, alphaBIF , UManProp , LambdaUProp, UBifMode )
         TIC = CPUTime() - TIC
         WRITE(TIMEFILEFICH,*) 'Step ', Step, '- BIFCMserieprop',TIC             
         CALL CalCondLimCoeff( UBifMode, 0.0_dp , 1 , FlowSolution_Init)  
         ! 0 : put a given vector for both Constraints and Load
         ! 1 : multiply Both contraints and load by a scalar
         ! 2 : Only Constraints (CALIMP EVE)        
         
!###############################################################BifMode     Phi             
!          BifMode  = VBifMode         ! Delta Uc <=> Only V no Pressure inside (why??)
         BifMode  = UBifMode         ! Delta Uc <=> Only V no Pressure inside (why??)
         XPDD=DOT_PRODUCT(BifMode,BifMode)
         XNBifM=DSQRT(XPDD)
         WRITE(*,*) "VERIF - || BifMode || = ",XNBifM
         WRITE(146,*) "VERIF - || BifMode || = ",XNBifM
         CALL FLUSH(146)  
!###############################################################BifMode

!###############################################################DERIVATION AT CRITICAL POINT          
! ________               .__               __  .__               
! \______ \   ___________|__|__  _______ _/  |_|__| ____   ____  
!  |    |  \_/ __ \_  __ \  \  \/ /\__  \\   __\  |/  _ \ /    \ 
!  |    `   \  ___/|  | \/  |\   /  / __ \|  | |  (  <_> )   |  \
! /_______  /\___  >__|  |__| \_/  (____  /__| |__|\____/|___|  /
!         \/     \/                     \/                    \/     
! Flawless series derivated at critical point
! Uder is either : W for Tangent at BifPoint,      Lder != 0
!                : Phi for tangent at Limit Point, Lder  = 0
         CALL ANMDerivatePolySer(Uder,Lder,UManProp,LambdaUProp,OMan-1,alphaBIFSENSPARCOURS)
         ! Der = W or Der= Phi for Mizushima cas??
         
         
!###############################################################SOLUTION AT CRITICAL POINT 
!   _________      .__          __  .__               
!  /   _____/ ____ |  |  __ ___/  |_|__| ____   ____  
!  \_____  \ /  _ \|  | |  |  \   __\  |/  _ \ /    \ 
!  /        (  <_> )  |_|  |  /|  | |  (  <_> )   |  \
! /_______  /\____/|____/____/ |__| |__|\____/|___|  /
!         \/                                       \/ 
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"           
         write(*,*) " - -  COMPUTE NEW SOLUTION AT CRITICAL POINT " 
         write(*,*) "*-*-*-   ANMSolPol UManProp alphaBIF"
! - - - - - - : COMPUTE NEW SOLUTION AT CRITICAL POINT
!  USING PRESSURE
         UBif  = U0
         LUBif = lambda0
         CALL ANMSolPol( UBif, LUBif, UManProp, LambdaUProp, OMan-1, alphaBIF )
         CALL CalCondLimCoeff( UBif, LUBif , 1 , FlowSolution_Init)  
          
! ! ! !  USING NO PRESSURE         
! ! !          VBif = U0
! ! !          LVBif = lambda0
! ! ! !!!!!!! TEST 2015 03 30 : Ubif ou UPBif with pressure?         
! ! !          CALL ANMSolPol( VBif, LVBif, VManProp, LambdaVProp, OMan-1, alphaBIF )
! ! !          CALL CalCondLimCoeff( VBif, LVBif , 1 , FlowSolution_Init)       
! ! !          
!!!!!!! TEST 2015 03 30 : Ubif ou UPBif with pressure?         
         write(*,*) "*-*-*-   LUBif = ",LUBif         
         write(146,*) "*-*-*-   LUBif = ",LUBif 
!          write(*,*) "*-*-*-   LVBif = ",LVBif         
!          write(146,*) "*-*-*-   LVBif = ",LVBif          
         CALL FLUSH(146) 
         write(*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"            
!###############################################################SOLUTION AT CRITICAL POINT                 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       TIC = CPUTime()
       FlowSolution = -52.0_dp
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )  
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ResultOutputSolver',TIC       
       FlowSolution = UBif
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )
       
!        FlowSolution = -53.0_dp
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )        
       FlowSolution = BifMode
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )         
!           WRITE(*,*) "VERIF - || UBif || : ",DSQRT(DOT_PRODUCT(UBif,UBif))
!        FlowSolution = -54.0_dp
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )        
       FlowSolution = Uder
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )      
       
              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
         

!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!                            
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!                            
!   _________.__                     .__                
!  /   _____/|__| ____    ____  __ __|  | _____ _______ 
!  \_____  \ |  |/    \  / ___\|  |  \  | \__  \\_  __ \
!  /        \|  |   |  \/ /_/  >  |  /  |__/ __ \|  | \/
! /_______  /|__|___|  /\___  /|____/|____(____  /__|   
!         \/         \//_____/                 \/       
! ___________                                    __   
! \__    ___/____    ____    ____   ____   _____/  |_ 
!   |    |  \__  \  /    \  / ___\_/ __ \ /    \   __\
!   |    |   / __ \|   |  \/ /_/  >  ___/|   |  \  |  
!   |____|  (____  /___|  /\___  / \___  >___|  /__|  
!                \/     \//_____/      \/     \/      
! ________                              __                
! \_____  \ ______   ________________ _/  |_  ___________ 
!  /   |   \\____ \_/ __ \_  __ \__  \\   __\/  _ \_  __ \
! /    |    \  |_> >  ___/|  | \// __ \|  | (  <_> )  | \/
! \_______  /   __/ \___  >__|  (____  /__|  \____/|__|   
!         \/|__|        \/           \/                   

! -------- FREE MUMPS         
         WRITE(6,*) "Branches : HMumps_Free_mumpsIDL"
         CALL HMumps_Free_mumpsIDL(StiffMatrix)         
! ! ! !          IF(ASSOCIATED(StiffMatrix)) THEN
! ! ! !           write(6,*) "COUCOU CALL FreeMatrix(StiffMatrix)"
! ! ! !           CALL FLUSH(6)          
! ! ! !           CALL FreeMatrix(StiffMatrix)
! ! ! !           
! ! ! !           write(6,*) "StiffMatrix => NULL()"          
! ! ! !           CALL FLUSH(6)          
! ! ! !           StiffMatrix => NULL()
! ! ! !          END IF            
! - - - - - - - - - - - - DEALLOCATE - - - - - - - - - - - - - -  
         WRITE(6,*)  '- Ltc  : Assembly'
         CALL FLUSH(6)
         CALL OperatorsLtF( StiffMatrix, ForceVector, UBif, NSDOFs,               &
                            MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
                            MeshSol, DensitySol, LocalNodes, dt, Transient,       &
                            PseudoCompressibilityScale, TempSol ,TempPerm,        &
                            Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &           
                            CompressibilityModel,DivDiscretization,               &
                            GradPDiscretization,NewtonLinearization,              &
                            Gravity                                               &                           
                            )  
! SINGULAR TANGENT OPERATOR : CHECK RANK WITH MUMPS V5 AND GET NULL BASIS
         FACTOTODO=.TRUE.
         CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO )       
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! Test 20150213                            
          CALL CalCondLimCoeff( ForceVector , 0.0_dp, 2, FlowSolution_Init)  
! 0 : put a given vector for both Constraints and Load
! 1 : multiply Both contraints and load by a scalar
! 2 : Only Constraints (CALIMP EVE)          
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
         CALL CalCondLimCoeff( UBif, 0._dp, 0, FlowSolution_Init )      
         

         
!          
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>      
! ! AJOUT 20161017 - Corrrection de Lt c pour voir si ca permet de meilleurs tangentes
! !------------------------------------------------------------------------------
! !      R(U0,L0) =  L(U0) + Q(U0,U0) - L0 * F  ! F = 0
! !------------------------------------------------------------------------------
!        ResVec = 0.0_dp
!        CALL MatrixVectorMultiply(StiffMatrix,UBif,ResVec)
!        WRITE(*,*)  '|RRB         ||Lt X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))  
!        
! 
!        IF (Stabilize) THEN
!               CALL  ANMQABSTAB( abeQAB , UBif, UBif, NSDOFs, FlowPerm, &
!                          FQManTemp, &
!                          Uelex , Ueley , Uelez, &
!                          Velex , Veley , Velez, &                           
!                          DensityTMP,Material,NodalMuTMP)
!               WRITE(*,*)  '|RRB  test||Qstab(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))   
!               ResVec = ResVec - abeQAB
!        ELSE
!               CALL  ANMQAB( abeQAA , UBif, UBif, NSDOFs, FlowPerm, &
!                          FQManTemp, &
!                          Uelex , Ueley , Uelez, &
!                          Velex , Veley , Velez, &                           
!                          DensityTMP,Material)
!               WRITE(*,*)  '|RRB         ||Q(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAA,abeQAA ))         
!               ResVec = ResVec - abeQAA
!        ENDIF
! ! ! ! ! !      R = [ L(U0) + Q(U0,U0) + Q(U0,U0) ] - Q(U0,U0) - L0 * F
! ! ! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
! ! ! ! !        ManualResidualU0j = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
! ! ! ! !        WRITE(*,*)  '|RRB ||L(U0) + Q(U0,U0)|| =',ManualResidualU0j
! ! ! ! ! !      
!        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
!        ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
!        WRITE(*,*)  '|RESIDUAL as ||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab
!        WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,ManualResidualU0jstab
!        
!        RIKSCORRECTION=.FALSE.
!        IF(ManualResidualU0jstab.GT.PrecRIKS) THEN
!          IF (CallNewtonCorr.EQV..TRUE.) Then 
!            NEWTCORRECTION=.TRUE.
!          ELSE         
!            RIKSCORRECTION=.TRUE.
!          ENDIF
!        ENDIF  
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
!        IF (NEWTCORRECTION) THEN
!         NEWTITER=0
!         DO WHILE((ManualResidualU0jstab.GE.PrecRIKS).AND.(NEWTITER.LT.NEWTMAX))
!          NEWTITER=NEWTITER+1
!          WRITE(6,*) "----------------------------------------------------------------"
!          WRITE(6,*) "-                NEWTON CORRECTION ",NEWTITER
!          WRITE(6,*) "- PrecRIKS     = ",PrecRIKS
!          WRITE(6,*) "- ResidualSTAB = ",ManualResidualU0jstab
!          WRITE(6,*) "----------------------------------------------------------------"
! 
! !        CORRECTOR K Dubar = -R0   Dubar <=> ResVecTMP pour ne pas declarer un tab en plus
!          ResVecTMP=0.0_dp
!          CALL CalCondLimCoeff( ResVec , 0._dp, 1, FlowSolution_Init )
!          CALL HMumps_SolveLinearSystem( StiffMatrix, ResVecTMP, -ResVec, Solver , MUMPSFICH)
!          CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )
!          
!          WRITE(6,*) ' - NEWTON : ||Du||=',DSQRT(dot_product(ResVecTMP,ResVecTMP))
! 
!          CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )   
!          UBif = UBif + ResVecTMP
!          CALL CalCondLimCoeff( UBif ,   0._dp, 0, lambda0 * FlowSolution_Init )
!          
! 
! ! FlowSolution = FlowSolution + CORRECTION
!          CALL OperatorsLtF( StiffMatrix, ForceVector, UBif, NSDOFs,       &
!                           MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
!                           MeshSol, DensitySol, LocalNodes, dt, Transient,       &
!                           PseudoCompressibilityScale, TempSol ,TempPerm,        &
!                           Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &           
!                           CompressibilityModel,DivDiscretization,               &
!                           GradPDiscretization,NewtonLinearization,              &
!                           Gravity                                               &
!                         )
!          FACTOTODO=.TRUE.
!          ResVec = 0.0_dp
!          CALL MatrixVectorMultiply(StiffMatrix,UBif,ResVec)
!          WRITE(6,*)  '|NEWTON         ||Lt X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
!        
!          CALL  ANMQABSTAB( abeQAB , UBif, UBif, NSDOFs, FlowPerm, &
!                          FQManTemp, &
!                          Uelex , Ueley , Uelez, &
!                          Velex , Veley , Velez, &                           
!                          DensityTMP,Material,NodalMuTMP)
!          WRITE(6,*)  '|NEWTON   ||Qstab(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))     
!          ResVec = ResVec - abeQAB
!          CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
!          ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
!          WRITE(6,*)  '|NEWTON RESSTAB||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab      
!          WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,ManualResidualU0jstab
!          
!          CALL FLUSH(6)       
!          CALL CalCondLimCoeff( UBif , 0._dp, 0, FlowSolution_Init )
!          CALL DefaultDirichletBCs()   
!        END DO
!        WRITE(*,*) "---------------------------------------------------------------------"
!        
!        
! ! TODO ->  ON PEUT FAIRE LE TEST DE LA DIMENSION DU NOYAU ICI, CAR LtC déjà Facto!!
!        
!        
!        
! !LtC post Newton  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - -
! 
!     ! -------- FREE MUMPS         
!          WRITE(6,*) "Branches : HMumps_Free_mumpsIDL"
!          CALL HMumps_Free_mumpsIDL(StiffMatrix)         
! ! ! ! !          IF(ASSOCIATED(StiffMatrix)) THEN
! ! ! ! !           write(6,*) "COUCOU CALL FreeMatrix(StiffMatrix)"
! ! ! ! !           CALL FLUSH(6)          
! ! ! ! !           CALL FreeMatrix(StiffMatrix)
! ! ! ! !           
! ! ! ! !           write(6,*) "StiffMatrix => NULL()"          
! ! ! ! !           CALL FLUSH(6)          
! ! ! ! !           StiffMatrix => NULL()
! ! ! ! !          END IF            
! ! - - - - - - - - - - - - DEALLOCATE - - - - - - - - - - - - - -  
!          WRITE(6,*)  '- Ltc  : Assembly'
!          CALL FLUSH(6)
!          CALL OperatorsLtF( StiffMatrix, ForceVector, UBif, NSDOFs,               &
!                             MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
!                             MeshSol, DensitySol, LocalNodes, dt, Transient,       &
!                             PseudoCompressibilityScale, TempSol ,TempPerm,        &
!                             Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &           
!                             CompressibilityModel,DivDiscretization,               &
!                             GradPDiscretization,NewtonLinearization,              &
!                             Gravity                                               &                           
!                             )  
! ! SINGULAR TANGENT OPERATOR : CHECK RANK WITH MUMPS V5 AND GET NULL BASIS
!          FACTOTODO=.TRUE.
!          CALL ListAddLogical( Solver % Values, 'Linear System Refactorize', FACTOTODO )       
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! ! Test 20150213                            
!           CALL CalCondLimCoeff( ForceVector , 0.0_dp, 2, FlowSolution_Init)  
! ! 0 : put a given vector for both Constraints and Load
! ! 1 : multiply Both contraints and load by a scalar
! ! 2 : Only Constraints (CALIMP EVE)          
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
!          CALL CalCondLimCoeff( UBif, 0._dp, 0, FlowSolution_Init )      
! !LtC post Newton  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - -
! !LtC post Newton  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - - -  - - - -
! 
!           
!        
!        
!        ENDIF      
! !FIN  AJOUT 20161017 - Corrrection de Lt c pour voir si ca permet de meilleurs tangentes  
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
!          
!          
!          
         
         
         
         
         
         
         

!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!                            
! VERIFICATION
         CALL MatrixVectorMultiply(StiffMatrix,BifMode,Vtmp)
         LtB = DSQRT( DOT_PRODUCT( Vtmp,Vtmp ))
         WRITE(*,*)  'VERIF -||Ltc X BifMode||=',LtB                       
         WRITE(146,*)  'VERIF -||Ltc X BifMode||=',LtB                             
         
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!    
         WRITE(6,*)  '- Ltc  : DefaultDirichletBCs'
         CALL FLUSH(6)
         CALL DefaultDirichletBCs()       
         WRITE(6,*)  '- Ltc  : DefaultDirichletBCs : DONE'
         CALL FLUSH(6)
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!                            
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!                            


         


         
!###############################################################LIMIT POINT? 
! .__  .__        .__  __                  .__        __   
! |  | |__| _____ |__|/  |_  ______   ____ |__| _____/  |_ 
! |  | |  |/     \|  \   __\ \____ \ /  _ \|  |/    \   __\
! |  |_|  |  Y Y  \  ||  |   |  |_> >  <_> )  |   |  \  |  
! |____/__|__|_|  /__||__|   |   __/ \____/|__|___|  /__|  
!               \/           |__|                  \/   
! LIMIT POINT?       
!   UDER VS BifMode
!   Lder = 0?
!   CHECK abeCoC = <abePsiL,  Q(BifMode,BifMode)> !=0  for a  Limit point

         WRITE(6,*) " - - - - - - - - - - - - - - - - - - - - - - -"
         WRITE(6,*) " - Test POINT LIMITE : Lder_crit=",Lder
         WRITE(6,*) " - Test POINT LIMITE : UBifMode=Uder?"
      
         XPDD=DOT_PRODUCT(Uder,Uder)
         XNBifM=DSQRT(XPDD)
         Vtmp = BifMode - Uder*DOT_PRODUCT(BifMode,Uder)/XPDD
         XPDD=DOT_PRODUCT(Vtmp,Vtmp)
         XNBifM=DSQRT(XPDD)        
         WRITE(6,*) " - Test POINT LIMITE : BifMode - Proj Uder BifMode =",XNBifM         
         
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  OUTPUT   Uder     
!        FlowSolution = -56
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        FlowSolution = Uder
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )         
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   


! VERIF p32 MEI 2000
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !          
! ! ! CHECK abeCoC = <abePsiL,  Q(BifMode,BifMode)> !=0  for a  Limit point
! !           CALL ANMQAB( abeQBB , BifMode,BifMode, NSDOFs, FlowPerm, &
! !                          FQManTemp, &
! !                          Uelex , Ueley , Uelez, &
! !                          Velex , Veley , Velez, &                              
! !                          DensityTMP,Material)
! !           abeCoC = DOT_PRODUCT(abePsiL,abeQBB)
! !           WRITE(*,*) "VERIF - C = <LM,Q(BifMode,BifMode)> =",abeCoC                         
! !           WRITE(146,*) "VERIF - C = <LM,Q(BifMode,BifMode)> =",abeCoC
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !          
          
          
          
         IF (DABS(Lder).LE.1E-6) THEN
           write(6,*) "Point limite : 1 TANGENTE! FAIRE ATTENTION"
            ! Tangent = Mode
            ! One PostBif branch
!          CALL SENSCO(T(IT(18)),T(IT(42)),IP,APOSITIF,RAYON,NPOI,DA,NDL)
           WRITE(6,*) " - - - - - - - - - - - - - - - - - - - - - - -"         
           CALL FLUSH(6)         
           STOP 'LIMIT POINT NET YET IMPLEMENTED'
           
!###############################################################LIMIT POINT

         ELSE         
         
!###############################################################BIFURCATION POINT
! __________.__  _____                            __  .__               
! \______   \__|/ ____\_ _________   ____ _____ _/  |_|__| ____   ____  
!  |    |  _/  \   __\  |  \_  __ \_/ ___\\__  \\   __\  |/  _ \ /    \ 
!  |    |   \  ||  | |  |  /|  | \/\  \___ / __ \|  | |  (  <_> )   |  \
!  |______  /__||__| |____/ |__|    \___  >____  /__| |__|\____/|___|  /
!         \/                            \/     \/                    \/ 
!               .__        __   
! ______   ____ |__| _____/  |_ 
! \____ \ /  _ \|  |/    \   __\
! |  |_> >  <_> )  |   |  \  |  
! |   __/ \____/|__|___|  /__|  
! |__|                  \/   
     
        

!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!       
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*! 
!  _   _ _    _ __  __ ____  ______ _____       ____  ______ 
! | \ | | |  | |  \/  |  _ \|  ____|  __ \     / __ \|  ____|
! |  \| | |  | | \  / | |_) | |__  | |__) |   | |  | | |__   
! | . ` | |  | | |\/| |  _ <|  __| |  _  /    | |  | |  __|  
! | |\  | |__| | |  | | |_) | |____| | \ \    | |__| | |     
! |_| \_|\____/|_|  |_|____/|______|_|  \_\    \____/|_|     
!  _______       _   _  _____ ______ _   _ _______ _____ 
! |__   __|/\   | \ | |/ ____|  ____| \ | |__   __/ ____|
!    | |  /  \  |  \| | |  __| |__  |  \| |  | | | (___  
!    | | / /\ \ | . ` | | |_ |  __| | . ` |  | |  \___ \ 
!    | |/ ____ \| |\  | |__| | |____| |\  |  | |  ____) |
!    |_/_/    \_\_| \_|\_____|______|_| \_|  |_| |_____/ 
! CHECK IT WITH MUMPS V5!!
        
! NB TANGENTS: + 2 if Bif Simple
!              + 3 case of SqSq mulitple easy or 5 degen?


! - - - - - - : BIRUCATED BRANCHES : case of simple bifurcation
!  -   -   -  : A) TANGENTs == ABE with W, LeftMode, BifMode
!      -      : A-1 Assemble Lt at Critical solution  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! Test 20150213
    

       NBTANGENTS = 2
       If ( .NOT.Allocated(abeUT)) THEN
         ALLOCATE( abeUT( NDL , NBTANGENTS ),  abeLT( NBTANGENTS )  )
       ENDIF
       abeUT = 0.0_dp
       abeLT = 0.0_dp
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!       
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!       



        TIC = CPUTime()
!///////////////////////////////////////////////////////////////////////////////////////
!         Solve the Collection-system with MUMPS
          StiffAugW => StiffMatrix % CollectionMatrix
          IF(.NOT.ASSOCIATED(StiffAugW)) THEN
!             write(*,*) 'COUCOU AllocateMatrix'
            StiffAugW => AllocateMatrix()
            StiffAugW % FORMAT = MATRIX_LIST
          ELSE
            DEALLOCATE(StiffAugW % RHS)
            StiffAugW % Values = 0.0_dp
          END IF
! How many left modes?  Get it with Mumps? Delta Uc mutliple? 
! How many ABE? 2 with two leftmodes, or 1 with two bif modes?



!      COMPUTE Psi and W using one factorization : Help of MumpsV5 ICNTL(9)
       CALL BIFStaDimKer1ComputePsiW( StiffMatrix, ForceVector, BifMode,  &
                                 normMAN, NSDOFs, Solver,MUMPSFICH,           &
                                 StiffAugW, abePsiL, abeW )
       TIC = CPUTime() - TIC
       WRITE(TIMEFILEFICH,*) 'Step ', Step, '- BIFStaDimKer1ComputePsiW',TIC                                     
!///////////////////////////////////////////////////////////////////////////////////////


!---------------------------------------------------------------------------------------
!  _           __ _   __  __           _      
! | |         / _| | |  \/  |         | |     
! | |     ___| |_| |_| \  / | ___   __| | ___ 
! | |    / _ \  _| __| |\/| |/ _ \ / _` |/ _ \
! | |___|  __/ | | |_| |  | | (_) | (_| |  __/
! |______\___|_|  \__|_|  |_|\___/ \__,_|\___|
!      -      : A-3 LeftMode as Lt(Transpose==TRUE) PSI = 0 and <Ut1orth,PSI>=1      
!---------------------------------------------------------------------------------------


          NLeftMode = DSQRT(dot_product(abePsiL,abePsiL))
!           WRITE(*,*) "VERIF - <LeftMode,LeftMode> =",dot_product(abePsiL,abePsiL)          
!           WRITE(146,*) "VERIF - <LeftMode,LeftMode> =",dot_product(abePsiL,abePsiL)          
          WRITE(*,*) "VERIF - || LeftMode || = ",NLeftMode                                          
          WRITE(146,*) "VERIF - || LeftMode || = ",NLeftMode       
          
          PLB = DOT_PRODUCT(abePsiL,BifMode)
          WRITE(*,*) "VERIF - <LeftMode,BifMode> = 1? : ",PLB
          WRITE(146,*) "VERIF - <LeftMode,BifMode> = 1? : ",PLB

         WRITE(6,*)  '- LeftMode  : ResultOutputSolver'
         CALL FLUSH(6)            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
       FlowSolution = -57.0_dp
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient ) 
       FlowSolution = abePsiL
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
         WRITE(6,*)  '- LeftMode  : ResultOutputSolver : DONE'
         CALL FLUSH(6)            

!---------------------------------------------------------------------------------------
! __          __
! \ \        / /
!  \ \  /\  / / 
!   \ \/  \/ /  
!    \  /\  /   
!     \/  \/                                
!      -      : A-2 W as Lt W = F and <Ut1orth,W>=0
! Solve with constraint as lagrangian multipliers
! CollectionMatrix [Stiff, ModeBif]
! CollectionVector {F,0}
        CALL CalCondLimCoeff( abeW , 0.0_dp, 1, FlowSolution_Init)                                     
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          NTMP1 = dot_product(Uder,abeW)
          NTMP1 = NTMP1/dot_product(Uder,Uder)
          WRITE(*,*) "VERIF - <Uder,W>/<Uder,Uder> =",NTMP1
          ResVecTMP =NTMP1*Uder-abeW 
          NTMP1 = DSQRT( DOT_PRODUCT(ResVecTMP,ResVecTMP))
          WRITE(*,*) "VERIF - ||<Uder,W>/<Uder,Uder> * Uder - W||=",NTMP1
          
        
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!    
          LMW=dot_product(abePsiL,abeW)
          WRITE(*,*) "VERIF - <LeftMode,W> =",LMW
          WRITE(146,*) "VERIF - <LeftMode,W> =",LMW
! Attention ici Lt est avec Dirichlet! ca peut jouer sur resultat... test avec W et calcond?          
        CALL MatrixVectorMultiply(StiffMatrix,abeW,Vtmp)
        NTMP1 = DSQRT( DOT_PRODUCT( Vtmp,Vtmp )) 
        WRITE(*,*)  'VERIF -||Ltc X abeW||=',  NTMP1
        WRITE(146,*)  'VERIF -||Ltc X abeW||=',NTMP1
        ResVecTMP=0.0_dp
        ResVecTMP = Vtmp - ForceVector
        NTMP1 = DSQRT( DOT_PRODUCT( ResVecTMP,ResVecTMP ))  
        WRITE(*,*)  '||Ltc W - Faug||=',NTMP1  
        WRITE(146,*)  '||Ltc W - Faug||=',NTMP1       

        BMW = DOT_PRODUCT(abeW,BifMode)
        WRITE(*,*) "VERIF - <W,BifMode> = 0? : ",BMW
        WRITE(146,*) "VERIF - <W,BifMode> = 0? : ",BMW
        XPWW=dot_product(abeW,abeW)
        WRITE(*,*) "VERIF - <W,W> =",XPWW
        WRITE(146,*) "VERIF - <W,W> =",XPWW
        NW = DSQRT(XPWW)
        WRITE(*,*) "VERIF - || W || = ",NW
        WRITE(146,*) "VERIF - || W || = ",NW
        CALL FLUSH(146)
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
       FlowSolution = abeW
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!!*!*!*!*!*!          

       FlowSolution = -1.1_dp
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!---------------------------------------------------------------------------------------




!---------------------------------------------------------------------------------------
!           ____  ______ 
!     /\   |  _ \|  ____|
!    /  \  | |_) | |__   
!   / /\ \ |  _ <|  __|  
!  / ____ \| |_) | |____ 
! /_/    \_\____/|______|
!      -      : A-3 ABE: second degree " TPsiL *Q(Ut,Ut) = 0 "
!      -      : A-3-a : lambda1 +/-

! 
!  BIF BRIS SYM alors W et Phi sont les tangentes!!! on verfie que A et C sont nuls!!!
! 
! 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! abeCoA = abePsiLT X Q(abeW,abeW)
          TIC = CPUTime()
          CALL ANMQAB( abeQAA , abeW,abeW, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                         
                         DensityTMP,Material)
          TIC = CPUTime() - TIC
          WRITE(TIMEFILEFICH,*) 'Step ', Step, '- ANMQAB',TIC                          
          NTMP2 = DSQRT(dot_product(abeQAA,abeQAA)) 
          WRITE(*,*) "VERIF - || Q(W,W) || = ",NTMP2
          WRITE(146,*) "VERIF - || Q(W,W) || = ",NTMP2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!        FlowSolution = abeQAA
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          

          abeCoA = DOT_PRODUCT(abePsiL,abeQAA)
          WRITE(*,*) "VERIF - A = <LM,Q(W,W)> = =",abeCoA
          WRITE(146,*) "VERIF - A = <LM,Q(W,W)> = =",abeCoA

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !          
! abeCoB = abePsiLT X Q(abeW,BifMode) + abePsiL X Q(BifMode,abeW)
          CALL ANMQAB( abeQAB , abeW, BifMode, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                            
                         DensityTMP,Material)
          CALL ANMQAB( abeQBA , BifMode, abeW, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                               
                         DensityTMP,Material)
          abeCoB = DOT_PRODUCT(abePsiL,abeQAB)
          abeCoB = abeCoB + DOT_PRODUCT(abePsiL,abeQBA)
          WRITE(*,*) "VERIF - abeCoB =",abeCoB                         
          WRITE(146,*) "VERIF - abeCoB =",abeCoB                         
         
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !          
! abeCoC = abePsiLT X Q(BifMode,BifMode)
          CALL ANMQAB( abeQBB , BifMode,BifMode, NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)
          abeCoC = DOT_PRODUCT(abePsiL,abeQBB)
          WRITE(*,*) "VERIF - C = <LM,Q(BifMode,BifMode)> =",abeCoC                         
          WRITE(146,*) "VERIF - C = <LM,Q(BifMode,BifMode)> =",abeCoC                         
                         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!        FlowSolution = abeQBB
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                          
          NTMP2= DSQRT(dot_product(abeQBB,abeQBB))
          WRITE(*,*) "VERIF - || Q(BifMode,BifMode) || = ", NTMP2        
          WRITE(146,*) "VERIF - || Q(BifMode,BifMode) || = ",NTMP2        

! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -            
!  _____ _____ _____ _____ _____ _____ _____ _____ 
! |_   _|  _  |   | |   __|   __|   | |_   _|   __|
!   | | |     | | | |  |  |   __| | | | | | |__   |
!   |_| |__|__|_|___|_____|_____|_|___| |_| |_____|
! TEST PITCHFORK!!!!!          
! BIF SIMPLE
! SI <LM,Q(W,W)>=0 et        <LM,Q(BifMode,BifMode)> = 0 et  abeCoB!=0 alors PitchFork! 
! Uts  = L1 W
! Utas = Phi
! CF B. WERNER AND A. SPENCE 1984 "THE COMPUTATION OF SYMMETRY-BREAKING BIFURCATION POINTS"
! + Seydel récent
          ZPF=1e-8
          IF( (DABS(abeCoC).LE.ZPF).AND.(DABS(abeCoA).LE.ZPF).AND.(DABS(abeCoB).GT.ZPF))THEN
! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -            
! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -                      
! ! ! ! !  _____ _____ _______ _____ _    _ ______ ____  _____  _  __
! ! ! ! ! |  __ \_   _|__   __/ ____| |  | |  ____/ __ \|  __ \| |/ /
! ! ! ! ! | |__) || |    | | | |    | |__| | |__ | |  | | |__) | ' / 
! ! ! ! ! |  ___/ | |    | | | |    |  __  |  __|| |  | |  _  /|  <  
! ! ! ! ! | |    _| |_   | | | |____| |  | | |   | |__| | | \ \| . \ 
! ! ! ! ! |_|   |_____|  |_|  \_____|_|  |_|_|    \____/|_|  \_\_|\_\          
            WRITE(6,*) " - PITCHFORK BIFURCATION DETECTED - "
            WRITE(6,*) " - PITCHFORK <LM,Q(W,W)>                =",abeCoA
            WRITE(6,*) " - PITCHFORK <LM,Q(BM,BM)>              =",abeCoC                      
            WRITE(6,*) " - PITCHFORK <LM, Q(W,BM) + <LM,Q(BM,W) =",abeCoB
            PITCHFORKbifDETECTED=.TRUE.
           ! SYM : L1^2 <W,W> + L1^2 = 1
           abeLT(1) = 1._dp / DSQRT(XPWW +1)
           abeUT(:,1) = abeLT(1) * abeW 
           WRITE(*,*) "Utsym =",abeLT(1),"W+"
           
           ! ASYM : ||BM||=1
           abeLT(2) = 0.0_dp
           abeUT(:,2) =  BifMode
           WRITE(*,*) "Utasym = BifMode"
! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -                       
! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -            
          ELSE  ! sinon faire classique

          WRITE(*,*) " - - - - - - - - - - - - - - - -"                                   
          WRITE(146,*) " - - - - - - - - - - - - - - - -"                                   
          Discr = abeCoB*abeCoB - 4.0_dp * abeCoA * abeCoC
          WRITE(*,*) "VERIF - Discr =",Discr
          WRITE(146,*) "VERIF - Discr =",Discr
          
          DeuxA = 2.0_dp * abeCoA 
          WRITE(*,*) "VERIF - DeuxA =",DeuxA
          WRITE(146,*) "VERIF - DeuxA =",DeuxA
          
          XTMP = DSQRT(Discr)
          WRITE(*,*) "VERIF - DSQRT(Discr) =",XTMP
          WRITE(146,*) "VERIF - DSQRT(Discr) =",XTMP
          
          abeL1temp = ( -1.0_dp * abeCoB + XTMP ) / DeuxA
          abeL2temp = ( -1.0_dp * abeCoB - XTMP ) / DeuxA
!         L1 / E1
          WRITE(*,*) "VERIF - abeL1temp =",abeL1temp
          WRITE(146,*) "VERIF - abeL1temp =",abeL1temp
          WRITE(*,*) "VERIF - abeL2temp =",abeL2temp   
          WRITE(146,*) "VERIF - abeL2temp =",abeL2temp   
          WRITE(*,*) " - - - - - - - - - - - - - - - -"          
          WRITE(146,*) " - - - - - - - - - - - - - - - -"          
!      -      : A-3-b : Eta1
          XTMP = abeL1temp * abeL1temp
          ! VERSION PILOTAGE U,L EVE : 1/DSQRT(<W,> S1**2 + <DUc,DUc> + S1**2)  <u,u> + l^2  = 1
          abeE1 = 1.0_dp / DSQRT( XTMP * XPWW + XPDD + XTMP )
          abeLT(1) = abeL1temp * abeE1
          WRITE(*,*) "VERIF - abeL1 =",abeLT(1)
          WRITE(146,*) "VERIF - abeL1 =",abeLT(1)         
          WRITE(*,*) "VERIF - abeE1 eve=",abeE1
          WRITE(146,*) "VERIF - abeE1 eve=",abeE1
          
          XTMP = abeL2temp*abeL2temp
          abeE2 = 1.0_dp / DSQRT( XTMP*XPWW + XPDD + XTMP )          
          abeLT(2) = abeL2temp * abeE2          
          WRITE(*,*) "VERIF - abeL2 =",abeLT(2)
          WRITE(146,*) "VERIF - abeL2 =",abeLT(2)
          WRITE(*,*) "VERIF - abeE2 eve=",abeE2          
          WRITE(146,*) "VERIF - abeE2 eve=",abeE2          
!      -      : A-4 : Tangents

          WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - -"          
          WRITE(*,*) " - - - - - - - - TANGENTES - - - - - - - - - -"
          WRITE(146,*) " - - - - - - - - TANGENTES - - - - - - - - - -"
!         L1 E1          
          WRITE(*,*) "Ut1 =",abeLT(1),"W+",abeE1,"BifMode"
          WRITE(146,*) "Ut1 =",abeLT(1),"W+",abeE1,"BifMode"
          abeUT(:,1) = abeLT(1) * abeW + abeE1 * BifMode
          NTMP1=DSQRT(DOT_PRODUCT(abeUT(:,1),abeUT(:,1)))
          WRITE(*,*) " VERIF - || Ut1 ||= ",NTMP1
          WRITE(146,*) " VERIF - || Ut1 ||= ",NTMP1
          
          CALL ANMQAB( abeQBB , abeUT(:,1),abeUT(:,1), NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)
          NTMP1 = dot_product(abePsiL,abeQBB)   
          WRITE(*,*) 'Verif ABE < PSI,Q(UT,UT)>=',NTMP1 
          WRITE(146,*) 'Verif ABE < PSI,Q(UT,UT)>=',NTMP1  
                                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!        FlowSolution = abeT1
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = U0
!        FlowSol % Values => FlowSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         L2 E2
          WRITE(*,*) "Ut2 =",abeLT(2),"W+",abeE2,"BifMode"
          WRITE(146,*) "Ut2 =",abeLT(2),"W+",abeE2,"BifMode"

          abeUT(:,2) = abeLT(2) * abeW + abeE2 * BifMode
          NTMP1 = DSQRT(DOT_PRODUCT(abeUT(:,2),abeUT(:,2)))
          WRITE(*,*) " VERIF - || Ut2 ||= ",NTMP1
          WRITE(146,*) " VERIF - || Ut2 ||= ",NTMP1
          
          CALL ANMQAB( abeQBB , abeUT(:,2),abeUT(:,2), NSDOFs, FlowPerm, &
                         FQManTemp, &
                         Uelex , Ueley , Uelez, &
                         Velex , Veley , Velez, &                              
                         DensityTMP,Material)
          NTMP1=dot_product(abePsiL,abeQBB)  
          WRITE(*,*) 'Verif ABE < PSI,Q(UT,UT)>=',NTMP1         
          WRITE(146,*) 'Verif ABE < PSI,Q(UT,UT)>=',NTMP1       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!        FlowSolution = abeT2
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = U0
!        FlowSol % Values => FlowSolution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
          NTMP1=dot_product(abeUT(:,1),abeUT(:,2)) 
          WRITE(*,*) 'Verif TANGENTs <Ut1,Ut2>', NTMP1         
          WRITE(146,*) 'Verif TANGENTs <Ut1,Ut2>',NTMP1          
          WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - -"          
          WRITE(146,*) " - - - - - - - - - - - - - - - - - - - - - - -"    
          CALL FLUSH(146)
! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -          
         ENDIF ! END CASE PITCHFORK OR NOT

!---------------------------------------------------------------------------------------
! ______                      _               
! | ___ \                    | |              
! | |_/ /_ __ __ _ _ __   ___| |__   ___  ___ 
! | ___ \ '__/ _` | '_ \ / __| '_ \ / _ \/ __|
! | |_/ / | | (_| | | | | (__| | | |  __/\__ \
! \____/|_|  \__,_|_| |_|\___|_| |_|\___||___/
!
!  -   -   -  : B) Branches simple : 2 tangentes
! IF BIF SYM a=0, b!=0 et c=0 alors W et Phi sont les tangentes.
! Simplification du calcul.
!
!  abeTa : BranchePos a , BrancheNeg -a
!  abeTb : BranchePos a , BrancheNeg -a
!
! DO FOR TANGENTS: POS NEG
         WRITE(*,*) ""
         WRITE(*,*) ""         
         WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - -"          
         DO IT=1,NBTANGENTS
           WRITE(*,*) "|-|-|-| Branch for tangent ",IT
! DANGER DANGER ASTUCE POUR 'SANS RESTART' 3-IT=2, puis 1 pour avoir FB2 en sortie!!!!
! DANGER DANGER ASTUCE POUR 'SANS RESTART' 3-IT=2, puis 1 pour avoir FB2 en sortie!!!!
!            CALL  BIFs_Compute_serie_one_branch( UBranch, LBranch, abeUT(:,3-IT), abeLT(3-IT),   &
           IF (PITCHFORKbifDETECTED) THEN
! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -                       
             IF(IT==1) THEN
               SYMBRANCH=.TRUE.
             ELSE
               SYMBRANCH=.FALSE.
             ENDIF
             TIC = CPUTime()
             CALL  BIFs_Compute_serie_one_branch_SYM( UBranch, LBranch, abeUT(:,IT), abeLT(IT),   &
                                                  BifMode, abeW ,abePsiL, NDL, NSDOFs, OMan,  &
                                                  FQMan, FQManTemp, USAV, GradSAV,            &
                                                  Uelex, Ueley, Uelez, Velex, Veley, Velez ,  &
                                                  DensityTMP, Material,                       &
                                                  Solver, StiffMatrix,FlowSolution_Init,      &
                                                  MUMPSFICH, Vtmp,StiffAugW,                  &
                                                  abeQAA, abeQBB, FlowPerm,SYMBRANCH ,abeCoB)                             
! * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * -                                                              
             TIC = CPUTime() - TIC
             WRITE(TIMEFILEFICH,*) 'Step ', Step, '- BIFs_Compute_serie_one_branch_SYM',TIC
           ELSE
             TIC = CPUTime()
             CALL  BIFs_Compute_serie_one_branch( UBranch, LBranch, abeUT(:,IT), abeLT(IT),   &
                                                  BifMode, abeW ,abePsiL, NDL, NSDOFs, OMan,  &
                                                  FQMan, FQManTemp, USAV, GradSAV,            &
                                                  Uelex, Ueley, Uelez, Velex, Veley, Velez ,  &
                                                  DensityTMP, Material,                       &
                                                  Solver, StiffMatrix,FlowSolution_Init,      &
                                                  MUMPSFICH, Vtmp,StiffAugW,                  &
                                                  abeQAA, abeQBB, FlowPerm ) 
             TIC = CPUTime() - TIC
             WRITE(TIMEFILEFICH,*) 'Step ', Step, '- BIFs_Compute_serie_one_branch',TIC
           ENDIF
! DANGER DANGER ASTUCE POUR 'SANS RESTART' 3-IT=2, puis 1 pour avoir FB2 en sortie!!!!
! DANGER DANGER ASTUCE POUR 'SANS RESTART' 3-IT=2, puis 1 pour avoir FB2 en sortie!!!!
              
                   
                   
! RoV Serie AUG RoV Serie AUG RoV Serie AUG RoV Serie AUG RoV Serie AUG                   
!  -> SOL PAD, OU SOL POL sur Serie Ubranch avec a et -a
           amaxpolyBr = ANMseriesRoV( UBranch, OMan, TolMan )     
           WRITE(*,*) "|-|-|-| UBranch Range of Validity amaxpolyBr=",amaxpolyBr

! POS (U(a),L(a))
           WRITE(*,*) "|-|-|-| POSITIVE Branch for tangent ",IT
!            UTemp = UBif
           UTemp = UBif
           Ltemp = LUBif
! - - - - - - - - - - - - - - - - - - OUT - - - - - - - - - - - - - - - - - - - - - - - 
           CDBx = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 1) )
           CDBy = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 2) )
           CDBz = 0.0_dp
           if (dim == 3)  CDBz = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 3) )       
           WRITE(144,903)  Step,Ltemp,CDBx,CDBy,CDBz   
           CALL FLUSH(144)


! -- -- Branch POS -- -- Branch POS -- -- Branch POS -- -- Branch POS -- -- Branch POS
! TEST SOL POL
!             CALL ANMSolPol( UTemp, Ltemp, UBranch, LBranch, OMan, amaxpolyBr )
!             CALL CalCondLimCoeff( UTemp, Ltemp , 1 , FlowSolution_Init)
      
! TEST SOL POL MULIPLE NPOI       
           CALL ANMMultipleSolPol( UTemp, Ltemp, UBranch, LBranch, OMan, amaxpolyBr,      &
                                  NBSolByStep, NODEOUT, Step, FlowPerm, NSDOFs, dim)
           CALL CalCondLimCoeff( UTemp, Ltemp , 1 , FlowSolution_Init)
                                  
! - - - - - - - - - - - - - - - - - - OUT - - - - - - - - - - - - - - - - - - - - - - - 
           CDBx = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 1) )
           CDBy = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 2) )
           CDBz = 0.0_dp
           if (dim == 3)  CDBz = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 3) )       
           WRITE(144,903)  Step,Ltemp,CDBx,CDBy,CDBz       
           CALL FLUSH(144)
       
! - - - - - - - - - - - - - - - - - - FILE OUT - - - - - - - - - - - - - - - - - - - - -
! ! ! !            Binary  = .FALSE. ! FIRST, after BINARY TRUE IS BETTER FOR PRECISION
           SaveAll = .TRUE. 
           Solver % Mesh % SavesDone = 0  !  Trick to  enforce new file!!
           FlowSolution = UTemp
           FlowSol % Values => FlowSolution            
           SaveCount = SaveResult( 'BIF'//trim(str(BFN))//'BR'//trim(str(IT))//'P.restart', &
                                    Solver % Mesh,Step,Ltemp,Binary,SaveAll,FreeSurfaceFlag ) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
! RESTART SYM FORWARD
          IF ((SYMRESTART.AND.SYMBRANCH).AND.(Ltemp>LUBif)) THEN
              URestart = UTemp
              LRestart = Ltemp
          ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


! -- -- Branch NEG -- -- Branch NEG -- -- Branch NEG -- -- Branch NEG -- -- Branch NEG
! NEG (U(-a),L(-a))                                   
           WRITE(*,*) "|-|-|-| NEGATIVE Branch for tangent ",IT

            UTemp = UBif
            Ltemp = LUBif
! TEST SOL POL MULIPLE NPOI       
           CALL ANMMultipleSolPol( UTemp, Ltemp, UBranch, LBranch, OMan, -amaxpolyBr,      &
                                   NBSolByStep, NODEOUT, Step, FlowPerm, NSDOFs, dim)
           CALL CalCondLimCoeff( UTemp, Ltemp , 1 , FlowSolution_Init )
                                  
! - - - - - - - - - - - - - - - - - - OUT - - - - - - - - - - - - - - - - - - - - - - - 
       CDBx = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 1) )
       CDBy = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 2) )
       CDBz = 0.0_dp
       if (dim == 3)  CDBz = UTemp(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 3) )       
       WRITE(144,903)  Step,Ltemp,CDBx,CDBy,CDBz       
       CALL FLUSH(144)
! - - - - - - - - - - - - - - - - - - FILE OUT - - - - - - - - - - - - - - - - - - - - -
! ! ! !            Binary  = .FALSE. ! FIRST, after BINARY TRUE IS BETTER FOR PRECISION
           SaveAll = .TRUE. 
           Solver % Mesh % SavesDone = 0  !  Trick to  enforce new file!!
           
           FlowSolution = UTemp
           FlowSol % Values => FlowSolution  
           TIC = CPUTime()           
           SaveCount = SaveResult( 'BIF'//trim(str(BFN))//'BR'//trim(str(IT))//'N.restart', &
                                    Solver % Mesh,Step,Ltemp,Binary,SaveAll,FreeSurfaceFlag ) 
           TIC = CPUTime() - TIC
           WRITE(TIMEFILEFICH,*) 'Step ', Step, '- SaveResult',TIC           
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! ! ! !            FlowSolution = U0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
! Should be better to restart from a file...... waiting with a user choice at least, but for each bif the same choice...
          IF ((SYMRESTART.AND.SYMBRANCH).AND.(Ltemp>LUBif)) THEN
              URestart = UTemp
              LRestart = Ltemp
          ELSE IF ((SYMRESTART.EQV..FALSE.).AND.(SYMBRANCH.EQV..FALSE.)) THEN
              URestart = UTemp
              LRestart = Ltemp
          ENDIF
!TODO TODO : SAVE RESTART FILE AT Utemp LTemp          
! why not save FlowIni in a file as dirichlet_node_profile_save
         ENDDO

! - - - - - - - - - - - - DEALLOCATE - - - - - - - - - - - - - -
!          deallocate(StiffAugW % RHS)
!          deallocate(StiffAugW)
!          StiffAugW => Null()

! -------- FREE MUMPS
         WRITE(6,*) "Branches : HMumps_Free_mumpsIDL"
         CALL HMumps_Free_mumpsIDL(StiffAugW)         
         IF(ASSOCIATED(StiffAugW)) THEN
          write(6,*) "COUCOU CALL FreeMatrix(StiffAugW)"
          CALL FLUSH(6)          
          CALL FreeMatrix(StiffAugW)
          
          write(6,*) "StiffAugW => NULL()"          
          CALL FLUSH(6)          
          StiffAugW => NULL()
         END IF            
! - - - - - - - - - - - - DEALLOCATE - - - - - - - - - - - - - -
         
! CHECK A BRANCH TO SEE....
!          U0 = UTemp     
!          lambda0 = Ltemp
         ! Utemp is the last solution on a branch
         U0 = URestart     
         lambda0 = LRestart
         write(*,*) ' alphaBIFSENSPARCOURS=',alphaBIFSENSPARCOURS
         write(*,*) ' RESTART FROM BRANCH at LRestart=',LRestart
! ! ! ! !          BIFSTAGEOMDETECT = .FALSE. 
         alphaBIF = 0.0_dp   
!********************************************************************************         


! ! ! ! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !
! ! !        IF (alphaVBIF < 0) THEN
! ! !          Amax = ANMseriesRoV( UManProp, OMan-1, TolMan )
! ! !          write(*,*) ' alphaBIF=',alphaBIF         
! ! !          write(*,*) ' PF CritGeom : Bifurcation Behind : ANMSolPol Amax=',Amax
! ! !          CALL ANMSolPol( U0, lambda0, UManProp, LambdaUProp, OMan-1, Amax )               
! ! !          alphaBIF = 0.0_dp         
! ! !          
! ! ! ! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !
! ! !        ELSE   
! - - - - - - 

        
! ! !       ENDIF       
! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !
! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !! ! ! !

         ENDIF   !  CASE LP OR BP
       ENDIF     !  CASE CONT Classic or SINGULAR POINT

       
       CALL CalCondLimCoeff( U0, lambda0, 1, FlowSolution_Init )       
       FlowSolution = U0
       FlowSol % Values => FlowSolution
       
! ! ! ! ! !      Construction de la série INDICATEUR BIFSTA
! ! ! ! !        write(*,*) 'INDIC BifStead : ANMSolPol'
! ! ! ! !        CALL ANMSolPol( INDDU0, INDMu0, INDDUserie, INDMuserie, OMan, amax )               
! ! ! ! !        CALL CalCondLimCoeff( INDDU0, 0.0_dp, 1,FlowSolution_Init )

       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - - - - "
       WRITE(*,*) " - - - -      BIF STA INDICATOR  INFOS         - - - - "              
       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - - - - "       
       WRITE(*,*) " - amaxpolyUL    : ",amaxpolyUL
!        WRITE(*,*) " - amaxpolyINDDU : ",amaxpolyINDDU
       WRITE(*,*) " - amax          : ",amax 
         IF (amax<0) THEN
           WRITE(*,*) " - MOVING FORWARD IN Lambda ie <U1,U> < 0 "       
         ENDIF       
       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - - - - "       
       WRITE(*,*) " - Lambda0       : ",lambda0
       WRITE(*,*) " - Indicator Mu0 : ",INDMu0
       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - - - - "       
!        WRITE(*,*) " - Critr4DiagBif : ",CritrDiagBif
       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - - - - "
       WRITE(*,*) " - - - - - - - - - - - - - - - - - - - - - - - - - - - "
       
       write(*,*) 'sortie fichier 143'
!        WRITE(143,901)  Step,lambda0,amaxpolyUL,sumcrit(1),sumcrit(2),alphaVBIF,LVBif,       & 
!                        PPole,LambdaCritPad,amaxPad,L0Pad
       WRITE(143,901)  Step,lambda0,amaxpolyUL,sumcrit(1),sumcrit(2),calcalphabif(3),LVBif,       & 
                       PPole,LambdaCritPad,amaxPad,L0Pad
       CALL FLUSH(143)          
!------------------------------------------------------------------------------
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  
! ________  .__              __________.__  _____ 
! \______ \ |__|____     ____\______   \__|/ ____\
!  |    |  \|  \__  \   / ___\|    |  _/  \   __\ 
!  |    `   \  |/ __ \_/ /_/  >    |   \  ||  |   
! /_______  /__(____  /\___  /|______  /__||__|   
!         \/        \//_____/        \/      
! 2D NSDOF = 3 => 2 1 0       
! 3D NSDOF = 4 => 3 2 1 0       
! Exemple Uy = composante 2 : 2D => 1 = (3-2), 3D => 2 = (4-2)
! Exemple Ux = composante 1 : 2D => 1 = (3-1)+1, 3D => 2 = (4-1)+1
!        CritrDiagBif = FlowSolution(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - NdOutCOMP) )
       CDBx = FlowSolution(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 1) )
       CDBy = FlowSolution(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 2) )
       CDBz = 0.0_dp
       if (dim == 3)  CDBz = FlowSolution(NSDOFs*FlowPerm(NODEOUT) - (NSDOFs - 3) )       
       
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!------------------------------------------------------------------------------
       
       
       WRITE(144,903)  Step,lambda0,CDBx,CDBy,CDBz
       CALL FLUSH(144)
       
       
 


!------------------------------------------------------------------------------
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! ! ! ! ! ! ! RESIDUAL
! ! ! ! ! !        WRITE(*,*)  '| - - - - - COMPUTERESIDUALVECTOR  - - - - - - -'
! ! ! ! ! !        ResVec=0.0_dp      
! ! ! ! ! !        CALL PressureFree( U0, NoPressure, Utemp )
! ! ! ! ! !        CALL CalCondLimCoeff( Utemp, lambda0, 1, FlowSolution_Init )       
! ! ! ! ! ! 
! ! ! ! ! !        CALL HMFlowResidual( Utemp, lambda0, ResVec, ResVecTMP, Model, Model % Mesh , &
! ! ! ! ! !                             FlowPerm, NSDOFs, FNORMG , HkResidualNorm)
! ! ! ! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )
! ! ! ! ! !        
! ! ! ! ! !        HMResidualU0jp1 = DSQRT( DOT_PRODUCT( ResVec,ResVec ))    
! ! ! ! ! !        CALL PressureFree( ResVec, NoPressure, ResVec )
! ! ! ! ! !        
! ! ! ! ! !        HMResidualU0jp1NoC = DSQRT( DOT_PRODUCT( ResVec,ResVec ))       
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               ! ADD STUFF IN VTU FILES ! TEST               
! ! ! ! ! ! !           ResidualPointer => ResVec
! ! ! ! ! !            ResidualPointer => ResVec( 1: SIZE(ResVec) : NSDOFs )
! ! ! ! ! !            CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! ! !                   'Residual 1', 1,  ResidualPointer, FlowPerm) 
! ! ! ! ! ! 
! ! ! ! ! !            ResidualPointer => ResVec( 2: SIZE(ResVec) : NSDOFs )
! ! ! ! ! !            CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! ! !                   'Residual 2', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! ! ! 
! ! ! ! ! !            IF ( NSDOFs == 3 ) THEN
! ! ! ! ! !              ResidualPointer => ResVec( 3: SIZE(ResVec) : NSDOFs )
! ! ! ! ! !              CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! ! !                   'Residual Div', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! ! !            ELSE
! ! ! ! ! !              ResidualPointer => ResVec( 3: SIZE(ResVec) : NSDOFs )
! ! ! ! ! !              CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh,  Solver, &
! ! ! ! ! !                   'Residual 3', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! ! !              
! ! ! ! ! !              ResidualPointer => ResVec( 4: SIZE(ResVec) : NSDOFs )
! ! ! ! ! !              CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! ! !                   'Residual Div', 1,  ResidualPointer, FlowPerm ) 
! ! ! ! ! !            ENDIF
! ! ! ! ! !            ResidualPointer => ResVec
! ! ! ! ! !            CALL VariableAdd( Model % Mesh % Variables,  Model % Mesh, Solver, &
! ! ! ! ! !                 'Residual', NSDOFs, ResidualPointer ,FlowPerm )           
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
               !---------------------------------------------------------------------
                    
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$       
!        FlowSolution = -11._dp
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
       
!        FlowSolution = ResVec
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )  
!        CALL ListRemove( Solver % Values,'Output File Name')       
!        FlowSolution = U0
!        FlowSol % Values => FlowSolution
!        CALL ResultOutputSolver( Model,Solver,dt,Transient )         
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$       
! ! ! ! !        WRITE(*,*)  '|RRA1      ||NODALFORCES|| =',FNORMG
! ! ! ! !        WRITE(*,*)  '|RRA1     COMPUTED ResNorm =',HkResidualNorm       
! ! ! ! !        WRITE(*,*)  '|RRA1   ||HMFlowResidual|| =',HMResidualU0jp1
! ! ! ! !        WRITE(*,*)  '|RRA1||HMFlowResidualNoc|| =',HMResidualU0jp1NoC
! ! ! ! ! ! ! ! ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! ! !  CORRECTION WITH CLASSIC SOLVER:  PARAMETERS RED FROM SIF
! ! ! !        FlowSolution = U0       
! ! ! !        CALL CalCondLimCoeff( FlowSolution , 0.0_dp, 0, Lambda0 * FlowSolution_Init )        
! ! ! ! !        CALL CalCondLimCoeff( FlowSolution , 0.0_dp, 0,  FlowSolution_Init )        
! ! ! !        FlowSol % Values => FlowSolution  
! ! ! ! 
! ! ! ! !      4CORRECTION : Ucontour = L0 x Ud              
! ! ! !        CALL BoundaryConditionMultCOEFF(Model, Lambda0,.FALSE.)
! ! ! ! ! CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR       
! ! ! !        CALL FlowSolverCorrection( Model,Solver,dt,TransientSimulation,Lambda0)
! ! ! ! ! CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR       
! ! ! ! !        CALL BoundaryConditionMultCOEFF(Model, Lambda0,.TRUE.) !  RESET TO Ud
! ! ! !        CALL BoundaryConditionMultCOEFF(Model, 1._dp/Lambda0,.FALSE.) !  RESET TO Ud OR GET MAX
! ! ! ! ! SWAP FOR CORRECTED SOLUTIN
! ! ! !        FlowSolution => FlowSol % Values       
! ! ! !        U0 = FlowSolution       
! ! ! ! ! ! ! ! !------------------------------------------------------------------------------       
! ! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !        WRITE(*,*)  '| - - - - - COMPUTERESIDUALVECTOR AFTER CORRECTION- - - - - - -'
! ! ! ! ! 
! ! ! ! !        
! ! ! ! !        ResVec=0.0_dp
! ! ! ! !        CALL HMFlowResidual( U0, lambda0, ResVec, ResVecTMP, Model, Solver % Mesh , &
! ! ! ! !                             FlowPerm, NSDOFs, FNORMG , ResNorm)
! ! ! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )
! ! ! ! !        

       if (ResOUTornot) THEN
! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! ! ! ! !  _____        _____       __      _______ ________          __
! ! ! ! ! |  __ \ /\   |  __ \     /\ \    / /_   _|  ____\ \        / /
! ! ! ! ! | |__) /  \  | |__) |   /  \ \  / /  | | | |__   \ \  /\  / / 
! ! ! ! ! |  ___/ /\ \ |  _  /   / /\ \ \/ /   | | |  __|   \ \/  \/ /  
! ! ! ! ! | |  / ____ \| | \ \  / ____ \  /   _| |_| |____   \  /\  /   
! ! ! ! ! |_| /_/    \_\_|  \_\/_/    \_\/   |_____|______|   \/  \/                      
! ! ! ! !
! ! ! ! ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
! LAST VISU
       FlowSolution = U0
       FlowSol % Values => FlowSolution
       CALL ResultOutputSolver( Model,Solver,dt,Transient )    
! ! ! ! 
! ! ! ! ! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! ! ! ! !------------------------------------------------------------------------------
       ENDIF
        
! LAST VISU
       FlowSolution = U0
       FlowSol % Values => FlowSolution
! 4 RESTART       
! ! !        CALL CalCondLimCoeff( FlowSolution , 0.0_dp, 0, FlowSolution_Init )  ! 4 RESTART
       
! ! ! ! !        WRITE(*,*)  '|RRA2      ||NODALFORCES|| =',FNORMG
! ! ! ! !        WRITE(*,*)  '|RRB     COMPUTED ResNorm =',ResNorm       
! ! ! ! !        WRITE(*,*)  '|RRA2   ||HMFlowResidual|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))
! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !------------------------------------------------------------------------------

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
! ! ! ! !------------------------------------------------------------------------------
!        WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,HMResidualU0jp1,HkResidualNorm
       WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,ManualResidualU0jstab
       CALL FLUSH(RESIDUALFICH)       
       
       CALL FLUSH(TIMEFILEFICH)

       st = CPUTIme()-st
       totst = totst + st

       WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' TOTAL TIME:    (s)', st, totst
       CALL Info( 'ManFlowSolve', Message, Level=4 )
!      
       WRITE(*,*) '------------------------------------------------------------------------------'
       WRITE(*,*) ' STEP FINISHED : ', Step
       WRITE(*,*) '------------------------------------------------------------------------------'
       write(*,*) 'lambda      :', lambda0
       write(*,*) 'Norme U     :', L2UMan       
       Write(*,*) 'amaxpolyUL  :', amaxpolyUL
       Write(*,*) 'amaxPad     :', amaxPad
       write(*,*) 'CochMed     :', sumcrit(1),sumcrit(2),alphaVBIF
       WRITE(*,*) '------------------------------------------------------------------------------'
       WRITE(*,*) '- * - * - *- * - * - *- * - * - *- * - * - *- * - * - *- * - * - *- * - * - * '
       WRITE(*,*) '------------------------------------------------------------------------------'       

       n = NSDOFs * LocalNodes

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!     This hack is needed  cause of the fluctuating pressure levels
!------------------------------------------------------------------------------

!write(*,*) 'This hack is needed  cause of the fluctuating pressure levels'

       IF ( NonlinearRelax /= 1.0d0 ) THEN
         IF ( CompressibilityModel == Incompressible ) THEN
           s = FlowSolution(NSDOFs)
           FlowSolution(NSDOFs:n:NSDOFs) = FlowSolution(NSDOFs:n:NSDOFs)-s
           PSolution(NSDOFs:n:NSDOFs)=PSolution(NSDOFs:n:NSDOFs)-PSolution(NSDOFs)
         END IF

         FlowSolution(1:n) = (1-NonlinearRelax)*PSolution(1:n) + NonlinearRelax*FlowSolution(1:n)                    
       
         IF ( CompressibilityModel == Incompressible ) THEN
            FlowSolution(NSDOFs:n:NSDOFs)=FlowSolution(NSDOFs:n:NSDOFs)+s
         END IF

         RelaxBefore = GetLogical( Solver % Values, &
              'Nonlinear system Relaxation Before', gotIt )
         IF ( .NOT.gotIt .OR. RelaxBefore ) THEN
           CALL ListAddLogical( Solver % Values, 'Skip Compute Nonlinear Change', .FALSE. )
           CALL ComputeChange( Solver, .FALSE., n, FlowSolution, PSolution )
         END IF
       END IF
       RelativeChange = Solver % Variable % NonlinChange

!------------------------------------------------------------------------------

! !        WRITE( Message, * ) 'Result Norm     : ',Solver % Variable % Norm
! !        CALL Info( 'ManFlowSolve', Message, Level=4 )
! !        WRITE( Message, * ) 'Relative Change : ',RelativeChange
! !        CALL Info( 'ManFlowSolve', Message, Level=4 )
! ! 
! !        IF ( RelativeChange < NewtonTol .OR. &
! !              iter > NewtonIter ) NewtonLinearization = .TRUE.
! ! 
! !        IF ( RelativeChange < NonLinearTol .AND. Iter<NonlinearIter ) EXIT
! ! 
! ! !------------------------------------------------------------------------------
! ! !     If free surfaces in model, this will move the nodal points
! ! !------------------------------------------------------------------------------
! !        IF ( FreeSurfaceFlag ) THEN
! ! 
! !          IF ( RelativeChange < FreeSTol .OR. iter > FreeSIter ) ComputeFree = .TRUE.           
! ! 
! !          IF ( ComputeFree ) THEN
! !            Relaxation = GetCReal( Solver % Values, &
! !              'Free Surface Relaxation Factor', GotIt )
! ! 
! !            IF ( .NOT.GotIt ) Relaxation = 1.0d0
! ! 
! !            MBFlag = GetLogical(GetSolverParams(), 'Internal Move Boundary', GotIt)
! !            IF ( MBFlag .OR. .NOT. GotIt ) CALL MoveBoundary( Model, Relaxation )
! !          END IF
! !        END IF
!------------------------------------------------------------------------------

! 20161014 : sauvegarde du dernierPoints
! - - - - - - - - - - - - - - - - - - FILE OUT - - - - - - - - - - - - - - - - - - - - -
           Binary  = .TRUE. ! FIRST, after BINARY TRUE IS BETTER FOR PRECISION
           SaveAll = .TRUE. 
           Solver % Mesh % SavesDone = 0  !  Trick to  enforce new file!!          
           SaveCount = SaveResult( 'DernierPoint.restart',Solver % Mesh,Step,lambda0,Binary,SaveAll,FreeSurfaceFlag ) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     END DO ! ANMStep
!  ______ _   _ _____             _   _ __  __    _____ _               _                       
! |  ____| \ | |  __ \      /\   | \ | |  \/  |  / ____| |             | |                      
! | |__  |  \| | |  | |    /  \  |  \| | \  / | | (___ | |_ ___ _ __   | |     ___   ___  _ __  
! |  __| | . ` | |  | |   / /\ \ | . ` | |\/| |  \___ \| __/ _ \ '_ \  | |    / _ \ / _ \| '_ \ 
! | |____| |\  | |__| |  / ____ \| |\  | |  | |  ____) | ||  __/ |_) | | |___| (_) | (_) | |_) |
! |______|_| \_|_____/  /_/    \_\_| \_|_|  |_| |_____/ \__\___| .__/  |______\___/ \___/| .__/ 
!                                                              | |                       | |    
!                                                              |_|                       |_|      




    CLOSE(142)
    CLOSE(CRITGEOMFICH)
    
    
    WRITE(144,903)  Step,lambda0      
    
    CLOSE(BIFSTADIAGFICH)
    CLOSE(MUMPSFICH)
    CLOSE(ABETFICH)    
    CLOSE(PADEFICH)    
    CLOSE(RESIDUALFICH)    
    CLOSE(TIMEFILEFICH)    
    CLOSE(DBSKPFICH)    

    
!------------------------------------------------------------------------------    
!     CALL ListAddConstReal( Solver % Values, &
!         'Nonlinear System Relaxation Factor', NonlinearRelax )

!     IF (ListGetLogical(Solver % Values,'Adaptive Mesh Refinement',GotIt)) &
!       CALL RefineMesh( Model,Solver,FlowSolution,FlowPerm, &
!          FlowInsideResidual, FlowEdgeResidual, FlowBoundaryResidual ) 

!------------------------------------------------------------------------------
!     CALL CheckCircleBoundary()
!------------------------------------------------------------------------------

    IF ( UseLocalCoords ) THEN
       Coordinates = ModelCoords
       Model % DIMENSION = ModelDim
    END IF
    
   
    write (*,*) '\\\\\\\\\ END /////////'

! CRITBIF
  900 FORMAT( '#STEP', 7X, 'LAMBDA', 12X,'Amax',7X,'|',6X,'Sum1',12X,'Sum2',12X,'alphaBIF', &
              8X,'LcritCM',7X,'|',3X,'MinPole >0',8X,'LcritPad',8X,'AmxPad',8X,'L0Pad') 
  901 FORMAT( 1X, I4, 2X, G15.8, 2X, G15.8, 2X ,'|', G15.8, 2X, G15.8, 2X, G15.8,            &
              2X, G15.8,2X, '|', G15.8, 2X, G15.8, 2X, G15.8, 2X, G15.8) 
!   902 FORMAT( 1X,'STEP', 7X, 'LAMBDA', 12X,'Mu0', 12X,'Udiagbif')
!   903 FORMAT( 1X, I4, 2X, G15.8, 2X, G15.8, 2X, G15.8, 2X, G15.8, 2X, G15.8)  
! DIAG BIFSTA
  902 FORMAT( '#STEP', 7X, 'LAMBDA', 12X,'Udiagbif')
  903 FORMAT( 1X, I4, 2X, G15.8, 2X, G15.8, 2X, G15.8, 2X, G15.8)
! PADE
  904 FORMAT( '#STEP', 7X, 'LAMBDA', 12X,'Amaxpoly',7X,'AmxPade',7X,'|',5X,'POLE',5X,'LcritPad') 
  905 FORMAT( 1X,I4,2X,G15.8,2X,G15.8,2X,G15.8,2X,'|',2X,G15.8,2X,G15.8) 
! RESIDUAL
  906 FORMAT( '#STEP', 7X, 'LAMBDA', 12X,'||U0||', 12X,'||L(U0) + Q(U0,U0)||',7X, &
              '||LSTAB(U0) + QSTAB(U0,U0)||')   
!               'HMResidualU0j+1',7X,'|',7X,'Hk Residual Nurm U0j+1') 

  907 FORMAT( 1X,I4,2X, G15.8,2X, G15.8,2X,G15.8,2X, 2X, G15.8) 
              
! BifSta indicators (order 0 first to test) works using LT dirichlet and Norm condition
  908 FORMAT( 'STEP', 7X, 'LAMBDA', 12X,'Cadou2006',7X,'Vannucci',7X,'NormConditionAsCadou2012') 
  909 FORMAT( 1X,I4,2X,      G15.8,2X,      G15.8,2X,2X,    G15.8,2X,   G15.8) 




CONTAINS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CheckCircleBoundary()
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: x,y,phi,x0,y0,r
      LOGICAL :: GotIt
      INTEGER :: i,j,k,l
!------------------------------------------------------------------------------

      l = 0
      DO i=1,Model % NumberOfBCs
         IF ( .NOT.ListgetLogical( Model % BCs(i) % Values, &
                  'Circle Boundary', GotIt ) ) CYCLE

         x0 = ListGetConstReal( Model % BCs(i) % Values, 'Circle X', GotIt )
         IF ( .NOT. GotIt ) x0 = 0.0d0

         y0 = ListGetConstReal( Model % BCs(i) % Values, 'Circle Y', GotIt )
         IF ( .NOT. GotIt ) y0 = 0.0d0

         R  = ListGetConstReal( Model % BCs(i) % Values, 'Circle R', GotIt )
         IF ( .NOT. GotIt ) R = 1.0d0

         DO j=Solver % Mesh % NumberOfBulkElements+1, &
            Solver % Mesh % NumberOfBulkElements+ &
               Solver % Mesh % NumberOfBoundaryElements
            Element => Solver % Mesh % Elements(j)
            IF ( Element % BoundaryInfo % Constraint &
                 /= Model % BCs(i) % Tag ) CYCLE

            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes
            DO k=1,n
               x = Solver % Mesh % Nodes % x(NodeIndexes(k)) - x0
               y = Solver % Mesh % Nodes % y(NodeIndexes(k)) - y0

               phi = ATAN2( y,x )
               x = R * COS( phi ) 
               y = R * SIN( phi ) 

               Solver % Mesh % Nodes % x(NodeIndexes(k)) = x + x0
               Solver % Mesh % Nodes % y(NodeIndexes(k)) = y + y0
            END DO
            l = l + 1
        END DO
     END DO

     IF ( l > 0 ) THEN
        WRITE( Message, * ) 'Elements on Circle', l
        CALL Info( 'ManFlowSolve', Message, Level=6 )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CheckCircleBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE ANMPathFollowing
!------------------------------------------------------------------------------


!> Initialization of the main solver: AdvectionDiffusionSolver
!------------------------------------------------------------------------------
   SUBROUTINE ManFlowSolver_init( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
     TYPE(Model_t), TARGET :: Model    !< All model information (mesh, materials, BCs, etc...)
     REAL(KIND=dp) :: Timestep         !< Timestep size for time dependent simulations
     LOGICAL :: TransientSimulation    !< Steady state or transient simulation
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Params

     Params => GetSolverParams()
     CALL ListAddInteger( Params,'Time Derivative Order',1 )

   END SUBROUTINE ManFlowSolver_init





 
!> \} 
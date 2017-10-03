!> \ingroup IMNS
!> \}
! # YG@2016
! # yann.guevel@univ-ubs.fr
! # Part of ANM continuation/detection/branch-switching in ELMER FEM
! # Jean-Marc Cadou / Gregory girault / Yann Guevel
! # Institut de Recherche Dupuy de Lome - Universite Bretagne Sud


!------------------------------------------------------------------------------
MODULE DiscreteOperators

!------------------------------------------------------------------------------

  USE SolverUtils
  USE DefUtils
  !USE ANMToolBox ! COMMENTED ON 2016 03 21 : test compil modules
  USE FEMtoolBox    ! New module 20160324 : calcondlim pour le moment
  

!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
    

CONTAINS
!------------------------------------------------------------------------------
!>  CFD 
!>  OperatorsLtF : Lt = L() + Q(,U) + Q(U,)  et F avec PRCONDLIM
!>    |-> MANTangentOperatorCompose (MODULE ANM...)

!>  ANMPFrhsFQ   : CONTINUATION CLASSIQUE Second membre MAN comme la somme des termes croisés de séries 
!>    |-> MANFQelemOPTI

!>  ANMBranchrhsFQ : CONTINUATION POST BIFURCATION STATIONNAIRE Second membre MAN comme la somme des termes croisés de séries 
!>    |-> MANBranchFQelemOPTI


!>  ANMSBIrhsFnl : Indicateur cadou2006 : second membre MAN
!>    |-> ANMSBIrhsFnlElemOPTI


!>  ANMQAB       : Q(A,B) comme u.grad u
!>    |-> MANQABelem



!>  RESIDUAL : Lt(U) - Q(U,U) - Lambda F  , F=0, unless bodyforce...
!>    |-> RESIDUAL_ELEM : Loop over element + Elementary residual NS


!>  HMFlowResidual
!>    |-> HMLocalFlowInsideResidual

!------------------------------------------------------------------------------
!  _______       _   _  _____ ______ _   _ _______ 
! |__   __|/\   | \ | |/ ____|  ____| \ | |__   __|
!    | |  /  \  |  \| | |  __| |__  |  \| |  | |   
!    | | / /\ \ | . ` | | |_ |  __| | . ` |  | |   
!    | |/ ____ \| |\  | |__| | |____| |\  |  | |   
!    |_/_/    \_\_| \_|\_____|______|_| \_|  |_|                                                  
!  _____                                _      _          
! (  _  )                              ( )    (_ )        
! | (_) |  ___   ___    __    ___ ___  | |_    | |  _   _ 
! |  _  |/',__)/',__) /'__`\/' _ ` _ `\| '_`\  | | ( ) ( )
! | | | |\__, \\__, \(  ___/| ( ) ( ) || |_) ) | | | (_) |
! (_) (_)(____/(____/`\____)(_) (_) (_)(_,__/'(___)`\__, |
!                                                  ( )_| |
!                                                  `\___/'
!
!---------------------- Composition de l'opérateur Tangent Lt -----------------
!------------------------- Lt(*)=L(*) + Q(*,U0) + Q(U0,*) ---------------------
!---------------------- Second membre en consequence          -----------------

    SUBROUTINE OperatorsLtF( StiffMatrix, ForceVector, USOL, NSDOFs,               &
                             MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
                             MeshSol, DensitySol, LocalNodes, dt, Transient,       &
                             PseudoCompressibilityScale, TempSol ,TempPerm,        &
                             Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &
                             CompressibilityModel,DivDiscretization,               &
                             GradPDiscretization,NewtonLinearization,              &
                             Gravity                                               &                           
                             ) 
!------------------------------------------------------------------------------
     USE Differentials  
     USE MaterialModels
     
!------------------------------------------------------------------------------     
    USE NavierStokes
    USE NavierStokesGeneral
    USE NavierStokesCylindrical     
!------------------------------------------------------------------------------    
! MODULE FAIT MAISON     
!      USE ANMToolBox
!------------------------------------------------------------------------------    


!------------------------------------------------------------------------------
! ROUTINE ARGUMENTS    
     TYPE(Matrix_t),POINTER   :: StiffMatrix
     REAL(KIND=dp),POINTER    :: ForceVector(:),MeshVelocity(:),Temperature(:)
     REAL(KIND=dp)    :: USOL(:)
     INTEGER          :: NSDOFs, LocalNodes
     INTEGER, POINTER :: FlowPerm(:),TempPerm(:), MeshPerm(:)
     TYPE(Solver_t)   :: Solver
     TYPE(Model_t)    :: Model
     TYPE(Variable_t), POINTER :: MeshSol,DensitySol,TempSol
     REAL(KIND=dp)    :: dt,PseudoCompressibilityScale
     LOGICAL          :: Transient,DivDiscretization,GradPDiscretization,NewtonLinearization
     REAL(KIND=dp), POINTER :: TempPrev(:)
     INTEGER :: CompressibilityModel
!------------------------------------------------------------------------------
! LOCAL
      
    INTEGER :: t,k,n,nb,nd,i,istat,j
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
    INTEGER :: body_id,bf_id,eq_id
    TYPE(ValueList_t),POINTER :: Material, BC, BodyForce, Equation
    CHARACTER(LEN=MAX_NAME_LEN) :: CompressibilityFlag, StabilizeFlag, VarName
    INTEGER :: ModelCoords, ModelDim
    LOGICAL :: Stabilize, GotForceBC, GotIt, &
                  MBFlag, Convect  = .TRUE., NormalTangential, RelaxBefore, &
                  ComputeFree=.FALSE., &
                  Rotating    
     LOGICAL :: AllocationsDone = .FALSE., FreeSurfaceFlag, &
         PseudoPressureExists, PseudoCompressible, Bubbles, &
         Porous =.FALSE., PotentialForce=.FALSE., Hydrostatic=.FALSE., &
         MagneticForce =.FALSE., UseLocalCoords, PseudoPressureUpdate     
         
     REAL(KIND=dp), POINTER ::  gWork(:,:), LayerThickness(:), &
                                SurfaceRoughness(:)
     REAL(KIND=dp) :: Gravity(3),AngularVelocity(3)
     REAL(KIND=dp) :: ReferencePressure,SpecificHeatRatio,Tdiff
     
     REAL(KIND=DP), POINTER :: Pwrk(:,:,:)
     TYPE(Nodes_t) :: ElementNodes
     

     

!------------------------------------------------------------------------------
     

     REAL(KIND=dp),ALLOCATABLE :: MASS(:,:),STIFF(:,:), LoadVector(:,:), &
       Viscosity(:),FORCE(:), TimeForce(:), PrevDensity(:),Density(:),   &
       U(:),V(:),W(:),MU(:),MV(:),MW(:), Pressure(:),Alpha(:),Beta(:),   &
       ExtPressure(:),PrevPressure(:), HeatExpansionCoeff(:),            &
       ReferenceTemperature(:), Permeability(:),Mx(:),My(:),Mz(:),       &
       LocalTemperature(:), GasConstant(:), HeatCapacity(:),             &
       LocalTempPrev(:),SlipCoeff(:,:), PseudoCompressibility(:),        &
       PseudoPressure(:), Drag(:,:), PotentialField(:),    &
       PotentialCoefficient(:)

!      SAVE U,V,W,MASS,STIFF,LoadVector,Viscosity, TimeForce,FORCE,ElementNodes,  &
!        Alpha,Beta,ExtPressure,Pressure,PrevPressure, PrevDensity,Density,       &
!        AllocationsDone,LocalNodes, HeatExpansionCoeff,ReferenceTemperature,     &
!        Permeability,Mx,My,Mz,LayerThickness, SlipCoeff, SurfaceRoughness,       &
!        LocalTemperature, GasConstant, HeatCapacity, LocalTempPrev,MU,MV,MW,     &
!        PseudoCompressibilityScale, PseudoCompressibility, PseudoPressure,       &
!        PseudoPressureExists, PSolution, Drag, PotentialField, PotentialCoefficient, &
!        ComputeFree, Indexes     
     
!------------------------------------------------------------------------------
! LOCAL ELEMENTAIRE INIT ET ALLOCATE REMOVE FROM UPPER SUBR

       N = Solver % Mesh % MaxElementDOFs
       NSDOFs = Solver % Variable % DOFs
       
!------------------------------------------------------------------------------       
       IF( AllocationsDone ) THEN
          DEALLOCATE(                                &
               U,  V,  W,                            &
               MU, MV, MW,                           &
               Indexes,                              &
               Pressure,                             &
               PrevPressure,                         &
               PseudoCompressibility,                &
               PrevDensity,Density,                  &
               LayerThickness,                       &
               SurfaceRoughness,                     &
               Permeability,                         &
               Mx,My,Mz,                             &
               SlipCoeff, Drag,                      &
               TimeForce,FORCE, Viscosity,           &
               MASS,  STIFF,                         &
               HeatExpansionCoeff,                   &
               GasConstant, HeatCapacity,            &
               ReferenceTemperature,                 & 
               LocalTempPrev, LocalTemperature,      &
               PotentialField, PotentialCoefficient, &
               LoadVector, Alpha, Beta,   &
               ExtPressure,                          &     
               STAT=istat )
       END IF     
       ALLOCATE( U(N),  V(N),  W(N),                             &
                 MU(N), MV(N), MW(N),                            &
                 Indexes( N ),                                   &
                 Pressure( N ),                                  &
                 PrevPressure( N ),                              &
                 PseudoCompressibility( N ),                     &
                 PrevDensity(N),Density( N ),                    &
                 LayerThickness(N),                              &
                 SurfaceRoughness(N),                            &
                 Permeability(N),                                &
                 Mx(N),My(N),Mz(N),                              &
                 SlipCoeff(3,N), Drag(3,N),                      &
                 TimeForce( 2*NSDOFs*N ),                        &
                 FORCE( 2*NSDOFs*N ), Viscosity( N ),            &
                 MASS(  2*NSDOFs*N,2*NSDOFs*N ),                 &
                 STIFF( 2*NSDOFs*N,2*NSDOFs*N ),                 &
                 HeatExpansionCoeff(N),                          &
                 GasConstant( N ), HeatCapacity( N ),            &
                 ReferenceTemperature(N),                        & 
                 LocalTempPrev(N), LocalTemperature(N),          &
                 PotentialField( N ), PotentialCoefficient( N ), &
                 LoadVector( 4,N ), Alpha( N ), Beta( N ),       &
                 ExtPressure( N ),                               &
                 STAT=istat )
     
!------------------------------------------------------------------------------
!        write(*,*) 'InitializeToZero'
       StiffMatrix => Solver % Matrix
       ForceVector => Solver % Matrix % RHS
       CALL InitializeToZero( StiffMatrix, ForceVector )
       
       NULLIFY(Pwrk) 
       
!------------------------------------------------------------------------------
       bf_id   = -1
       body_id = -1

       CALL StartAdvanceOutput( 'OperatorsLtF', 'Assembly: ' )
       
!------------------------------------------------------------------------------       
! YANN 201501 : ON FORCE CONVECT
           Convect = .TRUE.
! YANN 201501 : ON FORCE CONVECT       
!------------------------------------------------------------------------------
!  ______ _                           _     _                       
! |  ____| |                         | |   | |                      
! | |__  | | ___ _ __ ___   ___ _ __ | |_  | |     ___   ___  _ __  
! |  __| | |/ _ \ '_ ` _ \ / _ \ '_ \| __| | |    / _ \ / _ \| '_ \ 
! | |____| |  __/ | | | | |  __/ | | | |_  | |___| (_) | (_) | |_) |
! |______|_|\___|_| |_| |_|\___|_| |_|\__| |______\___/ \___/| .__/ 
!                                                            | |    
!                                                            |_|    
!------ Boucle sur les éléments du domaine ------------------------------------
       DO t = 1,GetNOFActive()

         CALL AdvanceOutput( t,GetNOFActive() )

         Element => GetActiveElement(t)
         NodeIndexes => Element % NodeIndexes

!------------------------------------------------------------------------------
         IF ( Element % BodyId /= body_id ) THEN
           body_id = Element % BodyId

           eq_id = ListGetInteger( Model % Bodies(body_id) % Values,'Equation', &
                                     minv=1, maxv=Model % NumberOfEquations )
           Equation => Model % Equations(eq_id) % Values

           bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
                                     'Body Force', gotIt, 1, Model % NumberOfBodyForces )
           IF( bf_id > 0 ) THEN
             BodyForce => Model % BodyForces(bf_id) % Values
           END IF



           k = ListGetInteger( Model % Bodies(body_id) % Values, 'Material', &
                                  minv=1, maxv=Model % NumberOfMaterials )
           Material => Model % Materials(k) % Values
             
!------------------------------------------------------------------------------             
           CompressibilityFlag = ListGetString( Material, 'Compressibility Model', GotIt)
                                                  
           IF ( .NOT.GotIt ) CompressibilityModel = Incompressible
           PseudoCompressible = .FALSE.

           SELECT CASE( CompressibilityFlag )
             CASE( 'incompressible' )
               CompressibilityModel = Incompressible

             CASE( 'perfect gas', 'perfect gas equation 1' )
               CompressibilityModel = PerfectGas1

             CASE( 'thermal' )
               CompressibilityModel = Thermal

             CASE( 'user defined' )
               CompressibilityModel = UserDefined1

             CASE( 'pressure dependent' )
               CompressibilityModel = UserDefined2

             CASE( 'artificial compressible' )
               CompressibilityModel = Incompressible 
               PseudoCompressible = .TRUE.

             CASE DEFAULT
               CompressibilityModel = Incompressible
           END SELECT
!------------------------------------------------------------------------------
           Gotit = .FALSE.
           IF ( bf_id > 0 ) THEN
             MagneticForce = ListGetLogical( BodyForce,'Lorentz Force', gotIt )
             Hydrostatic = ListGetLogical( BodyForce,'Hydrostatic Pressure',gotIt )
           END IF
           IF ( .NOT. GotIt ) THEN
             Hydrostatic = ListGetLogical( Equation,'Hydrostatic Pressure',gotIt )
           END IF
!            WRITE(*,*) "SIF VERIF - Hydrostatic=",Hydrostatic
!------------------------------------------------------------------------------
           Rotating = .FALSE.
           IF( bf_id > 0 ) THEN
             gWork => ListGetConstRealArray( BodyForce,'Angular Velocity',GotIt)
             IF ( GotIt ) THEN
               IF( Coordinates == Cartesian ) THEN
                 AngularVelocity = gWork(1:3,1)
                 Rotating = .TRUE.
               ELSE
                 CALL Fatal('ManFlowSolve','Rotating coordinate implemented only for cartesian coodinates')
               END IF
             ELSE
               AngularVelocity = 0.0_dp
             END IF
           END IF
         END IF
!------------------------------------------------------------------------------

         n = GetElementNOFNodes()
         nb = GetElementNOFBDOFs()
         nd = GetElementDOFs( Indexes )
     
         CALL GetElementNodes( ElementNodes )
!------------------------------------------------------------------------------         
         SELECT CASE( NSDOFs )
           CASE(3)
             U(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd))-2)
             V(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd))-1)
             W(1:nd) = 0.0d0
           CASE(4)
             U(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd))-3)
             V(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd))-2)
             W(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd))-1)
         END SELECT
!------------------------------------------------------------------------------         
!          write(*,*) t,U,V
!------------------------------------------------------------------------------
         MU(1:nd) = 0.0d0
         MV(1:nd) = 0.0d0
         MW(1:nd) = 0.0d0
         IF ( ASSOCIATED( MeshVelocity ) ) THEN
            SELECT CASE( MeshSol % DOFs )
            CASE(2)
               IF ( ALL( MeshPerm( Indexes(1:nd) ) > 0 ) ) THEN
                  MU(1:nd) = MeshVelocity(2*MeshPerm(Indexes(1:nd))-1)
                  MV(1:nd) = MeshVelocity(2*MeshPerm(Indexes(1:nd))-0)
               END IF
            CASE(3)
               IF ( ALL( MeshPerm( NodeIndexes ) > 0 ) ) THEN
                  MU(1:nd) = MeshVelocity(3*MeshPerm(Indexes(1:nd))-2)
                  MV(1:nd) = MeshVelocity(3*MeshPerm(Indexes(1:nd))-1)
                  MW(1:nd) = MeshVelocity(3*MeshPerm(Indexes(1:nd))-0)
               END IF
            END SELECT
         END IF
!------------------------------------------------------------------------------         
         LocalTemperature = 0.0d0
         LocalTempPrev    = 0.0d0
         IF ( ASSOCIATED( TempSol ) ) THEN
            IF ( ALL( TempPerm(NodeIndexes) > 0 ) ) THEN
               LocalTemperature(1:nd) = Temperature( TempPerm(Indexes(1:nd)) )
               IF ( Transient .AND. CompressibilityModel /= Incompressible) THEN
                 LocalTempPrev(1:nd) = TempPrev( TempPerm(Indexes(1:nd)) )
               END IF
            END IF
         END IF
!------------------------------------------------------------------------------
         
         ReferencePressure = 0.0d0

         PrevDensity = 0.0d0
         Density = 0.0d0
!------------------------------------------------------------------------------
         SELECT CASE( CompressibilityModel )
!------------------------------------------------------------------------------
           CASE( Incompressible )
!------------------------------------------------------------------------------
             Pressure(1:nd) = USOL( NSDOFs*FlowPerm(Indexes(1:nd)) )
             Density(1:n) = GetReal( Material, 'Density' )

             IF(PseudoCompressible) THEN
               Pressure(1:n) = GetReal( Material,'Artificial Pressure', GotIt )
               IF(.NOT. GotIt) THEN
                 Pressure(1:nd) = PseudoPressure(FlowPerm(Indexes(1:nd))) 
               ELSE
                 Pressure(n+1:nd) = 0.0d0
               END IF
               PseudoCompressibility(1:n) = PseudoCompressibilityScale * &
                   GetReal(Material,'Artificial Compressibility', GotIt )
               IF(.NOT. gotIt) PseudoCompressibility(1:n) = 0.0d0
             END IF
!------------------------------------------------------------------------------
           CASE( PerfectGas1 )

              ! Use  ReferenceTemperature in .sif file for fixed temperature
              ! field. At the moment can not have both fixed T ideal gas and
              ! Boussinesq force:
              !-------------------------------------------------------------
              IF ( .NOT. ASSOCIATED( TempSol ) ) THEN
                 LocalTemperature(1:n) = GetReal( Material, &
                   'Reference Temperature' )
                 LocalTempPrev = LocalTemperature
              END IF

              HeatCapacity(1:n) = GetReal( Material, 'Heat Capacity', GotIt )


              ! Read Specific Heat Ratio:
              !--------------------------
              SpecificHeatRatio = ListGetConstReal( Material, &
                     'Specific Heat Ratio', GotIt )
              IF ( .NOT.GotIt ) SpecificHeatRatio = 5.d0/3.d0


              ! For an ideal gas, \gamma, c_p and R are really a constant
              ! GasConstant is an array only since HeatCapacity formally is
              !------------------------------------------------------------
              GasConstant(1:n) = ( SpecificHeatRatio - 1.d0 ) *  &
                   HeatCapacity(1:n) / SpecificHeatRatio


              ! For ideal gases take pressure deviation p_d as the
              ! dependent variable: p = p_0 + p_d
              ! Read p_0
              !---------------------------------------------------
              ReferencePressure = ListGetConstReal( Material, &
                      'Reference Pressure', GotIt )
              IF ( .NOT.GotIt ) ReferencePressure = 0.0d0

              Pressure(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd)))
              IF ( Transient ) THEN
                PrevPressure(1:nd) = Solver % Variable % PrevValues( &
                          NSDOFs*FlowPerm(Indexes(1:nd)),1 )
              END IF
              Density(1:n) = ( Pressure(1:n) + ReferencePressure ) / &
                 ( GasConstant(1:n) * LocalTemperature(1:n) )

           CASE( UserDefined1 )
             Pressure(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd)) )
             IF ( ASSOCIATED( DensitySol ) ) THEN
               Density(1:nd) = DensitySol % Values( DensitySol % Perm(Indexes(1:nd)) )
!                IF ( Transient ) THEN
!                   PrevDensity(1:nd) = DensitySol % PrevValues( &
!                        DensitySol % Perm(Indexes(1:nd)),1)
!                 END IF
             ELSE
               Density(1:n) = GetReal( Material,'Density' )
!                IF ( Transient ) THEN
!                  IF (.NOT.ALLOCATED(pDensity0))  THEN
!                    ALLOCATE(pDensity0(LocalNodes), &
!                             pDensity1(LocalNodes))
!                  END IF

!                  IF ( Timestep==1 ) &
!                      pDensity0(Indexes(1:n)) = Density(1:n)
!                  pDensity1(Indexes(1:n)) = Density(1:n)
!                  PrevDensity(1:n) = pDensity0(Indexes(1:n))
!                END IF
             END IF

           CASE( UserDefined2 )
             Density(1:n) = GetReal( Material,'Density' )
             Pressure(1:nd) = USOL(NSDOFs*FlowPerm(Indexes(1:nd)))

           CASE( Thermal )
             Pressure(1:n) = USOL(NSDOFs*FlowPerm(Indexes(1:nd)))

             HeatExpansionCoeff(1:n) = GetReal( Material, &
               'Heat Expansion Coefficient' )

             ReferenceTemperature(1:n) = GetReal( Material, &
               'Reference Temperature' )

             Density(1:n) = GetReal( Material,'Density' )

             IF( Transient ) THEN
               PrevDensity(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                  ( LocalTempPrev(1:n) - ReferenceTemperature(1:n) ) )
             END IF
             Density(1:n) = Density(1:n) * ( 1 - HeatExpansionCoeff(1:n)  * &
                  ( LocalTemperature(1:n) - ReferenceTemperature(1:n) ) )
!------------------------------------------------------------------------------
         END SELECT
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!        Read in porous media defs
!------------------------------------------------------------------------------
         Porous = ListGetLogical( Material,'Porous Media', GotIt)
         IF(Porous) THEN
           CALL GetRealArray( Material,  Pwrk,'Porous Resistivity',GotIt)
           
           IF( .NOT. GotIt ) THEN
             Drag( 1,1:n) = GetReal( Material,'Porous Resistivity 1',GotIt )
             Drag( 2,1:n) = GetReal( Material,'Porous Resistivity 2',GotIt ) 
             IF( NSDOFs -1 > 2 ) THEN
               Drag( 3,1:n) = GetReal( Material,'Porous Resistivity 3',GotIt ) 
             END IF
           ELSE IF ( SIZE(Pwrk,1) == 1 ) THEN
             DO i=1,NSDOFs-1
               Drag( i,1:n ) = Pwrk( 1,1,1:n )
             END DO
           ELSE 
             DO i=1,MIN(NSDOFs,SIZE(Pwrk,1))
               Drag(i,1:n) = Pwrk(i,1,1:n)
             END DO
           END IF
         END IF
!------------------------------------------------------------------------------
!        Viscosity = Laminar viscosity
!------------------------------------------------------------------------------
         Viscosity(1:n) = GetReal( Material,'Viscosity' )
!------------------------------------------------------------------------------
!        Set body forces, if any
!------------------------------------------------------------------------------
         LoadVector = 0.0D0

         IF ( bf_id > 0 ) THEN
           HeatExpansionCoeff   = 0.0D0
           ReferenceTemperature = 0.0D0

!------------------------------------------------------------------------------
!          Boussinesq body force & gravity
!------------------------------------------------------------------------------
           IF ( ListGetLogical( BodyForce,'Boussinesq',gotIt) ) THEN

             HeatExpansionCoeff(1:n) = GetReal( Material, &
                 'Heat Expansion Coefficient' )
             
             ReferenceTemperature(1:n) = GetReal( Material, &
                 'Reference Temperature' )

             DO i=1,n
               k = TempPerm(NodeIndexes(i))
               IF ( k > 0 ) THEN
                 IF ( Hydrostatic ) THEN
                   Tdiff = 1 - HeatExpansionCoeff(i) * &
                      (Temperature(k) - ReferenceTemperature(i))

                   IF ( Tdiff <= 0.0D0 ) THEN
                      CALL Warn( 'ManFlowSolve','Zero or negative density.' )
                   END IF
                 ELSE
                   Tdiff = -HeatExpansionCoeff(i) * &
                               (Temperature(k) - ReferenceTemperature(i))
                 END IF
  
                 LoadVector(1,i)   = Gravity(1) * Tdiff
                 LoadVector(2,i)   = Gravity(2) * Tdiff
                 IF ( NSDOFs > 3 ) THEN
                   LoadVector(3,i) = Gravity(3) * Tdiff
                 END IF
               END IF
             END DO
           ELSE IF ( Hydrostatic ) THEN
             LoadVector(1,1:n)   = Gravity(1)
             LoadVector(2,1:n)   = Gravity(2)
             IF ( NSDOFs > 3 ) LoadVector(3,1:n) = Gravity(3)
           END IF
!            WRITE(*,*) "SIF VERIF - LoadVectorX=",LoadVector(1,1:n)
!            WRITE(*,*) "SIF VERIF - LoadVectorY=",LoadVector(2,1:n)
!            WRITE(*,*) "SIF VERIF - LoadVectorZ=",LoadVector(3,1:n)
!------------------------------------------------------------------------------
           LoadVector(1,1:n) = LoadVector(1,1:n) + ListGetReal( BodyForce, &
               'Flow Bodyforce 1',n,NodeIndexes,gotIt )
           
           LoadVector(2,1:n) = LoadVector(2,1:n) + ListGetReal( BodyForce, &
               'Flow Bodyforce 2',n,NodeIndexes,gotIt )
           
           IF ( NSDOFs > 3 ) THEN
             LoadVector(3,1:n) = LoadVector(3,1:n) + ListGetReal( BodyForce, &
                 'Flow Bodyforce 3',n,NodeIndexes,gotIt )
           END IF

!------------------------------------------------------------------------------
           
           PotentialForce = ListGetLogical( BodyForce,'Potential Force',gotIt) 
           IF(PotentialForce) THEN
             PotentialField(1:n) = ListGetReal( BodyForce, &
                 'Potential Field',n,NodeIndexes)             
             PotentialCoefficient(1:n) = ListGetReal( BodyForce, &
                 'Potential Coefficient',n,NodeIndexes)
           END IF
        
!------------------------------------------------------------------------------
         END IF ! of body forces
!------------------------------------------------------------------------------
!
! NOTE: LoadVector is multiplied by density inside *Navier* routines
!
         IF ( Transient ) THEN
           SELECT CASE( CompressibilityModel )
           CASE( PerfectGas1 )
             IF ( ASSOCIATED( TempSol ) ) THEN
               DO i=1,n
                 k = TempPerm(NodeIndexes(i))
                 IF ( k > 0 ) THEN
                    LoadVector(NSDOFs,i) = LoadVector(NSDOFs,i) + &
                      ( Temperature(k) - TempPrev(k) ) / dt
                 END IF
               END DO
             END IF
           CASE( UserDefined1, Thermal )
              DO i=1,n
                LoadVector(NSDOFs,i) = LoadVector(NSDOFs,i) - &
                  ( Density(i) - PrevDensity(i) ) / (Density(i)*dt)
              END DO
           END SELECT
         END IF

!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
! CARTESIAN INCOMPRESSIBLE
!------------------------------------------------------------------------------
!---------------------- Composition de l'opérateur Tangent Lt -----------------
!------------------------- Lt(*)=L(*) + Q(*,U0) + Q(U0,*) ---------------------
!------------------------------------------------------------------------------

         CALL MANTangentOperatorCompose( MASS,STIFF,FORCE, LoadVector,               &
                   Viscosity,Density,U,V,W,MU,MV,MW,ReferencePressure+Pressure(1:n), &
                   LocalTemperature, Convect, StabilizeFlag, CompressibilityModel,   &
                   PseudoCompressible, PseudoCompressibility, GasConstant, Porous,   &
                   Drag, PotentialForce, PotentialField, PotentialCoefficient,       &
                   MagneticForce, Rotating, AngularVelocity, DivDiscretization,      &
                   GradPDiscretization, NewtonLinearization, Element,n,ElementNodes)
                   
!------------------------------------------------------------------------------
!        If time dependent simulation, add mass matrix to global 
!        matrix and global RHS vector
!          write(*,*) 'nb=',nb,Bubbles
!------------------------------------------------------------------------------
!          IF ( CompressibilityModel /= Incompressible .AND. &
!                  StabilizeFlag == 'stabilized' ) THEN
!             Bubbles = .TRUE.
!             StabilizeFlag = 'bubbles'
!             write(*,*) 'aa'
!          END IF
!          IF ( Element % TYPE % BasisFunctionDegree <= 1 .AND. &
!               StabilizeFlag == 'p2/p1' ) THEN
!             Bubbles = .TRUE.
!             StabilizeFlag = 'bubbles'
!             write(*,*) 'bb'            
!          END IF
!          write(*,*) 'nb=',nb,Bubbles
!          IF ( nb==0 .AND. Bubbles ) nb = n
!          write(*,*) 'nb=',nb
!          TimeForce = 0.0_dp
!          IF ( Transient ) THEN
!------------------------------------------------------------------------------
!          NOTE: the following will replace STIFF and FORCE
!          with the combined information
!------------------------------------------------------------------------------
!            write(*,*) 'Default1stOrderTime',t
!            CALL Default1stOrderTime( MASS, STIFF, FORCE )
!          END IF

!          IF ( nb > 0 ) THEN
!             write(*,*) 'NSCondensate',t,nd,nb,NSDOFs-1      
!             CALL NSCondensate( nd, nb, NSDOFs-1, STIFF, FORCE, TimeForce )
!          END IF
! 
! !------------------------------------------------------------------------------
!        Add local stiffness matrix and force vector to global matrix & vector
!------------------------------------------------------------------------------
!           write(*,*) 'DEBUG - STIFF ',STIFF
!           write(*,*) 'DEBUG - FORCE ',FORCE
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------      
       END DO ! Fin boucle sur les éléments du domaine
!------------------------------------------------------------------------------
!  ______ _   _ _____    ______ _                           _     _                       
! |  ____| \ | |  __ \  |  ____| |                         | |   | |                      
! | |__  |  \| | |  | | | |__  | | ___ _ __ ___   ___ _ __ | |_  | |     ___   ___  _ __  
! |  __| | . ` | |  | | |  __| | |/ _ \ '_ ` _ \ / _ \ '_ \| __| | |    / _ \ / _ \| '_ \ 
! | |____| |\  | |__| | | |____| |  __/ | | | | |  __/ | | | |_  | |___| (_) | (_) | |_) |
! |______|_| \_|_____/  |______|_|\___|_| |_| |_|\___|_| |_|\__| |______\___/ \___/| .__/ 
!                                                                                  | |    
!                                                                                  |_|  
!------------------------------------------------------------------------------
       CALL DefaultFinishBulkAssembly()

     
       CALL Info( 'OperatorsLtF', 'Assembly done', Level=4 )
!  ___    _   _  ___       _____                                _      _          
! (  _`\ ( ) ( )(  _`\    (  _  )                              ( )    (_ )        
! | (_(_)| `\| || | ) |   | (_) |  ___   ___    __    ___ ___  | |_    | |  _   _ 
! |  _)_ | , ` || | | )   |  _  |/',__)/',__) /'__`\/' _ ` _ `\| '_`\  | | ( ) ( )
! | (_( )| |`\ || |_) |   | | | |\__, \\__, \(  ___/| ( ) ( ) || |_) ) | | | (_) |
! (____/'(_) (_)(____/'   (_) (_)(____/(____/`\____)(_) (_) (_)(_,__/'(___)`\__, |
!                                                                          ( )_| |
!                                                                          `\___/'
!------------------------------------------------------------------------------

!  _____               _                _____           _ _ _   _             
! | __  |___ _ _ ___ _| |___ ___ _ _   |     |___ ___ _| |_| |_|_|___ ___ ___ 
! | __ -| . | | |   | . | .'|  _| | |  |   --| . |   | . | |  _| | . |   |_ -|
! |_____|___|___|_|_|___|__,|_| |_  |  |_____|___|_|_|___|_|_| |_|___|_|_|___|
!                               |___|                                         
!     Neumann & Newton boundary conditions
!------------------------------------------------------------------------------

       DO t = 1,GetNOFBoundaryElements()

         Element => GetBoundaryElement(t)
         IF ( .NOT. ActiveBoundaryElement() ) CYCLE

         n = GetElementNOFNodes()

!        The element type 101 (point element) can only be used
!        to set Dirichlet BCs, so skip \B4em at this stage.
         IF( .NOT. PossibleFluxElement(Element) ) CYCLE

         CALL GetElementNodes( ElementNodes )
         NodeIndexes => Element % NodeIndexes

         BC => GetBC()
         IF ( .NOT. ASSOCIATED(BC) ) CYCLE

!------------------------------------------------------------------------------
         GotForceBC = GetLogical( BC, 'Flow Force BC',gotIt )
         IF ( .NOT. gotIt ) GotForceBC = .TRUE.

         IF ( GotForceBC ) THEN
           LoadVector  = 0.0d0
           Alpha       = 0.0d0
           ExtPressure = 0.0d0
           Beta        = 0.0d0
           SlipCoeff   = 0.0d0
           STIFF = 0.0d0
           FORCE = 0.0d0

!------------------------------------------------------------------------------
!          (at the moment the following is done...)
!          BC: \tau \cdot n = \alpha n +  @\beta/@t + R_k u_k + F
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          normal force BC: \tau\cdot n = \alpha n
!------------------------------------------------------------------------------
           IF ( GetLogical( BC, 'Free Surface',gotIt) ) THEN
             Alpha(1:n) = GetReal( BC,'Surface Tension Coefficient', gotIt )
           END IF

           ExtPressure(1:n) = GetReal( BC, 'External Pressure', GotForceBC )
           IF(.NOT. GotForceBC) ExtPressure(1:n) = GetReal( BC, 'Normal Pressure', GotForceBC )
!------------------------------------------------------------------------------
!         tangential force BC:
!         \tau\cdot n = @\beta/@t (tangential derivative of something)
!------------------------------------------------------------------------------
              
           IF ( ASSOCIATED( TempSol ) ) THEN
             Beta(1:n) = GetReal( BC, &
                 'Surface Tension Expansion Coefficient',gotIt )

             IF ( gotIt ) THEN
               DO j=1,n
                 k = TempPerm( NodeIndexes(j) )
                 IF ( k>0 ) Beta(j) = 1.0_dp - Beta(j) * Temperature(k)
               END DO
               Beta(1:n) = Beta(1:n) * GetReal(BC, 'Surface Tension Coefficient' )
             ELSE
               Beta(1:n) = GetReal( BC,'Surface Tension Coefficient', gotIt ) 
             END IF
           END IF

!------------------------------------------------------------------------------
!         force in given direction BC: \tau\cdot n = F
!------------------------------------------------------------------------------

           LoadVector(1,1:n) =  GetReal( BC, 'Pressure 1', GotIt )
           LoadVector(2,1:n) =  GetReal( BC, 'Pressure 2', GotIt )
           LoadVector(3,1:n) =  GetReal( BC, 'Pressure 3', GotIt )
           LoadVector(4,1:n) =  GetReal( BC, 'Mass Flux', GotIt )

!------------------------------------------------------------------------------
!          slip boundary condition BC: \tau\cdot n = R_k u_k
!------------------------------------------------------------------------------

           SlipCoeff = 0.0d0
           SlipCoeff(1,1:n) =  GetReal( BC, 'Slip Coefficient 1',GotIt )
           SlipCoeff(2,1:n) =  GetReal( BC, 'Slip Coefficient 2',GotIt )
           SlipCoeff(3,1:n) =  GetReal( BC, 'Slip Coefficient 3',GotIt )

           NormalTangential = GetLogical( BC, &
                 'Normal-Tangential Velocity', GotIt )

           IF (.NOT.GotIt) THEN
             NormalTangential = GetLogical( BC, &
                    'Normal-Tangential '//GetVarName(Solver % Variable), GotIt )
           END IF
!------------------------------------------------------------------------------
! Cartesian
           CALL NavierStokesBoundary(  STIFF, FORCE, &
              LoadVector, Alpha, Beta, ExtPressure, SlipCoeff, NormalTangential,   &
                 Element, n, ElementNodes )
!------------------------------------------------------------------------------

           IF ( GetLogical( BC, 'Wall Law',GotIt ) ) THEN
!
             Density(1:n)   = GetParentMatProp( 'Density' )
             Viscosity(1:n) = GetParentMatProp( 'Viscosity' )
!
             LayerThickness(1:n) = GetReal( BC, 'Boundary Layer Thickness' )
             SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness',GotIt )
!
             SELECT CASE( NSDOFs )
               CASE(3)
                 U(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-2 )
                 V(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-1 )
                 W(1:n) = 0.0d0
             
               CASE(4)
                 U(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-3 )
                 V(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-2 )
                 W(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-1 )
             END SELECT
!
             CALL NavierStokesWallLaw( STIFF,FORCE,     &
               LayerThickness,SurfaceRoughness,Viscosity,Density,U,V,W, &
                      Element,n, ElementNodes )
!
           ELSE IF ( GetLogical( BC, 'VMS Wall', GotIt ) ) THEN
             Density(1:n)   = GetParentMatProp( 'Density' )
             Viscosity(1:n) = GetParentMatProp( 'Viscosity' )
!
             LayerThickness(1:n) = GetReal( BC, 'Boundary Layer Thickness', GotIt )
             SurfaceRoughness(1:n) = GetReal( BC, 'Surface Roughness',GotIt )
!
             SELECT CASE( NSDOFs )
               CASE(3)
                 U(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-2 )
                 V(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-1 )
                 W(1:n) = 0.0d0            
               CASE(4)
                 U(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-3 )
                 V(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-2 )
                 W(1:n) = USOL( NSDOFs*FlowPerm(NodeIndexes)-1 )
             END SELECT
             CALL VMSWalls( STIFF,FORCE,     &
                            LayerThickness,SurfaceRoughness,Viscosity,Density,U,V,W, &
                            Element,n, ElementNodes )

           END IF
!------------------------------------------------------------------------------
           IF ( Transient ) THEN
             MASS = 0.0d0
             CALL Default1stOrderTime( MASS, STIFF, FORCE )
           END IF
!------------------------------------------------------------------------------
!         Add local stiffness matrix and force vector to
!         global matrix & vector
!--------------------------------------------------------------------------
           CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
         END IF
       END DO
!------------------------------------------------------------------------------
       !
       ! IMPLEMENT NOSLIP WALL BC CODE:
       ! ------------------------------
       DO i=1,Model % NumberOFBCs
         BC => Model % BCs(i) % Values
         IF ( GetLogical(  BC, 'Noslip wall BC', gotit ) ) THEN
           IF ( VarName  == 'flow solution' ) THEN
             CALL ListAddConstReal( BC, 'Velocity 1', 0.0_dp )
             CALL ListAddConstReal( BC, 'Velocity 2', 0.0_dp )
             IF ( NSDOFs>3 ) CALL ListAddConstReal( BC, 'Velocity 3', 0.0_dp )
           ELSE
             DO j=1,NSDOFs-1
               CALL ListAddConstReal( BC, ComponentName( &
                  Solver % Variable % Name, j), 0.0_dp )
             END DO
           END IF
         END IF
       END DO

       CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
!# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
!TEST 20150317 - Remove Dirichlet from here, and put it outside to check residual without Diric on Stiff   
!        CALL DefaultDirichletBCs()
!# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

       
!!! POSSIBLE DE NE PAS REMETTRE LE PROFIL A CHAQUES FOIS DANS FlowSOLUTION
!!! IL FAUT UN OBJET Variable QUI CONTIENDRAI le FlowInit!
! ! !        CALL DefaultDirichletBCs(Solver,FlowSolution_Init)

!        CALL Info( 'OperatorsLtF', 'Dirichlet conditions done', Level=4 )
!  _____ _____ ____     _____               _                _____           _ _ _   _             
! |   __|   | |    \   | __  |___ _ _ ___ _| |___ ___ _ _   |     |___ ___ _| |_| |_|_|___ ___ ___ 
! |   __| | | |  |  |  | __ -| . | | |   | . | .'|  _| | |  |   --| . |   | . | |  _| | . |   |_ -|
! |_____|_|___|____/   |_____|___|___|_|_|___|__,|_| |_  |  |_____|___|_|_|___|_|_| |_|___|_|_|___|
!                                                    |___|              
!------------------------------------------------------------------------------

    END SUBROUTINE OperatorsLtF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! - - - - - - - - V2 STABILIZED Franca Test
! TANGENT OPERATOR
   SUBROUTINE MANTangentOperatorCompose  (                                            &
       MassMatrix, StiffMatrix, ForceVector, LoadVector, Nodalmu,        &
       Nodalrho, Ux, Uy, Uz, MUx, MUy, MUz, NodalPressure, NodalTemperature,&
       Convect, StabilizeFlag, Cmodel, &
       PseudoCompressible, NodalCompressibility, NodalGasC, Porous, NodalDrag,  &
       PotentialForce, PotentialField, PotentialCoefficient, MagneticForce,     &
       Rotating, Omega, &
       DivDiscretization,  gradPDiscretization, NewtonLinearization,            &
       Element, n, Nodes )

  USE DefUtils
  USE Differentials
  USE MaterialModels
        IMPLICIT NONE
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: MassMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL(KIND=dp) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT:
!
!  REAL(KIND=dp) :: Nodalmu(:)
!     INPUT: Nodal values for viscosity (i.e. if turbulence model or
!             power-law viscosity is used, the values vary in space)
!
!  REAL(KIND=dp) :: Nodalrho(:)
!     INPUT: Nodal values of density
!
!  REAL(KIND=dp) :: Ux(:),Uy(:),Uz(:)
!     INPUT: Nodal values of velocity components from previous iteration
!
!  REAL(KIND=dp) :: NodalPressure(:)
!     INPUT: Nodal values of total pressure from previous iteration
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilization be used ?
!
!  LOGICAL :: Compressible1, Compressible2
!     INPUT: Should compressible flow terms be added ?
!
!  LOGICAL :: PseudoCompressible
!     INPUT: Should artificial compressibility be added ?
!
!  REAL(KIND=dp) :: NodalCompressibility(:)
!     INPUT: Artificial compressibility for the nodes
!
!  LOGICAL :: MagneticForce
!      INPUT: Should Lorentz force for magnetohydrodynamics be included
!
!  LOGICAL :: Rotating
!      INPUT: Is the coordinate system rotating
!  
!  REAL(KIND=dp) :: Omega(:)
!      INPUT: If previous is True, components of angular velocity
!
!  LOGICAL :: NewtonLinearization
!      INPUT: Picard or Newton linearization of the convetion term ?
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

     REAL(KIND=dp),TARGET :: MassMatrix(:,:),StiffMatrix(:,:),ForceVector(:)
     REAL(KIND=dp), DIMENSION(:) :: Ux,Uy,Uz,MUx,MUy,MUz,Omega
     LOGICAL :: PseudoCompressible, Porous, &
         NewtonLinearization, Convect, DivDiscretization, gradPDiscretization, &
         PotentialForce, MagneticForce, Rotating
     CHARACTER(LEN=*) :: StabilizeFlag
     REAL(KIND=dp) :: Nodalmu(:),Nodalrho(:), &
       NodalPressure(:), LoadVector(:,:), NodalTemperature(:), &
       NodalCompressibility(:), NodalDrag(:,:), PotentialField(:), &
       PotentialCoefficient(:), NodalGasC(:)

     INTEGER :: n, Cmodel

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

     REAL(KIND=dp) :: LES(6,6,n)
     TYPE(Variable_t), POINTER :: LESVar
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp) :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)

     REAL(KIND=dp) :: Velo(3),UVelo(3),Grad(3,3),Force(4),Metric(3,3),Symb(3,3,3), Drag(3), Coord(3)

     REAL(KIND=dp), POINTER :: A(:,:),M(:,:),Load(:), Jac(:,:)
     REAL(KIND=dp) :: SU(n,4,4),SW(n,4,4),LrF(3), LSx(3,3), VC(6,6), B(6,3), G(3,6)

     REAL(KIND=dp) :: Lambda=1.0,Re,Tau,Delta,Re1,Re2,GasC, hScale

     REAL(KIND=dp) :: VNorm,hK,mK,mu,dmudx(3),Temperature, &
               drhodp,drhodp_n(n)

     REAL(KIND=dp) :: drhodx(3), rho, Pressure, dTemperaturedx(3), &
                      dPressuredx(3),dPrevPressuredx(3), Compress, masscoeff

     REAL(KIND=dp), TARGET :: JacM(8*n,8*n), SOL(8*n)

     REAL(KIND=dp) :: Tau_M, Tau_C, Gmat(3,3),Gvec(3), dt=0._dp, C1, NodalVelo(4,n), &
       RM(n,4,4),LC(3,n), gradP(n), PRM(4), GradNodal(n,4,3), PVelo(3), NodalPVelo(4,n), &
       Work(3,n), maxl, minl, Strain(3,3), muder0, muder, ViscConstantCondition

     INTEGER :: i,j,k,l,c,p,q,t,dim

     REAL(KIND=dp) :: s,u,v,w,volume

     INTEGER, POINTER :: EdgeMap(:,:)

     TYPE(ElementType_t), POINTER :: LinearType, SaveType
     INTEGER :: LinearBasis, LinearCode(3:8) = (/ 303,404,504,605,706,808 /)
  
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis, deg(100), Order
     REAL(KIND=dp), POINTER :: NodalC(:,:,:)
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, Transient, stat, Bubbles, PBubbles, Stabilize, &
                VMS, P2P1, Isotropic, drhodp_found, Compressible, ViscNewtonLin, &
                ViscNonnewtonian, LaplaceDiscretization

!------------------------------------------------------------------------------
!20150414 : picard makes pressure looks good with ANM, but Newton doesn't...?
     NewtonLinearization = .TRUE.
     Convect = .TRUE.
!------------------------------------------------------------------------------     

     dim = CoordinateSystemDimension()
     Params => GetSolverParams()


     hScale = GetCReal( Params, 'H scale', Found )
     IF ( .NOT. Found )  hScale = 1._dp

     c = dim + 1

     ForceVector = 0._dp
     MassMatrix  = 0._dp
     StiffMatrix = 0._dp
!     JacM = 1._dp
     IF ( NewtonLinearization ) JacM = 0._dp
     
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis = n
     Bubbles = .FALSE.

     Material => GetMaterial()
     Isotropic=GetString( Material, 'Viscosity Model',ViscNonnewtonian )/='anisotropic'
!      write(*,*) 'isotropic=',isotropic
!      write(*,*) 'ViscNonnewtonian=',ViscNonnewtonian
     IF( ViscNonnewtonian ) THEN
       ViscConstantCondition = GetCReal( Params, 'Newtonian Viscosity Condition',Found)
!       write(*,*) 'ViscConstantCondition =',ViscConstantCondition
       IF( Found .AND. ViscConstantCondition > 0.0_dp ) ViscNonnewtonian = .FALSE.
     END IF

     LaplaceDiscretization = GetLogical( Params,'Laplace Discretization', Found)

     Transient = GetString(GetSimulation(),'Simulation Type',Found)=='transient'

     Compressible = Cmodel /= Incompressible
     drhodp_found = .FALSE.
     drhodp_n(1:n) = GetReal( Material, 'drho/dp', drhodp_found )

! - - - - - - - - V2 STABILIZED Franca Test
     Stabilize = StabilizeFlag == 'stabilized'
! - - - - - - - - V2 STABILIZED Franca Test

     LinearBasis = 0


     
     IntegStuff = GaussPoints( Element )
     
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n


     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!------------------------------------------------------------------------------
     
     hK = element % hK*hscale
     mK = element % StabilizationMK

     IF ( Stabilize ) THEN
       dNodalBasisdx = 0._dp
       DO p=1,n
         u = Element % TYPE % NodeU(p)
         v = Element % TYPE % NodeV(p)
         w = Element % TYPE % NodeW(p)
         stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
         dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
       END DO
     ENDIF
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
     

     ViscNewtonLin = .FALSE.
     
!      ViscNewtonLin = .TRUE.
!      ViscNonnewtonian = .FALSE.
!      ViscNonnewtonian = .TRUE.
      
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
    DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
      
      
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------

      stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
              Basis, dBasisdx, Bubbles=Bubbles )

      s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!     rho at the integration point
!------------------------------------------------------------------------------
      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
      dPressuredx = 0._dp
      dTemperaturedx = 0._dp
      
!------------------------------------------------------------------------------
!     Velocity from previous iteration (relative to mesh velocity)
!     at the integration point
!------------------------------------------------------------------------------
!     U0
      Velo = 0.0_dp
      Velo(1) = SUM( (Ux(1:n)-MUx(1:n))*Basis(1:n) )
      Velo(2) = SUM( (Uy(1:n)-MUy(1:n))*Basis(1:n) )
      IF ( DIM > 2 ) Velo(3) = SUM( (Uz(1:n)-MUz(1:n))*Basis(1:n) )

!     Grad(U0) 
      Grad = 0.0_dp
      DO i=1,3
        Grad(1,i) = SUM( Ux(1:n) * dBasisdx(1:n,i) )
        Grad(2,i) = SUM( Uy(1:n) * dBasisdx(1:n,i) )
        IF ( DIM > 2 ) Grad(3,i) = SUM( Uz(1:n) * dBasisdx(1:n,i) )
      END DO
!------------------------------------------------------------------------------
!     Force at integration point
!------------------------------------------------------------------------------
!     LoadVector 
      Force = 0.0_dp
      DO i=1,c
        Force(i) = SUM( LoadVector(i,1:n) * Basis(1:n) )
!        write(*,*) 'Force(',i,')=',Force(i)
      END DO
      
! - - - -  - - - -  - - - -  - - - -  - - - -       
      IF ( MagneticForce ) THEN
         LrF = LorentzForce( Element,Nodes,u,v,w,n )
         Force(1:DIM) = Force(1:DIM) + Lrf(1:DIM) / rho
      END IF
! - - - -  - - - -  - - - -  - - - -  - - - -       

!------------------------------------------------------------------------------
!     Additional forces due to gradient forces (electrokinetic flow) and
!     viscous drag in porous media.
!------------------------------------------------------------------------------

!       IF ( Convect .AND. NewtonLinearization ) THEN
!      write(*,*) 'Convect .AND. NewtonLinearization=', Convect .AND. NewtonLinearization
!         Uvelo = 0.0_dp
!         Uvelo(1) = SUM( Basis(1:n) * Ux(1:n) )
!         Uvelo(2) = SUM( Basis(1:n) * Uy(1:n) )
!         IF ( DIM > 2 ) Uvelo(3) = SUM( Basis(1:n) * Uz(1:n) )

        
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 20150324  - Test Fv = F OK
! 20150324  - Test Fv = F -L(U0) - Q(U0,U0) pas OK car Contour à 1, mais interieur >1
! Profile vitesse quasi ok... mais countour non ok...
!------------------------------------------------------------------------------
! ! ! ! ! AVANT : IF CONVECT AND N LINEAR Q(U0,U0)
! ! ! ! !         DO i=1,dim
! ! ! ! !           DO j=1,dim
! ! ! ! !             Force(i) = Force(i) + Grad(i,j) * Uvelo(j) 
! ! ! ! !           END DO
! ! ! ! !         END DO
!       U0xGrad.U0 ajouté au second membre (c'est le mal !) ===> tout seul mais avec LtU0 ca passe mieux
! 20150320 - REMODIF 
! ! !        DO i=1,dim
! ! !          DO j=1,dim
! ! ! ! 20150324
! ! ! ! Fv = F - ( L(U0) + Q(U0,U0) )
! ! ! ! -( L(U0) + Q(U0,U0) ) = -( StiffMatrix X U0            - Q(U0,U0)  )
! ! ! !                       = -( L(U0) + Q(U0,U0) + Q(U0,U0) - Q(U0,U0)   )
! ! ! !                       = - MATMUL(STIFF, U0) + Q(U0,U0)
! ! ! ! test Lin classique (marche pas avec ANM) TOUT SEUL, TEST EN COURS POUR Fv EVE LIKE
! ! !            Force(i) = Force(i) + Uvelo(j) * Grad(i,j)  ! Inclure avec MATMUL STIFF SOL!!!
! ! ! !            Force(i) = Force(i) - Uvelo(j) * Grad(i,j)    ! test TOUT SEUL
! ! !          END DO
! ! !        END DO
! 20150320 - REMODIF        
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!       END IF
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!     Effective viscosity & derivatives at integration point
!------------------------------------------------------------------------------
!       write(*,*) 'avant mu, ViscNewtonLin=',ViscNewtonLin
      IF ( isotropic ) THEN
        mu = SUM( Nodalmu(1:n) * Basis(1:n) )
!       write(*,*) 'isotropic mu =',mu

        IF( ViscNonnewtonian ) THEN 
        ! ViscNonnewtonian = .false. --> fluide non newtionien

          mu = EffectiveViscosity( mu, rho, Ux, Uy, Uz, &
              Element, Nodes, n, n, u, v, w,  muder0 )
          
          ViscNewtonLin = NewtonLinearization .AND. muder0/=0._dp
          IF ( ViscNewtonLin )  Strain = (Grad+TRANSPOSE(Grad))/2

        END IF
        

      ELSE
        DO i=1,6
        DO j=1,6
          VC(i,j) = SUM(NodalC(i,j,1:n)*Basis(1:n))
        END DO
        END DO
      END IF
!       write(*,*) 'apres avoir recupere mu, ViscNewtonLin=',ViscNewtonLin


! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
      IF ( Stabilize ) THEN
        DO i=1,3
          dmudx(i) = SUM( Nodalmu(1:n)*dBasisdx(1:n,i) )
        END DO
!------------------------------------------------------------------------------
!       Stabilization parameters Tau & Delta
!------------------------------------------------------------------------------
        IF ( Convect ) THEN
           VNorm = MAX( SQRT( SUM(Velo(1:DIM)**2) ), 1.0d-12 )
           Re = MIN( 1.0d0, rho * mK * hK * VNorm / (4 * mu) )
 
           Tau = hK * Re / (2 * rho * VNorm)
           Delta = rho * Lambda * Re * hK * VNorm
        ELSE
           Delta = 0._dp
           Tau   = mK * hK**2  / ( 8 * mu )
        END IF
!------------------------------------------------------------------------------
!       SU will contain residual of ns-equations (except for the time derivative
!       and force terms). SW will contain the weight function values.
!------------------------------------------------------------------------------
        SU(1:n,:,:) = 0.0D0
        SW(1:n,:,:) = 0.0D0

        DO p=1,N
          DO i=1,DIM
            SU(p,i,c) = SU(p,i,c) + dBasisdx(p,i)
            IF(Porous) THEN
              SU(p,i,i) = SU(p,i,i) + mu * Drag(i) * Basis(p)
            END IF

            IF ( Convect ) THEN
              DO j=1,DIM
                SU(p,i,i) = SU(p,i,i) + rho * dBasisdx(p,j) * Velo(j)
              END DO
            END IF

            DO j=1,DIM
              SU(p,i,i) = SU(p,i,i) - dmudx(j) * dBasisdx(p,j)
              SU(p,i,j) = SU(p,i,j) - dmudx(j) * dBasisdx(p,i)
              SU(p,i,i) = SU(p,i,i) - mu * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
              SU(p,i,j) = SU(p,i,j) - mu * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
            END DO

            IF ( Convect .AND. NewtonLinearization ) THEN
              DO j=1,DIM
                SU(p,i,j) = SU(p,i,j) + rho * Grad(i,j) * Basis(p)
              END DO
            END IF
!
!------------------------------------------------------------------------------

            IF ( Convect ) THEN
              SW(p,c,i) = SW(p,c,i) + rho * dBasisdx(p,i)
              DO j=1,dim
                SW(p,i,i) = SW(p,i,i) + rho * dBasisdx(p,j) * Velo(j)
              END DO
            ELSE
              SW(p,c,i) = SW(p,c,i) + dBasisdx(p,i)
            END IF

            DO j=1,dim
              SW(p,i,i) = SW(p,i,i) - dmudx(j) * dBasisdx(p,j)
              SW(p,j,i) = SW(p,j,i) - dmudx(j) * dBasisdx(p,i)
              SW(p,i,i) = SW(p,i,i) - mu * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
              SW(p,j,i) = SW(p,j,i) - mu * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
            END DO
          END DO
        END DO
      ENDIF
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
      
      
      
!#ifdef LES
!#define LSDelta (SQRT(2.0d0)*Element % hK)
!#define LSGamma 6

!------------------------------------------------------------------------------
!    Loop over basis functions (of both unknowns and weights)
!------------------------------------------------------------------------------
     DO p=1,NBasis

     BaseP =  Basis(p)
 
     DO q=1,NBasis
         i = c*(p-1)
         j = c*(q-1)
         M => MassMatrix ( i+1:i+c,j+1:j+c )
         A => StiffMatrix( i+1:i+c,j+1:j+c )
         IF ( ViscNewtonLin ) Jac => JacM( i+1:i+c,j+1:j+c )
!------------------------------------------------------------------------------
!      First plain Navier-Stokes
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!      Mass matrix:
!------------------------------------------------------------------------------
       ! Momentum equations
       DO i=1,dim
         M(i,i) = M(i,i) + s*rho*Basis(q)*Basis(p)
       END DO
!------------------------------------------------------------------------------
!      Stiffness matrix:
!------------------------------

!------------------------------------------------------------------------------
!      Diffusive terms
!      Convection terms, Picard linearization
!------------------------------------------------------------------------------
!write(*,*) 'juste avant Picard ViscNewtonLin=',ViscNewtonLin
       IF ( ViscNewtonLin ) THEN
         DO i=1,dim

           muder = muder0 * 4 * SUM( Strain(i,:) * dBasisdx(q,:) )
           !do j=1,size(muder)

           !end do
!           muder = 1*4*SUM(Strain(i,:)*dBasisdx(q,:))
           DO j=1,dim
             Jac(j,i)=Jac(j,i) + s * 2 * muder * SUM( Strain(j,:) * dBasisdx(p,:) )
!             write(*,*) 'Jac(',j,',',i,')=', Jac(J,i)
           END DO
         END DO
       END IF


       !FlowSolver:
! Add new Navier-Stokes discretization option 'Laplace Discretization'. 
! The diffusion term discretized is \mu div(grad(u)) instead of the default 2\mu div(\epsilon(u)). 
! The resulting linear system is smaller (in terms of nonzeros) and usually easier to solve. 
! This has, however, an effect on the meaning of the various boundary conditions.
       DO i=1,dim
         DO j = 1,dim
           IF ( Isotropic ) THEN
             A(i,i) = A(i,i) + s * mu * dBasisdx(q,j) * dBasisdx(p,j)     ! Laplacien  Grad u : Grad v
             IF ( divDiscretization ) THEN             
                A(i,j) = A(i,j) + s * mu * dBasisdx(q,j) * dBasisdx(p,i)
             ELSE IF (.NOT.LaplaceDiscretization) THEN             
                A(i,j) = A(i,j) + s * mu * dBasisdx(q,i) * dBasisdx(p,j)  ! Grad u T : Grad v T
             END IF

!------------------------------------------------------------------------------
!  For compressible flow add grad((2/3) \mu div(u))
!------------------------------------------------------------------------------
             IF ( Compressible ) THEN
               A(i,j) = A(i,j) - s * ( 2._dp / 3._dp ) * mu * &
                          dBasisdx(q,j) * dBasisdx(p,i)
             END IF
           END IF ! boucle IF ( Istropic )


!    Lt(*)=L(*) + Q(*,U0) + Q(U0,*)
! ICI on connait Uvelo et on plaque dbasisdx <-> quasi grad(*)
! on a donc Q(U0,*)
           IF ( Convect ) THEN
              A(i,i) = A(i,i) + s * rho * dBasisdx(q,j) * Velo(j) * Basis(p)
           END IF
         END DO ! Fin de la boucle j
!------------------------------------------------------------------------------
 
         ! Pressure terms:
         ! --------------- 
         IF ( gradPDiscretization  ) THEN
            A(i,c) = A(i,c) + s * dBasisdx(q,i) * Basis(p)
         ELSE
            A(i,c) = A(i,c) - s * Basis(q) * dBasisdx(p,i)
         END IF


         ! Continuity equation:
         !---------------------
         IF ( gradPDiscretization ) THEN
           A(c,i) = A(c,i) - s * rho * Basis(q) * dBasisdx(p,i)
         ELSE
           IF ( Compressible .OR. Convect ) THEN
             A(c,i) = A(c,i) + s * rho * dBasisdx(q,i) * BaseP
           ELSE
             A(c,i) = A(c,i) + s * dBasisdx(q,i) * BaseP
           END IF
         END IF ! IF ( gradPDiscretization )
       END DO ! Boucle sur i = 1 , dim ( construction de A )

!------------------------------------------------------------------------------
!      Convection, Newton linearization
! ANM : PICARD MAKES PRESSURE LOOKS GOOD, BUT NEWTON PUT STRANGE PRESSURE INLET ONLY...
!------------------------------------------------------------------------------
!    Lt(*)=L(*) + Q(U0,*) + Q(*,U0)
! Ici on a Grad est calculé en U0 donc on a Q(*,U0) == * fois gradU0
       IF ( Convect .AND. NewtonLinearization ) THEN
         DO i=1,dim
           DO j=1,dim
              A(i,j) = A(i,j) + s * rho * Basis(q) * Grad(i,j) *  Basis(p)
           END DO
         END DO
       END IF
       
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test         
!------------------------------------------------------------------------------
!      Add stabilization...
!------------------------------------------------------------------------------
       IF ( Stabilize ) THEN 
          DO i=1,dim
             DO j=1,c
               M(j,i) = M(j,i) + s * Tau * rho * Basis(q) * SW(p,j,i)
             END DO

             DO j=1,dim
               A(j,i) = A(j,i) + s * Delta * dBasisdx(q,i) * dBasisdx(p,j)
             END DO
          END DO
!------------------------------------------------------------------------------          
!------------------------------------------------------------------------------          
          A = A + s * Tau * MATMUL(SW(p,1:c,1:dim),SU(q,1:dim,1:c))
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
          
       ENDIF
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test         
       
       
     END DO ! Boucle sur q=1 , Nbasis
     END DO ! Boucle sur p=1 , Nbasis

!------------------------------------------------------------------------------
!    The righthand side...
!------------------------------------------------------------------------------


     DO p=1,NBasis
       Load => ForceVector( c*(p-1)+1 : c*(p-1)+c )

       DO i=1,c
         Load(i) = Load(i) + s * rho * Force(i) * Basis(p)
!         write(*,*) 'load(',i,')=',load(i)
       END DO
       
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test                
       IF ( Stabilize ) THEN
         DO i=1,DIM
           DO j=1,c
             Load(j) = Load(j) + s * Tau * rho * Force(i) * SW(p,j,i)
           END DO
         END DO
       ENDIF
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test         
       
       
      
     END DO ! Fin boucle p=1,NBasis
   END DO ! Fin boucle t=1, N_Integ

! ! ! ! !    L(U) + Q(U,U) = F
! ! ! ! !    U = U0 + V
! ! ! ! ! StiffMatrix = L(V) + Q(U0,V) + Q(V,U0)
! ! ! ! !    L(V) + Q(U0,V) + Q(V,U0)  = Lt assemblé en U0
! ! ! ! ! 20150324
! ! ! ! ! L(U0) + Q(U0,U0) = StiffMatrix X U0 - Q(U0,U0)  = L(U0) + Q(U0,U0) + Q(U0,U0) - Q(U0,U0)
! ! ! ! !    Fv = F - L(U0) - Q(U0,U0)
     SOL=0._dp
     SOL(1:c*n:c) = Ux(1:n)
     SOL(2:c*n:c) = Uy(1:n)
     p = c*NBasis
     IF ( dim>2 ) SOL(3:c*n:c) = Uz(1:n)
!      ForceVector(1:p) = ForceVector(1:p) - MATMUL(StiffMatrix(1:p,1:p),SOL(1:p))

! TO SOLVE:
!    Lt ( V ) + Q(V,V) = Fv
!  V = sum a^i * W_i   MAN avec un lambda_i
   
!open(unit=255,file='forcevectorcompose_avantmaj.dat',status='unknown',form='formatted')
!    do i=1,size(forcevector)
!      write(255,*) 'ForceVector(',i,')=',ForceVector(i)
!    end do
!      write(*,*) "COUCOU avViscNewtonLin=",ViscNewtonLin
! ! ! ! !    IF ( ViscNewtonLin ) THEN
! ! ! ! ! !      write(*,*) "COUCOU apViscNewtonLin=",ViscNewtonLin
! ! ! ! ! !!!!!! ON N'EST JAMAIS ICI
! ! ! ! !      SOL=0._dp
! ! ! ! !      SOL(1:c*n:c) = Ux(1:n)
! ! ! ! !      SOL(2:c*n:c) = Uy(1:n)
! ! ! ! !      IF ( dim>2 ) SOL(3:c*n:c) = Uz(1:n)
! ! ! ! !      p = c*nBasis
! ! ! ! !      StiffMatrix(1:p,1:p) = StiffMatrix(1:p,1:p) + JacM(1:p,1:p)
! ! ! ! !      ForceVector(1:p) = ForceVector(1:p) + MATMUL(JacM(1:p,1:p),SOL(1:p))
! ! ! ! !    END IF



!open(unit=87,file='forcevectorcompose.dat',status='unknown',form='formatted')
!    do i=1,size(forcevector)
!      write(87,*) 'ForceVector(',i,')=',ForceVector(i)
!    end do
!------------------------------------------------------------------------------
 END SUBROUTINE MANTangentOperatorCompose

 
 
 
  
!--
! RHS


!                       __      __ 
!         /\ |\ ||\/|  |__)|__|(_  
!        /--\| \||  |  | \ |  |__) 
! FQ classique = Somme des Q croisés
  SUBROUTINE  ANMPFrhsFQ( FQMan, UMan, IO, NSDOFs, FlowPerm, &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez, &
                         Density,Material,FlowSolution_Init)
!  

    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
!

    IMPLICIT NONE     
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: FQMan(:),FlowSolution_Init(:)
    REAL(KIND=dp) :: UMan(:,:)    
    REAL(KIND=dp) :: USAV(:,:,:,:)
    REAL(KIND=dp) :: GradSAV( :,:,:,:,: )   
    INTEGER       :: NSDOFs, IO
    REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:)     
    TYPE(ValueList_t),POINTER :: Material    
    REAL(KIND=dp) :: Density(:)
    REAL(KIND=dp) :: FQManTemp(:)
    
    INTEGER, POINTER :: FlowPerm(:)
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    INTEGER :: t,n,nb,nd
    REAL(KIND=dp), POINTER :: FQlocal(:)   
    
    
    INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
                            
!------------------------------------------------------------------------------

    FQManTemp = 0.0_dp
    DO t = 1 , GetNOFActive()  ! Boucle ELEMENT
      CALL AdvanceOutput( t,GetNOFActive() )      
      Element => GetActiveElement(t)      
!       NodeIndexes => Element % NodeIndexes
      Indexes => Element % NodeIndexes
      
      n = GetElementNOFNodes(Element)      
      nb = GetElementNOFBDOFs(Element)       
      nd = GetElementDOFs( Indexes )
      CALL GetElementNodes( ElementNodes )
!
      Density(1:n) = GetReal( Material, 'Density' )

!------------------------------------------------------------------------------
!     Vitesse au point d'intégration de l'itération précédente
!------------------------------------------------------------------------------
      SELECT CASE( NSDOFs )
        CASE(3) ! 2D
          Uelex(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
          Ueley(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
          Uelez(1:n) = 0.0_dp
        CASE(4) ! 3D
          Uelex(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-3 , IO-1)
          Ueley(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
          Uelez(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
      END SELECT
!------------------------------------------------------------------------------  
      CALL MANFQelemOPTI( FQManTemp, FQlocal, Density,                     &
                          Element, n, IO, ElementNodes,                    &
                          USAV, GradSAV, Element % ElementIndex , NSDOFs,  &
                          Uelex, Ueley, Uelez )
!------------------------------------------------------------------------------                                       

      CALL UpdateGlobalForce( FQMan, FQManTemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
    END DO ! FIN Boucle ELEMENT
!------------------------------------------------------------------------------
!       Conditions de Diriclet modifiées pour le 2nd membre----> 0           
    CALL CalCondLimCoeff( FQMan, 0.0_dp, 1,FlowSolution_Init )
!                           
 END SUBROUTINE ANMPFrhsFQ



 SUBROUTINE  MANFQelemOPTI( FQManTemp, FQlocal, Nodalrho,     &
                            Element, n, io, ElementNodes,     &
                            USAV , GradSAV, numel , NSDOFs,   &
                            Uelex, Ueley, Uelez )
!
  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE Adaptive
  USE SolverUtils

        IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: FQManTemp(:) 
     REAL(KIND=dp), POINTER :: FQlocal(:)
     REAL(KIND=dp) :: Nodalrho(:)
  
     ! Composante de la vitesse MAN ordre r et p-r
     INTEGER :: io, n   ! p = ordremanencours
     ! OMan ordre de la série MAN, nmax le nombre max de degré de liberté
     INTEGER, POINTER :: Flowperm(:)
      
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: Element
!----     
     REAL(KIND=dp) :: USAV(:,:,:,:)
     REAL(KIND=dp) :: GradSAV( :,:,:,:,: )
     INTEGER       :: numel,NSDOFs
     REAL(KIND=dp) ,ALLOCATABLE, TARGET :: UMan(:,:)     
     REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:)     
!----   

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp) :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)


     REAL(KIND=dp) :: Ur(3), GradUpmr(3,3), FQelemtemp(4)

     INTEGER :: i,j,k,l,p,c,q,t,r,dim

     REAL(KIND=dp) :: s,u,v,w,volume
     REAL(KIND=dp) :: rho
    
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, stat
!------------------------------------------------------------------------------
     
     REAL(KIND=dp) :: tmp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!              write(*,*) 'io=',io
     dim = CoordinateSystemDimension()
     Params => GetSolverParams()

     c = dim + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!     
     FQManTemp = 0.0_dp 
!!!!!!!!!!!!!!!!!!!!!!!!!!     
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis = n

     IntegStuff = GaussPoints( Element )

     U_Integ => IntegStuff % u   
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------    
!              write(*,*) 'io=',io
!******************************************************************************
    DO t=1,N_Integ                        ! boucle sur les points d'intégration
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
              Basis, dBasisdx )
!------------------------------------------------------------------------------
      s = detJ * S_Integ(t)

      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
! 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! On sauvegarde U(element,gausspoint)  et Grad U(element,gausspoint) à l'ordre en cours
! en procédant ainsi on ne recalcule pas à chaques ordre, des termes déjà calculé aux ordres précédents.
!------------------------------------------------------------------------------
      USAV( numel , t , io-1 , 1 ) = SUM( Basis(1:n) * Uelex(1:n) )
      USAV( numel , t , io-1 , 2 ) = SUM( Basis(1:n) * Ueley(1:n) )
      IF (dim > 2) USAV( numel , t , io-1 , 3) = SUM( Basis(1:n) * Uelez(1:n) )
!------------------------------------------------------------------------------
! Puis le gradient
! Probleme : comment stocker : element, gausspoint, ordre, composante1, composante2
      DO j=1,3
         GradSAV( numel , t , io-1 , 1 , j ) = SUM( Uelex(1:n) * dBasisdx(1:n,j) )
         GradSAV( numel , t , io-1 , 2 , j ) = SUM( Ueley(1:n) * dBasisdx(1:n,j) )
         IF ( DIM > 2 ) GradSAV( numel , t , io-1 , 3 , j ) = SUM( Uelez(1:n) * dBasisdx(1:n,j) )
      END DO
!------------------------------------------------------------------------------
! Somme des termes croisées par composante spatiale
!------------------------------------------------------------------------------
      FQelemtemp = 0.0_dp
      DO r=1, io-1                  ! Boucle CROSS
       DO i = 1 ,dim
        DO j = 1 ,dim
          FQelemtemp(i) = FQelemtemp(i) - USAV( numel , t , r , j ) * GradSAV( numel , t , io-r , i , j )
        END DO
       END DO
      END DO ! FIN Boucle CROSS
      
!------------------------------------------------------------------------------
! Assemblage
      DO p = 1 , NBasis !boucle sur les fonctions de formes des fonctions tests Wi = Ni
        FQlocal => FQManTemp( c*(p-1)+1 : c*(p-1)+c  )
        DO i = 1, c
          FQlocal(i) = FQlocal(i) + s * rho * FQelemtemp(i) * Basis(p)
        END DO
      END DO ! fin boucle sur les noeuds

!------------------------------------------------------------------------------
   END DO                       ! fin de la boucle sur les points d'intégration
!******************************************************************************   

 END SUBROUTINE MANFQelemOPTI

!  __     __                  __      __ 
! |_ |\ ||  \   /\ |\ ||\/|  |__)|__|(_  
! |__| \||__/  /--\| \||  |  | \ |  |__)                                         
!------------------------------------------------------------------------------    
 
 
 
 
 
 
!------------------------------------------------------------------------------    
!------------------------------------------------------------------------------    
!  ___ ___    _   _  _  ___ _  _ 
! | _ ) _ \  /_\ | \| |/ __| || |
! | _ \   / / _ \| .` | (__| __ |
! |___/_|_\/_/ \_\_|\_|\___|_||_|
!                       __      __ 
!         /\ |\ ||\/|  |__)|__|(_  
!        /--\| \||  |  | \ |  |__) 
! FQ  = Somme des Q croisés ordre > 1 pour Calcul Branches
! Sum j=2, p-1 Q (Uj, U p+1-j )
  SUBROUTINE  ANMBranchrhsFQ( FQMan, UMan, IO, NSDOFs, FlowPerm, &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez, &
                         Density,Material,FlowSolution_Init)
!  

    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
!
    IMPLICIT NONE     
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: FQMan(:),FlowSolution_Init(:)
    REAL(KIND=dp) :: UMan(:,:)    
    REAL(KIND=dp) :: USAV(:,:,:,:)
    REAL(KIND=dp) :: GradSAV( :,:,:,:,: )   
    INTEGER       :: NSDOFs, IO
    REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:)     
    TYPE(ValueList_t),POINTER :: Material    
    REAL(KIND=dp) :: Density(:)
    REAL(KIND=dp) :: FQManTemp(:)
    
    INTEGER, POINTER :: FlowPerm(:)
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    INTEGER :: t,n,nb,nd
    REAL(KIND=dp), POINTER :: FQlocal(:)   
    
    
    INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
                            
!------------------------------------------------------------------------------

    FQManTemp = 0.0_dp
    DO t = 1 , GetNOFActive()  ! Boucle ELEMENT
      CALL AdvanceOutput( t,GetNOFActive() )      
      Element => GetActiveElement(t)      
!       NodeIndexes => Element % NodeIndexes
      Indexes => Element % NodeIndexes
      
      n = GetElementNOFNodes(Element)      
      nb = GetElementNOFBDOFs(Element)       
      nd = GetElementDOFs( Indexes )
      CALL GetElementNodes( ElementNodes )
!
      Density(1:n) = GetReal( Material, 'Density' )

!------------------------------------------------------------------------------
!     Vitesse au point d'intégration de l'itération précédente
!------------------------------------------------------------------------------
      SELECT CASE( NSDOFs )
        CASE(3) ! 2D
          Uelex(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
          Ueley(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
          Uelez(1:n) = 0.0_dp
        CASE(4) ! 3D
          Uelex(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-3 , IO-1)
          Ueley(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
          Uelez(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
      END SELECT
!------------------------------------------------------------------------------  
      CALL MANBranchFQelemOPTI( FQManTemp, FQlocal, Density,                     &
                          Element, n, IO, ElementNodes,                    &
                          USAV, GradSAV, Element % ElementIndex , NSDOFs,  &
                          Uelex, Ueley, Uelez )
!------------------------------------------------------------------------------                                       

      CALL UpdateGlobalForce( FQMan, FQManTemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
    END DO ! FIN Boucle ELEMENT
!------------------------------------------------------------------------------
!       Conditions de Diriclet modifiées pour le 2nd membre----> 0           
    CALL CalCondLimCoeff( FQMan, 0.0_dp, 1,FlowSolution_Init )
!                           
 END SUBROUTINE ANMBranchrhsFQ
!  ___ ___    _   _  _  ___ _  _ 
! | _ ) _ \  /_\ | \| |/ __| || |
! | _ \   / / _ \| .` | (__| __ |
! |___/_|_\/_/ \_\_|\_|\___|_||_|
! Sum j=2, p-1 Q (Uj, U p+1-j )  : positif car passé en neg pour le second membre
 SUBROUTINE  MANBranchFQelemOPTI( FQManTemp, FQlocal, Nodalrho,     &
                            Element, n, io, ElementNodes,     &
                            USAV , GradSAV, numel , NSDOFs,   &
                            Uelex, Ueley, Uelez )
! Q + Q + !
  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE Adaptive
  USE SolverUtils

        IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: FQManTemp(:) 
     REAL(KIND=dp), POINTER :: FQlocal(:)
     REAL(KIND=dp) :: Nodalrho(:)
  
     ! Composante de la vitesse MAN ordre r et p-r
     INTEGER :: io, n   ! p = ordremanencours
     ! OMan ordre de la série MAN, nmax le nombre max de degré de liberté
     INTEGER, POINTER :: Flowperm(:)
      
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: Element
!----     
     REAL(KIND=dp) :: USAV(:,:,:,:)
     REAL(KIND=dp) :: GradSAV( :,:,:,:,: )
     INTEGER       :: numel,NSDOFs
     REAL(KIND=dp) ,ALLOCATABLE, TARGET :: UMan(:,:)     
     REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:)     
!----   

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp) :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)


     REAL(KIND=dp) :: Ur(3), GradUpmr(3,3), FQelemtemp(4)

     INTEGER :: i,j,k,l,p,c,q,t,r,dim

     REAL(KIND=dp) :: s,u,v,w,volume
     REAL(KIND=dp) :: rho
    
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, stat
!------------------------------------------------------------------------------
     
     REAL(KIND=dp) :: tmp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!              write(*,*) 'io=',io
     dim = CoordinateSystemDimension()
     Params => GetSolverParams()

     c = dim + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!     
     FQManTemp = 0.0_dp 
!!!!!!!!!!!!!!!!!!!!!!!!!!     
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis = n

     IntegStuff = GaussPoints( Element )

     U_Integ => IntegStuff % u   
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------    
!******************************************************************************
    DO t=1,N_Integ                        ! boucle sur les points d'intégration
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
              Basis, dBasisdx )
!------------------------------------------------------------------------------
      s = detJ * S_Integ(t)

      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
! 
!------------------------------------------------------------------------------
!!!!!!!!!! Termes déjà caclulé pour le second membre Lt Vchap = FQ, donc ici on utilise
!!!!!!!!!! directement les termes stockés dans USAV ET GRADUSAV
!!!!!!!!!!
! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! On sauvegarde U(element,gausspoint)  et Grad U(element,gausspoint) à l'ordre en cours
! ! ! ! ! ! en procédant ainsi on ne recalcule pas à chaques ordre, des termes déjà calculé aux ordres précédents.
! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !       USAV( numel , t , io-1 , 1 ) = SUM( Basis(1:n) * Uelex(1:n) )
! ! ! ! !       USAV( numel , t , io-1 , 2 ) = SUM( Basis(1:n) * Ueley(1:n) )
! ! ! ! !       IF (dim > 2) USAV( numel , t , io-1 , 3) = SUM( Basis(1:n) * Uelez(1:n) )    
! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! Puis le gradient
! ! ! ! ! ! Probleme : comment stocker : element, gausspoint, ordre, composante1, composante2
! ! ! ! !       DO j=1,3
! ! ! ! !          GradSAV( numel , t , io-1 , 1 , j ) = SUM( Uelex(1:n) * dBasisdx(1:n,j) )
! ! ! ! !          GradSAV( numel , t , io-1 , 2 , j ) = SUM( Ueley(1:n) * dBasisdx(1:n,j) )
! ! ! ! !          IF ( DIM > 2 ) GradSAV( numel , t , io-1 , 3 , j ) = SUM( Uelez(1:n) * dBasisdx(1:n,j) )
! ! ! ! !       END DO
!------------------------------------------------------------------------------
! Somme des termes croisées par composante spatiale
!------------------------------------------------------------------------------
      FQelemtemp = 0.0_dp
!       DO r=1, io-1                  ! Boucle CROSS
      DO r=2, io-1                  ! Boucle CROSS
      DO i = 1 ,dim
        DO j = 1 ,dim
        !! Uj, Up+1-j : U2,U2 => U3,U2 + U2,U3 => U4,U2 + U3,U3 + U2,U4
          FQelemtemp(i) = FQelemtemp(i) + USAV( numel , t , r , j ) * GradSAV( numel , t , io+1-r , i , j )
        END DO
      END DO
      END DO ! FIN Boucle CROSS
      
!------------------------------------------------------------------------------
! Assemblage
      DO p = 1 , NBasis !boucle sur les fonctions de formes des fonctions tests Wi = Ni
        FQlocal => FQManTemp( c*(p-1)+1 : c*(p-1)+c  )
        DO i = 1, c
          FQlocal(i) = FQlocal(i) + s * rho * FQelemtemp(i) * Basis(p)
        END DO
      END DO ! fin boucle sur les noeuds

!------------------------------------------------------------------------------
   END DO                       ! fin de la boucle sur les points d'intégration
!******************************************************************************   

 END SUBROUTINE MANBranchFQelemOPTI
!  ___ ___    _   _  _  ___ _  _ 
! | _ ) _ \  /_\ | \| |/ __| || |
! | _ \   / / _ \| .` | (__| __ |
! |___/_|_\/_/ \_\_|\_|\___|_||_|
         
!  __     __                  __      __ 
! |_ |\ ||  \   /\ |\ ||\/|  |__)|__|(_  
! |__| \||__/  /--\| \||  |  | \ |  |__)                                         
!------------------------------------------------------------------------------    
  
 
 
 
 
 
 
 
 
 
 
 
!------------------------------------------------------------------------------    
!  _____  _    _  _____ 
! |  __ \| |  | |/ ____|
! | |__) | |__| | (___  
! |  _  /|  __  |\___ \ 
! | | \ \| |  | |____) |
! |_|  \_\_|  |_|_____/ 
!                       
!  _____ _   _ _____ _____ _____       _______ ____  _____  
! |_   _| \ | |  __ \_   _/ ____|   /\|__   __/ __ \|  __ \ 
!   | | |  \| | |  | || || |       /  \  | | | |  | | |__) |
!   | | | . ` | |  | || || |      / /\ \ | | | |  | |  _  / 
!  _| |_| |\  | |__| || || |____ / ____ \| | | |__| | | \ \ 
! |_____|_| \_|_____/_____\_____/_/    \_\_|  \____/|_|  \_\
! FQ Steady Bif Indicator
! Cadou 2006
  SUBROUTINE  ANMSBIrhsFnl( INDFNL, INDDU0, INDDUserie, IO, NSDOFs, FlowPerm, &
                         USAV, GradSAV,                                       &
                         INDDUSAV, GradINDDUSAV,                              &
                         INDFNLtemp, INDUelex , INDUeley , INDUelez,          &
                         Density,Material,NDL,FlowSolution_Init)
!  

    USE DefUtils
  USE Differentials
  USE MaterialModels
  USE Adaptive    
    USE SolverUtils
!
    IMPLICIT NONE     
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: INDFNL(:) , INDDU0(:),FlowSolution_Init(:)
    REAL(KIND=dp) :: INDDUserie(:,:)    
    REAL(KIND=dp) :: USAV(:,:,:,:)
    REAL(KIND=dp) :: GradSAV( :,:,:,:,: )  
    REAL(KIND=dp) :: INDDUSAV(:,:,0:,:)
    REAL(KIND=dp) :: GradINDDUSAV( :,:,0:,:,: )      
    INTEGER       :: NSDOFs, IO , IK, NDL
    REAL(KIND=dp) :: INDUelex(:), INDUeley(:), INDUelez(:)
    TYPE(ValueList_t),POINTER :: Material    
    REAL(KIND=dp) :: Density(:)
    REAL(KIND=dp) :: INDFNLtemp(:)
    
    INTEGER, POINTER :: FlowPerm(:)
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    INTEGER :: t,n,nb,nd
    REAL(KIND=dp), POINTER :: FQlocal(:)   
    
    
    INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
                            
!------------------------------------------------------------------------------
!          DO IK=1,NDL
!           if (INDDU0(IK).GT.1e-15) THEN
!              WRITE(*,*) "INDDU0(",IK,")=",INDDU0(IK)
!           ENDIF
!          ENDDO
         
    INDFNLtemp = 0.0_dp
    DO t = 1 , GetNOFActive()  ! Boucle ELEMENT
      CALL AdvanceOutput( t,GetNOFActive() )      
      Element => GetActiveElement(t)      
!       NodeIndexes => Element % NodeIndexes
      Indexes => Element % NodeIndexes
      
      n = GetElementNOFNodes(Element)      
      nb = GetElementNOFBDOFs(Element)       
      nd = GetElementDOFs( Indexes )
      CALL GetElementNodes( ElementNodes )
!
      Density(1:n) = GetReal( Material, 'Density' )

!------------------------------------------------------------------------------
!     Fnl_1 = - [ Q( INDU_0 , U1) + Q( U1 , INDU_0 ) ]
!=>   Fnl_k = - SUM_{i=1}^k [ Q( INDU_k-i , Ui) + Q( Ui , INDU_k-i ) ]       
! REMARQUES:
! USAV ET GRADUSAV deja en mémoire, pas besoin de reprendre les infos aux elements.
! Par contre il faut avoir INDDU et Grad INDDU, donc même technique en OPTIM
!------------------------------------------------------------------------------
      IF(IO==1) THEN         ! ON PREND Delta U_O en direct car n'est plus passé en tab(:,0)....
        SELECT CASE( NSDOFs )
          CASE(3) ! 2D
            INDUelex(1:n) = INDDU0( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
            INDUeley(1:n) = INDDU0( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
            INDUelez(1:n) = 0.0_dp          
          CASE(4) ! 3D
            INDUelex(1:n) = INDDU0( NSDOFs*FlowPerm(Indexes(1:nd))-3 )
            INDUeley(1:n) = INDDU0( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
            INDUelez(1:n) = INDDU0( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
        END SELECT
      ELSE
        SELECT CASE( NSDOFs )
          CASE(3) ! 2D
            INDUelex(1:n) = INDDUserie( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
            INDUeley(1:n) = INDDUserie( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
            INDUelez(1:n) = 0.0_dp          
          CASE(4) ! 3D
            INDUelex(1:n) = INDDUserie( NSDOFs*FlowPerm(Indexes(1:nd))-3 , IO-1)
            INDUeley(1:n) = INDDUserie( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
            INDUelez(1:n) = INDDUserie( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
        END SELECT      
      ENDIF
!------------------------------------------------------------------------------  
      CALL ANMSBIrhsFnlElemOPTI( INDFNLtemp, FQlocal, Density, Element, n, io,           & 
                                 ElementNodes, USAV, GradSAV, INDDUSAV , GradINDDUSAV,   &
                                 Element % ElementIndex , NSDOFs,                        &
                                 INDUelex, INDUeley, INDUelez )
!------------------------------------------------------------------------------                                       

      CALL UpdateGlobalForce( INDFNL, INDFNLtemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
    END DO ! FIN Boucle ELEMENT
!------------------------------------------------------------------------------
!       Conditions de Diriclet modifiées pour le 2nd membre----> 0           
    CALL CalCondLimCoeff( INDFNL, 0.0_dp, 1,FlowSolution_Init )
!                           
 END SUBROUTINE ANMSBIrhsFnl

!!!! INDICATOR  RHS !!!!

 SUBROUTINE  ANMSBIrhsFnlElemOPTI( FTemp, FQlocal, Nodalrho,     &
                            Element, n, io, ElementNodes,     &
                            USAV, GradSAV,INDDUSAV , GradINDDUSAV, numel , NSDOFs,   &
                            INDUelex, INDUeley, INDUelez )
!
  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE Adaptive
  USE SolverUtils

        IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: FTemp(:) 
     REAL(KIND=dp), POINTER :: FQlocal(:)
     REAL(KIND=dp) :: Nodalrho(:)
  
     ! Composante de la vitesse MAN ordre r et p-r
     INTEGER :: io, n   ! p = ordremanencours
     ! OMan ordre de la série MAN, nmax le nombre max de degré de liberté
     INTEGER, POINTER :: Flowperm(:)
      
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: Element
!----     
     REAL(KIND=dp) :: USAV(:,:,:,:)
     REAL(KIND=dp) :: GradSAV( :,:,:,:,: )
     REAL(KIND=dp) :: INDDUSAV(:,:,0:,:)
     REAL(KIND=dp) :: GradINDDUSAV( :,:,0:,:,: )     
     INTEGER       :: numel,NSDOFs
     REAL(KIND=dp) ,ALLOCATABLE, TARGET :: UMan(:,:)     
     REAL(KIND=dp) :: INDUelex(:), INDUeley(:), INDUelez(:)     
!----   

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp) :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)


     REAL(KIND=dp) :: Ur(3), GradUpmr(3,3), FQelemtemp(4)

     INTEGER :: i,j,k,l,p,c,q,t,r,dim

     REAL(KIND=dp) :: s,u,v,w,volume
     REAL(KIND=dp) :: rho
    
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, stat
!------------------------------------------------------------------------------
     
     REAL(KIND=dp) :: tmp
!------------------------------------------------------------------------------
!!!! INDICATOR  RHS !!!!
!------------------------------------------------------------------------------
!              write(*,*) 'io=',io
     dim = CoordinateSystemDimension()
     Params => GetSolverParams()

     c = dim + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!     
     FTemp = 0.0_dp 
!!!!!!!!!!!!!!!!!!!!!!!!!!     
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis = n

     IntegStuff = GaussPoints( Element )

     U_Integ => IntegStuff % u   
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------    
!              write(*,*) 'io=',io
!******************************************************************************
    DO t=1,N_Integ                        ! boucle sur les points d'intégration
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
              Basis, dBasisdx )
!------------------------------------------------------------------------------
      s = detJ * S_Integ(t)

      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
! 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! On sauvegarde U(element,gausspoint)  et Grad U(element,gausspoint) à l'ordre en cours
! en procédant ainsi on ne recalcule pas à chaques ordre, des termes déjà calculé aux ordres précédents.
!------------------------------------------------------------------------------
!       WRITE(*,*) " numel , t ,", io-1," : INDUelex",INDUelex
!       WRITE(*,*) " numel , t ,", io-1," : INDUeley",INDUeley      
      INDDUSAV( numel , t , io-1 , 1 ) = SUM( Basis(1:n) * INDUelex(1:n) )
      INDDUSAV( numel , t , io-1 , 2 ) = SUM( Basis(1:n) * INDUeley(1:n) )
      IF (dim > 2) INDDUSAV( numel , t , io-1 , 3) = SUM( Basis(1:n) * INDUelez(1:n) )    
!------------------------------------------------------------------------------
! Puis le gradient
! Probleme : comment stocker : element, gausspoint, ordre, composante1, composante2
      DO j=1,3
         GradINDDUSAV( numel , t , io-1 , 1 , j ) = SUM( INDUelex(1:n) * dBasisdx(1:n,j) )
         GradINDDUSAV( numel , t , io-1 , 2 , j ) = SUM( INDUeley(1:n) * dBasisdx(1:n,j) )
         IF ( DIM > 2 ) GradINDDUSAV( numel , t , io-1 , 3 , j ) = SUM( INDUelez(1:n) * dBasisdx(1:n,j) )
      END DO
!------------------------------------------------------------------------------
! Somme des termes croisées par composante spatiale
!        Fnl_k = - SUM_{i=1}^k [ Q( INDU_k-i , Ui) + Q( Ui , INDU_k-i ) ]         
!------------------------------------------------------------------------------
      FQelemtemp = 0.0_dp
      DO r=1, io                  ! r [1,N] et io [N-1,0]
      DO i = 1 ,dim
        DO j = 1 ,dim
          FQelemtemp(i) = FQelemtemp(i) - INDDUSAV( numel , t , io - r , j ) *      GradSAV( numel , t ,      r , i , j )                  
          FQelemtemp(i) = FQelemtemp(i) -     USAV( numel , t ,      r , j ) * GradINDDUSAV( numel , t , io - r , i , j )
        END DO
      END DO
      END DO ! FIN Boucle CROSS
      
!------------------------------------------------------------------------------
! Assemblage
      DO p = 1 , NBasis !boucle sur les noeuds
        FQlocal => FTemp( c*(p-1)+1 : c*(p-1)+c  )
        DO i = 1, c
          FQlocal(i) = FQlocal(i) + s * rho * FQelemtemp(i)*Basis(p)
        END DO
      END DO ! fin boucle sur les noeuds

!------------------------------------------------------------------------------
   END DO                       ! fin de la boucle sur les points d'intégration
!******************************************************************************   

 END SUBROUTINE ANMSBIrhsFnlElemOPTI

!------------------------------------------------------------------------------    
!  ______ _   _ _____  
! |  ____| \ | |  __ \ 
! | |__  |  \| | |  | |
! |  __| | . ` | |  | |
! | |____| |\  | |__| |
! |______|_| \_|_____/ 
!  _____  _    _  _____ 
! |  __ \| |  | |/ ____|
! | |__) | |__| | (___  
! |  _  /|  __  |\___ \ 
! | | \ \| |  | |____) |
! |_|  \_\_|  |_|_____/ 
!                       
!  _____ _   _ _____ _____ _____       _______ ____  _____  
! |_   _| \ | |  __ \_   _/ ____|   /\|__   __/ __ \|  __ \ 
!   | | |  \| | |  | || || |       /  \  | | | |  | | |__) |
!   | | | . ` | |  | || || |      / /\ \ | | | |  | |  _  / 
!  _| |_| |\  | |__| || || |____ / ____ \| | | |__| | | \ \ 
! |_____|_| \_|_____/_____\_____/_/    \_\_|  \____/|_|  \_\
                     
!------------------------------------------------------------------------------    
   

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

   
   
   
   
   
   
!  _   _ _         ____                       _             
! | \ | | |       / __ \                     | |            
! |  \| | |      | |  | |_ __   ___ _ __ __ _| |_ ___  _ __ 
! | . ` | |      | |  | | '_ \ / _ \ '__/ _` | __/ _ \| '__|
! | |\  | |____  | |__| | |_) |  __/ | | (_| | || (_) | |   
! |_| \_|______|  \____/| .__/ \___|_|  \__,_|\__\___/|_|   
!                       | |                                 
!                       |_|                                 
! Q(A,B) = A * Grad(B)
! utile pour Residu
! sert pour ABE : PsiT * Q( LSK , LSK ) = 0
   
  SUBROUTINE  ANMQAB( Q , A, B, NSDOFs, FlowPerm, &
                         QTemp, &
                         Aelex , Aeley , Aelez, &
                         Belex , Beley , Belez, &                         
                         Density,Material)
!  

    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
!     USE ANMToolBox
    
!
    IMPLICIT NONE     
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Q(:) , A(:) , B(:)
    INTEGER       :: NSDOFs
    REAL(KIND=dp) :: Aelex(:), Aeley(:), Aelez(:)     
    REAL(KIND=dp) :: Belex(:), Beley(:), Belez(:)    
    TYPE(ValueList_t),POINTER :: Material    
    REAL(KIND=dp) :: Density(:)
    REAL(KIND=dp) :: QTemp(:)
    
    INTEGER, POINTER :: FlowPerm(:)
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    INTEGER :: t,n,nb,nd
    REAL(KIND=dp), POINTER :: Qlocal(:)   
    
    
    INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
                            
!------------------------------------------------------------------------------

    QTemp = 0.0_dp
    Q = 0.0_dp
    DO t = 1 , GetNOFActive()  ! Boucle ELEMENT
      CALL AdvanceOutput( t,GetNOFActive() )      
      Element => GetActiveElement(t)      
!       NodeIndexes => Element % NodeIndexes
      Indexes => Element % NodeIndexes
      
      n = GetElementNOFNodes(Element)      
      nb = GetElementNOFBDOFs(Element)       
      nd = GetElementDOFs( Indexes )
      CALL GetElementNodes( ElementNodes )
!
      Density(1:n) = GetReal( Material, 'Density' )

!------------------------------------------------------------------------------
!     Vitesse au point d'intégration de l'itération précédente
!------------------------------------------------------------------------------
      SELECT CASE( NSDOFs )
        CASE(3) ! 2D
          Aelex(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Aeley(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
          Aelez(1:n) = 0.0_dp
          Belex(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Beley(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
          Belez(1:n) = 0.0_dp          
        CASE(4) ! 3D
          Aelex(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-3 )
          Aeley(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Aelez(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
          Belex(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-3 )
          Beley(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Belez(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-1 )          
      END SELECT
!------------------------------------------------------------------------------  
      CALL MANQABelem( QTemp, Qlocal, Density,           &
                       Element, n, ElementNodes,     &
                       Element % ElementIndex , NSDOFs,  &
                       Aelex, Aeley, Aelez ,             &
                       Belex, Beley, Belez )
!------------------------------------------------------------------------------                                       

      CALL UpdateGlobalForce( Q , QTemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
    END DO ! FIN Boucle ELEMENT
!------------------------------------------------------------------------------
!       Conditions de Diriclet modifiées pour le 2nd membre----> 0           
! !     CALL CalCondLimCoeff( Q , 0.0_dp, 1,FlowSolution_Init )
!                           
 END SUBROUTINE ANMQAB



 SUBROUTINE  MANQABelem( FQManTemp, FQlocal, Nodalrho,     &
                            Element, n, ElementNodes,     &
                            numel , NSDOFs,   &
                            Aelex, Aeley, Aelez ,             &
                            Belex, Beley, Belez )
!
  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE Adaptive
  USE SolverUtils

        IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: FQManTemp(:) 
     REAL(KIND=dp), POINTER :: FQlocal(:)
     REAL(KIND=dp) :: Nodalrho(:)
  
     INTEGER ::  n   
     ! OMan ordre de la série MAN, nmax le nombre max de degré de liberté
     INTEGER, POINTER :: Flowperm(:)
      
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: Element
!----     
     INTEGER       :: numel,NSDOFs
     REAL(KIND=dp) ,ALLOCATABLE, TARGET :: UMan(:,:)     
     REAL(KIND=dp) :: Aelex(:), Aeley(:), Aelez(:)     
     REAL(KIND=dp) :: Belex(:), Beley(:), Belez(:)     
!----   
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp) :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)

     REAL(KIND=dp) :: AL(3), gradBL(3,3), FQelemtemp(4)

     INTEGER :: i,j,k,l,p,c,q,t,r,dim

     REAL(KIND=dp) :: s,u,v,w,volume
     REAL(KIND=dp) :: rho
    
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, stat
!------------------------------------------------------------------------------     
     REAL(KIND=dp) :: tmp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     Params => GetSolverParams()

     c = dim + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!     
     FQManTemp = 0.0_dp 
!!!!!!!!!!!!!!!!!!!!!!!!!!     
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis = n

     IntegStuff = GaussPoints( Element )

     U_Integ => IntegStuff % u   
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------    
!******************************************************************************
    DO t=1,N_Integ                        ! boucle sur les points d'intégration
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
              Basis, dBasisdx )
!------------------------------------------------------------------------------
      s = detJ * S_Integ(t)

      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! On sauvegarde U(element,gausspoint)  et Grad U(element,gausspoint) à l'ordre en cours
! en procédant ainsi on ne recalcule pas à chaques ordre, des termes déjà calculé aux ordres précédents.
!------------------------------------------------------------------------------
      AL( 1 ) = SUM( Basis(1:n) * Aelex(1:n) )
      AL( 2 ) = SUM( Basis(1:n) * Aeley(1:n) )
      IF (dim > 2) AL( 3 ) = SUM( Basis(1:n) * Aelez(1:n) )    
!------------------------------------------------------------------------------
! Puis le gradient
! Probleme : comment stocker : element, gausspoint, ordre, composante1, composante2
      DO j=1,3
         gradBL( 1 , j ) = SUM( Belex(1:n) * dBasisdx(1:n,j) )
         gradBL( 2 , j ) = SUM( Beley(1:n) * dBasisdx(1:n,j) )
         IF ( DIM > 2 ) gradBL( 3 , j ) = SUM( Belez(1:n) * dBasisdx(1:n,j) )
      END DO
!------------------------------------------------------------------------------
! Somme des termes croisées par composante spatiale
!------------------------------------------------------------------------------
      FQelemtemp = 0.0_dp
      DO i = 1 ,dim
        DO j = 1 ,dim
          FQelemtemp(i) = FQelemtemp(i) + AL( j ) * gradBL( i , j )
        END DO
      END DO
      
!------------------------------------------------------------------------------
! Assemblage
      DO p = 1 , NBasis !boucle sur les noeuds
        FQlocal => FQManTemp( c*(p-1)+1 : c*(p-1)+c  )
        DO i = 1, c
          FQlocal(i) = FQlocal(i) + s * rho * FQelemtemp(i) * Basis(p)
        END DO
      END DO ! fin boucle sur les noeuds

!------------------------------------------------------------------------------
   END DO                       ! fin de la boucle sur les points d'intégration
!******************************************************************************   

 END SUBROUTINE MANQABelem


 
 
 
 
 
 
 
 
! Q(A,B) Version Stabilisée Franca et. al 1992 : ONLY Q(A,B) RHS IS CONCERNED HERE 
 
  SUBROUTINE  ANMQABSTAB( Q , A, B, NSDOFs, FlowPerm, &
                         QTemp, &
                         Aelex , Aeley , Aelez, &
                         Belex , Beley , Belez, &                         
                         Density,Material,Nodalmu)
!  

    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
!     USE ANMToolBox
    
!
    IMPLICIT NONE     
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Q(:) , A(:) , B(:)
    INTEGER       :: NSDOFs
    REAL(KIND=dp) :: Aelex(:), Aeley(:), Aelez(:)     
    REAL(KIND=dp) :: Belex(:), Beley(:), Belez(:)    
    TYPE(ValueList_t),POINTER :: Material    
    REAL(KIND=dp) :: Density(:),Nodalmu(:)
    REAL(KIND=dp) :: QTemp(:)
    
    INTEGER, POINTER :: FlowPerm(:)
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    INTEGER :: t,n,nb,nd
    REAL(KIND=dp), POINTER :: Qlocal(:)   
    
    
    INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
                            
!------------------------------------------------------------------------------

    QTemp = 0.0_dp
    Q = 0.0_dp
    DO t = 1 , GetNOFActive()  ! Boucle ELEMENT
      CALL AdvanceOutput( t,GetNOFActive() )      
      Element => GetActiveElement(t)      
!       NodeIndexes => Element % NodeIndexes
      Indexes => Element % NodeIndexes
      
      n = GetElementNOFNodes(Element)      
      nb = GetElementNOFBDOFs(Element)       
      nd = GetElementDOFs( Indexes )
      CALL GetElementNodes( ElementNodes )
!
      Density(1:n) = GetReal( Material, 'Density'   )
      Nodalmu(1:n) = GetReal( Material, 'Viscosity' )

!------------------------------------------------------------------------------
!     Vitesse au point d'intégration de l'itération précédente
!------------------------------------------------------------------------------
      SELECT CASE( NSDOFs )
        CASE(3) ! 2D
          Aelex(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Aeley(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
          Aelez(1:n) = 0.0_dp
          Belex(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Beley(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
          Belez(1:n) = 0.0_dp          
        CASE(4) ! 3D
          Aelex(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-3 )
          Aeley(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Aelez(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
          Belex(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-3 )
          Beley(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
          Belez(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-1 )          
      END SELECT
!------------------------------------------------------------------------------  
      CALL MANQABelemSTAB( QTemp, Qlocal, Density,           &
                       Element, n, ElementNodes,     &
                       Element % ElementIndex , NSDOFs,  &
                       Aelex, Aeley, Aelez ,             &
                       Belex, Beley, Belez,Nodalmu )
!------------------------------------------------------------------------------                                       

      CALL UpdateGlobalForce( Q , QTemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
    END DO ! FIN Boucle ELEMENT
!------------------------------------------------------------------------------
!       Conditions de Diriclet modifiées pour le 2nd membre----> 0           
! !     CALL CalCondLimCoeff( Q , 0.0_dp, 1,FlowSolution_Init )
!                           
 END SUBROUTINE ANMQABSTAB



 SUBROUTINE  MANQABelemSTAB( FQManTemp, FQlocal, Nodalrho,     &
                            Element, n, ElementNodes,     &
                            numel , NSDOFs,   &
                            Aelex, Aeley, Aelez ,             &
                            Belex, Beley, Belez ,Nodalmu)
!
  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE Adaptive
  USE SolverUtils

        IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: FQManTemp(:) 
     REAL(KIND=dp), POINTER :: FQlocal(:)
     REAL(KIND=dp) :: Nodalrho(:)
  
     INTEGER ::  n   
     ! OMan ordre de la série MAN, nmax le nombre max de degré de liberté
     INTEGER, POINTER :: Flowperm(:)
      
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: Element
!----     
     INTEGER       :: numel,NSDOFs
     REAL(KIND=dp) ,ALLOCATABLE, TARGET :: UMan(:,:)     
     REAL(KIND=dp) :: Aelex(:), Aeley(:), Aelez(:)     
     REAL(KIND=dp) :: Belex(:), Beley(:), Belez(:)     
!----   
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp) :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)

     REAL(KIND=dp) :: AL(3), gradBL(3,3), FQelemtemp(4)

     INTEGER :: i,j,k,l,p,c,q,t,r,dim

     REAL(KIND=dp) :: s,u,v,w,volume
     REAL(KIND=dp) :: rho
    
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, stat
!------------------------------------------------------------------------------     
     REAL(KIND=dp) :: tmp
!------------------------------------------------------------------------------     
! STABILIZE TEST     
     REAL(KIND=dp) :: VNorm,hK,mK,mu,dmudx(3),Re,Tau,Delta,Lambda,hScale
     REAL(KIND=dp) :: SW(n,4,4)
     REAL(KIND=dp) :: Nodalmu(:)
     LOGICAL ::  Stabilize
     
!------------------------------------------------------------------------------
     Stabilize = .TRUE.
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     Params => GetSolverParams()

     c = dim + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!     
     FQManTemp = 0.0_dp 
!!!!!!!!!!!!!!!!!!!!!!!!!!     
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis = n

     IntegStuff = GaussPoints( Element )

     U_Integ => IntegStuff % u   
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!------------------------------------------------------------------------------
     hScale = GetCReal( Params, 'H scale', Found )
     IF ( .NOT. Found )  hScale = 1._dp
     hK = Element % hK*hscale
     mK = Element % StabilizationMK

!      IF ( Stabilize ) THEN
       dNodalBasisdx = 0._dp
       DO p=1,n
         u = Element % TYPE % NodeU(p)
         v = Element % TYPE % NodeV(p)
         w = Element % TYPE % NodeW(p)
         stat = ElementInfo( Element, ElementNodes, u,v,w, detJ, Basis, dBasisdx )
         dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
       END DO
!      ENDIF
     
     
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------    
!******************************************************************************
    DO t=1,N_Integ                        ! boucle sur les points d'intégration
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
              Basis, dBasisdx )
!------------------------------------------------------------------------------
      s = detJ * S_Integ(t)

      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! On sauvegarde U(element,gausspoint)  et Grad U(element,gausspoint) à l'ordre en cours
! en procédant ainsi on ne recalcule pas à chaques ordre, des termes déjà calculé aux ordres précédents.
!------------------------------------------------------------------------------
      AL( 1 ) = SUM( Basis(1:n) * Aelex(1:n) )
      AL( 2 ) = SUM( Basis(1:n) * Aeley(1:n) )
      IF (dim > 2) AL( 3 ) = SUM( Basis(1:n) * Aelez(1:n) )    
!------------------------------------------------------------------------------
! Puis le gradient
! Probleme : comment stocker : element, gausspoint, ordre, composante1, composante2
      DO j=1,3
         gradBL( 1 , j ) = SUM( Belex(1:n) * dBasisdx(1:n,j) )
         gradBL( 2 , j ) = SUM( Beley(1:n) * dBasisdx(1:n,j) )
         IF ( DIM > 2 ) gradBL( 3 , j ) = SUM( Belez(1:n) * dBasisdx(1:n,j) )
      END DO
!------------------------------------------------------------------------------
! Somme des termes croisées par composante spatiale
!------------------------------------------------------------------------------
      FQelemtemp = 0.0_dp
      DO i = 1 ,dim
        DO j = 1 ,dim
          FQelemtemp(i) = FQelemtemp(i) + AL( j ) * gradBL( i , j )
        END DO
      END DO
      

      
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
!       IF ( Stabilize ) THEN
        DO i=1,3
          dmudx(i) = SUM( Nodalmu(1:n)*dBasisdx(1:n,i) )
        END DO
        mu = SUM( Nodalmu(1:n) * Basis(1:n) )
        
!------------------------------------------------------------------------------
!       Stabilization parameters Tau & Delta
!------------------------------------------------------------------------------
!         IF ( Convect ) THEN
           VNorm = MAX( SQRT( SUM(AL(1:DIM)**2) ), 1.0d-12 )
           Re = MIN( 1.0d0, rho * mK * hK * VNorm / (4 * mu) )
           Lambda = 1.
           Tau = hK * Re / (2 * rho * VNorm)
           Delta = rho * Lambda * Re * hK * VNorm
!         ELSE
!            Delta = 0._dp
!            Tau   = mK * hK**2  / ( 8 * mu )
!         END IF
!------------------------------------------------------------------------------
!       SU will contain residual of ns-equations (except for the time derivative
!       and force terms). SW will contain the weight function values.
!------------------------------------------------------------------------------
!         SU(1:n,:,:) = 0.0D0
        SW(1:n,:,:) = 0.0D0

        DO p=1,N
          DO i=1,DIM
! ! ! !             SU(p,i,c) = SU(p,i,c) + dBasisdx(p,i)
! ! ! !             IF(Porous) THEN
! ! ! !               SU(p,i,i) = SU(p,i,i) + mu * Drag(i) * Basis(p)
! ! ! !             END IF
! ! ! ! 
! ! ! !             IF ( Convect ) THEN
! ! ! !               DO j=1,DIM
! ! ! !                 SU(p,i,i) = SU(p,i,i) + rho * dBasisdx(p,j) * Velo(j)
! ! ! !               END DO
! ! ! !             END IF
! ! ! ! 
! ! ! !             DO j=1,DIM
! ! ! !               SU(p,i,i) = SU(p,i,i) - dmudx(j) * dBasisdx(p,j)
! ! ! !               SU(p,i,j) = SU(p,i,j) - dmudx(j) * dBasisdx(p,i)
! ! ! !               SU(p,i,i) = SU(p,i,i) - mu * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
! ! ! !               SU(p,i,j) = SU(p,i,j) - mu * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
! ! ! !             END DO
! ! ! ! 
! ! ! !             IF ( Convect .AND. NewtonLinearization ) THEN
! ! ! !               DO j=1,DIM
! ! ! !                 SU(p,i,j) = SU(p,i,j) + rho * Grad(i,j) * Basis(p)
! ! ! !               END DO
! ! ! !             END IF
!
!------------------------------------------------------------------------------

!             IF ( Convect ) THEN
              SW(p,c,i) = SW(p,c,i) + rho * dBasisdx(p,i)
              DO j=1,dim
                SW(p,i,i) = SW(p,i,i) + rho * dBasisdx(p,j) * AL(j)
!                 SW(p,i,i) = SW(p,i,i) + rho * Grad(i,j)* AL(j) 
              END DO
!             ELSE
!               SW(p,c,i) = SW(p,c,i) + dBasisdx(p,i)
!             END IF

            DO j=1,dim
              SW(p,i,i) = SW(p,i,i) - dmudx(j) * dBasisdx(p,j)
              SW(p,j,i) = SW(p,j,i) - dmudx(j) * dBasisdx(p,i)
              SW(p,i,i) = SW(p,i,i) - mu * SUM(dNodalBasisdx(p,1:n,j)*dBasisdx(1:n,j))
              SW(p,j,i) = SW(p,j,i) - mu * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
            END DO
          END DO
        END DO
!       ENDIF
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test     
! ! ! ! ! ! - - - - - - - - V2 STABILIZED Franca Test        
!------------------------------------------------------------------------------
! Assemblage
      DO p = 1 , NBasis !boucle sur les noeuds
        FQlocal => FQManTemp( c*(p-1)+1 : c*(p-1)+c  )
        DO i = 1, c
          FQlocal(i) = FQlocal(i) + s * rho * FQelemtemp(i) * Basis(p)
        END DO
        
!        IF ( Stabilize ) THEN
         DO i=1,DIM
           DO j=1,c
             FQlocal(j) = FQlocal(j) + s * Tau * rho * FQelemtemp(i) * SW(p,j,i)
           END DO
         END DO
!        ENDIF        
        
      END DO ! fin boucle sur les noeuds

!------------------------------------------------------------------------------
   END DO                       ! fin de la boucle sur les points d'intégration
!******************************************************************************   

 END SUBROUTINE MANQABelemSTAB 
 
 
 
 
 
 
 
 
 
 
 
 
 

! ! ! !  
! ! ! ! ! Q(A,B) Version Stabilisée Franca et. al 1992 : ONLY Q(A,B) RHS IS CONCERNED HERE 
! ! ! !  
! ! ! !   SUBROUTINE  ANMQABSTABautonom( Q , A, B, QTemp, Model, Solver)
! ! ! ! !  
! ! ! !     USE DefUtils
! ! ! !     USE Differentials
! ! ! !     USE MaterialModels
! ! ! !     USE Adaptive    
! ! ! !     USE SolverUtils
! ! ! ! !     USE ANMToolBox
! ! ! !     
! ! ! ! !
! ! ! !     IMPLICIT NONE     
! ! ! ! !------------------------------------------------------------------------------
! ! ! !     TYPE(Model_t) :: Model
! ! ! !     TYPE(Solver_t), TARGET :: Solver
! ! ! ! 
! ! ! !     REAL(KIND=dp) :: Q(:) , A(:) , B(:)
! ! ! !   
! ! ! !     
! ! ! !     REAL(KIND=dp) :: QTemp(:)
! ! ! !     
! ! ! !      
! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !    Local variables
! ! ! ! !------------------------------------------------------------------------------
! ! ! !     INTEGER       :: NSDOFs
! ! ! !     REAL(KIND=dp) :: Aelex(:), Aeley(:), Aelez(:)     
! ! ! !     REAL(KIND=dp) :: Belex(:), Beley(:), Belez(:)    
! ! ! !     TYPE(ValueList_t),POINTER :: Material    
! ! ! !     REAL(KIND=dp) :: Density(:),Nodalmu(:)  
! ! ! !     INTEGER :: t,n,nb,nd, k
! ! ! !     REAL(KIND=dp), POINTER :: Qlocal(:)   
! ! ! !     
! ! ! !     
! ! ! !     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
! ! ! !     TYPE(Nodes_t) :: ElementNodes
! ! ! !     TYPE(Element_t), POINTER :: Element
! ! ! !     
! ! ! !     INTEGER, POINTER :: FlowPerm(:)
! ! ! !     
! ! ! !                             
! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !        NDL       = StiffMatrix % NumberOfRows
! ! ! !     N         = Solver % Mesh % MaxElementDOFs     
! ! ! !     ALLOCATE( Aelex(N), Aeley(N), Aelez(N),              &
! ! ! !               Belex(N), Beley(N), Belez(N),              & 
! ! ! !               Density(N),Nodalmu(N)
! ! ! !               STAT=istat )
! ! ! !               
! ! ! !     NSDOFs    = Solver % Variable % DOFs
! ! ! !     FlowPerm  => Solver % Variable % Perm
! ! ! ! 
! ! ! !     QTemp = 0.0_dp
! ! ! !     Q = 0.0_dp
! ! ! !     DO t = 1 , GetNOFActive()  ! Boucle ELEMENT
! ! ! !       CALL AdvanceOutput( t,GetNOFActive() )      
! ! ! !       Element => GetActiveElement(t)      
! ! ! ! !       NodeIndexes => Element % NodeIndexes
! ! ! !       Indexes => Element % NodeIndexes
! ! ! !       
! ! ! !       n = GetElementNOFNodes(Element)      
! ! ! !       nb = GetElementNOFBDOFs(Element)       
! ! ! !       nd = GetElementDOFs( Indexes )
! ! ! !       CALL GetElementNodes( ElementNodes )
! ! ! ! 
! ! ! ! !!    Material parameters: density, viscosity, etc.
! ! ! ! !    ----------------------------------------------
! ! ! !       k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
! ! ! !                   minv=1, maxv=Model % NumberOfMaterials )
! ! ! ! 
! ! ! !       Material => Model % Materials(k) % Values
! ! ! !       Density(1:n) = GetReal( Material, 'Density'   )
! ! ! !       Nodalmu(1:n) = GetReal( Material, 'Viscosity' )
! ! ! ! !    ----------------------------------------------
! ! ! ! 
! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !     Vitesse au point d'intégration de l'itération précédente
! ! ! ! !------------------------------------------------------------------------------
! ! ! !       SELECT CASE( NSDOFs )
! ! ! !         CASE(3) ! 2D
! ! ! !           Aelex(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
! ! ! !           Aeley(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
! ! ! !           Aelez(1:n) = 0.0_dp
! ! ! !           Belex(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
! ! ! !           Beley(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
! ! ! !           Belez(1:n) = 0.0_dp          
! ! ! !         CASE(4) ! 3D
! ! ! !           Aelex(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-3 )
! ! ! !           Aeley(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
! ! ! !           Aelez(1:n) = A( NSDOFs*FlowPerm(Indexes(1:nd))-1 )
! ! ! !           Belex(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-3 )
! ! ! !           Beley(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-2 )
! ! ! !           Belez(1:n) = B( NSDOFs*FlowPerm(Indexes(1:nd))-1 )          
! ! ! !       END SELECT
! ! ! ! !------------------------------------------------------------------------------  
! ! ! !       CALL MANQABelemSTAB( QTemp, Qlocal, Density,           &
! ! ! !                        Element, n, ElementNodes,             &
! ! ! !                        Element % ElementIndex , NSDOFs,      &
! ! ! !                        Aelex, Aeley, Aelez ,                 &
! ! ! !                        Belex, Beley, Belez,Nodalmu )
! ! ! ! !------------------------------------------------------------------------------                                       
! ! ! ! 
! ! ! !       CALL UpdateGlobalForce( Q , QTemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
! ! ! !     END DO ! FIN Boucle ELEMENT
! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !       Conditions de Diriclet modifiées pour le 2nd membre----> 0           
! ! ! ! ! !     CALL CalCondLimCoeff( Q , 0.0_dp, 1,FlowSolution_Init )
! ! ! ! 
! ! ! !     DEALLOCATE( Aelex, Aeley, Aelez, Belex, Beley, Belez,    & 
! ! ! !                 Density,Nodalmu)
! ! ! ! !                           
! ! ! !  END SUBROUTINE ANMQABSTABautonom
! ! ! !  








!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! .______       _______     _______. __   _______   __    __       ___       __      
! |   _  \     |   ____|   /       ||  | |       \ |  |  |  |     /   \     |  |     
! |  |_)  |    |  |__     |   (----`|  | |  .--.  ||  |  |  |    /  ^  \    |  |     
! |      /     |   __|     \   \    |  | |  |  |  ||  |  |  |   /  /_\  \   |  |     
! |  |\  \----.|  |____.----)   |   |  | |  '--'  ||  `--'  |  /  _____  \  |  `----.
! | _| `._____||_______|_______/    |__| |_______/  \______/  /__/     \__\ |_______|

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! RESIDUAL Version : Operators :
!      R(U0,L0) =  L(U0) + Q(U0,U0) - L0 * F  ! F = 0
!      R = [ L(U0) + Q(U0,U0) + Q(U0,U0) ]_stab - Qstab(U0,U0) 
!      R = LtStab(U0) - Qstab(U0,U0) 
!
! OUT : 
!       - ResNorm = ||Lstab(U0) + Qstab(U0,U0)||
!       - ResVec  =   Lstab(U0) + Qstab(U0,U0)
!------------------------------------------------------------------------------
! ! !    SUBROUTINE OperatorRESIDUALstab( ResNorm, ResVec, StiffMatrix, FlowSolution, &
! ! !                                     abeQAB, FlowSolution_Init,Model,Solver )
! ! !      USE DefUtils
! ! !      USE SolverUtils       
! ! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ! ! ! - - - IN
! ! !      TYPE(Matrix_t),POINTER   :: StiffMatrix
! ! !      REAL(KIND=dp)    :: FlowSolution(:),FlowSolution_Init(:),abeQAB(:)
! ! !      TYPE(Model_t) :: Model
! ! !      TYPE(Solver_t), TARGET :: Solver     
! ! ! ! - - - OUT
! ! !      REAL(KIND=dp)    :: ResNorm
! ! !      REAL(KIND=dp)    :: ResVec(:)
! ! ! ! - - - LOCAL
! ! ! 
! ! ! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ! !      
! ! !        CALL MatrixVectorMultiply(StiffMatrix,FlowSolution,ResVec)
! ! !        WRITE(*,*)  ' - RESIDSTAB interm ||Ltstab X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))  
! ! ! 
! ! !        CALL  ANMQABSTABautonom( abeQAB , FlowSolution, FlowSolution,  FQManTemp, &      
! ! !                                 Model, Solver)
! ! !        WRITE(*,*)  ' - RESIDSTAB interm ||Qstab(U0,U0)||  =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))                        
! ! !               
! ! ! !      R = [ L(U0) + Q(U0,U0) + Q(U0,U0) ] - Q(U0,U0) - L0 * F     
! ! !        ResVec = ResVec - abeQAB
! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
! ! !        ResNorm = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
! ! !        WRITE(*,*)  ' - RESIDSTAB R0= ||Lstab(U0) + Qstab(U0,U0)|| =',ResNorm
! ! ! 
! ! !    END SUBROUTINE OperatorRESIDUALstab
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   

! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! !      
! ! ! ! ! AJOUT 20161017 - Corrrection de Lt c pour voir si ca permet de meilleurs tangentes
! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !      R(U0,L0) =  L(U0) + Q(U0,U0) - L0 * F  ! F = 0
! ! ! ! !------------------------------------------------------------------------------
! ! ! !        ResVec = 0.0_dp
! ! ! !        CALL MatrixVectorMultiply(StiffMatrix,UBif,ResVec)
! ! ! !        WRITE(*,*)  '|RRB         ||Lt X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))  
! ! ! !        
! ! ! ! 
! ! ! !        IF (Stabilize) THEN
! ! ! !               CALL  ANMQABSTAB( abeQAB , UBif, UBif, NSDOFs, FlowPerm, &
! ! ! !                          FQManTemp, &
! ! ! !                          Uelex , Ueley , Uelez, &
! ! ! !                          Velex , Veley , Velez, &                           
! ! ! !                          DensityTMP,Material,NodalMuTMP)
! ! ! !               WRITE(*,*)  '|RRB  test||Qstab(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))   
! ! ! !               ResVec = ResVec - abeQAB
! ! ! !        ELSE
! ! ! !               CALL  ANMQAB( abeQAA , UBif, UBif, NSDOFs, FlowPerm, &
! ! ! !                          FQManTemp, &
! ! ! !                          Uelex , Ueley , Uelez, &
! ! ! !                          Velex , Veley , Velez, &                           
! ! ! !                          DensityTMP,Material)
! ! ! !               WRITE(*,*)  '|RRB         ||Q(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAA,abeQAA ))         
! ! ! !               ResVec = ResVec - abeQAA
! ! ! !        ENDIF
! ! ! !        
! ! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
! ! ! !        ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
! ! ! !        WRITE(*,*)  '|RESIDUAL as ||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab
! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>











!------------------------------------------------------------------------------
! COMPUTE RESIDUAL OVER THE BODY ELELENTS
! No boundary, No Edge : TODO
!------------------------------------------------------------------------------
  SUBROUTINE HMFlowResidual( U0, lambda0, ResidualVector, Rtemp, Model, Mesh , FlowPerm,&
                             NSDOFs, FNORMG , ResidualNorm)
!>  NEXT SOLUTION IS FORMED : U0 L0 with the actual computed series and amax  
!    R(U0,Lambda0) =  L(U0) + Q(U0,U0) - Lambda0 * F
! - - - - - 
       USE DefUtils
! - - - - - 
       IMPLICIT NONE
!------------------------------------------------------------------------------
       REAL(KIND=dp)           :: ResidualVector(:) , U0(:) , lambda0
       REAL(KIND=dp)           :: Rtemp(:)
! - - - - -      
       TYPE(Model_t)           :: Model
       TYPE( Mesh_t ), POINTER :: Mesh
       INTEGER, POINTER        :: FlowPerm(:)
       INTEGER                 :: NSDOFs
! - - - - - 
       REAL(KIND=dp)           :: Fnorm, FNORMG , ResidualNorm, ResidualNormEle
!------------------------------------------------------------------------------
!      Local variables
!------------------------------------------------------------------------------
       INTEGER                  :: t,n,nb,nd     
       INTEGER, POINTER         :: NodeIndexes(:), Indexes(:)
       TYPE(Element_t), POINTER :: Element
       TYPE(Nodes_t)            :: ElementNodes
!------------------------------------------------------------------------------
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -     
       FNORMG = 0.0_dp
       ResidualNorm =  0.0_dp
!      Boucle sur les éléments du domaine 
       DO t = 1,GetNOFActive()
!        -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -       
         CALL AdvanceOutput( t,GetNOFActive() )
         Element => GetActiveElement(t)     
         Indexes => Element % NodeIndexes
         n = GetElementNOFNodes(Element)      
!        -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!        LOCAL RESIDUAL
         CALL HMLocalFlowInsideResidual( Model, Element, Mesh, U0, FlowPerm, Fnorm, Rtemp, & 
                                         lambda0, ResidualNormEle )
!        -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
!        ASSEMBLY         
         CALL UpdateGlobalForce( ResidualVector , Rtemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
         FNORMG = FNORMG + Fnorm
         ResidualNorm = ResidualNorm + ResidualNormEle
!        -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -         
       END DO ! Fin boucle sur les éléments du domaine
       ResidualNorm = DSQRT(ResidualNorm)
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   
! ! !        CALL CalCondLimCoeff( ResidualVector , 0.0_dp, 1, FlowSolution_Init )
!   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   
!------------------------------------------------------------------------------
    END SUBROUTINE HMFlowResidual
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
!    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    
!  ___  _    ___  _   _  ___  _  _  ___  _   ___ __ __
! | __|| |  | __|| \_/ || __|| \| ||_ _|/ \ | o \\ V /
! | _| | |_ | _| | \_/ || _| | \\ | | || o ||   / \ / 
! |___||___||___||_| |_||___||_|\_| |_||_n_||_|\\ |_| 
!    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    
! |  | /~~\ |\  /||~~  |\  /|  /\  |~~\ |~~  |~~\|~~/~~\~|~|~~\ |   | /\  |  
! |--||    || \/ ||--  | \/ | /__\ |   ||--  |__/|--'--. | |   ||   |/__\ |  
! |  | \__/ |    ||__  |    |/    \|__/ |__  |  \|__\__/_|_|__/  \_//    \|__
!    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    
!------------------------------------------------------------------------------
!> Compute the residual of the Navier-Stokes equation for the bulk elements.
!> Home made version : 
!>  -> compute Norm as proposed by Elmer
!>  --> Out the Residual Vector!!
!
! Quant == FlowSolution
!    R(U0,Lambda0) =  L(U0) + Q(U0,U0) - Lambda0 * F
!!!!!!!! Lambda0 * F !!!!!!!!!!
!------------------------------------------------------------------------------
   SUBROUTINE HMLocalFlowInsideResidual( Model, Element,  &
          Mesh, Quant, Perm, Fnorm, Rtemp , lambda0, ResidualNorm) 
!------------------------------------------------------------------------------
     USE DefUtils
!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     INTEGER :: Perm(:)
     REAL(KIND=dp) :: Quant(:), Indicator(2), FNorm, lambda0
     TYPE( Mesh_t ), POINTER    :: Mesh
     TYPE( Element_t ), POINTER :: Element
!------------------------------------------------------------------------------

     TYPE(Nodes_t) :: Nodes

     INTEGER :: i,j,k,l,m,n,t,DIM,DOFs, p,NBasis,c

     LOGICAL :: stat, GotIt, Compressible, Convect

     TYPE( Variable_t ), POINTER :: Var

     REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
     REAL(KIND=dp) :: Density, Viscosity,u, v, w, s, detJ
     REAL(KIND=dp) :: Residual(4), ResidualNorm, Area, ReferencePressure, dt

     REAL(KIND=dp), ALLOCATABLE :: NodalViscosity(:), NodalDensity(:), &
            Basis(:),  dBasisdx(:,:), ddBasisddx(:,:,:)
     REAL(KIND=dp),ALLOCATABLE :: Velocity(:,:), Pressure(:)
     REAL(KIND=dp),ALLOCATABLE :: PrevVelo(:,:), PrevPres(:)
     REAL(KIND=dp),ALLOCATABLE :: Temperature(:), NodalForce(:,:)
     REAL(KIND=dp),ALLOCATABLE :: HeatCapacity(:), ReferenceTemperature(:), &
                      HeatExpansionCoeff(:)

     REAL(KIND=dp) :: SpecificHeatRatio

     REAL(KIND=dp), POINTER :: Gravity(:,:)

     TYPE(ValueList_t), POINTER :: Material

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
!------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: Rtemp(:) 
     REAL(KIND=dp), POINTER :: Rlocal(:)
     REAL(KINd=dp) :: Temp
!------------------------------------------------------------------------------

!    Initialize:
!    -----------
     Indicator = 0.0d0
     FNorm = 0.0d0
     
     IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN

     Metric = 0.0d0
     DO i=1,3
        Metric(i,i) = 1.0d0
     END DO

     SELECT CASE( CurrentCoordinateSystem() )
        CASE( AxisSymmetric, CylindricSymmetric )
           DIM = 3
        CASE DEFAULT
           DIM = CoordinateSystemDimension()
     END SELECT

     DOFs = DIM + 1
     IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs-1
!
!    Element nodal points:
!    ---------------------
     n = Element % TYPE % NumberOfNodes

     ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     ALLOCATE( NodalViscosity(n), NodalDensity(n), Basis(n), dBasisdx(n,3), &
        ddBasisddx(n,3,3), Velocity(3,n), Pressure(n), PrevVelo(3,n),  &
        PrevPres(n), Temperature(n), NodalForce(4,n), HeatCapacity(n), &
        ReferenceTemperature(n), HeatExpansionCoeff(n) )

!
!    Material parameters: density, viscosity, etc.
!    ----------------------------------------------
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
                  minv=1, maxv=Model % NumberOfMaterials )

     Material => Model % Materials(k) % Values

     NodalDensity(1:n) = ListGetReal( &
         Material, 'Density', n, Element % NodeIndexes, GotIt )

     NodalViscosity(1:n) = ListGetReal( &
         Material, 'Viscosity', n, Element % NodeIndexes, GotIt )

     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Equation', &
                      minv=1, maxv=Model % NumberOfEquations   )

!      Convect = ListGetLogical( Model % Equations(k) % Values, &
!                    'NS Convect', GotIt )
!      IF ( .NOT. GotIt ) Convect = .TRUE.
       Convect = .TRUE.
       
!    Elementwise nodal solution:
!    ---------------------------
     Velocity = 0.0d0
     DO k=1,DOFs-1
        Velocity(k,1:n) = Quant( DOFs*Perm(Element % NodeIndexes)-DOFs+k )
     END DO
     Pressure(1:n) = Quant( DOFs*Perm(Element % NodeIndexes) )

!
!    Check for time dep.
!    -------------------
     PrevPres(1:n)     = Pressure(1:n)
     PrevVelo(1:3,1:n) = Velocity(1:3,1:n)

     dt = Model % Solver % dt

     IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
        Var => VariableGet( Model % Variables, 'Flow Solution', .TRUE. )

        PrevVelo = 0.0d0
        DO k=1,DOFs-1
           PrevVelo(k,1:n) = &
              Var % PrevValues(DOFs*Var % Perm(Element % NodeIndexes)-DOFs+k,1)
        END DO
        PrevPres(1:n)=Var % PrevValues(DOFs*Var % Perm(Element % NodeIndexes),1)
     END IF


!
!    Check for compressible flow equations:
!    --------------------------------------
     Compressible = .FALSE.

!      IF (  ListGetString( Material, 'Compressibility Model', GotIt ) == &
!                       'perfect gas equation 1' ) THEN
! 
!         Compressible = .TRUE.
! 
!         Var => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )
!         IF ( ASSOCIATED( Var ) ) THEN
!            Temperature(1:n) = &
!                Var % Values( Var % Perm(Element % NodeIndexes) )
!         ELSE
!            Temperature(1:n) = ListGetReal( Material, &
!                'Reference Temperature',n,Element % NodeIndexes )
!         END IF
! 
!         SpecificHeatRatio = ListGetConstReal( Material, &
!                   'Specific Heat Ratio' )
! 
!         ReferencePressure = ListGetConstReal( Material, &
!                    'Reference Pressure' )
! 
!         HeatCapacity(1:n) = ListGetReal( Material, &
!                       'Heat Capacity',n,Element % NodeIndexes )
! 
!         NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
!               ( (SpecificHeatRatio - 1) * HeatCapacity(1:n) * Temperature(1:n) )
!      END IF
!
!    Body Forces:
!    ------------
!
     k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, &
       'Body Force', GotIt, 1, Model % NumberOfBodyForces )

     NodalForce = 0.0d0

     IF ( GotIt .AND. k > 0  ) THEN
!
!       Boussinesq approximation of heat expansion for
!       incompressible flow equations:
!
!       Density for the force term equals to
!
!       \rho = rho_0 (1-\beta(T-T_0)),
!
!       where \beta is the  heat expansion  coefficient,
!       T temperature and \rho_0 and T_0 correspond to
!       stress free state. Otherwise density is assumed
!       constant.
!       ----------------------------------------------
        IF (ListGetLogical(Model % BodyForces(k) % Values,'Boussinesq',GotIt)) THEN

           Var => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )
           IF ( ASSOCIATED( Var ) ) THEN
              Temperature(1:n) = &
                  Var % Values( Var % Perm(Element % NodeIndexes) )

              HeatExpansionCoeff(1:n) = ListGetReal( Material, &
                 'Heat Expansion Coefficient',n,Element % NodeIndexes )

              ReferenceTemperature(1:n) = ListGetReal( Material, &
                 'Reference Temperature',n,Element % NodeIndexes )

              Gravity => ListGetConstRealArray( Model % Constants, &
                             'Gravity' )

              k = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Equation', &
                        minv=1, maxv=Model % NumberOfEquations )

              IF ( ListGetLogical( Model % Equations(k) % Values, &
                            'Hydrostatic Pressure', GotIt) ) THEN
                 DO i=1,DIM
                    NodalForce(i,1:n) = ( 1 - HeatExpansionCoeff(1:n) * &
                       ( Temperature(1:n) - ReferenceTemperature(1:n) ) ) * &
                            Gravity(i,1) * Gravity(4,1)
                 END DO
              ELSE
                 DO i=1,DIM
                    NodalForce(i,1:n) = ( -HeatExpansionCoeff(1:n) * &
                       ( Temperature(1:n) - ReferenceTemperature(1:n) ) ) * &
                            Gravity(i,1) * Gravity(4,1)
                 END DO
              END IF
           END IF
        END IF

!
!       Given external force:
!       ---------------------
        NodalForce(1,1:n) = NodalForce(1,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Flow BodyForce 1', &
                  n, Element % NodeIndexes, GotIt )

        NodalForce(2,1:n) = NodalForce(2,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Flow BodyForce 2', &
                  n, Element % NodeIndexes, GotIt )

        NodalForce(3,1:n) = NodalForce(3,1:n) + ListGetReal( &
             Model % BodyForces(k) % Values, 'Flow BodyForce 3', &
                  n, Element % NodeIndexes, GotIt )
     END IF
!
!    Integrate square of residual over element:
!    ------------------------------------------
     ResidualNorm = 0.0d0
     Area = 0.0d0
!    Local Residual Vector     
     Rtemp = 0.0_dp
     NBasis = GetElementNOFNodes(Element)      
     
     IntegStuff = GaussPoints( Element )

     DO t=1,IntegStuff % n
        u = IntegStuff % u(t)
        v = IntegStuff % v(t)
        w = IntegStuff % w(t)

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
            Basis, dBasisdx, ddBasisddx, .TRUE. )

        IF ( CurrentCoordinateSystem() == Cartesian ) THEN
           s = IntegStuff % s(t) * detJ
        ELSE
           u = SUM( Basis(1:n) * Nodes % x(1:n) )
           v = SUM( Basis(1:n) * Nodes % y(1:n) )
           w = SUM( Basis(1:n) * Nodes % z(1:n) )

           CALL CoordinateSystemInfo( Metric, SqrtMetric, &
                      Symb, dSymb, u, v, w )
           s = IntegStuff % s(t) * detJ * SqrtMetric
        END IF

        Density   = SUM( NodalDensity(1:n)   * Basis(1:n) )
        Viscosity = SUM( NodalViscosity(1:n) * Basis(1:n) )
!
!       Residual of the navier-stokes equations:
!       or more generally:
!       ----------------------------------------------------------
        Residual = 0.0d0
        DO i=1,DIM
! 
           Temp = 0.0_dp
!          given force:
!          -------------
           Residual(i) = -Density * SUM( NodalForce(i,1:n) * Basis(1:n) )
! !            write (*,*) "INRES : given force              - ",t," -",i,"= ",Residual(i) - Temp
! !            Temp = Residual(i)
           
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!!!!!!!!   F = - rho NodalForce
!!!!!!!!   R = - Lambda0 * F
! ! ! !            Residual(i) = lambda0 * (-Density * SUM( NodalForce(i,1:n) * Basis(1:n) ))
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
!            IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!             + grad(p):
!             ----------
              Residual(i) = Residual(i) + SUM( Pressure(1:n) * dBasisdx(1:n,i) )
!            write (*,*) "INRES : grad(p)                  - ",t," -",i,"= ",Residual(i) - Temp
           Temp = Residual(i)
           

              DO j=1,DIM
!
!                - 2 ( \mu \epsilon^{ij} )_{,j}:
!                -------------------------------
                 Residual(i) = Residual(i) - Viscosity * &
                     SUM( Velocity(i,1:n) * ddBasisddx(1:n,j,j) )        ! laplace= div(grad U) Ui,jj

! !                  Residual(i) = Residual(i) - &
! !                       SUM( NodalViscosity(1:n) * dBasisdx(1:n,j) ) * &   ! si viscosité non constante
! !                           SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )    

                  Residual(i) = Residual(i) - Viscosity * &
                      SUM( Velocity(j,1:n) * ddBasisddx(1:n,i,j) )       ! div((grad U)T) ij Uj,ij

! !                   Residual(i) = Residual(i) - &
! !                       SUM( NodalViscosity(1:n) * dBasisdx(1:n,j) ) * &   ! si viscosité non constante
! !                           SUM( Velocity(j,1:n) * dBasisdx(1:n,i) )

!                   IF ( Compressible ) THEN
! !
! !                    + (2/3) grad(\mu div(u)):
! !                    -------------------------
!                      Residual(i) = Residual(i) + &
!                         Viscosity * ( 2.0d0 / 3.0d0 ) * &
!                            SUM( Velocity(j,1:n) * ddBasisddx(1:n,j,i) )
! 
!                      Residual(i) = Residual(i) + &
!                          SUM( NodalViscosity(1:n) * dBasisdx(1:n,i) ) * &
!                              SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )
! 
!                   END IF
              END DO
!            write (*,*) "INRES : -2(\mu\epsilon^{ij})_{,j}- ",t," -",i,"= ",Residual(i) - Temp
           Temp = Residual(i)
           

!               IF ( Convect ) THEN
!
!                + \rho * (@u/@t + u.grad(u)):
!                -----------------------------
! ! ! ! !                  Residual(i) = Residual(i) + Density *  &
! ! ! ! !                      SUM((Velocity(i,1:n)-PrevVelo(i,1:n))*Basis(1:n)) / dt

                 DO j=1,DIM
                    Residual(i) = Residual(i) + &
                        Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
                                  SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )   ! Ui Ni Ui Nipwd,j
                 END DO
!            write (*,*) "INRES : \rho u.grad(u)           - ",t," -",i,"= ",Residual(i) - Temp
           Temp = Residual(i)
           
                 
!               END IF
!            ELSE
! !             + g^{ij}p_{,j}:
! !             ---------------
!               DO j=1,DIM
!                  Residual(i) = Residual(i) + Metric(i,j) * &
!                       SUM( Pressure(1:n) * dBasisdx(1:n,i) )
!               END DO
! 
! !             - g^{jk} (\mu u^i_{,k})_{,j}):
! !             ------------------------------
!               DO j=1,DIM
!                  DO k=1,DIM
!                     Residual(i) = Residual(i) -   &
!                          Metric(j,k) * Viscosity * &
!                          SUM( Velocity(i,1:n) * ddBasisddx(1:n,j,k) )
! 
!                     DO l=1,DIM
!                        Residual(i) = Residual(i) +  &
!                             Metric(j,k) * Viscosity * Symb(j,k,l) * &
!                             SUM( Velocity(i,1:n) * dBasisdx(1:n,l) )
! 
!                        Residual(i) = Residual(i) -  &
!                             Metric(j,k) * Viscosity * Symb(l,j,i) * &
!                             SUM( Velocity(l,1:n) * dBasisdx(1:n,k) )
! 
!                        Residual(i) = Residual(i) -  &
!                             Metric(j,k) * Viscosity * Symb(l,k,i) * &
!                             SUM( Velocity(l,1:n) * dBasisdx(1:n,j) )
! 
!                        Residual(i) = Residual(i) -  &
!                             Metric(j,k) * Viscosity * dSymb(l,j,i,k) * &
!                             SUM( Velocity(l,1:n) * Basis(1:n) )
! 
!                        DO m=1,DIM
!                           Residual(i) = Residual(i) - Metric(j,k) * Viscosity *&
!                                   Symb(m,k,i) * Symb(l,j,m) * &
!                                         SUM( Velocity(l,1:n) * Basis(1:n) )
! 
!                           Residual(i) = Residual(i) + Metric(j,k) * Viscosity *&
!                                   Symb(j,k,m) * Symb(l,m,i) * &
!                                         SUM( Velocity(l,1:n) * Basis(1:n) )
!                        END DO
!                     END DO
!                  END DO
!               END DO
! 
! !             - g^{ik} (\mu u^j_{,k})_{,j}):
! !             ------------------------------
!               DO j=1,DIM
!                  DO k=1,DIM
!                     Residual(i) = Residual(i) -   &
!                          Metric(i,k) * Viscosity * &
!                          SUM( Velocity(j,1:n) * ddBasisddx(1:n,j,k) )
! 
!                     DO l=1,DIM
!                        Residual(i) = Residual(i) +  &
!                             Metric(i,k) * Viscosity * Symb(j,k,l) * &
!                             SUM( Velocity(j,1:n) * dBasisdx(1:n,l) )
! 
!                        Residual(i) = Residual(i) -  &
!                             Metric(i,k) * Viscosity * Symb(l,j,j) * &
!                             SUM( Velocity(l,1:n) * dBasisdx(1:n,k) )
! 
!                        Residual(i) = Residual(i) -  &
!                             Metric(i,k) * Viscosity * Symb(l,k,j) * &
!                             SUM( Velocity(l,1:n) * dBasisdx(1:n,j) )
! 
!                        Residual(i) = Residual(i) -  &
!                             Metric(i,k) * Viscosity * dSymb(l,j,j,k) * &
!                             SUM( Velocity(l,1:n) * Basis(1:n) )
! 
!                        DO m=1,DIM
!                           Residual(i) = Residual(i) - Metric(i,k) * Viscosity *&
!                                   Symb(m,k,j) * Symb(l,j,m) * &
!                                         SUM( Velocity(l,1:n) * Basis(1:n) )
! 
!                           Residual(i) = Residual(i) + Metric(i,k) * Viscosity *&
!                                   Symb(j,k,m) * Symb(l,m,j) * &
!                                         SUM( Velocity(l,1:n) * Basis(1:n) )
!                        END DO
!                     END DO
!                  END DO
!               END DO
! 
!               IF ( Convect ) THEN
! !
! !                + \rho * (@u/@t + u^j u^i_{,j}):
! !                --------------------------------
!                  Residual(i) = Residual(i) + Density *  &
!                      SUM((Velocity(i,1:n)-PrevVelo(i,1:n))*Basis(1:n)) / dt
! 
!                  DO j=1,DIM
!                     Residual(i) = Residual(i) + &
!                          Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
!                          SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )
! 
!                     DO k=1,DIM
!                        Residual(i) = Residual(i) + &
!                             Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
!                             Symb(j,k,i) * SUM( Velocity(k,1:n) * Basis(1:n) )
!                     END DO
!                  END DO
!               END IF
!            END IF
!            write (*,*) "INRES : Residual      - ",t," -",i,"=",Residual(i)
        END DO

!
!       Continuity equation:
!       --------------------
!         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!
!          + \rho * div(u):
!          ----------------
           DO j=1,DIM
              Residual(DIM+1) = Residual(DIM+1) + &
                   Density * SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )    ! Ui,i
           END DO
!            write (*,*) "INRES : \rho * div(u)            - ",t," -",i,"=",Residual(i) - Temp
           Temp = 0.0_dp
           

!            IF ( Compressible ) THEN
! !
! !             + u.grad(\rho):
! !             ----------------
!               DO j=1,DIM
!                  Residual(DIM+1) = Residual(DIM+1) + &
!                       SUM( Velocity(j,1:n) * Basis(1:n) ) *  &
!                            SUM( NodalDensity(1:n) * dBasisdx(1:n,j) ) 
!               END DO
!            END IF
!         ELSE
! !
! !          + \rho * u^j_{,j}:
! !          ------------------
!            DO j=1,DIM
!               Residual(DIM+1) = Residual(DIM+1) + &
!                    Density * SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )
! 
!               DO k=1,DIM
!                  Residual(DIM+1) = Residual(DIM+1) + Density * &
!                       Symb(k,j,j) * SUM( Velocity(k,1:n) * Basis(1:n) )
!               END DO
!            END DO
! 
!            IF ( Compressible ) THEN
! !
! !             + u^j \rho_{,j}:
! !             ----------------
!               DO j=1,DIM
!                  Residual(DIM+1) = Residual(DIM+1) + &
!                       SUM( Velocity(j,1:n) * Basis(1:n) ) *  &
!                       SUM( NodalDensity(1:n) * dBasisdx(1:n,j) ) 
!               END DO
!            END IF
!         END IF

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        
        DO i=1,DIM
           FNorm = FNorm + s * (Density * SUM(NodalForce(i,1:n)*Basis(1:n))**2)
        END DO 
        Area = Area + s
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
!         
           ResidualNorm = ResidualNorm + &
             s * (Element % hK**2 * SUM(Residual(1:dim)**2) + Residual(dim+1)**2 )
!              
!         ELSE
!            CALL InvertMatrix( Metric,3 )
!            DO i=1,dim
!               DO j=1,dim
!                  ResidualNorm = ResidualNorm + &
!                     s * Element % hK **2 * Metric(i,j) * Residual(i) * Residual(j)
!               END DO
!            END DO
!            ResidualNorm = ResidualNorm + s * Residual(dim+1)**2
!         END IF

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                                  _     _       
!     /\                          | |   | |      
!    /  \   ___ ___  ___ _ __ ___ | |__ | |_   _ 
!   / /\ \ / __/ __|/ _ \ '_ ` _ \| '_ \| | | | |
!  / ____ \\__ \__ \  __/ | | | | | |_) | | |_| |
! /_/    \_\___/___/\___|_| |_| |_|_.__/|_|\__, |
!                                           __/ |
!                                          |___/ 
           
           Temp = DSQRT(SUM(Residual(1:dim)**2) + Residual(dim+1)**2)
!            write (*,*) "INRES : ||Residual||  - ",t," -",Temp
!            write (*,*) "INRES :  s            - ",t," -",s

      c = dim + 1
      DO p = 1 , NBasis                                       ! boucle sur les noeuds
        Rlocal => Rtemp( c*(p-1)+1 : c*(p-1)+c  )             !   Pointeur
        DO i = 1, dim + 1                                     !   Boucle sur les DOFs
!         ICI  s = GaussWeight * detJ        
          Rlocal(i) = Rlocal(i) + s * Residual(i) * Basis(p)    !     wi * Ri * Ni ??????
        END DO                                                !   fin Boucle sur les DOFs
      END DO                                                  ! fin boucle sur les noeuds 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -        

!------------------------------------------------------------------------------
     END DO                       ! fin de la boucle sur les points d'intégration
!******************************************************************************   

     FNorm = Area * FNorm
     Indicator = ResidualNorm
     


     DEALLOCATE( NodalViscosity, NodalDensity, Basis, dBasisdx,           &
        ddBasisddx, Velocity, Pressure, PrevVelo, PrevPres, Temperature,   &
        NodalForce, HeatCapacity, ReferenceTemperature, HeatExpansionCoeff,&
        Nodes % x, Nodes % y, Nodes % z )
!------------------------------------------------------------------------------
  END SUBROUTINE HMLocalFlowInsideResidual
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------





! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! !  _____           _     _             _ 
! ! ! ! ! ! ! ! |  __ \         (_)   | |           | |
! ! ! ! ! ! ! ! | |__) |___  ___ _  __| |_   _  __ _| |
! ! ! ! ! ! ! ! |  _  // _ \/ __| |/ _` | | | |/ _` | |
! ! ! ! ! ! ! ! | | \ \  __/\__ \ | (_| | |_| | (_| | |
! ! ! ! ! ! ! ! |_|  \_\___||___/_|\__,_|\__,_|\__,_|_|
! ! ! ! ! ! ! !                                        
! ! ! ! ! ! ! ! BY
! ! ! ! ! ! ! !  ______ _                           _    
! ! ! ! ! ! ! ! |  ____| |                         | |   
! ! ! ! ! ! ! ! | |__  | | ___ _ __ ___   ___ _ __ | |_  
! ! ! ! ! ! ! ! |  __| | |/ _ \ '_ ` _ \ / _ \ '_ \| __| 
! ! ! ! ! ! ! ! | |____| |  __/ | | | | |  __/ | | | |_  
! ! ! ! ! ! ! ! |______|_|\___|_| |_| |_|\___|_| |_|\__|       
! ! ! ! ! ! ! !                                        
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! !> Compute the residual of the Navier-Stokes equation for the boundary elements.
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !   FUNCTION FlowBoundaryResidual( Model, Edge, Mesh, &
! ! ! ! ! ! !         Quant, Perm, Gnorm ) RESULT( Indicator )
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !      USE DefUtils
! ! ! ! ! ! !      IMPLICIT NONE
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !      TYPE(Model_t) :: Model
! ! ! ! ! ! !      INTEGER :: Perm(:)
! ! ! ! ! ! !      REAL(KIND=dp) :: Quant(:), Indicator(2), Gnorm
! ! ! ! ! ! !      TYPE( Mesh_t ), POINTER    :: Mesh
! ! ! ! ! ! !      TYPE( Element_t ), POINTER :: Edge
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(Nodes_t) :: Nodes, EdgeNodes
! ! ! ! ! ! !      TYPE(Element_t), POINTER :: Element
! ! ! ! ! ! ! 
! ! ! ! ! ! !      INTEGER :: i,k,n,l,t,bc,DIM,DOFs,Pn,En
! ! ! ! ! ! !      LOGICAL :: stat, GotIt, Compressible
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp) :: Grad(3,3), Grad1(3,3), Stress(3,3), Normal(3), ForceSolved(3), &
! ! ! ! ! ! !                       EdgeLength,u, v, w, s, detJ
! ! ! ! ! ! !      REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE :: EdgeBasis(:), dEdgeBasisdx(:,:), Basis(:),dBasisdx(:,:)
! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE ::  x(:), y(:), z(:), ExtPressure(:)
! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE :: Temperature(:), Tension(:), SlipCoeff(:,:)
! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE :: Velocity(:,:), Pressure(:), Force(:,:), NodalViscosity(:)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp) :: Residual(3), ResidualNorm, Viscosity, Slip, Dir(3)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(Variable_t), POINTER  :: TempSol
! ! ! ! ! ! !      TYPE(ValueList_t), POINTER :: Material
! ! ! ! ! ! !      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !    Initialize:
! ! ! ! ! ! ! !    -----------
! ! ! ! ! ! !      Indicator = 0.0d0
! ! ! ! ! ! !      Gnorm     = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Metric = 0.0d0
! ! ! ! ! ! !      DO i=1,3
! ! ! ! ! ! !         Metric(i,i) = 1.0d0
! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !      SELECT CASE( CurrentCoordinateSystem() )
! ! ! ! ! ! !         CASE( AxisSymmetric, CylindricSymmetric )
! ! ! ! ! ! !            DIM = 3
! ! ! ! ! ! !         CASE DEFAULT
! ! ! ! ! ! !            DIM = CoordinateSystemDimension()
! ! ! ! ! ! !      END SELECT
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DOFs = DIM + 1
! ! ! ! ! ! !      IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs-1
! ! ! ! ! ! ! !    
! ! ! ! ! ! ! !    --------------------------------------------------
! ! ! ! ! ! !      Element => Edge % BoundaryInfo % Left
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IF ( .NOT. ASSOCIATED( Element ) ) THEN
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Element => Edge % BoundaryInfo % Right
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ELSE IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) THEN
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Element => Edge % BoundaryInfo % Right
! ! ! ! ! ! ! 
! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IF ( .NOT. ASSOCIATED( Element ) ) RETURN
! ! ! ! ! ! !      IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN
! ! ! ! ! ! ! 
! ! ! ! ! ! !      En = Edge % TYPE % NumberOfNodes
! ! ! ! ! ! !      Pn = Element % TYPE % NumberOfNodes
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
! ! ! ! ! ! !      EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
! ! ! ! ! ! !      EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ALLOCATE( Nodes % x(Pn), Nodes % y(Pn), Nodes % z(Pn) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
! ! ! ! ! ! !      Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
! ! ! ! ! ! !      Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ALLOCATE( EdgeBasis(En), dEdgeBasisdx(En,3), Basis(Pn), dBasisdx(Pn,3), &
! ! ! ! ! ! !       x(En), y(En), z(En), ExtPressure(En), Temperature(Pn), Tension(En),    &
! ! ! ! ! ! !       SlipCoeff(3,En), Velocity(3,Pn), Pressure(Pn), Force(3,En), NodalViscosity(En) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DO l = 1,En
! ! ! ! ! ! !        DO k = 1,Pn
! ! ! ! ! ! !           IF ( Edge % NodeIndexes(l) == Element % NodeIndexes(k) ) THEN
! ! ! ! ! ! !              x(l) = Element % TYPE % NodeU(k)
! ! ! ! ! ! !              y(l) = Element % TYPE % NodeV(k)
! ! ! ! ! ! !              z(l) = Element % TYPE % NodeW(k)
! ! ! ! ! ! !              EXIT
! ! ! ! ! ! !           END IF
! ! ! ! ! ! !        END DO
! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Integrate square of residual over boundary element:
! ! ! ! ! ! ! !    ---------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Indicator    = 0.0d0
! ! ! ! ! ! !      EdgeLength   = 0.0d0
! ! ! ! ! ! !      ResidualNorm = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DO bc=1,Model % NumberOfBCs
! ! ! ! ! ! !         IF ( Edge % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !       IF ( .NOT. ListGetLogical( Model % BCs(bc) % Values, &
! ! ! ! ! ! ! !                 'Flow Force BC', gotIt ) ) CYCLE
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Get material parameters:
! ! ! ! ! ! ! !       ------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! !         k = ListGetInteger(Model % Bodies(Element % BodyId) % Values,'Material', &
! ! ! ! ! ! !                      minv=1, maxv=Model % NumberOfMaterials )
! ! ! ! ! ! !         Material => Model % Materials(k) % Values
! ! ! ! ! ! ! 
! ! ! ! ! ! !         NodalViscosity(1:En) = ListGetReal( Material, &
! ! ! ! ! ! !                  'Viscosity', En, Edge % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Compressible = .FALSE.
! ! ! ! ! ! !         IF ( ListGetString( Material, 'Compressibility Model', GotIt ) == &
! ! ! ! ! ! !                'perfect gas equation 1' ) Compressible = .TRUE.
! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !       Given traction:
! ! ! ! ! ! ! !       ---------------
! ! ! ! ! ! !         Force = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Force(1,1:En) = ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !             'Pressure 1', En, Edge % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Force(2,1:En) = ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !             'Pressure 2', En, Edge % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Force(3,1:En) = ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !             'Pressure 3', En, Edge % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Force in normal direction:
! ! ! ! ! ! ! !       ---------------------------
! ! ! ! ! ! !         ExtPressure(1:En) = ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !           'External Pressure', En, Edge % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Slip BC condition:
! ! ! ! ! ! ! !       ------------------
! ! ! ! ! ! !         SlipCoeff = 0.0d0
! ! ! ! ! ! !         SlipCoeff(1,1:En) =  ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !              'Slip Coefficient 1',En,Edge % NodeIndexes,GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         SlipCoeff(2,1:En) =  ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !              'Slip Coefficient 2',En,Edge % NodeIndexes,GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         SlipCoeff(3,1:En) =  ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !              'Slip Coefficient 3',En,Edge % NodeIndexes,GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Surface tension induced by temperature gradient (or otherwise):
! ! ! ! ! ! ! !       ---------------------------------------------------------------
! ! ! ! ! ! !         TempSol => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         IF ( ASSOCIATED( TempSol ) ) THEN
! ! ! ! ! ! !           Tension(1:En) = ListGetReal( Model % BCs(bc) % Values, &
! ! ! ! ! ! !            'Surface Tension Expansion Coefficient',En,Edge % NodeIndexes,gotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( gotIt ) THEN
! ! ! ! ! ! !               DO n=1,En
! ! ! ! ! ! !                  k = TempSol % Perm( Edge % NodeIndexes(n) )
! ! ! ! ! ! !                  IF (k>0) Tension(n) = 1.0d0 - Tension(n) * TempSol % Values(k)
! ! ! ! ! ! !               END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !               Tension(1:En) = Tension(1:En) * ListGetReal( &
! ! ! ! ! ! !                  Model % BCs(bc) % Values,'Surface Tension Coefficient', &
! ! ! ! ! ! !                                En, Edge % NodeIndexes ) 
! ! ! ! ! ! !            ELSE
! ! ! ! ! ! !               Tension(1:En) = ListGetReal( &
! ! ! ! ! ! !                   Model % BCs(bc) % Values,'Surface Tension Coefficient', &
! ! ! ! ! ! !                          En, Edge % NodeIndexes,gotIt ) 
! ! ! ! ! ! !            END IF
! ! ! ! ! ! !         ELSE
! ! ! ! ! ! !            Tension(1:En) = ListGetReal( &
! ! ! ! ! ! !                Model % BCs(bc) % Values,'Surface Tension Coefficient', &
! ! ! ! ! ! !                       En, Edge % NodeIndexes,gotIt ) 
! ! ! ! ! ! !         END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       If dirichlet BC for velocity in any direction given,
! ! ! ! ! ! ! !       nullify force in that directon:
! ! ! ! ! ! ! !       ------------------------------------------------------------------
! ! ! ! ! ! !         Dir = 1
! ! ! ! ! ! !         s = ListGetConstReal( Model % BCs(bc) % Values, 'Velocity 1', GotIt )
! ! ! ! ! ! !         IF ( GotIt ) Dir(1) = 0
! ! ! ! ! ! ! 
! ! ! ! ! ! !         s = ListGetConstReal( Model % BCs(bc) % Values, 'Velocity 2', GotIt )
! ! ! ! ! ! !         IF ( GotIt ) Dir(2) = 0
! ! ! ! ! ! ! 
! ! ! ! ! ! !         s = ListGetConstReal( Model % BCs(bc) % Values, 'Velocity 3', GotIt )
! ! ! ! ! ! !         IF ( GotIt ) Dir(3) = 0
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Elementwise nodal solution:
! ! ! ! ! ! ! !       ---------------------------
! ! ! ! ! ! !         Velocity = 0.0d0
! ! ! ! ! ! !         DO k=1,DOFs-1
! ! ! ! ! ! !            Velocity(k,1:Pn) = Quant(DOFs*Perm(Element % NodeIndexes)-DOFs + k)
! ! ! ! ! ! !         END DO
! ! ! ! ! ! !         Pressure(1:Pn) = Quant( DOFs*Perm(Element % NodeIndexes) )
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !       do the integration:
! ! ! ! ! ! ! !       -------------------
! ! ! ! ! ! !         EdgeLength   = 0.0d0
! ! ! ! ! ! !         ResidualNorm = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !         IntegStuff = GaussPoints( Edge )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         DO t=1,IntegStuff % n
! ! ! ! ! ! !            u = IntegStuff % u(t)
! ! ! ! ! ! !            v = IntegStuff % v(t)
! ! ! ! ! ! !            w = IntegStuff % w(t)
! ! ! ! ! ! ! 
! ! ! ! ! ! !            stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
! ! ! ! ! ! !                EdgeBasis, dEdgeBasisdx )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !               s = IntegStuff % s(t) * detJ
! ! ! ! ! ! !            ELSE
! ! ! ! ! ! !               u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
! ! ! ! ! ! !               v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
! ! ! ! ! ! !               w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
! ! ! ! ! ! !       
! ! ! ! ! ! !               CALL CoordinateSystemInfo( Metric, SqrtMetric, &
! ! ! ! ! ! !                           Symb, dSymb, u, v, w )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               s = IntegStuff % s(t) * detJ * SqrtMetric
! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Normal = NormalVector( Edge, EdgeNodes, u, v, .TRUE. )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            u = SUM( EdgeBasis(1:En) * x(1:En) )
! ! ! ! ! ! !            v = SUM( EdgeBasis(1:En) * y(1:En) )
! ! ! ! ! ! !            w = SUM( EdgeBasis(1:En) * z(1:En) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
! ! ! ! ! ! !               Basis, dBasisdx )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Viscosity = SUM( NodalViscosity(1:En) * EdgeBasis(1:En) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Residual = 0.0d0
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          Given force at the integration point:
! ! ! ! ! ! ! !          -------------------------------------
! ! ! ! ! ! !            Residual = Residual + MATMUL( Force(:,1:En), EdgeBasis(1:En) ) - &
! ! ! ! ! ! !                  SUM( ExtPressure(1:En) * EdgeBasis(1:En) ) * Normal
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          Slip velocity BC:
! ! ! ! ! ! ! !          -----------------
! ! ! ! ! ! !            DO i=1,DIM
! ! ! ! ! ! !               Slip = SUM( SlipCoeff(i,1:En) * EdgeBasis(1:En) )
! ! ! ! ! ! !               Residual(i) = Residual(i) - &
! ! ! ! ! ! !                    Slip * SUM( Velocity(i,1:Pn) * Basis(1:Pn) )
! ! ! ! ! ! !            END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          Tangential tension force:
! ! ! ! ! ! ! !          -------------------------
! ! ! ! ! ! !            DO i=1,DIM
! ! ! ! ! ! !               Residual(i) = Residual(i) + &
! ! ! ! ! ! !                    SUM( dEdgeBasisdx(1:En,i) * Tension(1:En) )
! ! ! ! ! ! !            END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          Force given by the computed solution:
! ! ! ! ! ! ! !          -------------------------------------
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          Stress tensor on the boundary:
! ! ! ! ! ! ! !          ------------------------------
! ! ! ! ! ! !            Grad = MATMUL( Velocity(:,1:Pn), dBasisdx(1:Pn,:) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
! ! ! ! ! ! !               Grad1 = Grad
! ! ! ! ! ! !               DO i=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     DO l=1,DIM
! ! ! ! ! ! !                        Grad1(i,k) = Grad1(i,k) - &
! ! ! ! ! ! !                           Symb(k,l,i) * SUM ( Velocity(l,1:Pn) * Basis(1:Pn) )
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !               Grad = 0.0d0
! ! ! ! ! ! !               DO i=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     DO l=1,DIM
! ! ! ! ! ! !                        Grad(i,k) = Grad(i,k) + Metric(k,l) * Grad1(i,l)
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Stress = Viscosity * ( Grad + TRANSPOSE(Grad) )
! ! ! ! ! ! !            Stress = Stress - Metric * SUM( Pressure(1:Pn) * Basis(1:Pn) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( Compressible ) THEN
! ! ! ! ! ! !               IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !                  DO i=1,DIM
! ! ! ! ! ! !                     DO k=1,DIM
! ! ! ! ! ! !                        Stress(i,i) = Stress(i,i) - &
! ! ! ! ! ! !                            (2.0d0/3.0d0) * Viscosity * Grad(k,k)
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               ELSE
! ! ! ! ! ! !                  DO i=1,DIM
! ! ! ! ! ! !                     DO k=1,DIM
! ! ! ! ! ! !                        DO l=1,DIM
! ! ! ! ! ! !                           Stress(i,k) = Stress(i,k) - &
! ! ! ! ! ! !                              Metric(i,k) * (2.0d0/3.0d0) * Viscosity * Grad(l,l)
! ! ! ! ! ! !                        END DO
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END IF
! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !            ForceSolved = MATMUL(Stress,Normal)
! ! ! ! ! ! !            Residual = Residual - ForceSolved * Dir
! ! ! ! ! ! ! 
! ! ! ! ! ! !            EdgeLength = EdgeLength + s
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !               Gnorm = Gnorm + s * SUM( ForceSolved**2 )
! ! ! ! ! ! !               ResidualNorm = ResidualNorm + s * SUM( Residual(1:DIM) ** 2 )
! ! ! ! ! ! !            ELSE
! ! ! ! ! ! !               CALL InvertMatrix( Metric,3 )
! ! ! ! ! ! !               DO i=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     ResidualNorm = ResidualNorm + &
! ! ! ! ! ! !                             s * Metric(i,k) * Residual(i) * Residual(k)
! ! ! ! ! ! !                     Gnorm = GNorm + s * Metric(i,k) * &
! ! ! ! ! ! !                                         ForceSolved(i) * ForceSolved(k)
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END IF
! ! ! ! ! ! !         END DO
! ! ! ! ! ! !         EXIT
! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IF ( CoordinateSystemDimension() == 3 ) EdgeLength = SQRT(EdgeLength)
! ! ! ! ! ! !      Indicator = EdgeLength * ResidualNorm
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
! ! ! ! ! ! !      DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DEALLOCATE( EdgeBasis, dEdgeBasisdx, Basis, dBasisdx, x, y, z,   &
! ! ! ! ! ! !       ExtPressure, Temperature, Tension,SlipCoeff, Velocity, Pressure, &
! ! ! ! ! ! !       Force, NodalViscosity )
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !   END FUNCTION FlowBoundaryResidual
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    
! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! !> Compute the residual of the Navier-Stokes equation for the edge elements.
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !   FUNCTION FlowEdgeResidual( Model,Edge,Mesh,Quant,Perm ) RESULT( Indicator )
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !      USE DefUtils
! ! ! ! ! ! !      IMPLICIT NONE
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(Model_t) :: Model
! ! ! ! ! ! !      INTEGER :: Perm(:)
! ! ! ! ! ! !      REAL(KIND=dp) :: Quant(:), Indicator(2)
! ! ! ! ! ! !      TYPE( Mesh_t ), POINTER    :: Mesh
! ! ! ! ! ! !      TYPE( Element_t ), POINTER :: Edge
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(Nodes_t) :: Nodes, EdgeNodes
! ! ! ! ! ! !      TYPE(Element_t), POINTER :: Element
! ! ! ! ! ! ! 
! ! ! ! ! ! !      INTEGER :: i,j,k,l,n,t,DIM,DOFs,En,Pn
! ! ! ! ! ! !      LOGICAL :: stat, GotIt
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
! ! ! ! ! ! !      REAL(KIND=dp) :: Stress(3,3,2), Jump(3), Viscosity
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp) :: Grad(3,3), Grad1(3,3), Normal(3)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE :: NodalViscosity(:), x(:), y(:), z(:), &
! ! ! ! ! ! !                EdgeBasis(:), Basis(:), dBasisdx(:,:)
! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE :: Velocity(:,:), Pressure(:)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp) :: ResidualNorm, EdgeLength, u, v, w, s, detJ
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !    Initialize:
! ! ! ! ! ! ! !    -----------
! ! ! ! ! ! !      SELECT CASE( CurrentCoordinateSystem() )
! ! ! ! ! ! !         CASE( AxisSymmetric, CylindricSymmetric )
! ! ! ! ! ! !            DIM = 3
! ! ! ! ! ! !         CASE DEFAULT
! ! ! ! ! ! !            DIM = CoordinateSystemDimension()
! ! ! ! ! ! !      END SELECT
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DOFs = DIM + 1
! ! ! ! ! ! !      IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs - 1
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Metric = 0.0d0
! ! ! ! ! ! !      DO i = 1,3
! ! ! ! ! ! !         Metric(i,i) = 1.0d0
! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Grad = 0.0d0
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    ---------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Element => Edge % BoundaryInfo % Left
! ! ! ! ! ! !      n = Element % TYPE % NumberOfNodes
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Element => Edge % BoundaryInfo % Right
! ! ! ! ! ! !      n = MAX( n, Element % TYPE % NumberOfNodes )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      En = Edge % TYPE % NumberOfNodes
! ! ! ! ! ! !      ALLOCATE( EdgeNodes % x(En), EdgeNodes % y(En), EdgeNodes % z(En) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      EdgeNodes % x = Mesh % Nodes % x(Edge % NodeIndexes)
! ! ! ! ! ! !      EdgeNodes % y = Mesh % Nodes % y(Edge % NodeIndexes)
! ! ! ! ! ! !      EdgeNodes % z = Mesh % Nodes % z(Edge % NodeIndexes)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ALLOCATE( NodalViscosity(En), x(En), y(En), z(En), EdgeBasis(En), &
! ! ! ! ! ! !            Basis(n), dBasisdx(n,3), Velocity(3,n), Pressure(n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !    Integrate square of jump over edge:
! ! ! ! ! ! ! !    ------------------------------------
! ! ! ! ! ! !      ResidualNorm = 0.0d0
! ! ! ! ! ! !      EdgeLength   = 0.0d0
! ! ! ! ! ! !      Indicator    = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IntegStuff = GaussPoints( Edge )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DO t=1,IntegStuff % n
! ! ! ! ! ! ! 
! ! ! ! ! ! !         u = IntegStuff % u(t)
! ! ! ! ! ! !         v = IntegStuff % v(t)
! ! ! ! ! ! !         w = IntegStuff % w(t)
! ! ! ! ! ! ! 
! ! ! ! ! ! !         stat = ElementInfo( Edge, EdgeNodes, u, v, w, detJ, &
! ! ! ! ! ! !              EdgeBasis, dBasisdx )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Normal = NormalVector( Edge, EdgeNodes, u, v, .FALSE. )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !            s = IntegStuff % s(t) * detJ
! ! ! ! ! ! !         ELSE
! ! ! ! ! ! !            u = SUM( EdgeBasis(1:En) * EdgeNodes % x(1:En) )
! ! ! ! ! ! !            v = SUM( EdgeBasis(1:En) * EdgeNodes % y(1:En) )
! ! ! ! ! ! !            w = SUM( EdgeBasis(1:En) * EdgeNodes % z(1:En) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            CALL CoordinateSystemInfo( Metric, SqrtMetric, &
! ! ! ! ! ! !                        Symb, dSymb, u, v, w )
! ! ! ! ! ! !            s = IntegStuff % s(t) * detJ * SqrtMetric
! ! ! ! ! ! !         END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Stress = 0.0d0
! ! ! ! ! ! !         DO i = 1,2
! ! ! ! ! ! !            IF ( i==1 ) THEN
! ! ! ! ! ! !               Element => Edge % BoundaryInfo % Left
! ! ! ! ! ! !            ELSE
! ! ! ! ! ! !               Element => Edge % BoundaryInfo % Right
! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) CYCLE
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Pn = Element % TYPE % NumberOfNodes
! ! ! ! ! ! !            Nodes % x(1:Pn) = Mesh % Nodes % x(Element % NodeIndexes)
! ! ! ! ! ! !            Nodes % y(1:Pn) = Mesh % Nodes % y(Element % NodeIndexes)
! ! ! ! ! ! !            Nodes % z(1:Pn) = Mesh % Nodes % z(Element % NodeIndexes)
! ! ! ! ! ! ! 
! ! ! ! ! ! !            DO j = 1,En
! ! ! ! ! ! !               DO k = 1,Pn
! ! ! ! ! ! !                  IF ( Edge % NodeIndexes(j) == Element % NodeIndexes(k) ) THEN
! ! ! ! ! ! !                     x(j) = Element % TYPE % NodeU(k)
! ! ! ! ! ! !                     y(j) = Element % TYPE % NodeV(k)
! ! ! ! ! ! !                     z(j) = Element % TYPE % NodeW(k)
! ! ! ! ! ! !                     EXIT
! ! ! ! ! ! !                  END IF
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !            u = SUM( EdgeBasis(1:En) * x(1:En) )
! ! ! ! ! ! !            v = SUM( EdgeBasis(1:En) * y(1:En) )
! ! ! ! ! ! !            w = SUM( EdgeBasis(1:En) * z(1:En) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
! ! ! ! ! ! !                Basis, dBasisdx )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            k = ListGetInteger( Model % Bodies( Element % BodyId) % Values, 'Material', &
! ! ! ! ! ! !                             minv=1, maxv=Model % NumberOfMaterials )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            NodalViscosity(1:En) = ListGetReal( &
! ! ! ! ! ! !                Model % Materials(k) % Values, 'Viscosity', &
! ! ! ! ! ! !                     En, Edge % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Viscosity = SUM( NodalViscosity(1:En) * EdgeBasis(1:En) )
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          Elementwise nodal solution:
! ! ! ! ! ! ! !          ---------------------------
! ! ! ! ! ! !            Velocity = 0.0d0
! ! ! ! ! ! !            DO k=1,DOFs-1
! ! ! ! ! ! !               Velocity(k,1:Pn) = Quant(DOFs*Perm(Element % NodeIndexes)-DOFs+k)
! ! ! ! ! ! !            END DO
! ! ! ! ! ! !            Pressure(1:Pn) = Quant( DOFs*Perm(Element % NodeIndexes) )
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          Stress tensor on the edge:
! ! ! ! ! ! ! !          --------------------------
! ! ! ! ! ! !            Grad = MATMUL( Velocity(:,1:Pn), dBasisdx(1:Pn,:) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
! ! ! ! ! ! !               Grad1 = Grad
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     DO l=1,DIM
! ! ! ! ! ! !                        Grad1(j,k) = Grad1(j,k) - &
! ! ! ! ! ! !                           Symb(k,l,j) * SUM ( Velocity(l,1:Pn) * Basis(1:Pn) )
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !               Grad = 0.0d0
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     DO l=1,DIM
! ! ! ! ! ! !                        Grad(j,k) = Grad(j,k) + Metric(k,l) * Grad1(j,l)
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Stress(:,:,i) = Viscosity * ( Grad + TRANSPOSE(Grad) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  Stress(j,j,i) = Stress(j,j,i) - SUM( Pressure(1:Pn) * Basis(1:Pn))
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     Stress(j,j,i) = Stress(j,j,i) - (2.0d0/3.0d0)*Viscosity*Grad(k,k)
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            ELSE
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     Stress(j,k,i) = Stress(j,k,i) - &
! ! ! ! ! ! !                            Metric(j,k) * SUM( Pressure(1:Pn) * Basis(1:Pn) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                     DO l=1,DIM
! ! ! ! ! ! !                        Stress(j,k,i) = Stress(j,k,i) - &
! ! ! ! ! ! !                            Metric(j,k) * (2.0d0/3.0d0) * Viscosity * Grad(l,l)
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !         END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !         EdgeLength = EdgeLength + s
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Jump = MATMUL( ( Stress(:,:,1) - Stress(:,:,2)), Normal )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !            ResidualNorm = ResidualNorm + s * SUM( Jump(1:DIM) ** 2 )
! ! ! ! ! ! !         ELSE
! ! ! ! ! ! !            CALL InvertMatrix( Metric,3 )
! ! ! ! ! ! !            DO i=1,DIM
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  ResidualNorm = ResidualNorm + s*Metric(i,j)*Jump(i)*Jump(j)
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END DO
! ! ! ! ! ! !         END IF
! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Indicator = EdgeLength * ResidualNorm
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z)
! ! ! ! ! ! !      DEALLOCATE( EdgeNodes % x, EdgeNodes % y, EdgeNodes % z)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DEALLOCATE( NodalViscosity, x, y, z, EdgeBasis, &
! ! ! ! ! ! !            Basis, dBasisdx, Velocity, Pressure )
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !   END FUNCTION FlowEdgeResidual
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! !> Compute the residual of the Navier-Stokes equation for the bulk elements.
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !    FUNCTION FlowInsideResidual( Model, Element,  &
! ! ! ! ! ! !           Mesh, Quant, Perm, Fnorm ) RESULT( Indicator )
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !      USE DefUtils
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !      IMPLICIT NONE
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !      TYPE(Model_t) :: Model
! ! ! ! ! ! !      INTEGER :: Perm(:)
! ! ! ! ! ! !      REAL(KIND=dp) :: Quant(:), Indicator(2), FNorm
! ! ! ! ! ! !      TYPE( Mesh_t ), POINTER    :: Mesh
! ! ! ! ! ! !      TYPE( Element_t ), POINTER :: Element
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(Nodes_t) :: Nodes
! ! ! ! ! ! ! 
! ! ! ! ! ! !      INTEGER :: i,j,k,l,m,n,t,DIM,DOFs
! ! ! ! ! ! ! 
! ! ! ! ! ! !      LOGICAL :: stat, GotIt, Compressible, Convect
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE( Variable_t ), POINTER :: Var
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp) :: SqrtMetric, Metric(3,3), Symb(3,3,3), dSymb(3,3,3,3)
! ! ! ! ! ! !      REAL(KIND=dp) :: Density, Viscosity,u, v, w, s, detJ
! ! ! ! ! ! !      REAL(KIND=dp) :: Residual(4), ResidualNorm, Area, ReferencePressure, dt
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp), ALLOCATABLE :: NodalViscosity(:), NodalDensity(:), &
! ! ! ! ! ! !             Basis(:),  dBasisdx(:,:), ddBasisddx(:,:,:)
! ! ! ! ! ! !      REAL(KIND=dp),ALLOCATABLE :: Velocity(:,:), Pressure(:)
! ! ! ! ! ! !      REAL(KIND=dp),ALLOCATABLE :: PrevVelo(:,:), PrevPres(:)
! ! ! ! ! ! !      REAL(KIND=dp),ALLOCATABLE :: Temperature(:), NodalForce(:,:)
! ! ! ! ! ! !      REAL(KIND=dp),ALLOCATABLE :: HeatCapacity(:), ReferenceTemperature(:), &
! ! ! ! ! ! !                       HeatExpansionCoeff(:)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp) :: SpecificHeatRatio
! ! ! ! ! ! ! 
! ! ! ! ! ! !      REAL(KIND=dp), POINTER :: Gravity(:,:)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(ValueList_t), POINTER :: Material
! ! ! ! ! ! ! 
! ! ! ! ! ! !      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !    Initialize:
! ! ! ! ! ! ! !    -----------
! ! ! ! ! ! !      Indicator = 0.0d0
! ! ! ! ! ! !      FNorm = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IF ( ANY( Perm( Element % NodeIndexes ) <= 0 ) ) RETURN
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Metric = 0.0d0
! ! ! ! ! ! !      DO i=1,3
! ! ! ! ! ! !         Metric(i,i) = 1.0d0
! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !      SELECT CASE( CurrentCoordinateSystem() )
! ! ! ! ! ! !         CASE( AxisSymmetric, CylindricSymmetric )
! ! ! ! ! ! !            DIM = 3
! ! ! ! ! ! !         CASE DEFAULT
! ! ! ! ! ! !            DIM = CoordinateSystemDimension()
! ! ! ! ! ! !      END SELECT
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DOFs = DIM + 1
! ! ! ! ! ! !      IF ( CurrentCoordinateSystem() == AxisSymmetric ) DOFs = DOFs-1
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Element nodal points:
! ! ! ! ! ! ! !    ---------------------
! ! ! ! ! ! !      n = Element % TYPE % NumberOfNodes
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
! ! ! ! ! ! !      Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
! ! ! ! ! ! !      Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
! ! ! ! ! ! !      Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      ALLOCATE( NodalViscosity(n), NodalDensity(n), Basis(n), dBasisdx(n,3), &
! ! ! ! ! ! !         ddBasisddx(n,3,3), Velocity(3,n), Pressure(n), PrevVelo(3,n),  &
! ! ! ! ! ! !         PrevPres(n), Temperature(n), NodalForce(4,n), HeatCapacity(n), &
! ! ! ! ! ! !         ReferenceTemperature(n), HeatExpansionCoeff(n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Material parameters: density, viscosity, etc.
! ! ! ! ! ! ! !    ----------------------------------------------
! ! ! ! ! ! !      k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Material', &
! ! ! ! ! ! !                   minv=1, maxv=Model % NumberOfMaterials )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Material => Model % Materials(k) % Values
! ! ! ! ! ! ! 
! ! ! ! ! ! !      NodalDensity(1:n) = ListGetReal( &
! ! ! ! ! ! !          Material, 'Density', n, Element % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      NodalViscosity(1:n) = ListGetReal( &
! ! ! ! ! ! !          Material, 'Viscosity', n, Element % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      k = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Equation', &
! ! ! ! ! ! !                       minv=1, maxv=Model % NumberOfEquations   )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      Convect = ListGetLogical( Model % Equations(k) % Values, &
! ! ! ! ! ! !                    'NS Convect', GotIt )
! ! ! ! ! ! !      IF ( .NOT. GotIt ) Convect = .TRUE.
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Elementwise nodal solution:
! ! ! ! ! ! ! !    ---------------------------
! ! ! ! ! ! !      Velocity = 0.0d0
! ! ! ! ! ! !      DO k=1,DOFs-1
! ! ! ! ! ! !         Velocity(k,1:n) = Quant( DOFs*Perm(Element % NodeIndexes)-DOFs+k )
! ! ! ! ! ! !      END DO
! ! ! ! ! ! !      Pressure(1:n) = Quant( DOFs*Perm(Element % NodeIndexes) )
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Check for time dep.
! ! ! ! ! ! ! !    -------------------
! ! ! ! ! ! !      PrevPres(1:n)     = Pressure(1:n)
! ! ! ! ! ! !      PrevVelo(1:3,1:n) = Velocity(1:3,1:n)
! ! ! ! ! ! ! 
! ! ! ! ! ! !      dt = Model % Solver % dt
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IF ( ListGetString( Model % Simulation, 'Simulation Type') == 'transient' ) THEN
! ! ! ! ! ! !         Var => VariableGet( Model % Variables, 'Flow Solution', .TRUE. )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         PrevVelo = 0.0d0
! ! ! ! ! ! !         DO k=1,DOFs-1
! ! ! ! ! ! !            PrevVelo(k,1:n) = &
! ! ! ! ! ! !               Var % PrevValues(DOFs*Var % Perm(Element % NodeIndexes)-DOFs+k,1)
! ! ! ! ! ! !         END DO
! ! ! ! ! ! !         PrevPres(1:n)=Var % PrevValues(DOFs*Var % Perm(Element % NodeIndexes),1)
! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Check for compressible flow equations:
! ! ! ! ! ! ! !    --------------------------------------
! ! ! ! ! ! !      Compressible = .FALSE.
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IF (  ListGetString( Material, 'Compressibility Model', GotIt ) == &
! ! ! ! ! ! !                       'perfect gas equation 1' ) THEN
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Compressible = .TRUE.
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Var => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )
! ! ! ! ! ! !         IF ( ASSOCIATED( Var ) ) THEN
! ! ! ! ! ! !            Temperature(1:n) = &
! ! ! ! ! ! !                Var % Values( Var % Perm(Element % NodeIndexes) )
! ! ! ! ! ! !         ELSE
! ! ! ! ! ! !            Temperature(1:n) = ListGetReal( Material, &
! ! ! ! ! ! !                'Reference Temperature',n,Element % NodeIndexes )
! ! ! ! ! ! !         END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !         SpecificHeatRatio = ListGetConstReal( Material, &
! ! ! ! ! ! !                   'Specific Heat Ratio' )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         ReferencePressure = ListGetConstReal( Material, &
! ! ! ! ! ! !                    'Reference Pressure' )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         HeatCapacity(1:n) = ListGetReal( Material, &
! ! ! ! ! ! !                       'Heat Capacity',n,Element % NodeIndexes )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         NodalDensity(1:n) =  (Pressure(1:n) + ReferencePressure) * SpecificHeatRatio / &
! ! ! ! ! ! !               ( (SpecificHeatRatio - 1) * HeatCapacity(1:n) * Temperature(1:n) )
! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Body Forces:
! ! ! ! ! ! ! !    ------------
! ! ! ! ! ! ! !
! ! ! ! ! ! !      k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, &
! ! ! ! ! ! !        'Body Force', GotIt, 1, Model % NumberOfBodyForces )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      NodalForce = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IF ( GotIt .AND. k > 0  ) THEN
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Boussinesq approximation of heat expansion for
! ! ! ! ! ! ! !       incompressible flow equations:
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Density for the force term equals to
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       \rho = rho_0 (1-\beta(T-T_0)),
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       where \beta is the  heat expansion  coefficient,
! ! ! ! ! ! ! !       T temperature and \rho_0 and T_0 correspond to
! ! ! ! ! ! ! !       stress free state. Otherwise density is assumed
! ! ! ! ! ! ! !       constant.
! ! ! ! ! ! ! !       ----------------------------------------------
! ! ! ! ! ! !         IF (ListGetLogical(Model % BodyForces(k) % Values,'Boussinesq',GotIt)) THEN
! ! ! ! ! ! ! 
! ! ! ! ! ! !            Var => VariableGet( Mesh % Variables, 'Temperature', .TRUE. )
! ! ! ! ! ! !            IF ( ASSOCIATED( Var ) ) THEN
! ! ! ! ! ! !               Temperature(1:n) = &
! ! ! ! ! ! !                   Var % Values( Var % Perm(Element % NodeIndexes) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               HeatExpansionCoeff(1:n) = ListGetReal( Material, &
! ! ! ! ! ! !                  'Heat Expansion Coefficient',n,Element % NodeIndexes )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               ReferenceTemperature(1:n) = ListGetReal( Material, &
! ! ! ! ! ! !                  'Reference Temperature',n,Element % NodeIndexes )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               Gravity => ListGetConstRealArray( Model % Constants, &
! ! ! ! ! ! !                              'Gravity' )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               k = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Equation', &
! ! ! ! ! ! !                         minv=1, maxv=Model % NumberOfEquations )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               IF ( ListGetLogical( Model % Equations(k) % Values, &
! ! ! ! ! ! !                             'Hydrostatic Pressure', GotIt) ) THEN
! ! ! ! ! ! !                  DO i=1,DIM
! ! ! ! ! ! !                     NodalForce(i,1:n) = ( 1 - HeatExpansionCoeff(1:n) * &
! ! ! ! ! ! !                        ( Temperature(1:n) - ReferenceTemperature(1:n) ) ) * &
! ! ! ! ! ! !                             Gravity(i,1) * Gravity(4,1)
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               ELSE
! ! ! ! ! ! !                  DO i=1,DIM
! ! ! ! ! ! !                     NodalForce(i,1:n) = ( -HeatExpansionCoeff(1:n) * &
! ! ! ! ! ! !                        ( Temperature(1:n) - ReferenceTemperature(1:n) ) ) * &
! ! ! ! ! ! !                             Gravity(i,1) * Gravity(4,1)
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END IF
! ! ! ! ! ! !            END IF
! ! ! ! ! ! !         END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Given external force:
! ! ! ! ! ! ! !       ---------------------
! ! ! ! ! ! !         NodalForce(1,1:n) = NodalForce(1,1:n) + ListGetReal( &
! ! ! ! ! ! !              Model % BodyForces(k) % Values, 'Flow BodyForce 1', &
! ! ! ! ! ! !                   n, Element % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         NodalForce(2,1:n) = NodalForce(2,1:n) + ListGetReal( &
! ! ! ! ! ! !              Model % BodyForces(k) % Values, 'Flow BodyForce 2', &
! ! ! ! ! ! !                   n, Element % NodeIndexes, GotIt )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         NodalForce(3,1:n) = NodalForce(3,1:n) + ListGetReal( &
! ! ! ! ! ! !              Model % BodyForces(k) % Values, 'Flow BodyForce 3', &
! ! ! ! ! ! !                   n, Element % NodeIndexes, GotIt )
! ! ! ! ! ! !      END IF
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !    Integrate square of residual over element:
! ! ! ! ! ! ! !    ------------------------------------------
! ! ! ! ! ! !      ResidualNorm = 0.0d0
! ! ! ! ! ! !      Area = 0.0d0
! ! ! ! ! ! ! 
! ! ! ! ! ! !      IntegStuff = GaussPoints( Element )
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DO t=1,IntegStuff % n
! ! ! ! ! ! !         u = IntegStuff % u(t)
! ! ! ! ! ! !         v = IntegStuff % v(t)
! ! ! ! ! ! !         w = IntegStuff % w(t)
! ! ! ! ! ! ! 
! ! ! ! ! ! !         stat = ElementInfo( Element, Nodes, u, v, w, detJ, &
! ! ! ! ! ! !             Basis, dBasisdx, ddBasisddx, .TRUE. )
! ! ! ! ! ! ! 
! ! ! ! ! ! !         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !            s = IntegStuff % s(t) * detJ
! ! ! ! ! ! !         ELSE
! ! ! ! ! ! !            u = SUM( Basis(1:n) * Nodes % x(1:n) )
! ! ! ! ! ! !            v = SUM( Basis(1:n) * Nodes % y(1:n) )
! ! ! ! ! ! !            w = SUM( Basis(1:n) * Nodes % z(1:n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !            CALL CoordinateSystemInfo( Metric, SqrtMetric, &
! ! ! ! ! ! !                       Symb, dSymb, u, v, w )
! ! ! ! ! ! !            s = IntegStuff % s(t) * detJ * SqrtMetric
! ! ! ! ! ! !         END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !         Density   = SUM( NodalDensity(1:n)   * Basis(1:n) )
! ! ! ! ! ! !         Viscosity = SUM( NodalViscosity(1:n) * Basis(1:n) )
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Residual of the navier-stokes equations:
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       or more generally:
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       ----------------------------------------------------------
! ! ! ! ! ! ! !
! ! ! ! ! ! !         Residual = 0.0d0
! ! ! ! ! ! !         DO i=1,DIM
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          given force:
! ! ! ! ! ! ! !          -------------
! ! ! ! ! ! !            Residual(i) = -Density * SUM( NodalForce(i,1:n) * Basis(1:n) )
! ! ! ! ! ! !  
! ! ! ! ! ! !            IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! ! !             + grad(p):
! ! ! ! ! ! ! !             ----------
! ! ! ! ! ! !               Residual(i) = Residual(i) + SUM( Pressure(1:n) * dBasisdx(1:n,i) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !                - 2 ( \mu \epsilon^{ij} )_{,j}:
! ! ! ! ! ! ! !                -------------------------------
! ! ! ! ! ! !                  Residual(i) = Residual(i) - Viscosity * &
! ! ! ! ! ! !                      SUM( Velocity(i,1:n) * ddBasisddx(1:n,j,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                  Residual(i) = Residual(i) - &
! ! ! ! ! ! !                       SUM( NodalViscosity(1:n) * dBasisdx(1:n,j) ) * &
! ! ! ! ! ! !                           SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                   Residual(i) = Residual(i) - Viscosity * &
! ! ! ! ! ! !                       SUM( Velocity(j,1:n) * ddBasisddx(1:n,i,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                   Residual(i) = Residual(i) - &
! ! ! ! ! ! !                       SUM( NodalViscosity(1:n) * dBasisdx(1:n,j) ) * &
! ! ! ! ! ! !                           SUM( Velocity(j,1:n) * dBasisdx(1:n,i) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                   IF ( Compressible ) THEN
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !                    + (2/3) grad(\mu div(u)):
! ! ! ! ! ! ! !                    -------------------------
! ! ! ! ! ! !                      Residual(i) = Residual(i) + &
! ! ! ! ! ! !                         Viscosity * ( 2.0d0 / 3.0d0 ) * &
! ! ! ! ! ! !                            SUM( Velocity(j,1:n) * ddBasisddx(1:n,j,i) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                      Residual(i) = Residual(i) + &
! ! ! ! ! ! !                          SUM( NodalViscosity(1:n) * dBasisdx(1:n,i) ) * &
! ! ! ! ! ! !                              SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                   END IF
! ! ! ! ! ! !               END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !               IF ( Convect ) THEN
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !                + \rho * (@u/@t + u.grad(u)):
! ! ! ! ! ! ! !                -----------------------------
! ! ! ! ! ! !                  Residual(i) = Residual(i) + Density *  &
! ! ! ! ! ! !                      SUM((Velocity(i,1:n)-PrevVelo(i,1:n))*Basis(1:n)) / dt
! ! ! ! ! ! ! 
! ! ! ! ! ! !                  DO j=1,DIM
! ! ! ! ! ! !                     Residual(i) = Residual(i) + &
! ! ! ! ! ! !                         Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
! ! ! ! ! ! !                             SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END IF
! ! ! ! ! ! !            ELSE
! ! ! ! ! ! ! !             + g^{ij}p_{,j}:
! ! ! ! ! ! ! !             ---------------
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  Residual(i) = Residual(i) + Metric(i,j) * &
! ! ! ! ! ! !                       SUM( Pressure(1:n) * dBasisdx(1:n,i) )
! ! ! ! ! ! !               END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !             - g^{jk} (\mu u^i_{,k})_{,j}):
! ! ! ! ! ! ! !             ------------------------------
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     Residual(i) = Residual(i) -   &
! ! ! ! ! ! !                          Metric(j,k) * Viscosity * &
! ! ! ! ! ! !                          SUM( Velocity(i,1:n) * ddBasisddx(1:n,j,k) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                     DO l=1,DIM
! ! ! ! ! ! !                        Residual(i) = Residual(i) +  &
! ! ! ! ! ! !                             Metric(j,k) * Viscosity * Symb(j,k,l) * &
! ! ! ! ! ! !                             SUM( Velocity(i,1:n) * dBasisdx(1:n,l) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        Residual(i) = Residual(i) -  &
! ! ! ! ! ! !                             Metric(j,k) * Viscosity * Symb(l,j,i) * &
! ! ! ! ! ! !                             SUM( Velocity(l,1:n) * dBasisdx(1:n,k) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        Residual(i) = Residual(i) -  &
! ! ! ! ! ! !                             Metric(j,k) * Viscosity * Symb(l,k,i) * &
! ! ! ! ! ! !                             SUM( Velocity(l,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        Residual(i) = Residual(i) -  &
! ! ! ! ! ! !                             Metric(j,k) * Viscosity * dSymb(l,j,i,k) * &
! ! ! ! ! ! !                             SUM( Velocity(l,1:n) * Basis(1:n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        DO m=1,DIM
! ! ! ! ! ! !                           Residual(i) = Residual(i) - Metric(j,k) * Viscosity *&
! ! ! ! ! ! !                                   Symb(m,k,i) * Symb(l,j,m) * &
! ! ! ! ! ! !                                         SUM( Velocity(l,1:n) * Basis(1:n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                           Residual(i) = Residual(i) + Metric(j,k) * Viscosity *&
! ! ! ! ! ! !                                   Symb(j,k,m) * Symb(l,m,i) * &
! ! ! ! ! ! !                                         SUM( Velocity(l,1:n) * Basis(1:n) )
! ! ! ! ! ! !                        END DO
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !             - g^{ik} (\mu u^j_{,k})_{,j}):
! ! ! ! ! ! ! !             ------------------------------
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  DO k=1,DIM
! ! ! ! ! ! !                     Residual(i) = Residual(i) -   &
! ! ! ! ! ! !                          Metric(i,k) * Viscosity * &
! ! ! ! ! ! !                          SUM( Velocity(j,1:n) * ddBasisddx(1:n,j,k) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                     DO l=1,DIM
! ! ! ! ! ! !                        Residual(i) = Residual(i) +  &
! ! ! ! ! ! !                             Metric(i,k) * Viscosity * Symb(j,k,l) * &
! ! ! ! ! ! !                             SUM( Velocity(j,1:n) * dBasisdx(1:n,l) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        Residual(i) = Residual(i) -  &
! ! ! ! ! ! !                             Metric(i,k) * Viscosity * Symb(l,j,j) * &
! ! ! ! ! ! !                             SUM( Velocity(l,1:n) * dBasisdx(1:n,k) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        Residual(i) = Residual(i) -  &
! ! ! ! ! ! !                             Metric(i,k) * Viscosity * Symb(l,k,j) * &
! ! ! ! ! ! !                             SUM( Velocity(l,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        Residual(i) = Residual(i) -  &
! ! ! ! ! ! !                             Metric(i,k) * Viscosity * dSymb(l,j,j,k) * &
! ! ! ! ! ! !                             SUM( Velocity(l,1:n) * Basis(1:n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                        DO m=1,DIM
! ! ! ! ! ! !                           Residual(i) = Residual(i) - Metric(i,k) * Viscosity *&
! ! ! ! ! ! !                                   Symb(m,k,j) * Symb(l,j,m) * &
! ! ! ! ! ! !                                         SUM( Velocity(l,1:n) * Basis(1:n) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                           Residual(i) = Residual(i) + Metric(i,k) * Viscosity *&
! ! ! ! ! ! !                                   Symb(j,k,m) * Symb(l,m,j) * &
! ! ! ! ! ! !                                         SUM( Velocity(l,1:n) * Basis(1:n) )
! ! ! ! ! ! !                        END DO
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !               IF ( Convect ) THEN
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !                + \rho * (@u/@t + u^j u^i_{,j}):
! ! ! ! ! ! ! !                --------------------------------
! ! ! ! ! ! !                  Residual(i) = Residual(i) + Density *  &
! ! ! ! ! ! !                      SUM((Velocity(i,1:n)-PrevVelo(i,1:n))*Basis(1:n)) / dt
! ! ! ! ! ! ! 
! ! ! ! ! ! !                  DO j=1,DIM
! ! ! ! ! ! !                     Residual(i) = Residual(i) + &
! ! ! ! ! ! !                          Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
! ! ! ! ! ! !                          SUM( Velocity(i,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !                     DO k=1,DIM
! ! ! ! ! ! !                        Residual(i) = Residual(i) + &
! ! ! ! ! ! !                             Density * SUM( Velocity(j,1:n) * Basis(1:n) ) * &
! ! ! ! ! ! !                             Symb(j,k,i) * SUM( Velocity(k,1:n) * Basis(1:n) )
! ! ! ! ! ! !                     END DO
! ! ! ! ! ! !                  END DO
! ! ! ! ! ! !               END IF
! ! ! ! ! ! !            END IF
! ! ! ! ! ! !         END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !       Continuity equation:
! ! ! ! ! ! ! !       --------------------
! ! ! ! ! ! !         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          + \rho * div(u):
! ! ! ! ! ! ! !          ----------------
! ! ! ! ! ! !            DO j=1,DIM
! ! ! ! ! ! !               Residual(DIM+1) = Residual(DIM+1) + &
! ! ! ! ! ! !                    Density * SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! !            END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( Compressible ) THEN
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !             + u.grad(\rho):
! ! ! ! ! ! ! !             ----------------
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  Residual(DIM+1) = Residual(DIM+1) + &
! ! ! ! ! ! !                       SUM( Velocity(j,1:n) * Basis(1:n) ) *  &
! ! ! ! ! ! !                            SUM( NodalDensity(1:n) * dBasisdx(1:n,j) ) 
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END IF
! ! ! ! ! ! !         ELSE
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !          + \rho * u^j_{,j}:
! ! ! ! ! ! ! !          ------------------
! ! ! ! ! ! !            DO j=1,DIM
! ! ! ! ! ! !               Residual(DIM+1) = Residual(DIM+1) + &
! ! ! ! ! ! !                    Density * SUM( Velocity(j,1:n) * dBasisdx(1:n,j) )
! ! ! ! ! ! ! 
! ! ! ! ! ! !               DO k=1,DIM
! ! ! ! ! ! !                  Residual(DIM+1) = Residual(DIM+1) + Density * &
! ! ! ! ! ! !                       Symb(k,j,j) * SUM( Velocity(k,1:n) * Basis(1:n) )
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! !            IF ( Compressible ) THEN
! ! ! ! ! ! ! !
! ! ! ! ! ! ! !             + u^j \rho_{,j}:
! ! ! ! ! ! ! !             ----------------
! ! ! ! ! ! !               DO j=1,DIM
! ! ! ! ! ! !                  Residual(DIM+1) = Residual(DIM+1) + &
! ! ! ! ! ! !                       SUM( Velocity(j,1:n) * Basis(1:n) ) *  &
! ! ! ! ! ! !                       SUM( NodalDensity(1:n) * dBasisdx(1:n,j) ) 
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END IF
! ! ! ! ! ! !         END IF
! ! ! ! ! ! ! 
! ! ! ! ! ! !         DO i=1,DIM
! ! ! ! ! ! !            FNorm = FNorm + s * (Density * SUM(NodalForce(i,1:n)*Basis(1:n))**2)
! ! ! ! ! ! !         END DO 
! ! ! ! ! ! !         Area = Area + s
! ! ! ! ! ! ! 
! ! ! ! ! ! !         IF ( CurrentCoordinateSystem() == Cartesian ) THEN
! ! ! ! ! ! !            ResidualNorm = ResidualNorm + &
! ! ! ! ! ! !              s * (Element % hK**2 * SUM(Residual(1:dim)**2) + Residual(dim+1)**2 )
! ! ! ! ! ! !         ELSE
! ! ! ! ! ! !            CALL InvertMatrix( Metric,3 )
! ! ! ! ! ! !            DO i=1,dim
! ! ! ! ! ! !               DO j=1,dim
! ! ! ! ! ! !                  ResidualNorm = ResidualNorm + &
! ! ! ! ! ! !                     s * Element % hK **2 * Metric(i,j) * Residual(i) * Residual(j)
! ! ! ! ! ! !               END DO
! ! ! ! ! ! !            END DO
! ! ! ! ! ! !            ResidualNorm = ResidualNorm + s * Residual(dim+1)**2
! ! ! ! ! ! !         END IF
! ! ! ! ! ! !      END DO
! ! ! ! ! ! ! 
! ! ! ! ! ! ! !    FNorm = Area * FNorm
! ! ! ! ! ! !      Indicator = ResidualNorm
! ! ! ! ! ! ! 
! ! ! ! ! ! !      DEALLOCATE( NodalViscosity, NodalDensity, Basis, dBasisdx,           &
! ! ! ! ! ! !         ddBasisddx, Velocity, Pressure, PrevVelo, PrevPres, Temperature,   &
! ! ! ! ! ! !         NodalForce, HeatCapacity, ReferenceTemperature, HeatExpansionCoeff,&
! ! ! ! ! ! !         Nodes % x, Nodes % y, Nodes % z )
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! !   END FUNCTION FlowInsideResidual
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! ! ! 

!  ______ _   _ _____  
! |  ____| \ | |  __ \ 
! | |__  |  \| | |  | |
! |  __| | . ` | |  | |
! | |____| |\  | |__| |
! |______|_| \_|_____/    
!  _   _ _         ____                       _             
! | \ | | |       / __ \                     | |            
! |  \| | |      | |  | |_ __   ___ _ __ __ _| |_ ___  _ __ 
! | . ` | |      | |  | | '_ \ / _ \ '__/ _` | __/ _ \| '__|
! | |\  | |____  | |__| | |_) |  __/ | | (_| | || (_) | |   
! |_| \_|______|  \____/| .__/ \___|_|  \__,_|\__\___/|_|   
!                       | |                                 
!                       |_|                                 
!    
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

















! DEBUG TEST
  SUBROUTINE  ANMPFrhsFQDEBUGLOOP( FQMan, UMan, IO, NSDOFs, FlowPerm, &
                         USAV, GradSAV, &
                         FQManTemp, Uelex , Ueley , Uelez, &
                         Density,Material,FlowSolution_Init, &
                         DEBTU,DEBTGU,DEBTSSS,DEBASS,DEBINIel)   
!  

    USE DefUtils
    USE Differentials
    USE MaterialModels
    USE Adaptive    
    USE SolverUtils
!

    IMPLICIT NONE     
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: FQMan(:),FlowSolution_Init(:)
    REAL(KIND=dp) :: UMan(:,:)    
    REAL(KIND=dp) :: USAV(:,:,:,:)
    REAL(KIND=dp) :: GradSAV( :,:,:,:,: )   
    INTEGER       :: NSDOFs, IO
    REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:)     
    TYPE(ValueList_t),POINTER :: Material    
    REAL(KIND=dp) :: Density(:)
    REAL(KIND=dp) :: FQManTemp(:)
    
    INTEGER, POINTER :: FlowPerm(:)
     
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
    INTEGER :: t,n,nb,nd
    REAL(KIND=dp), POINTER :: FQlocal(:)   
    
    
    INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) ::  TIC, TOTupdate,TOTFQelem,TOTinit
    REAL(KIND=dp) ::  DEBTU,DEBTGU,DEBTSSS,DEBASS,DEBINIel
!------------------------------------------------------------------------------
    TOTinit = 0.0_dp
    TOTupdate = 0.0_dp
    TOTFQelem = 0.0_dp    
    FQManTemp = 0.0_dp
    

    
    DO t = 1 , GetNOFActive()  ! Boucle ELEMENT
      TIC = CPUTime()    
      CALL AdvanceOutput( t,GetNOFActive() )      
      Element => GetActiveElement(t)      
!       NodeIndexes => Element % NodeIndexes
      Indexes => Element % NodeIndexes
      
      n = GetElementNOFNodes(Element)      
      nb = GetElementNOFBDOFs(Element)       
      nd = GetElementDOFs( Indexes )
      CALL GetElementNodes( ElementNodes )
!
      Density(1:n) = GetReal( Material, 'Density' )

!------------------------------------------------------------------------------
!     Vitesse au point d'intégration de l'itération précédente
!------------------------------------------------------------------------------
      SELECT CASE( NSDOFs )
        CASE(3) ! 2D
          Uelex(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
          Ueley(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
          Uelez(1:n) = 0.0_dp
        CASE(4) ! 3D
          Uelex(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-3 , IO-1)
          Ueley(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-2 , IO-1)
          Uelez(1:n) = UMan( NSDOFs*FlowPerm(Indexes(1:nd))-1 , IO-1)
      END SELECT
         TIC = CPUTime() - TIC
         TOTinit = TOTinit + TIC
!          WRITE(6,*) 'Element ', t, '- ANMPFrhsFQDEBUGLOOP init',TIC      
!------------------------------------------------------------------------------  
         TIC = CPUTime()
      CALL MANFQelemOPTIDEBUGLOOP( FQManTemp, FQlocal, Density,            &
                          Element, n, IO, ElementNodes,                    &
                          USAV, GradSAV, Element % ElementIndex , NSDOFs,  &
                          Uelex, Ueley, Uelez , &
                         DEBTU,DEBTGU,DEBTSSS,DEBASS,DEBINIel)      
!------------------------------------------------------------------------------                                       
         TIC = CPUTime() - TIC
!          WRITE(6,*) 'Element ', t, '- MANFQelemOPTIDEBUGLOOP',TIC
         TOTFQelem = TOTFQelem + TIC
         
         TIC = CPUTime()         
      CALL UpdateGlobalForce( FQMan, FQManTemp, n, NSDOFs, FlowPerm(Indexes(1:n)), UElement=Element )
         TIC = CPUTime() - TIC
!          WRITE(6,*) 'Element ', t, '- UpdateGlobalForce',TIC      
         TOTupdate = TOTupdate + TIC
    END DO ! FIN Boucle ELEMENT
    
    WRITE(*,*) " TOTAL UpdateGlobalForce= ",TOTupdate
    WRITE(*,*) " TOTAL TOTinit= ",TOTinit
    WRITE(*,*) " TOTAL MANFQelemOPTIDEBUGLOOP= ",TOTFQelem
!------------------------------------------------------------------------------
!       Conditions de Diriclet modifiées pour le 2nd membre----> 0           
    CALL CalCondLimCoeff( FQMan, 0.0_dp, 1,FlowSolution_Init )
!                           
 END SUBROUTINE ANMPFrhsFQDEBUGLOOP



 SUBROUTINE  MANFQelemOPTIDEBUGLOOP( FQManTemp, FQlocal, Nodalrho,     &
                            Element, n, io, ElementNodes,     &
                            USAV , GradSAV, numel , NSDOFs,   &
                            Uelex, Ueley, Uelez , &
                            DEBTU,DEBTGU,DEBTSSS,DEBASS,DEBINIel) 
!
  USE DefUtils
  USE Differentials
  USE MaterialModels
  USE Adaptive
  USE SolverUtils

        IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
     REAL(KIND=dp), TARGET :: FQManTemp(:) 
     REAL(KIND=dp), POINTER :: FQlocal(:)
     REAL(KIND=dp) :: Nodalrho(:)
  
     ! Composante de la vitesse MAN ordre r et p-r
     INTEGER :: io, n   ! p = ordremanencours
     ! OMan ordre de la série MAN, nmax le nombre max de degré de liberté
     INTEGER, POINTER :: Flowperm(:)
      
     INTEGER, POINTER :: NodeIndexes(:), Indexes(:)
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t), POINTER :: Element
!----     
     REAL(KIND=dp) :: USAV(:,:,:,:)
     REAL(KIND=dp) :: GradSAV( :,:,:,:,: )
     INTEGER       :: numel,NSDOFs
     REAL(KIND=dp) ,ALLOCATABLE, TARGET :: UMan(:,:)     
     REAL(KIND=dp) :: Uelex(:), Ueley(:), Uelez(:)     
!----   

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(2*n),dBasisdx(2*n,3),ddBasisddx(n,3,3)
     REAL(KIND=dp) :: detJ, NodalBasis(n), dLBasisdx(n,3)
     REAL(KIND=dp) :: PBasis(n), pdBasisdx(n,3),BaseP, dNodalBasisdx(n,n,3)


     REAL(KIND=dp) :: Ur(3), GradUpmr(3,3), FQelemtemp(4)

     INTEGER :: i,j,k,l,p,c,q,t,r,dim

     REAL(KIND=dp) :: s,u,v,w,volume
     REAL(KIND=dp) :: rho
    
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ, NBasis

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
     TYPE(ValueList_t), POINTER :: Material, Params

     LOGICAL :: Found, stat
!------------------------------------------------------------------------------
     
     REAL(KIND=dp) :: tmp
     ! DEBUG LOOP 
     REAL(KIND=dp) :: UTMP, GTMP, TIC
    REAL(KIND=dp) ::  DEBTU,DEBTGU,DEBTSSS,DEBASS,DEBINIel
     INTEGER :: itmp
    
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!              write(*,*) 'io=',io
     dim = CoordinateSystemDimension()
     Params => GetSolverParams()

     c = dim + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!     
     FQManTemp = 0.0_dp 
!!!!!!!!!!!!!!!!!!!!!!!!!!     
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     NBasis = n

     IntegStuff = GaussPoints( Element )

     U_Integ => IntegStuff % u   
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------    
!              write(*,*) 'io=',io
!******************************************************************************

    DO t=1,N_Integ                        ! boucle sur les points d'intégration
      TIC = CPUTime()
    
      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element, ElementNodes, u, v, w, detJ, &
              Basis, dBasisdx )
!------------------------------------------------------------------------------
      s = detJ * S_Integ(t)

      rho  = SUM( Nodalrho(1:n)*Basis(1:n) )
      
      TIC = CPUTime() - TIC      
      DEBINIel = DEBINIel + TIC      
! 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! On sauvegarde U(element,gausspoint)  et Grad U(element,gausspoint) à l'ordre en cours
! en procédant ainsi on ne recalcule pas à chaques ordre, des termes déjà calculé aux ordres précédents.
!------------------------------------------------------------------------------
      TIC = CPUTime()

      USAV( numel , t , io-1 , 1 ) = SUM( Basis(1:n) * Uelex(1:n) )
      USAV( numel , t , io-1 , 2 ) = SUM( Basis(1:n) * Ueley(1:n) )
      IF (dim > 2) USAV( numel , t , io-1 , 3) = SUM( Basis(1:n) * Uelez(1:n) )
      TIC = CPUTime() - TIC      
      DEBTU = DEBTU + TIC
!       WRITE(6,*) 'DEBUG TIME USAV ',TIC
      
!------------------------------------------------------------------------------
! Puis le gradient
! Probleme : comment stocker : element, gausspoint, ordre, composante1, composante2
      TIC = CPUTime()

      DO j=1,3
         GradSAV( numel , t , io-1 , 1 , j ) = SUM( Uelex(1:n) * dBasisdx(1:n,j) )
         GradSAV( numel , t , io-1 , 2 , j ) = SUM( Ueley(1:n) * dBasisdx(1:n,j) )
         IF ( DIM > 2 ) GradSAV( numel , t , io-1 , 3 , j ) = SUM( Uelez(1:n) * dBasisdx(1:n,j) )
      END DO
      TIC = CPUTime() - TIC      
      DEBTGU = DEBTGU + TIC

!       WRITE(6,*) 'DEBUG TIME GradSAV- ',TIC
      
!------------------------------------------------------------------------------
! Somme des termes croisées par composante spatiale
!------------------------------------------------------------------------------
      FQelemtemp = 0.0_dp
      
!       TIC = CPUTime()
!       
!       DO r=1, io-1      ! Boucle CROSS
!        DO i = 1 ,dim
!         DO j = 1 ,dim
!           UTMP = USAV( numel , t , r , j )
!           GTMP = GradSAV( numel , t , io-r , i , j )
!           FQelemtemp(i) = FQelemtemp(i) -  UTMP * GTMP
!         END DO
!        END DO
!       END DO ! FIN Boucle CROSS
!       TIC = CPUTime() - TIC      
!       WRITE(6,*) 'DEBUG TIME MANFQelemOPTIDEBUGLOOP DO r DO i DO j : ',TIC

      TIC = CPUTime()
      DO j = 1 ,dim
       DO i = 1 ,dim
        DO r=1, io-1      ! Boucle CROSS
         
          itmp = io-r
          UTMP =    USAV( numel , t , r , j )
          GTMP = GradSAV( numel , t , itmp , i , j )
          FQelemtemp(i) = FQelemtemp(i) -  UTMP * GTMP
          
        END DO
       END DO
      END DO ! FIN Boucle CROSS
      TIC = CPUTime() - TIC 
      DEBTSSS = DEBTSSS + TIC
      
!       WRITE(6,*) 'DEBUG TIME MANFQelemOPTIDEBUGLOOP DO j DO i DO r : ',TIC
      
!------------------------------------------------------------------------------
! Assemblage
      TIC = CPUTime()

      DO p = 1 , NBasis !boucle sur les fonctions de formes des fonctions tests Wi = Ni
        FQlocal => FQManTemp( c*(p-1)+1 : c*(p-1)+c  )
        DO i = 1, c
          FQlocal(i) = FQlocal(i) + s * rho * FQelemtemp(i) * Basis(p)
        END DO
      END DO ! fin boucle sur les noeuds
      TIC = CPUTime() - TIC 
      DEBASS = DEBASS + TIC

!------------------------------------------------------------------------------
   END DO                       ! fin de la boucle sur les points d'intégration
!******************************************************************************   

 END SUBROUTINE MANFQelemOPTIDEBUGLOOP

END MODULE DiscreteOperators

!> \} IMNS
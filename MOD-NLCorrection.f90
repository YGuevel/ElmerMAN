!> \ingroup IMNS
!> \}

!> \ingroup IMNS

! # YG@2016
! # yann.guevel@univ-ubs.fr
! # Part of ANM continuation/detection/branch-switching in ELMER FEM
! # Jean-Marc Cadou / Gregory girault / Yann Guevel
! # Institut de Recherche Dupuy de Lome - Universite Bretagne Sud


!------------------------------------------------------------------------------
!>   
!>  Cadou2010 : Reduced Newton + Modified Newton (using K*=[Lt0]_stepbefore first
!>  
!>  
!>
!> TODO : 
!>  
!------------------------------------------------------------------------------
MODULE NLCORRECTION

!------------------------------------------------------------------------------

  USE SolverUtils
  USE DefUtils
! ------------------------------  
  USE DiscreteOperators
  USE HomeMadeSolvers

! ------------------------------

!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
    
CONTAINS






!------------------------------------------------------------------------------------
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
     
! RIKS Order1 : CF Vannucci + pilotage en Dl CF Cadou2010 part 2.4!


!------------------------------------------------------------------------------------





!------------------------------------------------------------------------------------
! linear iterative solver involving a reduced subspace technique
! Reduced Newton + Modified Newton
! Base : UOrthog
! K* : choix Lt0 step before
! K  : assemblage avec solution à corriger
!
! Algo:
!  - Faire un pas avec Lt0j en U0_j
!                   + garder K* = par exemple Lt0j
!                   + garder Série Uorthog pour operateur prolongation P
!  - Nouveau pas en U0_{j+1}:
!    + Assembler Ltemp avec U0_{j+1}
!    + Calcul résidu en U0,L0 : R0 = || L(U0) + Q(U0,U0) - L0F||
!    + Boucle de correction Tant que R0 > ToleranceResidu
!           - operateur reduit:
!           -- Cinter = tP Ltemp
!           -- ksmall = tP Ltemp P   avec Ltemp non triangulé!
!           -- fsmall = tP F   : ici F = R0 ou -R0
!           Step1 : Reduced Newton-Raphson
!           --> COMPUTE dq :
!                          k dq = f - CUcor_{i-1}  avec Ucor_0=0 
!           Step2 : Modified full Newton-Raphson with K* = Lt0j, one may choose K* as ILU with elmer (TODO)
!           --> COMPUTE V :
!                       Kprec V = F - Ktmp(Pdq + U_{i-1})  avec F = R0 ou -RO 
!           Step3 : U_i = U_{i-1} + Pdq + V
!           
!           => COMPUTE NEW RESIDUAL

! IN : Kprec_LU, Knew_nonLU, Uorthog, U_0



!------------------------------------------------------------------------------------






































! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! 
! ! ! ! ! ! NEWTON : Need Residual
! ! ! ! ! 
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! !      
! ! ! ! ! ! AJOUT 20161017 - Corrrection de Lt c pour voir si ca permet de meilleurs tangentes
! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! ! !      R(U0,L0) =  L(U0) + Q(U0,U0) - L0 * F  ! F = 0
! ! ! ! ! !------------------------------------------------------------------------------
! ! ! ! !        ResVec = 0.0_dp
! ! ! ! !        CALL MatrixVectorMultiply(StiffMatrix,UBif,ResVec)
! ! ! ! !        WRITE(*,*)  '|RRB         ||Lt X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))  
! ! ! ! !        
! ! ! ! ! 
! ! ! ! !        IF (Stabilize) THEN
! ! ! ! !               CALL  ANMQABSTAB( abeQAB , UBif, UBif, NSDOFs, FlowPerm, &
! ! ! ! !                          FQManTemp, &
! ! ! ! !                          Uelex , Ueley , Uelez, &
! ! ! ! !                          Velex , Veley , Velez, &                           
! ! ! ! !                          DensityTMP,Material,NodalMuTMP)
! ! ! ! !               WRITE(*,*)  '|RRB  test||Qstab(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))   
! ! ! ! !               ResVec = ResVec - abeQAB
! ! ! ! !        ELSE
! ! ! ! !               CALL  ANMQAB( abeQAA , UBif, UBif, NSDOFs, FlowPerm, &
! ! ! ! !                          FQManTemp, &
! ! ! ! !                          Uelex , Ueley , Uelez, &
! ! ! ! !                          Velex , Veley , Velez, &                           
! ! ! ! !                          DensityTMP,Material)
! ! ! ! !               WRITE(*,*)  '|RRB         ||Q(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAA,abeQAA ))         
! ! ! ! !               ResVec = ResVec - abeQAA
! ! ! ! !        ENDIF
! ! ! ! !        
! ! ! ! !        CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
! ! ! ! !        ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
! ! ! ! !        WRITE(*,*)  '|RESIDUAL as ||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! 
! ! ! ! !        
! ! ! ! !        
! ! ! ! !        
! ! ! ! !        
! ! ! ! !        
! ! ! ! !        
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! 
! ! ! ! !       IF (NEWTCORRECTION) THEN
! ! ! ! !         NEWTITER=0
! ! ! ! !         DO WHILE((ManualResidualU0jstab.GE.PrecRIKS).AND.(NEWTITER.LT.NEWTMAX))
! ! ! ! !          NEWTITER=NEWTITER+1
! ! ! ! !          WRITE(6,*) "----------------------------------------------------------------"
! ! ! ! !          WRITE(6,*) "-                NEWTON CORRECTION ",NEWTITER
! ! ! ! !          WRITE(6,*) "- PrecRIKS     = ",PrecRIKS
! ! ! ! !          WRITE(6,*) "- ResidualSTAB = ",ManualResidualU0jstab
! ! ! ! !          WRITE(6,*) "----------------------------------------------------------------"
! ! ! ! ! 
! ! ! ! ! !        CORRECTOR K Dubar = -R0   Dubar <=> ResVecTMP pour ne pas declarer un tab en plus
! ! ! ! !          ResVecTMP=0.0_dp
! ! ! ! !          CALL CalCondLimCoeff( ResVec , 0._dp, 1, FlowSolution_Init )
! ! ! ! !          CALL HMumps_SolveLinearSystem( StiffMatrix, ResVecTMP, -ResVec, Solver , MUMPSFICH)
! ! ! ! !          CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )
! ! ! ! !          
! ! ! ! !          WRITE(6,*) ' - NEWTON : ||Du||=',DSQRT(dot_product(ResVecTMP,ResVecTMP))
! ! ! ! ! 
! ! ! ! !          CALL CalCondLimCoeff( ResVecTMP , 0._dp, 1, FlowSolution_Init )   
! ! ! ! !          UBif = UBif + ResVecTMP
! ! ! ! !          CALL CalCondLimCoeff( UBif ,   0._dp, 0, lambda0 * FlowSolution_Init )
! ! ! ! !          
! ! ! ! ! 
! ! ! ! ! ! FlowSolution = FlowSolution + CORRECTION
! ! ! ! !          CALL OperatorsLtF( StiffMatrix, ForceVector, UBif, NSDOFs,       &
! ! ! ! !                           MeshVelocity, FlowPerm, MeshPerm, Solver, Model ,     &
! ! ! ! !                           MeshSol, DensitySol, LocalNodes, dt, Transient,       &
! ! ! ! !                           PseudoCompressibilityScale, TempSol ,TempPerm,        &
! ! ! ! !                           Temperature,TempPrev,Bubbles,Stabilize,StabilizeFlag, &           
! ! ! ! !                           CompressibilityModel,DivDiscretization,               &
! ! ! ! !                           GradPDiscretization,NewtonLinearization,              &
! ! ! ! !                           Gravity                                               &
! ! ! ! !                         )
! ! ! ! !          FACTOTODO=.TRUE.
! ! ! ! !          ResVec = 0.0_dp
! ! ! ! !          CALL MatrixVectorMultiply(StiffMatrix,UBif,ResVec)
! ! ! ! !          WRITE(6,*)  '|NEWTON         ||Lt X U0)|| =',DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
! ! ! ! !        
! ! ! ! !          CALL  ANMQABSTAB( abeQAB , UBif, UBif, NSDOFs, FlowPerm, &
! ! ! ! !                          FQManTemp, &
! ! ! ! !                          Uelex , Ueley , Uelez, &
! ! ! ! !                          Velex , Veley , Velez, &                           
! ! ! ! !                          DensityTMP,Material,NodalMuTMP)
! ! ! ! !          WRITE(6,*)  '|NEWTON   ||Qstab(U0,U0)|| =',DSQRT( DOT_PRODUCT( abeQAB,abeQAB ))     
! ! ! ! !          ResVec = ResVec - abeQAB
! ! ! ! !          CALL CalCondLimCoeff( ResVec , 0.0_dp, 1, FlowSolution_Init )       
! ! ! ! !          ManualResidualU0jstab = DSQRT( DOT_PRODUCT( ResVec,ResVec ))           
! ! ! ! !          WRITE(6,*)  '|NEWTON RESSTAB||L(U0) + Q(U0,U0)|| =',ManualResidualU0jstab      
! ! ! ! !          WRITE(RESIDUALFICH,907)  Step,lambda0,NORMUO,ManualResidualU0j,ManualResidualU0jstab
! ! ! ! !          
! ! ! ! !          CALL FLUSH(6)       
! ! ! ! !          CALL CalCondLimCoeff( UBif , 0._dp, 0, FlowSolution_Init )
! ! ! ! !          CALL DefaultDirichletBCs()   
! ! ! ! !        END DO
! ! ! ! !        WRITE(*,*) "---------------------------------------------------------------------"
! ! ! ! !        ENDIF      
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! ! ! ! ! ! !<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>!<!><!><!><!><!><!><!><!><!><!><!><!>
! ! ! ! ! 








END MODULE NLCORRECTION
!> \} IMNS
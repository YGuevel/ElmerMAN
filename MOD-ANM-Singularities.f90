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
MODULE ANMSingularities

!------------------------------------------------------------------------------

  USE SolverUtils
  USE DefUtils

  USE FEMtoolBox    

  
!------------------------------------------------------------------------------
  IMPLICIT NONE
!------------------------------------------------------------------------------
    

CONTAINS





!------------------------------------------------------------------------------
! Compute ratios from a series in order to plot Domb Sykes
! Analyse to see influences of singularities
!------------------------------------------------------------------------------
 SUBROUTINE ANMDombSykesPlot( Lambda0, Lambda, UMan, ANMorder,IFICH ) 
!------------------------------------------------------------------------------
      IMPLICIT NONE
!------------------------------------------------------------------------------    
      REAL(KIND=dp)                 :: Lambda0
      REAL(KIND=dp), DIMENSION(:)   :: Lambda      
      REAL(KIND=dp), DIMENSION(:,:) :: UMan
      INTEGER                       :: ANMorder,IFICH
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: D(ANMorder) , Rat(ANMorder)
      REAL(KIND=dp) :: Num,Denum,tmp
      INTEGER       :: I,J,K
!------------------------------------------------------------------------------    
!     Init
      D   = 0.0_dp
      Rat = 0.0_dp
      
      WRITE(IFICH,*) "Lambda0=",Lambda0
      
!     Compute ratios as in CochMed13      
!     Dk = Ck/Ck-1 = <Xk,Xk> / <Xk-1,Xk>, with X = {U,L}
      DO K=2,ANMorder
        Num   = DOT_PRODUCT( UMan(:,K)   , UMan(:,K) ) + Lambda(k)   * Lambda(k)
        Denum = DOT_PRODUCT( UMan(:,K-1) , UMan(:,K) ) + Lambda(k-1) * Lambda(k)
        D(k) = Num/Denum
      ENDDO

!     Rk = Dk/Dk-1 = (1 - (1+m)/(k) ) / (1 - (1+m)/(k-1) ) 
      DO K=3,ANMorder
        Rat(k) = D(k) / D(k-1)
      ENDDO      
      
!     m if detection ok on last ANM orders
      WRITE(IFICH,*) 2,D(2),0.0_dp,0.0_dp

      DO K=3,ANMorder
        ! m deduced from Rk = Dk/Dk-1  = Ck/Ck-1 / Ck-1/Ck-2 != <Xk,Xk> / <Xk-2,Xk-1>
        tmp   = Rat(k) - 1.0_dp
        Num   = ( (k-2) * k * tmp ) + 1.0_dp
        Denum = (         k * tmp ) + 1.0_dp
        tmp   = Num/denum
        WRITE(IFICH,*) k,D(k),Rat(k),tmp
      ENDDO      
      
      CALL FLUSH(IFICH)


 END SUBROUTINE ANMDombSykesPlot












END MODULE ANMSingularities

!> \} IMNS
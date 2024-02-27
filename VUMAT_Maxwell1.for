C     ***************************************************************
C     Isotropic linear viscoelastic generalized Maxwell model VUMAT subroutine
C     author: Zhu Zhenyu (Institute of Advanced Structure Technology, Beijing Institute of Technology)
C     ***************************************************************

C     should change to vumat function
      subroutine vumat(
C     Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C     Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
      include 'mpif.h'

C     ****************************************************** 
C     Define Material Properties
C     ******************************************************

      DIMENSION PROPS(NPROPS), STATEV(NSTATV), PROPSV(NPROPSV),  
     1  STRAIN(NTENS), DSTRAIN(NTENS), STRESS(NTENS),  
     2  DFGRD0(3,3), DFGRD1(3,3)  
     3  DIMENSIONRELAX(NPROPS), RELAXTIME(NPROPS)  
C     
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0, four=4.0d0)

C     ******************************************************
C     Elastic Modulus and Poisson's Ratio  
C     ******************************************************

      E=props(1) ! Elastic modulus   of spring 0
      aNu=props(2) ! Poisson's ratio
      NMAX=NPROPS/two ! Number of Maxwell elements 
C     Compute elastic constants  
      G=E/(two*(one+aNu))  
      K=E/(three*(one-two*aNu))  

C     ******************************************************
C     Initialize relaxation modulus and relaxation time arrays    
C     ******************************************************

      DO I = 1,NMAX 
      RELAX(I) = PROPS(I + 2)  ! Relaxation modulus  
      RELAXTIME(I) = PROPS(I + NMAX + 1)  ! Relaxation time  
      END DO  
C     Compute elastic stiffness matrix  
      C11 = K + (four/three)*G  
      C12 = K - (two/three)*G  
      C44 = G

C     ******************************************************
C     End Pre work
C     ******************************************************

C     ******************************************************
C     Initialize stress and internal variables     
C     ******************************************************
 
      IF (KSTEP == 1 .AND. KINC == 1) THEN  
      DO I = 1, NTENS  
          STRESS(I) = 0.0  
      END DO  
      DO I = 1, NSTATV  
          STATEV(I) = 0.0  
      END DO  
      END IF  

C     ******************************************************
C     Compute viscous stress increment  
C     ******************************************************

      VISCOUS_STRESS = 0.0  
      DO I = 1, NMAX  
      RELAX_TERM = RELAX(I) * DSTRAIN(1) * EXP(-DTIME / RELAXTIME(I))  
      VISCOUS_STRESS = VISCOUS_STRESS + RELAX_TERM  
      END DO  

C     ******************************************************
C     Compute total stress increment  
C     ******************************************************

      DO I = 1, NTENS  
      DSTRESS(I) = C11 * DSTRAIN(I) + C12 * SUM(DSTRAIN) * IF(I .LE. NDI ,1.0 , 0.0)  
      IF (I .GT. NDI) THEN  
          DSTRESS(I) = C44 * DSTRAIN(I)  
      END IF  
      END DO  
      DSTRESS(1) = DSTRESS(1) + VISCOUS_STRESS 

C     ******************************************************
C     Update stress and internal variables  
C  
      DO I = 1, NTENS  
      STRESS(I) = STRESS(I) + DSTRESS(I)  
      END DO  
C  
      RETURN  
END
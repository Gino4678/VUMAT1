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
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0)

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
      RELAX(I)=PROPS(I + 2)  ! Relaxation modulus  
      RELAXTIME(I)=PROPS(I + NMAX + 1)  ! Relaxation time  
      END DO  

C     ******************************************************
C     End: Pre work
C     ******************************************************
c **************************************************
c * Start: 1st elastic increment
c **************************************************
      if (totalTime.eq.zero) then
      ! if (totalTIme.ge.zero) then
            do km=1, nblock
c volume strain increment
                  dPhi=strainInc(km,1)+strainInc(km,2)
     &                      +strainInc(km,3)
c update normal stress component
                  stressNew(km,1)=stressOld(km,1)
     &                  +GBar0Two*strainInc(km,1)+aLam0*dPhi
                  stressNew(km,2)=stressOld(km,2)
     &                  +GBar0Two*strainInc(km,2)+aLam0*dPhi
                  stressNew(km,3)=stressOld(km,3)
     &                  +GBar0Two*strainInc(km,3)+aLam0*dPhi
c update shear stress component
                  stressNew(km,4)=stressOld(km,4)
     &                  +GBar0Two*strainInc(km,4)
                  stressNew(km,5)=stressOld(km,5)
     &                  +GBar0Two*strainInc(km,5)
                  stressNew(km,6)=stressOld(km,6)
     &                  +GBar0Two*strainInc(km,6)
            end do
c **************************************************
c * End: 1st elastic increment
c **************************************************
      else
c **************************************************
c * Start: following increments
c **************************************************
            do km=1, nblock
                  ! if (totalTime.gt.8.0) then
                  ! write(*,*) 'Please input an integer:'
                  ! read(*,*) tmpread
                  ! end if
c read the visco strain 1 at start of this increment
                  do i=1, 6
                        e1Old(i)=stateOld(km,i)
                  end do
                  phi1Old=stateOld(km,7)
c get strain increment
                  do i=1, 6
                        dStrain(i)=strainInc(km,i)
                  end do
                  call TensorDecomp(dStrain, dE, dPhi)
c calculate stress increment due to visco 1
                  pre1=one-exp(-dt/GTauBar)
                  pre2=GTauBar/dt
                  do i=1, 6
                        dE1(i)=pre1*(pre2*dE(i)-e1Old(i))
                  end do
                  pre3=one-exp(-dt/aKTauBar)
                  pre4=aKTauBar/dt
                  dPhi1=pre3*(pre4*dPhi-phi1Old)
c update stress
                  do i=1, 6
                        dS(i)=GBarInfTwo*dE(i)+GBar1Two*dE1(i)
                  end do
                  dP=aKBarInf*dPhi+aKBar1*dPhi1
                  do i=1, 3
                        stressNew(km,i)=stressOld(km,i)+dS(i)+dP
                  end do
                  do i=4, 6
                        stressNew(km,i)=stressOld(km,i)+dS(i)
                  end do
c update state variables
                  do i=1, 6
                        stateNew(km,i)=e1Old(i)+dE1(i)
                  end do
                  stateNew(km,7)=phi1Old+dPhi1
            end do
c **************************************************
c * End: following increments
c **************************************************
      end if
c
      return
      end
c
      subroutine TensorDecomp(a,aDev,aTrace)
c
c     param:
c     a: an symmetric tensor in Voigt notation
c
c     return:
c     aDev: devatoric part of a in Voigt notation
c     aTrace: trace of a
c
      include 'vaba_param.inc'
c
      dimension a(6), aDev(6)
c
      parameter three=3.0d0
c
      aTrace=a(1)+a(2)+a(3)
c
      do i=1, 3
            aDev(i)=a(i)-aTrace/three
      end do
      do i=4, 6
            aDev(i)=a(i)
      end do
c
      return
      end      
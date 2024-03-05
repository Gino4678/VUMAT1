C     ***************************************************************
C     Isotropic linear viscoelastic generalized Maxwell model 
C     VUMAT subroutine
C     Author: Zhu Zhenyu (Institute of Advanced Structure Technology,
C                         Beijing Institute of Technology)
C     ***************************************************************
c should change to vumat function
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
      include 'mpif.h'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      ! dimension viscoStrain1Old(6), devViscoStrain1Old(6)
      ! dimension dDevViscoStrain1(6)
      ! dimension dStrain(6), dDevStrain(6)
      ! dimension dDevViscoStress(6), dDevElasStress(6)
      dimension e1Old(6),e2Old(6)
      dimension dStrain(6)
      dimension dE(6)   ! Partial strain increment
      dimension dE1(6), dE2(6)
      dimension dS(6)   ! Partial stress increment
c
      character*80 cmname
c
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0)
c ******************************************************
c * Start: Pre work
c ******************************************************
c
c read material model constants from props
      E0=props(1) ! young's modulus of maxwell model
      aNu=props(2) ! Poisson's ratio

      ag1=props(3) ! g1
      ak1=ag1 ! k1
      ag2=props(4) ! g2
      ak2=ag2 ! k2
      tau1=props(5) ! tau1
      tau2=props(6) ! tau2
! possion's ratio which is persumed to be constant
! thus the variation of bulk modulus is the same as the shear modulus
c
      EBarInf=E0/(one+ag1+ag2) ! long term young's modulus
      GBarInf=EBarInf/(two*(one+aNu)) ! long term shear modulus
      aKBarInf=EBarInf/(three*(one-two*aNu)) ! long term bulk modulus
      aLam0=(E0*aNu)/((one+aNu)*(one-two*aNu)) !instantaneous lame constant
      
      GBar1=ag1*EBarInf/(two*(one+aNu)) ! visco 1 shear modulus
      aKBar1=ag1*EBarInf/(three*(one-two*aNu)) ! visco 1 bulk modulus

      GBar2=ag2*EBarInf/(two*(one+aNu)) ! visco 2 shear modulus
      aKBar2=ag2*EBarInf/(three*(one-two*aNu)) ! visco 2 bulk modulus
c
c ******************************************************
c * End: Pre work
c ******************************************************
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
c read the visco strain 1 at start of this increment
                  do i=1, 6
                        e1Old(i)=stateOld(km,i)
                  end do
                  phi1Old=stateOld(km,7)
c read the visco strain 2 at start of this increment
                  do i=1, 6
                        e2Old(i)=stateOld(km,i+7)
                  end do
                  phi2Old=stateOld(km,14)
c get strain increment
                  do i=1, 6
                        dStrain(i)=strainInc(km,i)
                  end do
                  call TensorDecomp(dStrain, dE, dPhi)
c calculate stress increment due to visco 1
                  pre1=one-exp(-dt/tau1)
                  pre2=tau1/dt
                  do i=1, 6
                        dE1(i)=pre1*(pre2*dE(i)-e1Old(i))
                  end do
                  pre3=one-exp(-dt/tau1)
                  pre4=tau1/dt
                  dPhi1=pre3*(pre4*dPhi-phi1Old)
c calculate stress increment due to visco 2
                  pre5=one-exp(-dt/tau2)
                  pre6=tau2/dt
                  do i=1, 6
                        dE2(i)=pre1*(pre2*dE(i)-e2Old(i))
                  end do
                  pre7=one-exp(-dt/tau2)
                  pre8=tau2/dt
                  dPhi2=pre3*(pre4*dPhi-phi2Old)
c update stress
                  do i=1, 6
                        dS(i)=GBarInf*Two*dE(i)+GBar1*Two*dE1(i)
     1                                         +GBar2*Two*dE2(i)
                  end do
                  dP=aKBarInf*dPhi+aKBar1*dPhi1+aKBar2*dPhi2
                  do i=1, 3
                        stressNew(km,i)=stressOld(km,i)+dS(i)+dP
                  end do
                  do i=4, 6
                        stressNew(km,i)=stressOld(km,i)+dS(i)
                  end do
c update state variables
                  do i=1, 6
                        stateNew(km,i)=e1Old(i)+dE1(i)
                        stateNew(km,i+7)=e2Old(i)+dE2(i)
                  end do
                  stateNew(km,7)=phi1Old+dPhi1
                  stateNew(km,14)=phi2Old+dPhi2
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
c     a: an symmetric tensor 
c
c     return:
c     aDev: devatoric part
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
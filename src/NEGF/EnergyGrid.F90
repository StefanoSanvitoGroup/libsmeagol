! 
! Copyright (c) Authors:
! Ivan Rungger and Stefano Sanvito
! Trinity College Dublin, Ireland
! October 2008 
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! THE SUBROUTINES
!                   ENERGYGRID_SELFENERGIES_REAL2,
!                   ENERGYGRID_SELFENERGIES_IMAG2,
!                   ENERGYGRID_TRANSM,
!                   ENERGIES_IMAG,
!                   ENERGY_POINTS2,
!                   ALLOCATE_SIGMA_SINGLE,
!                   DEALLOCATE_ENERGYGRID_SELFENERGIES2,
!                   DEALLOCATE_ENERGYGRID_SELFENERGIES,
!                   READ_SELFENERGIES,
!                   WRITE_SELFENERGIES,
!                   CHECK_FILE,
!                   ALLOCATE_SIGMA2,
!                   WRITE_SINGLE_SELFENERGY,
!                   READSINGLE_SELFENERGY,
!                   SINGLE_SELFENERGY,
!                   ALL_SELFENERGIES,
!                   ENERGYGRID_ADAPT2,
!                   RECURSIVE_ENERGYGRID,
!                   REALLOC_GLOBAL,
!                   CALCULATE_WEIGHTS,
!                   REDISTRIBUTE_CONST,
!                   ADD_ENERGYPOINTS,
!                   REALLOC_MYDEPTH,
!                   DUPLICATE_ENERGYGRID,
!                   REALLOCATE_ENERGYGRID,
!                   DEALLOCATE_SELFENERGIES,
!                   GET_SELFENERGIES,
!                   WRITE_SELFENERGIES_TRANSMISSION,
!                   READ_SELFENERGIES_TRANSMISSION,
!                   SELFENERGY_K,
!                   HSFOURIER,
!                   SURFACEGF,
!                   EXTRACTSIGMA,
!                   CALCULATE_AND_WRITE_SELFENERGIES  
! AND
! THE MODULE
!                   MENERGYGRID  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
module mEnergyGrid

  use mConstants
  use mTypes
    
  implicit none

  private


  type(EnergyGridType), public :: ERealGrid
  type(EnergyGridType), public :: EImagGrid
  type(EnergyGridType), public :: ETransmGrid

!  type(EnergyPoint), allocatable, save, public :: ErealGlobal(:)


  real(kdp), parameter :: ToleranceT=12.0_kdp 

  public :: energygrid_selfenergies_real2
  public :: energygrid_selfenergies_imag2
  public :: energygrid_transm
  public :: get_selfenergies
  public :: deallocate_selfenergies
  public :: deallocate_energygrid_selfenergies
  public :: deallocate_energygrid_selfenergies2
  public :: allocate_sigma_single
  public :: extractsigma

  contains




  subroutine energygrid_selfenergies_real2(label,ef_lead,v,Temp,iter,node,nnodes,n1,nspin,Delta,ik,hgeneral,sgeneral,istep, inicoor, vinitial,vfinal,nEneRealBase,nk,nleadslr,leadsdim,storesigmai,LeadsVoltageShift,sigmatodiski)

    use negfmod
    use mMPI_NEGF
    
    character(LEN=*), intent(in) ::  label
    logical, intent(in) :: sigmatodiski
    integer, intent(in) ::  iter,node,nnodes,n1,nspin,ik,nk,istep,nEneRealBase,nleadslr,leadsdim(nleadslr),storesigmai, inicoor
    real(kdp), intent(in) ::  ef_lead,v,Temp,Delta,vinitial,vfinal,LeadsVoltageShift(nleadslr)
    type(matrixTypeGeneral), intent(in) :: hgeneral(nspin),sgeneral
    real(kdp) einitial,efinal

    integer ispin,ie,ikl,il,nene,neneglobal
    character fsigma*250,char_bias*7
    logical,save ::  fsigma_exists
    logical, external :: leqi


    ERealGrid%nLeads=nleadslr
    if(.not.allocated(ERealGrid%leadsTotalDim))allocate(ERealGrid%leadsTotalDim(ERealGrid%nLeads))

    if(nspin<=2)then
      ERealGrid%leadsTotalDim=leadsdim
      ERealGrid%nSpin=nspin
    else
      ERealGrid%leadsTotalDim=leadsdim * 2
      ERealGrid%nSpin=1
    endif
  
    ERealGrid%v=v
    ERealGrid%nk=nk
    ERealGrid%deltasigma=deltaimag
    ERealGrid%sLabel=trim(label)
    ERealGrid%InfoSigma=storesigmai
    ERealGrid%SigmaSuffix=".EREAL_SIGMA"

    If(setene)then
      einitial= setemin/13.6057D0 + ef_lead
      efinal= setemax/13.6057D0 + ef_lead
    else
      einitial = Ef_Lead-abs(V)*0.5_kdp-ToleranceT*Temp
      efinal = Ef_Lead+abs(V)*0.5_kdp+ToleranceT*Temp
    endif

    if (V == 0_kdp) then
      ERealGrid%nEnergies=0
    else

      IF (leqi(gridmethod,"adaptivegrid")) then
        ERealGrid%GridType=1
      else
        ERealGrid%GridType=0
      endif

      if(ERealGrid%GridType==1)then
        if(ITER.EQ.1.and.node.eq.0 .and. istep .eq. inicoor )  write(6,'(a)') "energygrid_selfenergies_real : Using adaptive grid"

        ERealGrid%nk=1

        ERealGrid%nEnergies=INT((ABS(einitial-efinal)/(deltaini/deltatode)) /Nnodes)
        If (ERealGrid%nEnergies.EQ.0) ERealGrid%nEnergies=1
        ERealGrid%nEnergiesGlobal=ERealGrid%nEnergies*Nnodes

        if (allocated(ERealGrid%e)) deallocate(ERealGrid%e) !!
        allocate(ERealGrid%e(ERealGrid%nEnergies))
        if (allocated(ERealGrid%w)) deallocate(ERealGrid%w) !!
        allocate(ERealGrid%w(ERealGrid%nEnergies))
        if (allocated(ERealGrid%ig)) deallocate(ERealGrid%ig) !!
        allocate(ERealGrid%ig(ERealGrid%nEnergies))
        nene=ERealGrid%nEnergies
        neneglobal=ERealGrid%nEnergiesGlobal
        call energygrid_adapt2(N1,NSPIN, einitial,efinal,Delta, V,nene,neneglobal, ik,hgeneral,sgeneral,nk,storesigmai,nleadslr,LeadsVoltageShift)

      else

!        if(ITER.EQ.1.and.(node.eq.0))  write(6,'(a)') "energygrid_selfenergies_real : Using Gaussian grid"

!        write(12347,*)"istep=",istep,iter,ik
        if (ITER.EQ.1.AND.ik.EQ.1 .and. istep .eq. inicoor) THEN

          ERealGrid%nEnergies=NINT(((abs(v)+2_kdp*ToleranceT*Temp)/(MAX(abs(vinitial),abs(Vfinal))+2_kdp*ToleranceT*Temp)) *nEneRealBase/Nnodes)
          If (ERealGrid%nEnergies.EQ.0.and.nEneRealBase.gt.0) ERealGrid%nEnergies=1
          ERealGrid%nEnergiesGlobal=ERealGrid%nEnergies*Nnodes

          if(ERealGrid%nEnergies>0)then

            if (allocated(ERealGrid%e)) deallocate(ERealGrid%e) !!
            allocate(ERealGrid%e(ERealGrid%nEnergies))
            if (allocated(ERealGrid%w)) deallocate(ERealGrid%w) !!
            allocate(ERealGrid%w(ERealGrid%nEnergies))
            if (allocated(ERealGrid%ig)) deallocate(ERealGrid%ig) !!
            allocate(ERealGrid%ig(ERealGrid%nEnergies))
!            allocate(ErealGlobal(ERealGrid%nEnergiesGlobal))

            call energy_points2(einitial,efinal, V,ERealGrid%nEnergies,ERealGrid%e,ERealGrid%w,ERealGrid%ig,storesigmai)

            if(mynode_inverse==0) call allocate_sigma2(storesigmai,LeadsVoltageShift,ERealGrid)
          endif

        endif

        if (mynode_inverse==0.and.ITER.EQ.1 .and. istep .eq. inicoor) THEN
!          write(12347,*)"istep=",istep,v,iter,ik,ERealGrid%nEnergies,nspin,ERealGrid%nLeads

          write(char_bias,'(f7.4)')v*13.6057_kdp
          fsigma=trim(label)//'.V_'//TRIM(ADJUSTL(char_bias))//'.EREAL_SIGMA'

          if(ik==1)then
            call check_file(fsigma,fsigma_exists)
          endif

          if(storesigmai==0.or.storesigmai==2)then
            if(.not.fsigma_exists)then
              if(storesigmai==0)then
                call all_selfenergies(ik,1,ERealGrid%nEnergies,1,ERealGrid%nspin,ERealGrid,0)
                if(sigmatodiski) call write_selfenergies(ik,fsigma,ERealGrid)
              else
                call calculate_and_write_selfenergies(ik,fsigma,ERealGrid)
              endif
            else
              if(storesigmai==0)then
                call read_selfenergies(ik,fsigma,ERealGrid)
              endif
            endif
          endif

        endif
      endif

    endif 

  end subroutine energygrid_selfenergies_real2


  subroutine energygrid_selfenergies_imag2(label,iter,ik,nenerg1,nenerg2, npoles, istep, inicoor, nspin,V,T,EFermi,r0,eb,nleadslr,leadsdim,nk,storesigmai,LeadsVoltageShift,sigmatodisk)


    use mMPI_NEGF

    implicit none

    logical, intent(in) :: sigmatodisk
    character(LEN=*), intent(in) ::  label
    integer, intent(in) ::  nleadslr,leadsdim(nleadslr),nk,storesigmai
    real(kdp), intent(in) ::  LeadsVoltageShift(nleadslr)

    integer iter,ik,nenerg1,nenerg2,npoles, NenergImNode,nspin
    double precision V,T,EFermi,r0,eb
    integer ikl,il,ie,ispin, istep, inicoor
    integer nline2,ncirc2,npoles2,nenesum,deltane,nadd,mpierror

    character fsigma*250,char_bias*7
    logical,save ::  fsigma_exists
    logical, external :: leqi

    EImagGrid%nLeads=nleadslr
    if(.not.allocated(EImagGrid%leadsTotalDim))allocate(EImagGrid%leadsTotalDim(EImagGrid%nLeads))

    if(nspin<=2)then
      EImagGrid%leadsTotalDim=leadsdim
      EImagGrid%nSpin=nspin
    else
      EImagGrid%leadsTotalDim=leadsdim * 2
      EImagGrid%nSpin=1
    endif
  
    EImagGrid%v=v
    EImagGrid%deltasigma=0.0_kdp
    EImagGrid%nk=nk
    EImagGrid%GridType=0
    EImagGrid%sLabel=trim(label)
    EImagGrid%InfoSigma=storesigmai
    EImagGrid%SigmaSuffix=".EIMAG_SIGMA"


    IF (ITER .EQ. 1 .AND. ik .EQ. 1 .and. istep .eq. inicoor) THEN


      nline2=nenerg1
      ncirc2=nenerg2
      npoles2=npoles

      nenesum=nline2+ncirc2+2*npoles2
      if(nenesum.ne.0)then
!        write(*,*)"nheads=",nheads
        deltane=MOD(nenesum,nheads)
        if(deltane.ne.0)deltane=nheads-deltane
        nadd=deltane/4
        nline2=nline2+nadd
        ncirc2=ncirc2+nadd
        npoles2=npoles2+nadd
        nenesum=nline2+ncirc2+2*npoles2
        deltane=MOD(nenesum,nheads)
        if(deltane.ne.0)deltane=nheads-deltane
        ncirc2=ncirc2+deltane
        nenesum=nline2+ncirc2+2*npoles2
        deltane=MOD(nenesum,nheads)
      endif

      nenerg1=nline2
      nenerg2=ncirc2
      npoles=npoles2

      EImagGrid%nEnergies=(nline2+ncirc2+2*npoles2)/nheads
      EImagGrid%nEnergiesGlobal=EImagGrid%nEnergies*nheads

      if (allocated(EImagGrid%e)) deallocate(EImagGrid%e) !!
      allocate(EImagGrid%e(EImagGrid%nEnergies))
      if (allocated(EImagGrid%w)) deallocate(EImagGrid%w) !!
      allocate(EImagGrid%w(EImagGrid%nEnergies))
      if (allocated(EImagGrid%ig)) deallocate(EImagGrid%ig) !!
      allocate(EImagGrid%ig(EImagGrid%nEnergies))
      call energies_imag(nline2,ncirc2,npoles2,EImagGrid%nEnergies,EImagGrid%nEnergies,R0,EB,EFermi,V,T,EImagGrid%e,EImagGrid%w,EImagGrid%ig,storesigmai)

      if(mynode_inverse==0) call allocate_sigma2(storesigmai,LeadsVoltageShift,EImagGrid)
    ENDIF

    if (mynode_inverse==0.and.ITER.EQ.1 .and. istep .eq. inicoor) THEN

      write(char_bias,'(f7.4)')v*13.6057_kdp
      fsigma=trim(label)//'.V_'//TRIM(ADJUSTL(char_bias))//'.EIMAG_SIGMA'
!
      if(ik==1)then
        call check_file(fsigma,fsigma_exists)
      endif

      if(storesigmai==0.or.storesigmai==2)then
        if(.not.fsigma_exists)then
          if(storesigmai==0)then
            call all_selfenergies(ik,1,EImagGrid%nEnergies,1,EImagGrid%nspin,EImagGrid,0)
            if(sigmatodisk)call write_selfenergies(ik,fsigma,EImagGrid)
          else
            call calculate_and_write_selfenergies(ik,fsigma,EImagGrid)
          endif
        else
          if(storesigmai==0)then
            call read_selfenergies(ik,fsigma,EImagGrid)
          endif
        endif
      endif

    endif
 



  end subroutine energygrid_selfenergies_imag2



  subroutine energygrid_transm(label,nspin,V,T,EFermi,nleadslr,leadsdim,nk,LeadsVoltageShift,deltaene,deltaenes,deltaenebig)

    use mMPI_NEGF
    use negfmod, only: deauto,trcde,trcef,tenergi,tenergf,deltaimag,nenet,storesigma,RytoeV

    implicit none

    character(LEN=*), intent(in) ::  label
    integer, intent(in) ::  nleadslr,leadsdim(nleadslr),nk,nspin
    real(kdp), intent(in) ::  LeadsVoltageShift(nleadslr),EFermi,v,t
    real(kdp), intent(out) ::  deltaene,deltaenes,deltaenebig

    real(kdp) kbtoltrans,EnerInitial,EnerFinal

    integer ie,ntglobal,nt

    ETransmGrid%nLeads=nleadslr
    if(.not.allocated(ETransmGrid%leadsTotalDim))allocate(ETransmGrid%leadsTotalDim(ETransmGrid%nLeads))
    ETransmGrid%leadsTotalDim=leadsdim
  
    ETransmGrid%v=v
    ETransmGrid%nSpin=nspin
    ETransmGrid%deltasigma=deltaimag
    ETransmGrid%nk=nk
    ETransmGrid%GridType=0
    ETransmGrid%sLabel=trim(label)
    ETransmGrid%GridType=storesigma
    ETransmGrid%SigmaSuffix=".TRC_SIGMA"


    if(deauto)then
      if(trcde.le.0D0)trcde=1D-3
      kbtoltrans=7D0
      EnerInitial=EFermi-0.5D0 * abs(v) - kbtoltrans * T
      EnerFinal=EFermi+0.5D0 * abs(v) + kbtoltrans * T
      ntglobal=INT((EnerFinal-EnerInitial- abs(v))/trcde)
!      write(*,*)"deltaemin=",2D0 * kbtoltrans * T,ntglobal
    elseif(trcef)then
      EnerInitial=tenergi+EFermi
      EnerFinal=tenergf+EFermi
      ntglobal=NeneT
    else
      EnerInitial=tenergi
      EnerFinal=tenergf
      ntglobal=NeneT
    endif

    if(trcde.gt.0D0)then
      ntglobal=INT((EnerFinal-EnerInitial)/trcde)
    endif
    if(myhead .EQ.0) then
      write(*,*)"Calculating transmission coefficient"
!      write(*,*)"eif=",EnerInitial,EnerFinal,trcde,ntglobal,EFermi,T
    endif

    nt=ntglobal/nheads
    if(nheads*nt < ntglobal)nt=nt+1
    ntglobal=nheads*nt

    ETransmGrid%nEnergies=nt
    ETransmGrid%nEnergiesGlobal=ETransmGrid%nEnergies*nheads

    allocate(ETransmGrid%e(ETransmGrid%nEnergies))
    allocate(ETransmGrid%w(ETransmGrid%nEnergies))
    allocate(ETransmGrid%ig(ETransmGrid%nEnergies))

    IF (ETransmGrid%nEnergiesGlobal.EQ.1.or.(EnerFinal-EnerInitial==0)) THEN
      deltaenes=0.0_kdp
      deltaenebig=deltaenes * nheads
      deltaene=deltaenes
      do ie=1,ETransmGrid%nEnergies
        ETransmGrid%e(ie)=EnerInitial + (ie-1)*deltaenebig + myhead * deltaenes
        ETransmGrid%w(ie)=1.0D0/(RytoeV * (ETransmGrid%nEnergiesGlobal))
!xxx: check the rytoev, can we remove it? Ideally we should keep everything in Ry
        ETransmGrid%ig(ie)=1
      enddo
    ELSE
      deltaenes=ABS(EnerFinal-EnerInitial)/(1D0 * (ETransmGrid%nEnergiesGlobal-1))
      deltaenebig=deltaenes * nheads
      deltaene=deltaenes
      do ie=1,ETransmGrid%nEnergies
        ETransmGrid%e(ie)=EnerInitial + (ie-1)*deltaenebig + myhead * deltaenes
        ETransmGrid%w(ie)=deltaenes
!        ETransmGrid%w(ie)=deltaenes/(RytoeV * ABS(EnerFinal-EnerInitial))
!xxx: check the whether we should rescale by energy window sometimes
        ETransmGrid%ig(ie)=1
      enddo
      if(myhead==0)ETransmGrid%w(1)=ETransmGrid%w(1)*0.5D0
      if(myhead==nheads-1)ETransmGrid%w(ETransmGrid%nEnergies)=ETransmGrid%w(ETransmGrid%nEnergies)*0.5D0
!      if(myhead==0)write(12347,*)"minene=",ETransmGrid%e(1)
!      if(myhead==nheads-1)write(12347,*)"maxene=",ETransmGrid%e(ETransmGrid%nEnergies)
    ENDIF

    if(myhead==0)then
      allocate(ETransmGrid%eGlobal(ETransmGrid%nEnergiesGlobal))

      IF (ETransmGrid%nEnergiesGlobal.EQ.1) THEN
        ETransmGrid%eGlobal(1)=EnerInitial
      ELSE
        do ie=1,ETransmGrid%nEnergiesGlobal
          ETransmGrid%eGlobal(ie)=EnerInitial + (ie-1)*deltaenes
        enddo
      ENDIF
    endif




  end subroutine energygrid_transm




  SUBROUTINE energies_imag(Nenerg1,Nenerg2,NPOLES,NenergIm, NenergImNode,R0,EB,Ef_fermi,V,T,Ei,W,iGlobal,infosigma)


  use mMPI_NEGF
  include "const2.h"

  integer, intent(in) :: Nenerg1,Nenerg2,NPOLES,NenergIm,NenergImNode,infosigma
  real(kdp), intent(in) ::  eb,Ef_fermi,v,t
  real(kdp), intent(inout) ::  r0
  complex(kdp), intent(inout):: ei(nenergim),w(nenergim)
  integer     , intent(inout):: iGlobal(nenergim)

  real(kdp) :: Gamma,r1,delta,theta0
  integer nenerg,i,IPOLES
  real(kdp):: X1(Nenerg1),WX1(Nenerg1),X2(Nenerg2),WX2(Nenerg2)
  DOUBLE COMPLEX, DIMENSION (Nenerg1+Nenerg2+2*NPOLES) :: Eiaux,Waux
! ***********************************************************************
!  Calculate the radius of the semicircle and its centre
!
  Gamma=20.D0*T
  Nenerg=Nenerg1+Nenerg2
  Delta=2*PI*NPOLES*T
  R1=((Ef_fermi-DABS(V/2.D0)-Gamma-EB)**2+Delta**2)/ (2.D0*(Ef_fermi-DABS(V/2.D0)-Gamma-EB))
  R0=EB+R1
  theta0=DATAN(Delta/(Ef_fermi-DABS(V/2.D0)-gamma-R0))


  CALL GAULEG(Ef_fermi+DABS(V/2.D0)+Gamma, Ef_fermi-DABS(V/2.D0)-Gamma,X1,WX1,Nenerg1)
  DO I=1,Nenerg1
   Eiaux(I)=DCMPLX(X1(I))+zi*Delta
   Waux(I)=DCMPLX(WX1(I))
  ENDDO

! *********************************************************************
! Calculate the points of integration for the semicircle in the complex
! plane. Uses a Gauss-Legendre algorithm
!
  CALL GAULEG(theta0,PI,X2,WX2,Nenerg2)
  DO I=1,Nenerg2
   Eiaux(Nenerg1+I)=DCMPLX(R0)+DCMPLX(R1)*CDEXP(zi*DCMPLX(X2(I)))
   Waux(Nenerg1+I)=DCMPLX(WX2(I))
  ENDDO

! *********************************************************************
! Calculate the Matsubara frequencies (poles of the Fermi-Dirac
! distribution) and the corresponding residues
!
  DO IPOLES=1,NPOLES
   Eiaux(NEnerg+IPOLES)=zi*(2*(IPOLES-1)+1)*PI*T +  Ef_fermi + V/2.D0
   Waux(NEnerg+IPOLES)=2*zi*PI*T
   Eiaux(NEnerg+NPOLES+IPOLES)=  zi*(2*(IPOLES-1)+1)*PI*T + Ef_fermi - V/2.D0
   Waux(NEnerg+NPOLES+IPOLES)=Waux(NEnerg+IPOLES)
  ENDDO


  if(infosigma==0.or.infosigma==1)then
    DO i=1,NenergImNode
     Ei(i)=Eiaux(myhead * NenergImNode + i)
     W(i)=Waux(myhead * NenergImNode + i)
     iGlobal(i)=myhead * NenergImNode + i
    ENDDO 


  else
    DO i=1,NenergImNode
     Ei(i)=Eiaux(myhead + (i-1) * nheads+1)
     W(i)=Waux(myhead + (i-1) * nheads+1)
     iGlobal(i)=myhead + (i-1) * nheads+1
    ENDDO 
  endif

  END SUBROUTINE energies_imag


  SUBROUTINE energy_points2(EnergI,EnergF, V,nene,energies,weights,iGlobals,infosigma)

    use mMPI_NEGF

    IMPLICIT NONE
     
    include "const2.h"
      
    integer, intent(in):: nene,infosigma
    DOUBLE PRECISION, intent(in):: EnergI,EnergF,V
    complex(kdp), intent(out):: energies(nene),weights(nene)
    integer     ,  intent(out):: iGlobals(nene)

    DOUBLE PRECISION,ALLOCATABLE,DIMENSION (:) ::  EX,EW
    INTEGER i,ihead,iglobal

    ALLOCATE(EX(ERealGrid%nEnergiesGlobal),EW(ERealGrid%nEnergiesGlobal))
    CALL GAULEG(EnergI,EnergF,EX,EW, ERealGrid%nEnergiesGlobal)


    if(infosigma==0.or.infosigma==1)then
      DO i=1,nene
        iglobal=myhead * nene+i
        energies(i)=EX(iglobal)
        weights(i)=EW(iglobal)
        iGlobals(i)=iglobal
      ENDDO
    else
      DO i=1,nene
        iglobal=myhead + (i-1) * nheads+1
        energies(i)=EX(iglobal)
        weights(i)=EW(iglobal)
        iGlobals(i)=iglobal
      ENDDO
    endif

    DEALLOCATE(EX,EW)

  END SUBROUTINE energy_points2

  subroutine allocate_sigma_single(sigma,n,mynode,storesigmai,energy,vshift,deltaimag)

    integer, intent(in) :: n,mynode,storesigmai
    real(kdp), intent(in) :: vshift,deltaimag
    complex(kdp), intent(in) :: energy
    type(SelfEnergyType), intent(inout) :: sigma
  
    allocate(sigma%sigma(n,n))
    sigma%n=n
    sigma%node=mynode
    sigma%InfoSigma=storesigmai
    sigma%e=energy-vshift+(0.0_kdp,1.0_kdp)*deltaimag
  
  end subroutine allocate_sigma_single

  subroutine deallocate_energygrid_selfenergies2(energygrid1)

    type(EnergyGridType), intent(inout) :: energygrid1
  
    integer il
  
    if(allocated(energygrid1%sigma))then
      do il=1,energygrid1%nLeads
        if(allocated(energygrid1%sigma(il,1,1,1)%sigma)) deallocate(energygrid1%sigma(il,1,1,1)%sigma)
      enddo
      deallocate(energygrid1%sigma)
    endif
  
    if(allocated(energygrid1%e))deallocate(energygrid1%e)
    if(allocated(energygrid1%w))deallocate(energygrid1%w)
    if(allocated(energygrid1%ig))deallocate(energygrid1%ig)
    if(allocated(energygrid1%eGlobal))deallocate(energygrid1%eGlobal)

  end subroutine deallocate_energygrid_selfenergies2


  subroutine deallocate_energygrid_selfenergies(energygrid1)

    type(EnergyGridType), intent(inout) :: energygrid1
  
    integer ik,ispin,ie,il
  
    if(allocated(energygrid1%sigma))then
      if(energygrid1%sigma(1,1,1,1)%InfoSigma==0)then
        do ik=1,energygrid1%nk
          do ispin=1,energygrid1%nspin
            do ie=1,energygrid1%nEnergies
              do il=1,energygrid1%nLeads
                deallocate(energygrid1%sigma(il,ie,ispin,ik)%sigma)
              enddo
            enddo
          enddo
        enddo
      endif
      deallocate(energygrid1%sigma)
    endif
  
    deallocate(energygrid1%e)
    deallocate(energygrid1%w)
    deallocate(energygrid1%ig)

  end subroutine deallocate_energygrid_selfenergies


  subroutine read_selfenergies(ik,fsigma,energygrid1)

    use mMPI_NEGF

    integer, intent(in) :: ik
    character(LEN=*), intent(in) ::  fsigma
    type(EnergyGridType), intent(inout) :: energygrid1

    integer ispin,ie,il,ihead,bytes_sigma,sigma_io
    integer MPIerror,i1,i2
    type(SelfEnergyType), allocatable :: sigma_buffer(:)
    integer rec1
#ifdef MPI
    INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

    if(myhead==0)then

      sigma_io=22349

      allocate(sigma_buffer(energygrid1%nLeads))
      do il=1,energygrid1%nLeads
        call allocate_sigma_single(sigma_buffer(il),energygrid1%leadsTotalDim(il),myhead,0,(0.0_kdp,0.0_kdp),0.0_kdp,0.0_kdp)
      enddo

      bytes_sigma=0
      do il=1,energygrid1%nLeads
        bytes_sigma=bytes_sigma + 16 *  sigma_buffer(il)%n**2
      enddo

      OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='OLD', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)

      rec1=energygrid1%nEnergiesGlobal*energygrid1%nspin*(ik-1)
      do ie=1,energygrid1%nEnergies
        do ispin=1,energygrid1%nspin
          rec1=rec1+1
          read(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

          do il=1,energygrid1%nLeads
            energygrid1%sigma(il,ie,ispin,ik)%sigma=sigma_buffer(il)%sigma
          enddo
        enddo
      enddo


#ifdef MPI
      do ihead=1,nheads-1

        do ie=1,energygrid1%nEnergies
          do ispin=1,energygrid1%nspin

            rec1=rec1+1
            read(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

            do il=1,energygrid1%nLeads
              call MPI_SEND(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex,ihead,1,inverseheads_comm,MPIerror)
            enddo

          enddo
        enddo
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)

      enddo
#endif

      close(sigma_io)

      do il=1,energygrid1%nLeads
        deallocate(sigma_buffer(il)%sigma)
      enddo
      deallocate(sigma_buffer)

#ifdef MPI
    else

      allocate(sigma_buffer(energygrid1%nLeads))
      do il=1,energygrid1%nLeads
        call allocate_sigma_single(sigma_buffer(il),energygrid1%leadsTotalDim(il),myhead,0,(0.0_kdp,0.0_kdp),0.0_kdp,0.0_kdp)
      enddo

      do ihead=1,nheads-1
        if(myhead.eq.ihead)then

          do ie=1,energygrid1%nEnergies
            do ispin=1,energygrid1%nspin
              do il=1,energygrid1%nLeads
                call MPI_RECV(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex, 0, 1, inverseheads_comm, istatus, MPIerror)
                energygrid1%sigma(il,ie,ispin,ik)%sigma=sigma_buffer(il)%sigma
              enddo
            enddo
          enddo

        endif
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo

      do il=1,energygrid1%nLeads
        deallocate(sigma_buffer(il)%sigma)
      enddo
      deallocate(sigma_buffer)

#endif

    endif



  end SUBROUTINE read_selfenergies


  subroutine write_selfenergies(ik,fsigma,energygrid1)

    use mMPI_NEGF

    integer, intent(in) :: ik
    character(LEN=*), intent(in) ::  fsigma
    type(EnergyGridType), intent(in) :: energygrid1

    integer ispin,ie,il,ihead,bytes_sigma,sigma_io,rec1
    integer MPIerror,i1,i2
    type(SelfEnergyType), allocatable :: sigma_buffer(:)
#ifdef MPI
    INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

    if(myhead==0)then

      allocate(sigma_buffer(energygrid1%nLeads))
      do il=1,energygrid1%nLeads
        call allocate_sigma_single(sigma_buffer(il),energygrid1%leadsTotalDim(il),myhead,0,(0.0_kdp,0.0_kdp),0.0_kdp,0.0_kdp)
      enddo

      bytes_sigma=0
      do il=1,energygrid1%nLeads
        bytes_sigma=bytes_sigma + 16 *  sigma_buffer(il)%n**2
      enddo


      sigma_io=22349
      if(ik==1)then
        OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='UNKNOWN', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)
      else
        OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='OLD', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)
      endif

      rec1=energygrid1%nEnergiesGlobal*energygrid1%nspin*(ik-1)
      do ie=1,energygrid1%nEnergies
        do ispin=1,energygrid1%nspin

          rec1=rec1+1
          write(sigma_io,REC=rec1)(energygrid1%sigma(il,ie,ispin,ik)%sigma,il=1,energygrid1%nLeads)

        enddo
      enddo



#ifdef MPI
      do ihead=1,nheads-1

        do ie=1,energygrid1%nEnergies
          do ispin=1,energygrid1%nspin

            do il=1,energygrid1%nLeads
              call MPI_RECV(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex, ihead, 1, inverseheads_comm, istatus, MPIerror)
            enddo
            rec1=rec1+1
            write(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

          enddo
        enddo

        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo
#endif
      close(sigma_io)

      do il=1,energygrid1%nLeads
        deallocate(sigma_buffer(il)%sigma)
      enddo
      deallocate(sigma_buffer)


#ifdef MPI
    else


      do ihead=1,nheads-1
        if(myhead.eq.ihead)then

          do ie=1,energygrid1%nEnergies
            do ispin=1,energygrid1%nspin
              do il=1,energygrid1%nLeads
                call MPI_SEND(energygrid1%sigma(il,ie,ispin,ik)%sigma,energygrid1%sigma(il,ie,ispin,ik)%n**2,DAT_dcomplex,0,1,inverseheads_comm,MPIerror)
              enddo
            enddo
          enddo

        endif
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo


#endif
    endif

  end SUBROUTINE write_selfenergies


  subroutine check_file(fname,fileexists)
   
    use mMPI_NEGF

    character(LEN=*), intent(in) ::  fname
    logical, intent(out) ::  fileexists
    integer   iostatv,MPIerror,iostatv2

    if(myhead==0)then
      open(unit=22349,file=fname,STATUS='OLD',iostat=iostatv)
      if (iostatv .eq. 0 ) then
        fileexists=.true.
      else
        fileexists=.false.
      endif
      close(unit=22349,iostat=iostatv2)
!      write(*,*)"iostatv=",iostatv,fileexists
    endif
#ifdef MPI
    Call MPI_Bcast(fileexists,1,MPI_LOGICAL,0,inverseheads_comm,MPIerror)
#endif
    

  end subroutine check_file


  subroutine allocate_sigma2(storesigmai,LeadsVoltageShift,energygrid1)
  
    use mMPI_NEGF
  
    implicit none
    
    type(EnergyGridType), intent(inout) :: energygrid1
    integer, intent(in) :: storesigmai
    real(kdp), intent(in) :: LeadsVoltageShift(energygrid1%nLeads)
    integer ikl,ispin,ie,il,n
  
    if (allocated(energygrid1%sigma)) deallocate(energygrid1%sigma)
    allocate(energygrid1%sigma(energygrid1%nLeads,energygrid1%nEnergies,energygrid1%nspin,energygrid1%nk))

    do ikl=1,energygrid1%nk
      do ispin=1,energygrid1%nspin
        do ie=1,energygrid1%nEnergies
          do il=1,energygrid1%nLeads
            if(il==1)then
              energygrid1%sigma(il,ie,ispin,ikl)%side='L'
            elseif(il==2)then
              energygrid1%sigma(il,ie,ispin,ikl)%side='R'
            endif
 
            n=energygrid1%leadsTotalDim(il)
            energygrid1%sigma(il,ie,ispin,ikl)%n=n
            energygrid1%sigma(il,ie,ispin,ikl)%node=myhead
            energygrid1%sigma(il,ie,ispin,ikl)%InfoSigma=storesigmai
            energygrid1%sigma(il,ie,ispin,ikl)%e=energygrid1%e(ie)-LeadsVoltageShift(il)+(0.0_kdp,1.0_kdp)*energygrid1%deltasigma
  
            if(storesigmai==0) allocate(energygrid1%sigma(il,ie,ispin,ikl)%sigma(n,n))

          enddo
        enddo
      enddo
    enddo
  
  end subroutine allocate_sigma2

  subroutine write_single_selfenergy(sigma_buffer,ik,ispin,ie,energygrid1)

    use mMPI_NEGF

    integer, intent(in) :: ik,ispin,ie
    type(EnergyGridType), intent(in) :: energygrid1
    type(SelfEnergyType), intent(inout) :: sigma_buffer(energygrid1%nLeads)


    integer il,ihead,bytes_sigma,sigma_io
    integer MPIerror,i1,i2
    integer rec1
#ifdef MPI
    INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif
    character fsigma*250,char_bias*7


    write(char_bias,'(f7.4)')energygrid1%v*13.6057_kdp
    fsigma=trim(energygrid1%sLabel)//'.V_'//TRIM(ADJUSTL(char_bias))//trim(energygrid1%SigmaSuffix)

    if(myhead==0)then

      sigma_io=22349

      bytes_sigma=0
      do il=1,energygrid1%nLeads
        bytes_sigma=bytes_sigma + 16 *  sigma_buffer(il)%n**2
      enddo

      if(ik==1)then
        OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='UNKNOWN', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)
      else
        OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='OLD', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)
      endif


!***********************
      rec1=energygrid1%nEnergiesGlobal*energygrid1%nspin*(ik-1)
      rec1=rec1+(ie-1) * nheads *energygrid1%nspin+ispin

      write(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

#ifdef MPI
      do ihead=1,nheads-1

        do il=1,energygrid1%nLeads
          call MPI_RECV(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex, ihead, 1, inverseheads_comm, istatus, MPIerror)
        enddo

        rec1=rec1+energygrid1%nspin
        write(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

        CALL MPI_BARRIER(inverseheads_comm, MPIerror)

      enddo
#endif

      
      close(sigma_io)

#ifdef MPI
    else

      do ihead=1,nheads-1
        if(myhead.eq.ihead)then

          do il=1,energygrid1%nLeads
            call MPI_SEND(sigma_buffer(il)%sigma,sigma_buffer%n**2,DAT_dcomplex,0,1,inverseheads_comm,MPIerror)
          enddo

        endif
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo

#endif
    endif

  end SUBROUTINE write_single_selfenergy


  subroutine readsingle_selfenergy(sigma_buffer,ik,ispin,ie,energygrid1)

    use mMPI_NEGF

    integer, intent(in) :: ik,ispin,ie
    type(EnergyGridType), intent(inout) :: energygrid1
    type(SelfEnergyType), intent(inout) :: sigma_buffer(energygrid1%nLeads)


    integer il,ihead,bytes_sigma,sigma_io
    integer MPIerror,i1,i2
    integer rec1
#ifdef MPI
    INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif
    character fsigma*250,char_bias*7


    write(char_bias,'(f7.4)')energygrid1%v*13.6057_kdp
    fsigma=trim(energygrid1%sLabel)//'.V_'//TRIM(ADJUSTL(char_bias))//trim(energygrid1%SigmaSuffix)

    if(myhead==0)then

      sigma_io=22349

      bytes_sigma=0
      do il=1,energygrid1%nLeads
        bytes_sigma=bytes_sigma + 16 *  sigma_buffer(il)%n**2
      enddo

      OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='OLD', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)

!***********************
      rec1=energygrid1%nEnergiesGlobal*energygrid1%nspin*(ik-1)
      rec1=rec1+ie * nheads *energygrid1%nspin+ispin

#ifdef MPI
      do ihead=nheads-1,1,-1

        rec1=rec1-energygrid1%nspin
        read(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

        do il=1,energygrid1%nLeads
          call MPI_SEND(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex,ihead,1,inverseheads_comm,MPIerror)
        enddo

        CALL MPI_BARRIER(inverseheads_comm, MPIerror)

      enddo
#endif

      rec1=rec1-energygrid1%nspin
      read(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

      close(sigma_io)

#ifdef MPI
    else

      do ihead=nheads-1,1,-1
        if(myhead.eq.ihead)then

          do il=1,energygrid1%nLeads
            call MPI_RECV(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex, 0, 1, inverseheads_comm, istatus, MPIerror)
          enddo

        endif
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo

#endif
    endif

  end SUBROUTINE readsingle_selfenergy


  subroutine single_selfenergy(sigmai,n,il,ik,ispin,ie,energygrid1)

    use negfmod, only : emtimings
    use sigma, only: h0_l,h1_l,s0_l,s1_l,h0_r,h1_r,s0_r,s1_r
    use mSelfenergies, only : SelfEnergyGeneral

    integer, intent(in) :: il,ik,ispin,ie,n
    complex(kdp),intent(out) :: sigmai(n,n)
    type(EnergyGridType), intent(in) :: energygrid1

    integer nchan
    integer*4:: sc_0,sc_1,sc_r,sc_m

    if(emtimings)CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)

    if(energygrid1%sigma(il,ie,ispin,ik)%Side=='L')then
      call SelfEnergyGeneral(energygrid1%sigma(il,ie,ispin,ik)%Side,n ,energygrid1%sigma(il,ie,ispin,ik)%e,H0_L(:,:,ispin), H1_L(:,:,ispin),S0_L(:,:),S1_L(:,:),sigmai,nchan,0.0_kdp,.true.)
    else
      call SelfEnergyGeneral(energygrid1%sigma(il,ie,ispin,ik)%Side,n ,energygrid1%sigma(il,ie,ispin,ik)%e,H0_R(:,:,ispin), H1_R(:,:,ispin),S0_R(:,:),S1_R(:,:),sigmai,nchan,0.0_kdp,.true.)
    endif

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'tsigma1',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif

  end SUBROUTINE single_selfenergy


  subroutine all_selfenergies(ik,istart,iend,ispinstart,ispinend,energygrid1,ispinHin)

    use negfmod, only : emtimings
    use sigma, only: h0_l,h1_l,s0_l,s1_l,h0_r,h1_r,s0_r,s1_r
    use mSelfenergies, only : SelfEnergyGeneral

    integer, intent(in) :: ik,istart,iend,ispinstart,ispinend,ispinHin
    type(EnergyGridType), intent(inout) :: energygrid1

    integer ispin,ie,nchan,ispinh
    integer*4:: sc_0,sc_1,sc_r,sc_m
    logical :: DoFourier 

    do ispin=ispinstart,ispinend
      if(ispinHin==0)then
        ispinh=ispin
      else
        ispinh=ispinHin
      endif
      DoFourier=.true.
      do ie=istart,iend
        if(emtimings)CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)

        call SelfEnergyGeneral('L',energygrid1%sigma(1,ie,ispin,ik)%n ,energygrid1%sigma(1,ie,ispin,ik)%e,H0_L(:,:,ispinh), H1_L(:,:,ispinh),S0_L(:,:),S1_L(:,:),energygrid1%sigma(1,ie,ispin,ik)%sigma ,nchan,0.0_kdp,DoFourier)
        DoFourier=.false.
        
        if(ispinHin.ne.0)energygrid1%sigma(1,ie,ispin,ik)%nchannels=nchan

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')'tsigma2',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif
      enddo

      DoFourier=.true.
      do ie=istart,iend
        if(emtimings)CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)

        call SelfEnergyGeneral('R',energygrid1%sigma(2,ie,ispin,ik)%n ,energygrid1%sigma(2,ie,ispin,ik)%e,H0_R(:,:,ispinh), H1_R(:,:,ispinh),S0_R(:,:),S1_R(:,:),energygrid1%sigma(2,ie,ispin,ik)%sigma ,nchan,0.0_kdp,DoFourier)
        DoFourier=.false.
        
        if(ispinHin.ne.0)energygrid1%sigma(2,ie,ispin,ik)%nchannels=nchan

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')'tsigma3',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif
      enddo

    enddo

  end SUBROUTINE all_selfenergies

   
  subroutine energygrid_adapt2(N1,NSPIN, energi,energf,Delta, V,Nenerg_div,Nenerg_div_nodes, ik,hgeneral,sgeneral,nk,storesigmai,nleadslr,LeadsVoltageShift)


! **********************************************************************
! Written by Ivan Rungger, Chaitanya Das Pemmaraju and 
! Alexandre Reily Rocha, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *********************************************************************


      use mMPI_NEGF

      use sigma
      use negfmod
      use global_meshvar
      use mTypes
           
      implicit none 
      
      integer, parameter :: maxsize=100 !maximum number of refinement steps, replace with input option, replace with maxdepth
      
      integer, intent(in):: nleadslr,storesigmai
      real(kdp), intent(in) :: LeadsVoltageShift(nleadslr)
      integer :: N1,NL,NR,NSPIN,Nenerg_div, nenerg_div_nodes,I,ik,N_IVPOINTS, idepth,isize,Nenerg_new,Nenerg_total, Nenerg_div_start,Nenerg_div_end,ii,jj,ind,nk
      
      integer, dimension (nheads+1,maxsize) :: Nenerg_pernode
      
      double precision, dimension (2*maxsize) :: Energyranges

      double precision :: EnergI,EnergF,V,Delta,dE,deltaout
      double precision de_ini_used, deltaini_used
      type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
      integer iglobal
      
#ifdef MPI
      INTEGER :: MPIerror
#endif
 
! **************************************************************
! If the self-energies have already been allocated, we deallocate
! them before-hand

      de_ini_used=ABS(EnergI-EnergF)/(1D0 * nenerg_div_nodes)
      deltaini_used=de_ini_used * deltatode
!      if(myhead.eq.0) then
!        write(*,*)"nenerg_div=",nenerg_div,nenerg_div_nodes, de_ini_used,deltaini_used,deltaini,deltatode
!        write(*,*)"energy range=",setene,energi,energf
!      endif


! initial grid----------------------------
      N_IVPOINTS=0

      if (allocated(EXG)) deallocate (EXG,CONSTG)
      allocate(EXG(Nenerg_div_nodes),CONSTG(Nenerg_div_nodes))
      CONSTG=0.d0

      dE=(EnergF-EnergI)/(Nenerg_div_nodes-1)
      DO I=1,Nenerg_div_nodes
        EXG(I)=EnergI+(EnergF-EnergI)*(I-1)/(Nenerg_div_nodes-1)
      ENDDO

      NEnerg_pernode=0
      NEnerg_pernode(1:nheads,1)=NEnerg_div
      Nenerg_pernode(nheads+1,1)=1

      DO i=1,Nenerg_div
        iglobal=myhead * Nenerg_div+i
        ERealGrid%e(i)=EXG(iglobal)
        ERealGrid%w(i)=constg(iglobal)
        ERealGrid%ig(i)=iglobal
      ENDDO
! end setting up initial grid
      
      call allocate_sigma2(storesigmai,LeadsVoltageShift,ERealGrid)

      allocate(mydepth(ERealGrid%nEnergies)) 
      mydepth(:)=1
  
!---------initialize variables--------------------------
      idepth=1
      isize=1
      Nenerg_div_start=1
      Nenerg_div_end=Nenerg_div
      Nenerg_new=Nenerg_div_nodes
      Nenerg_total=Nenerg_div_nodes
      Delta=deltaini_used

      Energyranges(:)=0.0d0
      Energyranges(1)=EnergI
      Energyranges(2)=EnergF
       

      nl=ERealGrid%leadsTotalDim(1)
      nr=ERealGrid%leadsTotalDim(2)
      write(12347,*)"eneindex0=",Nenerg_div_start,Nenerg_div_end,nl,nr,n1
!---------find energy mesh iteratively--------------------------
      call recursive_energygrid(N1,NL,NR,NSPIN, IDepth,isize,Delta,V,dE,Nenerg_div_start,  Nenerg_div_end,Nenerg_new,Nenerg_total, ik,deltaimag, Nenerg_pernode,Energyranges,critam,hgeneral,sgeneral,storesigmai,LeadsVoltageShift,nleadslr)


      Nenerg_div_nodes=Nenerg_total
      Nenerg_div=SUM(Nenerg_pernode(myhead+1,1:isize))
      deltaout=delta
 
      if(delta/1000D0.gt.deltamin)then
        Delta=Delta/1000D0
      else
        Delta=deltamin
      endif
      if(myhead .eq. 0) write(6,*) "adaptive:  Final parameters", " ik=",ik,Nenerg_div_nodes,Nenerg_div,deltaout,delta,ERealGrid%nEnergies,ERealGrid%nEnergiesGlobal

      call calculate_weights(Nenerg_div_nodes,  maxsize,Energyranges)     


      call redistribute_const(maxsize, isize,Nenerg_pernode)         

      deallocate(mydepth)

         

      end subroutine energygrid_adapt2

! *************************************************************************

  subroutine recursive_energygrid(N1,NL,NR,NSPIN, IDepth,isize,Delta,V,dE,Nenerg_div_start, Nenerg_div_end,Nenerg_new,Nenerg_total, ik,deltaimag, Nenerg_pernode,Energyranges,criteria,hgeneral,sgeneral,storesigmai,LeadsVoltageShift,nleadslr)

! *************************************************************************
! subroutine to calculate the adaptive mesh iteratively
!
! Written by Ivan Rungger, Chaitanya Das Pemmaraju and 
! Alexandre Reily Rocha, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *********************************************************************


      use mMPI_NEGF

      use sigma
      use global_meshvar
      use negfmod,only:ikpmod,deltamin,em_iscf,negfon,maxdepth, ndivisions,inversion_solver
      use mTypes
      use mMatrixUtil
      use mONInterface
            
      implicit none

      include "const2.h"
      
      integer, intent(in) :: storesigmai
      integer, intent(in):: nleadslr
      real(kdp), intent(in) :: LeadsVoltageShift(nleadslr)

      integer, parameter :: maxsize=100
      
      double precision :: criteria

      integer :: N1,MPIerror,idepth,Nenerg_div_start, Nenerg_div_end,I,ISPIN,NSPIN,NL,NR,II,ik, INFO,JJ,Nenerg_new,Nenerg_total,NPeaks, Idepth_aux,isize, lwork,Nenerg_old

 
      integer, dimension (:), allocatable :: IPIV 
      
      integer, dimension (nheads+1,maxsize) :: Nenerg_pernode
      
      integer, dimension (:), allocatable :: INDEX
      
      integer, dimension (:), allocatable :: INDEXG

      double precision :: Delta,Delta_aux,dE,dE_new,V,Ei,deltaimag

      double precision, dimension (2*maxsize) :: Energyranges
           
      double complex gammal(nl,nl),gammar(nr,nr)

      double complex, allocatable, dimension (:) :: WORK
      
      double complex  totaldos(NSPIN)
      double precision  totaldos2(nspin)
      double precision  totalpdos(n1)
      type(ioType) :: io

!--------------------
      INTEGER::NENERG_LOC,tot_nploc, tot_nploc_new,irec
      double complex, allocatable::TEMPEXG(:)
      type(matrixTypeGeneral) :: gfmat,gfout
      type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
      double complex gf_iter1l(n1,nl),gf_iter2l(n1,nl), gf_iter1r(n1,nr),gf_iter2r(n1,nr)
      integer gfmattype,nnz,i2,i1,ind,ind2,ind3
      integer  nnzrow(n1)
      double complex gfadd
      double complex gfu
      double complex gf_diag(n1)
      double complex, allocatable :: EXG_local(:)

!--------------------           
      lwork=10 * N1

      if(.not.negfon)then
        gfmattype=0 !for dense GF matrix
      else
        gfmattype=2 !for sparse GF matrix
      endif

      do irec=1,maxdepth+1


      NENERG_LOC=Nenerg_div_end-Nenerg_div_start+1

      if(myhead .eq. 0) write(*,*)"adaptive:",delta,de,idepth,maxdepth, Nenerg_div_start,Nenerg_div_end, Nenerg_div_end-Nenerg_div_start

!      write(*,*)"getting self-energies"
      call all_selfenergies(1,Nenerg_div_start,Nenerg_div_end,1,ERealGrid%nspin,ERealGrid,0)
!      write(*,*)"got self-energies"


      if(allocated(INDEXG)) then
        deallocate(INDEXG,TEMPEXG)
      endif

      allocate(INDEXG(Nenerg_new))
      ALLOCATE(TEMPEXG(Nenerg_new))


      indexG=0
      if (Nenerg_div_start.le.Nenerg_div_end) then
        allocate(IPIV(N1),WORK(lwork))
        allocate(index(Nenerg_div_end-Nenerg_div_start+1))
        index=0
      endif
        
      if(idepth.gt.maxdepth.or.(delta.le.deltamin)) then
!       write(6,*) "adaptive: reached maximum depth - returning",
!     &    idepth,maxdepth,delta,deltamin
        if(allocated(IPIV))then
          deallocate(IPIV,WORK)
        endif            
        exit 
      endif
      
      
!---Calculate the total DOS to identify peaks-------------------------!      
      DO I=Nenerg_div_start,Nenerg_div_end

        Ei=ERealGrid%e(i)
        
        totaldos=0.0D0

        DO ISPIN=1,NSPIN

          if(.not.negfon)then
            nnzrow=0
            nnz=n1*n1
          else
           call findnnzgf2(nnz,nnzrow,n1,n1,nl,nr, ERealGrid%sigma(1,i,ispin,1)%sigma, ERealGrid%sigma(2,i,ispin,1)%sigma, hgeneral(ispin)) 
          endif

          call AllocateMatrixGeneral(n1,n1,nnz,gfmattype,gfmat, "adaptivegrid", io)
          
         write(12347,*)"eneindex=",i,ei,Nenerg_div_start,Nenerg_div_end,nl,nr,n1

         call setgfelementsgeneral_nc(Ei+zi*Delta,nspin,ispin,gfmat,nnz,n1,nl,nr, ERealGrid%sigma(1,i,ispin,1)%sigma, ERealGrid%sigma(2,i,ispin,1)%sigma,hgeneral,sgeneral)

          if(.not.negfon)then

!            write(*,*)"calculating inverse"

            CALL ZGETRF(N1,N1,gfmat%matdense%a,N1,IPIV,INFO)
            CALL ZGETRI(N1,gfmat%matdense%a,N1,IPIV,WORK,lwork,INFO)
            gf_iter1l=gfmat%matdense%a(:,1:nl)
            gf_iter1r=gfmat%matdense%a(:,n1-nr+1:n1)

!            write(*,*)"done inverting",i,myhead

          else

!            write(*,*)"calculating inverseon"

            call AllocateMatrixGeneral(n1,nl+nr,n1*(nl+nr),0,gfout, "recursive_energygrid", io)

!            call invertdiagonalandcolumnsONGeneral(N1,gfmat,nl,nr,gfout)
            call InvertONGeneral(N1,gfmat,nl,nr,gfout,3,inversion_solver)
            gf_iter1l=gfout%matdense%a(:,1:nl)
            gf_iter1r=gfout%matdense%a(:,nl+1:nl+nr)

            call DestroyMatrixGeneral(gfout,"adaptivegrid",io)

            write(*,*)"done invertingon",i,myhead

          endif

          gammal=ERealGrid%sigma(1,i,ispin,1)%sigma
          gammal=zi*(gammal-DCONJG(TRANSPOSE(gammal)))
          gammar=ERealGrid%sigma(2,i,ispin,1)%sigma
          gammar=zi*(gammar-DCONJG(TRANSPOSE(gammar)))

          do ii=1,n1
            gf_diag(ii)=0D0
            do ind=sgeneral%matSparse%q(ii), sgeneral%matSparse%q(ii+1)-1
              if(gfmattype.eq.0)then
                gfu=zi * (gfmat%matdense%a(sgeneral%matSparse%j(ind),ii)- DCONJG(gfmat%matdense%a(ii,sgeneral%matSparse%j(ind))))
              else

                call inddensetoindsparsegeneral(sgeneral%matSparse%j(ind),ii,ind2,gfmat)
                call inddensetoindsparsegeneral(ii,sgeneral%matSparse%j(ind),ind3,gfmat)
                if(ind2.ne.0.and.ind3.ne.0)then
                  gfu=zi * (gfmat%matSparse%b(ind2)- DCONJG(gfmat%matSparse%b(ind3)))
                else
                  gfu=0D0
                endif

              endif
              gf_diag(ii)=gf_diag(ii)+sgeneral%matSparse%b(ind)*gfu
            enddo
          enddo


          CALL ZGEMM('N','N',N1,NL,NL,(1.D0,0.D0),GF_iter1l,N1, gammal,NL,(0.D0,0.D0),GF_iter2l,N1)
          
          do ii=1,n1
            do ind=sgeneral%matSparse%q(ii), sgeneral%matSparse%q(ii+1)-1
              jj=sgeneral%matSparse%j(ind)

              gfadd=0D0
              do i1=1,nl
                gfadd=gfadd+gf_iter2l(ii,i1)* DCONJG(gf_iter1l(jj,i1))
              enddo

              gf_diag(ii)=gf_diag(ii)-sgeneral%matSparse%b(ind)*gfadd
             
            enddo
          enddo


          CALL ZGEMM('N','N',N1,NR,NR,(1.D0,0.D0),GF_iter1r,N1, gammar,NR,(0.D0,0.D0),GF_iter2r,N1)
          
          do ii=1,n1
            do ind=sgeneral%matSparse%q(ii), sgeneral%matSparse%q(ii+1)-1
              jj=sgeneral%matSparse%j(ind)

              gfadd=0D0
              do i1=1,nr
                gfadd=gfadd+gf_iter2r(ii,i1)* DCONJG(gf_iter1r(jj,i1))
              enddo

              gf_diag(ii)=gf_diag(ii)-sgeneral%matSparse%b(ind)*gfadd
             
            enddo
          enddo

          Do II=1,N1
            if(orbital_BS(II)) then
              totaldos(ISPIN)=totaldos(ISPIN)+gf_diag(ii)
            endif
          enddo

          call DestroyMatrixGeneral(gfmat,"adaptivegrid",io)

          if (DREAL(totaldos(ISPIN)).gt.(1.0d0/(criteria*Delta))) then
            INDEX(I-Nenerg_div_start+1)=1
          endif
        ENDDO

        if(NSPIN.eq.2) write(12346,*)"tdos=",ei,DREAL(totaldos(1)),DREAL(totaldos(2)),DIMAG(totaldos(2)),INDEX(I-Nenerg_div_start+1)/(criteria * Delta),'iscf=',em_iscf,'ik=',ik 
        if(NSPIN.eq.1)  write(12346,*)"tdos=",ei,DREAL(totaldos(1)), DIMAG(totaldos(1)), INDEX(I-Nenerg_div_start+1)/(criteria * Delta), 'iscf=',em_iscf,'ik=',ik 
      ENDDO            

      if (allocated(IPIV))then
        deallocate(IPIV,WORK)
      endif            


#ifdef MPI
      CALL MPI_GATHER(INDEX(1),NENERG_LOC, MPI_INTEGER,INDEXG(1),NENERG_LOC,MPI_INTEGER,0, inverseheads_comm,MPIerror)
      CALL MPI_GATHER(ERealGrid%e(Nenerg_div_start),NENERG_LOC, DAT_dcomplex,TEMPEXG(1),NENERG_LOC, DAT_dcomplex,0,inverseheads_comm,MPIerror)

      call mpi_bcast(INDEXG,Nenerg_new,MPI_INTEGER,0, inverseheads_comm,MPIError)
      call mpi_bcast(TEMPEXG(1),Nenerg_new,DAT_dcomplex,0,inverseheads_comm,MPIError)
#else
      INDEXG=index
      TEMPEXG(1:NENERG_LOC)=ERealGrid%e(Nenerg_div_start:Nenerg_div_start+NENERG_LOC-1)
#endif

      IF(ALLOCATED(INDEX))DEALLOCATE(INDEX)

      if (maxval(INDEXG).eq.0) then
        exit
      endif

      Nenerg_old=Nenerg_new
 


      call add_energypoints(Nenerg_old,indexg,tempexg,npeaks, ndivisions,nheads,Nenerg_new,nenerg_total,isize,maxsize, Nenerg_pernode,NENERG_LOC,tot_nploc,tot_nploc_new,myhead,de_new)

      idepth=idepth+1

      allocate(EXG_local(NENERG_LOC))
#ifdef MPI
      CALL MPI_SCATTER(EXG(Nenerg_total+1),NENERG_LOC,DAT_dcomplex,EXG_local(1),NENERG_LOC,DAT_dcomplex,0,inverseheads_comm,MPIerror)
#else
      EXG_local(1:NENERG_LOC)=exg(Nenerg_total+1:Nenerg_total+NENERG_LOC)
#endif
      call reallocate_energygrid(nspin,tot_nploc,tot_nploc_new,nenerg_total,nenerg_loc,EXG_local,v,idepth,ERealGrid,storesigmai,LeadsVoltageShift)
      deallocate(EXG_local)

      call realloc_mydepth(NL,NR,nspin,tot_nploc, tot_nploc_new,nenerg_total,nenerg_loc,v,idepth)


      Nenerg_div_start=tot_nploc+1
      Nenerg_div_end=tot_nploc_new
      Delta=Delta/(dE/dE_New)
      Nenerg_total=Nenerg_new+Nenerg_total        
      de=dE_New

      Energyranges(2*isize-1)=DREAL(EXG(Nenerg_total-Nenerg_new+1))
      Energyranges(2*isize)  =DREAL(EXG(Nenerg_total))

      idepth_aux=idepth
      Delta_aux=delta


      end do


      deallocate(indexG)
      DEALLOCATE(TEMPEXG)
      RETURN
      
  end subroutine recursive_energygrid



  SUBROUTINE realloc_global(Nenerg_total,Nenerg_new_tot)
       
! *************************************************************************
! Written by Chaitanya Das Pemmaraju, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: pemmaras@tcd.ie
! *********************************************************************

       use global_meshvar
       
       implicit none
       
       integer :: Nenerg_total,Nenerg_new_tot
       
       double complex, allocatable, dimension (:) :: aux
       
       if (Nenerg_new_tot.gt.Nenerg_total) then
        if (allocated(CONSTG)) then
   allocate(aux(Nenerg_total))
         
   aux=CONSTG
   deallocate(CONSTG)
   allocate(CONSTG(Nenerg_new_tot))
   CONSTG=(0.d0,0.d0)
   CONSTG(1:Nenerg_total)=aux
   
   aux=EXG
   deallocate(EXG)
   allocate(EXG(Nenerg_new_tot))
   EXG=(0.d0,0.d0)
   EXG(1:Nenerg_total)=aux
  else
   allocate(CONSTG(Nenerg_new_tot))
   allocate(EXG(Nenerg_new_tot))
   CONSTG=0.d0
   EXG=0.d0
  endif
       endif
       
  end subroutine realloc_global


  subroutine calculate_weights(Nenerg_div_nodes, maxsize,Energyranges)

! *************************************************************************
! Written by Chaitanya Das Pemmaraju, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: pemmaras@tcd.ie
! *********************************************************************

       use global_meshvar  
  
  implicit none 
  
  integer :: Nenerg_div_nodes, maxsize,I,jj
     
  integer, dimension (Nenerg_div_nodes) :: INDEX_aux
  
  double precision :: aux1,hl,hr
  
  double precision, dimension (Nenerg_div_nodes) :: EXG_aux,CONST_aux
  
  double precision, dimension (2*maxsize) :: Energyranges
  
  integer, dimension (2*maxsize) :: IZ
  
  
  EXG_aux=DREAL(EXG)
  DO I=1,Nenerg_div_nodes
   INDEX_aux(I)=I
  ENDDO
  
  caLL SSORT (EXG_aux, INDEX_aux, Nenerg_div_nodes, 2)
  
  IZ=0
  CALL SSORT (Energyranges, IZ, maxsize*2, 1)
  
  CONST_aux=0.d0
  aux1=0.0d0
  if(Nenerg_div_nodes .eq. 1) then
    CONST_aux=1.0d0
    GoTo 111
  endif

  if(mod(Nenerg_div_nodes,2) .eq. 0) then
    DO I=2, Nenerg_div_nodes-2, 2
      hl=DABS(EXG_aux(I)-EXG_aux(I-1))
      hr=DABS(EXG_aux(I+1)-EXG_aux(I))
      if(dabs(hr-hl) .lt. 1.0D-12) then
        CONST_aux(I-1)=CONST_aux(I-1)+hl/3.0d0
        CONST_aux(I)=CONST_aux(I)+(4.0d0*hl/3.0d0)
        CONST_aux(I+1)=CONST_aux(I+1)+hl/3.0d0
      else
        CONST_aux(I-1)=CONST_aux(I-1)+hl/2.0d0
        CONST_aux(I)=CONST_aux(I)+hl/2.0d0+hr/2.0d0
        CONST_aux(I+1)=CONST_aux(I+1)+hr/2.0d0
      endif
      
    ENDDO
    hl=DABS(EXG_aux(Nenerg_div_nodes)- EXG_aux(Nenerg_div_nodes-1))
    CONST_aux(Nenerg_div_nodes-1)=CONST_aux(Nenerg_div_nodes-1)+hl/2.0d0
    CONST_aux(Nenerg_div_nodes)= CONST_aux(Nenerg_div_nodes)+hl/2.0d0
  else 
    DO I=2, Nenerg_div_nodes-1, 2
      hl=DABS(EXG_aux(I)-EXG_aux(I-1))
      hr=DABS(EXG_aux(I+1)-EXG_aux(I))
      if(dabs(hr-hl) .lt. 1.0D-12) then
        CONST_aux(I-1)=CONST_aux(I-1)+hl/3.0d0
        CONST_aux(I)=CONST_aux(I)+(4.0d0*hl/3.0d0)
        CONST_aux(I+1)=CONST_aux(I+1)+hl/3.0d0
      else
        CONST_aux(I-1)=CONST_aux(I-1)+hl/2.0d0
        CONST_aux(I)=CONST_aux(I)+hl/2.0d0+hr/2.0d0
        CONST_aux(I+1)=CONST_aux(I+1)+hr/2.0d0
      endif
    ENDDO
  endif
  
 111  continue

  aux1=0.0d0
  DO JJ=1,Nenerg_div_nodes
   CONSTG(INDEX_aux(JJ))=DCMPLX(CONST_aux(JJ))
   aux1=aux1+CONSTG(INDEX_aux(JJ))
  ENDDO

!$$$  write(6,*) "Sum of Weights =", aux1
!$$$  write(6,*) "Width of interval =",
!$$$     $             EXG_aux(Nenerg_div_nodes)-EXG_aux(1) 
  

  end subroutine calculate_weights
! ****************************************************************************************
  SUBROUTINE redistribute_const(maxsize, Nsize,Nenerg_pernode)

! *************************************************************************
! Written by Chaitanya Das Pemmaraju, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: pemmaras@tcd.ie
! *********************************************************************


      use mMPI_NEGF

      use sigma
      use global_meshvar

      implicit none 

      integer :: nsize,maxsize,I, MPIerror,gptr,lptr,npts
      integer, dimension (nheads+1,maxsize) :: Nenerg_pernode

  
      gptr=1
      lptr=1
    
      do i=1, nsize
    
        npts=Nenerg_pernode(myhead+1,i)
    
#ifdef MPI
        CALL MPI_SCATTER(CONSTG(gptr), npts,DAT_dcomplex, ERealGrid%w(lptr), npts,DAT_dcomplex,0, inverseheads_comm,MPIerror)          
#else
        ERealGrid%w(lptr:lptr+npts-1)=CONSTG(gptr:gptr+npts-1)  
#endif
    
        gptr=gptr+SUM(Nenerg_pernode(1:nheads,i))
        lptr=lptr+npts
    
      enddo


       
  end subroutine redistribute_const


  subroutine add_energypoints(ne_old,indexg,tempexg,npeaks, ndivisions,nnodes,ne_new,nenerg_total,isize,maxsize, Nenerg_pernode,NENERG_LOC,tot_nploc,tot_nploc_new,mynode,de_new)

! *************************************************************************
! Written by Chaitanya Das Pemmaraju, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: pemmaras@tcd.ie
! *********************************************************************

      use global_meshvar
      implicit none

      integer npeaks,ne_old,indexg(ne_old),i,ipeak,tot_npts, tot_npts_new,maxsize,tot_nploc,tot_nploc_new, ndivisions,nnodes,ne_new,ii,nenerg_total,isize,Nenerg_pernode(Nnodes+1,maxsize),NENERG_LOC, mynode
      double complex tempexg(ne_old)
      integer,allocatable::Npts(:),ilft(:),irgt(:)
      double precision, allocatable::Elft(:),Ergt(:)
      double precision::Ewindow,de_new,ei

      NPEAKS=0
      if (INDEXG(1).eq.1) then
        NPeaks=1
      else
        NPeaks=0
      endif

      DO I=2,ne_old
        if ((INDEXG(I).EQ.1).AND.(INDEXG(I-1).EQ.0)) then
          NPeaks=NPeaks+1
        endif
      ENDDO  


      allocate(elft(npeaks))
      allocate(ergt(npeaks))
      allocate(ilft(npeaks))
      allocate(irgt(npeaks))
      allocate(npts(npeaks))

      IF(indexG(1) .eq. 1) then
        ipeak=1
        elft(1)=tempexg(1)
        ilft(1)=1
      else
        ipeak=0
      endif

      DO I=2,ne_old
        if ((INDEXG(I).EQ.1).AND.(INDEXG(I-1).EQ.0)) then
          iPeak=iPeak+1
          elft(ipeak)=tempexg(I-1)
          ilft(ipeak)=I-1
        endif
        if((INDEXG(I).EQ.0) .AND. (INDEXG(I-1).EQ.1)) then
          ergt(ipeak)=tempexg(I)
          irgt(ipeak)=I
        endif
      ENDDO  

      if(indexG(ne_old) .eq.1) then
        ergt(ipeak)=tempexg(ne_old)
        irgt(ipeak)=ne_old
      endif  

      tot_npts=0
      Ewindow=0.0d0
      do ipeak=1, npeaks
        npts(ipeak)=irgt(ipeak)-ilft(ipeak)+1
        tot_npts=tot_npts+npts(ipeak)
        Ewindow=Ewindow+ergt(ipeak)-elft(ipeak)
      enddo
      tot_npts_new=tot_npts*ndivisions
      if(MOD(tot_npts_new,Nnodes) .NE. 0) tot_npts_new=(tot_npts_new/Nnodes+1)*Nnodes
      dE_New=Ewindow/(tot_npts_new+1)

      ne_new=tot_npts_new

      call realloc_global(Nenerg_total,ne_new+Nenerg_total)

      II=1
      DO ipeak=1, Npeaks
        II=II-1
        Ei=elft(ipeak)+dE_New/Sqrt(17.0d0)
        DO
          II=II+1
          Ei=Ei+dE_New
          IF(Ei .ge. ergt(ipeak) .or. II.gt.ne_new) exit
          EXG(II+Nenerg_total)=Ei
        ENDDO
      ENDDO

      DO WHILE(II .LE. ne_new) 
        EXG(II+Nenerg_total)=Ei
        II=II+1
        Ei=Ei+dE_New
      ENDDO

      deallocate(elft,ergt,ilft,irgt,npts)

      isize=isize+1 !replece with idepth?
      Nenerg_pernode(1:Nnodes,isize)=ne_new/Nnodes
      NENERG_LOC=ne_new/Nnodes

      tot_nploc=SUM(Nenerg_pernode(Mynode+1,1:isize-1))
      tot_nploc_new=SUM(Nenerg_pernode(Mynode+1,1:isize))


  end subroutine add_energypoints
  

  subroutine realloc_mydepth(NL,NR,nspin,tot_nploc,tot_nploc_new,nenerg_total,nenerg_loc,v,idepth)
       
      use global_meshvar
      use sigma,only: mydepth
#ifdef MPI
      use mMPI_NEGF
#endif

      implicit none
     

      integer NL,NR,tot_nploc,nspin,tot_nploc_new,nenerg_total, nenerg_loc,MPIerror,idepth
      integer,allocatable::dpth_tmp(:)
      double precision V

      allocate(dpth_tmp(tot_nploc))
      dpth_tmp=mydepth
      deallocate(mydepth)
      allocate(mydepth(tot_nploc_new))
      mydepth(1:tot_nploc)=dpth_tmp(1:tot_nploc)
      mydepth(tot_nploc+1:tot_nploc_new)=idepth
      deallocate(dpth_tmp)
        

  end subroutine realloc_mydepth


  subroutine duplicate_energygrid(energygrid1,energygrid2,storesigmai,LeadsVoltageShift)

    use mMPI_NEGF
    
    type(EnergyGridType), intent(in) :: energygrid1
    type(EnergyGridType), intent(out) :: energygrid2
    integer, intent(in) :: storesigmai
    real(kdp), intent(in) :: LeadsVoltageShift(energygrid1%nLeads)

    integer ie,ispin,ik,il

    energygrid2%nEnergies=energygrid1%nEnergies
    energygrid2%nEnergiesGlobal=energygrid1%nEnergiesGlobal
    energygrid2%nSpin=energygrid1%nSpin
    energygrid2%nk=energygrid1%nk
    energygrid2%v=energygrid1%v
    energygrid2%nLeads=energygrid1%nLeads
    energygrid2%GridType=energygrid1%GridType

    allocate(energygrid2%e(energygrid2%nEnergies))
    allocate(energygrid2%w(energygrid2%nEnergies))
    allocate(energygrid2%ig(energygrid2%nEnergies))
    energygrid2%e=energygrid1%e
    energygrid2%w=energygrid1%w
    energygrid2%ig=energygrid1%ig

    allocate(energygrid2%leadsTotalDim(energygrid2%nLeads))
    energygrid2%leadsTotalDim=energygrid1%leadsTotalDim

    call allocate_sigma2(storesigmai,LeadsVoltageShift,energygrid2)

    do ik=1,energygrid2%nk
      do ispin=1,energygrid2%nspin
        do ie=1,energygrid2%nEnergies
          do il=1,energygrid2%nLeads
            energygrid2%sigma(il,ie,ispin,ik)=energygrid1%sigma(il,ie,ispin,ik)
          enddo
        enddo
      enddo
    enddo

  end subroutine duplicate_energygrid



  subroutine reallocate_energygrid(nspin,nold,nnew,ntotal,nlocal,exg_local,v,idepth,energygrid1,storesigmai,LeadsVoltageShift)

    use negfmod
    use mMPI_NEGF
    
    type(EnergyGridType),intent(inout) :: energygrid1
    integer, intent(in) :: nspin,nold,nnew,ntotal,nlocal,idepth
    real(kdp), intent(in) :: v
    integer, intent(in) :: storesigmai
    real(kdp), intent(in) :: LeadsVoltageShift(energygrid1%nLeads)
    complex(kdp), intent(in) :: exg_local(nlocal)

    type(EnergyGridType) :: energygrid_buffer
    integer ie,ispin,ik,il


    call duplicate_energygrid(energygrid1,energygrid_buffer,storesigmai,LeadsVoltageShift)

    deallocate(energygrid1%e)
    deallocate(energygrid1%w)
    deallocate(energygrid1%ig)

    do ik=1,energygrid1%nk
      do ispin=1,energygrid1%nspin
        do ie=1,energygrid1%nEnergies
          do il=1,energygrid1%nLeads
            deallocate(energygrid1%sigma(il,ie,ispin,ik)%sigma)
          enddo
        enddo
      enddo
    enddo
    deallocate(energygrid1%sigma)



    energygrid1%nEnergies=nnew
    energygrid1%nEnergiesGlobal=ntotal

    allocate(energygrid1%e(energygrid1%nEnergies))
    allocate(energygrid1%w(energygrid1%nEnergies))
    allocate(energygrid1%ig(energygrid1%nEnergies))
    energygrid1%e(1:nold)=energygrid_buffer%e
    energygrid1%w(1:nold)=energygrid_buffer%w
    energygrid1%ig(1:nold)=energygrid_buffer%ig
    energygrid1%e(nold+1:nnew)=exg_local(:)

    call allocate_sigma2(storesigmai,LeadsVoltageShift,energygrid1)

    do ik=1,energygrid1%nk
      do ispin=1,energygrid1%nspin
        do ie=1,nold
          do il=1,energygrid1%nLeads
            energygrid1%sigma(il,ie,ispin,ik)%sigma=energygrid_buffer%sigma(il,ie,ispin,ik)%sigma
          enddo
        enddo
      enddo
    enddo

    deallocate(energygrid_buffer%e)
    deallocate(energygrid_buffer%w)
    deallocate(energygrid_buffer%ig)

    do ik=1,energygrid_buffer%nk
      do ispin=1,energygrid_buffer%nspin
        do ie=1,energygrid_buffer%nEnergies
          do il=1,energygrid_buffer%nLeads
            deallocate(energygrid_buffer%sigma(il,ie,ispin,ik)%sigma)
          enddo
        enddo
      enddo
    enddo
    deallocate(energygrid_buffer%sigma)

  end subroutine reallocate_energygrid


  subroutine deallocate_selfenergies(ie,ispin,ik,energygrid1)

    integer, intent(in) :: ie,ispin,ik
    type(EnergyGridType), intent(inout) :: energygrid1

    integer il

    do il=1,energygrid1%nLeads
      deallocate(energygrid1%sigma(il,1,1,1)%sigma)
    enddo
    deallocate(energygrid1%sigma)

  end subroutine deallocate_selfenergies


  subroutine get_selfenergies(ie,ispin,ik,energygrid1,v,storesigmai,LeadsVoltageShift,sigmatodisk)

    use mMPI_NEGF

    logical, intent(in) :: sigmatodisk
    integer, intent(in) :: ie,ispin,ik,storesigmai
    real(kdp), intent(in) :: v
    type(EnergyGridType), intent(inout) :: energygrid1
    real(kdp), intent(in) ::  LeadsVoltageShift(energygrid1%nLeads)

    character fsigma*250,char_bias*7
    logical,save ::  fsigma_exists
    integer il

    allocate(energygrid1%sigma(energygrid1%nLeads,1,1,1))
    do il=1,energygrid1%nLeads
      call allocate_sigma_single(energygrid1%sigma(il,1,1,1),energygrid1%leadsTotalDim(il),myhead,storesigmai,energygrid1%e(ie),LeadsVoltageShift(il),energygrid1%deltasigma)
    enddo

    write(char_bias,'(f7.4)')v*13.6057_kdp
    fsigma=trim(energygrid1%sLabel)//'.V_'//TRIM(ADJUSTL(char_bias))//trim(energygrid1%SigmaSuffix)

    if(ik==1.and.ie==1.and.ispin==1)then
      call check_file(fsigma,fsigma_exists)
    endif

    if(.not.fsigma_exists)then
      call all_selfenergies(1,1,1,1,1,energygrid1,ispin)
      if(sigmatodisk)call write_selfenergies_transmission(ik,ispin,ie,fsigma,ETransmGrid)
    else
      call read_selfenergies_transmission(ik,ispin,ie,fsigma,ETransmGrid)
    endif



  end subroutine get_selfenergies


  subroutine write_selfenergies_transmission(ik,ispin,ie,fsigma,energygrid1)

    use mMPI_NEGF

    integer, intent(in) :: ik,ispin,ie
    character(LEN=*), intent(in) ::  fsigma
    type(EnergyGridType), intent(in) :: energygrid1

    integer il,ihead,bytes_sigma,sigma_io,rec1
    integer MPIerror,i1,i2
    type(SelfEnergyType), allocatable :: sigma_buffer(:)
#ifdef MPI
    INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

    if(myhead==0)then

      allocate(sigma_buffer(energygrid1%nLeads))
      do il=1,energygrid1%nLeads
        call allocate_sigma_single(sigma_buffer(il),energygrid1%leadsTotalDim(il),myhead,0,(0.0_kdp,0.0_kdp),0.0_kdp,0.0_kdp)
      enddo

      bytes_sigma=0
      do il=1,energygrid1%nLeads
        bytes_sigma=bytes_sigma + 16 *  sigma_buffer(il)%n**2
      enddo

      sigma_io=22349
      if(ik==1)then
        OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='UNKNOWN', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)
      else
        OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='OLD', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)
      endif

      rec1= (ie-1) * nheads  + ETransmGrid%nEnergiesGlobal  * (ispin-1) + ETransmGrid%nEnergiesGlobal  * ETransmGrid%nspin * (ik-1)  
      rec1=rec1+1
      write(sigma_io,REC=rec1)(energygrid1%sigma(il,1,1,1)%sigma,il=1,energygrid1%nLeads)

#ifdef MPI
      do ihead=1,nheads-1
        do il=1,energygrid1%nLeads
          call MPI_RECV(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex, ihead, 1, inverseheads_comm, istatus, MPIerror)
        enddo
        rec1=rec1+1
        write(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo
#endif

      close(sigma_io)

      do il=1,energygrid1%nLeads
        deallocate(sigma_buffer(il)%sigma)
      enddo
      deallocate(sigma_buffer)


#ifdef MPI
    else


      do ihead=1,nheads-1
        if(myhead.eq.ihead)then

          do il=1,energygrid1%nLeads
            call MPI_SEND(energygrid1%sigma(il,1,1,1)%sigma,energygrid1%sigma(il,1,1,1)%n**2,DAT_dcomplex,0,1,inverseheads_comm,MPIerror)
          enddo

        endif
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo

#endif

    endif

  end SUBROUTINE write_selfenergies_transmission

  subroutine read_selfenergies_transmission(ik,ispin,ie,fsigma,energygrid1)

    use mMPI_NEGF

    integer, intent(in) :: ik,ispin,ie
    character(LEN=*), intent(in) ::  fsigma
    type(EnergyGridType), intent(inout) :: energygrid1

    integer il,ihead,bytes_sigma,sigma_io
    integer MPIerror,i1,i2
    type(SelfEnergyType), allocatable :: sigma_buffer(:)
    integer rec1
#ifdef MPI
    INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

    if(myhead==0)then

      sigma_io=22349

      allocate(sigma_buffer(energygrid1%nLeads))
      do il=1,energygrid1%nLeads
        call allocate_sigma_single(sigma_buffer(il),energygrid1%leadsTotalDim(il),myhead,0,(0.0_kdp,0.0_kdp),0.0_kdp,0.0_kdp)
      enddo

      bytes_sigma=0
      do il=1,energygrid1%nLeads
        bytes_sigma=bytes_sigma + 16 *  sigma_buffer(il)%n**2
      enddo

      OPEN(UNIT=sigma_io,FILE=fsigma,STATUS='OLD', FORM='UNFORMATTED',ACCESS='DIRECT', RECL=bytes_sigma/4)

!***********************
      rec1= (ie-1) * nheads  + ETransmGrid%nEnergiesGlobal  * (ispin-1) + ETransmGrid%nEnergiesGlobal  * ETransmGrid%nspin * (ik-1)  

      rec1=rec1+1
      read(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)
 
      do il=1,energygrid1%nLeads
        energygrid1%sigma(il,1,1,1)%sigma=sigma_buffer(il)%sigma
      enddo


#ifdef MPI
      do ihead=1,nheads-1

        rec1=rec1+1
        read(sigma_io,REC=rec1)(sigma_buffer(il)%sigma,il=1,energygrid1%nLeads)

        do il=1,energygrid1%nLeads
          call MPI_SEND(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex,ihead,1,inverseheads_comm,MPIerror)
        enddo

        CALL MPI_BARRIER(inverseheads_comm, MPIerror)

      enddo
#endif

      close(sigma_io)

      do il=1,energygrid1%nLeads
        deallocate(sigma_buffer(il)%sigma)
      enddo
      deallocate(sigma_buffer)

#ifdef MPI
    else

      allocate(sigma_buffer(energygrid1%nLeads))
      do il=1,energygrid1%nLeads
        call allocate_sigma_single(sigma_buffer(il),energygrid1%leadsTotalDim(il),myhead,0,(0.0_kdp,0.0_kdp),0.0_kdp,0.0_kdp)
      enddo

      do ihead=1,nheads-1
        if(myhead.eq.ihead)then

         do il=1,energygrid1%nLeads
           call MPI_RECV(sigma_buffer(il)%sigma,sigma_buffer(il)%n**2,DAT_dcomplex, 0, 1, inverseheads_comm, istatus, MPIerror)
           energygrid1%sigma(il,1,1,1)%sigma=sigma_buffer(il)%sigma
         enddo

        endif
        CALL MPI_BARRIER(inverseheads_comm, MPIerror)
      enddo

      do il=1,energygrid1%nLeads
        deallocate(sigma_buffer(il)%sigma)
      enddo
      deallocate(sigma_buffer)

#endif

    endif

  end SUBROUTINE read_selfenergies_transmission


  
  subroutine extractsigma(energygrid1,sigma,ie,ispin,ik)

    integer, intent(in) :: ie,ispin,ik
    type(EnergyGridType), intent(inout) :: energygrid1
    type(SelfEnergyType), intent(inout) :: sigma(energygrid1%nLeads)

    integer il,iksigma

    if((energygrid1%GridType.eq.1))then
      iksigma=1
    else
      iksigma=ik
    endif

    if(sigma(1)%InfoSigma==0)then
      do il=1,energygrid1%nLeads
        sigma(il)%sigma=energygrid1%sigma(il,ie,ispin,iksigma)%sigma
      enddo
    elseif(sigma(1)%InfoSigma==1)then
      do il=1,energygrid1%nLeads
        call single_selfenergy(sigma(il)%sigma,sigma(il)%n,il,iksigma,ispin,ie,energygrid1)
      enddo
    elseif(sigma(1)%InfoSigma==2)then
      call readsingle_selfenergy(sigma,iksigma,ispin,ie,energygrid1)
    endif

  end subroutine extractsigma

  subroutine calculate_and_write_selfenergies(ik,fsigma,energygrid1)

    use mMPI_NEGF

    integer, intent(in) :: ik
    character(LEN=*), intent(in) ::  fsigma
    type(EnergyGridType), intent(in) :: energygrid1

    integer ispin,ie,il
    type(SelfEnergyType), allocatable :: sigmaleads(:)

    allocate(sigmaleads(energygrid1%nLeads))
    do il=1,energygrid1%nLeads
      call allocate_sigma_single(sigmaleads(il), energygrid1%leadsTotalDim(il),myhead, energygrid1%InfoSigma,(0.0D0,0.0D0),0.0D0,0.0D0)
    enddo

    do ie=1,energygrid1%nEnergies
      do ispin=1,energygrid1%nspin
        do il=1,energygrid1%nLeads
          call single_selfenergy(sigmaleads(il)%sigma,sigmaleads(il)%n,il,ik,ispin,ie,energygrid1)
        enddo
        call write_single_selfenergy(sigmaleads,ik,ispin,ie,energygrid1)
      enddo
    enddo

    do il=1,energygrid1%nLeads
      deallocate(sigmaleads(il)%sigma)
    enddo
    deallocate(sigmaleads)

 
  end subroutine calculate_and_write_selfenergies


end module mEnergyGrid



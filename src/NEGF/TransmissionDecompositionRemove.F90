module mTransmissionDecompositionRemove
  
  use mTypes

  implicit none

  private
  
  public :: ComputeVg
  public :: ConvertExpandFourierSstates
  
  integer, public                       :: nTstatesS,nTstatesL,nTstatesLlg,nTstatesR,nTstatesRlg
  integer, public                       :: nopenLrg,nopenRrg,nopenLlg,nopenRlg
  double complex, allocatable, public   :: Smat(:,:)
  double complex, allocatable, public   :: utildeR(:,:),uAtildeL(:,:), utildeLlg(:,:),uAtildeRlg(:,:)
  double complex, allocatable, public   :: uR(:,:),uAL(:,:),uLlg(:,:), uARlg(:,:),zAL(:), zR(:),zLlg(:),zARlg(:)
  double complex, allocatable, public   :: kRrg(:),kALrg(:),kLlg(:),kARlg(:)
  double complex, allocatable, public   :: VgLlg(:), VgLrg(:), VgRrg(:), VgRlg(:)
  double precision, allocatable, public :: MuR(:,:),MuAL(:,:),MuLlg(:,:), MuARlg(:,:)
  double complex, allocatable, public   :: ALlglg(:,:),ALrglg(:,:), ARrgrg(:,:),ARlgrg(:,:)

  type(FourierScatteringStates), public :: FourierPhiS(2) ! element 1 is for the left electrode, element 2 for the right electrode

  contains

  subroutine ComputeVg(side,n,ene,S0,S1,km1,k1)

    use negfmod, only: ikpmod,TransmissionMatrixPDOS,ef_em,TransmissionMatrixPDOSNWrite,TransmissionMatrixiSetPhase,EM_NSPINBlocks
    use mComputeULR, only : FindNth

    integer, intent(in) :: n
    double complex, intent(in) :: ene
    CHARACTER(LEN=1),intent(in) :: side
    DOUBLE COMPLEX, intent(in) :: S0(n,n),s1(n,n)
    DOUBLE COMPLEX, intent(in) :: km1(n,n),k1(n,n)

    integer i,j,i1
    double complex vri(n),vrj(n)
    DOUBLE COMPLEX :: s01k(n,n),k1t(n,n),vg,vrbuf(n),expiphi
    DOUBLE COMPLEX :: svec(4),sbuf(4),ei
    double precision pdostolerance,alpha,kappa
    integer iSetPhase
    double COMPLEX, parameter :: ii=(0.0D0,1.0D0)

    ei=ene-ef_em

    iSetPhase=TransmissionMatrixiSetPhase
!    if(iSetPhase>0)then
!      write(12347,*)"The following element of the wave function is set to be real:",iSetPhase
!    endif

    if(side.eq.'L')then

      if(allocated(VgLrg))then
        deallocate(VgLrg,MuAL,VgLlg,MuLlg)
        allocate(zAL(nTstatesL),zLlg(nTstatesLlg))
        zAL=exp(ii * kALrg)
        zLlg=exp(ii * kLlg)
      endif
!start computing VgLrg
      allocate(VgLrg(nTstatesL))
      if(EM_NSPINBlocks==4)then
        allocate(MuAL(4,nTstatesL))
      else
        allocate(MuAL(4,1))
      endif
      do i=1,nTstatesL
        k1t=(k1 * zAL(i) - km1 / zAL(i))
        vri=uAL(:,i)
        
        s01k=s0+s1 * zAL(i) + DCONJG(TRANSPOSE(s1))/ zAL(i)
        vrbuf=matmul(s01k,vri)
        vg=0.0D0
        do j=1,n
          vg=vg+ DCONJG(vri(j)) * vrbuf(j)
        enddo
        vg=sqrt(vg)
!        write(12347,*)"vgnorm=",DREAL(ene),DREAL(vg),DIMAG(vg)
        vri=vri/vg
        
        if(iSetPhase>0)then
          if(abs(vri(iSetPhase)).ne.0.0D0) then
            expiphi=abs(vri(iSetPhase))/vri(iSetPhase)
          else
            expiphi=1.0D0
            write(12347,*)"First element of wave function is exactly 0, not setting phase shift"
          endif
        else
          expiphi=1.0D0
        endif


!        call alphaofz2(alpha,kappa,expiphi)
!        write(12347,*)"expiphiLrg=",13.6056981D0 * dreal(ei),i,alpha,kappa,abs(expiphi),DREAL(expiphi),DIMAG(expiphi),ikpmod
        vri=vri*expiphi

        uAL(:,i)=vri
        uAtildeL(:,i)=uAtildeL(:,i)*vg*expiphi

        if(EM_NSPINBlocks==4)then
          vrbuf=vrbuf*expiphi/vg
          svec=0.0D0

          if(TransmissionMatrixPDOS)then
            call FindNth(TransmissionMatrixPDOSNWrite,n,vri,vrbuf,pdostolerance)
          endif

          do j=1,n,2
            sbuf(1)= (DCONJG(vri(j)) * vrbuf(j)+DCONJG(vri(j+1)) * vrbuf(j+1))
            sbuf(2)= (DCONJG(vri(j+1)) * vrbuf(j)+DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(3)=(0.0D0,1.0D0)*(DCONJG(vri(j+1)) * vrbuf(j)-DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(4)= (DCONJG(vri(j)) * vrbuf(j)-DCONJG(vri(j+1)) * vrbuf(j+1))

            svec(1)=svec(1)+ sbuf(1)
            svec(2)=svec(2)+ sbuf(2)
            svec(3)=svec(3)+ sbuf(3)
            svec(4)=svec(4)+ sbuf(4)
            if(TransmissionMatrixPDOS.and.abs(sbuf(1)).ge.pdostolerance)then
              write(12347,*)"PDOS_L_rg=",13.6056981D0 * dreal(ei),i,(j+1)/2,dreal(sbuf(1)),dreal(sbuf(2)),dreal(sbuf(3)),dreal(sbuf(4)),1,ikpmod
            endif
          enddo
          MuAL(:,i)=DREAL(svec)
!          write(12347,*)"svec1=",DREAL(ene),DREAL(svec(1)),DIMAG(svec(1)),i
!          write(12347,*)"svec2=",DREAL(ene),DREAL(svec(2)),DIMAG(svec(2)),i
!          write(12347,*)"svec3=",DREAL(ene),DREAL(svec(3)),DIMAG(svec(3)),i
!          write(12347,*)"svec4=",DREAL(ene),DREAL(svec(4)),DIMAG(svec(4)),i
        endif
 
        vrbuf=matmul(k1t,vri)
        vg=0D0
        do j=1,n
          vg=vg+(0D0,1D0) *  DCONJG(vri(j)) * vrbuf(j)
        enddo
!        vg=vg* 13.6057d0 ! this transforms the units of the group velocity from Ry to eV
!        write(12347,*)"VgLrg=",DREAL(ene),DREAL(vg),DIMAG(vg)
        VgLrg(i)=vg
      enddo
      

!end computing VgLrg

!start computing VgLlg
      allocate(VgLlg(nTstatesLlg))
      if(EM_NSPINBlocks==4)then
        allocate(MuLlg(4,nTstatesLlg))
      else
        allocate(MuLlg(4,1))
      endif
      do i=1,nTstatesLlg
        k1t=(k1 * zLlg(i) - km1 / zLlg(i))
        vri=uLlg(:,i)
        
        s01k=s0+s1 * zLlg(i) + DCONJG(TRANSPOSE(s1))/ zLlg(i)
        vrbuf=matmul(s01k,vri)
        vg=0.0D0
        do j=1,n
          vg=vg+ DCONJG(vri(j)) * vrbuf(j)
        enddo
        vg=sqrt(vg)
!        write(12347,*)"vgnorm2=",DREAL(ene),DREAL(vg),DIMAG(vg)

        vri=vri/vg

        if(iSetPhase>0)then
          if(abs(vri(iSetPhase)).ne.0.0D0) then
            expiphi=abs(vri(iSetPhase))/vri(iSetPhase)
          else
            expiphi=1.0D0
            write(12347,*)"First element of wave function is exactly 0, not setting phase shift"
          endif
        else
          expiphi=1.0D0
        endif




!        call alphaofz2(alpha,kappa,expiphi)
!        write(12347,*)"expiphiLlg=",13.6056981D0 * dreal(ei),i,alpha,kappa,DREAL(expiphi),DIMAG(expiphi),abs(expiphi),ikpmod
        vri=vri*expiphi

        uLlg(:,i)=vri
        utildeLlg(:,i)=utildeLlg(:,i)*vg*expiphi


        if(EM_NSPINBlocks==4)then
          vrbuf=vrbuf*expiphi/vg
          svec=0.0D0

          if(TransmissionMatrixPDOS)then
            call FindNth(TransmissionMatrixPDOSNWrite,n,vri,vrbuf,pdostolerance)
          endif

          do j=1,n,2
            sbuf(1)= (DCONJG(vri(j)) * vrbuf(j)+DCONJG(vri(j+1)) * vrbuf(j+1))
            sbuf(2)= (DCONJG(vri(j+1)) * vrbuf(j)+DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(3)=(0.0D0,1.0D0)*(DCONJG(vri(j+1)) * vrbuf(j)-DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(4)= (DCONJG(vri(j)) * vrbuf(j)-DCONJG(vri(j+1)) * vrbuf(j+1))

            svec(1)=svec(1)+ sbuf(1)
            svec(2)=svec(2)+ sbuf(2)
            svec(3)=svec(3)+ sbuf(3)
            svec(4)=svec(4)+ sbuf(4)
            if(TransmissionMatrixPDOS.and.abs(sbuf(1)).ge.pdostolerance)then
              write(12347,*)"PDOS_L_lg=",13.6056981D0 * dreal(ei),i,(j+1)/2,dreal(sbuf(1)),dreal(sbuf(2)),dreal(sbuf(3)),dreal(sbuf(4)),1,ikpmod
            endif
          enddo
          MuLlg(:,i)=DREAL(svec)
        endif
  
        vrbuf=matmul(k1t,vri)
        vg=0D0
        do j=1,n
          vg=vg+(0D0,1D0) *  DCONJG(vri(j)) * vrbuf(j)
        enddo
!        vg=vg* 13.6057d0 ! this transforms the units of the group velocity from Ry to eV
!        write(12347,*)"VgLlg=",DREAL(ene),DREAL(vg),DIMAG(vg)
        VgLlg(i)=vg
      enddo
      

!end computing VgLlg

      allocate(ALrglg(nTstatesL,nTstatesLlg))

        
      do i=1,nTstatesL
        vri=uAL(:,i)

        do j=1,nTstatesLlg

          vrj=uLlg(:,j)
          s01k=s0+0.5D0* s1 * (zLlg(j)+zAL(i)) + 0.5D0 * DCONJG(TRANSPOSE(s1))* (1.0D0/ zLlg(j)+1.0D0/ zAL(i))
          vrbuf=matmul(s01k,vrj)
          vg=0.0D0
          do i1=1,n
            vg=vg+ DCONJG(vri(i1)) * vrbuf(i1)
          enddo
          ALrglg(i,j)=vg

!          call alphaofz2(alpha,kappa,ALrglg(i,j))
!          write(12347,*)"ALrglg=",13.6056981D0 * dreal(ei),i,j,abs(ALrglg(i,j)),alpha,dreal(ALrglg(i,j)),dimag(ALrglg(i,j)),ikpmod

        enddo
      enddo

      if(.false.)then

        allocate(ALlglg(nTstatesLlg,nTstatesLlg))

        do i=1,nTstatesLlg
          vri=uLlg(:,i)

          do j=1,nTstatesLlg

            vrj=uLlg(:,j)
            s01k=s0+0.5D0* s1 * (zLlg(i)+zLlg(j)) +0.5D0 *  DCONJG(TRANSPOSE(s1))* (1.0D0/ zLlg(j)+1.0D0/ zLLg(i))
            vrbuf=matmul(s01k,vrj)
            vg=0.0D0
            do i1=1,n
              vg=vg+ DCONJG(vri(i1)) * vrbuf(i1)
            enddo
            ALlglg(i,j)=vg

!            call alphaofz2(alpha,kappa,ALlglg(i,j))
!            write(12347,*)"ALlglg=",13.6056981D0 * dreal(ei),i,j,abs(ALlglg(i,j)),alpha,dreal(ALlglg(i,j)),dimag(ALlglg(i,j)),ikpmod

          enddo
        enddo


        deallocate(ALlglg)
      endif




      deallocate(zLlg)
      deallocate(zAL)

    else


      if(allocated(VgRlg))then
        deallocate(VgRlg,MuARlg,VgRrg,MuR)
        allocate(zARlg(nTstatesRlg),zR(nTstatesR))
        zARlg=exp(-ii * kARlg)
        zR=exp(-ii * kRrg)
      endif
!start computing VgRrg
      allocate(VgRrg(nTstatesR))
      if(EM_NSPINBlocks==4)then
        allocate(MuR(4,nTstatesR))
      else
        allocate(MuR(4,1))
      endif
      do i=1,nTstatesR
        k1t=(k1 / zR(i) - km1 * zR(i))
        vri=uR(:,i)
        
        s01k=s0+s1 / zR(i) + DCONJG(TRANSPOSE(s1))* zR(i)
        vrbuf=matmul(s01k,vri)
        vg=0.0D0
        do j=1,n
          vg=vg+ DCONJG(vri(j)) * vrbuf(j)
        enddo
        vg=sqrt(vg)
!        write(12347,*)"vgnorm3=",DREAL(ene),DREAL(vg),DIMAG(vg)
        vri=vri/vg


        if(iSetPhase>0)then
          if(abs(vri(iSetPhase)).ne.0.0D0) then
            expiphi=abs(vri(iSetPhase))/vri(iSetPhase)
          else
            expiphi=1.0D0
            write(12347,*)"First element of wave function is exactly 0, not setting phase shift"
          endif
        else
          expiphi=1.0D0
        endif






!        call alphaofz2(alpha,kappa,expiphi)
!        write(12347,*)"expiphiRrg=",13.6056981D0 * dreal(ei),i,alpha,kappa,DREAL(expiphi),DIMAG(expiphi),abs(expiphi),ikpmod
        vri=vri*expiphi

        uR(:,i)=vri
        utildeR(:,i)=utildeR(:,i)*vg*expiphi

        if(EM_NSPINBlocks==4)then
          vrbuf=vrbuf*expiphi/vg
          svec=0.0D0

          if(TransmissionMatrixPDOS)then
            call FindNth(TransmissionMatrixPDOSNWrite,n,vri,vrbuf,pdostolerance)
          endif

          do j=1,n,2
            sbuf(1)= (DCONJG(vri(j)) * vrbuf(j)+DCONJG(vri(j+1)) * vrbuf(j+1))
            sbuf(2)= (DCONJG(vri(j+1)) * vrbuf(j)+DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(3)=(0.0D0,1.0D0)*(DCONJG(vri(j+1)) * vrbuf(j)-DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(4)= (DCONJG(vri(j)) * vrbuf(j)-DCONJG(vri(j+1)) * vrbuf(j+1))

            svec(1)=svec(1)+ sbuf(1)
            svec(2)=svec(2)+ sbuf(2)
            svec(3)=svec(3)+ sbuf(3)
            svec(4)=svec(4)+ sbuf(4)
            if(TransmissionMatrixPDOS.and.abs(sbuf(1)).ge.pdostolerance)then
              write(12347,*)"PDOS_R_rg=",13.6056981D0 * dreal(ei),i,(j+1)/2,dreal(sbuf(1)),dreal(sbuf(2)),dreal(sbuf(3)),dreal(sbuf(4)),1,ikpmod
            endif
          enddo
          MuR(:,i)=DREAL(svec)
        endif
  
        vrbuf=matmul(k1t,vri)
        vg=0D0
        do j=1,n
          vg=vg+(0D0,1D0) *  DCONJG(vri(j)) * vrbuf(j)
        enddo
!        vg=vg* 13.6057d0 ! this transforms the units of the group velocity from Ry to eV
!        write(12347,*)"VgRrg=",DREAL(ene),DREAL(vg),DIMAG(vg)
        VgRrg(i)=vg
      enddo
      

!end computing VgRrg


!start computing VgRlg
      allocate(VgRlg(nTstatesRlg))
      if(EM_NSPINBlocks==4)then
        allocate(MuARlg(4,nTstatesRlg))
      else
        allocate(MuARlg(4,1))
      endif
      do i=1,nTstatesRlg
        k1t=(k1 / zARlg(i) - km1 * zARlg(i))
        vri=uARlg(:,i)
        
        s01k=s0+s1 / zARlg(i) + DCONJG(TRANSPOSE(s1))* zARlg(i)
        vrbuf=matmul(s01k,vri)
        vg=0.0D0
        do j=1,n
          vg=vg+ DCONJG(vri(j)) * vrbuf(j)
        enddo
        vg=sqrt(vg)
!        write(12347,*)"vgnorm3=",DREAL(ene),DREAL(vg),DIMAG(vg)
        vri=vri/vg


        if(iSetPhase>0)then
          if(abs(vri(iSetPhase)).ne.0.0D0) then
            expiphi=abs(vri(iSetPhase))/vri(iSetPhase)
          else
            expiphi=1.0D0
            write(12347,*)"First element of wave function is exactly 0, not setting phase shift"
          endif
        else
          expiphi=1.0D0
        endif

!        call alphaofz2(alpha,kappa,expiphi)
!        write(12347,*)"expiphiRrg=",13.6056981D0 * dreal(ei),i,alpha,kappa,DREAL(expiphi),DIMAG(expiphi),abs(expiphi),ikpmod
        vri=vri*expiphi

        uARlg(:,i)=vri
        uAtildeRlg(:,i)=uAtildeRlg(:,i)*vg*expiphi

        if(EM_NSPINBlocks==4)then
          vrbuf=vrbuf*expiphi/vg
          svec=0.0D0

          if(TransmissionMatrixPDOS)then
            call FindNth(TransmissionMatrixPDOSNWrite,n,vri,vrbuf,pdostolerance)
          endif

          do j=1,n,2
            sbuf(1)= (DCONJG(vri(j)) * vrbuf(j)+DCONJG(vri(j+1)) * vrbuf(j+1))
            sbuf(2)= (DCONJG(vri(j+1)) * vrbuf(j)+DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(3)=(0.0D0,1.0D0)*(DCONJG(vri(j+1)) * vrbuf(j)-DCONJG(vri(j)) * vrbuf(j+1))
            sbuf(4)= (DCONJG(vri(j)) * vrbuf(j)-DCONJG(vri(j+1)) * vrbuf(j+1))

            svec(1)=svec(1)+ sbuf(1)
            svec(2)=svec(2)+ sbuf(2)
            svec(3)=svec(3)+ sbuf(3)
            svec(4)=svec(4)+ sbuf(4)
            if(TransmissionMatrixPDOS.and.abs(sbuf(1)).ge.pdostolerance)then
              write(12347,*)"PDOS_R_lg=",13.6056981D0 * dreal(ei),i,(j+1)/2,dreal(sbuf(1)),dreal(sbuf(2)),dreal(sbuf(3)),dreal(sbuf(4)),1,ikpmod
            endif
          enddo
          MuARlg(:,i)=DREAL(svec)
        endif
  
        vrbuf=matmul(k1t,vri)
        vg=0D0
        do j=1,n
          vg=vg+(0D0,1D0) *  DCONJG(vri(j)) * vrbuf(j)
        enddo
!        vg=vg* 13.6057d0 ! this transforms the units of the group velocity from Ry to eV
!        write(12347,*)"VgRlg=",DREAL(ene),DREAL(vg),DIMAG(vg)
        VgRlg(i)=vg
      enddo
      

!end computing VgRlg

      allocate(ARlgrg(nTstatesRlg,nTstatesR))
        
      do i=1,nTstatesRlg
        vri=uARlg(:,i)

        do j=1,nTstatesR

          vrj=uR(:,j)
          s01k=s0+0.5D0* s1 * (1.0D0/zR(j)+1.0D0/zARlg(i)) + 0.5D0 * DCONJG(TRANSPOSE(s1))* (zR(j)+zARlg(i))
          vrbuf=matmul(s01k,vrj)
          vg=0.0D0
          do i1=1,n
            vg=vg+ DCONJG(vri(i1)) * vrbuf(i1)
          enddo
          ARlgrg(i,j)=vg

!          call alphaofz2(alpha,kappa,ARlgrg(i,j))
!          write(12347,*)"ARlgrg=",13.6056981D0 * dreal(ei),i,j,abs(ARlgrg(i,j)),alpha,dreal(ARlgrg(i,j)),dimag(ARlgrg(i,j)),ikpmod

        enddo
      enddo


      if(.false.)then
        allocate(ARrgrg(nTstatesR,nTstatesR))
  
        do i=1,nTstatesR
          vri=uR(:,i)

          do j=1,nTstatesR

            vrj=uR(:,j)
            s01k=s0+0.5D0* s1 * (1.0D0/zR(j)+1.0D0/zR(i)) + 0.5D0 * DCONJG(TRANSPOSE(s1))* (zR(j)+zR(i))
            vrbuf=matmul(s01k,vrj)
            vg=0.0D0
            do i1=1,n
              vg=vg+ DCONJG(vri(i1)) * vrbuf(i1)
            enddo
            ARrgrg(i,j)=vg

!            call alphaofz2(alpha,kappa,ARrgrg(i,j))
!            write(12347,*)"ARrgrg=",13.6056981D0 * dreal(ei),i,j,abs(ARrgrg(i,j)),alpha,dreal(ARrgrg(i,j)),dimag(ARrgrg(i,j)),ikpmod

          enddo
        enddo



        deallocate(ARrgrg)
      endif


      deallocate(zARlg)
      deallocate(zR)



    endif


  end subroutine ComputeVg

  SUBROUTINE ConvertExpandFourierSstates(il,FourierPhiS)

  use negfmod, only: EM_NSPINBlocks

  type(FourierScatteringStates), intent(inout) :: FourierPhiS
  integer, intent(in)         :: il

  integer ikx,iky,ndivxy(2),ind1,ind2,n,i1

  ndivxy=FourierPhiS%ndivxy
  n=FourierPhiS%n
  
  if(il==1)then
    nTstatesL=0
    nTstatesLlg=0
    do ikx=1,ndivxy(1)
      do iky=1,ndivxy(2)
        nTstatesL=nTstatesL+FourierPhiS%FourierSstates(ikx,iky)%nTstates(1)
        nTstatesLlg=nTstatesLlg+FourierPhiS%FourierSstates(ikx,iky)%nTstates(2)
      enddo
    enddo

!    write(12347,*)"nTstatesL,lg=",nTstatesL,nTstatesLlg
    allocate(uAL(n,nTstatesL))
    allocate(uAtildeL(n,nTstatesL))
    allocate(kALrg(nTstatesL))
    allocate(VgLrg(nTstatesL))
    if(EM_NSPINBlocks==4)then
      allocate(MuAL(4,nTstatesL))
    else
      allocate(MuAL(4,1))
    endif

    allocate(uLlg(n,nTstatesLlg))
    allocate(utildeLlg(n,nTstatesLlg))
    allocate(kLlg(nTstatesLlg))
    allocate(VgLlg(nTstatesLlg))
    if(EM_NSPINBlocks==4)then
      allocate(MuLlg(4,nTstatesLlg))
    else
      allocate(MuLlg(4,1))
    endif

    ind1=0
    ind2=0

    do ikx=1,ndivxy(1)
      do iky=1,ndivxy(2)
        do i1=1,FourierPhiS%FourierSstates(ikx,iky)%nTstates(1)
          ind1=ind1+1

          uAL(:,ind1)=FourierPhiS%Sstates(ikx,iky)%phi_in(:,i1)
          uAtildeL(:,ind1)=FourierPhiS%Sstates(ikx,iky)%phit_in(:,i1)
          kALrg(ind1)=FourierPhiS%Sstates(ikx,iky)%k_in(i1)
          VgLrg(ind1)=FourierPhiS%Sstates(ikx,iky)%v_in(i1)
          if(EM_NSPINBlocks==4)then
            MuAL(:,ind1)=FourierPhiS%Sstates(ikx,iky)%mu_in(:,i1)
          endif

        enddo

        do i1=1,FourierPhiS%FourierSstates(ikx,iky)%nTstates(2)
          ind2=ind2+1

          uLlg(:,ind2)=FourierPhiS%Sstates(ikx,iky)%phi_out(:,i1)
          utildeLlg(:,ind2)=FourierPhiS%Sstates(ikx,iky)%phit_out(:,i1)
          kLlg(ind2)=FourierPhiS%Sstates(ikx,iky)%k_out(i1)
          VgLlg(ind2)=FourierPhiS%Sstates(ikx,iky)%v_out(i1)
          if(EM_NSPINBlocks==4)then
            MuLlg(:,ind2)=FourierPhiS%Sstates(ikx,iky)%mu_out(:,i1)
          endif
      
        enddo
      enddo
    enddo
  endif

  if(il==2)then
    nTstatesR=0
    nTstatesRlg=0
    do ikx=1,ndivxy(1)
      do iky=1,ndivxy(2)
        nTstatesR=nTstatesR+FourierPhiS%FourierSstates(ikx,iky)%nTstates(2)
        nTstatesRlg=nTstatesrlg+FourierPhiS%FourierSstates(ikx,iky)%nTstates(1)
      enddo
    enddo
!    write(12347,*)"nTstatesR,lg=",nTstatesR,nTstatesRlg

    allocate(uARlg(n,nTstatesRlg))
    allocate(uAtildeRlg(n,nTstatesRlg))
    allocate(kARlg(nTstatesRlg))
    allocate(VgRlg(nTstatesRlg))
    if(EM_NSPINBlocks==4)then
      allocate(MuARlg(4,nTstatesRlg))
    else
      allocate(MuARlg(4,1))
    endif

    allocate(uR(n,nTstatesR))
    allocate(utildeR(n,nTstatesR))
    allocate(kRrg(nTstatesR))
    allocate(VgRrg(nTstatesR))
    if(EM_NSPINBlocks==4)then
      allocate(MuR(4,nTstatesR))
    else
      allocate(MuR(4,1))
    endif

    ind1=0
    ind2=0

    do ikx=1,ndivxy(1)
      do iky=1,ndivxy(2)
        do i1=1,FourierPhiS%FourierSstates(ikx,iky)%nTstates(2)
          ind1=ind1+1

          uARlg(:,ind1)=FourierPhiS%Sstates(ikx,iky)%phi_in(:,i1)
          uAtildeRlg(:,ind1)=FourierPhiS%Sstates(ikx,iky)%phit_in(:,i1)
          kARlg(ind1)=FourierPhiS%Sstates(ikx,iky)%k_in(i1)
          VgRlg(ind1)=FourierPhiS%Sstates(ikx,iky)%v_in(i1)
          if(EM_NSPINBlocks==4)then
            MuARlg(:,ind1)=FourierPhiS%Sstates(ikx,iky)%mu_in(:,i1)
          endif

        enddo

        do i1=1,FourierPhiS%FourierSstates(ikx,iky)%nTstates(1)
          ind2=ind2+1

          uR(:,ind2)=FourierPhiS%Sstates(ikx,iky)%phi_out(:,i1)
          utildeR(:,ind2)=FourierPhiS%Sstates(ikx,iky)%phit_out(:,i1)
          kRrg(ind2)=FourierPhiS%Sstates(ikx,iky)%k_out(i1)
          VgRrg(ind2)=FourierPhiS%Sstates(ikx,iky)%v_out(i1)
          if(EM_NSPINBlocks==4)then
            MuR(:,ind2)=FourierPhiS%Sstates(ikx,iky)%mu_out(:,i1)
          endif
      
        enddo
      enddo
    enddo
  endif

  end SUBROUTINE ConvertExpandFourierSstates




end module mTransmissionDecompositionRemove


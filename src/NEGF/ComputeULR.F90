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
!                   COMPUTEULR,
!                   COMPUTEULRSVD,
!                   ALPHAOFZ2,
!                   FINDNTH  
! AND
! THE MODULE
!                   MCOMPUTEULR  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

module mComputeULR

  use mTypes
  implicit none

  private
  
  public :: ComputeULR
  public :: InitPhiS
  public :: SetKPhiS
  public :: ComputePhiSVD
  public :: ComputeVgMuPhi
  public :: SortSstates
  public :: FindNth
  public :: FreePhiS
  type(FourierScatteringStates), public :: PhiS(2) ! element 1 is for the left electrode, element 2 for the right electrode,in general make it allocatable for more leads
  
  contains
      
  subroutine ComputeULR(side,n,zvr,zvl,ene,vrg,vrgb,vlg,vlgb,sigma,tolki)

! *****************************************************************
! Written by Ivan Rungger, October 2013
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    integer, intent(in)          :: n
    double complex, intent(in)   :: zvr(n),zvl(n),ene
    double complex, intent(in)   :: vrg(n,n),vrgb(n,n),vlg(n,n),vlgb(n,n),sigma(n,n)
    double precision, intent(in) :: tolki
    CHARACTER(LEN=1), intent(in) :: side

    If (side .eq. 'L') then
      call FreeSstates(PhiS(1)%FourierSstates(PhiS(1)%ik(1),PhiS(1)%ik(2)))
      call ComputePhiS(PhiS(1)%FourierSstates(PhiS(1)%ik(1),PhiS(1)%ik(2)),side,n,zvr,zvl,ene,vrg,vrgb,vlg,vlgb,sigma,tolki,1)
      call ComputePhiS(PhiS(1)%FourierSstates(PhiS(1)%ik(1),PhiS(1)%ik(2)),side,n,zvr,zvl,ene,vrg,vrgb,vlg,vlgb,sigma,tolki,2)
    else
      call FreeSstates(PhiS(2)%FourierSstates(PhiS(2)%ik(1),PhiS(2)%ik(2)))
      call ComputePhiS(PhiS(2)%FourierSstates(PhiS(2)%ik(1),PhiS(2)%ik(2)),side,n,1.0D0/zvl,1.0D0/zvr,ene,vlg,vlgb,vrg,vrgb,sigma,tolki,1)
      call ComputePhiS(PhiS(2)%FourierSstates(PhiS(2)%ik(1),PhiS(2)%ik(2)),side,n,1.0D0/zvl,1.0D0/zvr,ene,vlg,vlgb,vrg,vrgb,sigma,tolki,2)
    endif

  end subroutine ComputeULR

  SUBROUTINE InitPhiS(Phi,side,ndivxy,n,e)

  integer, intent(in)           :: ndivxy(2)
  CHARACTER(LEN=1),intent(in)   :: SIDE
  integer, intent(in)           :: n
  double complex, intent(in)    :: e
  type(FourierScatteringStates), intent(out) :: Phi

  Phi%n=n
  Phi%ns=n/(ndivxy(1)*ndivxy(2))
  Phi%ndivxy(1)=ndivxy(1)
  Phi%ndivxy(2)=ndivxy(2)
  Phi%e=e
  Phi%side=side

  Phi%ik(1)=1
  Phi%ik(2)=1

  allocate(Phi%kxy(2,ndivxy(1),ndivxy(2)))
  allocate(Phi%FourierSstates(ndivxy(1),ndivxy(2)))

  end SUBROUTINE InitPhiS


  SUBROUTINE FreePhiS(Phi)

  type(FourierScatteringStates), intent(inout) :: Phi

  integer ikx,iky

  do ikx=1,Phi%ndivxy(1)
    do iky=1,Phi%ndivxy(2)
      call FreeSstates(Phi%FourierSstates(ikx,iky))
    enddo
  enddo

  Phi%n=0
  Phi%ns=0
  Phi%ndivxy(1)=0
  Phi%ndivxy(2)=0
  Phi%e=0.0D0
  Phi%side='n'

  Phi%ik(1)=0
  Phi%ik(2)=0

  deallocate(Phi%kxy)
  deallocate(Phi%FourierSstates)

  end SUBROUTINE FreePhiS

  SUBROUTINE FreeSstates(Phi)

  type(ScatteringStates), intent(inout) :: Phi

  Phi%nTstates=0
  Phi%n=0

  Phi%e=0.0D0
  Phi%side='n'

  if(allocated(Phi%phi_in))deallocate(Phi%phi_in)
  if(allocated(Phi%phit_in))deallocate(Phi%phit_in)
  if(allocated(Phi%k_in))deallocate(Phi%k_in)
  if(allocated(Phi%v_in))deallocate(Phi%v_in)
  if(allocated(Phi%z_in))deallocate(Phi%z_in)
  if(allocated(Phi%mu_in))deallocate(Phi%mu_in)
  if(allocated(Phi%phiAll_in))deallocate(Phi%phiAll_in)
  if(allocated(Phi%phitAll_in))deallocate(Phi%phitAll_in)

  if(allocated(Phi%phi_out))deallocate(Phi%phi_out)
  if(allocated(Phi%phit_out))deallocate(Phi%phit_out)
  if(allocated(Phi%k_out))deallocate(Phi%k_out)
  if(allocated(Phi%v_out))deallocate(Phi%v_out)
  if(allocated(Phi%z_out))deallocate(Phi%z_out)
  if(allocated(Phi%mu_out))deallocate(Phi%mu_out)
  if(allocated(Phi%phiAll_out))deallocate(Phi%phiAll_out)
  if(allocated(Phi%phitAll_out))deallocate(Phi%phitAll_out)

  end SUBROUTINE FreeSstates


  SUBROUTINE SetKPhiS(Phi,ndivxy,k)

  integer, intent(in)           :: ndivxy(2)
  double precision, intent(in)  :: k(2,ndivxy(1),ndivxy(2))
  type(FourierScatteringStates), intent(inout) :: Phi

  Phi%kxy=k

  end SUBROUTINE SetKPhiS



  subroutine ComputePhiS(Phi,side,n,zvr,zvl,ene,vrg,vrgb,vlg,vlgb,sigma,tol,InOut)

    integer, intent(in)          :: n
    double precision, intent(in) :: tol
    double complex, intent(in)   :: vrg(n,n),vrgb(n,n),vlg(n,n),vlgb(n,n),sigma(n,n)
    double complex, intent(in)   :: zvr(n),zvl(n),ene
    CHARACTER(LEN=1), intent(in) :: side
    type(ScatteringStates), intent(inout) :: Phi
    integer, intent(in)          :: InOut ! 1 is for incoming states, 2 is four outgoing states

    integer IPIV(n),info, il,nopen,ncl
    integer iprop(n) 
    double precision alphak,lz
    double complex, allocatable  ::  work(:)

    Phi%n=n
    Phi%e=ene
    Phi%side=side

    if(inout==1)then
      allocate(Phi%phiAll_in(n,n))
    else
      allocate(Phi%phiAll_out(n,n))
    endif
    nopen=0
    ncl=0
    iprop=0
    do il=1,n
      if(inout==1)then
        call alphaofz2(alphak,lz,zvr(il))
      else
        call alphaofz2(alphak,lz,zvl(il))
      endif

      if(abs(lz)<tol)then
        nopen=nopen+1
        ncl=ncl+1
        iprop(ncl)=il
      endif

    enddo
! here we can set nTstates=n if we want all states for testing
! in principle we can choose nTstates any value between
! nopen and n
    if(.true.)then
      Phi%nTstates(InOut)=nopen
    else
      Phi%nTstates(InOut)=n
    endif

    if(InOut==1)then
      allocate(Phi%phi_in(n,Phi%nTstates(InOut)))
      allocate(Phi%z_in(Phi%nTstates(InOut)))
      allocate(Phi%k_in(Phi%nTstates(InOut)))
      call set_values_Phi(Phi%nTstates(InOut),Phi%phi_in,Phi%phiAll_in,Phi%z_in,Phi%k_in,    vrg,vlg,zvr,zvl,vrg,vlg,zvr,zvl,n,nopen,iprop,tol)
    else
      allocate(Phi%phi_out(n,Phi%nTstates(InOut)))
      allocate(Phi%z_out(Phi%nTstates(InOut)))
      allocate(Phi%k_out(Phi%nTstates(InOut)))
      call set_values_Phi(Phi%nTstates(InOut),Phi%phi_out,Phi%phiAll_out,Phi%z_out,Phi%k_out,vlg,vrg,zvl,zvr,vrg,vlg,zvr,zvl,n,nopen,iprop,tol)
    endif

    allocate(work(n**2))
    if(InOut==1)then
      CALL ZGETRF(n,n,Phi%phiAll_in,n,IPIV,INFO)
      CALL ZGETRI(n,Phi%phiAll_in,n,IPIV,work,n**2,INFO)
      Phi%phiAll_in=DCONJG(transpose(Phi%phiAll_in))

      allocate(Phi%phit_in(n,Phi%nTstates(InOut)))
      Phi%phit_in=Phi%phiAll_in(:,1:Phi%nTstates(InOut))
      deallocate(Phi%phiAll_in)
    else
      CALL ZGETRF(n,n,Phi%phiAll_out,n,IPIV,INFO)
      CALL ZGETRI(n,Phi%phiAll_out,n,IPIV,work,n**2,INFO)
      Phi%phiAll_out=DCONJG(transpose(Phi%phiAll_out))

      allocate(Phi%phit_out(n,Phi%nTstates(InOut)))
      Phi%phit_out=Phi%phiAll_out(:,1:Phi%nTstates(InOut))
      deallocate(Phi%phiAll_out)
    endif
    deallocate(work)

  end subroutine ComputePhiS

  subroutine set_values_Phi(nTstates,phi,phiAll,z,k,vrg,vlg,zvr,zvl,vrg2,vlg2,zvr2,zvl2,n,nopen,iprop,tol)

    integer, intent(in)          :: nTstates
    integer, intent(in)          :: nopen
    integer, intent(in)          :: n
    integer, intent(in)          :: iprop(n) 
    double precision, intent(in) :: tol
    double complex, intent(in)   :: vrg(n,n),vlg(n,n)
    double complex, intent(in)   :: zvr(n),zvl(n)
    double complex, intent(in)   :: vrg2(n,n),vlg2(n,n)
    double complex, intent(in)   :: zvr2(n),zvl2(n)

    double complex, intent(out) ::  phi(n,nTstates) 
    double complex, intent(out) ::  k(nTstates) 
    double complex, intent(out) ::  z(nTstates) 
    double complex, intent(out) ::  phiAll(n,n) 
    double precision alphak,lz
    integer il,ncl
    double complex zvu

    do il=1,nopen
      phi(:,il)=vrg(:,iprop(il))
      phiAll(:,il)=vrg(:,iprop(il))
      z(il)=zvr(iprop(il))
      call alphaofz2(alphak,lz,z(il))
      k(il)=alphak+(0.0D0,1.0D0)*lz
    enddo

    ncl=nopen
    do il=1,n

      zvu=zvl2(il)
      call alphaofz2(alphak,lz,zvu)

      if(abs(lz)>=tol)then
        ncl=ncl+1
        if(ncl<=nTstates) then
          phi(:,ncl)=vlg2(:,il)
          z(ncl)=zvl2(il)
          k(ncl)=alphak+(0.0D0,1.0D0)*lz
        endif
        phiAll(:,ncl)=vlg2(:,il)
      endif

    enddo

  end subroutine set_values_Phi


  subroutine ComputePhiSVD(side,n,m,matdi,matc,vpm1,uvt,Phi,InOut)

    use mTypes

    integer, intent(in)         :: n,m,InOut
    CHARACTER(LEN=1),intent(in) :: side
    DOUBLE COMPLEX, intent(in)  :: matc(n-m,m),matdi(n-m,n-m)
    DOUBLE COMPLEX, intent(in)  :: vpm1(n-m,m)
    DOUBLE COMPLEX, intent(in)  :: uvt(n,n)
    type(ScatteringStates), intent(inout) :: Phi

    double complex, allocatable :: vttoppart(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: buf(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: qu(:,:),vrlt4(:,:),qprime(:,:)
    integer i1,i2
    integer nTstates

    
    nTstates=Phi%nTstates(InOut)

    Phi%n=n
!start computing phi_tilde
    allocate(buf(m,nTstates))
    if(InOut==1)then
      buf=Phi%phit_in
      deallocate(Phi%phit_in)
      allocate(Phi%phit_in(n,nTstates))
    else
      buf=Phi%phit_out
      deallocate(Phi%phit_out)
      allocate(Phi%phit_out(n,nTstates))
    endif

    allocate(vttoppart(n,m))
    if(side.eq.'L')then
      vttoppart=DCONJG(transpose(uvt(1:m,:)))
    else
      vttoppart=uvt(:,1:m)
    endif

    if(InOut==1)then
      Phi%phit_in=matmul(vttoppart,buf)
    else
      Phi%phit_out=matmul(vttoppart,buf)
    endif
    deallocate(buf)
    deallocate(vttoppart)
!end computing phi_tilde
      
!start computing phi
    allocate(vrlt4(m,nTstates))
    allocate(qu(n-m,nTstates))


!note that Phi%z_in=e^(+i k_n) for both left and right side
    do i1=1,nTstates
      if(side.eq.'L')then
        if(InOut==1)then
          vrlt4(:,i1)=Phi%phi_in(:,i1) * Phi%z_in(i1)
        else
          vrlt4(:,i1)=Phi%phi_out(:,i1) * Phi%z_out(i1)
        endif
      else
        if(InOut==1)then
          vrlt4(:,i1)=Phi%phi_in(:,i1) / Phi%z_in(i1)
        else
          vrlt4(:,i1)=Phi%phi_out(:,i1) / Phi%z_out(i1)
        endif
      endif
    enddo

    if(InOut==1)then
      qu=-matmul(matdi,matmul(vpm1,vrlt4)+matmul(matc,Phi%phi_in))
    else
      qu=-matmul(matdi,matmul(vpm1,vrlt4)+matmul(matc,Phi%phi_out))
    endif
    deallocate(vrlt4)


    allocate(qprime(n,nTstates))

    if(InOut==1)then
      qprime(1:m,:)=Phi%phi_in
      deallocate(Phi%phi_in)
    else
      qprime(1:m,:)=Phi%phi_out
      deallocate(Phi%phi_out)
    endif

    qprime(m+1:n,:)=qu
    deallocate(qu)

    if(InOut==1)then
      allocate(Phi%phi_in(n,nTstates))
      if(side.eq.'L')then
        Phi%phi_in=matmul(DCONJG(TRANSPOSE(uvt)),qprime)
      else
        Phi%phi_in=matmul(uvt,qprime)
      endif
    else
      allocate(Phi%phi_out(n,nTstates))
      if(side.eq.'L')then
        Phi%phi_out=matmul(DCONJG(TRANSPOSE(uvt)),qprime)
      else
        Phi%phi_out=matmul(uvt,qprime)
      endif
    endif
    deallocate(qprime)
!end computing phi

  end subroutine ComputePhiSVD


  subroutine ComputeVgMuPhiGet(Phi,side,n,ene,S0,S1,km1,k1,InOut,nam)

    use negfmod, only: ikpmod,TransmissionMatrixPDOS,ef_em,TransmissionMatrixPDOSNWrite,TransmissionMatrixiSetPhase,EM_NSPINBlocks

    integer, intent(in)                   :: n
    integer, intent(in)                   :: InOut
    double complex, intent(in)            :: ene
    CHARACTER(LEN=1),intent(in)           :: side
    DOUBLE COMPLEX, intent(in)            :: S0(n,n),s1(n,n)
    DOUBLE COMPLEX, intent(in)            :: km1(n,n),k1(n,n)
    CHARACTER(LEN=*), intent(in)          :: nam
    type(ScatteringStates), intent(inout) :: Phi


    integer nTstates
    integer i,j
    double complex vri(n)
    DOUBLE COMPLEX :: s01k(n,n),k1t(n,n),vg,vrbuf(n),expiphi
    DOUBLE COMPLEX :: svec(4),sbuf(4)
    double precision pdostolerance
    integer iSetPhase
    double complex z(Phi%nTstates(InOut))

    nTstates=Phi%nTstates(InOut)

    iSetPhase=TransmissionMatrixiSetPhase
!    if(iSetPhase>0)then
!      write(12347,*)"The following element of the wave function is set to be real:",iSetPhase
!    endif

    if(InOut==1)then
      z=Phi%z_in
      allocate(Phi%v_in(nTstates))
      if(EM_NSPINBlocks==4)then
        allocate(Phi%mu_in(4,nTstates))
      else
        allocate(Phi%mu_in(4,1))
      endif
    else
      z=Phi%z_out
      allocate(Phi%v_out(nTstates))
      if(EM_NSPINBlocks==4)then
        allocate(Phi%mu_out(4,nTstates))
      else
        allocate(Phi%mu_out(4,1))
      endif
    endif

    do i=1,nTstates
      if(InOut==1)then
        vri=Phi%phi_in(:,i)
      else
        vri=Phi%phi_out(:,i)
      endif
      
      s01k=s0+s1 * z(i) + DCONJG(TRANSPOSE(s1))/ z(i)
      vrbuf=matmul(s01k,vri)
      vg=0.0D0
      do j=1,n
        vg=vg+ DCONJG(vri(j)) * vrbuf(j)
      enddo
      vg=sqrt(vg)
      vri=vri/vg
      
      if(iSetPhase>0)then
        if(abs(vri(iSetPhase)).ne.0.0D0) then
          expiphi=abs(vri(iSetPhase))/vri(iSetPhase)
        else
          expiphi=1.0D0
!          write(12347,*)"First element of wave function is exactly 0, not setting phase shift"
        endif
      else
        expiphi=1.0D0
      endif
      vri=vri*expiphi

      if(InOut==1)then
        Phi%phi_in(:,i)=vri
        Phi%phit_in(:,i)=Phi%phit_in(:,i)*vg*expiphi
      else
        Phi%phi_out(:,i)=vri
        Phi%phit_out(:,i)=Phi%phit_out(:,i)*vg*expiphi
      endif

      if(EM_NSPINBlocks==4)then
!this is only correct if EM.OrderN T
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
            write(12347,*)nam,13.6056981D0 * dreal(ene-ef_em),i,(j+1)/2,dreal(sbuf(1)),dreal(sbuf(2)),dreal(sbuf(3)),dreal(sbuf(4)),1,ikpmod
          endif
        enddo
        if(InOut==1)then
          Phi%mu_in(:,i)=DREAL(svec)
        else
          Phi%mu_out(:,i)=DREAL(svec)
        endif
!xxx noncollinear version still to check
      endif
 
      k1t=k1 * z(i) - km1 / z(i)
      vrbuf=matmul(k1t,vri)
      vg=0D0
      do j=1,n
        vg=vg+(0D0,1D0) *  DCONJG(vri(j)) * vrbuf(j)
      enddo
!     the units of the group velocity are in Ry
!     vg=vg* 13.6057d0 ! this transforms the units of the group velocity from Ry to eV
      if(InOut==1)then
        Phi%v_in(i)=vg
      else
        Phi%v_out(i)=vg
      endif

    enddo

  end subroutine ComputeVgMuPhiGet


  subroutine ComputeVgMuPhi(side,n,ene,S0,S1,km1,k1)

    use negfmod, only: ikpmod,TransmissionMatrixPDOS,ef_em,TransmissionMatrixPDOSNWrite,TransmissionMatrixiSetPhase,EM_NSPINBlocks

    integer, intent(in)         :: n
    double complex, intent(in)  :: ene
    CHARACTER(len=1),intent(in) :: side
    double complex, intent(in)  :: S0(n,n),s1(n,n)
    double complex, intent(in)  :: km1(n,n),k1(n,n)

    If (side .eq. 'L') then
       call ComputeVgMuPhiGet(PhiS(1)%FourierSstates(PhiS(1)%ik(1),PhiS(1)%ik(2)),side,n,ene,S0,S1,km1,k1,1,"PDOS_L_rg")
       call ComputeVgMuPhiGet(PhiS(1)%FourierSstates(PhiS(1)%ik(1),PhiS(1)%ik(2)),side,n,ene,S0,S1,km1,k1,2,"PDOS_L_lg")
    else
       call ComputeVgMuPhiGet(PhiS(2)%FourierSstates(PhiS(2)%ik(1),PhiS(2)%ik(2)),side,n,ene,S0,S1,km1,k1,1,"PDOS_R_lg")
       call ComputeVgMuPhiGet(PhiS(2)%FourierSstates(PhiS(2)%ik(1),PhiS(2)%ik(2)),side,n,ene,S0,S1,km1,k1,2,"PDOS_R_rg")
    endif

  end subroutine ComputeVgMuPhi


  subroutine FindNth(k,n,vri,vrbuf,pdostolerance)

  integer, intent(in)           :: k,n
  double complex, intent(in)    :: vri(n)
  double complex, intent(in)    :: vrbuf(n)
  double precision, intent(out) :: pdostolerance

  integer i1,i2,nh,k2
  DOUBLE precision, allocatable :: sbuflist(:)
  integer maxindex
  double precision maxvalue

  nh=n/2
  k2=k
  if(k>n)k2=nh
  if(k<0)k2=nh

  allocate(sbuflist(nh))

  do i1=1,n,2
    sbuflist((i1+1)/2)= abs(DCONJG(vri(i1)) * vrbuf(i1)+DCONJG(vri(i1+1)) * vrbuf(i1+1))
  enddo

  do i1=1,k2
    maxindex=i1
    maxvalue=sbuflist(i1)
    do i2=i1+1,nh
      if(sbuflist(i2)>maxvalue)then
        maxvalue=sbuflist(i2)
        maxindex=i2
      endif
    enddo
    maxvalue=sbuflist(i1)
    sbuflist(i1)=sbuflist(maxindex)
    sbuflist(maxindex)=maxvalue
  enddo

  pdostolerance=sbuflist(k2)

!  do i1=1,nh
!    if(sbuflist(i1)>=pdostolerance)then
!      write(12347,*)"sbuflist G =",i1,sbuflist(i1),pdostolerance
!    else
!      write(12347,*)"sbuflist S =",i1,sbuflist(i1),pdostolerance
!    endif
!  enddo

  deallocate(sbuflist)

  end subroutine FindNth

  subroutine SortSstates(e,nl,n,k,phi,phit,v,mu,nspinblocks)

  integer, intent(in)             :: n,nl,nspinblocks
  double complex, intent(in)      :: e
  double complex, intent(inout)   :: k(n)
  double complex, intent(inout)   :: v(n)
  double precision, intent(inout) :: mu(4,n)
  double complex, intent(inout)   :: phi(nl,n)
  double complex, intent(inout)   :: phit(nl,n)

  integer i1,i2
  integer, allocatable :: ind(:),ind2(:)
  double precision, allocatable :: klist(:)
  double complex, allocatable :: phibuf(:)
  double precision  :: mubuf(4)
  double complex buf

  allocate(klist(n))
  allocate(ind(n),ind2(n))
  klist=k
  do i1=1,n
   ind(i1)=i1
  enddo

  call ssort(klist, ind, n, 2)

  i2=0
  do i1=1,n
    i2=i2+1
    ind2(ind(i1))=i2
  enddo

  allocate(phibuf(nl))
  do i1=1,n
    buf=k(i1)
    k(i1)=k(ind(i1))
    k(ind(i1))=buf

    buf=v(i1)
    v(i1)=v(ind(i1))
    v(ind(i1))=buf

    if(nspinblocks==4)then
      mubuf=mu(:,i1)
      mu(:,i1)=mu(:,ind(i1))
      mu(:,ind(i1))=mubuf
    endif

    phibuf=phi(:,i1)
    phi(:,i1)=phi(:,ind(i1))
    phi(:,ind(i1))=phibuf

    phibuf=phit(:,i1)
    phit(:,i1)=phit(:,ind(i1))
    phit(:,ind(i1))=phibuf

    ind2(ind(i1))=ind2(i1)
    ind(ind2(i1))=ind(i1)
  enddo

!!!! normalizing first element to 1; uncomment for testing
!!!  do i1=1,n
!!!    buf=phi(1,i1)
!!!    phi(:,i1)=phi(:,i1)/buf
!!!    phit(:,i1)=phit(:,i1)*DCONJG(buf)
!!!  enddo

  deallocate(phibuf,klist,ind,ind2)
 
  end subroutine SortSstates


  subroutine alphaofz2(alpha,kappa,zev)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************


    DOUBLE COMPLEX, intent(in)    ::  zev
    DOUBLE PRECISION, intent(out) ::  alpha,kappa
    DOUBLE PRECISION  lz,zr,zi,alpha0,pi

    pi=3.141592653589793D0

    lz=CDABS(zev)
    zr=DREAL(zev)
    zi=DIMAG(zev)
    if(zr.ge.0D0.and.zi.ge.0D0) then
      alpha0=0D0
    elseif(zr.ge.0D0.and.zi.lt.0D0) then
      zr=-DIMAG(zev)
      zi=DREAL(zev)
      alpha0=-pi/2D0
    elseif(zr.lt.0D0.and.zi.ge.0D0) then
      zr=DIMAG(zev)
      zi=-DREAL(zev)
      alpha0=pi/2D0
    elseif(zr.lt.0D0.and.zi.lt.0D0) then
      zr=-DREAL(zev)
      zi=-DIMAG(zev)
      alpha0=-pi
    endif

    alpha=ATAN(zi/zr)+alpha0
!    alpha=ATAN((zi/lz)/(zr/lz))+alpha0
    kappa=DLOG(lz)

  end subroutine alphaofz2

end module mComputeULR


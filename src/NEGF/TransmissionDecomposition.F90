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
!                   TRANSMISSIONSPINCOMPONENTS_GENERAL,
!                   TRANSMISSIONSPINCOMPONENTS,
!                   TRCMAT4,
!                   TRCMAT_2TERMS,
!                   TRCMAT,
!                   TRCMAT3,
!                   TRCMAT2,
!                   TRANSMISSIONCHANNELSDECOMPOSITION,
!                   WRITEWFS,
!                   NORMALIZEPSI,
!                   SHIFTRIGHTSELFENERGY,
!                   SENEMATMUL,
!                   FINDNRN,
!                   TMATRIX,
!                   ALPHAOFZ2  
! AND
! THE MODULE
!                   MTRANSMISSIONDECOMPOSITION  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
module mTransmissionDecomposition
  
  use mTypes

  implicit none

  private
  
  public :: TransmissionSpinComponents_general
  public :: TransmissionChannelsDecomposition
  public :: ShiftRightSelfenergy
  public :: FindnrN
  public :: Tmatrix

  contains

  subroutine TransmissionSpinComponents_general(GfLR,n1,nl,nr,sigmal,sigmar,nspin,Tr4G,nTr4G,ene)

  use mTypes
  use mMatrixUtil

  implicit none

      
  integer, intent(in) :: n1,nl,nr,nspin,nTr4G
  DOUBLE COMPLEX,intent(in) :: sigmal(nl,nl,nspin)
  DOUBLE COMPLEX,intent(in) :: sigmar(nr,nr,nspin)
  DOUBLE COMPLEX, intent(in) :: GfLR(n1,n1,nspin)
  double complex, intent(out) :: Tr4G(nTr4G)
  double precision, intent(in) :: ene

  DOUBLE COMPLEX, ALLOCATABLE :: gammalDOWN(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gammarDOWN(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: GfLRDOWN(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: GfLRUP(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gammalUP(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gammarUP(:,:)
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)

  type(matrixtype), ALLOCATABLE :: gammas(:)
  type(matrixtype), ALLOCATABLE :: GFs(:)

  integer i,i1,i2,i3,i4
  type(ioType) :: io

  io%isDebug=.false.


  allocate(GfLRUP(n1,n1))
  allocate(gammalUP(nl,nl))
  allocate(gammarUP(nr,nr))
  GfLRUP=GfLR(:,:,1)
  gammalUP=zi*(sigmal(:,:,1)-DCONJG(TRANSPOSE(sigmal(:,:,1))))
  gammarUP=zi*(sigmar(:,:,1)-DCONJG(TRANSPOSE(sigmar(:,:,1))))

  allocate(GfLRDOWN(n1,n1))
  allocate(gammalDOWN(nl,nl))
  allocate(gammarDOWN(nr,nr))
  GfLRDOWN=GfLR(:,:,2)
  gammalDOWN=zi*(sigmal(:,:,2)-DCONJG(TRANSPOSE(sigmal(:,:,2))))
  gammarDOWN=zi*(sigmar(:,:,2)-DCONJG(TRANSPOSE(sigmar(:,:,2))))


  allocate(gammas(4))
  allocate(GFs(2))

  call AllocateMatrix(nl,nl, gammas(1),"TSg",io)
  call AllocateMatrix(nl,nl, gammas(2),"TSg",io)
  call AllocateMatrix(nr,nr, gammas(3),"TSg",io)
  call AllocateMatrix(nr,nr, gammas(4),"TSg",io)

  call AllocateMatrix(n1,n1, GFs(1),"TSg",io)
  call AllocateMatrix(n1,n1, GFs(2),"TSg",io)


  gammas(1)%a=0.5D0 * (gammalUP+gammalDOWN)
  gammas(2)%a=0.5D0 * (gammalUP-gammalDOWN)
  gammas(3)%a=0.5D0 * (gammarUP+gammarDOWN)
  gammas(4)%a=0.5D0 * (gammarUP-gammarDOWN)

  GFs(1)%a=0.5D0 * (GfLRUP+GfLRDOWN)
  GFs(2)%a=0.5D0 * (GfLRUP-GfLRDOWN)

!  gfs(1)%a=0.0D0
!  gfs(2)%a=0.0D0
!
!  do i=1,n1
!    gfs(1)%a(i,i)=1.0D0
!    gfs(2)%a(i,i)=1.0D0
!  enddo


  Tr4G=0.0D0

  call trcmat(gammalUP,GfLRUP(n1-nr+1:n1,1:nl),gammarUP,GfLRUP(n1-nr+1:n1,1:nl),nl,nr,Tr4G(nTr4G-1))
  call trcmat(gammalDOWN,GfLRDOWN(n1-nr+1:n1,1:nl),gammarDOWN,GfLRDOWN(n1-nr+1:n1,1:nl),nl,nr,Tr4G(nTr4G))

  deallocate(GfLRUP)
  deallocate(gammalUP)
  deallocate(gammarUP)
  deallocate(GfLRDOWN)
  deallocate(gammalDOWN)
  deallocate(gammarDOWN)


  if(.true.)then

!Important: we assume that Gamma_L is not spin-polarized, so that we can skip the multiplications involving Gamma_L^S

    i=0
    do i1=1,2
      do i3=3,4
        do i2=1,2
          do i4=1,2
            i=i+1
            if(i1==2)cycle
            call trcmat(gammas(i1)%a,GFs(i2)%a(n1-nr+1:n1,1:nl),gammas(i3)%a,GFs(i4)%a(n1-nr+1:n1,1:nl),nl,nr,Tr4G(i))
!            write(12347,*)"Tr4G_i=",ene,i,DREAL(Tr4G(i)),DIMAG(Tr4G(i))
          enddo
        enddo
      enddo
    enddo

    do i1=3,4
      do i3=1,2
        do i2=1,2
          do i4=1,2
            i=i+1
            if(i3==2)cycle
            call trcmat(gammas(i1)%a,GFs(i2)%a(1:nl,n1-nr+1:n1),gammas(i3)%a,GFs(i4)%a(1:nl,n1-nr+1:n1),nr,nl,Tr4G(i))
          enddo
        enddo
      enddo
    enddo

    do i1=1,2
      do i3=1,2
        do i2=1,2
          do i4=1,2
            i=i+1
            if(i1==2)cycle
            if(i3==2)cycle
            call trcmat(gammas(i1)%a,GFs(i2)%a(1:nl,1:nl),gammas(i3)%a,GFs(i4)%a(1:nl,1:nl),nl,nl,Tr4G(i))
          enddo
        enddo
      enddo
    enddo

    do i1=3,4
      do i3=3,4
        do i2=1,2
          do i4=1,2
            i=i+1
            call trcmat(gammas(i1)%a,GFs(i2)%a(n1-nr+1:n1,n1-nr+1:n1),gammas(i3)%a,GFs(i4)%a(n1-nr+1:n1,n1-nr+1:n1),nr,nr,Tr4G(i))
          enddo
        enddo
      enddo
    enddo

  endif

  if(.true.)then

    do i1=1,2
      do i2=1,2
        i=i+1
        if(i1==2)cycle
        call trcmat_2Terms(gammas(i1)%a,GFs(i2)%a(1:nl,1:nl),nl,Tr4G(i))
      enddo
    enddo

    do i1=3,4
      do i2=1,2
        i=i+1
        call trcmat_2Terms(gammas(i1)%a,GFs(i2)%a(n1-nr+1:n1,n1-nr+1:n1),nr,Tr4G(i))
      enddo
    enddo

  endif


  do i=1,4
    call DestroyMatrix(gammas(i),"TSg",io)
  enddo
  do i=1,2
    call DestroyMatrix(GFs(i),"TSg",io)
  enddo

  deallocate(gammas,GFs)


  end subroutine TransmissionSpinComponents_general


  subroutine TransmissionSpinComponents(GfLR,nl,nr,sigmal,sigmar,nspin,ispin,tud,ene)

  implicit none

      
  integer, intent(in) :: nl,nr,nspin,ispin
  DOUBLE COMPLEX,intent(in) :: sigmal(nl,nl)
  DOUBLE COMPLEX,intent(in) :: sigmar(nr,nr)
  DOUBLE COMPLEX, intent(in) :: GfLR(nr,nl)
  double complex, intent(out) :: tud
  double precision, intent(in) :: ene

  DOUBLE COMPLEX, ALLOCATABLE :: gammalDOWN(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gammarDOWN(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: GfLRDOWN(:,:)
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  DOUBLE COMPLEX, ALLOCATABLE, SAVE :: GfLRUP(:,:)
  DOUBLE COMPLEX, ALLOCATABLE, SAVE :: gammalUP(:,:)
  DOUBLE COMPLEX, ALLOCATABLE, SAVE :: gammarUP(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: g0(:,:),gS(:,:)
  double complex tuu,tdd,tdu,t0,ts
  logical gett0s

  integer ii,jj,ind,i


  gett0s=.false.

  if(ispin==1)then
    allocate(GfLRUP(nr,nl))
    allocate(gammalUP(nl,nl))
    allocate(gammarUP(nr,nr))
    GfLRUP=GfLR
    gammalUP=zi*(sigmal-DCONJG(TRANSPOSE(sigmal)))
    gammarUP=zi*(sigmar-DCONJG(TRANSPOSE(sigmar)))
    return
  elseif(ispin==2)then
    allocate(GfLRDOWN(nr,nl))
    allocate(gammalDOWN(nl,nl))
    allocate(gammarDOWN(nr,nr))
    GfLRDOWN=GfLR
    gammalDOWN=zi*(sigmal-DCONJG(TRANSPOSE(sigmal)))
    gammarDOWN=zi*(sigmar-DCONJG(TRANSPOSE(sigmar)))
  endif

  
  if(gett0s)then
    allocate(g0(nr,nl),gS(nr,nl))
    g0=0.5D0 * (GfLRUP+GfLRDOWN)
    gS=0.5D0 * (GfLRUP-GfLRDOWN)
    call trcmat(gammalUP,g0,gammarDOWN,g0,nl,nr,t0)
    call trcmat(gammalUP,gS,gammarDOWN,gS,nl,nr,tS)
    tud=dreal(t0)+ zi * dreal(tS)
    deallocate(g0,gS)
  else
    call trcmat(gammalUP,GfLRUP,gammarDOWN,GfLRDOWN,nl,nr,tud)
  endif

!  write(12347,*)"tsplitr=",ene,dreal(tuu),dreal(tdd),dreal(tud),dreal(tdu)
!  write(12347,*)"tsplitc=",ene,dimag(tuu),dimag(tdd),dimag(tud),dimag(tdu)

  deallocate(GfLRUP)
  deallocate(gammalUP)
  deallocate(gammarUP)
  deallocate(GfLRDOWN)
  deallocate(gammalDOWN)
  deallocate(gammarDOWN)

  end subroutine TransmissionSpinComponents


  subroutine trcmat4(gf1Column,gf2Column,gamma1,g1,gamma2,g2,n1,nl,nr,imin,imax,ene,ntrcchannels,tchannels,psil,psir,outputwfs)
!calculates Trace[gamma1 g1^dagger gamma2 g1] for lead 1
! and
!calculates Trace[gamma2 g2^dagger gamma1 g2] for lead 2

    implicit none

    integer, intent(in) :: nl,nr,n1,imin,imax,ntrcchannels
    DOUBLE COMPLEX,intent(inout) :: gamma1(nl,nl)
    DOUBLE COMPLEX,intent(inout) :: gamma2(nr,nr)
    DOUBLE COMPLEX,intent(in) :: gf1Column(n1,nl)
    DOUBLE COMPLEX,intent(in) :: gf2Column(n1,nr)
    DOUBLE COMPLEX, intent(in) :: g1(nr,nl)
    DOUBLE COMPLEX, intent(in) :: g2(nl,nr)
    double precision, intent(in) :: ene
    double precision, intent(out) :: tchannels(ntrcchannels)
    DOUBLE COMPLEX, intent(out) :: psil(n1,ntrcchannels)
    DOUBLE COMPLEX, intent(out) :: psir(n1,ntrcchannels)
    logical, intent(in) :: outputwfs
     
    DOUBLE COMPLEX, ALLOCATABLE :: tau(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: GammalGdagger(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: tmatl(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: tmatr(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: vr(:,:),vl(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: eigl(:)
    DOUBLE PRECISION, ALLOCATABLE :: eigr(:)
    DOUBLE COMPLEX, ALLOCATABLE :: gamma1po(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: gamma2po(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: sqrtgamma1(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: sqrtgamma2(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: u(:,:),vt(:,:),work(:)
    DOUBLE precision, ALLOCATABLE :: sv(:)
    integer info,i,j,i2,i1,iminl,iminr,imaxl,imaxr
    integer nmin,lwork
    double precision, allocatable :: rwork(:)
    DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0


    ALLOCATE(vl(nl,nl))
    ALLOCATE(eigl(nl))
    call geigenvalueshermitian(gamma1,eigl,nl,vl,info)

    do i=1,nl
      eigl(i)=sqrt(abs(eigl(i)))
    enddo

    ALLOCATE(sqrtgamma1(nl,nl))
    sqrtgamma1=0.0D0
    do i1=1,nl
      do j=1,nl
        do i=1,nl
           sqrtgamma1(i,j)=sqrtgamma1(i,j) + vl(i,i1) * eigl(i1) *  DCONJG(vl(j,i1))  
         enddo
      enddo
    enddo
!    write(12347,*)"deltagsgsg=",maxval(abs(gamma1-matmul(sqrtgamma1,sqrtgamma1)))
    deallocate(vl,eigl)

    ALLOCATE(vr(nr,nr))
    ALLOCATE(eigr(nr))
    call geigenvalueshermitian(gamma2,eigr,nr,vr,info)
!    do i=1,nr
!      if(abs(eigr(i))>1d-13) write(12347,*)"eiggammar=",ene,eigr(i),i
!    enddo

    do i=1,nr
      eigr(i)=sqrt(abs(eigr(i)))
    enddo

    ALLOCATE(sqrtgamma2(nr,nr))
    sqrtgamma2=0.0D0
    do i1=1,nr
      do j=1,nr
        do i=1,nr
           sqrtgamma2(i,j)=sqrtgamma2(i,j) + vr(i,i1) * eigr(i1) *  DCONJG(vr(j,i1))  
         enddo
      enddo
    enddo
    deallocate(vr,eigr)


    ALLOCATE(GammalGdagger(nl,nr))
    CALL ZGEMM('N','C',nl,nr,nl,(1.D0,0.D0),sqrtgamma1,nl,g1,nr,(0.D0,0.D0),GammalGdagger,nl)
    ALLOCATE(tau(nl,nr))
    CALL ZGEMM('N','N',nl,nr,nr,(1.D0,0.D0),GammalGdagger,nl,sqrtgamma2,nr,(0.D0,0.D0),tau,nl)
    deallocate(GammalGdagger)


    nmin=min(nl,nr)
    allocate(u(nl,nmin))
    allocate(vt(1,1))
    allocate(sv(nmin))
    lwork=40 * (nl+nr)
    allocate(work(lwork))
    allocate(rwork(5 * nmin))

    if(outputwfs)then
      CALL ZGESVD( 'S', 'N', nl, nr, tau, nl,sv, u, nl, vt,1, work, lwork, rwork, INFO )
    else
      CALL ZGESVD( 'N', 'N', nl, nr, tau, nl,sv, u, nl, vt,1, work, lwork, rwork, INFO )
    endif
    deallocate(tau,work,rwork)


    iminl=max(1,imin)
    imaxl=min(nmin,imax)
    do i=1,imaxl-iminl+1
      tchannels(i)=sv(i)*sv(i)
    enddo
    deallocate(sv)

    if(.not.outputwfs)then
      deallocate(u,vt)
      deallocate(sqrtgamma2,sqrtgamma1)

      return
    endif


!this part is only executed if wave functions are calculated
    ALLOCATE(tmatl(nl,nmin))
    tmatl=0.0D0
    do j=iminl,imaxl
      do i=1,nl
        do i2=1,nl
          tmatl(i,j)=tmatl(i,j)+sqrtgamma1(i,i2)  * u(i2,j) 
        enddo
      enddo
    enddo
    deallocate(u,vt)
    
!    write(12347,*)"nmin=",nl,nr,nmin
    ALLOCATE(gamma1po(nl,nl))
!    gamma1po=matmul(tmatl,dconjg(transpose(tmatl)))
    CALL ZGEMM('N','C',nl,nl,nmin,(1.D0,0.D0),tmatl,nl,tmatl,nl,(0.D0,0.D0),gamma1po,nl)
    gamma1=gamma1po
    deallocate(gamma1po)

!    write(12347,*)"deltagamma1=",maxval(abs(gamma1-gamma1po))
!    write(12347,*)"deltagsgsg2=",maxval(abs(gamma1-matmul(matmul(sqrtgamma1,u),matmul(dconjg(transpose(u)),sqrtgamma1))))

    ALLOCATE(GammalGdagger(nr,nl))
    CALL ZGEMM('N','C',nr,nl,nr,(1.D0,0.D0),sqrtgamma2,nr,g2,nl,(0.D0,0.D0),GammalGdagger,nr)
    ALLOCATE(tau(nr,nl))
    CALL ZGEMM('N','N',nr,nl,nl,(1.D0,0.D0),GammalGdagger,nr,sqrtgamma1,nl,(0.D0,0.D0),tau,nr)
    deallocate(GammalGdagger)


    nmin=min(nl,nr)
    allocate(u(nr,nmin))
    allocate(vt(1,1))
    allocate(sv(nmin))
    lwork=40 * (nl+nr)
    allocate(work(lwork))
    allocate(rwork(5 * nmin))

    CALL ZGESVD( 'S', 'N', nr, nl, tau, nr,sv, u, nr, vt,1, work, lwork, rwork, INFO )
    deallocate(tau,work,rwork)

    ALLOCATE(tmatr(nr,nmin))

    iminr=max(1,imin)
    imaxr=min(nmin,imax)
    tmatr=0.0D0
    do j=iminr,imaxr
      do i=1,nr
        do i2=1,nr
          tmatr(i,j)=tmatr(i,j)+sqrtgamma2(i,i2)  * u(i2,j)
        enddo
      enddo
    enddo
    deallocate(u,sv,vt)

    deallocate(sqrtgamma1)
    deallocate(sqrtgamma2)
    

!    write(12347,*)"nmin=",nr,nr,nmin
    ALLOCATE(gamma2po(nr,nr))
    CALL ZGEMM('N','C',nr,nr,nmin,(1.D0,0.D0),tmatr,nr,tmatr,nr,(0.D0,0.D0),gamma2po,nr)
!    write(12347,*)"deltagamma2=",maxval(abs(gamma2-gamma2po))
    gamma2=gamma2po
    deallocate(gamma2po)

    psil=matmul(gf1Column,tmatl(:,iminl:imaxl))/sqrt(2.0D0 * pi)
    psir=matmul(gf2Column,tmatr(:,iminr:imaxr))/sqrt(2.0D0 * pi)

    deallocate(tmatl,tmatr)


  end subroutine trcmat4


  subroutine trcmat_2Terms(gamma1,g1,n,t)
    implicit none

    integer, intent(in) :: n
    DOUBLE COMPLEX,intent(in) :: gamma1(n,n)
    DOUBLE COMPLEX,intent(in) :: g1(n,n)
    double complex, intent(out) :: t

    integer i,j

    t=0.0D0
    do j=1,n
      do i=1,n
        t=t+DCONJG(gamma1(i,j))*g1(i,j)
      enddo
    enddo

  end subroutine trcmat_2Terms



  subroutine trcmat(gamma1,g1,gamma2,g2,nl,nr,t)
!calculated Trace[gamma1 g1^dagger gamma2 g1]
    implicit none

    integer, intent(in) :: nl,nr
    DOUBLE COMPLEX,intent(in) :: gamma1(nl,nl)
    DOUBLE COMPLEX,intent(in) :: gamma2(nr,nr)
    DOUBLE COMPLEX, intent(in) :: g1(nr,nl)
    DOUBLE COMPLEX, intent(in) :: g2(nr,nl)
    double complex, intent(out) :: t

    DOUBLE COMPLEX, ALLOCATABLE :: GammarG(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: GGammal(:,:)
    integer i,j


    ALLOCATE(GGammal(nr,nl))
    ALLOCATE(GammarG(nr,nl))

    CALL ZGEMM('N','N',nr,nl,nl,(1.D0,0.D0),g1,nr,gamma1,nl,(0.D0,0.D0),GGammal,nr)
    CALL ZGEMM('N','N',nr,nl,nr,(1.D0,0.D0),gamma2,nr,g2,nr,(0.D0,0.D0),GammarG,nr)

    t=0.0D0
    do j=1,nl
      do i=1,nr
        t=t+DCONJG(GGammal(i,j))*GammarG(i,j)
      enddo
    enddo

    deallocate(GammarG,GGammal)

  end subroutine trcmat


  subroutine trcmat3(gamma1,g1,gamma2,g2,nl,nr,t)
!calculated Trace[gamma1 g1^dagger gamma2 g1]
    implicit none

    integer, intent(in) :: nl,nr
    DOUBLE COMPLEX,intent(in) :: gamma1(nl,nl)
    DOUBLE COMPLEX,intent(in) :: gamma2(nr,nr)
    DOUBLE COMPLEX, intent(in) :: g1(nr,nl)
    DOUBLE COMPLEX, intent(in) :: g2(nr,nl)
    double complex, intent(out) :: t

    DOUBLE COMPLEX, ALLOCATABLE :: GammalGdagger(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: GammarG(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: tmat(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: GammarGT(:,:)
    integer i,j
!    double complex :: t2


    ALLOCATE(GammalGdagger(nl,nr))
    ALLOCATE(GammarG(nr,nl))
!    ALLOCATE(tmat(nl,nl))

    CALL ZGEMM('N','C',nl,nr,nl,(1.D0,0.D0),gamma1,nl,g1,nr,(0.D0,0.D0),GammalGdagger,nl)
    CALL ZGEMM('N','N',nr,nl,nr,(1.D0,0.D0),gamma2,nr,g2,nr,(0.D0,0.D0),GammarG,nr)
!    CALL ZGEMM('N','N',nl,nl,nr,(1.D0,0.D0),GammalGdagger,nl,GammarG,nr,(0.D0,0.D0),tmat,nl)
!    t2=0.0D0
!    do i=1,nl
!      t2=t2+tmat(i,i)
!    enddo
    t=0.0D0
    do j=1,nr
      do i=1,nl
        t=t+GammalGdagger(i,j)*GammarG(j,i)
      enddo
    enddo
!    write(12347,*)"delta_t=",abs(t2-t)


!    deallocate(tmat)
    deallocate(GammarG,GammalGdagger)

  end subroutine trcmat3

  subroutine trcmat2(gamma1,g1,gamma2,g2,nl,nr,t)
!calculated Trace[gamma1 g1^dagger gamma2 g1]
    implicit none

    integer, intent(in) :: nl,nr
    DOUBLE COMPLEX,intent(in) :: gamma1(nl,nl)
    DOUBLE COMPLEX,intent(in) :: gamma2(nr,nr)
    DOUBLE COMPLEX, intent(in) :: g1(nr,nl)
    DOUBLE COMPLEX, intent(in) :: g2(nr,nl)
    double complex, intent(out) :: t

    DOUBLE COMPLEX, ALLOCATABLE :: GammalGdagger(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: GammarG(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: tmat(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: GammarGT(:,:)
    integer i,j
!    double complex :: t2


    ALLOCATE(GammalGdagger(nl,nr))
    ALLOCATE(GammarG(nr,nl))
!    ALLOCATE(tmat(nl,nl))

    CALL ZGEMM('N','C',nl,nr,nl,(1.D0,0.D0),gamma1,nl,g1,nr,(0.D0,0.D0),GammalGdagger,nl)
    CALL ZGEMM('N','N',nr,nl,nr,(1.D0,0.D0),gamma2,nr,g2,nr,(0.D0,0.D0),GammarG,nr)
!    CALL ZGEMM('N','N',nl,nl,nr,(1.D0,0.D0),GammalGdagger,nl,GammarG,nr,(0.D0,0.D0),tmat,nl)
!    t2=0.0D0
!    do i=1,nl
!      t2=t2+tmat(i,i)
!    enddo
    ALLOCATE(GammarGT(nl,nr))

    do i=1,nl
      do j=1,nr
        GammarGT(i,j)=GammarG(j,i)
      enddo
    enddo

    t=0.0D0
    do j=1,nr
      do i=1,nl
        t=t+GammalGdagger(i,j)*GammarGT(i,j)
      enddo
    enddo
!    write(12347,*)"delta_t=",abs(t2-t)


!    deallocate(tmat)
    deallocate(GammarG,GammalGdagger,GammarGT)

  end subroutine trcmat2



  subroutine TransmissionChannelsDecomposition(gfLColumn,gfRColumn,GfLR,GfRL,n1,nl,nr,sigmal,sigmar,nspin,ispin,ene,ef,gamma1op,gamma2op,imin,imax,ntrcchannels,tchannels,ik,ie_global,nk,iv,slabel,kpoint,outputwfs)

  use mTypes
  implicit none

      
  logical, intent(in) :: outputwfs
  integer, intent(in) :: n1,nl,nr,nspin,ispin,imin,imax,ntrcchannels,ik,iv,nk,ie_global
  CHARACTER, intent(in) :: slabel*20
  DOUBLE precision,intent(in) :: kpoint(2)
  DOUBLE COMPLEX,intent(in) :: sigmal(nl,nl)
  DOUBLE COMPLEX,intent(in) :: sigmar(nr,nr)
  DOUBLE COMPLEX,intent(in) :: gfLColumn(n1,nl)
  DOUBLE COMPLEX,intent(in) :: gfRColumn(n1,nr)
  DOUBLE COMPLEX, intent(in) :: GfLR(nr,nl)
  DOUBLE COMPLEX, intent(in) :: GfRL(nl,nr)
  double precision, intent(in) :: ene,ef
  DOUBLE COMPLEX,intent(out) :: gamma1op(nl,nl)
  DOUBLE COMPLEX,intent(out) :: gamma2op(nr,nr)
  double precision, intent(out) :: tchannels(ntrcchannels)

  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  DOUBLE COMPLEX, allocatable :: psil(:,:)
  DOUBLE COMPLEX, allocatable :: psir(:,:)
  double precision e

  integer ii,jj,ind,i


  if(outputwfs)then
    allocate(psil(n1,ntrcchannels))
    allocate(psir(n1,ntrcchannels))
  else
    allocate(psil(1,1))
    allocate(psir(1,1))
  endif

  gamma1op=zi*(sigmal-DCONJG(TRANSPOSE(sigmal)))
  gamma2op=zi*(sigmar-DCONJG(TRANSPOSE(sigmar)))
  
  call trcmat4(gfLColumn,gfRColumn,gamma1op,GfLR,gamma2op,GfRL,n1,nl,nr,imin,imax,ene,ntrcchannels,tchannels,psil,psir,outputwfs)

  if(outputwfs)then
!    write(12347,*)"maxval_psi12",maxval(abs(psil-psil2))
    e=(ene-ef)*13.6056981D0
    call  writewfs(iv,slabel,ik,ie_global,e,nk,n1,ntrcchannels,psil,ispin,tchannels,kpoint,nspin,"_L")
    call  writewfs(iv,slabel,ik,ie_global,e,nk,n1,ntrcchannels,psir,ispin,tchannels,kpoint,nspin,"_R")
  endif

  deallocate(psil)
  deallocate(psir)

  end subroutine TransmissionChannelsDecomposition




  SUBROUTINE writewfs(iv,slabel,ik,ie_global,e,nk,N1,nbs,psi,ispin,eig,kpoint,nspin,nam)
  
    use mMPI_NEGF

    implicit none

    integer, intent(in) :: nbs,ispin,ik,n1,iv,nk,nspin,ie_global
    double precision, intent(in) :: e
    DOUBLE COMPLEX, intent(in) :: psi(N1,nbs)
    DOUBLE precision,intent(in) :: eig(nbs)
    DOUBLE precision,intent(in) :: kpoint(2)
    CHARACTER(LEN=*), intent(in) :: nam

    character*30  charnode
    character fname*43
    DOUBLE COMPLEX :: psin(N1),psint(N1)
    double complex enec
    integer i,j
    double precision kpointmod(3)
    CHARACTER :: chivv*8
    CHARACTER :: slabel*20, paste*100, filetr*100
    DOUBLE PRECISION :: ef,v
    CHARACTER :: iec*50
    CHARACTER :: ec*50
    CHARACTER :: ikc*24
    character pasbias*55
    integer iufiletr,is


!    call io_assign(iufiletr)

    write( iec, '(i50)')   ie_global
    write( ec, '(f12.3)')  e
    write( ikc, '(i24 )' ) ik
    write( chivv, '(i7 )' ) iv
    chivv = pasbias(chivv, '.')
    filetr=TRIM(ADJUSTL(TRIM(ADJUSTL(slabel))//"-ie"//TRIM(ADJUSTL(iec))//"-e_"//TRIM(ADJUSTL(ec))//"eV-k_"//TRIM(ADJUSTL(ikc))// TRIM(ADJUSTL(nam))))
    if(ispin==1)then
      filetr = paste(filetr,'_UP')
    elseif(ispin==2)then
      filetr = paste(filetr,'_DOWN')
    endif
    filetr = paste(filetr,'.WFS')
    filetr = paste( chivv, filetr )

    fname=TRIM(ADJUSTL(filetr))


    kpointmod(1)=kpoint(1)
    kpointmod(2)=kpoint(2)
    kpointmod(3)=0.0D0

    iufiletr=12348
    OPEN(UNIT=iufiletr,FILE=fname, form='unformatted',status='unknown')

     write(iufiletr) 1
     write(iufiletr) nspin
     write(iufiletr) n1

     endfile (iufiletr)
     backspace (iufiletr)
     close (iufiletr)


      do is=1,nspin
        open(iufiletr, file=fname, form='unformatted', position='append', status='old' )
        write(iufiletr) ik,kpointmod(1),kpointmod(2),kpointmod(3)
        write(iufiletr) is
        write(iufiletr) nbs

        do i=1,nbs

        psin=psi(:,i)
        psint=psi(:,i)
        call NormalizePsi(n1,psin)

        enec=eig(i)
        call writephibin(N1,psin,enec,i)

        enddo
      enddo

    close(iufiletr)


    write(12347,*)"start of wf output",ispin,ik,nbs
    do is=1,nspin
      do i=1,nbs
          write(12347,*)"energy=",eig(i)
          write(12347,*)
          write(12347,'(a72)') ' ***********************************************************************'
          write(12347,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik
          write(12347,'(a22,2x,i6)') 'Spin component = ',ispin
          write(12347,'(a22,2x,i6)') 'Num. wavefunctions = ',i

          psin=psi(:,i)
          psint=psi(:,i)
          call NormalizePsi(n1,psin)

          enec=eig(i)
          call writephiascii(N1,psin,psint,enec,i)

      enddo
    enddo
    write(12347,*)"end of wf output",ispin,ik



  end SUBROUTINE writewfs


  subroutine NormalizePsi(n,psi)

  implicit none

  integer, intent(in) :: n
  double complex, intent (inout):: psi(n)
  integer i
  double complex norm,maxpsi

  norm=0.0D0
  maxpsi=0.0D0
  do i=1,n
    norm=norm+abs(psi(i))**2
    if(abs(psi(i))>abs(maxpsi))then
      maxpsi=psi(i)
    endif
  enddo
  norm=sqrt(norm)

  if(norm.ne.0)then
    psi=psi/norm
    maxpsi=maxpsi/norm
    psi=(psi/maxpsi)* abs(maxpsi)
  endif

  

   
  end subroutine NormalizePsi


  subroutine ShiftRightSelfenergy(n1,nl,nr,nsplit,KC,KmC,KTip,SigmaRN,nrN,work,ispin)

  implicit none

      
  integer, intent(in) :: n1,nl,nr,nsplit,ispin
  integer, intent(in) :: nrN
  DOUBLE COMPLEX, intent(out) :: SigmaRN(nrN,nrN,2)
!  DOUBLE COMPLEX :: SigmaRN2(nrN,nrN)
!  DOUBLE COMPLEX :: GfLN(nsplit-1,nsplit-1)
  DOUBLE COMPLEX, intent(inout) :: work(n1*n1)

  double complex,intent(inout) ::Ktip(n1-nsplit+1,n1-nsplit+1,2)
  double complex,intent(inout) ::KC(nrN,n1-nsplit+1,2)
  double complex,intent(inout) ::KmC(n1-nsplit+1,nrN,2)

  integer nt,nb,info,i,is
  integer, allocatable :: IPIV(:)

  nt=n1-nsplit+1
  nb=nsplit-1

  allocate(IPIV(nt))
  CALL ZGETRF(nt,nt,Ktip(:,:,ispin),nt,IPIV,   INFO)
  CALL ZGETRI(nt,Ktip(:,:,ispin),  nt,IPIV,WORK,n1**2,INFO)
  deallocate(IPIV)

  if(ispin==1)return


  if(.true.)then
    KC(:,:,1)=0.5D0*(KC(:,:,1)+KC(:,:,2))
    KmC(:,:,1)=0.5D0*(KmC(:,:,1)+KmC(:,:,2))
    KC(:,:,2)=KC(:,:,1)
    KmC(:,:,2)=KmC(:,:,1)
  endif

  do is=1,2
    call SeneMatMul(nrN,nt,SigmaRN(:,:,is),KC(:,:,is),Ktip(:,:,is),KmC(:,:,is))
  enddo
!  if(.false.)then
!    GfLN=kmat(1:nb,1:nb)
!    GfLN(nb-nrN+1:nb,nb-nrN+1:nb)=GfLN(nb-nrN+1:nb,nb-nrN+1:nb)-SigmaRN
!    allocate(IPIV(nb))
!    CALL ZGETRF(nb,nb,GfLN,nb,IPIV,   INFO)
!    CALL ZGETRI(nb,GfLN,  nb,IPIV,WORK,n1**2,INFO)
!    deallocate(IPIV)
!    write(12347,*)"deltag_b= ",maxval(abs(Gf(1:nb,1:nb)-GfLN))
!  endif


  end subroutine ShiftRightSelfenergy



  subroutine SeneMatMul(n1,n2,sigma,kc,k2,kmc)

    implicit none

    integer, intent(in) :: n1,n2
    DOUBLE COMPLEX,intent(in) :: kc(n1,n2)
    DOUBLE COMPLEX,intent(in) :: kmc(n2,n1)
    DOUBLE COMPLEX,intent(in) :: k2(n2,n2)
    double complex, intent(out) :: sigma(n1,n1)

    DOUBLE COMPLEX, ALLOCATABLE :: GammalGdagger(:,:)

    ALLOCATE(GammalGdagger(n1,n2))

    CALL ZGEMM('N','N',n1,n2,n2,(1.D0,0.D0),kc,n1,k2,n2,(0.D0,0.D0),GammalGdagger,n1)
    CALL ZGEMM('N','N',n1,n1,n2,(1.D0,0.D0),GammalGdagger,n1,kmc,n2,(0.D0,0.D0),sigma,n1)

    deallocate(GammalGdagger)

  end subroutine SeneMatMul





  subroutine FindnrN(KmC,nt,nb,nrN)

  implicit NONE

  integer, intent (in) :: nb,nt
  integer, intent (out) :: nrN
  double complex, intent (in) :: KmC(nt,nb)
  double precision tol

  integer i1,i2

  nrN=0
  tol=0.0D0

  outloop: do i2=1,nb
    do i1=1,nt
!      write(12347,*)"KmC=", i1,i2,abs(KmC(i1,i2))
      if(abs(KmC(i1,i2))>tol)then
        nrN=nb-i2+1
        exit outloop
      endif
    enddo
  enddo outloop
!  write(12347,*)"nrN2=",nb,i2,nrN

  end subroutine FindnrN

  subroutine Tmatrix(e,ef,nl,nr,gfRL,gfLL,gfLR,gfRR,sigmal,sigmar,ispin,nspin)

  use negfmod, only: ikpmod,kpointmod,TransmissionMatrixWFS,EM_NSPINBlocks, ndivxyNL
  use mSigmaFourier, only : ExpandSstates
  use mComputeULR, only : SortSstates,PhiS,FreePhiS
  use sigma, only: h0_l,h1_l,s0_l,s1_l,h0_r,h1_r,s0_r,s1_r
  use mTransmissionDecompositionRemove

  implicit none

  integer, intent(in) :: nl,nr,ispin,nspin
  double complex, intent(in) :: e
  double precision, intent(in) :: ef
  double complex, intent(in) :: gfRL(nr,nl)
  double complex, intent(in) :: gfLL(nl,nl)
  double complex, intent(in) :: gfLR(nl,nr)
  double complex, intent(in) :: gfRR(nr,nr)
  double complex, intent(in) :: sigmal(nl,nl),sigmar(nr,nr)

  double complex, ALLOCATABLE :: tm(:,:),tmL(:),tmR(:)
  double complex, ALLOCATABLE :: tmRL_R(:),tmRl_L(:)
  double complex, ALLOCATABLE :: rmRR_R(:),rmRR_Rlg(:)

  double complex, ALLOCATABLE :: tmRL(:,:)
  double complex, ALLOCATABLE :: rmRR(:,:),rm2RR(:,:),rm1RR(:,:)
  double complex, ALLOCATABLE :: rm(:,:),rmL(:),rmLlg(:),rm2(:,:),rm1(:,:)
  double complex, ALLOCATABLE :: SSdagger(:,:)
  double complex, allocatable :: gammal(:,:),gammar(:,:)
  double complex, allocatable :: velocityl(:,:),velocityr(:,:)
  double complex, allocatable :: velocityLrg(:),velocityRrg(:),velocityLlg(:),velocityRlg(:)
  double complex, allocatable :: mat1(:,:),mat2(:,:)
  double complex :: ttotal,rtotal,ttotalR,rtotallg
  double precision alpha, kappa,alpha2
  integer i1,i2,i3
  double complex :: ene

  ene=e-ef

  call ExpandSstates(PhiS(1),'L',PhiS(1)%ndivxy,PhiS(1)%n,PhiS(1)%ns,PhiS(1)%e,PhiS(1)%kxy)
  call FreePhiS(PhiS(1))
  call ExpandSstates(PhiS(2),'R',PhiS(2)%ndivxy,PhiS(2)%n,PhiS(2)%ns,PhiS(2)%e,PhiS(2)%kxy)
  call FreePhiS(PhiS(2))


  if(nTstatesL>0)then
    call SortSstates(e,nl,nTstatesL,kALrg,uAL,uAtildeL,VgLrg,MuAL,EM_NSPINBlocks)
  endif
  if(nTstatesLlg> 0 ) then
    call SortSstates(e,nl,nTstatesLlg,kLlg,uLlg,utildeLlg,VgLlg,MuLlg,EM_NSPINBlocks)
  endif
  if(nTstatesR>0)then
    call SortSstates(e,nr,nTstatesR,kARlg,uARlg,uAtildeRlg,VgRlg,MuARlg,EM_NSPINBlocks)
  endif
  if(nTstatesRlg>0)then
    call SortSstates(e,nr,nTstatesRlg,kRrg,uR,utildeR,VgRrg,MuR,EM_NSPINBlocks)
  endif

  call ComputeVg('L',nl,e,S0_L,S1_L,DCONJG(TRANSPOSE(H1_L(:,:,ispin))) - e*DCONJG(TRANSPOSE(S1_L)),H1_L(:,:,ispin) - e*S1_L)
  call ComputeVg('R',nr,e,S0_R,S1_R,DCONJG(TRANSPOSE(H1_R(:,:,ispin))) - e*DCONJG(TRANSPOSE(S1_R)),H1_R(:,:,ispin) - e*S1_R)


  if(TransmissionMatrixWFS)then
    do i1=1,nTstatesL
      write(12347,*)"k_L_rg=",13.6056981D0 * dreal(ene),i1,dreal(kALrg(i1)),dimag(kALrg(i1)),ispin,ikpmod
      do i2=1,nl
        write(12347,*)"phiA_L_rg=",13.6056981D0 * dreal(ene),i1,i2,dreal(uAL(i2,i1)),dimag(uAL(i2,i1)),sign(abs(uAL(i2,i1)),dreal(uAL(i2,i1))),ikpmod
      enddo

      do i2=1,nl
        write(12347,*)"phiAt_L_rg=",13.6056981D0 * dreal(ene),i1,i2,dreal(uAtildeL(i2,i1)),dimag(uAtildeL(i2,i1)),sign(abs(uAtildeL(i2,i1)),dreal(uAtildeL(i2,i1))),ikpmod
      enddo

    enddo

    do i1=1,nTstatesR
      write(12347,*)"k_R_rg=",13.6056981D0 * dreal(ene),i1,dreal(kRrg(i1)),dimag(kRrg(i1)),ispin,ikpmod
      do i2=1,nr
        write(12347,*)"phi_R_rg=",13.6056981D0 * dreal(ene),i1,i2,dreal(uR(i2,i1)),dimag(uR(i2,i1)),sign(abs(uR(i2,i1)),dreal(uR(i2,i1))),ispin,ikpmod
      enddo
      do i2=1,nr
        write(12347,*)"phit_R_rg=",13.6056981D0 * dreal(ene),i1,i2,dreal(utildeR(i2,i1)),dimag(utildeR(i2,i1)),sign(abs(utildeR(i2,i1)),dreal(utildeR(i2,i1))),ispin,ikpmod
      enddo

    enddo
    do i1=1,nTstatesLlg
      write(12347,*)"k_L_lg=",13.6056981D0 * dreal(ene),i1,dreal(kLlg(i1)),dimag(kLlg(i1)),ispin,ikpmod
      do i2=1,nl
        write(12347,*)"phi_L_lg=",13.6056981D0 * dreal(ene),i1,i2,dreal(uLlg(i2,i1)),dimag(uLlg(i2,i1)),sign(abs(uLlg(i2,i1)),dreal(uLlg(i2,i1))),ispin,ikpmod
      enddo
      do i2=1,nl
        write(12347,*)"phit_L_lg=",13.6056981D0 * dreal(ene),i1,i2,dreal(utildeLlg(i2,i1)),dimag(utildeLlg(i2,i1)),sign(abs(utildeLlg(i2,i1)),dreal(utildeLlg(i2,i1))),ispin,ikpmod
      enddo

    enddo

  endif

  allocate(velocityLrg(nTstatesL),velocityRrg(nTstatesR),velocityLlg(nTstatesLlg),velocityRlg(nTstatesRlg))


  allocate(gammal(nl,nl))
  gammal=(0.0D0,1.0D0)*(sigmal-DCONJG(transpose(sigmal)))

  allocate(mat1(nl,nTstatesL))
  mat1=matmul(gammal,uAL)

  velocityLrg=0.0D0
  do i1=1,nTstatesL
    do i2=1,nl
      velocityLrg(i1)=velocityLrg(i1)+DCONJG(uAL(i2,i1))*mat1(i2,i1)
    enddo
  enddo
  deallocate(mat1)

  allocate(mat1(nl,nTstatesLlg))
  mat1=matmul(gammal,uLlg)

  velocityLlg=0.0D0
  do i1=1,nTstatesLlg
    do i2=1,nl
      velocityLlg(i1)=velocityLlg(i1)+DCONJG(uLlg(i2,i1))*mat1(i2,i1)
    enddo
  enddo
  deallocate(mat1)
!note that velocityLlg is actually the absolute value of the velocity for left going states

  deallocate(gammal)

  allocate(mat2(nr,nTstatesR))
  allocate(gammar(nr,nr))
  gammar=(0.0D0,1.0D0)*(sigmar-DCONJG(transpose(sigmar)))
  mat2=matmul(gammar,uR)

  velocityRrg=0.0D0
  do i1=1,nTstatesR
    do i2=1,nr
      velocityRrg(i1)=velocityRrg(i1)+DCONJG(uR(i2,i1))*mat2(i2,i1)
    enddo
  enddo
  deallocate(mat2)


  allocate(mat2(nr,nTstatesRlg))
  mat2=matmul(gammar,uARlg)

  velocityRlg=0.0D0
  do i1=1,nTstatesRlg
    do i2=1,nr
      velocityRlg(i1)=velocityRlg(i1)+DCONJG(uARlg(i2,i1))*mat2(i2,i1)
    enddo
  enddo
!note that velocityRlg is actually the absolute value of the velocity for left going states
  deallocate(mat2)

  deallocate(gammar)

  do i1=1,nTstatesL
    write(12347,*)"VgL_rg=",13.6056981D0 * dreal(ene),i1,dreal(VgLrg(i1)),dreal(velocityLrg(i1)),ispin,ikpmod
  enddo
  do i1=1,nTstatesLlg
    write(12347,*)"VgL_lg=",13.6056981D0 * dreal(ene),i1,dreal(-VgLlg(i1)),dreal(velocityLlg(i1)),ispin,ikpmod
  enddo
  do i1=1,nTstatesR
    write(12347,*)"VgR_rg=",13.6056981D0 * dreal(ene),i1,dreal(VgRrg(i1)),dreal(velocityRrg(i1)),ispin,ikpmod
  enddo
  do i1=1,nTstatesRlg
    write(12347,*)"VgR_lg=",13.6056981D0 * dreal(ene),i1,dreal(-VgRlg(i1)),dreal(velocityRlg(i1)),ispin,ikpmod
  enddo



  do i1=1,nTstatesL
    if(abs(VgLrg(i1)-velocityLrg(i1))>1.0d-7) write(12347,*)"Delta_VgL_rg=",13.6056981D0 * dreal(ene),i1,dreal(VgLrg(i1)-velocityLrg(i1)),dimag(VgLrg(i1)-velocityLrg(i1)),ispin,ikpmod
  enddo
  do i1=1,nTstatesLlg
    if(abs(-VgLlg(i1)-velocityLlg(i1))>1.0d-7) write(12347,*)"Delta_VgL_lg=",13.6056981D0 * dreal(ene),i1,dreal(-VgLlg(i1)-velocityLlg(i1)),dimag(-VgLlg(i1)-velocityLlg(i1)),ispin,ikpmod
  enddo
  do i1=1,nTstatesR
    if(abs(VgRrg(i1)-velocityRrg(i1))>1.0d-7) write(12347,*)"Delta_VgR_rg=",13.6056981D0 * dreal(ene),i1,dreal(VgRrg(i1)-velocityRrg(i1)),dimag(VgRrg(i1)-velocityRrg(i1)),ispin,ikpmod
  enddo
  do i1=1,nTstatesRlg
    if(abs(-VgRlg(i1)-velocityRlg(i1))>1.0d-7) write(12347,*)"Delta_VgR_lg=",13.6056981D0 * dreal(ene),i1,dreal(-VgRlg(i1)-velocityRlg(i1)),dimag(-VgRlg(i1)-velocityRlg(i1)),ispin,ikpmod
  enddo

!start T and R for states incoming from the left

  allocate(tm(nTstatesR,nTstatesL))
  allocate(rm(nTstatesLlg,nTstatesL))
  allocate(rm1(nTstatesLlg,nTstatesL))
  allocate(rm2(nTstatesLlg,nTstatesL))
!tm and rm can be multiplied by an arbitrary phase, but it needs to be the same for both
  tm=(0.0D0,1.0D0)*matmul(matmul(DCONJG(TRANSPOSE(utildeR)),gfRL),uAtildeL)
  rm1=(0.0D0,1.0D0)*matmul(matmul(DCONJG(TRANSPOSE(utildeLlg)),gfLL),uAtildeL)
  rm2=matmul(DCONJG(TRANSPOSE(utildeLlg)),uAL)

  do i2=1,nTstatesL
    do i1=1,nTstatesR
      tm(i1,i2)=sqrt(velocityRrg(i1)*velocityLrg(i2)) * tm(i1,i2)
    enddo
  enddo

  do i2=1,nTstatesL
    do i1=1,nTstatesLlg
      rm1(i1,i2)=sqrt(velocityLlg(i1)*velocityLrg(i2)) * rm1(i1,i2)
      rm2(i1,i2)=sqrt((velocityLlg(i1)/velocityLrg(i2))) * rm2(i1,i2)
    enddo
  enddo

  rm=rm1-rm2


  if(.true.)then

    do i2=1,nTstatesL
      do i1=1,nTstatesR
!transmission matrix output
!optional : convert k to 1/angstrom
        call alphaofz2(alpha,kappa,tm(i1,i2))

        if(EM_NSPINBlocks.ne.4)then
          write(12347,*)"tkLR_iLout_jLin=",13.6056981D0 * dreal(ene),i1,i2,abs(tm(i1,i2)),alpha,dreal(kRrg(i1))-dreal(kALrg(i2)),dreal(kRrg(i1)),dreal(kALrg(i2)),dreal(VgRrg(i1)),dreal(VgLrg(i2)),ispin,ikpmod
        else
          write(12347,*)"tkLR_iLout_jLin=",13.6056981D0 * dreal(ene),i1,i2,abs(tm(i1,i2)),alpha,dreal(kRrg(i1))-dreal(kALrg(i2)),dreal(kRrg(i1)),dreal(kALrg(i2)),dreal(VgRrg(i1)),dreal(VgLrg(i2)),ispin,ikpmod,MuR(1,i1),MuR(2,i1),MuR(3,i1),MuR(4,i1),MuAL(1,i2),MuAL(2,i2),MuAL(3,i2),MuAL(4,i2)
        endif
      enddo
    enddo

    do i2=1,nTstatesL
      do i1=1,nTstatesLlg
!reflection matrix output
!optional : convert k to 1/angstrom
        call alphaofz2(alpha,kappa,rm(i1,i2))
        call alphaofz2(alpha2,kappa,ALrglg(i2,i1))
        if(EM_NSPINBlocks.ne.4)then
          write(12347,*)"rkLL_iLout_jLin=",13.6056981D0 * dreal(ene),i1,i2,abs(rm(i1,i2)),alpha,dreal(kLlg(i1))-dreal(kALrg(i2)),dreal(kLlg(i1)),dreal(kALrg(i2)),dreal(VgLlg(i1)),dreal(VgLrg(i2)),abs(ALrglg(i2,i1)),alpha2,ispin,ikpmod
        else
          write(12347,*)"rkLL_iLout_jLin=",13.6056981D0 * dreal(ene),i1,i2,abs(rm(i1,i2)),alpha,dreal(kLlg(i1))-dreal(kALrg(i2)),dreal(kLlg(i1)),dreal(kALrg(i2)),dreal(VgLlg(i1)),dreal(VgLrg(i2)),abs(ALrglg(i2,i1)),alpha2,ispin,ikpmod,MuLlg(1,i1),MuLlg(2,i1),MuLlg(3,i1),MuLlg(4,i1),MuAL(1,i2),MuAL(2,i2),MuAL(3,i2),MuAL(4,i2)
        endif

      enddo
    enddo
    if(allocated(ALrglg))deallocate(ALrglg) 

  endif


!end T and R for states incoming from the left

!start T and R for states incoming from the right

  allocate(tmRL(nTstatesLlg,nTstatesRlg))
  allocate(rmRR(nTstatesR,nTstatesRlg))
  allocate(rm1RR(nTstatesR,nTstatesRlg))
  allocate(rm2RR(nTstatesR,nTstatesRlg))

!tmRL and rmRR can be multiplied by an arbitrary phase, but it needs to be the same for both
  tmRL=(0.0D0,1.0D0)*matmul(matmul(DCONJG(TRANSPOSE(utildeLlg)),gfLR),uAtildeRlg)
  rm1RR=(0.0D0,1.0D0)*matmul(matmul(DCONJG(TRANSPOSE(utildeR)),gfRR),uAtildeRlg)
  rm2RR=matmul(DCONJG(TRANSPOSE(utildeR)),uARlg)

  do i2=1,nTstatesRlg
    do i1=1,nTstatesLlg
      tmRL(i1,i2)=sqrt(velocityLlg(i1)*velocityRlg(i2)) * tmRL(i1,i2)
    enddo
  enddo

  do i2=1,nTstatesRlg
    do i1=1,nTstatesR
      rm1RR(i1,i2)=sqrt(velocityRrg(i1)*velocityRlg(i2)) * rm1RR(i1,i2)
      rm2RR(i1,i2)=sqrt((velocityRrg(i1)/velocityRlg(i2))) * rm2RR(i1,i2)
    enddo
  enddo

  rmRR=rm1RR-rm2RR

  if(.true.)then

    do i2=1,nTstatesRlg
      do i1=1,nTstatesLlg
        call alphaofz2(alpha,kappa,tmRL(i1,i2))

        if(EM_NSPINBlocks.ne.4)then
          write(12347,*)"tkRL_iLout_jRin=",13.6056981D0 * dreal(ene),i1,i2,abs(tmRL(i1,i2)),alpha,dreal(kLlg(i1))-dreal(kARlg(i2)),dreal(kLlg(i1)),dreal(kARlg(i2)),dreal(VgLlg(i1)),dreal(VgRlg(i2)),ispin,ikpmod
        else
          write(12347,*)"tkRL_iLout_jRin=",13.6056981D0 * dreal(ene),i1,i2,abs(tmRL(i1,i2)),alpha,dreal(kLlg(i1))-dreal(kARlg(i2)),dreal(kLlg(i1)),dreal(kARlg(i2)),dreal(VgLlg(i1)),dreal(VgRlg(i2)),ispin,ikpmod,MuLlg(1,i1),MuLlg(2,i1),MuLlg(3,i1),MuLlg(4,i1),MuARlg(1,i2),MuARlg(2,i2),MuARlg(3,i2),MuARlg(4,i2)
        endif
      enddo
    enddo

    do i2=1,nTstatesRlg
      do i1=1,nTstatesR
        call alphaofz2(alpha,kappa,rmRR(i1,i2))
        call alphaofz2(alpha2,kappa,ARlgrg(i2,i1))
        if(EM_NSPINBlocks.ne.4)then
          write(12347,*)"rkRR_iRout_jRin=",13.6056981D0 * dreal(ene),i1,i2,abs(rmRR(i1,i2)),alpha,dreal(kRrg(i1))-dreal(kARlg(i2)),dreal(kRrg(i1)),dreal(kARlg(i2)),dreal(VgRrg(i1)),dreal(VgRlg(i2)),abs(ARlgrg(i2,i1)),alpha2,ispin,ikpmod
        else
          write(12347,*)"rkRR_iRout_jRin=",13.6056981D0 * dreal(ene),i1,i2,abs(rmRR(i1,i2)),alpha,dreal(kRrg(i1))-dreal(kARlg(i2)),dreal(kRrg(i1)),dreal(kARlg(i2)),dreal(VgRrg(i1)),dreal(VgRlg(i2)),abs(ARlgrg(i2,i1)),alpha2,ispin,ikpmod,MuR(1,i1),MuR(2,i1),MuR(3,i1),MuR(4,i1),MuARlg(1,i2),MuARlg(2,i2),MuARlg(3,i2),MuARlg(4,i2)
        endif

      enddo
    enddo
    if(allocated(ARlgrg))deallocate(ARlgrg)

  endif

!end T and R for states incoming from the right

  deallocate(velocityLrg,velocityRrg,velocityLlg,velocityRlg)
  deallocate(VgLrg,VgLlg,VgRrg,VgRlg)

  deallocate(MuAL,MuR,MuLlg,MuARlg)

  deallocate(rm1,rm2)
  deallocate(utildeR, uAtildeL,utildeLlg,uAtildeRlg)
  deallocate(uAL,uR,uLlg,uARlg)
  deallocate(kRrg,kALrg,kLlg,kARlg)


!start calculation and output of partial sums of T matrix squared and total transmission


  allocate(tmL(nTstatesL),tmR(nTstatesR))
  allocate(rmL(nTstatesL),rmLlg(nTstatesLlg))
  call TR_Total("Lrg",nTstatesL, nTstatesR, nTstatesLlg, tm,rm,tmL,tmR,rmL,rmLlg,ttotal,rtotal,ttotalR,rtotallg,ene,ispin,ikpmod)
  allocate(tmRL_R(nTstatesRlg),tmRL_L(nTstatesLlg))
  allocate(rmRR_R(nTstatesR),rmRR_Rlg(nTstatesRlg))
  call TR_Total("Rlg",nTstatesRlg, nTstatesLlg, nTstatesR, tmRL,rmRR,tmRL_R,tmRL_L,rmRR_R,rmRR_Rlg,ttotal,rtotal,ttotalR,rtotallg,ene,ispin,ikpmod)

  nTstatesS=nTstatesLlg+nTstatesR
  allocate(Smat(nTstatesS,nTstatesS))
  Smat(1:nTstatesLlg,1:nTstatesLlg)=rm
  Smat(1:nTstatesLlg,nTstatesLlg+1:nTstatesS)=tmRL
  Smat(nTstatesLlg+1:nTstatesS,1:nTstatesLlg)=tm
  Smat(nTstatesLlg+1:nTstatesS,nTstatesLlg+1:nTstatesS)=rmRR 

  if(.true.)then
    allocate(SSdagger(nTstatesS,nTstatesS))
    SSdagger=matmul(Smat,DCONJG(TRANSPOSE(Smat)))
    do i1=1,nTstatesS
      SSdagger(i1,i1)=SSdagger(i1,i1)-1.0D0
    enddo
    call writemat8(dreal(ene),SSdagger,nTstatesS,nTstatesS,1.0d-8,"SSdagger-1")
    deallocate(SSdagger)
  endif

  if(.true.)then
!No other operations are currenlty performed on Smat
!xxx: this was calling writemat8b, changed it to 8; 8b was printing transposed matrix
    call writemat8(dreal(ene),Smat,nTstatesS,nTstatesS,0.0d-8,"Smat")
  endif
  deallocate(Smat)

  deallocate(tm)
  deallocate(rm)
  deallocate(tmL,tmR,rmL,rmLlg)

  deallocate(tmRL)
  deallocate(rmRR)
  deallocate(tmRL_R,tmRL_L,rmRR_R,rmRR_Rlg)

  end subroutine Tmatrix


  subroutine alphaofz2(alpha,kappa,zev)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************


    DOUBLE PRECISION   alpha,kappa,lz,zr,zi,alpha0,pi
    DOUBLE COMPLEX     zev

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


  subroutine TR_Total(prefix,n_in, n_out, n_back, tm,rm,tm_in,tm_out,rm_in,rm_back,ttotal_in,rtotal_in,ttotal_out,rtotal_back,ene,ispin,ikpmod)

  integer, intent(in)         :: n_in, n_back, n_out,ispin,ikpmod
  double complex, intent(in)  :: tm(n_out,n_in)
  double complex, intent(in)  :: rm(n_back,n_in)
  double complex, intent(in)  :: ene
  character(LEN=*)            :: prefix

  double complex, intent(out) :: tm_in(n_in),tm_out(n_out)
  double complex, intent(out) :: rm_in(n_in),rm_back(n_back)
  double complex, intent(out) :: ttotal_in,rtotal_in,ttotal_out,rtotal_back

  integer i1,i2

  tm_in=0.0D0
  do i2=1,n_in
    do i1=1,n_out
      tm_in(i2)=tm_in(i2)+abs(tm(i1,i2))**2
    enddo
  enddo

!!  tm_out=0.0D0
!!  do i2=1,n_in
!!    do i1=1,n_out
!!      tm_out(i1)=tm_out(i1)+abs(tm(i1,i2))**2
!!    enddo
!!  enddo

  rm_in=0.0D0
  do i2=1,n_in
    do i1=1,n_back
      rm_in(i2)=rm_in(i2)+abs(rm(i1,i2))**2
    enddo
  enddo

!!  rm_back=0.0D0
!!  do i2=1,n_in
!!    do i1=1,n_back
!!      rm_back(i1)=rm_back(i1)+abs(rm(i1,i2))**2
!!    enddo
!!  enddo

  tm_out=0.0D0
  do i2=1,n_out
    do i1=1,n_in
      tm_out(i2)=tm_out(i2)+abs(tm(i2,i1))**2
    enddo
  enddo

  rm_back=0.0D0
  do i2=1,n_back
    do i1=1,n_in
      rm_back(i2)=rm_back(i2)+abs(rm(i2,i1))**2
    enddo
  enddo



  ttotal_in=0.0D0
  rtotal_in=0.0D0
  do i2=1,n_in
    ttotal_in=ttotal_in+tm_in(i2)
    rtotal_in=rtotal_in+rm_in(i2)
  enddo

  rtotal_back=0.0D0
  do i2=1,n_back
    rtotal_back=rtotal_back+rm_back(i2)
  enddo

  ttotal_out=0.0D0
  do i2=1,n_out
    ttotal_out=ttotal_out+tm_out(i2)
  enddo


!output of transmission

  do i2=1,n_in
    write(12347,*)trim(adjustl(prefix)),": Tm_in,Rm_in,Tm_in+Rm_in=",13.6056981D0 * dreal(ene),i2,dreal(tm_in(i2)),dimag(tm_in(i2)),dreal(rm_in(i2)),dimag(rm_in(i2)),dreal(tm_in(i2)+rm_in(i2)),ispin,ikpmod
  enddo
  write(12347,*)trim(adjustl(prefix)),": T_in,R_in,T_in+R_in=",13.6056981D0 * dreal(ene),dreal(ttotal_in),dreal(rtotal_in),dreal(ttotal_in+rtotal_in),ispin,ikpmod

!!  do i2=1,n_back
!!    write(12347,*)trim(adjustl(prefix))," Rm_back=",13.6056981D0 * dreal(ene),i2,dreal(rm_back(i2)),dimag(rm_back(i2)),ispin,ikpmod
!!  enddo
!!  write(12347,*)trim(adjustl(prefix))," R_back,R_back-R_in=",13.6056981D0 * dreal(ene),dreal(rtotal_back),dreal(rtotal_back-rtotal_in),ispin,ikpmod
!!
!!  do i2=1,n_out
!!    write(12347,*)trim(adjustl(prefix))," Tm_out=",13.6056981D0 * dreal(ene),i2,dreal(tm_out(i2)),dimag(tm_out(i2)),ispin,ikpmod
!!  enddo
!!  write(12347,*)trim(adjustl(prefix))," T_out,T_out-T_in=",13.6056981D0 * dreal(ene),dreal(ttotal_out),dreal(ttotal_out-ttotal_in),ispin,ikpmod
 
  end subroutine TR_Total

end module mTransmissionDecomposition

!xxx print out difference in sigma, gamma as well; check difference in first elements of phi (not multiplied with eik)

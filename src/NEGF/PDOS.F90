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
!                   EM_DOS,
!                   EM_DOS_GENERAL,
!                   EM_DOS_NC,
!                   EM_DOS_NC2,
!                   SET_S_HALF,
!                   SET_GF_0,
!                   GETRHOPDOS_NC_noON  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

SUBROUTINE em_dos(n,gf,empdos, emdostotk,empdostotk,sgeneral)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
  use mTypes

  implicit none
  integer n,ii,jj,ind,ind2,ind3
  type(matrixTypeGeneral) :: sgeneral
  type(matrixTypeGeneral) :: gf
  double complex gfu
  double complex gf_diag(n),gf_diag2(n)
  double precision pi,totpdos
  double precision emdostotk,empdostotk(n)
  logical empdos

  pi=3.141592653589793D0

  do ii=1,n
    gf_diag(ii)=0D0
    do ind=sgeneral%matSparse%q(ii),  sgeneral%matSparse%q(ii+1)-1
      if(gf%mattype.eq.0)then
!        write(*,*)"using dense matrix"
        gfu=((0.0D0,0.5D0)/ pi) * (gf%matdense%a(sgeneral%matSparse%j(ind),ii)-DCONJG(gf%matdense%a(ii,sgeneral%matSparse%j(ind))))
      elseif(gf%mattype.eq.2)then
!        write(*,*)"using sparse matrix"
        call inddensetoindsparsegeneral(sgeneral%matSparse%j(ind),ii,ind2,gf)
        call inddensetoindsparsegeneral(ii,sgeneral%matSparse%j(ind),ind3,gf)
        if(ind2.ne.0.and.ind3.ne.0)then
          gfu=((0.0D0,0.5D0)/ pi) * (gf%matSparse%b(ind2)- DCONJG(gf%matSparse%b(ind3)))
        else
          gfu=0D0
        endif
      endif
      gf_diag(ii)=gf_diag(ii)+sgeneral%matSparse%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,n
    if(empdos)then
      empdostotk(ii)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk=totpdos/13.6057D0


end SUBROUTINE em_dos




SUBROUTINE em_dos_general(n,nuo,NSpinBlocks,NspinComplexMatrix,gf,empdos, emdostotk,empdostotk,sgeneral)

  use mTypes

  implicit none

  integer, intent(in) :: n,NSpinBlocks,nuo,NspinComplexMatrix
  type(matrixTypeGeneral), intent(inout) :: gf
  logical, intent(in) :: empdos
  type(matrixTypeGeneral), intent(inout) :: sgeneral
  double precision, intent(out) :: emdostotk(NSpinBlocks),empdostotk(nuo,NSpinBlocks)
  integer nmat

  if(NSpinBlocks<=3)then
    call em_dos(nuo,gf,empdos,emdostotk,empdostotk,sgeneral)
  else
    if(NspinComplexMatrix==4)then
      call em_dos_nc(n,NSpinBlocks,gf,empdos,emdostotk,empdostotk,sgeneral)
    else

!      call PrintSparse2DenseReorderedNC(sgeneral,n/2, 1, NSpinBlocks, NspinComplexMatrix,2,"sp_2")
!      call PrintSparse2DenseReorderedNC(gf,n/2, NspinComplexMatrix, NSpinBlocks, NspinComplexMatrix,1,"gp_2")
      call em_dos_nc2(n,NSpinBlocks,NspinComplexMatrix,gf%matSparse,empdos,emdostotk,empdostotk,sgeneral%matSparse)
    endif
  endif



end SUBROUTINE em_dos_general


SUBROUTINE em_dos_nc(n,NSpinBlocks,gf,empdos, emdostotk,empdostotk,sgeneral)

  use mTypes

  implicit none

  integer, intent(in) :: n,NSpinBlocks
  type(matrixTypeGeneral), intent(in) :: gf
  logical, intent(in) :: empdos
  type(matrixTypeGeneral), intent(in) :: sgeneral
  double precision, intent(out) :: emdostotk(NSpinBlocks),empdostotk(n/2,NSpinBlocks)

  double complex, allocatable :: GF_0(:,:),GF_x(:,:),GF_y(:,:),GF_z(:,:)
  double complex, allocatable :: A_12(:,:)
  double complex gfu
  double complex gf_diag(n/2)
  double precision pi,totpdos
  integer nhalf,ii,ind

  pi=3.141592653589793D0

  nhalf=n/2


  allocate(gf_0(nhalf,nhalf))
  gf_0=gf%matdense%a(1:nhalf,1:nhalf)+gf%matdense%a(nhalf+1:n,nhalf+1:n)

  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=sgeneral%matSparse%q(ii),  sgeneral%matSparse%q(ii+1)-1
      gfu=((0.0D0,0.5D0)/ pi) * (gf_0(sgeneral%matSparse%j(ind),ii)-DCONJG(gf_0(ii,sgeneral%matSparse%j(ind))))
      gf_diag(ii)=gf_diag(ii)+sgeneral%matSparse%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,1)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(1)=totpdos/13.6057D0

  deallocate(gf_0)


  allocate(gf_z(nhalf,nhalf))
  gf_z=gf%matdense%a(1:nhalf,1:nhalf)-gf%matdense%a(nhalf+1:n,nhalf+1:n)

  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=sgeneral%matSparse%q(ii),  sgeneral%matSparse%q(ii+1)-1
      gfu=((0.0D0,0.5D0)/ pi) * (gf_z(sgeneral%matSparse%j(ind),ii)-DCONJG(gf_z(ii,sgeneral%matSparse%j(ind))))
      gf_diag(ii)=gf_diag(ii)+sgeneral%matSparse%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,4)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(4)=totpdos/13.6057D0

  deallocate(gf_z)


  allocate(gf_x(nhalf,nhalf))
  allocate(a_12(nhalf,nhalf))
  a_12=((0.0D0,0.5D0)/pi) * (gf%matdense%a(1:nhalf,nhalf+1:n)-DCONJG(transpose(gf%matdense%a(nhalf+1:n,1:nhalf))))
  gf_x=a_12+DCONJG(transpose(a_12))

  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=sgeneral%matSparse%q(ii),  sgeneral%matSparse%q(ii+1)-1
      gfu=gf_x(sgeneral%matSparse%j(ind),ii)
      gf_diag(ii)=gf_diag(ii)+sgeneral%matSparse%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,2)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(2)=totpdos/13.6057D0

  deallocate(gf_x)
  deallocate(a_12)


  allocate(gf_y(nhalf,nhalf))
  allocate(a_12(nhalf,nhalf))
  a_12=((0.0D0,0.5D0)/pi) * (gf%matdense%a(1:nhalf,nhalf+1:n)-DCONJG(transpose(gf%matdense%a(nhalf+1:n,1:nhalf))))
  gf_y=(0.0D0,1.0D0) * (a_12-DCONJG(transpose(a_12)))

  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=sgeneral%matSparse%q(ii),  sgeneral%matSparse%q(ii+1)-1
      gfu=gf_y(sgeneral%matSparse%j(ind),ii)
      gf_diag(ii)=gf_diag(ii)+sgeneral%matSparse%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,3)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(3)=totpdos/13.6057D0

  deallocate(gf_y)
  deallocate(a_12)



end SUBROUTINE em_dos_nc


SUBROUTINE em_dos_nc2(n,NSpinBlocks,NspinComplexMatrix,gf,empdos, emdostotk,empdostotk,s)

  use mTypes
  use mMatrixUtil

  implicit none

  integer, intent(in) :: n,NSpinBlocks,NspinComplexMatrix
  logical, intent(in) :: empdos
  type(matrixSparseType), intent(in) :: gf
  type(matrixSparseType), intent(in) :: s
  double precision, intent(inout) :: emdostotk(NSpinBlocks),empdostotk(n/2,NSpinBlocks)

  type(matrixSparseType) :: Shalf
  type(matrixTypeGeneral) :: gf_2(4)
  type(matrixSparseType) :: GF_0,GF_x,GF_y,GF_z,A_12
  type(matrixSparseType) :: GF_11,GF_12,GF_21,GF_22
  type(matrixSparseType) :: GF_11d,GF_12d,GF_21d,GF_22d
  double complex gfu
  double complex gf_diag(n/2)
  double precision pi,totpdos
  integer nhalf,ii,ind
  type(ioType) :: io

  io%isDebug=.false.

  pi=3.141592653589793D0

  nhalf=n/2

  call set_s_half(nhalf,s,Shalf)
!"  call PrintSparse2DenseReorderedNC(Shalf,nhalf, 1, NSpinBlocks, 4,2,"sp_1")


!!  call set_s_half(nhalf,s,Shalf_general)
!!  call set_gf_0(nhalf,gf_2(1)%matsparse,gf_2(3)%matsparse,gf_2(4)%matsparse,gf_2(2)%matsparse,gf,Shalf_general)
!!  call PrintSparse2DenseReorderedNC(gf_2,nhalf, 4, NSpinBlocks, 4,1,"gp_1")
  call set_gf_0(nhalf,gf_11,gf_12,gf_21,gf_22,gf,Shalf)

  call AllocateMatrixCRS(nhalf,nhalf,gf_11d,gf_11%nnz,"em_dos_nc2",io)
  call mathermitianCRS(gf_11,gf_11d)
  call AllocateMatrixCRS(nhalf,nhalf,gf_22d,gf_22%nnz,"em_dos_nc2",io)
  call mathermitianCRS(gf_22,gf_22d)

  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=Shalf%q(ii),  Shalf%q(ii+1)-1
      gfu=((0.0D0,0.5D0)/ pi) * (gf_11%b(ind)-gf_11d%b(ind)+gf_22%b(ind)-gf_22d%b(ind))

      gf_diag(ii)=gf_diag(ii)+Shalf%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,1)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(1)=totpdos/13.6057D0

  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=Shalf%q(ii),  Shalf%q(ii+1)-1
      gfu=((0.0D0,0.5D0)/ pi) * (gf_11%b(ind)-gf_11d%b(ind)-gf_22%b(ind)+gf_22d%b(ind))

      gf_diag(ii)=gf_diag(ii)+Shalf%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,4)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(4)=totpdos/13.6057D0

  call AllocateMatrixCRS(nhalf,nhalf,gf_12d,gf_12%nnz,"em_dos_nc2",io)
  call mathermitianCRS(gf_12,gf_12d)
  call AllocateMatrixCRS(nhalf,nhalf,gf_21d,gf_21%nnz,"em_dos_nc2",io)
  call mathermitianCRS(gf_21,gf_21d)


  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=Shalf%q(ii),  Shalf%q(ii+1)-1
      gfu=((0.0D0,0.5D0)/ pi) * (gf_12%b(ind)-gf_12d%b(ind)+gf_21%b(ind)-gf_21d%b(ind))
!      write(12347,*)"gfu=",ii,Shalf%j(ind),gfu,gf_12%b(ind),gf_12d%b(ind),gf_21%b(ind),gf_21d%b(ind)

      gf_diag(ii)=gf_diag(ii)+Shalf%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,2)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(2)=totpdos/13.6057D0


!!!  a_12=((0.0D0,0.5D0)/pi) * (gf%matdense%a(1:nhalf,nhalf+1:n)-DCONJG(transpose(gf%matdense%a(nhalf+1:n,1:nhalf))))
!!!  gf_x=a_12+DCONJG(transpose(a_12))

  do ii=1,nhalf
    gf_diag(ii)=0D0
    do ind=Shalf%q(ii),  Shalf%q(ii+1)-1
      gfu=(0.0D0,1.0D0) * ((0.0D0,0.5D0)/ pi) * (gf_12%b(ind)+gf_12d%b(ind)-gf_21%b(ind)-gf_21d%b(ind))
!      write(12347,*)"gfu=",ii,Shalf%j(ind),gfu,gf_12%b(ind),gf_12d%b(ind),gf_21%b(ind),gf_21d%b(ind)

      gf_diag(ii)=gf_diag(ii)+Shalf%b(ind)*gfu
    enddo
  enddo

  totpdos=0D0
  do ii=1,nhalf
    if(empdos)then
      empdostotk(ii,3)=DREAL(gf_diag(ii))/13.6057D0
    endif
    totpdos=totpdos+DREAL(gf_diag(ii))
  enddo
  emdostotk(3)=totpdos/13.6057D0


  call DestroyMatrixSparse(Shalf,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_11,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_12,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_21,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_22,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_11d,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_12d,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_21d,"em_dos_nc2",io)
  call DestroyMatrixSparse(gf_22d,"em_dos_nc2",io)

!----->  a_12=((0.0D0,0.5D0)/pi) * (gf%matdense%a(1:nhalf,nhalf+1:n)-DCONJG(transpose(gf%matdense%a(nhalf+1:n,1:nhalf))))
!----->  gf_y=(0.0D0,1.0D0) * (a_12-DCONJG(transpose(a_12)))
!----->


end SUBROUTINE em_dos_nc2


SUBROUTINE  set_s_half(nhalf,s,Shalf)

  use mTypes
  use mMatrixUtil

  implicit none

  integer, intent(in) :: nhalf
  type(matrixSparseType), intent(in) :: s
  type(matrixSparseType), intent(out) :: Shalf

  integer i,j,ind,nnz,ind2
  type(ioType) :: io

!  write(12347,*)"nnz=",s%nnz,s%nnz/4,nhalf
  call AllocateMatrixCRS(nhalf,nhalf,Shalf,s%nnz/4,"em_dos_nc2",io)

  ind2=1
  Shalf%q(1)=1
  do i=1,nhalf
    do ind=s%q(2*i),s%q(2*i+1)-1
      if(mod(s%j(ind),2)==0)then
        Shalf%j(ind2)=s%j(ind)/2
!        write(12347,*)"sj=",Shalf%j(ind2),s%j(ind),ind2
        Shalf%b(ind2)=s%b(ind)
        ind2=ind2+1
      endif
    enddo
    Shalf%q(i+1)=ind2
!    write(12347,*)"si=",i,nhalf,Shalf%q(i+1),((s%q(2*i+1)-1)/4)+1
  enddo

end SUBROUTINE set_s_half


  SUBROUTINE set_gf_0(nhalf,gf_11,gf_12,gf_21,gf_22,gf,s)

  use mTypes
  use mMatrixUtil

  implicit none

  integer, intent(in) :: nhalf
  type(matrixSparseType), intent(in) :: s
  type(matrixSparseType), intent(in) :: gf
  type(matrixSparseType), intent(out) :: GF_11,gf_12,gf_21,gf_22

  integer i,j,ind,nnz,ind2,i2,i3,iup,idown,jup,jdown,iPup,iPdown
  type(ioType) :: io
  integer w(nhalf)

  call duplicatematrix(s,gf_11,.false.,io)
  call duplicatematrix(s,gf_12,.false.,io)
  call duplicatematrix(s,gf_21,.false.,io)
  call duplicatematrix(s,gf_22,.false.,io)

  gf_11%b=0.0D0
  gf_12%b=0.0D0
  gf_21%b=0.0D0
  gf_22%b=0.0D0

  w=0
  do i=1,nhalf
    do ind=s%q(i),s%q(i+1)-1
      w(s%j(ind))=ind
    enddo
    
    iPup=2*i-1
    iPdown=2*i

    do ind=gf%q(iPup),gf%q(iPup+1)-1
      if(mod(gf%j(ind),2)==1)then
        j=(gf%j(ind)+1)/2
        if(w(j).ne.0)gf_11%b(w(j))=gf%b(ind)
      else
        j=(gf%j(ind))/2
        if(w(j).ne.0)gf_12%b(w(j))=gf%b(ind)
      endif
    enddo

    do ind=gf%q(iPdown),gf%q(iPdown+1)-1
      if(mod(gf%j(ind),2)==1)then
        j=(gf%j(ind)+1)/2
        if(w(j).ne.0)gf_21%b(w(j))=gf%b(ind)
      else
        j=(gf%j(ind))/2
        if(w(j).ne.0)gf_22%b(w(j))=gf%b(ind)
      endif
    enddo

    do ind=s%q(i),s%q(i+1)-1
      w(s%j(ind))=0
    enddo
  enddo

end SUBROUTINE set_gf_0




!SUBROUTINE  set_gf_0(nhalf,gf_0,gf,s)
!
!  use mTypes
!
!  implicit none
!
!  integer, intent(in) :: nhalf
!  type(matrixSparseType), intent(in) :: gf
!  type(matrixSparseType), intent(in) :: s
!  type(matrixSparseType), intent(out) :: GF_0
!
!  integer i,j,ind
!  type(ioType) :: io
!
!  call AllocateMatrixCRS(nhalf,nhalf,gf_0,s%nnz,"em_dos_nc2",io)
!  gf_0%j=s%j
!  gf_0%q=s%q
!  gf_0%b=0.0D0
!
!
!
!  gf_0=gf%matdense%a(1:nhalf,1:nhalf)+gf%matdense%a(nhalf+1:n,nhalf+1:n)
!
!
!
!
!end SUBROUTINE set_gf_0
  SUBROUTINE GetRhoPDOS_nc_noON(rhogeneral,sgeneral,nspin,empdos,emdostotk,RhoPDOS)

    use mTypes
    use mConstants
    use mMatrixUtil, only:WriteMatrixSparse,DuplicateMatrix,DestroyMatrixSparse,Row2Cols,Col2Rows

    implicit none
    integer, intent(in) :: nspin
    type(matrixTypeGeneral), intent(in) :: rhogeneral(nspin),sgeneral
    real(kdp), intent(out) :: emdostotk(nspin)
    logical, intent(in)    :: empdos
    real(kdp), intent(out) :: RhoPDOS(sgeneral%irows,nspin)

    type(MatrixSparseType) :: rho,rhot,s
    type(ioType) :: io
    integer i,n,ind

    n=sgeneral%irows


    call DuplicateMatrix(sgeneral%matsparse,rho,.false.,io)
    call DuplicateMatrix(sgeneral%matsparse,rhot,.false.,io)
    call DuplicateMatrix(sgeneral%matsparse,s,.true.,io)

    s%b=DCONJG(s%b)


    if(empdos) RhoPDOS=0.0D0
    emdostotk=0.0D0

    rho%b=0.5D0*(rhogeneral(1)%matsparse%b+rhogeneral(2)%matsparse%b)

    if(empdos)then
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          RhoPDOS(i,1)=RhoPDOS(i,1)+rho%b(ind)*s%b(ind)
        enddo
        emdostotk(1)=emdostotk(1)+RhoPDOS(i,1)
      enddo
    else
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          emdostotk(1)=emdostotk(1)+rho%b(ind)*s%b(ind)
        enddo
      enddo
    endif

    rho%b=0.5D0*(rhogeneral(1)%matsparse%b-rhogeneral(2)%matsparse%b)

    if(empdos)then
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          RhoPDOS(i,4)=RhoPDOS(i,4)+rho%b(ind)*s%b(ind)
        enddo
        emdostotk(4)=emdostotk(4)+RhoPDOS(i,4)
      enddo
    else
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          emdostotk(4)=emdostotk(4)+rho%b(ind)*s%b(ind)
        enddo
      enddo
    endif

    rho%b=rhogeneral(3)%matsparse%b
    call mathermitianCRS(rho,rhot)

    call Row2Cols(rho,io)
    call Col2Rows(rho,io)
    call Row2Cols(rhot,io)
    call Col2Rows(rhot,io)
    call Row2Cols(s,io)
    call Col2Rows(s,io)

    rho%b=0.5D0*(rho%b+rhot%b)
    if(empdos)then
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          RhoPDOS(i,2)=RhoPDOS(i,2)+rho%b(ind)*s%b(ind)
        enddo
        emdostotk(2)=emdostotk(2)+RhoPDOS(i,2)
      enddo
    else
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          emdostotk(2)=emdostotk(2)+rho%b(ind)*s%b(ind)
        enddo
      enddo
    endif

    call DestroyMatrixSparse(rho,"GetRhoPDOS_nc_noON",io)
    call DuplicateMatrix(sgeneral%matsparse,rho,.false.,io)

    rho%b=rhogeneral(3)%matsparse%b
    call mathermitianCRS(rho,rhot)

    call Row2Cols(rho,io)
    call Col2Rows(rho,io)
    call Row2Cols(rhot,io)
    call Col2Rows(rhot,io)

    rho%b=(0.0D0,0.5D0)*(rho%b-rhot%b)
    if(empdos)then
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          RhoPDOS(i,3)=RhoPDOS(i,3)+rho%b(ind)*s%b(ind)
        enddo
        emdostotk(3)=emdostotk(3)+RhoPDOS(i,3)
      enddo
    else
      do i=1,n
        do ind=sgeneral%matsparse%q(i),sgeneral%matsparse%q(i+1)-1
          emdostotk(3)=emdostotk(3)+rho%b(ind)*s%b(ind)
        enddo
      enddo
    endif

!Check factor 2
    RhoPDOS=2.0D0 * RhoPDOS/13.6057D0
    emdostotk=2.0D0 * emdostotk/13.6057D0

    call DestroyMatrixSparse(rho,"GetRhoPDOS_nc_noON",io)
    call DestroyMatrixSparse(rhot,"GetRhoPDOS_nc_noON",io)
    call DestroyMatrixSparse(s,"GetRhoPDOS_nc_noON",io)


 end SUBROUTINE GetRhoPDOS_nc_noON



SUBROUTINE em_dos_SingleLead_general(n,nuo,NSpinBlocks,NspinComplexMatrix,empdos, emdostotk,empdostotk,sgeneral,rhogeneral)

  use mTypes

  implicit none

  integer, intent(in) :: n,NSpinBlocks,nuo,NspinComplexMatrix
  logical, intent(in) :: empdos
  type(matrixTypeGeneral), intent(in) :: sgeneral
  type(matrixTypeGeneral), intent(in) :: rhogeneral(NspinComplexMatrix)
  double precision, intent(out) :: emdostotk(NSpinBlocks),empdostotk(nuo,NSpinBlocks)

!  if(NSpinBlocks<=3)then
!    for collinear spins this still needs to be implemented
  if(NSpinBlocks>3)then
    if(NspinComplexMatrix==4)then
      call GetRhoPDOS_nc_noON(rhogeneral,sgeneral,NspinBlocks,empdos,emdostotk,empdostotk)
    else
      call GetRhoPDOS_nc_ON(rhogeneral(1)%matsparse,sgeneral%matsparse,NspinBlocks,empdos,emdostotk,empdostotk)
    endif
  endif


end SUBROUTINE em_dos_SingleLead_general



  SUBROUTINE GetRhoPDOS_nc_ON(rhosparse,ssparse,nspin,empdos,emdostotk,RhoPDOS)

    use mTypes
    use mConstants
    use mMatrixUtil, only:WriteMatrixSparse,DuplicateMatrix,DestroyMatrixSparse,Row2Cols,Col2Rows

    implicit none

    integer, intent(in):: nspin
    type(matrixSparseType), intent(in) :: rhosparse,ssparse
    real(kdp), intent(out) :: RhoPDOS(ssparse%irows/2,nspin)
    real(kdp), intent(out) :: emdostotk(nspin)
    logical, intent(in)    :: empdos

    integer i,n,ind,n2,i2

    n2=ssparse%irows
    n=ssparse%irows/2
!    write(12347,*)"n=",n,n2

    if(empdos)then
      RhoPDOS=0.0D0

      do i2=1,n2,2
        i=(i2+1)/2
        do ind=ssparse%q(i2),ssparse%q(i2+1)-1
          if(mod(rhosparse%j(ind),2)==1)then
            RhoPDOS(i,1)=RhoPDOS(i,1)+rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            RhoPDOS(i,4)=RhoPDOS(i,4)+rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            RhoPDOS(i,2)=RhoPDOS(i,2)+rhosparse%b(ind+1)*DCONJG(ssparse%b(ind))
            RhoPDOS(i,3)=RhoPDOS(i,3)-rhosparse%b(ind+1)*DCONJG(ssparse%b(ind))
          endif
        enddo
      enddo

      do i2=2,n2,2
        i=i2/2
        do ind=ssparse%q(i2),ssparse%q(i2+1)-1
          if(mod(rhosparse%j(ind),2)==0)then
            RhoPDOS(i,1)=RhoPDOS(i,1)+rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            RhoPDOS(i,4)=RhoPDOS(i,4)-rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            RhoPDOS(i,2)=RhoPDOS(i,2)+rhosparse%b(ind-1)*DCONJG(ssparse%b(ind))
            RhoPDOS(i,3)=RhoPDOS(i,3)+rhosparse%b(ind-1)*DCONJG(ssparse%b(ind))
          endif
        enddo
      enddo

      RhoPDOS=RhoPDOS/13.6057D0

      emdostotk=0.0D0
      do i=1,n
        emdostotk(:)=emdostotk(:)+RhoPDOS(i,:)
      enddo

    else

      emdostotk=0.0D0

      do i2=1,n2,2
        i=(i2+1)/2
        do ind=ssparse%q(i2),ssparse%q(i2+1)-1
          if(mod(rhosparse%j(ind),2)==1)then
            emdostotk(1)=emdostotk(1)+rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            emdostotk(4)=emdostotk(4)+rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            emdostotk(2)=emdostotk(2)+rhosparse%b(ind+1)*DCONJG(ssparse%b(ind))
            emdostotk(3)=emdostotk(3)-rhosparse%b(ind+1)*DCONJG(ssparse%b(ind))
          endif
        enddo
      enddo

      do i2=2,n2,2
        i=i2/2
        do ind=ssparse%q(i2),ssparse%q(i2+1)-1
          if(mod(rhosparse%j(ind),2)==0)then
            emdostotk(1)=emdostotk(1)+rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            emdostotk(4)=emdostotk(4)-rhosparse%b(ind)*DCONJG(ssparse%b(ind))
            emdostotk(2)=emdostotk(2)+rhosparse%b(ind-1)*DCONJG(ssparse%b(ind))
            emdostotk(3)=emdostotk(3)+rhosparse%b(ind-1)*DCONJG(ssparse%b(ind))
          endif
        enddo
      enddo

      emdostotk=emdostotk/13.6057D0

    endif

 end SUBROUTINE GetRhoPDOS_nc_ON



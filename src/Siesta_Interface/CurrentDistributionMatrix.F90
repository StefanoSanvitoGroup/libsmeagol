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
! THE SUBROUTINE
!                   CURRENTDISTRIBUTIONMATRIX  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
module mCurrDistTotal

implicit none

public CurrentDistributionMatrix
public CurrentDistributionMatrix_Vectors

contains

subroutine CurrentDistributionMatrix_Vectors(Jrho,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,indxuo,no,v)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF
  use negfmod, only: em_isa,em_iaorb,em_iphorb,em_nau,em_nas,em_nso,BohrToAng
#ifdef MPI
  use parallel
#endif

  implicit none

  integer, intent(in) :: maxnh,NspinComplexMatrix,NspinRealInputMatrix,n1,n1loc


  double precision, intent(in) :: Jrho(maxnh,NspinRealInputMatrix),xij(3,maxnh),v

  integer, intent(in) :: indxuo(no),numh(n1loc),listhptr(n1loc),listh(maxnh)

  integer io,j,no,ind,jo,iuo,juo,i,iio,ii,ispin,indstart

  integer ind2,is,maxnelerow,ind4
  double precision kxij,kpoint(3)
  type(matrixSparseType) :: jmat(NspinComplexMatrix)
  double precision, allocatable :: dx(:,:)
  double precision, allocatable :: Normdx(:)
  double precision, allocatable :: Currdx(:,:,:)
  double complex, allocatable :: hrow(:,:),srow(:)
  double complex, allocatable :: hrowbuffer(:)
  integer buffersize
  double complex, allocatable :: bbuf(:)
  double precision, allocatable :: dxbuf(:,:)
  integer, allocatable :: jbuf(:)
  
  integer  mynode,nnodes,bnode,MPIerror,neletotal,ia,ja,ia_old
  integer, allocatable :: nstartrow(:),nelerow(:),listj(:),listjval(:),nelenz(:)
  type(ioType) :: iout
#ifdef MPI
  INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

#ifdef MPI
  call MPI_Comm_Rank(negfo_comm,myNode,MPIerror)
  call MPI_Comm_Size(negfo_comm,NNodes,MPIerror)
#else
  MyNode = 0
  NNodes = 1
#endif



  iout%isDebug=.false.

  
!   Find matrix for each k-point
  allocate(nelerow(em_nau),nelenz(em_nas),listj(em_nas),listjval(em_nas))

! Find number of non-zero elements in each row of the matrix
  nelerow=0
  nelenz=0
  listj=0
  ia_old=0
  ind2=0
  do io = 1,n1
    ia=em_iaorb(io)

    if(ia.ne.ia_old)then
      do j=1,ind2
        nelenz(listj(j))=0
      enddo
      ind2=0
      ia_old=ia
    else
      do j=1,ind2
        nelenz(listj(j))=1
      enddo
    endif

    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then

      call GlobalToLocalOrb(io,MyNode,NNodes,iio)

      do j = 1,numh(iio)
        ind = listhptr(iio) + j
        jo = listh(ind)
        ja=em_iaorb(jo)

        if(nelenz(ja)==0)then
          ind2=ind2+1
          listj(ind2)=ja
          nelenz(ja)=1
        endif

      enddo

      nelerow(ia)=ind2

    endif

#ifdef MPI
    call MPI_Bcast(ind2,1,MPI_integer,BNode,negfo_comm,MPIerror)
    call MPI_Bcast(nelerow(ia),1,MPI_integer,BNode,negfo_comm,MPIerror)
    call MPI_Bcast(listj,ind2,MPI_integer,BNode,negfo_comm,MPIerror)
#endif

  enddo

!  do ia=1,em_nau
!    write(12347,*)"atomneighbors=",ia,nelerow(ia)
!  enddo

  neletotal=0
  do ia = 1,em_nau
    neletotal=neletotal+nelerow(ia)
  enddo
!  write(12347,*)"neletotal=",neletotal

  if(mynode_inverse == 0)then
    do ispin=1,NspinComplexMatrix
      call AllocateMatrixCRS(em_nau,em_nas,jmat(ispin),neletotal,"jmat",iout)
    enddo

    do ispin=1,NspinComplexMatrix
      jmat(ispin)%q(1)=1
      do ii=1,em_nau
        jmat(ispin)%q(ii+1)=jmat(ispin)%q(ii)+nelerow(ii)
      enddo
    enddo

  endif


! Find number of non-zero elements in each row of the matrix
  nelerow=0
  nelenz=0
  listj=0
  listjval=0
  ia_old=0
  ind2=0
  do ispin=1,NspinComplexMatrix
    jmat(ispin)%b=0.0D0
    jmat(ispin)%j=0
  enddo

  allocate(dx(3,jmat(1)%nnz))
  dx=0.0D0

  do io = 1,n1
    ia=em_iaorb(io)

    if(ia.ne.ia_old)then
      do j=1,ind2
        nelenz(listj(j))=0
      enddo
      ind2=0
      ia_old=ia
      indstart=jmat(1)%q(ia)-1
    else
      do j=1,ind2
        nelenz(listj(j))=listjval(j)
      enddo
    endif

    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then

      call GlobalToLocalOrb(io,MyNode,NNodes,iio)

!xxx: this does n really work for noncollinear spins?
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
        jo = listh(ind)
        ja=em_iaorb(jo)


        if(nelenz(ja)==0)then
          ind2=ind2+1
          listj(ind2)=ja
          listjval(ind2)=ind2
          nelenz(ja)=ind2
          do ispin=1,NspinComplexMatrix
            jmat(ispin)%j(indstart+ind2)=ja
          enddo
          dx(:,indstart+ind2)=xij(:,ind)
!          write(12347,*)"dx=",indstart+ind2,ind,ia,ja,dx(:,indstart+ind2)
        endif

        do ispin=1,NspinComplexMatrix
          jmat(ispin)%b(indstart+nelenz(ja))=jmat(ispin)%b(indstart+nelenz(ja))+Jrho(ind,ispin)
        enddo
      enddo

      nelerow(ia)=ind2

    endif

#ifdef MPI
    call MPI_Bcast(ind2,1,MPI_integer,BNode,negfo_comm,MPIerror)
    call MPI_Bcast(nelerow(ia),1,MPI_integer,BNode,negfo_comm,MPIerror)
    call MPI_Bcast(listj,ind2,MPI_integer,BNode,negfo_comm,MPIerror)
    call MPI_Bcast(listjval,ind2,MPI_integer,BNode,negfo_comm,MPIerror)
#endif

  enddo

!  do ia=1,em_nau
!    write(12347,*)"atomneighbors=",ia,nelerow(ia)
!  enddo

!  neletotal=0
!  do ia = 1,em_nau
!    neletotal=neletotal+nelerow(ia)
!    write(12347,*)"nelerow",nelerow(ia)
!  enddo
!  write(12347,*)"neletotal2=",neletotal


!  write(*,*)"exiting ",mynode
!  call MPI_Finalize( MPIerror )
!  stop


#ifdef MPI
  if(mynode_inverse == 0)then
 
    do ispin=1,NspinComplexMatrix
#ifdef NoMPIInPlace
      allocate(bbuf(jmat(1)%nnz))
      bbuf=jmat(ispin)%b
      call MPI_ALLREDUCE(bbuf,jmat(ispin)%b,jmat(1)%nnz,DAT_dcomplex,MPI_SUM,negfo_comm,MPIerror)
      deallocate(bbuf)
      allocate(jbuf(jmat(1)%nnz))
      jbuf=jmat(ispin)%j
      call MPI_ALLREDUCE(jbuf,jmat(ispin)%j,jmat(1)%nnz,MPI_integer,MPI_SUM,negfo_comm,MPIerror)
      deallocate(jbuf)
#else
      call MPI_ALLREDUCE(MPI_IN_PLACE,jmat(ispin)%b,jmat(1)%nnz,DAT_dcomplex,MPI_SUM,negfo_comm,MPIerror)
      call MPI_ALLREDUCE(MPI_IN_PLACE,jmat(ispin)%j,jmat(1)%nnz,MPI_integer,MPI_SUM,negfo_comm,MPIerror)
#endif
    enddo
#ifdef NoMPIInPlace
    allocate(dxbuf(3,jmat(1)%nnz))
    dxbuf=dx
    call MPI_ALLREDUCE(dxbuf,dx,jmat(1)%nnz * 3,DAT_double,MPI_SUM,negfo_comm,MPIerror)
    deallocate(dxbuf)
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,dx,jmat(1)%nnz * 3,DAT_double,MPI_SUM,negfo_comm,MPIerror)
#endif
  endif
#endif

  dx=dx * BohrToAng

  if(mynode == 0) then
    write(12347,*)"Start of current distribution matrix for bias (eV)",v * 13.6057D0
    do ispin=1,NspinComplexMatrix
      write(12347,*)"  spin index",ispin
      write(12347,*)"  atom_I,atom_J,Re(jmat),Im(jmat)"
      call PrintMatrixCRS3VectorsDouble(jmat(ispin),dx,3,jmat(ispin)%nnz,"jmat_dx",iout)
      write(12347,*)
    enddo
    write(12347,*)"End of current distribution matrix for bias (eV)",v * 13.6057D0
    write(12347,*)
  endif


  deallocate(nelerow,nelenz,listj,listjval)


  allocate(Normdx(jmat(1)%nnz))
  do ind=1,jmat(1)%nnz
    if(norm2(dx(:,ind)).ne.0.0D0)then
      Normdx(ind)=2.0D0 * norm2(dx(:,ind))
    else
      Normdx(ind)=2.0D0
    endif
  enddo



  if(mynode == 0)then
    allocate(Currdx(NspinComplexMatrix,3,jmat(1)%iRows))
    Currdx=0.0D0

    do ispin=1,NspinComplexMatrix
      do ii=1,jmat(ispin)%iRows
        do ind=jmat(ispin)%q(ii),jmat(ispin)%q(ii+1)-1
          Currdx(ispin,:,ii)=Currdx(ispin,:,ii)+DREAL(jmat(ispin)%b(ind)) * (dx(:,ind)/Normdx(ind))
        enddo
      enddo
    enddo

    write(12347,*)"Start of current distribution vectors for bias (eV)",v * 13.6057D0
    do ispin=1,NspinComplexMatrix
      write(12347,*)"  spin index",ispin
      write(12347,*)"  atom_I,current_x,current_y,current_z"
      do ii=1,jmat(ispin)%iRows
        write(12347,*)"current_vector=",ii,Currdx(ispin,1,ii),Currdx(ispin,2,ii),Currdx(ispin,3,ii)
      enddo
    enddo
    write(12347,*)"End of current distribution vectors for bias (eV)",v * 13.6057D0

    deallocate(Currdx)
  endif


  if(mynode_inverse == 0)then
    do ispin=1,NspinComplexMatrix
      call DestroyMatrixCRS(jmat(ispin),"jmat",iout)
    enddo
  endif
  deallocate(dx,Normdx)



end subroutine CurrentDistributionMatrix_Vectors

subroutine CurrentDistributionMatrix(DMImag,OmegaImag,H,S,maxnh,numh,listhptr,listh,n1loc,xij,xo,NspinRealInputMatrix,NspinComplexMatrix,N1,indxuo,no,v,PrintVectors)

  use mMatrixUtil
  use mTypes
  use mConstants
  use mMPI_NEGF
  use negfmod, only: em_isa,em_iaorb,em_iphorb,em_nau,em_nas,em_nso,BohrToAng, em_nflux,em_nbfluxStart

#ifdef MPI
  use parallel
#endif
!  use mIO_sxdhe, only : WriteRealMatPFormatted

  implicit none

  integer, intent(in) :: maxnh,NspinComplexMatrix,NspinRealInputMatrix,n1,n1loc
  logical, intent(in) :: PrintVectors


  double precision, intent(in) :: DMImag(maxnh,NspinRealInputMatrix),xij(3,maxnh),v
  double precision, intent(in) :: xo(3,n1)
  double precision, intent(in) :: OmegaImag(maxnh,NspinRealInputMatrix)
  double precision, intent(in) :: H(maxnh,NspinRealInputMatrix)
  double precision, intent(in) :: S(maxnh)

  double precision, allocatable :: Jmat(:,:)

  integer, intent(in) :: indxuo(no),numh(n1loc),listhptr(n1loc),listh(maxnh)

  integer io,j,no,ind,jo,iuo,juo,i,iio,ii,ispin,indstart,js

  integer ind2,is,maxnelerow,ind4
  double precision kxij,kpoint(3)
  double precision, allocatable :: dx(:,:)
  double precision, allocatable :: Currdx(:,:,:)
  double complex, allocatable :: hrow(:,:),srow(:)
  double complex, allocatable :: hrowbuffer(:)
  integer buffersize
  double precision Zsplit,dy,dz,ly
  real(kdp), allocatable :: JLR(:)
  real(kdp), allocatable :: JLRbuf(:)

  real(kdp), allocatable :: H0(:),Hx(:),Hy(:),Hz(:)
  real(kdp), allocatable :: ImRho0(:),ImRhox(:),ImRhoy(:),ImRhoz(:)
  real(kdp), allocatable :: ImOmega0(:),ImOmegax(:),ImOmegay(:),ImOmegaz(:)
  
  integer  mynode,nnodes,bnode,MPIerror,neletotal,ia,ja,ia_old,iu
  integer, allocatable :: nstartrow(:),nelerow(:),listj(:),listjval(:),nelenz(:)
  type(ioType) :: iout
  real(kdp) PI
  integer nl,nr
  real(kdp), parameter :: rydberg_to_ev=13.6056981D0
  logical, save :: FirstCall = .true.
#ifdef MPI
  INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

#ifdef MPI
  call MPI_Comm_Rank(negfo_comm,myNode,MPIerror)
  call MPI_Comm_Size(negfo_comm,NNodes,MPIerror)
#else
  MyNode = 0
  NNodes = 1
#endif

  PI=2.0_kdp * acos(0.0_kdp)

  iout%isDebug=.false.


  allocate(Jmat(maxnh,NspinRealInputMatrix))
  jmat=0.0D0
  if(NspinRealInputMatrix==1)then

    do ind=1,maxnh
      jmat(ind,1)=4.0_kdp* PI*(H(ind,1)*DMImag(ind,1)-S(ind)*OmegaImag(ind,1))
    enddo

  elseif(NspinRealInputMatrix==2)then


    allocate(H0(maxnh),Hz(maxnh))
    allocate(ImRho0(maxnh),ImRhoz(maxnh))
    allocate(ImOmega0(maxnh),ImOmegaz(maxnh))

    H0=0.5D0*(H(:,1)+H(:,2))
    Hz=0.5D0*(H(:,1)-H(:,2))

    ImRho0=0.5D0*(DMImag(:,1)+DMImag(:,2))
    ImRhoz=0.5D0*(DMImag(:,1)-DMImag(:,2))

    ImOmega0=0.5D0*(OmegaImag(:,1)+OmegaImag(:,2))
    ImOmegaz=0.5D0*(OmegaImag(:,1)-OmegaImag(:,2))

!    call printhsd(H0,maxnh,numh,listhptr,listh,n1loc,n1,  indxuo,no,"h0_t")
!!!    call WriteRealMatPFormatted( MyNode,Nnodes,maxnh, n1loc, n1,1, numh, listhptr, listh, H0,12347,"h0=")
!!!    call WriteRealMatPFormatted( MyNode,Nnodes,maxnh, n1loc, n1,1, numh, listhptr, listh, Hz,12347,"hz=")
!!!    call WriteRealMatPFormatted( MyNode,Nnodes,maxnh, n1loc, n1,1, numh, listhptr, listh, S,12347,"s=")
!!!    call WriteRealMatPFormatted( MyNode,Nnodes,maxnh, n1loc, n1,1, numh, listhptr, listh, ImRho0,12347,"imrho0=")
!!!    call WriteRealMatPFormatted( MyNode,Nnodes,maxnh, n1loc, n1,1, numh, listhptr, listh, ImRhoz,12347,"imrhoz=")
!!!    call WriteRealMatPFormatted( MyNode,Nnodes,maxnh, n1loc, n1,1, numh, listhptr, listh, ImOmega0,12347,"imomega0=")
!!!    call WriteRealMatPFormatted( MyNode,Nnodes,maxnh, n1loc, n1,1, numh, listhptr, listh, ImOmegaz,12347,"imomegaz=")


    do ind=1,maxnh
      jmat(ind,1)=8.0_kdp* PI*(H0(ind)*ImRho0(ind)-S(ind)*ImOmega0(ind)+Hz(ind)*ImRhoz(ind))
      jmat(ind,2)=8.0_kdp* PI*(Hz(ind)*ImRho0(ind)+H0(ind)*ImRhoz(ind)-S(ind)*ImOmegaz(ind))
    enddo

    deallocate(H0,Hz,ImRho0,ImRhoz,ImOmega0,ImOmegaz)

  else

    allocate(H0(maxnh),Hx(maxnh),Hy(maxnh),Hz(maxnh))
    allocate(ImRho0(maxnh),ImRhox(maxnh),ImRhoy(maxnh),ImRhoz(maxnh))
    allocate(ImOmega0(maxnh),ImOmegax(maxnh),ImOmegay(maxnh),ImOmegaz(maxnh))

    H0=0.5D0*(H(:,1)+H(:,2))
    Hx=H(:,3)
    Hy=-H(:,4)
    Hz=0.5D0*(H(:,1)-H(:,2))

    ImRho0=0.5D0*(DMImag(:,1)+DMImag(:,2))
    ImRhox=DMImag(:,3)
    ImRhoy=DMImag(:,4)
    ImRhoz=0.5D0*(DMImag(:,1)-DMImag(:,2))

    ImOmega0=0.5D0*(OmegaImag(:,1)+OmegaImag(:,2))
    ImOmegax=OmegaImag(:,3)
    ImOmegay=OmegaImag(:,4)
    ImOmegaz=0.5D0*(OmegaImag(:,1)-OmegaImag(:,2))

    do ind=1,maxnh

      jmat(ind,1)=8.0_kdp* PI*(H0(ind)*ImRho0(ind)-S(ind)*ImOmega0(ind)+Hx(ind)*ImRhox(ind)+Hy(ind)*ImRhoy(ind)+Hz(ind)*ImRhoz(ind))

      jmat(ind,2)=8.0_kdp* PI*(Hx(ind)*ImRho0(ind)+H0(ind)*ImRhox(ind)-S(ind)*ImOmegax(ind))
      jmat(ind,3)=8.0_kdp* PI*(Hy(ind)*ImRho0(ind)+H0(ind)*ImRhoy(ind)-S(ind)*ImOmegay(ind))
      jmat(ind,4)=8.0_kdp* PI*(Hz(ind)*ImRho0(ind)+H0(ind)*ImRhoz(ind)-S(ind)*ImOmegaz(ind))

    enddo

    deallocate(H0,Hx,Hy,Hz,ImRho0,ImRhox,ImRhoy,ImRhoz,ImOmega0,ImOmegax,ImOmegay,ImOmegaz)
  endif

  if(MyNode==0)then
    call io_assign(iu)
    IF (FirstCall) THEN
      OPEN(UNIT=iu,FILE="STT_V.dat",RECL=1000000,status='REPLACE')
      FirstCall=.false.
    ELSE
      OPEN(UNIT=iu,FILE="STT_V.dat",RECL=1000000,POSITION='APPEND')
    ENDIF
  endif

  allocate(JLR(NspinRealInputMatrix))
  do ii=1,em_nflux
    Zsplit=xo(3,em_nbfluxStart(ii))-0.5_kdp
!    write(12347,*)"Zsplit=",Zsplit,em_nflux,ii

    JLR=0.0_kdp

    do i=1,n1
!      write(12347,*)"xo_i=",xo(3,i)
      if(xo(3,i)>=Zsplit)cycle !only orbitals with xo(3,1) < Zsplit contribute to the first summation for the flux from left to right across Zsplit

      call WhichNodeOrb(i,NNodes,BNode)

      if (MyNode.eq.BNode) then
        call GlobalToLocalOrb(i,MyNode,NNodes,iio)
        do j = 1,numh(iio)
          ind = listhptr(iio) + j
          js=listh(ind)
          juo = indxuo(js)
!          ly=5.2912327691999996
!          dy=-xo(2,i)+xo(2,js)
!          if(dy>4.0D0*ly)then
!            dy=dy-7.0D0*ly
!          endif
          if(xo(3,i)+xij(3,ind)>=Zsplit)then
!            write(12347,*)"xij=",i,juo,xij(1,ind),xij(2,ind),xij(3,ind),xo(3,i)+xij(3,ind),jmat(ind,1),ii
!!            nl=45
!!            nr=45
!!            if(i<=nl.and.juo>=n1-nr+1)then
!!              write(12347,*)"sskipping xij=",i,js,juo,xij(1,ind),xij(2,ind),xij(3,ind),xo(3,i)+xij(3,ind),jmat(ind,1),ii
!!            else
            JLR(:)=JLR(:)+jmat(ind,:)
!!            endif
          endif
        enddo
      endif
    enddo

#ifdef MPI
#ifdef NoMPIInPlace
    allocate(JLRbuf(NspinRealInputMatrix))
    JLRbuf=JLR
    call MPI_ALLREDUCE(JLRbuf(1),JLR(1),NspinRealInputMatrix,DAT_double,MPI_SUM,negfo_comm,MPIerror)
    deallocate(JLRbuf)
#else
    call MPI_ALLREDUCE(MPI_IN_PLACE,JLR(1),NspinRealInputMatrix,DAT_double,MPI_SUM,negfo_comm,MPIerror)
#endif
#endif

!    if(NspinRealInputMatrix==1)then
!      write(12346,*)"JLR=  ",ii,JLR(1),v*rydberg_to_ev
!    elseif(NspinRealInputMatrix==2)then
!      write(12346,*)"JLR=  ",ii,JLR(1),JLR(2),v*rydberg_to_ev
!    elseif(NspinRealInputMatrix==4)then
!      write(12346,*)"JLR=  ",ii,JLR(1),JLR(2),JLR(3),JLR(4),v*rydberg_to_ev
!    else
!      write(12346,*)"NspinRealInputMatrix needs to be 1, 2, or 4"
!      call stopnegf
!    endif
   
    if(mynode.eq.0)then
      if(NspinRealInputMatrix==1)then
        write(iu,*)ii,JLR(1),em_nbfluxStart(ii),v*rydberg_to_ev
      elseif(NspinRealInputMatrix==2)then
        write(iu,*)ii,JLR(1),JLR(2),em_nbfluxStart(ii),v*rydberg_to_ev
      elseif(NspinRealInputMatrix==4)then
        write(iu,*)ii,JLR(1),JLR(2),JLR(3),JLR(4),em_nbfluxStart(ii),v*rydberg_to_ev
      else
        write(iu,*)"NspinRealInputMatrix needs to be 1, 2, or 4"
        call stopnegf
      endif
    endif

  enddo

  if(mynode==0)call io_close(iu)

  if(PrintVectors)then
    call CurrentDistributionMatrix_Vectors(jmat,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,indxuo,no,v)
  endif

  deallocate(JLR)
  deallocate(Jmat)


end subroutine CurrentDistributionMatrix


end module mCurrDistTotal

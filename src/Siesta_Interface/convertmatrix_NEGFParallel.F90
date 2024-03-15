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
!                   CONVERTMATRIXSIESTATOSMEAGOLPARALLEL2  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine convertmatrixsiestatosmeagolparallel2(H,S,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,indxuo,no,kpoint,hgeneral,sgeneral,rhogeneral,ematgeneral,emforces,maxnelerow,nProcs,iProc)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF
#ifdef MPI
  use parallel
#endif

  implicit none

  logical, intent(in) :: emforces
  integer, intent(in)::nProcs,iProc
  integer maxnh,NspinRealInputMatrix,NspinComplexMatrix,n1,n1loc
  double precision H(maxnh,NspinRealInputMatrix),S(maxnh),xij(3,maxnh)

  integer io,j,no,ind,jo,iuo,juo,indxuo(no),i,numh(n1loc),listhptr(n1loc),listh(maxnh),iio,ii,ispin
  integer ind2,is,maxnelerow,ind3,ind4
  double precision kxij,kpoint(3)
  type(matrixTypeGeneral) :: rhogeneral(NspinComplexMatrix),ematgeneral(NspinComplexMatrix)
  type(matrixTypeGeneral) :: hgeneral(NspinComplexMatrix),sgeneral
  double complex, allocatable :: hrow(:,:),srow(:)
  

  double complex sckxij
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  integer  mynode,nnodes,bnode,MPIerror,neletotal
  integer, allocatable :: nstartrow(:),nelerow(:),listj(:),nelenz(:)
  type(ioType) :: iout
  integer n1l,i1,i2,nele,nnz,istart,iend,nrows,neletarget,jproc,mpi_group

#ifdef MPI
  call MPI_Comm_Rank(negf_comm,myNode,MPIerror)
  call MPI_Comm_Size(negf_comm,NNodes,MPIerror)
#else
  MyNode = 0
  NNodes = 1
#endif

  iout%isDebug=.false.

! Find matrix for each k-point
  allocate(nelerow(n1),listj(n1),nelenz(n1))

! Find number of non-zero elements in each row of the matrix
  listj=0
  nelenz=0
  nelerow=0
  do io = 1,n1
    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then
      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      ind2=0
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
        juo = indxuo(listh(ind))
!        if(nelenz(juo)==0)then
        if(nelenz(juo)==0.and.((h(ind,1).ne.0.0D0).or.(s(ind).ne.0.0D0)))then
          ind2=ind2+1
          listj(ind2)=juo
          nelenz(juo)=1
        endif
      enddo
      nelerow(io)=ind2

      do j=1,ind2
        nelenz(listj(j))=0
      enddo

    endif
#ifdef MPI
      call MPI_Bcast(nelerow(io),1,MPI_integer,BNode,negf_comm,MPIerror)
#endif
  enddo
!  do i=1,n1
!    write(12347,*)"nelerowi=",i,nelerow(i),mynode
!  enddo
  maxnelerow=maxval(nelerow)

  deallocate(listj)

  nnz=sum(nelerow)
  neletarget=nnz/nProcs
!  write(12347,*)"nele=",nnz,neletarget,nProcs


  nele=0
  istart=1
  jproc=0
  iend=-1
  do io = 1,n1
    nele=nele+nelerow(io)
!    write(12347,*)"neleloop=",nnz,nele,nelerow(io),neletarget,nProcs,n1
    if(nele >= neletarget.or.io==n1)then
!      write(12347,*)"neleloopexit=",nnz,nele,nelerow(io),neletarget,nProcs,n1
      iend=io
      if(jproc==iproc)exit
      istart=io+1
      nele=0
      jproc=jproc+1
    endif
  enddo
  nrows=iend-istart+1
!  write(12347,*)"nele2=",istart,iend,nrows,nele

  allocate(nstartrow(nrows))
  neletotal=0
  nstartrow=0
  do io = 1,nrows
    nstartrow(io)=neletotal+1
    neletotal=neletotal+nelerow(io+istart-1)
  enddo

  mpi_group=-1
  do ispin=1,NspinComplexMatrix
    call AllocateMatrixGeneral(mpi_group,inverse_comm,nnodes_inverse,mynode_inverse,n1,n1,nrows,n1,1,istart,neletotal,3,hgeneral(ispin),"convertmatrix",iout)
  enddo
  call AllocateMatrixGeneral(mpi_group,inverse_comm,nnodes_inverse,mynode_inverse,n1,n1,nrows,n1,1,istart,neletotal,3,sgeneral,"convertmatrix",iout)

  do ispin=1,NspinComplexMatrix
    do ii=1,nrows
      hgeneral(ispin)%matSparseP%matSparse%q(ii)=nstartrow(ii)
    enddo
    hgeneral(ispin)%matSparseP%matSparse%q(nrows+1)=neletotal+1
  enddo

  do ii=1,nrows
    sgeneral%matSparseP%matSparse%q(ii)=nstartrow(ii)
  enddo
  sgeneral%matSparseP%matSparse%q(nrows+1)=neletotal+1

  allocate(hrow(maxnelerow,NspinComplexMatrix),srow(maxnelerow))
  allocate(listj(maxnelerow))
  listj=0
  nelenz=0
  ind3=1
  srow=0.0D0
  hrow=0.0D0
  do io = 1,n1
    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then
      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      ind2=0
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
        juo = indxuo(listh(ind))
!        if(nelenz(juo)==0)then
        if(nelenz(juo)==0.and.((h(ind,1).ne.0.0D0).or.(s(ind).ne.0.0D0)))then
          ind2=ind2+1
          listj(ind2)=juo
          nelenz(juo)=ind2
        endif
        if(nelenz(juo)==0)cycle
        ind4=nelenz(juo)

        kxij = kpoint(1) * xij(1,ind) +  kpoint(2) * xij(2,ind) 

!       note: changed sign of the complex part due to a change of the indices
        sckxij = cos(kxij) + zi*sin(kxij)

        do ispin=1,NspinComplexMatrix
          hrow(ind4,ispin)= hrow(ind4,ispin)+h(ind,ispin)*sckxij
        enddo
        srow(ind4)= srow(ind4)+s(ind)*sckxij

      enddo


      do j=1,ind2
        nelenz(listj(j))=0
      enddo

    endif


#ifdef MPI
    call MPI_Bcast(ind2,1,MPI_integer,BNode,negf_comm,MPIerror)
    call MPI_Bcast(listj(1),nelerow(io),MPI_integer,BNode,negf_comm,MPIerror)
    do ispin=1,NspinComplexMatrix
      call MPI_Bcast(hrow(1,ispin),nelerow(io),DAT_dcomplex,BNode,negf_comm,MPIerror)
    enddo
    call MPI_Bcast(srow(1),nelerow(io),DAT_dcomplex,BNode,negf_comm,MPIerror)
#endif

    if(io<=iend.and.io>=istart)then
      do ispin=1,NspinComplexMatrix
        hgeneral(ispin)%matSparseP%matSparse%b(ind3:ind3+nelerow(io)-1)=hrow(1:nelerow(io),ispin)
        hgeneral(ispin)%matSparseP%matSparse%j(ind3:ind3+nelerow(io)-1)=listj(1:nelerow(io))
      enddo
      sgeneral%matSparseP%matSparse%b(ind3:ind3+nelerow(io)-1)=srow(1:nelerow(io))
      sgeneral%matSparseP%matSparse%j(ind3:ind3+nelerow(io)-1)=listj(1:nelerow(io))
      ind3=ind3+nelerow(io)
    endif

    do j=1,ind2
      do ispin=1,NspinComplexMatrix
        hrow(j,ispin)=0.0D0
      enddo
      srow(j)=0.0D0
    enddo

  enddo
  deallocate(listj)

  deallocate(hrow,srow)

  deallocate(nelerow,nstartrow,nelenz)

  do ispin=1,NspinComplexMatrix
    call AllocateMatrixGeneral(mpi_group,inverse_comm,nnodes_inverse,mynode_inverse,n1,n1,nrows,n1,1,istart,neletotal,3,rhogeneral(ispin),"convertmatrix",iout)
  enddo
  do ispin=1,NspinComplexMatrix
    do ii=1,nrows
      rhogeneral(ispin)%matSparseP%matSparse%q(ii)=hgeneral(ispin)%matSparseP%matSparse%q(ii)
    enddo
    rhogeneral(ispin)%matSparseP%matSparse%q(nrows+1)=hgeneral(ispin)%matSparseP%matSparse%q(nrows+1)
    rhogeneral(ispin)%matSparseP%matSparse%j(:)=hgeneral(ispin)%matSparseP%matSparse%j(:)
    rhogeneral(ispin)%matSparseP%matSparse%b(:)=0.0D0
  enddo


  if(emforces)then
    do ispin=1,NspinComplexMatrix
      call AllocateMatrixGeneral(mpi_group,inverse_comm,nnodes_inverse,mynode_inverse,n1,n1,nrows,n1,1,istart,neletotal,3,ematgeneral(ispin),"convertmatrix",iout)
    enddo
    do ispin=1,NspinComplexMatrix
      do ii=1,nrows
        ematgeneral(ispin)%matSparseP%matSparse%q(ii)=hgeneral(ispin)%matSparseP%matSparse%q(ii)
      enddo
      ematgeneral(ispin)%matSparseP%matSparse%q(nrows+1)=hgeneral(ispin)%matSparseP%matSparse%q(nrows+1)
      ematgeneral(ispin)%matSparseP%matSparse%j(:)=hgeneral(ispin)%matSparseP%matSparse%j(:)
      ematgeneral(ispin)%matSparseP%matSparse%b(:)=0.0D0
    enddo
  endif


end subroutine convertmatrixsiestatosmeagolparallel2


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
!                   CONVERTMATRIXSIESTATOSMEAGOLPARALLELK2,
!                   SETISTARTIEND  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine convertmatrixsiestatosmeagolparallelK2(H,S,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,indxuo,no,hgeneral,sgeneral,rhogeneral,ematgeneral,xijk,emforces,maxnelerow,nProcs,iProc,HSblocksize)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF
  use negfmod, only : emtimings
#ifdef MPI
  use parallel
#endif

  implicit none

  logical, intent(in) :: emforces
  integer, intent(in)::nProcs,iProc
  integer, intent(out) :: HSblocksize
  integer maxnh,NspinRealInputMatrix,NspinComplexMatrix,n1,n1loc
  double precision H(maxnh,NspinRealInputMatrix),S(maxnh),xij(3,maxnh)

  integer io,j,no,ind,jo,iuo,juo,indxuo(no),i,numh(n1loc),listhptr(n1loc),listh(maxnh),iio,ii,ispin
  integer ind2,is,maxnelerow,ind3,ind4
  type(matrixTypeGeneral) :: rhogeneral(NspinRealInputMatrix),ematgeneral(NspinRealInputMatrix)
  type(matrixTypeGeneral) :: hgeneral(NspinRealInputMatrix),sgeneral
  type(matrixTypeGeneral) :: xijK(3)
  double complex, allocatable :: hrow(:,:),srow(:),xrow(:,:)
  

  double complex sckxij
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  integer  mynode,nnodes,bnode,MPIerror,neletotal
  integer, allocatable :: nstartrow(:),nelerow(:),listj(:)
  integer, allocatable :: noderowR(:),ilocalR(:)
  integer, allocatable :: noderowRbuf(:)
  integer, allocatable :: noderowS(:),ilocalS(:)
  type(ioType) :: iout
  integer n1l,i1,i2,istart,iend,nrows,mpi_group,nloc
  integer  istart2,iend2,rem,i0
  integer*4:: sc_0,sc_1,sc_r,sc_m
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


  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif




  iout%isDebug=.false.


  do io = 1,n1
    call WhichNodeOrb(io,NNodes,BNode)
    if(BNode==1)then
      HSblocksize=io-1
      exit
    endif
  enddo
!  write(12347,*)"HSblocksize=",HSblocksize


! Find matrix for each k-point
  allocate(nelerow(n1))

! Find number of non-zero elements in each row of the matrix
  nelerow=0
  do io = 1,n1
    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then
      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      nelerow(io)=numh(iio)
    endif
#ifdef MPI
    call MPI_Bcast(nelerow(io),1,MPI_integer,BNode,negfo_comm,MPIerror)
#endif
  enddo
!  do i=1,n1
!    write(12347,*)"nelerowi=",i,nelerow(i),mynode
!  enddo
  maxnelerow=maxval(nelerow)

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)') 't_cm_nelerow',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif

  call setIstartIend(istart,iend,iproc,nprocs,nelerow,n1)
  nrows=iend-istart+1

!  write(12347,*)"istart=",iproc,istart,iend,nrows

  allocate(nstartrow(nrows))
  neletotal=0
  nstartrow=0
  do io = 1,nrows
    nstartrow(io)=neletotal+1
    neletotal=neletotal+nelerow(io+istart-1)
  enddo

  mpi_group=-1
  do ispin=1,NspinRealInputMatrix
    call AllocateMatrixGeneral(mpi_group,negf_comm,nnodes_negf,mynode_negf,n1,n1,nrows,n1,1,istart,neletotal,3,hgeneral(ispin),"convertmatrix_k",iout)
  enddo
  do ispin=1,3
    call AllocateMatrixGeneral(mpi_group,negf_comm,nnodes_negf,mynode_negf,n1,n1,nrows,n1,1,istart,neletotal,3,xijK(ispin),"convertmatrix_k",iout)
  enddo
  call AllocateMatrixGeneral(mpi_group,negf_comm,nnodes_negf,mynode_negf,n1,n1,nrows,n1,1,istart,neletotal,3,sgeneral,"convertmatrix_k",iout)

  do ispin=1,NspinRealInputMatrix
    do ii=1,nrows
      hgeneral(ispin)%matSparseP%matSparse%q(ii)=nstartrow(ii)
    enddo
    hgeneral(ispin)%matSparseP%matSparse%q(nrows+1)=neletotal+1
  enddo

  do ii=1,nrows
    sgeneral%matSparseP%matSparse%q(ii)=nstartrow(ii)
  enddo
  sgeneral%matSparseP%matSparse%q(nrows+1)=neletotal+1

  do ispin=1,3
    do ii=1,nrows
      xijK(ispin)%matSparseP%matSparse%q(ii)=nstartrow(ii)
    enddo
    xijK(ispin)%matSparseP%matSparse%q(nrows+1)=neletotal+1
  enddo


  if(emtimings)then
    call system_clock(sc_1,sc_r,sc_m)
    write(12347,'(a,f12.6)') 't_cm_allocations',(sc_1-sc_0)*1.0d0/sc_r
    call system_clock(sc_0,sc_r,sc_m)
  endif



!new stuff


  
!   Find matrix for each k-point
  allocate(noderowR(n1),ilocalR(n1))

  noderowR=0
  ilocalR=0
!  write(12347,*)"ihorz=",sgeneral%matSparseP%matSparse%iVert,sgeneral%matSparseP%matSparse%iHorz,sgeneral%matSparseP%matSparse%iRows,sgeneral%matSparseP%matSparse%iCols






  nloc=sgeneral%matSparseP%matSparse%iRows
  istart2=0
  iend2=0

  i0=nloc/nnodes_groupk
  rem=mod(nloc,nnodes_groupk)
  if(mynode_groupk<rem.and.rem>0)then
    istart2=mynode_groupk * (i0+1)
    iend2=istart2+i0
  else
    istart2=rem * (i0+1) + (mynode_groupk-rem) * i0
    iend2=istart2+i0-1
  endif
!  write(12347,*)"istart,iend=",istart2,iend2,i0,rem,mynode_groupk,nnodes_groupk,nloc
  if(mynode_groupk==0)then
    istart2=0
    iend2=nloc-1
  else
    istart2=nloc+1
    iend2=-1
  endif

!  call stopnegf
  do i=sgeneral%matSparseP%matSparse%iVert+istart2,sgeneral%matSparseP%matSparse%iVert+iend2
    noderowR(i)=mynode
    ilocalR(i)=i-sgeneral%matSparseP%matSparse%iVert+1
!    write(12347,*)"noderowR(i)=",i,noderowR(i),ilocalR(i),mynode_groupk,nnodes_groupk
  enddo

#ifdef MPI
#ifdef NoMPIInPlace
  allocate(noderowRbuf(n1))
  noderowRbuf=noderowR
  call MPI_ALLREDUCE(noderowRbuf,noderowR,n1,MPI_integer,MPI_SUM,negfo_comm,MPIerror)
  deallocate(noderowRbuf)
#else
  call MPI_ALLREDUCE(MPI_IN_PLACE,noderowR,n1,MPI_integer,MPI_SUM,negfo_comm,MPIerror)
#endif
#endif

!  do i=1,n1
!    write(12347,*)"noderowR_all(i)=",i,noderowR(i)
!  enddo

  allocate(noderowS(n1),ilocalS(n1))
  noderowS=-1
  ilocalS=-1

  do io = 1,n1

    call WhichNodeOrb(io,NNodes,BNode)
    noderowS(io)=BNode
    if (MyNode.eq.BNode) then
      call GlobalToLocalOrb(io,MyNode,NNodes,ilocalS(io))
    endif

  enddo


!  do i=1,n1
!    write(12347,*)"noderowSR_all(i)=",i,noderowS(i),noderowR(i),ilocalS(i),ilocalR(i)
!  enddo





  if(emtimings)then
    call system_clock(sc_1,sc_r,sc_m)
    write(12347,'(a,f12.6)') 't_cm_e0',(sc_1-sc_0)*1.0d0/sc_r
    call system_clock(sc_0,sc_r,sc_m)
  endif



  allocate(hrow(maxnelerow,NspinRealInputMatrix),srow(maxnelerow),xrow(3,maxnelerow))
  allocate(listj(maxnelerow))
  listj=0
  ind3=1
  srow=0.0D0
  hrow=0.0D0
  xrow=0.0D0
  do io = 1,n1
 
!    if(emtimings)then
!      call system_clock(sc_0,sc_r,sc_m)
!       write(12347,*) io
!    endif


    BNode=noderowS(io)
    if (MyNode.eq.BNode) then
      iio=ilocalS(io)
      do ispin=1,NspinRealInputMatrix
        hrow(1:nelerow(io),ispin)= h(listhptr(iio)+1:listhptr(iio)+nelerow(io),ispin)
      enddo
      srow(1:nelerow(io))= s(listhptr(iio)+1:listhptr(iio)+nelerow(io))
      do ispin=1,3
        xrow(ispin,1:nelerow(io))= xij(ispin,listhptr(iio)+1:listhptr(iio)+nelerow(io))
      enddo
      listj(1:nelerow(io))= listh(listhptr(iio)+1:listhptr(iio)+nelerow(io))

 
!      if(emtimings)then
!        call system_clock(sc_1,sc_r,sc_m)
!        write(12347,'(a,f12.6)') 't_cm_e1',(sc_1-sc_0)*1.0d0/sc_r
!        call system_clock(sc_0,sc_r,sc_m)
!      endif



    endif
  

#ifdef MPI
    if(noderowR(io).ne.noderowS(io))then
 
!      if(emtimings)then
!        call system_clock(sc_0,sc_r,sc_m)
!      endif


      if (MyNode==BNode) then
        call MPI_SEND(listj(1),nelerow(io),MPI_integer, noderowR(io),1,negfo_comm,MPIerror)
        do ispin=1,NspinRealInputMatrix
          call MPI_SEND(hrow(1,ispin),nelerow(io),DAT_dcomplex, noderowR(io),1,negfo_comm,MPIerror)
        enddo
        call MPI_SEND(srow(1),nelerow(io),DAT_dcomplex, noderowR(io),1,negfo_comm,MPIerror)
        call MPI_SEND(xrow(1,1),3 * nelerow(io),DAT_dcomplex, noderowR(io),1,negfo_comm,MPIerror)

 
!        if(emtimings)then
!          call system_clock(sc_1,sc_r,sc_m)
!          write(12347,'(a,f12.6)') 't_cm_e2s',(sc_1-sc_0)*1.0d0/sc_r
!          call system_clock(sc_0,sc_r,sc_m)
!        endif



      elseif(MyNode==noderowR(io))then
        call MPI_RECV(listj(1),nelerow(io),MPI_integer, BNode, 1, negfo_comm, istatus, MPIerror)
        do ispin=1,NspinRealInputMatrix
          call MPI_RECV(hrow(1,ispin),nelerow(io),DAT_dcomplex, BNode, 1, negfo_comm, istatus, MPIerror)
        enddo
        call MPI_RECV(srow(1),nelerow(io),DAT_dcomplex, BNode, 1, negfo_comm, istatus, MPIerror)
        call MPI_RECV(xrow(1,1),3 * nelerow(io),DAT_dcomplex, BNode, 1, negfo_comm, istatus, MPIerror)

 
!        if(emtimings)then
!          call system_clock(sc_1,sc_r,sc_m)
!          write(12347,'(a,f12.6)') 't_cm_e2r',(sc_1-sc_0)*1.0d0/sc_r
!          call system_clock(sc_0,sc_r,sc_m)
!        endif
!



      endif

    endif
#endif


    if(MyNode==noderowR(io))then
      ind3=sgeneral%matSparseP%matSparse%q(ilocalR(io))
      do ispin=1,NspinRealInputMatrix
        hgeneral(ispin)%matSparseP%matSparse%b(ind3:ind3+nelerow(io)-1)=hrow(1:nelerow(io),ispin)
        hgeneral(ispin)%matSparseP%matSparse%j(ind3:ind3+nelerow(io)-1)=listj(1:nelerow(io))
      enddo
      sgeneral%matSparseP%matSparse%b(ind3:ind3+nelerow(io)-1)=srow(1:nelerow(io))
      sgeneral%matSparseP%matSparse%j(ind3:ind3+nelerow(io)-1)=listj(1:nelerow(io))
      do ispin=1,3
        xijK(ispin)%matSparseP%matSparse%b(ind3:ind3+nelerow(io)-1)=xrow(ispin,1:nelerow(io))
        xijK(ispin)%matSparseP%matSparse%j(ind3:ind3+nelerow(io)-1)=listj(1:nelerow(io))
      enddo
    endif


  enddo
  deallocate(listj)

  deallocate(hrow,srow,xrow)

  deallocate(nelerow,nstartrow)

#ifdef MPI
 
        if(emtimings)then
          call system_clock(sc_0,sc_r,sc_m)
        endif
!--->    do io=1,nnodes_groupk
!      do i=sgeneral%matSparseP%matSparse%iVert+istart2,sgeneral%matSparseP%matSparse%iVert+iend2
!        noderowR(i)=mynode
!        ilocalR(i)=i-sgeneral%matSparseP%matSparse%iVert+1
!        write(12347,*)"noderowR(i)=",i,noderowR(i),ilocalR(i),mynode_groupk,nnodes_groupk
!      enddo
      
      do ispin=1,NspinRealInputMatrix
        call MPI_Bcast(hgeneral(ispin)%matSparseP%matSparse%b(1),neletotal,DAT_dcomplex,0,groupk_comm,MPIerror)
        call MPI_Bcast(hgeneral(ispin)%matSparseP%matSparse%j(1),neletotal,MPI_integer,0,groupk_comm,MPIerror)
      enddo

      call MPI_Bcast(sgeneral%matSparseP%matSparse%j(1),neletotal,MPI_integer,0,groupk_comm,MPIerror)
      call MPI_Bcast(sgeneral%matSparseP%matSparse%b(1),neletotal,DAT_dcomplex,0,groupk_comm,MPIerror)

      do ispin=1,3
        call MPI_Bcast(xijK(ispin)%matSparseP%matSparse%b(1),neletotal,DAT_dcomplex,0,groupk_comm,MPIerror)
        call MPI_Bcast(xijK(ispin)%matSparseP%matSparse%j(1),neletotal,MPI_integer,0,groupk_comm,MPIerror)
      enddo
!--->    enddo
 
        if(emtimings)then
          call system_clock(sc_1,sc_r,sc_m)
          write(12347,'(a,f12.6)') 't_cm_bc',(sc_1-sc_0)*1.0d0/sc_r
          call system_clock(sc_0,sc_r,sc_m)
        endif

#endif

  deallocate(noderowR,ilocalR)
  deallocate(noderowS,ilocalS)


! end new stuff


  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)') 't_cm_elements',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif


  do ispin=1,NspinRealInputMatrix
    call AllocateMatrixGeneral(mpi_group,negf_comm,nnodes_negf,mynode_negf,n1,n1,nrows,n1,1,istart,neletotal,3,rhogeneral(ispin),"convertmatrix_k",iout)
  enddo
  do ispin=1,NspinRealInputMatrix
    do ii=1,nrows
      rhogeneral(ispin)%matSparseP%matSparse%q(ii)=hgeneral(ispin)%matSparseP%matSparse%q(ii)
    enddo
    rhogeneral(ispin)%matSparseP%matSparse%q(nrows+1)=hgeneral(ispin)%matSparseP%matSparse%q(nrows+1)
    rhogeneral(ispin)%matSparseP%matSparse%j(:)=hgeneral(ispin)%matSparseP%matSparse%j(:)
    rhogeneral(ispin)%matSparseP%matSparse%b(:)=0.0D0
  enddo


  if(emforces)then
    do ispin=1,NspinRealInputMatrix
      call AllocateMatrixGeneral(mpi_group,negf_comm,nnodes_negf,mynode_negf,n1,n1,nrows,n1,1,istart,neletotal,3,ematgeneral(ispin),"convertmatrix",iout)
    enddo
    do ispin=1,NspinRealInputMatrix
      do ii=1,nrows
        ematgeneral(ispin)%matSparseP%matSparse%q(ii)=hgeneral(ispin)%matSparseP%matSparse%q(ii)
      enddo
      ematgeneral(ispin)%matSparseP%matSparse%q(nrows+1)=hgeneral(ispin)%matSparseP%matSparse%q(nrows+1)
      ematgeneral(ispin)%matSparseP%matSparse%j(:)=hgeneral(ispin)%matSparseP%matSparse%j(:)
      ematgeneral(ispin)%matSparseP%matSparse%b(:)=0.0D0
    enddo
  endif


  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)') 't_cm_allocate',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif


!  call stopnegf

end subroutine convertmatrixsiestatosmeagolparallelK2



  subroutine setIstartIend(istart,iend,iproc,nprocs,nelerow,n1)
  implicit none

  integer, intent(out):: istart,iend
  integer, intent(in)::  iproc,nprocs,n1,nelerow(n1)

  integer io,nele,nnz,neletarget,jproc,iostart,ioend
  integer nnzv(n1+1)
  integer, allocatable :: nrows(:)



!  write(12347,*)"nele=",nnz,neletarget,nProcs,iproc
  nnzv(n1+1)=0
  do io=n1,1,-1
    nnzv(io)=nnzv(io+1)+nelerow(io)
  enddo


  nele=0
  jproc=0
  iostart=1
  ioend=n1-1
  do jproc=0,nProcs-1

    if(jproc==iproc)istart=iostart
    if(jproc==nprocs-1)exit
    if(iostart>n1)exit
      
    nele=0
!    nnz=sum(nelerow(iostart:n1))
    nnz=nnzv(iostart)
    neletarget=nnz/(nProcs-jproc)
!    write(12347,*)"nnz=",nnz,neletarget,jproc,nprocs
    do io = iostart,ioend
      nele=nele+nelerow(io)
      if(nele +nelerow(io+1)/2>= neletarget)then
        if(jproc==iproc)then
          iend=io
        endif
        iostart=io+1
        exit
      endif
    enddo
  enddo
  if(iproc==nprocs-1) iend=n1

!  iostart=n1+1 ! this line is just for testing

  if(iostart>n1)then
    allocate(nrows(nprocs))

    nrows=0
    jproc=0
    do io = 1,n1
      nrows(jproc+1)=nrows(jproc+1)+1
      jproc=jproc+1
      if(jproc==nprocs)jproc=0
    enddo

!    do jproc=0,nprocs-1
!      write(12347,*)"nrows=",nrows(jproc+1)
!    enddo
!    write(12347,*)"totrows=",sum(nrows)
    
    io=1
    do jproc=0,nprocs-1
      if(iproc==jproc)then
        istart=io
        iend=istart+nrows(jproc+1)-1
      endif
      io=io+nrows(jproc+1)
    enddo

    deallocate(nrows)
  endif

!  nnz=nnzv(istart)-nnzv(iend+1)
!  write(12347,*)"istartend=",iproc,istart,iend,iend-istart+1,nnz


  end subroutine setIstartIend

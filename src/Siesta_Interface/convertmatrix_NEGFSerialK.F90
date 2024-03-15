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
!                   CONVERTMATRIXSIESTATOSMEAGOLPARALLEL_2K  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

subroutine convertmatrixsiestatosmeagolparallel_2k(sgeneralp,H,S,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,nmat,indxuo,no,kpoint,hgeneral,sgeneral,rhogeneral,ematgeneral,emforces,maxnelerow)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF
  use negfmod, only : emtimings

  implicit none

  logical, intent(in) :: emforces
  integer maxnh,NspinComplexMatrix,NspinRealInputMatrix,n1,n1loc,nmat
  integer nspinMin
  double precision H(maxnh,NspinRealInputMatrix),S(maxnh),xij(3,maxnh)

  integer io,j,no,ind,jo,iuo,juo,indxuo(no),i,numh(n1loc),listhptr(n1loc),listh(maxnh),iio,ii,ispin,io2
  integer ind2,is,maxnelerow,ind4,ind5
  double precision kxij,kpoint(3)
  type(matrixTypeGeneral) :: rhogeneral(NspinComplexMatrix),ematgeneral(NspinComplexMatrix)
  type(matrixTypeGeneral) :: hgeneral(NspinComplexMatrix),sgeneral
  type(matrixSparsePType), intent(in) :: sgeneralp
  type(matrixTypeGeneral) :: hbuf
  double complex, allocatable :: hrow(:,:),srow(:)
  double complex, allocatable :: hrowbuffer(:)
  integer buffersize
  

  double complex sckxij
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  integer  mynode,nnodes,bnode,MPIerror,neletotal
  integer, allocatable :: nstartrow(:),nelerow(:),listj(:),nelenz(:),listj2(:),noderow(:),ilocal(:)
  integer, allocatable :: noderowbuf(:)
  integer*4:: sc_0,sc_1,sc_r,sc_m
  type(ioType) :: iout
#ifdef MPI
  INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

#ifdef MPI
  call MPI_Comm_Rank(negf_comm,myNode,MPIerror)
  call MPI_Comm_Size(negf_comm,NNodes,MPIerror)
#else
  MyNode = 0
  NNodes = 1
#endif

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif




  iout%isDebug=.false.

  nspinMin=NspinComplexMatrix
  if(NspinComplexMatrix>1)nspinMin=2

!   Find matrix for each k-point
  allocate(nelerow(n1),listj(n1),nelenz(n1),noderow(n1),ilocal(n1))

  noderow=0
  ilocal=0

!  write(12347,*)"ihorz=",sgeneralp%matSparse%iVert,sgeneralp%matSparse%iHorz,sgeneralp%matSparse%iRows,sgeneralp%matSparse%iCols,mynode,nnodes

!xxx here we can set bnode from siesta, so that noderow is general
  do i=sgeneralp%matSparse%iVert,sgeneralp%matSparse%iVert+sgeneralp%matSparse%iRows-1
    noderow(i)=mynode
    ilocal(i)=i-sgeneralp%matSparse%iVert+1
!    write(12347,*)"noderow(i)=",i,noderow(i)
  enddo

#ifdef MPI
#ifdef NoMPIInPlace
  allocate(noderowbuf(n1))
  noderowbuf=noderow
  call MPI_ALLREDUCE(noderowbuf,noderow,n1,MPI_integer,MPI_SUM,negf_comm,MPIerror)
  deallocate(noderowbuf)
#else
  call MPI_ALLREDUCE(MPI_IN_PLACE,noderow,n1,MPI_integer,MPI_SUM,negf_comm,MPIerror)
#endif
#endif

!  do i=1,n1
!    write(12347,*)"noderow_all(i)=",i,noderow(i)
!  enddo


! Find number of non-zero elements in each row of the matrix
  listj=0
  nelenz=0
  nelerow=0
  do io = 1,n1
!    call WhichNodeOrb(io,NNodes,BNode)
    bnode=noderow(io)
    if (MyNode.eq.BNode) then
!      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      iio=ilocal(io)
      ind2=0
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
        juo = indxuo(listh(ind))
        if(nelenz(juo)==0)then
!        if(nelenz(juo)==0.and.((h(ind,1).ne.0.0D0).or.(s(ind).ne.0.0D0)))then
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

  deallocate(listj)
  allocate(nstartrow(nmat))

  neletotal=0
  nstartrow=0
  if(nmat==n1)then
    do io = 1,n1
      nstartrow(io)=neletotal+1
      neletotal=neletotal+nelerow(io)
    enddo
  else
    nelerow=4*nelerow
    do io = 1,n1
      nstartrow(2*io-1)=neletotal+1
      neletotal=neletotal+nelerow(io)/2
      nstartrow(2*io)=neletotal+1
      neletotal=neletotal+nelerow(io)/2
    enddo
  endif
  maxnelerow=maxval(nelerow) ! 1 row stored in hrow


  if(mynode_inverse == 0)then
    do ispin=1,NspinComplexMatrix
      call AllocateMatrixGeneral(nmat,nmat,neletotal,2, hgeneral(ispin),"convertmatrix",iout)
    enddo
    call AllocateMatrixGeneral(nmat,nmat,neletotal,2, sgeneral,"convertmatrix",iout)
  endif


  if(mynode_inverse == 0)then
    do ispin=1,NspinComplexMatrix
      do ii=1,nmat
        hgeneral(ispin)%matSparse%q(ii)=nstartrow(ii)
      enddo
      hgeneral(ispin)%matSparse%q(nmat+1)=neletotal+1
    enddo

    do ii=1,nmat
      sgeneral%matSparse%q(ii)=nstartrow(ii)
    enddo
    sgeneral%matSparse%q(nmat+1)=neletotal+1

  endif

  allocate(hrow(maxnelerow,NspinComplexMatrix),srow(maxnelerow))
  allocate(listj(maxnelerow))

#ifdef MPI
!    buffersize=maxnelerow * nspin*sizeof(sckxij) + MPI_BSEND_OVERHEAD
  buffersize=maxnelerow + (MPI_BSEND_OVERHEAD/sizeof(sckxij))+1
!  write(12347,*)"buffersize1=",buffersize,maxnelerow,MPI_BSEND_OVERHEAD,sizeof(sckxij)
  allocate(hrowbuffer(buffersize))
  buffersize=buffersize * sizeof(sckxij)
!  write(12347,*)"buffersize2=",buffersize
#endif

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_cm_setup',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif


  if(n1.ne.nmat) allocate(listj2(maxnelerow))
  listj=0
  nelenz=0
  srow=0.0D0
  hrow=0.0D0
  do io = 1,n1

#ifdef MPI
        CALL MPI_BARRIER(negf_comm, MPIerror)
#endif

!    call WhichNodeOrb(io,NNodes,BNode)
    BNode=noderow(io)
    if (MyNode.eq.BNode) then
!      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      iio=ilocal(io)
      ind2=0
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
        juo = indxuo(listh(ind))
        if(nelenz(juo)==0)then
!        if(nelenz(juo)==0.and.((h(ind,1).ne.0.0D0).or.(s(ind).ne.0.0D0)))then
          ind2=ind2+1
          listj(ind2)=juo
          nelenz(juo)=ind2
        endif
        if(nelenz(juo)==0)cycle

        kxij = kpoint(1) * xij(1,ind) +  kpoint(2) * xij(2,ind) 

!       note: changed sign of the complex part due to a change of the indices
        sckxij = cos(kxij) + zi*sin(kxij)


        if(nmat==n1)then
          ind4=nelenz(juo)
          do ispin=1,nspinMin
            hrow(ind4,ispin)= hrow(ind4,ispin)+h(ind,ispin)*sckxij
          enddo
  
          if(NspinRealInputMatrix > 2)then
            hrow(ind4,3)= hrow(ind4,3)+(h(ind,3)+zi * h(ind,4))*sckxij
            if(NspinRealInputMatrix > 4)then
              hrow(ind4,1)= hrow(ind4,1)+zi * h(ind,5)*sckxij
              hrow(ind4,2)= hrow(ind4,2)+zi * h(ind,6)*sckxij
            endif
          endif
  
          srow(ind4)= srow(ind4)+s(ind)*sckxij
        else
          ind4=2*nelenz(juo)
          ind5=ind4+nelerow(io)/2

          hrow(ind4-1,1)= hrow(ind4-1,1)+h(ind,1)*sckxij
          hrow(ind4,1)= hrow(ind4,1)+(h(ind,3)+zi * h(ind,4))*sckxij

          hrow(ind5,1)= hrow(ind5,1)+h(ind,2)*sckxij

          if(NspinRealInputMatrix > 4)then
            hrow(ind4-1,1)= hrow(ind4-1,1)+zi * h(ind,5)*sckxij
            hrow(ind5,1)= hrow(ind5,1)+zi * h(ind,6)*sckxij
          endif
  
          srow(ind4-1)= srow(ind4-1)+s(ind)*sckxij
          srow(ind5)= srow(ind5)+s(ind)*sckxij

        endif

      enddo

      do j=1,ind2
        nelenz(listj(j))=0
      enddo
      
!      write(12347,*)"ind2,nele=",ind2,nelerow(io)/4,nelerow(io)

      if(nmat.ne.n1)then
        do j=1,ind2
          listj2(2*j-1)=2*listj(j)-1
          listj2(2*j)=2*listj(j)

          listj2(2*ind2+2*j-1)=2*listj(j)-1
          listj2(2*ind2+2*j)=2*listj(j)
        enddo
        listj(1:4*ind2)=listj2(1:4*ind2)
      endif

    endif

#ifdef MPI
    if (MyNode.eq.BNode) then
      if(mynode .ne. 0 )then
        do ispin=1,NspinComplexMatrix
          call MPI_SEND(listj(1),nelerow(io),MPI_integer,0,1,negf_comm,MPIerror)


!            call MPI_SEND(hrow(1,ispin),nelerow(io),DAT_dcomplex,0,1,negf_comm,MPIerror)
            call MPI_BUFFER_ATTACH(hrowbuffer,buffersize,MPIerror)
            call MPI_BSEND(hrow(1,ispin),nelerow(io),DAT_dcomplex,0,1,negf_comm,MPIerror)
            call MPI_BUFFER_DETACH(hrowbuffer,buffersize,MPIerror)

        enddo
        call MPI_SEND(listj(1),nelerow(io),MPI_integer,0,1,negf_comm,MPIerror)
!            call MPI_SEND(srow(1),nelerow(io),DAT_dcomplex,0,1,negf_comm,MPIerror)
        call MPI_BUFFER_ATTACH(hrowbuffer,buffersize,MPIerror)
        call MPI_BSEND(srow(1),nelerow(io),DAT_dcomplex,0,1,negf_comm,MPIerror)
        call MPI_BUFFER_DETACH(hrowbuffer,buffersize,MPIerror)
      endif
    endif
#endif

    if(n1==nmat)then
      io2=io
    else
      io2=2*io-1
    endif
    if(mynode == 0.and. BNode .ne. 0)then

#ifdef MPI
      do ispin=1,NspinComplexMatrix
        call MPI_RECV(hgeneral(ispin)%matSparse%j(nstartrow(io2)),nelerow(io),MPI_integer, BNode, 1, negf_comm, istatus, MPIerror)
        call MPI_RECV(hgeneral(ispin)%matSparse%b(nstartrow(io2)), nelerow(io),DAT_dcomplex, BNode, 1, negf_comm, istatus, MPIerror)
      enddo
      call MPI_RECV(sgeneral%matSparse%j(nstartrow(io2)),nelerow(io),MPI_integer, BNode, 1, negf_comm, istatus, MPIerror)
      call MPI_RECV(sgeneral%matSparse%b(nstartrow(io2)), nelerow(io),DAT_dcomplex, BNode, 1, negf_comm, istatus, MPIerror)
#endif
    elseif(mynode == 0.and. BNode .eq. 0)then
      do ispin=1,NspinComplexMatrix
        hgeneral(ispin)%matSparse%j(nstartrow(io2):nstartrow(io2)+nelerow(io)-1)=listj(1:nelerow(io))
        hgeneral(ispin)%matSparse%b(nstartrow(io2):nstartrow(io2)+nelerow(io)-1)=hrow(1:nelerow(io),ispin)
      enddo
      sgeneral%matSparse%j(nstartrow(io2):nstartrow(io2)+nelerow(io)-1)=listj(1:nelerow(io))
      sgeneral%matSparse%b(nstartrow(io2):nstartrow(io2)+nelerow(io)-1)=srow(1:nelerow(io))
    endif

    do j=1,nelerow(io)
      do ispin=1,NspinComplexMatrix
        hrow(j,ispin)=0.0D0
      enddo
      srow(j)=0.0D0
    enddo

  enddo
  deallocate(hrow,srow)
#ifdef MPI
  deallocate(hrowbuffer)
#endif
  deallocate(listj)
  if(n1.ne.nmat) deallocate(listj2)

  deallocate(nelerow,nstartrow,nelenz)

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_cm_convert',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif


  
  if(n1.ne.nmat.and.mynode==0)then
    call AllocateMatrixGeneral(nmat,nmat,hgeneral(1)%matSparse%nnz,2,hbuf,"convertmatrix",iout)

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_cm_convertnc_alloc',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif



    call mathermitianCRS(hgeneral(1)%MatSparse,hbuf%MatSparse)

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_cm_convertnc_math1',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif


    call mathermitianCRS(hbuf%MatSparse,hgeneral(1)%MatSparse)

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_cm_convertnc_math2',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif


!    call Row2Cols(hgeneral(1)%matSparse,iout)
!    call Col2Rows(hgeneral(1)%matSparse,iout)
!    call Row2Cols(sgeneral%matSparse,iout)
!    call Col2Rows(sgeneral%matSparse,iout)

    do i=1,nmat
      do ind=hgeneral(1)%matSparse%q(i),hgeneral(1)%matSparse%q(i+1)-1
        j=hgeneral(1)%matSparse%j(ind)
!        write(12347,*)"infh0",hgeneral(1)%matSparse%b(ind),ind,i,j,hbuf%matSparse%j(ind),hgeneral(1)%matSparse%i(ind)
        if(.false.)then ! xxx we should make sure H  is hermitian (and S as well)
          if(mod(i,2)==mod(j,2))then
!            write(12347,*)"infh1",hgeneral(1)%matSparse%b(ind),ind,i,j,hbuf%matSparse%j(ind)!hgeneral(1)%matSparse%i(ind)
            hgeneral(1)%matSparse%b(ind)=0.5D0*(hgeneral(1)%matSparse%b(ind)+hbuf%matSparse%b(ind))
          endif
        endif
        if((mod(i,2)==0).and.(mod(j,2)==1))then
!          write(12347,*)"infh2",hgeneral(1)%matSparse%b(ind),ind,i,j,hbuf%matSparse%j(ind)!hgeneral(1)%matSparse%i(ind)
!at this stage i is even and j is odd, we set H=H^dagger for down-up elements
          hgeneral(1)%matSparse%b(ind)=hbuf%matSparse%b(ind)
        endif
      enddo
    enddo

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_cm_convertnc_loop',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif



    call mathermitianCRS(sgeneral%matSparse,hbuf%matSparse)

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_cm_convertnc_math3',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif


    call mathermitianCRS(hbuf%matSparse,sgeneral%matSparse)
 
    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_cm_convertnc_math4',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif

  
    call DestroyMatrixGeneral(hbuf,"convertmatrix",iout)

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_cm_convertnc2',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif




  endif

#ifdef MPI
  if(mynode_inverse == 0)then
    call MPI_Bcast(sgeneral%matSparse%j(1),neletotal,MPI_integer,0,inverseheads_comm,MPIerror)
    do is = 1,NspinComplexMatrix
      hgeneral(is)%matSparse%j=sgeneral%matSparse%j
    enddo

    call MPI_Bcast(sgeneral%matSparse%b(1),2*neletotal,DAT_double,0,inverseheads_comm,MPIerror)
    do is = 1,NspinComplexMatrix
      call MPI_Bcast(hgeneral(is)%matSparse%b(1),neletotal,DAT_dcomplex,0,inverseheads_comm,MPIerror)
    enddo
  endif
#endif

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_cm_bcast',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif



  if(mynode_inverse == 0)then
    do ispin=1,NspinComplexMatrix
      call AllocateMatrixGeneral(nmat,nmat, sgeneral%matSparse%nnz,2, rhogeneral(ispin),"convertmatrix",iout)
    enddo
    do ispin=1,NspinComplexMatrix
      rhogeneral(ispin)%matSparse%q(:)=sgeneral%matSparse%q(:)
      rhogeneral(ispin)%matSparse%j(:)=sgeneral%matSparse%j(:)
      rhogeneral(ispin)%matSparse%b(:)=0.0D0
    enddo

    if(emforces)then
      do ispin=1,NspinComplexMatrix
        call AllocateMatrixGeneral(nmat,nmat, sgeneral%matSparse%nnz,2, ematgeneral(ispin),"convertmatrix",iout)
      enddo
      do ispin=1,NspinComplexMatrix
        ematgeneral(ispin)%matSparse%q(:)=sgeneral%matSparse%q(:)
        ematgeneral(ispin)%matSparse%j(:)=sgeneral%matSparse%j(:)
        ematgeneral(ispin)%matSparse%b(:)=0.0D0
      enddo
    endif

  endif

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_cm_rhoemat',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif



end subroutine convertmatrixsiestatosmeagolparallel_2k



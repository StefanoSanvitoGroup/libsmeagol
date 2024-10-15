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
!                   WRITE_SXDHE,
!                   READ_SXDHE,
!                   WRITEREALMATP,
!                   READREALMATP,
!                   PRINTHSD,
!                   PRINTHSDX  
! AND
! THE MODULE
!                   MIO_SXDHE  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
MODULE mIO_sxdhe 

  use precision
  use parallel
  use ionew
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  CONTAINS

!
subroutine write_sxdhe( maxnd, nbasis, nspin, numd, listdptr, listd, s, xij,dm,h,omega, trlabel)


  character paste*33, trlabel*(*)
  integer   maxnd, nbasis, nspin, isave
  integer   listd(maxnd), numd(nbasis), listdptr(nbasis)
  real(dp)    s(maxnd)
  real(dp)    xij(3,maxnd)
  real(dp)    dm(maxnd,nspin)
  real(dp)    H(maxnd,nspin)
  real(dp)    omega(maxnd,nspin)
  real(dp), allocatable ::    xijt(:,:)

! Internal variables and arrays
  character fname*33, sname*30
  logical   exist1, exist2, exist3, frstme, Steps
  integer   im, is, unit1, unit2, m, nb, ndmax, ns, iu
  integer   Node, Nodes, nbasistot, ml, ndmaxg
  integer, dimension(:), allocatable, save :: numdg
  real(dp), dimension(:), allocatable, save :: buffer
#ifdef MPI
  integer   MPIerror, Request, Status(MPI_Status_Size)
  integer   BNode
  integer, dimension(:), allocatable, save :: ibuffer
#endif
  external          chkdim, paste, timer, memory

  save      frstme, fname
  data      frstme /.true./

! Get the Node number
#ifdef MPI
  call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
  call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
  Node = 0
  Nodes = 1
#endif
! Find file name


  fname="file.DHE"

! Find total number of basis functions over all Nodes
#ifdef MPI
  if (Nodes.gt.1) then
    call MPI_AllReduce(nbasis,nbasistot,1,MPI_integer,MPI_sum,MPI_Comm_World,MPIerror)
  else
    nbasistot = nbasis
  endif
#else
  nbasistot = nbasis
#endif


! Allocate local buffer array for globalised numd
  allocate(numdg(nbasistot))


  if (Node.eq.0) then
    call io_assign(unit1)
    open( unit1, file=fname, form='unformatted', status='unknown' )
    rewind(unit1)
    write(unit1) nbasistot, nspin
  endif
  
  ! Create globalised numd
  do m = 1,nbasistot
#ifdef MPI
    call WhichNodeOrb(m,Nodes,BNode)
    if (BNode.eq.0.and.Node.eq.BNode) then
      call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
      ml = m
#endif
      numdg(m) = numd(ml)
#ifdef MPI
    elseif (Node.eq.BNode) then
      call GlobalToLocalOrb(m,Node,Nodes,ml)
      call MPI_ISend(numd(ml),1,MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
      call MPI_Wait(Request,Status,MPIerror)
    elseif (Node.eq.0) then
      call MPI_IRecv(numdg(m),1,MPI_integer, BNode,1,MPI_Comm_World,Request,MPIerror)
      call MPI_Wait(Request,Status,MPIerror)
    endif
    if (BNode.ne.0) then
      call MPI_Barrier(MPI_Comm_World,MPIerror)
    endif
#endif
  enddo
  
  ! Write out numd array
  if (Node.eq.0) then
    ndmaxg = 0
    do m = 1,nbasistot
      ndmaxg = max(ndmaxg,numdg(m))
    enddo
    write(unit1) (numdg(m),m=1,nbasistot)
#ifdef MPI
    allocate(ibuffer(ndmaxg))
  else
#endif
    ndmaxg=1
  endif
  
  allocate(buffer(ndmaxg))

  ! Write out listd array
  do m = 1,nbasistot
#ifdef MPI
    call WhichNodeOrb(m,Nodes,BNode)
    if (BNode.eq.0.and.Node.eq.BNode) then
      call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
      ml = m
#endif
      write(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
#ifdef MPI
    elseif (Node.eq.0) then
      call MPI_IRecv(ibuffer,numdg(m),MPI_integer,BNode,1, MPI_Comm_World,Request,MPIerror)
      call MPI_Wait(Request,Status,MPIerror)
    elseif (Node.eq.BNode) then
      call GlobalToLocalOrb(m,Node,Nodes,ml)
      call MPI_ISend(listd(listdptr(ml)+1),numd(ml),MPI_integer,  0,1,MPI_Comm_World,Request,MPIerror)
      call MPI_Wait(Request,Status,MPIerror)
    endif
    if (BNode.ne.0) then
      call MPI_Barrier(MPI_Comm_World,MPIerror)
      if (Node.eq.0) then
        write(unit1) (ibuffer(im),im=1,numdg(m))
      endif
    endif
#endif
  enddo
  
#ifdef MPI
  if (Node.eq.0) then
    call memory('D','I',size(ibuffer),'iodm')
    deallocate(ibuffer)
  endif
#endif
  
  call WriteRealMatP( node,nodes,maxnd, nbasis, nbasistot,1, numd, listdptr, listd, S,unit1,ndmaxg,buffer,numdg)
 
  allocate(xijt(maxnd,3))
!  xijt=transpose(xij)
  xijt(:,1)=xij(1,:)
  xijt(:,2)=xij(2,:)
  xijt(:,3)=xij(3,:)
  call WriteRealMatP( node,nodes,maxnd, nbasis, nbasistot,3, numd, listdptr, listd, xijt,unit1,ndmaxg,buffer,numdg)
  deallocate(xijt)
   
  call WriteRealMatP( node,nodes,maxnd, nbasis, nbasistot,nspin, numd, listdptr, listd, dm,unit1,ndmaxg,buffer,numdg)

  call WriteRealMatP( node,nodes,maxnd, nbasis, nbasistot,nspin, numd, listdptr, listd, H,unit1,ndmaxg,buffer,numdg)

  call WriteRealMatP( node,nodes,maxnd, nbasis, nbasistot,nspin, numd, listdptr, listd, omega,unit1,ndmaxg,buffer,numdg)
   

  if (Node.eq.0) then
    call io_close(unit1)
  endif
  deallocate(buffer)
  
  
  deallocate(numdg)

end subroutine write_sxdhe


!
subroutine read_sxdhe( maxnd, nbasis, nspin, numd, listdptr, listd, s, xij,dm,h,omega, trlabel)


  character paste*33, trlabel*(*)
  integer   maxnd, nbasis, nspin, isave
  integer   listd(maxnd), numd(nbasis), listdptr(nbasis)
  real(dp)    s(maxnd)
  real(dp)    xij(3,maxnd)
  real(dp)    dm(maxnd,nspin)
  real(dp)    H(maxnd,nspin)
  real(dp)    omega(maxnd,nspin)
  real(dp), allocatable ::    xijt(:,:)

! Internal variables and arrays
  character fname*33, sname*30
  logical   exist1, exist2, exist3, frstme, Steps
  integer   im, is, unit1, unit2, m, nb, ndmax, ns, iu
  integer   Node, Nodes, nbasistot, ml, ndmaxg
  integer, dimension(:), allocatable, save :: numdg
  real(dp), dimension(:), allocatable, save :: buffer
#ifdef MPI
  integer   MPIerror, Request, Status(MPI_Status_Size)
  integer   BNode
  integer, dimension(:), allocatable, save :: ibuffer
#endif
  external          chkdim, paste, timer, memory

  save      frstme, fname
  data      frstme /.true./

! Get the Node number
#ifdef MPI
  call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
  call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
  Node = 0
  Nodes = 1
#endif
! Find file name


  fname="file.DHE"

! Find total number of basis functions over all Nodes
#ifdef MPI
  if (Nodes.gt.1) then
    call MPI_AllReduce(nbasis,nbasistot,1,MPI_integer,MPI_sum,MPI_Comm_World,MPIerror)
  else
    nbasistot = nbasis
  endif
#else
  nbasistot = nbasis
#endif


! Allocate local buffer array for globalised numd
  allocate(numdg(nbasistot))

  if (Node.eq.0) then
    inquire (file=fname,         exist=exist3)
  endif

#ifdef MPI
  call MPI_Bcast(exist3,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
  
  if (Node.eq.0) then
    write(6,'(/,a)') 'iodm: Reading Density Matrix from files', fname
    call io_assign(unit1)
    open( unit1, file=fname, form='unformatted', status='unknown' )
    rewind(unit1)
    read(unit1) nb, ns
  endif
  
! Communicate the values to all Nodes and adjust to allow for
! distributed memory before checking the dimensions
#ifdef MPI
  call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
  
! Check dimensions
  call chkdim( 'iodm', 'nbasis', nbasistot, nb, 0 )
  if(nspin.lt.4)then
    call chkdim( 'iodm', 'nspin',  nspin,  ns, 0 )
  else
    call chkdim( 'iodm', 'nspin',  nspin,  ns, 1 )
  endif
  
  
  
  if (Node.eq.0) then
    read(unit1) (numdg(m),m=1,nbasistot)
  endif
#ifdef MPI
  call MPI_Bcast(numdg,nbasistot,MPI_integer,0,MPI_Comm_World, MPIerror)
#endif
  
! Convert global numd pointer to local form and generate listdptr
  ndmax = 0
  do m = 1,nbasis
    call LocalToGlobalOrb(m,Node,Nodes,ml)
    numd(m) = numdg(ml)
    ndmax = ndmax + numd(m)
    if (m .eq. 1) then
      listdptr(1) = 0
    else
      listdptr(m) = listdptr(m-1) + numd(m-1)
    endif
  enddo
  ndmaxg = 0
  do m = 1,nbasistot
    ndmaxg = max(ndmaxg,numdg(m))
  enddo
  
! Check size of first dimension of dm
  call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )
  
#ifdef MPI
! Create buffer arrays for transfering density matrix between nodes and lists
  allocate(buffer(ndmaxg))
  call memory('A','D',ndmaxg,'iodm')
  allocate(ibuffer(ndmaxg))
  call memory('A','I',ndmaxg,'iodm')
#else
  allocate(buffer(1))
#endif
  
  do m = 1,nbasistot
#ifdef MPI
    call WhichNodeOrb(m,Nodes,BNode)
    if (BNode.eq.0.and.Node.eq.BNode) then
      call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
      ml = m
#endif
      read(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
#ifdef MPI
    elseif (Node.eq.0) then
      read(unit1) (ibuffer(im),im=1,numdg(m))
      call MPI_ISend(ibuffer,numdg(m),MPI_integer, BNode,1,MPI_Comm_World,Request,MPIerror)
      call MPI_Wait(Request,Status,MPIerror)
    elseif (Node.eq.BNode) then
      call GlobalToLocalOrb(m,Node,Nodes,ml)
      call MPI_IRecv(listd(listdptr(ml)+1),numd(ml),MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
      call MPI_Wait(Request,Status,MPIerror)
    endif
    if (BNode.ne.0) then
      call MPI_Barrier(MPI_Comm_World,MPIerror)
    endif
#endif
  enddo
  
#ifdef MPI
  call memory('D','I',size(ibuffer),'iodm')
  deallocate(ibuffer)
#endif
  

  call ReadRealMatP( node,nodes,maxnd, nbasis, nbasistot,ns, numd, listdptr, S,unit1,ndmaxg,buffer,numdg)
  
! start reading xij
  allocate(xijt(maxnd,3))
  xijt=10.0D0
  call ReadRealMatP( node,nodes,maxnd, nbasis, nbasistot,3, numd, listdptr,xijt,unit1,ndmaxg,buffer,numdg)
!  xij=transpose(xijt)
  xij(1,:)=xijt(:,1)
  xij(2,:)=xijt(:,2)
  xij(3,:)=xijt(:,3)
  deallocate(xijt)
! end reading xij
  

  call ReadRealMatP( node,nodes,maxnd, nbasis, nbasistot,ns, numd, listdptr, dm,unit1,ndmaxg,buffer,numdg)

  call ReadRealMatP( node,nodes,maxnd, nbasis, nbasistot,ns, numd, listdptr, H,unit1,ndmaxg,buffer,numdg)

  call ReadRealMatP( node,nodes,maxnd, nbasis, nbasistot,ns, numd, listdptr, omega,unit1,ndmaxg,buffer,numdg)


  
  deallocate(buffer)
  if (Node.eq.0) then
    call io_close(unit1)
  endif
  
  
    ! Deallocate local buffer array for globalised numd
  call memory('D','I',size(numdg),'iodm')
  deallocate(numdg)

end subroutine read_sxdhe

subroutine WriteRealMatPFormatted( node,nodes,maxnd, nbasis, nbasistot,nspin, numd, listdptr, listd, mat,output_unit,nam)


  integer,intent(in)::   maxnd, nbasis, nspin,node,nodes,nbasistot,output_unit
  integer, intent(in)::   listd(maxnd), numd(nbasis), listdptr(nbasis)
  real(dp), intent(in)::    mat(maxnd,nspin)
  CHARACTER(LEN=*) :: nam

#ifdef MPI
  integer   MPIerror, Request, istatus(MPI_Status_Size)
  integer   BNode
#endif
  integer is,m,ml,im
  integer ndmaxg
  integer, allocatable :: numdg(:)
  real(dp), allocatable::   buffer(:)
  integer, allocatable ::   listdbuf(:)
  integer maxnumd,maxnumdin,numdval,ind,j,i

  real(dp), allocatable :: matdense(:,:)

  allocate(numdg(nbasistot))

  maxnumdin=maxval(numd)
#ifdef MPI
  CALL MPI_REDUCE(maxnumdin,maxnumd,1, MPI_INTEGER,MPI_MAX,0, MPI_Comm_World,MPIerror)
#else
  maxnumd = maxnumdin
#endif

  if(node==0) allocate(buffer(maxnumd),listdbuf(maxnumd),matdense(nbasistot,nbasistot))

  do is=1,nspin
    if(node==0)matdense=0.0_dp
    do m=1,nbasistot
#ifdef MPI
      call WhichNodeOrb(m,Nodes,BNode)
      if (BNode.eq.0.and.Node.eq.BNode) then
        call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
        ml = m
#endif
        numdval=numd(ml)
        listdbuf(1:numd(ml))=listd(listdptr(ml)+1:listdptr(ml)+numd(ml))
        buffer(1:numd(ml))=mat(listdptr(ml)+1:listdptr(ml)+numd(ml),is)
#ifdef MPI
      elseif (Node.eq.0) then
        call MPI_RECV(numdval, 1,mpi_integer,BNode,6,MPI_Comm_World,istatus,MPIerror)
        call MPI_IRecv(listdbuf,numdval,MPI_INTEGER, BNode,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,istatus,MPIerror)

        call MPI_IRecv(buffer,numdval,DAT_double, BNode,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,istatus,MPIerror)
      elseif (Node.eq.BNode) then
        call GlobalToLocalOrb(m,Node,Nodes,ml)
        call MPI_SEND(numd(ml), 1,mpi_integer,0,6,MPI_Comm_World,MPIerror)

        call MPI_ISend(listd(listdptr(ml)+1),numd(ml),MPI_INTEGER,0,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,istatus,MPIerror)

        call MPI_ISend(mat(listdptr(ml)+1,is),numd(ml),DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,istatus,MPIerror)
      endif
      call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif

      if (Node.eq.0) then
        do ind=1,numdval
          matdense(m,listdbuf(ind))=buffer(ind)
          write(output_unit,*)nam,m,listdbuf(ind),buffer(ind),is
        enddo
      endif


    enddo
    if (Node.eq.0) then
      do i=1,nbasistot
        do j=1,nbasistot
          write(output_unit,*)"dense:",nam,i,j,matdense(i,j),is
        enddo
      enddo
    endif


  enddo

  if(Node==0)deallocate(buffer,listdbuf,matdense)

end subroutine WriteRealMatPFormatted


subroutine WriteRealMatP( node,nodes,maxnd, nbasis, nbasistot,nspin, numd, listdptr, listd, mat,output_unit,ndmaxg,buffer,numdg)



  integer,intent(in)::   maxnd, nbasis, nspin,node,nodes,nbasistot,output_unit,ndmaxg
  integer, intent(in)::   listd(maxnd), numd(nbasis), listdptr(nbasis),numdg(nbasistot)
  real(dp), intent(in)::    mat(maxnd,nspin)
  real(dp), intent(inout)::    buffer(ndmaxg)

#ifdef MPI
  integer   MPIerror, Request, Status(MPI_Status_Size)
  integer   BNode
#endif
  integer is,m,ml,im


  do is=1,nspin
    do m=1,nbasistot
#ifdef MPI
      call WhichNodeOrb(m,Nodes,BNode)
      if (BNode.eq.0.and.Node.eq.BNode) then
        call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
        ml = m
#endif
        write(output_unit) (mat(listdptr(ml)+im,is),im=1,numd(ml))
#ifdef MPI
      elseif (Node.eq.0) then
        call MPI_IRecv(buffer,numdg(m),DAT_double, BNode,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,Status,MPIerror)
      elseif (Node.eq.BNode) then
        call GlobalToLocalOrb(m,Node,Nodes,ml)
        call MPI_ISend(mat(listdptr(ml)+1,is),numd(ml),DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,Status,MPIerror)
      endif
      if (BNode.ne.0) then
        call MPI_Barrier(MPI_Comm_World,MPIerror)
        if (Node.eq.0) then
          write(output_unit) (buffer(im),im=1,numdg(m))
        endif
      endif
#endif
    enddo
  enddo


end subroutine WriteRealMatP




subroutine ReadRealMatP( node,nodes,maxnd, nbasis, nbasistot,nspin, numd, listdptr, mat,output_unit,ndmaxg,buffer,numdg)


  integer,intent(in)::   maxnd,nbasis,nspin,node,nodes,nbasistot,output_unit,ndmaxg
  integer, intent(in)::   numd(nbasis), listdptr(nbasis),numdg(nbasistot)
  real(dp), intent(out)::    mat(maxnd,nspin)
  real(dp), intent(inout)::    buffer(ndmaxg)

#ifdef MPI
  integer   MPIerror, Request, Status(MPI_Status_Size)
  integer   BNode
#endif
  integer is,m,ml,im


  do is = 1,nspin
    do m = 1,nbasistot
#ifdef MPI
      call WhichNodeOrb(m,Nodes,BNode)
      if (BNode.eq.0.and.Node.eq.BNode) then
        call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
        ml = m
#endif
        read(output_unit) (mat(listdptr(ml)+im,is),im=1,numd(ml))
#ifdef MPI
      elseif (Node.eq.0) then
        read(output_unit) (buffer(im),im=1,numdg(m))
        call MPI_ISend(buffer,numdg(m),DAT_double, BNode,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,Status,MPIerror)
      elseif (Node.eq.BNode) then
        call GlobalToLocalOrb(m,Node,Nodes,ml)
        call MPI_IRecv(mat(listdptr(ml)+1,is),numd(ml), DAT_double,0,1,MPI_Comm_World,Request,MPIerror)
        call MPI_Wait(Request,Status,MPIerror)
      endif
      if (BNode.ne.0) then
        call MPI_Barrier(MPI_Comm_World,MPIerror)
      endif
#endif
    enddo
  enddo


end subroutine ReadRealMatP



subroutine printhsd(S,maxnh,numh,listhptr,listh,n1loc,N1,indxuo,no,nam)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF
#ifdef MPI
  use parallel
#endif

  implicit none


  integer maxnh,NspinComplexMatrix,NspinRealInputMatrix,n1,n1loc
  integer nspinMin
  double precision S(maxnh)

  integer io,j,no,ind,jo,iuo,juo,indxuo(no),i,numh(n1loc),listhptr(n1loc),listh(maxnh),iio,ii,ispin
  integer ind2,is,maxnelerow,ind4,mynode,nnodes,mpierror,bnode
  CHARACTER(LEN=*) :: nam

#ifdef MPI
  call MPI_Comm_Rank(negfo_comm,myNode,MPIerror)
  call MPI_Comm_Size(negfo_comm,NNodes,MPIerror)
#else
  MyNode = 0
  NNodes = 1
#endif




  do io = 1,n1

#ifdef MPI
    CALL MPI_BARRIER(negfo_comm, MPIerror)
#endif

    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then
      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
!        juo = indxuo(listh(ind))
        write(12347,*)nam,io,listh(ind),s(ind)
      enddo
    endif

  enddo



end subroutine printhsd





subroutine printhsdx(S,xijo,maxnh,numh,listhptr,listh,n1loc,N1,indxuo,no,nam)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF
#ifdef MPI
  use parallel
#endif

  implicit none


  integer maxnh,NspinComplexMatrix,NspinRealInputMatrix,n1,n1loc
  integer nspinMin
  double precision S(maxnh)
  double precision xijo(3,maxnh)

  integer io,j,no,ind,jo,iuo,juo,indxuo(no),i,numh(n1loc),listhptr(n1loc),listh(maxnh),iio,ii,ispin
  integer ind2,is,maxnelerow,ind4,mynode,nnodes,mpierror,bnode
  CHARACTER(LEN=*) :: nam

#ifdef MPI
  call MPI_Comm_Rank(negfo_comm,myNode,MPIerror)
  call MPI_Comm_Size(negfo_comm,NNodes,MPIerror)
#else
  MyNode = 0
  NNodes = 1
#endif




  do io = 1,n1

#ifdef MPI
    CALL MPI_BARRIER(negfo_comm, MPIerror)
#endif

    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then
      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
!        juo = indxuo(listh(ind))
        write(12347,'(a,2i6,4E20.5)')nam,io,listh(ind),s(ind),xijo(1,ind),xijo(2,ind),xijo(3,ind)
      enddo
    endif

  enddo






end subroutine printhsdx




END MODULE mIO_sxdhe

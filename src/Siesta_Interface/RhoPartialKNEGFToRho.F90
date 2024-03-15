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
!                   RHONEGFK_TO_RHOSIESTA  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

subroutine rhonegfk_to_rhosiesta(maxnelerow,NspinRealInputMatrix,rhogeneral,Dnew,maxnh,no,nuo,nuotot, listh, numh, listhptr,indxuo)

      use mTypes
      use mMPI_NEGF
#ifdef MPI
      use parallel
#endif

      implicit none

      integer maxnh,NspinRealInputMatrix
      type(matrixTypeGeneral) :: rhogeneral(NspinRealInputMatrix)
      double precision Dnew(maxnh,NspinRealInputMatrix)
      integer maxnelerow
      integer no,nuo,nuotot,listh(maxnh), numh(nuo), listhptr(nuo),indxuo(no)
      double precision kpoint(3),xij(3,maxnh),wk

      double complex, allocatable :: rhorow(:,:)
      integer io,BNode,MPIerror,nelerow,ind,iuo,juo,ispin,j,jo,iio,ind2,i
      integer, allocatable :: ind2ofj(:)
      integer, allocatable :: noderow(:),ilocal(:)
      integer, allocatable :: noderowbuf(:)
#ifdef MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

!   Find matrix for each k-point
  allocate(noderow(nuotot),ilocal(nuotot))

  noderow=0
  ilocal=0
!  write(12347,*)"ihorz=",rhogeneral(1)%matsparseP%matSparse%iVert,rhogeneral(1)%matsparseP%matSparse%iHorz,rhogeneral(1)%matsparseP%matSparse%iRows,rhogeneral(1)%matsparseP%matSparse%iCols

  if(mynode_groupk==0)then
    do i=rhogeneral(1)%matsparseP%matSparse%iVert,rhogeneral(1)%matsparseP%matSparse%iVert+rhogeneral(1)%matsparseP%matSparse%iRows-1
      noderow(i)=mynode_negfo
      ilocal(i)=i-rhogeneral(1)%matsparseP%matSparse%iVert+1
!      write(12347,*)"noderow(i)=",i,noderow(i)
    enddo
  endif

#ifdef MPI
#ifdef NoMPIInPlace
  allocate(noderowbuf(nuotot))
  noderowbuf=noderow
  call MPI_ALLREDUCE(noderowbuf,noderow,nuotot,MPI_integer,MPI_SUM,negf_comm,MPIerror)
  deallocate(noderowbuf)
#else
  call MPI_ALLREDUCE(MPI_IN_PLACE,noderow,nuotot,MPI_integer,MPI_SUM,negfo_comm,MPIerror)
#endif
#endif

!  do i=1,nuotot
!    write(12347,*)"noderowr_all(i)=",i,noderow(i)
!  enddo


       allocate(rhorow(maxnelerow,NspinRealInputMatrix))

       do io = 1,nuotot
         call WhichNodeOrb(io,nnodes_negfo,BNode)
#ifdef MPI
         CALL MPI_BARRIER(negfo_comm, MPIerror)
#endif
         if(mynode_negfo.eq.BNode)then

           if(mynode_negfo.ne.noderow(io))then
#ifdef MPI
             call MPI_RECV(nelerow,1,MPI_integer,noderow(io),1,negfo_comm,istatus,MPIerror)
             do ispin=1,NspinRealInputMatrix
               call MPI_RECV(rhorow(1,ispin),nelerow,DAT_dcomplex,noderow(io),1, negfo_comm,istatus,MPIerror)
             enddo
#endif
           else
             iio=ilocal(io)
             nelerow=rhogeneral(1)%matSparseP%matSparse%q(iio+1)-rhogeneral(1)%matSparseP%matSparse%q(iio)
             do ispin=1,NspinRealInputMatrix
               rhorow(1:nelerow,ispin)=rhogeneral(ispin)%matSparseP%matSparse%b(rhogeneral(1)%matSparseP%matSparse%q(iio): rhogeneral(1)%matSparseP%matSparse%q(iio+1)-1)
             enddo
           endif

           call GlobalToLocalOrb(io,mynode_negfo,nnodes_negfo,iio)

!           write(12347,*)"numhi=",iio,nelerow,numh(iio)

           do ispin=1,NspinRealInputMatrix
             Dnew(listhptr(iio)+1:listhptr(iio)+numh(iio),ispin)= rhorow(1:numh(iio),ispin)
           enddo

#ifdef MPI
         elseif(mynode_negfo.eq.noderow(io))then

           iio=ilocal(io)
           nelerow=rhogeneral(1)%matSparseP%matSparse%q(iio+1)-rhogeneral(1)%matSparseP%matSparse%q(iio)
           call MPI_SEND(nelerow,1,MPI_integer, BNode,1,negfo_comm,MPIerror)
           do ispin=1,NspinRealInputMatrix
             call MPI_SEND(rhogeneral(ispin)%matSparseP%matSparse%b(rhogeneral(1)%matSparseP%matSparse%q(iio)),nelerow,DAT_dcomplex,BNode,1,negfo_comm,MPIerror)
           enddo

#endif
         endif
#ifdef MPI
        CALL MPI_BARRIER(negfo_comm, MPIerror)
#endif
       enddo
       deallocate(rhorow)





       deallocate(noderow,ilocal)

!       call writematreal(Dnew,maxnh,NspinComplexMatrix,"dsck")
!       call MPI_Finalize( MPIerror )
!       stop

 end subroutine rhonegfk_to_rhosiesta



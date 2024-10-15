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
!                   CREATE_COMMUNICATORS_NEGF,
!                   DESTROY_COMMUNICATORS_NEGF,
!                   WHICHNODEORB,
!                   GLOBALTOLOCALORB  
! AND
! THE MODULE
!                   MMPI_NEGF  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
module mMPI_NEGF

#ifdef MPI
 use mpi_siesta
#endif

 implicit none
 private

#ifdef MPI
 public :: MPI_integer
 public :: MPI_logical
 public :: MPI_character
 public :: MPI_SUM
#ifndef NoMPIInPlace
 public :: mpi_in_place
#endif
 public :: MPI_STATUS_SIZE
 public :: MPI_BSEND_OVERHEAD
 public :: DAT_double
 public :: DAT_dcomplex
#endif

 public :: create_communicators_negf
 public :: destroy_communicators_negf
 integer, public :: negfo_comm,nnodes_negfo,mynode_negfo
! integer, public :: inverseheads_kimag_comm,myhead_i,nheads_i
! integer, public :: inverseheads_kreal_comm,myhead_r,nheads_r
! integer, public :: inverseheads_transm_comm,myhead_t,nheads_t
 integer, public :: inverseheads_comm,myhead,nheads
 integer, public :: inverse_comm,mynode_inverse,nnodes_inverse
 integer, public :: negf_comm,mynode_negf,nnodes_negf
 integer, public :: groupk_comm,mynode_groupk,nnodes_groupk
! integer, public :: sigma_comm,mynode_sigma,nnodes_sigma
! integer, public :: sigmaheads_comm,mynodeheads_sigma,nnodesheads_sigma
! integer, public :: kpts_comm,mynode_kpts,nnodes_kpts
#ifndef MPI
  public whichnodeorb
  public globaltolocalorb
#endif

contains


 subroutine create_communicators_negf(parent_comm,nprocs_inverse,NParallelK)


  implicit none

  integer, intent(in):: parent_comm,nprocs_inverse,NParallelK

  integer group_world,zeronode, inverse_group,k_group
  integer, allocatable :: members_inverse(:)
  integer, allocatable :: members_K(:),groups_K(:)
  integer MPIerror,ii,nprocs_k

#ifdef MPI

  call MPI_COMM_DUP(parent_comm, negfo_comm, MPIerror)
  call MPI_COMM_SIZE(negfo_comm,nnodes_negfo,MPIerror)
  call MPI_COMM_RANK(negfo_comm,mynode_negfo,MPIerror)


  call MPI_Comm_group(negfo_comm,group_world,MPIerror)


  if(mod(nnodes_negfo,NParallelK).ne.0)then
    if(mynode_negfo==0)then
      write(*,*)"The total number of MPI processes must be an integer multiple of EM.ParallelOverKNum."
      write(*,*)"Please change either the number of MPI processes or the value of EM.ParallelOverKNum."
    endif
    call stopnegf
  endif

  nprocs_k=nnodes_negfo/NParallelK
!  write(12347,*)"nprocs_k=",nprocs_k,nnodes_negfo,NParallelK
  allocate(members_K(nprocs_k))
  zeronode=mynode_negfo-mod(mynode_negfo,nprocs_k)
!  write(12347,*)"zeronode_k=",zeronode
  do ii=1,nprocs_k
    members_K(ii)=zeronode+ii-1
  enddo
  call MPI_Group_incl(group_world, nprocs_k, members_K,k_group, MPIerror)
  deallocate(members_K)

  call MPI_COMM_CREATE(negfo_comm, k_group, negf_comm, MPIerror)
  call MPI_Group_free(k_group,MPIerror)

  CALL MPI_COMM_RANK(negf_comm,mynode_negf,MPIerror)
  CALL MPI_COMM_SIZE(negf_comm,nnodes_negf,MPIerror)

!  write(12347,*)"negf_comm_info=",mynode_negf,nnodes_negf,mynode_negfo,nnodes_negfo


  allocate(groups_K(NParallelK))
  zeronode=mod(mynode_negfo,nprocs_k)
!  write(*,*)"zeronode_group_k=",zeronode
  do ii=1,NParallelK
!    write(12347,*)"ii=",ii
    groups_K(ii)=zeronode+(ii-1)*nprocs_k
  enddo
!"  do ii=1,NParallelK
!"    write(12347,*)"groups_K(ii)=",ii,groups_K(ii)
!"  enddo
  call MPI_Group_incl(group_world, NParallelK, groups_K,k_group, MPIerror)
  deallocate(groups_K)


  call MPI_COMM_CREATE(negfo_comm, k_group, groupk_comm, MPIerror)
  call MPI_Group_free(k_group,MPIerror)

  CALL MPI_COMM_RANK(groupk_comm,mynode_groupk,MPIerror)
  CALL MPI_COMM_SIZE(groupk_comm,nnodes_groupk,MPIerror)

!  write(12347,*)"groupk_comm_info=",mynode_groupk,nnodes_groupk,mynode_negfo,nnodes_negfo

  call MPI_Group_free(group_world,MPIerror)


!  write(*,*)"nprocs_inverse=",nprocs_inverse

  call MPI_Comm_group(negf_comm,group_world,MPIerror)
  allocate(members_inverse(nprocs_inverse))
  zeronode=mynode_negf-mod(mynode_negf,nprocs_inverse)
  do ii=1,nprocs_inverse
    members_inverse(ii)=zeronode+ii-1
  enddo
  call MPI_Group_incl(group_world, nprocs_inverse, members_inverse,inverse_group, MPIerror)
  deallocate(members_inverse)

  call MPI_COMM_CREATE(negf_comm, inverse_group, inverse_comm, MPIerror)
  call MPI_Group_free(inverse_group,MPIerror)

  CALL MPI_COMM_RANK(inverse_comm,mynode_inverse,MPIerror)
  CALL MPI_COMM_SIZE(inverse_comm,nnodes_inverse,MPIerror)




  if(mynode_inverse.eq.0)then
    call MPI_COMM_SPLIT(negf_comm, 1 ,mynode_negf/nprocs_inverse, inverseheads_comm, MPIerror)
    CALL MPI_COMM_RANK(inverseheads_comm,myhead,MPIerror)
    CALL MPI_COMM_SIZE(inverseheads_comm,nheads,MPIerror)
  else
    call MPI_COMM_SPLIT(negf_comm, MPI_UNDEFINED,mynode_negf/nprocs_inverse, inverseheads_comm, MPIerror)
  endif

  call MPI_BCAST(myhead,1,MPI_INTEGER,0,inverse_comm,MPIerror)
  call MPI_BCAST(nheads,1,MPI_INTEGER,0,inverse_comm,MPIerror)


  call MPI_Group_free(group_world,MPIerror)


#else
  myhead=0
  nheads=1
  mynode_negf=0
  nnodes_negf=1
  mynode_negfo=0
  nnodes_negfo=1
  mynode_inverse=0
  nnodes_inverse=1
!  write(*,*)"nheads=",myhead,nheads
#endif

 end subroutine create_communicators_negf




 subroutine destroy_communicators_negf()

#ifdef MPI
  integer MPIerror

  if(mynode_inverse.eq.0)then
    call MPI_Comm_free(inverseheads_comm,MPIerror)
  endif
  call MPI_Comm_free(inverse_comm,MPIerror)
  call MPI_Comm_free(negf_comm,MPIerror)
  call MPI_Comm_free(negfo_comm,MPIerror)
  call MPI_Comm_free(groupk_comm,MPIerror)
#endif


 end subroutine destroy_communicators_negf

#ifndef MPI
  subroutine whichnodeorb(nlocal,nnodes,nodeglobal)
    integer, intent(in) :: nlocal,nnodes
    integer, intent(out) :: nodeglobal

    nodeglobal=0
  end subroutine whichnodeorb

  subroutine globaltolocalorb(nglobal,mynode,nnodes,nlocal)
    integer, intent(in) :: nglobal,nnodes,mynode
    integer, intent(out) :: nlocal

    nlocal=nglobal

  end subroutine globaltolocalorb
#endif

end module mMPI_NEGF


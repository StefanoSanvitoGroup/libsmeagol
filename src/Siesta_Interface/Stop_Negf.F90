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
!                   STOPNEGF,
!                   STOP_NOT_IMPLEMENTED  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
  subroutine stopnegf

  use mMPI_NEGF
  use negfmod, only : outinfo
  integer MPIerror
 
  if (outinfo) write(12347,*)"exiting"
#ifdef MPI
!        call MPI_Barrier(MPI_Comm_World,MPIerror)
  call MPI_Finalize( MPIerror )
#endif
  stop

  end subroutine stopnegf


  SUBROUTINE stop_not_implemented(message)


#ifdef MPI
  use mpi_siesta
#endif
  IMPLICIT NONE
  CHARACTER(LEN=*) :: message
  integer mynode
#ifdef MPI
  INTEGER :: MPIerror
#endif
 
#ifdef MPI
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,mynode,MPIerror)
#else
  MyNode=0
#endif

  if(mynode.eq.0)then
    write(*,*)message
    write(*,*)"stopping program"
  endif

#ifdef MPI
  call MPI_Barrier(MPI_Comm_World,MPIerror)
  call MPI_Finalize( MPIerror )
#endif
  stop

  END SUBROUTINE stop_not_implemented



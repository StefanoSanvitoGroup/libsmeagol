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
!                   RHOG_REDUCE  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

      subroutine rhoG_reduce(N1,NL,NR,NSPIN,NspinComplexMatrix,
     .    rhobs_general, rhogeneral)
 
! *****************************************************************
! Written by Alexandre Reily Rocha and Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: rochaa@tcd.ie, reilya@if.usp.br, runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:      October 2008
! *****************************************************************

      
      use negfmod
      use mTypes
      USE mMPI_NEGF

      implicit none
      integer, intent(in) :: NspinComplexMatrix
      integer :: N1,NL,NR,NSPIN,II,JJ,ISPIN,
     &    ind,j
      type(matrixTypeGeneral) :: rhobs_general(NspinComplexMatrix)
      type(matrixTypeGeneral) :: rhogeneral(NspinComplexMatrix)
       
      double complex, parameter :: zi=(0.0d0,1.0d0)
      integer nspinMin,i,ind1,ind2,nsend,nmax
      integer*4:: sc_0,sc_1,sc_r,sc_m

#ifdef MPI
      double complex, allocatable::rhobuffsparse(:)
!      double complex::rhobuffsparse(rhogeneral(1)%matSparse%nnz)
      INTEGER :: MPIerror
#endif
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


#ifdef MPI
      DO ISPIN=1,NspinComplexMatrix
!        rhobuffsparse=0.0d0
        nmax=0
        do i=1,rhogeneral(ispin)%matSparse%iRows
          if(rhogeneral(ispin)%matSparse%q(i+1)-
     &        rhogeneral(ispin)%matSparse%q(i).gt.nmax)then
            nmax=rhogeneral(ispin)%matSparse%q(i+1)-
     &          rhogeneral(ispin)%matSparse%q(i)
          endif
        enddo

!        write(12347,*)"rho",nmax
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_rho1',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif




        if(.true.)then
          allocate(rhobuffsparse(nmax))
          do i=1,rhogeneral(ispin)%matSparse%iRows
            ind1=rhogeneral(ispin)%matSparse%q(i)
            ind2=rhogeneral(ispin)%matSparse%q(i+1)
            nsend=ind2-ind1
            call MPI_Barrier(inverseheads_comm,MPIerror)
            CALL MPI_REDUCE(rhogeneral(ispin)%matSparse%b(ind1),
     &          rhobuffsparse(1),
     &          nsend, DAT_dcomplex,MPI_SUM,
     &          0,inverseheads_comm,MPIerror)
            if(myhead .eq. 0)rhogeneral(ispin)%matSparse%b(ind1:ind2-1)=
     &        rhobuffsparse(1:nsend)
          enddo
          deallocate(rhobuffsparse)
        else

          allocate(rhobuffsparse(rhogeneral(ispin)%matSparse%nnz))
          CALL MPI_REDUCE(rhogeneral(ispin)%matSparse%b(1),
     &        rhobuffsparse(1),
     &        rhogeneral(ispin)%matSparse%nnz, DAT_dcomplex,MPI_SUM,
     &        0,inverseheads_comm,MPIerror)
          if(myhead .eq. 0)rhogeneral(ispin)%matSparse%b(:) =
     &          rhobuffsparse(:)
          deallocate(rhobuffsparse)


        endif
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_rho2',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

      ENDDO
#endif
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_rho3',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif

          

      IF (myhead.EQ.0) THEN

        if(bs_add.and.(bs_method.eq.1).and.
     &       allocated(rhobs_general(1)%matSparse%b))then
          do ispin=1,NspinComplexMatrix
            rhogeneral(ispin)%matSparse%b(:)=
     &          rhogeneral(ispin)%matSparse%b(:)+
     &          rhobs_general(ispin)%matSparse%b(:)
          enddo
        endif

        if(bs_add.and.(bs_method.eq.0).and.
     &       allocated(rhobs_delta))then
          do ii=1,n1
            do ind=rhogeneral(1)%matSparse%q(ii),
     &          rhogeneral(1)%matSparse%q(ii+1)-1
              j=rhogeneral(1)%matSparse%j(ind)
              do ispin=1,NspinComplexMatrix
                rhogeneral(ispin)%matSparse%b(ind)=
     &              rhogeneral(ispin)%matSparse%b(ind)+
     &              rhobs_delta(ii,j,ispin)
              enddo
            enddo
          enddo
        endif

      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_rho4',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif

 

      if(myhead.eq.0)then

!        write(12347,*)"rho_reducing3"
        nspinMin=nspin
        if(nspin>1)nspinMin=2
        do ispin=1,NspinComplexMatrix
          rhogeneral(ispin)%matSparse%b(:)=2.D0/REAL(nspinMin)*
     &        rhogeneral(ispin)%matSparse%b(:)
        enddo
!        write(12347,*)"rho_reducing4"
      endif
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_rho5',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


#ifdef MPI
      do ispin=1,NspinComplexMatrix
!        write(12347,*)"rho_reducing5"
        call MPI_Bcast(rhogeneral(ispin)%matSparse%b(1),
     &      rhogeneral(1)%matSparse%nnz,DAT_dcomplex,0,
     &      inverseheads_comm,MPIerror)
!        write(12347,*)"rho_reducing6"
      enddo
#endif
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_rho6',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


      end subroutine rhoG_reduce

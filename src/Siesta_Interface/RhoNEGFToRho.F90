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
!                   RHONEGF_TO_RHOSIESTA  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine rhonegf_to_rhosiesta(maxnelerow,NspinRealInputMatrix,NspinComplexMatrix,rhogeneral,Dnew,maxnh,no,nuo,nuotot,nmat, listh, numh, listhptr,indxuo,kpoint,xij,wk,timereversal)

      use mTypes
      use mMPI_NEGF
#ifdef MPI
      use parallel
#endif

      implicit none

      integer, intent(in) :: maxnh,NspinRealInputMatrix,NspinComplexMatrix,nmat
      integer :: nspinMin
      type(matrixTypeGeneral) :: rhogeneral(NspinComplexMatrix)
      double precision Dnew(maxnh,NspinRealInputMatrix)
      integer maxnelerow
      integer no,nuo,nuotot,listh(maxnh), numh(nuo), listhptr(nuo),indxuo(no)
      double precision kpoint(3),xij(3,maxnh),wk
      logical, intent(in):: timereversal

      double complex, allocatable :: rhorow(:,:)
      integer, allocatable :: jrow(:)
      integer io,BNode,MPIerror,nelerow,ind,iuo,juo,ispin,j,jo,iio,ind2
      integer, allocatable :: ind2ofj(:)
      integer istart,iend,ioprimeUp,ioprimeDown,joprimeUp,joprimeDown,ind2uu,ind2ud,ind2du,ind2dd
#ifdef MPI
      integer istatus(MPI_STATUS_SIZE)
#endif
      double complex  sckxij
      double complex, parameter :: zi=(0.D0,1.D0)
      double precision kxij
      integer*4:: sc_0,sc_1,sc_r,sc_m

!      write(12347,*)"nmat=",nmat,nuotot,NspinComplexMatrix,NspinRealInputMatrix


      nspinMin=NspinComplexMatrix
      if(NspinComplexMatrix>2.and.nmat==nuotot)nspinMin=2
!      write(12347,*)"nspinMin=",nspinMin,nmat,nuotot,NspinComplexMatrix


       allocate(rhorow(maxnelerow,NspinComplexMatrix))
       allocate(jrow(maxnelerow))
       allocate(ind2ofj(2*nmat))
       ind2ofj=0
       do io = 1,nuotot

         call WhichNodeOrb(io,nnodes_negf,BNode)

         if(mynode_negf.eq.BNode)then

           if(nmat==nuotot)then
             istart= rhogeneral(1)%matSparse%q(io)
             iend  = rhogeneral(1)%matSparse%q(io+1)-1
           else
             ioprimeUp=2*io-1
             ioprimeDown=2*io
             istart= rhogeneral(1)%matSparse%q(ioprimeUp)
             iend  = rhogeneral(1)%matSparse%q(ioprimeUp+2)-1
           endif
           nelerow=iend-istart+1
!          write(12347,*)"nelerow=",istart,iend,nelerow,maxnelerow,nmat
           
           do ispin=1,NspinComplexMatrix
             rhorow(1:nelerow,ispin)=rhogeneral(ispin)%matSparse%b(istart:iend)
           enddo
           jrow(1:nelerow)=rhogeneral(1)%matSparse%j(istart:iend)

           if(nmat==nuotot)then
             do ind2=1,nelerow
               ind2ofj(jrow(ind2))=ind2
             enddo
           else
             do ind2=1,nelerow/2
               ind2ofj(jrow(ind2))=ind2
!               write(12347,*)"ind2ofju=",jrow(ind2),jrow(ind2)+nmat,ind2
             enddo
             do ind2=nelerow/2+1,nelerow
               ind2ofj(jrow(ind2)+nmat)=ind2
!               write(12347,*)"ind2ofjd=",jrow(ind2)+nmat,jrow(ind2),ind2
             enddo
           endif

           call GlobalToLocalOrb(io,mynode_negf,nnodes_negf,iio)
           do j = 1,numh(iio)
             
             ind = listhptr(iio) + j
             jo = listh(ind)
             iuo = indxuo(io)
             juo = indxuo(jo)
!             write(12347,*)"ind=",iio,io,iuo,jo,juo,j,ind
!             call inddensetoindsparserow(juo,ind2,jrow,nelerow)
!             write(12347,*)"ind2ofjuo=",iuo,juo,ind2,ind2ofj(juo),ind2-ind2ofj(juo)
             kxij = kpoint(1) * xij(1,ind) + kpoint(2) * xij(2,ind) 
             sckxij = cos(kxij) - zi*sin(kxij)


             if(nmat==nuotot)then
               ind2=ind2ofj(juo)
               if(ind2.eq.0)then
                 write(12347,*)"ind2zero",iuo,juo
               endif

!               write(12347,*)"sckxij=",dreal(sckxij),dimag(sckxij)
               do ispin=1,nspinMin
                 Dnew(ind,ispin)= Dnew(ind,ispin) + wk* dreal(rhorow(ind2,ispin)*sckxij)
               enddo

               if(NspinRealInputMatrix > 2)then
                 if(timereversal)then
                   Dnew(ind,3)= Dnew(ind,3) + wk* dreal(rhorow(ind2,3)*sckxij+DCONJG(rhorow(ind2,4))*DCONJG(sckxij)) * 0.5D0
                   Dnew(ind,4)= Dnew(ind,4) + wk* dimag(rhorow(ind2,3)*sckxij+DCONJG(rhorow(ind2,4))*DCONJG(sckxij)) * 0.5D0
                 else
!this gives t  he same result as for timereversal .true., but it is left separate for clarity for the moment
!the results   for +k and -k are different in general, so we need to use both +k and -k
                   Dnew(ind,3)= Dnew(ind,3) + wk* dreal(rhorow(ind2,3)*sckxij+rhorow(ind2,4)*sckxij) * 0.5D0
                   Dnew(ind,4)= Dnew(ind,4) + wk* dimag(rhorow(ind2,3)*sckxij-rhorow(ind2,4)*sckxij) * 0.5D0
                 endif

                 if(NspinRealInputMatrix > 4)then
                   Dnew(ind,5)= Dnew(ind,5) + wk* dimag(rhorow(ind2,1)*sckxij)
                   Dnew(ind,6)= Dnew(ind,6) + wk* dimag(rhorow(ind2,2)*sckxij)
                 endif
               endif
             else
               joprimeUp=2*juo-1
               joprimeDown=2*juo
               ind2uu=ind2ofj(joprimeUp)
               ind2ud=ind2ofj(joprimeUp+1)
               ind2du=ind2ofj(joprimeUp+nmat)
               ind2dd=ind2ofj(joprimeUp+nmat+1)

               Dnew(ind,1)= Dnew(ind,1) + wk* dreal(rhorow(ind2uu,1)*sckxij)
               Dnew(ind,2)= Dnew(ind,2) + wk* dreal(rhorow(ind2dd,1)*sckxij)

               if(timereversal)then
                 Dnew(ind,3)= Dnew(ind,3) + wk* dreal(rhorow(ind2ud,1)*sckxij+DCONJG(rhorow(ind2du,1))*DCONJG(sckxij)) * 0.5D0
                 Dnew(ind,4)= Dnew(ind,4) + wk* dimag(rhorow(ind2ud,1)*sckxij+DCONJG(rhorow(ind2du,1))*DCONJG(sckxij)) * 0.5D0
               else
                 Dnew(ind,3)= Dnew(ind,3) + wk* dreal(rhorow(ind2ud,1)*sckxij+rhorow(ind2du,1)*sckxij) * 0.5D0
                 Dnew(ind,4)= Dnew(ind,4) + wk* dimag(rhorow(ind2ud,1)*sckxij-rhorow(ind2du,1)*sckxij) * 0.5D0
               endif

               if(NspinRealInputMatrix > 4)then
                 Dnew(ind,5)= Dnew(ind,5) + wk* dimag(rhorow(ind2uu,1)*sckxij)
                 Dnew(ind,6)= Dnew(ind,6) + wk* dimag(rhorow(ind2dd,1)*sckxij)
               endif
             endif
             
!             write(12347,*)"dnew=",iuo,juo,dnew(ind,1),dnew(ind,2),dnew(ind,3),dnew(ind,4),dnew(ind,5),dnew(ind,6)

           enddo

           do ind2=1,nelerow
             ind2ofj(jrow(ind2))=0
           enddo

        endif
!#ifdef MPI
!        CALL MPI_BARRIER(negf_comm, MPIerror)
!#endif
       enddo
       deallocate(rhorow,jrow,ind2ofj)

!       call writematreal(Dnew,maxnh,NspinComplexMatrix,"dsck")
!       call MPI_Finalize( MPIerror )
!       stop

 end subroutine rhonegf_to_rhosiesta






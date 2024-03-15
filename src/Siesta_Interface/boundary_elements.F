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
!                   BOUNDARY_ELEMENTS_MPI  
! AND
! THE MODULE
!                   SET_RHOBD  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
      MODULE  set_rhobd

      use negfmod
      use mMPI_NEGF
#ifdef MPI
      use parallel
#endif
    
      implicit none

      CONTAINS

      subroutine boundary_elements_mpi(neql,maxnhl,iequivL,maxnhg,
     .    listhg,listig,listhptrg,nuotot,dnew,DMbulkL,
     .    NspinRealInputMatrix,
     .    maxnh,numh,nuo,listhptr,listh,DMlead,iBias,flag)
 
      implicit none
      integer neql,maxnhl,nuotot,iequivL(maxnhL,2),maxnhg,nuo,
     .    listhg(maxnhg),listig(maxnhg),listhptrg(nuotot),
     .    NspinRealInputMatrix,maxnh,
     .    numh(nuo),listhptr(nuo),listh(maxnh)
      double precision DMbulkL(maxnhL,NspinRealInputMatrix),
     .    dnew(maxnh,NspinRealInputMatrix),dnewele
      integer ind,mpierror,iuo,jo,bnode,i,j,ind2,jo2,iio,
     .    ispin
!zrx,
      integer :: iopt
      double precision, optional :: DMlead(maxnhL,NspinRealInputMatrix),
     .    iBias
      character(len=*), optional :: flag
!zrx,


#ifdef MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

  
!zrx,
      if(present(flag)) then
        if(flag.eq.'L') then
          iopt=2
        else if(flag.eq.'R') then
          iopt=3
        endif
      else
         iopt=1
      endif 



      do ind=1,NeqL

!       call findij(iequivL(ind,1),nuotot,maxnhg,listhg,listhptrg,iuo,jo)
!       write(12347,*)"findij=",iuo,jo,listig(iequivL(ind,1)),
!     .     listhg(iequivL(ind,1))
       iuo=listig(iequivL(ind,1))
       jo=listhg(iequivL(ind,1))
       call WhichNodeOrb(iuo,nnodes_negfo,BNode)
       if(mynode_negfo.eq.BNode)then
          call GlobalToLocalOrb(iuo,mynode_negfo,nnodes_negfo,iio)
          do j = 1,numh(iio)
            ind2 = listhptr(iio) + j
            jo2 = listh(ind2)
            if(jo2.eq.jo)exit
          enddo
          IF(iopt.eq.1) THEN
             Dnew(ind2,:)=DMbulkL(iequivL(ind,2),:)
          ELSE if(iopt.eq.2) then
             Dnew(ind2,:)=DMbulkL(iequivL(ind,2),:)+
     .   (iBias/2.d0)*DMlead(iequivL(ind,2),:)
          ELSE if(iopt.eq.3) then
             Dnew(ind2,:)=DMbulkL(iequivL(ind,2),:)-
     .   (iBias/2.d0)*DMlead(iequivL(ind,2),:)
          ENDIF
        endif
      enddo
!zrx,

      end subroutine boundary_elements_mpi

      END MODULE set_rhobd

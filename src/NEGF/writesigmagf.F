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
!                   WRITESIGMAINFO,
!                   WRITESIGMA  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

      subroutine writesigmainfo(slabel,n1,nl,nr,Nnodes,NeneT_div,
     &    NSPIN,nk,ei,tenergi,deltaenebig,deltaenes,v,ef)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************


      integer    iufilesigma,n1,nl,nr,Nnodes,NeneT_div,NSPIN,nk,nodei,i
      double precision tenergi,deltaenebig,deltaenes,v,ef
      double complex  EiCompR, EiCompL,ei
      CHARACTER(LEN=*) :: slabel
      CHARACTER :: filesigmainfo*55
      
        
      filesigmainfo=TRIM(ADJUSTL(slabel)) // ".TRC.sigmainfo"

      call io_assign(iufilesigma)
      OPEN(UNIT=iufilesigma,FILE=filesigmainfo,RECL=1000000,
     &     status='unknown')

      write(iufilesigma,*)"# n1,nl,nr,number of nodes, energy points per
     & node, number of spins, number of k-points"
      write(iufilesigma,*)n1,nl,nr,Nnodes,NeneT_div,NSPIN,nk
      write(iufilesigma,*)"#node index, local energy index, energy,
     &    energy of left selfenergy, energy of right selfenergy"
      do nodei=0,Nnodes-1
        DO I=1,NeneT_div
          IF (NeneT.EQ.1) THEN
            Ei=TEnergI
          ELSE 
            Ei=TEnergI + (i-1)*deltaenebig + nodei * deltaenes
          ENDIF
          EiCompL=Ei-V/2.0D0
          EiCompR=Ei+V/2.0D0
          write(iufilesigma,*)nodei,i,(ei-ef)*
     &        13.6056981D0,(DREAL(EiCompL)-ef)*
     &        13.6056981D0,(DREAL(EiCompR)-ef)*13.6056981D0
        enddo
      enddo
      call io_close(iufilesigma)

      end subroutine writesigmainfo


      subroutine writesigma(slabel,ispin,mynode,ik,i,n1,nl,nr,Ei,ef,
     &    sigmal,sigmar,gf_iter)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      CHARACTER(LEN=*) :: slabel
      CHARACTER :: filetrksigmal*55
      character*20   charnode
      integer ISPIN,mynode,iufilesigma,ik,i,i1,i2,nl,nr,n1
      double precision ef
      double complex ei
      double complex sigmal(nl,nl),sigmar(nr,nr),gf_iter(n1,n1)

      write(charnode,'(i10)')MyNode
      IF (ISPIN .EQ. 1)
     &  filetrksigmal=TRIM(ADJUSTL(slabel)) // ".TRC.sigmal.k.up_"
      IF (ISPIN .EQ. 2)
     &  filetrksigmal=TRIM(ADJUSTL(slabel)) // ".TRC.sigmal.k.down_"

      filetrksigmal=TRIM(ADJUSTL(filetrksigmal)) // 
     &      TRIM(ADJUSTL(charnode))
 
!      write(*,*)"writesigma,outfilename",TRIM(filetrksigmal),MyNode
      call io_assign(iufilesigma)
      IF ((ik .EQ. 1).AND.(I.EQ.1)) THEN
       OPEN(UNIT=iufilesigma,FILE=filetrksigmal,RECL=1000000,
     &     status='unknown')
      ELSE
       OPEN(UNIT=iufilesigma,FILE=filetrksigmal,RECL=1000000,
     &     POSITION='append')
      ENDIF

      do i1=1,nL
        do i2=1,nL
          WRITE(iufilesigma,*)ispin,ik,(Ei-ef)* 13.6056981D0,
     &        DREAL(SigmaL(i1,i2)),
     &        DIMAG(SigmaL(i1,i2))
        enddo
      enddo

      do i1=1,nr
        do i2=1,nr
          WRITE(iufilesigma,*)ispin,ik,(Ei-ef)* 13.6056981D0,
     &        DREAL(SigmaR(i1,i2)),
     &        DIMAG(SigmaR(i1,i2))
        enddo
      enddo

      do i1=1,n1
        do i2=1,n1
          WRITE(iufilesigma,*)ispin,ik,(Ei-ef)* 13.6056981D0,
     &        DREAL(GF_iter(i1,i2)),
     &        DIMAG(GF_iter(i1,i2))
        enddo
      enddo


      call io_close(iufilesigma)

      end subroutine writesigma

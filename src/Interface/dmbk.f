 
! Copyright (c) Smeagol Authors:
! A. R. Rocha, V. Garcia-Suarez, S. Bailey, C. J. Lambert, J. Ferrer and
! S. Sanvito 2003-2005
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
! SMEAGOL IS DISTRIBUTED ONLY THROUGH THE OFICIAL WEBSITE (www.smeagol.tcd.ie)
! UPON COMPLETION OF THE "SMEAGOL ACADEMIC LICENSE" .
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
	SUBROUTINE DMBK(Nbasis,ndmax,nspin,slabel,
     &   numd,listd,listdptr,DM, iopt)
C *********************************************************************
C Subroutine to read a saved density matrix for the left or right leads
C
C Written by Alexandre Reily Rocha, June 2003
C Computational Spintronics Group
C Trinity College Dublin
C e-mail: rochaa@tcd.ie
C ***************************** HISTORY ***********************************
C Original version:	June 2003
C ********************* INPUTS ****************************************	
C INTEGER NBasis		: Number of basis orbitals
C INTEGER ndmax			: Maximum number of non-zero elements of
C				  the density matrix
C INTEGER nspin			: Number of spin components
C CHAR*20 slabel		: Name of the file (.DM) to be read
C INTEGER numd(Nbasis)		: Number of non-zero elements per orbital
C INTEGER listd(ndmax)		: Nonzero Hamiltonian-matrix element  
C                          	  column indexes for each matrix row
C INTEGER listdptr(Nbasis)	: Pointer to each row (-1) of the
C                          	  Hamiltonian matrix
C ******************** OUTPUTS ****************************************
C DOUBLE DM(ndmax,nspin)	: Leads density matrix
C *********************************************************************

	
      use ionew
	IMPLICIT NONE
	
	INTEGER :: nb,ns,Nbasis,ndmax,nspin
	DOUBLE PRECISION:: DM(ndmax,nspin)
	INTEGER :: listd(ndmax)
	INTEGER, DIMENSION (Nbasis) :: numd,listdptr,numdg
	INTEGER :: m,is,ml,ndmaxg,im
	CHARACTER slabel*20, paste*25
	EXTERNAL paste
        integer iounit
!zrx,
       integer :: iopt

        call io_assign(iounit)
        if(iopt.eq.1) then
          OPEN(UNIT=iounit,file=paste(slabel,'.DM'),FORM='unformatted',
     &            STATUS='OLD')
        else if(iopt.eq.2) then
          OPEN(UNIT=iounit,file=paste(slabel,'.EDM'),FORM='unformatted',
     &            STATUS='OLD')
        endif
        REWIND(iounit)
        READ(iounit) nb,ns


        READ(iounit) (numdg(m),m=1,Nbasis)

        DO m = 1,Nbasis
	 ml=m
         numd(m) = numdg(ml)
         IF (m .eq. 1) then
          listdptr(1) = 0
         ELSE
          listdptr(m) = listdptr(m-1) + numd(m-1)
         ENDIF
        ENDDO
        ndmaxg = 0
        DO m = 1,Nbasis
         ndmaxg = max(ndmaxg,numdg(m))
        ENDDO

        DO m = 1,Nbasis
         ml = m
         READ(iounit) (listd(listdptr(ml)+im),im=1,numd(ml))
        ENDDO

        DO is = 1,nspin
         DO m = 1,Nbasis
          ml = m
          READ(iounit) (DM(listdptr(ml)+im,is),im=1,numd(ml))
         ENDDO
        ENDDO

        call io_close(iounit)

	END SUBROUTINE DMBK

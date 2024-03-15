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
!                   SET_BOUNDARYELEMENTSON  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
SUBROUTINE set_boundaryelementson(n_replace_L,n_replace_R,NSPIN,V,N1,NL,H0_L,S0_L,H1_L,S1_L, NR,H0_R,S0_R,H1_R,S1_R,hgeneral,sgeneral,Set_HBoundary_Leads,Set_HLR_Zero,tol)


! ********************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:	October 2007
! ********************************************************************
  use negfmod, only : nlSetZeroRatio,nrSetZeroRatio, outinfo
  use mTypes

  IMPLICIT NONE
  INTEGER,intent(in) :: n_replace_L,n_replace_R
  INTEGER         NSPIN,N1,NL,NR,ISPIN,I,NSLICES,II,JJ,ind,i1,i2
  DOUBLE COMPLEX  H0_L(NL,NL,NSPIN),H1_L(NL,NL,NSPIN)
  DOUBLE COMPLEX  S0_L(NL,NL),S1_L(NL,NL)
  DOUBLE COMPLEX  H0_R(NR,NR,NSPIN),H1_R(NR,NR,NSPIN)
  DOUBLE COMPLEX  S0_R(NR,NR),S1_R(NR,NR)
  DOUBLE PRECISION V
  type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
  integer nspinMin
  integer ii0,jj0,istart,iend,istart_1,iend_1,istart_2,iend_2
  integer nlSetZero,nrSetZero

  logical,intent(in)::Set_HBoundary_Leads,Set_HLR_Zero
  double precision, intent(in) :: tol
 
  logical debugoutput

  debugoutput=.false.


  nlSetZero=nlSetZeroRatio * nl
  nrSetZero=nrSetZeroRatio * nr
  if(nlSetZeroRatio.ne.1.0D0.or.nrSetZeroRatio.ne.1.0D0 .and. outinfo) write(12347,*)"nlSetZero,nrSetZero=",nlSetZero,nrSetZero,nl,nr,nlSetZeroRatio,nrSetZeroRatio

  IF ((NL .NE. N1) .AND. Set_HLR_Zero) THEN 
    do ii=1,nlSetZero
      do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
        if(sgeneral%matSparse%j(ind).gt.(n1-nrSetZero))then
          do ispin=1,nspin
            hgeneral(ispin)%matSparse%b(ind)=0D0
          enddo
          sgeneral%matSparse%b(ind)=0D0
        endif
      enddo
    enddo
  ENDIF

  IF ((NR .NE. N1) .AND. Set_HLR_Zero) THEN
    do ii=n1-nrSetZero+1,n1
      do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
        if(sgeneral%matSparse%j(ind).le.nlSetZero)then
          do ispin=1,nspin
            hgeneral(ispin)%matSparse%b(ind)=0D0
          enddo
          sgeneral%matSparse%b(ind)=0D0
        endif
      enddo
    enddo
  ENDIF


  nspinMin=nspin
  if(nspin>2)nspinMin=2

  if(Set_HBoundary_Leads)then

    do i1=0,n_replace_L-1
      if(debugoutput .and. outinfo) write(12347,*)"replace_0L",i1

      istart=i1*nl
      iend=(i1+1)*nl

      do ii=istart+1,iend
        do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
          if((sgeneral%matSparse%j(ind).gt.istart).and.(sgeneral%matSparse%j(ind).le.iend))then
            ii0=ii-istart
            jj0=sgeneral%matSparse%j(ind)-istart
 
            do ispin=1,nspinMin
              hgeneral(ispin)%matSparse%b(ind)=H0_L(ii0,jj0,ispin)+V/2.D0*S0_L(ii0,jj0)
            enddo
            if(nspin>2) hgeneral(3)%matSparse%b(ind)=H0_L(ii0,jj0,3)
 
            sgeneral%matSparse%b(ind)=S0_L(ii0,jj0)
          endif
        enddo
      enddo

    enddo

    do i1=0,n_replace_R-1
      if(debugoutput .and. outinfo) write(12347,*)"replace_0R",i1

      istart=n1-(i1+1)*nr
      iend=n1-i1*nr

      do ii=istart+1,iend
        do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
          if((sgeneral%matSparse%j(ind).gt.istart).and.(sgeneral%matSparse%j(ind).le.iend))then
            ii0=ii-istart
            jj0=sgeneral%matSparse%j(ind)-istart
 
            do ispin=1,nspinMin
              hgeneral(ispin)%matSparse%b(ind)=H0_R(ii0,jj0,ispin)-V/2.D0*S0_R(ii0,jj0)
            enddo
            if(nspin>2) hgeneral(3)%matSparse%b(ind)=H0_R(ii0,jj0,3)
 
            sgeneral%matSparse%b(ind)=S0_R(ii0,jj0)
          endif
        enddo
      enddo

    enddo

    do i1=1,n_replace_L-1
      if(debugoutput .and. outinfo) write(12347,*)"replace_1L",i1

      istart_1=(i1-1)*nl
      iend_1=i1*nl
      istart_2=i1*nl
      iend_2=(i1+1)*nl

      do ii=istart_1+1,iend_1
        do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
          if((sgeneral%matSparse%j(ind).gt.istart_2).and.(sgeneral%matSparse%j(ind).le.iend_2))then
            ii0=ii-istart_1
            jj0=sgeneral%matSparse%j(ind)-istart_2
 
            do ispin=1,nspinMin
              hgeneral(ispin)%matSparse%b(ind)=H1_L(ii0,jj0,ispin)+V/2.D0*S1_L(ii0,jj0)
            enddo
            if(nspin>2) hgeneral(3)%matSparse%b(ind)=H1_L(ii0,jj0,3)
 
            sgeneral%matSparse%b(ind)=S1_L(ii0,jj0)
          endif
        enddo
      enddo

    enddo

    do i1=1,n_replace_L-1
      if(debugoutput .and. outinfo) write(12347,*)"replace_m1L",i1
      istart_1=i1*nl
      iend_1=(i1+1)*nl
      istart_2=(i1-1)*nl
      iend_2=i1*nl

      do ii=istart_1+1,iend_1
        do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
          if((sgeneral%matSparse%j(ind).gt.istart_2).and.(sgeneral%matSparse%j(ind).le.iend_2))then
            ii0=ii-istart_1
            jj0=sgeneral%matSparse%j(ind)-istart_2
 
            do ispin=1,nspinMin
              hgeneral(ispin)%matSparse%b(ind)=DCONJG(H1_L(jj0,ii0,ispin))+V/2.D0*DCONJG(S1_L(jj0,ii0))
            enddo
            if(nspin>2) hgeneral(3)%matSparse%b(ind)=DCONJG(H1_L(jj0,ii0,4))
 
            sgeneral%matSparse%b(ind)=DCONJG(S1_L(jj0,ii0))
          endif
        enddo
      enddo





    enddo



    do i1=1,n_replace_R-1
      if(debugoutput .and. outinfo) write(12347,*)"replace_1R",i1

      istart_1=n1-(i1+1)*nr
      iend_1=n1-i1*nr
      istart_2=n1-i1*nr
      iend_2=n1-(i1-1)*nr

      do ii=istart_1+1,iend_1
        do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
          if((sgeneral%matSparse%j(ind).gt.istart_2).and.(sgeneral%matSparse%j(ind).le.iend_2))then
            ii0=ii-istart_1
            jj0=sgeneral%matSparse%j(ind)-istart_2
 
            do ispin=1,nspinMin
              hgeneral(ispin)%matSparse%b(ind)=H1_R(ii0,jj0,ispin)+V/2.D0*S1_R(ii0,jj0)
            enddo
            if(nspin>2) hgeneral(3)%matSparse%b(ind)=H1_R(ii0,jj0,3)
 
            sgeneral%matSparse%b(ind)=S1_R(ii0,jj0)
          endif
        enddo
      enddo

    enddo


    do i1=1,n_replace_R-1
      if(debugoutput .and. outinfo) write(12347,*)"replace_m1R",i1

      istart_1=n1-i1*nr
      iend_1=n1-(i1-1)*nr
      istart_2=n1-(i1+1)*nr
      iend_2=n1-i1*nr

      do ii=istart_1+1,iend_1
        do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
          if((sgeneral%matSparse%j(ind).gt.istart_2).and.(sgeneral%matSparse%j(ind).le.iend_2))then
            ii0=ii-istart_1
            jj0=sgeneral%matSparse%j(ind)-istart_2
 
            do ispin=1,nspinMin
              hgeneral(ispin)%matSparse%b(ind)=DCONJG(H1_R(jj0,ii0,ispin))+V/2.D0*DCONJG(S1_R(jj0,ii0))
            enddo
            if(nspin>2) hgeneral(3)%matSparse%b(ind)=DCONJG(H1_R(jj0,ii0,4))
 
            sgeneral%matSparse%b(ind)=DCONJG(S1_R(jj0,ii0))
          endif
        enddo
      enddo

    enddo





  else
      if (outinfo) write(12347,*)"not setting boundary elements of H from leads"
  endif


  DO ind=1,sgeneral%matSparse%nnz
    DO ISPIN=1,NSPIN
      if(abs(hgeneral(ISPIN)%matSparse%b(ind)).lt.tol) hgeneral(ispin)%matSparse%b(ind)=0D0
    ENDDO
    if(abs(sgeneral%matSparse%b(ind)).lt.tol) sgeneral%matSparse%b(ind)=0D0
  ENDDO


end SUBROUTINE set_boundaryelementson

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
!                   SET_BOUNDARYELEMENTSON2  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
SUBROUTINE set_boundaryelementson2(n_replace_L,n_replace_R,V,N1,NL,NR,hgeneral,sgeneral,Set_HBoundary_Leads,Set_HLR_Zero,tol)


! ********************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:	October 2007
! ********************************************************************
  use sigma, only:h0_l,s0_l,h1_l,s1_l,h0_r,s0_r,h1_r,s1_r
  use negfmod, only : nlSetZeroRatio,nrSetZeroRatio, outinfo
  use mTypes

  IMPLICIT NONE
  INTEGER,intent(in) :: n_replace_L,n_replace_R
  INTEGER         N1,NL,NR,I,NSLICES,II,JJ,ind,i1,i2
  DOUBLE PRECISION V
  type(matrixTypeGeneral) :: hgeneral,sgeneral
  integer ii0,jj0,istart,iend,istart_1,iend_1,istart_2,iend_2

  logical,intent(in)::Set_HBoundary_Leads,Set_HLR_Zero
  double precision, intent(in) :: tol
  integer nlSetZero,nrSetZero
 
  logical debugoutput

  debugoutput=.false.

  if (outinfo) write(12347,*)"nl,nr=",nl,nr,n1

  nlSetZero=nlSetZeroRatio * nl
  nrSetZero=nrSetZeroRatio * nr
  if((nlSetZeroRatio.ne.1.0D0.or.nrSetZeroRatio.ne.1.0D0) .and. outinfo) write(12347,*)"nlSetZero,nrSetZero=",nlSetZero,nrSetZero,nl,nr,nlSetZeroRatio,nrSetZeroRatio

  IF ((NL .NE. N1) .AND. Set_HLR_Zero) THEN 
    do ii=1,nlSetZero
      do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
        if(sgeneral%matSparse%j(ind).gt.(n1-nrSetZero))then
          hgeneral%matSparse%b(ind)=0D0
          sgeneral%matSparse%b(ind)=0D0
        endif
      enddo
    enddo
  ENDIF

  IF ((NR .NE. N1) .AND. Set_HLR_Zero) THEN
    do ii=n1-nrSetZero+1,n1
      do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
        if(sgeneral%matSparse%j(ind).le.nlSetZero)then
          hgeneral%matSparse%b(ind)=0D0
          sgeneral%matSparse%b(ind)=0D0
        endif
      enddo
    enddo
  ENDIF



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
 
            if(mod(ii0,2)==mod(jj0,2))then
              hgeneral%matSparse%b(ind)=H0_L(ii0,jj0,1)+V/2.D0*S0_L(ii0,jj0)
            else
              hgeneral%matSparse%b(ind)=H0_L(ii0,jj0,1)
            endif

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

            if(mod(ii0,2)==mod(jj0,2))then
              hgeneral%matSparse%b(ind)=H0_R(ii0,jj0,1)-V/2.D0*S0_R(ii0,jj0)
            else
              hgeneral%matSparse%b(ind)=H0_R(ii0,jj0,1)
            endif

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
  
            if(mod(ii0,2)==mod(jj0,2))then
              hgeneral%matSparse%b(ind)=H1_L(ii0,jj0,1)+V/2.D0*S1_L(ii0,jj0)
            else
              hgeneral%matSparse%b(ind)=H1_L(ii0,jj0,1)
            endif

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
 
            if(mod(ii0,2)==mod(jj0,2))then
              hgeneral%matSparse%b(ind)=H1_L(ii0,jj0,1)+V/2.D0*S1_L(ii0,jj0)
            else
              hgeneral%matSparse%b(ind)=H1_L(ii0,jj0,1)
            endif

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

            if(mod(ii0,2)==mod(jj0,2))then
              hgeneral%matSparse%b(ind)=H1_r(ii0,jj0,1)-V/2.D0*S1_r(ii0,jj0)
            else
              hgeneral%matSparse%b(ind)=H1_r(ii0,jj0,1)
            endif

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
 
            if(mod(ii0,2)==mod(jj0,2))then
              hgeneral%matSparse%b(ind)=H1_r(ii0,jj0,1)-V/2.D0*S1_r(ii0,jj0)
            else
              hgeneral%matSparse%b(ind)=H1_r(ii0,jj0,1)
            endif
 
            sgeneral%matSparse%b(ind)=DCONJG(S1_R(jj0,ii0))
          endif
        enddo
      enddo

    enddo

  else
      if (outinfo) write(12347,*)"not setting boundary elements of H from leads"
  endif


  DO ind=1,sgeneral%matSparse%nnz
    if(abs(hgeneral%matSparse%b(ind)).lt.tol)then 
      hgeneral%matSparse%b(ind)=0D0
    endif
    if(abs(sgeneral%matSparse%b(ind)).lt.tol) sgeneral%matSparse%b(ind)=0D0  
  ENDDO


end SUBROUTINE set_boundaryelementson2

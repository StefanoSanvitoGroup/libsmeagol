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
!                   SET_BOUNDARYELEMENTSGENERAL,
!                   SET_BOUNDARYELEMENTSONP  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!


SUBROUTINE set_boundaryelementsgeneral(n_replace_L,n_replace_R,NSPIN,NspinComplexMatrix, V,N1,NL,H0_L,S0_L,H1_L,S1_L, NR,H0_R,S0_R,H1_R,S1_R,hgeneral,sgeneral,Set_HBoundary_Leads,Set_HLR_Zero,tol)


! ********************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:	October 2007
! ********************************************************************
  use mTypes

  IMPLICIT NONE
  INTEGER,intent(in) :: n_replace_L,n_replace_R,NspinComplexMatrix
  INTEGER         NSPIN,N1,NL,NR,ISPIN,I,II,JJ,ind
  DOUBLE COMPLEX  H0_L(NL,NL,NSPIN),H1_L(NL,NL,NSPIN)
  DOUBLE COMPLEX  S0_L(NL,NL),S1_L(NL,NL)
  DOUBLE COMPLEX  H0_R(NR,NR,NSPIN),H1_R(NR,NR,NSPIN)
  DOUBLE COMPLEX  S0_R(NR,NR),S1_R(NR,NR)
  DOUBLE PRECISION V
  type(matrixTypeGeneral) :: hgeneral(NspinComplexMatrix),sgeneral
  double precision, intent(in) :: tol

  logical,intent(in)::Set_HBoundary_Leads,Set_HLR_Zero

  if(sgeneral%mattype==2)then
    if(nspin.eq.NspinComplexMatrix)then
!      write(12347,*)"calling set_boundaryelementson"
      call set_boundaryelementson(n_replace_L,n_replace_R,NSPIN,V,N1,NL,H0_L,S0_L,H1_L,S1_L, NR,H0_R,S0_R,H1_R,S1_R,hgeneral,sgeneral,Set_HBoundary_Leads,Set_HLR_Zero,tol)
    else
!      write(12347,*)"calling set_boundaryelementson2"
      call set_boundaryelementson2(n_replace_L,n_replace_R, V,2*N1,2*NL,2* NR,hgeneral(1),sgeneral,Set_HBoundary_Leads,Set_HLR_Zero,tol)
    endif
  elseif(sgeneral%mattype==3)then
!    write(12347,*)"calling set_boundaryelementsonP"
    call set_boundaryelementsonP(n_replace_L,NSPIN, V,N1,NL,H0_L,S0_L,H1_L,S1_L, NR,H0_R,S0_R,H1_R,S1_R,hgeneral,sgeneral,Set_HBoundary_Leads,Set_HLR_Zero,tol)
  endif


end SUBROUTINE set_boundaryelementsgeneral


SUBROUTINE set_boundaryelementsonP(NSLICES,NSPIN, V,N1,NL,H0_L,S0_L,H1_L,S1_L, NR,H0_R,S0_R,H1_R,S1_R,hgeneral,sgeneral,Set_HBoundary_Leads,Set_HLR_Zero,tol)


! ********************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:	October 2007
! ********************************************************************
  use mTypes

  IMPLICIT NONE
  INTEGER         NSPIN,N1,NL,NR,ISPIN,I,NSLICES,II,JJ,ind
  DOUBLE COMPLEX  H0_L(NL,NL,NSPIN),H1_L(NL,NL,NSPIN)
  DOUBLE COMPLEX  S0_L(NL,NL),S1_L(NL,NL)
  DOUBLE COMPLEX  H0_R(NR,NR,NSPIN),H1_R(NR,NR,NSPIN)
  DOUBLE COMPLEX  S0_R(NR,NR),S1_R(NR,NR)
  DOUBLE PRECISION V
  type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral

  logical,intent(in)::Set_HBoundary_Leads,Set_HLR_Zero
  double precision, intent(in) :: tol
  integer         iend,istart


  istart=sgeneral%matSparseP%matSparse%iVert
  iend=istart+sgeneral%matSparseP%matSparse%iRows-1

  IF ((NL .NE. N1) .AND. Set_HLR_Zero) THEN 
    do ii=1,nl
      if(ii<=iend.and.ii>=istart)then
        do ind=sgeneral%matSparseP%matSparse%q(ii-istart+1),sgeneral%matSparseP%matSparse%q(ii-istart+2)-1
          if(sgeneral%matSparseP%matSparse%j(ind).gt.(n1-nl))then
            do ispin=1,nspin
              hgeneral(ispin)%matSparseP%matSparse%b(ind)=0D0
            enddo
            sgeneral%matSparseP%matSparse%b(ind)=0D0
          endif
        enddo
      endif
    enddo
  ENDIF

  IF ((NR .NE. N1) .AND. Set_HLR_Zero) THEN
    do ii=n1-nr+1,n1
      if(ii<=iend.and.ii>=istart)then
        do ind=sgeneral%matSparseP%matSparse%q(ii-istart+1),sgeneral%matSparseP%matSparse%q(ii-istart+2)-1
          if(sgeneral%matSparseP%matSparse%j(ind).le.nr)then
            do ispin=1,nspin
              hgeneral(ispin)%matSparseP%matSparse%b(ind)=0D0
            enddo
            sgeneral%matSparseP%matSparse%b(ind)=0D0
          endif
        enddo
      endif
    enddo
  ENDIF



  if(Set_HBoundary_Leads)then
    do ii=1,nl
      if(ii<=iend.and.ii>=istart)then
        do ind=sgeneral%matSparseP%matSparse%q(ii-istart+1),sgeneral%matSparseP%matSparse%q(ii-istart+2)-1
          if(sgeneral%matSparseP%matSparse%j(ind).le.nl)then

            do ispin=1,nspin
              hgeneral(ispin)%matSparseP%matSparse%b(ind)=H0_L(ii,sgeneral%matSparseP%matSparse%j(ind),ispin)+V/2.D0*S0_L(ii,sgeneral%matSparseP%matSparse%j(ind))
            enddo
            sgeneral%matSparseP%matSparse%b(ind)=S0_L(ii,sgeneral%matSparseP%matSparse%j(ind))
          endif
        enddo
      endif
    enddo

         
    do ii=n1-nr+1,n1
      if(ii<=iend.and.ii>=istart)then
        do ind=sgeneral%matSparseP%matSparse%q(ii-istart+1),sgeneral%matSparseP%matSparse%q(ii-istart+2)-1
          if(sgeneral%matSparseP%matSparse%j(ind).ge.n1-nr+1)then
            do ispin=1,nspin
              hgeneral(ispin)%matSparseP%matSparse%b(ind)=H0_R(ii-n1+nr,sgeneral%matSparseP%matSparse%j(ind)-n1+nr,ispin)-V/2.D0*S0_R(ii-n1+nr,sgeneral%matSparseP%matSparse%j(ind)-n1+nr)
            enddo
            sgeneral%matSparseP%matSparse%b(ind)=S0_R(ii-n1+nr,sgeneral%matSparseP%matSparse%j(ind)-n1+nr)
          endif
        enddo
      endif
    enddo
  endif


  DO ind=1,sgeneral%matSparseP%matSparse%nnz
    if(abs(hgeneral(1)%matSparseP%matSparse%b(ind)).lt.tol)then
      DO ISPIN=1,NSPIN
        hgeneral(ispin)%matSparseP%matSparse%b(ind)=0D0
      ENDDO
      sgeneral%matSparseP%matSparse%b(ind)=0D0
    endif
  ENDDO



end SUBROUTINE set_boundaryelementsonP


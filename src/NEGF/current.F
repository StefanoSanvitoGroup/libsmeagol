! 
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
	subroutine current(I,nk,N1,NL,NR,NSPIN,ISPIN,
     $      ik,Ei,
     &      V,T,Ef_Lead,
     &      CONST,Delta,
     $      wk,Ic,tk,hgeneral,sgeneral)

! *****************************************************************
! Calculates the current after a self consistent
! DFT calculation of a nanoscopic system coupled to
! charge reservoirs. 
!
! Written by Alexandre Reily Rocha and Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: rochaa@tcd.ie, reilya@if.usp.br, runggeri@tcd.ie
! *****************************************************************

      use negfmod, only: gamma_negf
      use mTypes
      use mEnergyGrid

	IMPLICIT NONE
	
	include "const2.h"
	
	INTEGER NSPIN,nk,NL,NR,ISPIN,N1,II,
     .   JJ,I,INFO,J,L,ik,ind,il
	
	double precision :: Ei,Ef_Lead,V,T,Delta,wk
	
	DOUBLE COMPLEX :: fL,fR,CONST

	INTEGER, DIMENSION (N1) :: IPIV

	DOUBLE PRECISION, DIMENSION (NSPIN) :: Ic,tk
	
	DOUBLE COMPLEX, DIMENSION (NL) :: Ic_aux
	
	DOUBLE COMPLEX, DIMENSION (N1,N1) :: GF_iter,GF_iter_dag
	
	
	
	DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: 
     .   Gamma1_aux,Gamma2_aux,Gamma_L,Gamma_R
	
	DOUBLE COMPLEX, DIMENSION (N1*N1) :: WORK
      type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
      type(SelfEnergyType), allocatable :: sigmaleads(:)
	
		
	allocate(Gamma1_aux(NL,NR),Gamma2_aux(NR,NL),
     &   Gamma_L(NL,NL),Gamma_R(NR,NR))
		
	IF ((Ei-Ef_Lead-e*V/2.D0)/(kB*T) .GT. 40.D0) THEN
	 fL=0.D0
	ELSE
	 fL=1.D0/(1+DEXP((Ei-Ef_Lead-e*V/2.D0)/(kB*T)))
	ENDIF

	IF ((Ei-Ef_Lead+e*V/2.D0)/(kB*T) .GT. 40.D0) THEN
	 fR=0.D0
	ELSE
	 fR=1.D0/(1+DEXP((Ei-Ef_Lead+e*V/2.D0)/(kB*T)))
	ENDIF

      tk(ispin)=0D0
!	DO ISPIN=1,NSPIN
       GF_iter=0D0
       do ii=1,n1
         do ind=sgeneral%matSparse%q(ii),
     .       sgeneral%matSparse%q(ii+1)-1
           GF_iter(ii,sgeneral%matSparse%j(ind))=
     .         (Ei+zi*Delta)*sgeneral%matSparse%b(ind)-
     .         hgeneral(ispin)%matSparse%b(ind)
         enddo
       enddo

!       GF_iter=(Ei+zi*Delta)*S0_Chain-
!     &    H_Chain(:,:,ISPIN)

         allocate(sigmaleads(ERealGrid%nLeads))
         do il=1,ERealGrid%nLeads
           call allocate_sigma_single(sigmaleads(il),
     &         ERealGrid%leadsTotalDim(il),0,
     &            ERealGrid%InfoSigma,
     &         (0.0D0,0.0D0),0.0D0,0.0D0)
         enddo

         call extractsigma(ERealGrid,sigmaleads,i,ispin,ik)


         GF_iter(1:NL,1:NL)=GF_iter(1:NL,1:NL)-
     &       sigmaleads(1)%sigma 
         GF_iter(N1-NR+1:N1,N1-NR+1:N1)=
     &   GF_iter(N1-NR+1:N1,N1-NR+1:N1)-
     &       sigmaleads(2)%sigma 

#ifdef ESSL
      if (.false.) then
#else
      if (gamma_negf) then
#endif
   	  CALL ZSYTRF('L',N1,GF_iter,N1,IPIV,WORK,N1**2,INFO)
	  CALL ZSYTRI('L',N1,GF_iter,N1,IPIV,WORK(1:2*N1),INFO)       
	 else
          CALL ZGETRF(N1,N1,GF_iter,N1,IPIV,INFO)
          CALL ZGETRI(N1,GF_iter,N1,IPIV,WORK,N1**2,INFO)
         endif
        	     
	 DO II=1,N1
	  DO JJ=1,N1
	   GF_iter_dag(II,JJ)=DCONJG(GF_iter(JJ,II))
	  ENDDO
	 ENDDO

        Gamma_L=zi*(sigmaleads(1)%sigma
     &      -DCONJG(Transpose(sigmaleads(1)%sigma)))
        Gamma_R=zi*(sigmaleads(2)%sigma
     &      -DCONJG(Transpose(sigmaleads(2)%sigma)))
          do il=1,ERealGrid%nLeads
            deallocate(sigmaleads(il)%sigma)
          enddo
          deallocate(sigmaleads)


	 CALL ZHEMM('L','U',NL,NR,(1.D0,0.D0),Gamma_L,NL,
     &    GF_iter_dag(1:NL,N1-NR+1:N1),NL,(0.D0,0.D0),Gamma1_aux,NL)
	 CALL ZHEMM('L','U',NR,NL,(1.D0,0.D0),Gamma_R,NR,
     &        GF_iter(N1-NR+1:N1,1:NL),NR,(0.D0,0.D0),Gamma2_aux,NR)

	 DO J=1,NL
	  Ic_aux(J)=0.D0
	  DO L=1,NR
	   Ic_aux(J)=Ic_aux(J)+Gamma1_aux(J,L)*Gamma2_aux(L,J)
	  ENDDO
	 ENDDO    

	 DO J=1,NL
	  Ic(ISPIN)=Ic(ISPIN)+wk*eh_const*
     &     DREAL(CONST*Ic_aux(J)*(fL-fR))
	 ENDDO
	 DO J=1,NL
	  tk(ISPIN)=tk(ISPIN)+DREAL(Ic_aux(J))
	 ENDDO
!	ENDDO

        deallocate(Gamma1_aux,Gamma2_aux,Gamma_L,Gamma_R)
       
	end subroutine current

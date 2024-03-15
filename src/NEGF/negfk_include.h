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
! THE CONTENTS
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
C ********* Variables *************************************
	
	IMPLICIT NONE

C ********** Parameters ***********************************
! Add variable inicoor, idyn, tmdskip, tmdsampling
! Meilin Bai, Dec 2012
!
	INTEGER :: N1,nsc(2),
     &   nspinL,maxnhL,nspinR,maxnhR,
     &   Nenerg,Nenerg1,Nenerg2,NPOLES,ITER,istep,ICASE,
     &   inicoor, idyn, tmdskip, tmdsampling 
	INTEGER, SAVE :: NL,NR
        INTEGER, SAVE :: N_IVPOINTS
	LOGICAL, SAVE :: PERIODIC
	
	INCLUDE "const2.h"

C ******** Labels ***********************************************************
	CHARACTER slabel*20,slabelL*20,slabelR*20,paste*25,curfile*25
	
C ********** Energy Related Variables  ********************

C Energy
	DOUBLE PRECISION :: Ei,Delta,EB
	DOUBLE PRECISION,SAVE :: Ef_Lead,R0,
     &                      Ef_LeadL,Ef_LeadR
	DOUBLE PRECISION :: Vini,Vfinal
C Temperature
	DOUBLE PRECISION :: T	
C Potencial bias
        DOUBLE PRECISION :: V,wk
        DOUBLE PRECISION, DIMENSION (2), SAVE :: Ic

	
	INTEGER :: IV
	INTEGER :: I,ISPIN,iuc
	INTEGER :: ik,nk,iksigma
	DOUBLE PRECISION, DIMENSION (3) :: kpoint
	INTEGER :: NSLICES
	LOGICAL :: CONVERGED,CALCTRANSM

C ******* Parallel Variables ******************************
        INTEGER :: Nnodes,MyNode


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
! THE MODULE
!                   SIGMA  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
      module sigma

!*****************************************************************
! Module that defines the real-axis self-energies
!
! Written by Alexandre Reily Rocha, Chaitanya Das Pemmaraju and
! Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: rochaa@tcd.ie, reilya@if.usp.br, runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version: October 2008
! *****************************************************************

      IMPLICIT NONE
    
      integer, dimension (:) , allocatable, save :: mydepth
      
! ********** Leads H0,H1,S0,S1 ****************************
      integer SigmaNl,SigmaNr
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: H0_L(:,:,:),H1_L(:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: S0_L(:,:),S1_L(:,:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: H0_R(:,:,:),H1_R(:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: S0_R(:,:),S1_R(:,:)


      integer SigmaNl2,SigmaNr2
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: H0_L2(:,:,:),H1_L2(:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: S0_L2(:,:),S1_L2(:,:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: H0_R2(:,:,:),H1_R2(:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: S0_R2(:,:),S1_R2(:,:)

! ********** Decimation Variables *************************
      INTEGER, ALLOCATABLE,DIMENSION (:,:), SAVE :: numberL,numberR

!---------Variables for Self Energy Interpolation----------------
      
      integer, save:: sigma_method
!------Fixed grid parameters
      logical,save::first_grid
      integer,save :: Tot_Gp,Gp_node
      double precision, save :: EngI1, EngF1, FixdE,Vb
      double complex, allocatable, save :: FixEG(:)
      integer, allocatable, save :: GpToNode(:),GpToLp(:)
!-----File I/O Variables----
      INTEGER,save:: Sig_RLen,Sig_RNum,fu_sl,fu_sr
      CHARACTER*10,save::nodenum
      CHARACTER*15,save::fn_sl,fn_sr
!     LOGICAL,SAVE::firstL,firstR
!     DATA firstL/.true./,firstR/.true./
      Character*8::SigSave
      LOGICAL  SigInter 
      
!----------------------------------------------------------------


      end module sigma
      

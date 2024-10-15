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
!                   NEGF_INTERFACE,
!                   SET_KPOINTS_NODE,
!                   OUTPUT_TOTALCHARGE,
!                   SET_BOUNDARY_ELEMENTS,
!                   RHO_BOUNDARYL2,
!                   RHO_BOUNDARYR2,
!                   BOUNDARY_ELEMENTS_MPI2,
!                   PRINTHSD  
! AND
! THE MODULE
!                   MNEGF_INTERFACE  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
MODULE mNegf_Interface 

  use mConstants
  use mTypes
  use negfmod
  use mMatrixUtil
  use mMPI_NEGF
  use set_rhobd
  use mCurrDistTotal, only : CurrentDistributionMatrix

  implicit none

  CONTAINS



  subroutine Negf_Interface(H, S, DM,Omega,xij, no_s, no_u,no_u_node, xa,na_u,na_s,NspinRealInputMatrix, maxnh, numh, listhptr,  listh,nkpts, kpoint, weight_k, indxuo, iaorb,MD_Step, inicoor, &
  IV_step,Vb ,SCF_step,Last_SCF_step,StructureChanged,&
   temp, nsc, slabel, NEnergR, NEnergIC, NEnergIL, NPoles,Delta, EnergLB, VInitial, VFinal, SpinCL, NSlices, TrCoeff, CalcIETS, ldos, idyn, tmdskip, tmdsampling)
  !


  use ScissorOperator, only : sco,scoSCF,SCOSetHamiltonianBlock,SCOApplyK
  use mIO_sxdhe, only : WriteRealMatPFormatted

  integer, intent(in) :: no_s,no_u,no_u_node,NspinRealInputMatrix,maxnh,nkpts,IV_step,SCF_step,MD_step, inicoor, idyn, tmdskip, tmdsampling
  integer, intent(in) :: numh(no_u_node),listhptr(no_u_node),listh(maxnh),indxuo(no_s),iaorb(no_s)
  
  logical, intent(in) :: Last_SCF_Step,StructureChanged

  real(kdp),intent(in) :: H(maxnh,NspinRealInputMatrix)
  real(kdp),intent(in) :: S(maxnh)
  real(kdp),intent(inout) :: DM(maxnh,NspinRealInputMatrix)
  real(kdp),intent(inout) :: Omega(maxnh,NspinRealInputMatrix)
  real(kdp),intent(in) :: xij(3,NspinRealInputMatrix)
  real(kdp),intent(in) :: kpoint(3,nkpts),weight_k(nkpts)
  real(kdp),intent(in) :: Vb
  integer, intent(in)  :: na_u,na_s
  real(kdp),intent(in) :: xa(3,na_s)


  real(kdp),allocatable :: DMImag(:,:),OmegaImag(:,:)

  real(kdp) :: weight_k_use(nkpts)
  integer  BNode, ie, ierror, io, iio, ispin,ispincompare
  integer MPIerror,is
#ifdef MPI
  integer istatus(MPI_STATUS_SIZE)
#endif
  integer, save :: maxnhg
  integer, dimension(:), allocatable, save :: numhg,listhptrg,listhg
  integer, dimension(:), allocatable, save :: listig

  integer  nsc(2),  NEnergR, NEnergIC, NEnergIL, Npoles,  SpinCL, NSlices,  node, Nodes,nkeff,nkP,ikP
  double precision  temp, Delta, EnergLB, VInitial, VFinal  
  double precision, allocatable :: Hk(:,:),Sk(:),xijKa(:,:),dnewk(:,:)
  character*20  slabel
  logical  TrCoeff,ldos, CalcIETS
  type(ioType) :: iout
  integer, allocatable :: numhk(:),listhptrk(:)

  integer iuo, j, juo, jo, jo2,ind, ind2,ik, i, iu1, iu2, nspinL, nspinR, nn, juoL, jL, indL, juoR, jR, indR, NspinComplexMatrix,NspinBlocks
  integer, save::  nuoL, nuoR, maxnhL, maxnhR,NeqR,NeqL,NeqOffL,NeqOffR
  integer, allocatable, save:: numdL(:), listdptrL(:), listdL(:), numdR(:), listdptrR(:), listdR(:), iequivL(:,:),iequivR(:,:),iequivOffL(:,:),iequivOffR(:,:)
  double precision kxij
  double precision, allocatable, save:: DMbulkL(:,:), DMbulkR(:,:)
  double complex sckxij
  double complex, parameter :: ii= (0.d0, 1.d0)
  character*20 slabelL, slabelR
  external negfk, memory, io_assign, io_close, dmbk
  type(matrixTypeGeneral), allocatable :: rhogeneralp(:)
  type(matrixTypeGeneral), allocatable :: ematgeneralp(:)
  type(matrixTypeGeneral), allocatable :: hgeneralp(:)
  type(matrixTypeGeneral) :: sgeneralp
  type(matrixTypeGeneral), allocatable :: hgeneralpK(:)
  type(matrixTypeGeneral), allocatable :: rhogeneralpK(:)
  type(matrixTypeGeneral), allocatable :: ematgeneralpK(:)
  type(matrixTypeGeneral) :: sgeneralpK
  type(matrixTypeGeneral) :: mat1
  type(matrixTypeGeneral), allocatable :: xijK(:)
  type(matrixSparseType) :: deltas
  integer maxnelerow,i1,HSblocksize,ikstart,ikend,ikstep,nprocskp
  integer maxnelerowk
  integer nmat
  integer*4:: sc_0,sc_1,sc_r,sc_m
  integer*4:: sc_0b
  integer*4:: sc_0k
  double precision qsol,buffer1
  real(kdp),allocatable :: xo(:,:)

!zrx,
  double precision, allocatable, save:: EDMbulkL(:,:), EDMbulkR(:,:)
!zrx,
  character*20 istring

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    CALL SYSTEM_CLOCK(sc_0b,sc_r,sc_m)
  endif

  em_Last_SCF_Step=Last_SCF_Step

  iout%isDebug=.false.
  NspinComplexMatrix=NspinRealInputMatrix
  NspinBlocks=NspinRealInputMatrix
  nmat=no_u
  if(NspinRealInputMatrix >= 4)then
    NspinBlocks=4
    if(.not.negfon)then
      NspinComplexMatrix=4
    else
      NspinComplexMatrix=1
      nmat=2 * no_u
    endif
  endif
  EM_NSPINBlocks=NspinBlocks

  allocate(rhogeneralp(NspinComplexMatrix))
  allocate(ematgeneralp(NspinComplexMatrix))
  allocate(hgeneralp(NspinComplexMatrix))
  if(NParallelK>1) then
    allocate(hgeneralpK(NspinRealInputMatrix))
    allocate(rhogeneralpK(NspinRealInputMatrix))
    allocate(ematgeneralpK(NspinRealInputMatrix))
    allocate(xijK(3))
  endif




  if(NParallelK>1)then
    call convertmatrixsiestatosmeagolparallelK2(H,S,maxnh,numh, listhptr,listh,no_u_node,xij,NspinRealInputMatrix, NspinComplexMatrix,no_u, indxuo,no_s,hgeneralpK,sgeneralpK,rhogeneralpK, ematgeneralpK,xijK,emforces,   maxnelerowk,nnodes_negf,mynode_negf,HSblocksize)
  endif

  call set_kpoints_node(nnodes_negfo,mynode_negfo,onekp,nkpts,NParallelK,nkP,ikstart,ikend,ikstep)


  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_ni_setup',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    CALL SYSTEM_CLOCK(sc_0k,sc_r,sc_m)
  endif


!  write(12347,*)"allocating dmimag"
  if(curr_dist)then
    allocate(DMImag(maxnh,NspinRealInputMatrix))
    allocate(OmegaImag(maxnh,NspinRealInputMatrix))
    DMImag = 0.d0
    OmegaImag = 0.d0

    allocate(xo(3,no_u))
    call Xa2Xo(xa,na_s,xo,no_u,no_s,iaorb)
!    do i=1,no_u
!      write(12347,*)"xorb=",i,xo(1,i),xo(2,i),xo(3,i)
!    enddo
  endif
 
  
!  call WriteRealMatPFormatted( mynode_negfo,nnodes_negfo,maxnh, no_u_node, no_u,NspinRealInputMatrix, numh, listhptr, listh, H,12347,"h0_ni")

  DM = 0.d0
  Omega = 0.d0
  ikP=0

  do ik = ikstart, ikend,ikstep

    ikP=ikP+1
    weight_k_use(:)=weight_k(:)
    if(onekp.eq.1)weight_k_use(1)=1D0
    kpointmod=kpoint(1:3,ik)
    ikpmod=ik
    ikpmodK=ikP
    wkmod=weight_k_use(ik)

    if(NParallelK>1)then
      allocate(Hk(sgeneralpK%matSparseP%matSparse%nnz, NspinRealInputMatrix),Sk(sgeneralpK%matSparseP%matSparse%nnz))
      allocate(xijKa(3,sgeneralpK%matSparseP%matSparse%nnz))

      do ispin=1,NspinRealInputMatrix
        Hk(:,ispin)=hgeneralpK(ispin)%matSparseP%matSparse%b(:)
      enddo
      Sk(:)=sgeneralpK%matSparseP%matSparse%b(:)
      do ispin=1,3
        xijKa(ispin,:)=xijK(ispin)%matSparseP%matSparse%b(:)
      enddo

      allocate(numhk(sgeneralpK%matSparseP%matSparse%iRows))
      allocate(listhptrk(sgeneralpK%matSparseP%matSparse%iRows))
      numhk=0
      do i=1,sgeneralpK%matSparseP%matSparse%iRows
        numhk(i)=sgeneralpK%matSparseP%matSparse%q(i+1)- sgeneralpK%matSparseP%matSparse%q(i)
        listhptrk(i)=sgeneralpK%matSparseP%matSparse%q(i)-1
      enddo

      call convertmatrixsiestatosmeagolparallel_2k(sgeneralpK%matSparseP,Hk,Sk, sgeneralpK%matSparseP%matSparse%nnz, numhk, listhptrk,  sgeneralpK%matSparseP%matSparse%j,  sgeneralpK%matSparseP%matSparse%irows,   xijKa,NspinRealInputMatrix,  NspinComplexMatrix,no_u,  nmat,   indxuo,no_s,kpoint(:,ik),hgeneralp,sgeneralp,rhogeneralp,  ematgeneralp,emforces, maxnelerow,nnodes_inverse,mynode_inverse,2) 
      deallocate(listhptrk,numhk)
      deallocate(Hk,Sk,xijKa)

    else
      call convertmatrixsiestatosmeagolgeneral(H,S,maxnh,numh, listhptr,listh,no_u_node,xij,NspinRealInputMatrix,  NspinComplexMatrix,no_u, nmat, indxuo,no_s,kpoint(:,ik),hgeneralp,sgeneralp,rhogeneralp,  ematgeneralp,emforces, maxnelerow,nnodes_inverse,mynode_inverse,2)

!!!!      if(.false.)then
!!!!        call PrintSparse2DenseReorderedNC(hgeneralp,no_u, NspinComplexMatrix, NspinBlocks, NspinComplexMatrix,1,"hn_0")
!!!!        call PrintSparse2DenseReorderedNC(sgeneralp,no_u, 1, NspinBlocks, NspinComplexMatrix,2,"sn_0")
!!!!      endif

    endif


    if(nprocs_hs.ne.1)then
      call convertmatrixsiestatosmeagolgeneral(H,S,maxnh,numh, listhptr,listh,no_u_node,xij,NspinRealInputMatrix,  NspinComplexMatrix,no_u, nmat, indxuo,no_s,kpoint(:,ik),hgeneralp,sgeneralp,rhogeneralp, ematgeneralp,emforces, maxnelerow,nnodes_inverse,mynode_inverse,3)
      do ispin=1,NspinComplexMatrix
        hgeneralp(ispin)%mattype=2
      enddo
      sgeneralp%mattype=2

      call MaxDifferenceMatrixCRS_CRSP(sgeneralp,sgeneralp,mynode_inverse,"(sgeneralp) ",iout)

      do ispin=1,NspinComplexMatrix
        call MaxDifferenceMatrixCRS_CRSP(hgeneralp(ispin),  hgeneralp(ispin),mynode_inverse,"(hgeneraludp) ",iout)
      enddo
    endif


    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_ni_convertmatrix',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!  ScissorOperator !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!will not work in parallel over k mode, check mynode and nnodes, and then broadcasts inside SCOApplyK
    if (sco.and.(Last_SCF_Step.or.scoSCF)) then
      SCOSetHamiltonianBlock=.true.
      call SCOApplyK(NspinRealInputMatrix,maxnh,H,S,xij,kpoint(:,ik),mynode_negfo,nnodes_negfo,no_u,no_s,numh,listhptr,listh,indxuo,hgeneralp,sgeneralp,NspinComplexMatrix)
    else
      SCOSetHamiltonianBlock=.false.
    endif
!!!!!!!!!!!!!!!!!!!  ScissorOperator !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!    do ispin=1,NspinComplexMatrix
!!!      write(istring,'(i10)')ispin
!!!      call PrintMatrixCRS(hgeneralp(1)%matSparse,"h_"//TRIM(ADJUSTL(istring)),iout)
!!!    enddo

!!!!    if(.true.)then
!!!!      call PrintSparse2DenseReorderedNC(hgeneralp,no_u, NspinComplexMatrix, NspinBlocks, NspinComplexMatrix,1,"hn_0")
!!!!      call PrintSparse2DenseReorderedNC(sgeneralp,no_u, 1, NspinBlocks, NspinComplexMatrix,2,"sn_0")
!!!!    endif


    call negfk(no_u,NspinBlocks,NspinComplexMatrix,nsc,NEnergR,NEnergIL,  NEnergIC, Npoles,MD_step, inicoor, IV_step,SCF_step,ikP,SpinCL,VInitial, VFinal,temp,Delta,Vb,EnergLB,Last_SCF_step,NSlices,TrCoeff, CalcIETS, nkP,kpoint(1:3,ik),slabel,weight_k_use(ik), hgeneralp,sgeneralp, rhogeneralp,ematgeneralp,ldos, idyn, tmdskip, tmdsampling)

!!!!    if(.true.)then
!!!!      call PrintSparse2DenseReorderedNC(rhogeneralp,no_u, NspinComplexMatrix, NspinBlocks, NspinComplexMatrix,1,"rho_0")
!!!!    endif

      if(RhoSetZeroIfHZero)then
      if(NspinBlocks==4.and.NspinBlocks==NspinComplexMatrix)then
        do ispin=1,NspinComplexMatrix
          do i=1,rhogeneralp(ispin)%Matsparse%iRows
            do ind=rhogeneralp(ispin)%Matsparse%q(i),rhogeneralp(ispin)%Matsparse%q(i+1)-1

              ispincompare=ispin
              if((ispin)==4)ispincompare=3
              if(hgeneralp(ispincompare)%Matsparse%b(ind)==0.0D0)then
                if(rhogeneralp(ispin)%Matsparse%b(ind).ne.0.0D0)then
                  rhogeneralp(ispin)%Matsparse%b(ind)=0.0D0
!                   write(12347,*)"drho=",i,rhogeneralp(ispin)%Matsparse%j(ind),DREAL(rhogeneralp(ispin)%Matsparse%b(ind)),DIMAG(rhogeneralp(ispin)%Matsparse%b(ind))
                endif
              endif
            enddo
          enddo
        enddo
      endif
    endif

!!!!    do ispin=1,NspinComplexMatrix
!!!!      write(istring,'(i10)')ispin
!!!!      call PrintMatrixCRS(rhogeneralp(ispin)%matSparse,"r_"//TRIM(ADJUSTL(istring)),iout)
!!!!    enddo


    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_ni_negfk',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif

    if(mynode_inverse.eq.0)then
      do ispin=1,NspinComplexMatrix
        hgeneralp(ispin)%mattype=2
        call DestroyMatrixGeneral(hgeneralp(ispin),"negfk",iout)
      enddo
      sgeneralp%mattype=2
      call DestroyMatrixGeneral(sgeneralp,"negfk",iout)
    endif

    if(nprocs_hs.ne.1)then
      do ispin=1,NspinComplexMatrix
        hgeneralp(ispin)%mattype=3
        call DestroyMatrixGeneral(hgeneralp(ispin),"negfk",iout)
        hgeneralp(ispin)%mattype=2
      enddo
      sgeneralp%mattype=3
      call DestroyMatrixGeneral(sgeneralp,"negfk",iout)
      sgeneralp%mattype=2
    endif

    if(NParallelK>1)then

      allocate(xijKa(3,sgeneralpK%matSparseP%matSparse%nnz))
      do ispin=1,3
        xijKa(ispin,:)=xijK(ispin)%matSparseP%matSparse%b(:)
      enddo

      allocate(dnewK(sgeneralpK%matSparseP%matSparse%nnz, NspinRealInputMatrix))
      dnewk=0.0D0

      call rhoK_to_rho2(rhogeneralpK(1)%matSparseP%matSparse,maxnelerow,NspinRealInputMatrix,NspinComplexMatrix,rhogeneralp,dnewk,maxnh,no_s,no_u_node,no_u,nmat, listh, numh, listhptr,indxuo, kpoint(:,ik),xijKa,weight_k_use(ik),timereversal)
      do ispin=1,NspinRealInputMatrix
        rhogeneralpK(ispin)%matSparseP%matSparse%b(:)=rhogeneralpK(ispin)%matSparseP%matSparse%b(:)+dnewk(:,ispin)
      enddo

      if(emforces) then
        dnewk=0.0D0
        call rhoK_to_rho2(rhogeneralpK(1)%matSparseP%matSparse,maxnelerow,NspinRealInputMatrix,NspinComplexMatrix,ematgeneralp,dnewk,maxnh,no_s,no_u_node,no_u,nmat, listh, numh, listhptr,indxuo, kpoint(:,ik),xijKa,weight_k_use(ik),timereversal)
        do ispin=1,NspinRealInputMatrix
          ematgeneralpK(ispin)%matSparseP%matSparse%b(:)=ematgeneralpK(ispin)%matSparseP%matSparse%b(:)+dnewk(:,ispin)
        enddo
      endif

      deallocate(dnewk,xijKa)

    else

!!!!      if(.true.)then
!!!!        call PrintSparse2DenseReorderedNC(rhogeneralp,no_u, NspinComplexMatrix, NspinBlocks, NspinComplexMatrix,1,"rho_1")
!!!!      endif

      call rhonegf_to_rhosiesta(maxnelerow,NspinRealInputMatrix,  NspinComplexMatrix,rhogeneralp,DM,   maxnh,no_s,no_u_node,no_u,nmat, listh, numh, listhptr,indxuo, kpoint(:,ik),xij,weight_k_use(ik),timereversal)

      if(emforces)   call rhonegf_to_rhosiesta(maxnelerow,NspinRealInputMatrix,  NspinComplexMatrix,ematgeneralp,    Omega, maxnh,no_s,no_u_node,no_u,nmat, listh, numh, listhptr,indxuo,  kpoint(:,ik),xij,weight_k_use(ik),timereversal)

      if(curr_dist)then
        call rhonegf_to_rhosiestaImaginary(maxnelerow,NspinRealInputMatrix,  NspinComplexMatrix,rhogeneralp,DMImag,       maxnh,no_s,no_u_node,no_u,nmat, listh, numh, listhptr,indxuo, kpoint(:,ik),xij,weight_k_use(ik),timereversal)
        call rhonegf_to_rhosiestaImaginary(maxnelerow,NspinRealInputMatrix,  NspinComplexMatrix,ematgeneralp,OmegaImag,   maxnh,no_s,no_u_node,no_u,nmat, listh, numh, listhptr,indxuo, kpoint(:,ik),xij,weight_k_use(ik),timereversal)
      endif





    endif

    if(.false.)call output_totalcharge(S,DM,NspinRealInputMatrix,maxnh,negfo_comm,mynode_negfo,"q_EM1=")

    if(mynode_inverse.eq.0)then
      do ispin=1,NspinComplexMatrix
        call DestroyMatrixGeneral(rhogeneralp(ispin),"emtkon",iout)
        call DestroyMatrixGeneral(ematgeneralp(ispin),  "emtkon",iout)
      enddo
    endif

    if(nprocs_hs.ne.1)then
      do ispin=1,NspinComplexMatrix
        rhogeneralp(ispin)%mattype=3
        ematgeneralp(ispin)%mattype=3
        call DestroyMatrixGeneral(rhogeneralp(ispin),"emtkon",iout)
        call DestroyMatrixGeneral(ematgeneralp(ispin), "emtkon",iout)
         rhogeneralp(ispin)%mattype=2
        ematgeneralp(ispin)%mattype=2
      enddo
    endif

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'t_ni_convertback',(sc_1-sc_0)*1.0d0/sc_r
      write(12347,'(A,i6,f12.6)')'t_ni_ik',ik,(sc_1-sc_0k)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      CALL SYSTEM_CLOCK(sc_0k,sc_r,sc_m)
    endif



  enddo

  if(curr_dist)then
    call CurrentDistributionMatrix(DMImag,OmegaImag,H,S,maxnh,numh,listhptr,listh,  no_u_node,xij,xo,NspinRealInputMatrix,NspinComplexMatrix,no_u,  indxuo,no_s,Vb,.false.)
    deallocate(DMImag,OmegaImag)
    deallocate(xo)
  endif

  if(em_endcode1) call stopnegf

  if(NParallelK>1)then

    call merge_rhok(NspinRealInputMatrix,rhogeneralpK)
    if(emforces) call merge_rhok(NspinRealInputMatrix,ematgeneralpK)
  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_merge_rhok',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
endif

    call rhonegfk_to_rhosiesta(maxnelerowk,NspinRealInputMatrix,  rhogeneralpK,DM,maxnh,no_s,no_u_node,no_u, listh, numh, listhptr,indxuo)
  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_rhonegfk_to_rhosiesta',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
endif
    if(emforces)  call rhonegfk_to_rhosiesta(maxnelerowk,NspinRealInputMatrix, ematgeneralpK,Omega, maxnh,no_s,no_u_node,no_u, listh, numh, listhptr,indxuo)

  endif

  CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  call set_boundary_elements(S,DM,Omega,NspinRealInputMatrix,no_u,no_u_node,maxnh,numh,listhptr,listh,SCF_step,MD_step,inicoor, IV_step,StructureChanged,Last_SCF_step,nsc,Vb)

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_ni_setboundary',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif


!  call printhsd(DM(:,1),maxnh,numh,listhptr,listh,no_u_node,no_u,  indxuo,no_s,"dmi_1")


  deallocate(rhogeneralp,ematgeneralp,hgeneralp)


  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)')'t_ni_all',(sc_1-sc_0b)*1.0d0/sc_r
  endif


  end subroutine Negf_Interface



  subroutine set_kpoints_node(nnodes_negfo,mynode_negfo,onekp,nkpts,NParallelK,nkP,ikstart,ikend,ikstep)


    integer, intent(in) :: onekp,nkpts,nnodes_negfo,mynode_negfo,NParallelK
    integer, intent(out) :: nkP,ikstart,ikend,ikstep

    integer nkeff,nprocskp,ikP,ik

    if(onekp.eq.1)then
      nkeff=1
    else
      nkeff=nkpts
    endif

    nprocskp=nnodes_negfo/NParallelK
!    write(12347,*)"nprocskp=",nprocskp,nnodes_negfo,NParallelK
    ikstart=mynode_negfo/nprocskp+1
    ikend=nkeff
    ikstep=NParallelK

    ikP=0
    do ik = ikstart, ikend,ikstep
      ikP=ikP+1
    enddo
    nkP=ikP
!    write(12347,*)"ikstart=",ikstart,ikend,ikstep,nprocskp,nkP

  end subroutine set_kpoints_node


  subroutine output_totalcharge(S,DM,NspinRealInputMatrix,maxnh,negfo_comm,mynode_negfo,nam)

    integer, intent(in) :: maxnh,NspinRealInputMatrix,negfo_comm,mynode_negfo
    real(kdp),intent(in) :: S(maxnh)
    real(kdp),intent(inout) :: DM(maxnh,NspinRealInputMatrix)
    character(len=*), intent(in) :: nam


    integer MPIerror,i,ispin
    real(kdp) buffer1,qsol
    

    qsol = 0.0D0
    do ispin = 1,min(NspinRealInputMatrix,2)
      do i = 1,maxnh
        qsol = qsol + DM(i,ispin) * s(i)
      enddo
    enddo
#ifdef MPI
    call MPI_AllReduce(qsol,buffer1,1,DAT_double, MPI_sum,negfo_comm,MPIerror)
    qsol=buffer1
#endif

    if (mynode_negfo.eq.0) then
      write(*,*)nam,qsol
    endif
  end subroutine output_totalcharge






  subroutine set_boundary_elements(S,DM,Omega,NspinRealInputMatrix,&
  no_u,no_u_node,maxnh,numh,listhptr,listh,SCF_step,MD_step,inicoor,&
  IV_step,StructureChanged,Last_SCF_step,nsc,Vb)

! Add argument inicoor intead of 0
! Meilin Bai, Dec 2012

  integer, intent(in) :: no_u,maxnh,no_u_node,SCF_step,MD_step,inicoor,IV_step,NspinRealInputMatrix,nsc(2)
  integer, intent(in) :: numh(no_u_node),listhptr(no_u_node),listh(maxnh)
  logical, intent(in) :: Last_SCF_Step,StructureChanged
  real(kdp),intent(in) :: S(maxnh)
  real(kdp),intent(inout) :: DM(maxnh,NspinRealInputMatrix)
  real(kdp),intent(inout) :: Omega(maxnh,NspinRealInputMatrix)
  real(kdp),intent(in) :: Vb

  integer, save :: maxnhg
  integer, dimension(:), allocatable, save :: numhg,listhptrg,listhg,listig
  integer, save::  nuoL, nuoR, maxnhL, maxnhR,NeqR,NeqL,NeqOffL,NeqOffR
  double precision, allocatable, save:: DMbulkL(:,:), DMbulkR(:,:)
  double precision, allocatable, save:: EDMbulkL(:,:), EDMbulkR(:,:)
  integer, allocatable, save:: numdL(:), listdptrL(:), listdL(:), numdR(:), listdptrR(:), listdR(:), iequivL(:,:),iequivR(:,:),iequivOffL(:,:),iequivOffR(:,:)

  integer io,jo,iio,BNode,MPIerror,iu1,iu2,nn
  integer nspinL, nspinR
  character*20 slabelL, slabelR
  logical debug_q_em

  debug_q_em=.false.

  if(.true..or..not.allocated(numhg))then
    allocate(numhg(no_u))
    allocate(listhptrg(no_u))

    do io = 1,no_u
      call WhichNodeOrb(io,nnodes_negfo,BNode)
      if (mynode_negfo.eq.BNode) then
        call GlobalToLocalOrb(io,mynode_negfo,nnodes_negfo,iio)
        numhg(io) = numh(iio)
      endif
#ifdef MPI
      call MPI_Bcast(numhg(io),1,MPI_integer,BNode, negfo_comm,MPIerror)
#endif
    enddo

    listhptrg(1) = 0
    do io = 2,no_u
      listhptrg(io) = listhptrg(io-1) + numhg(io-1)
    enddo

    maxnhg = listhptrg(no_u)+numhg(no_u)



    allocate(listhg(maxnhg))
    allocate(listig(maxnhg))

    do io = 1,no_u
      listig(listhptrg(io)+1:listhptrg(io)+numhg(io)) = io
      call WhichNodeOrb(io,nnodes_negfo,BNode)
      if (mynode_negfo.eq.BNode) then
        call GlobalToLocalOrb(io,mynode_negfo,nnodes_negfo,iio)
        do jo = 1,numhg(io)
          listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
        enddo
      endif
#ifdef MPI
      call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer, BNode,negfo_comm,MPIerror)
#endif
    enddo
  endif



  if (mynode_negfo .EQ. 0) then
!from here everything is executed on mynode_negfo 0
! Add bulk density matrix: this is done on mynode_negfo 0, and then sent to the other nnodes_negfo

    if (SCF_step.eq.1 .and. MD_step.eq.inicoor .and. IV_step.eq.0) then
      call io_assign(iu1)
      open(iu1,file='bulklft.DAT',status='old')
      read(iu1,*) slabelL, nuoL, nspinL, maxnhL

      call io_assign(iu2)
      open(iu2,file='bulkrgt.DAT',status='old')
      read(iu2,*) slabelR, nuoR, nspinR, maxnhR

      allocate(DMbulkL(maxnhL,NspinRealInputMatrix))
      allocate(DMbulkR(maxnhR,NspinRealInputMatrix))

      allocate(numdL(nuoL),listdptrL(nuoL),listdL(maxnhL))
      allocate(numdR(nuoR),listdptrR(nuoR),listdR(maxnhR))

      DMbulkL = 0.d0
      call dmbk(nuoL,maxnhL,nspinL,slabelL,numdL,listdL,listdptrL,DMbulkL,1)
      DMbulkR = 0.d0
      call dmbk(nuoR,maxnhR,nspinR,slabelR,  numdR,listdR,listdptrR,DMbulkR,1)
      
      if(ThetaLeadsL.ne.0.0D0.or.PhiLeadsL.ne.0.0D0)then
        call rotateDM(DMbulkL,maxnhL,NspinRealInputMatrix, ThetaLeadsL,PhiLeadsL)
      endif
  
      if(ThetaLeadsR.ne.0.0D0.or.PhiLeadsR.ne.0.0D0)then
        call rotateDM(DMbulkR,maxnhR,NspinRealInputMatrix,  ThetaLeadsR,PhiLeadsR)
      endif

      if(emforces)then
        allocate(EDMbulkL(maxnhL,NspinRealInputMatrix))
        allocate(EDMbulkR(maxnhR,NspinRealInputMatrix))
        EDMbulkL = 0.d0
        call dmbk(nuoL,maxnhL,nspinL,slabelL, numdL,listdL,listdptrL,EDMbulkL,2)
        EDMbulkR = 0.d0
        call dmbk(nuoR,maxnhR,nspinR,slabelR, numdR,listdR,listdptrR,EDMbulkR,2)
      endif


      call io_close(iu1)
      call io_close(iu2)

      allocate(iequivL(maxnhL,2),iequivR(maxnhR,2))
      if (Set_RhoBoundaryOverlap_Leads) then
       allocate(iequivOffL(maxnhL,2),iequivOffR(maxnhR,2))
      else
       allocate(iequivOffL(1,1),iequivOffR(1,1))
      endif

    endif

    if ((SCF_step.eq.1).and.(StructureChanged)) then

      NeqL=0
      NeqOffL=0
      NeqR=0
      NeqOffR=0

      nn = nsc(1)*nsc(2)

      call rho_boundaryL2(nn,nuoL,maxnhL,listdptrL,listdL,numdL,    no_u, maxnhg,numhg,listhptrg,listhg,neqL,iequivL,nuoR,     Set_RhoBoundaryOverlap_Leads,NeqOffL,iequivOffL)

      call rho_boundaryR2(nn,nuoR,maxnhR,listdptrR,listdR,numdR,   no_u,   maxnhg,numhg,listhptrg,listhg,neqR,iequivR,nuoL,       Set_RhoBoundaryOverlap_Leads,NeqOffR,iequivOffR)

    endif

  endif


!up to here everything is executed on mynode_negfo 0, only mynode_negfo 0 has the full DM
#ifdef MPI
  if ((SCF_step.eq.1).and.(StructureChanged)) then

    call MPI_Bcast(NeqL,1,MPI_integer,0, negfo_comm,MPIerror)
    call MPI_Bcast(NeqR,1,MPI_integer,0, negfo_comm,MPIerror)
    if(Set_RhoBoundaryOverlap_Leads)then
      call MPI_Bcast(NeqOffL,1,MPI_integer,0, negfo_comm,MPIerror)
      call MPI_Bcast(NeqOffR,1,MPI_integer,0, negfo_comm,MPIerror)
    endif
    call MPI_Bcast(maxnhl,1,MPI_integer,0, negfo_comm,MPIerror)
    call MPI_Bcast(maxnhr,1,MPI_integer,0, negfo_comm,MPIerror)

    if(mynode_negfo.ne.0)then
      if(.not.allocated(iequivL))allocate(iequivL(maxnhL,2))
      if(.not.allocated(iequivR))allocate(iequivR(maxnhR,2))
      if(.not.allocated(DMbulkL))allocate(DMbulkL(maxnhl, NspinRealInputMatrix))
      if(.not.allocated(DMbulkR))allocate(DMbulkR(maxnhr, NspinRealInputMatrix))
      if(emforces)then
        if(.not.allocated(EDMbulkL))allocate(EDMbulkL(maxnhl,NspinRealInputMatrix))
        if(.not.allocated(EDMbulkR))allocate(EDMbulkR(maxnhr, NspinRealInputMatrix))
      endif
      if(Set_RhoBoundaryOverlap_Leads)then
        if(.not.allocated(iequivoffL))allocate(iequivoffL(maxnhL,2))
        if(.not.allocated(iequivoffR))allocate(iequivoffR(maxnhR,2))
      endif
    endif

    call MPI_Bcast(iequivL(1,1),maxnhL*2,MPI_integer,0, negfo_comm,MPIerror)
    call MPI_Bcast(iequivR(1,1),maxnhR*2,MPI_integer,0, negfo_comm,MPIerror)

    if(Set_RhoBoundaryOverlap_Leads)then
      call MPI_Bcast(iequivoffL(1,1),maxnhL*2,MPI_integer,0, negfo_comm,MPIerror)
      call MPI_Bcast(iequivoffR(1,1),maxnhR*2,MPI_integer,0,  negfo_comm,MPIerror)
    endif
  endif
#endif  
!zrx,because all variables have been declared using 'save'
!DO NOT put this block into IOnode block
  IF((IV_step.eq.0).and.(MD_step.eq.inicoor).and.(SCF_step.eq.1)) THEN
#ifdef MPI
#ifdef NODAT
     call MPI_Bcast(DMbulkL(1,1),maxnhL*NspinRealInputMatrix,  MPI_double_precision,0,negfo_comm,MPIerror)
     call MPI_Bcast(DMbulkR(1,1),maxnhR*NspinRealInputMatrix,  MPI_double_precision, 0,negfo_comm,MPIerror)
     if(emforces)then
       call MPI_Bcast(EDMbulkL(1,1),maxnhL*NspinRealInputMatrix,  MPI_double_precision,0,negfo_comm,MPIerror)
       call MPI_Bcast(EDMbulkR(1,1),maxnhR*NspinRealInputMatrix,  MPI_double_precision,0,negfo_comm,MPIerror)
     endif
#else
     call MPI_Bcast(DMbulkL(1,1),maxnhL*NspinRealInputMatrix, DAT_double,0,negfo_comm,MPIerror)
     call MPI_Bcast(DMbulkR(1,1),maxnhR*NspinRealInputMatrix, DAT_double, 0,negfo_comm,MPIerror)
     if(emforces)then
       call MPI_Bcast(EDMbulkL(1,1),maxnhL*NspinRealInputMatrix, DAT_double, 0,negfo_comm,MPIerror)
       call MPI_Bcast(EDMbulkR(1,1),maxnhR*NspinRealInputMatrix, DAT_double, 0,negfo_comm,MPIerror)
     endif
#endif
#endif
  ENDIF

  if(debug_q_em)call output_totalcharge(S,DM,NspinRealInputMatrix,maxnh,negfo_comm,mynode_negfo,"q_EM alone as output=")

  if(.not.Last_SCF_step)then
    if(Set_RhoBoundary_Leads)then
      call boundary_elements_mpi2(neql,maxnhl,iequivL,maxnhg, listhg,listig,listhptrg,no_u,DM,DMbulkL,NspinRealInputMatrix,   maxnh,numh,no_u_node,listhptr,listh)
      call boundary_elements_mpi2(neqr,maxnhr,iequivR,maxnhg, listhg,listig,listhptrg,no_u,DM,DMbulkR, NspinRealInputMatrix,  maxnh, numh,no_u_node,listhptr,listh)
      if(debug_q_em)call output_totalcharge(S,DM,NspinRealInputMatrix,maxnh,negfo_comm,mynode_negfo,"q_EM alone after setting leads layers=")
    endif



    if(Set_RhoBoundaryOverlap_Leads)then
      call boundary_elements_mpi2(neqoffl,maxnhl,iequivoffL,maxnhg, listhg,listig,listhptrg,no_u,DM,DMbulkL, NspinRealInputMatrix, maxnh, numh,no_u_node,listhptr,listh)
      call boundary_elements_mpi2(neqoffr,maxnhr,iequivoffR,maxnhg,  listhg,listig,listhptrg,no_u,DM,DMbulkR, NspinRealInputMatrix,  maxnh, numh,no_u_node,listhptr,listh)
      if(debug_q_em)call output_totalcharge(S,DM,NspinRealInputMatrix,maxnh,negfo_comm,mynode_negfo,"q_EM including overlap contribution=")
    endif
  endif

  
!       if(Last_SCF_step .and. SetEBD) then
   if(emforces) then
     if(Set_RhoBoundary_Leads)then
       call boundary_elements_mpi2(neql,maxnhl,iequivL,maxnhg, listhg,listig,listhptrg,no_u,Omega,EDMbulkL,  NspinRealInputMatrix,  maxnh, numh,no_u_node,listhptr,listh,DMbulkL,Vb,'L')
       call boundary_elements_mpi2(neqr,maxnhr,iequivR,maxnhg,  listhg,listig,listhptrg,no_u,Omega,EDMbulkR,   NspinRealInputMatrix,   maxnh, numh,no_u_node,listhptr,listh,DMbulkR,Vb,'R')
     endif
     if(Set_RhoBoundaryOverlap_Leads)then
      call boundary_elements_mpi2(neqoffl,maxnhl,iequivoffL,maxnhg, listhg,listig,listhptrg,no_u,Omega,EDMbulkL,  NspinRealInputMatrix,    maxnh,numh,no_u_node,listhptr,listh,DMbulkL,Vb,'L')
      call boundary_elements_mpi2(neqoffr,maxnhr,iequivoffR,maxnhg,   listhg,listig,listhptrg,no_u,Omega,EDMbulkR,    NspinRealInputMatrix,    maxnh,numh,no_u_node,listhptr,listh,DMbulkR,Vb,'R')
     endif
   endif


  if(.true..or..false.)then ! these need to be recalculated if structure of matrices changes, e.g. due to structural relaxations
    deallocate(numhg)
    deallocate(listhptrg)
    deallocate(listhg)
    deallocate(listig)
  endif


  end subroutine set_boundary_elements

!
subroutine rho_boundaryL2(nn,nuoL,maxnhL,listdptrL,listdL,numdL,nuotot,maxnhg,numhg,listhptrg,listhg,neqL,iequivL,nuoR,Set_RhoBoundaryOverlap_Leads,NeqOffL,iequivOffL)

    implicit none

    integer nuoL,nuotot,maxnhg,maxnhL,numhg(nuotot),listhptrg(nuotot),listhg(maxnhg),listdptrL(nuoL),listdL(maxnhL),neqL,iequivL(maxnhL,2),nuoR,NeqOffL,iequivOffL(maxnhL,2),numdL(nuoL)
    logical Set_RhoBoundaryOverlap_Leads
    integer iuo,i,j,ind,juo,nn,indL,juoL,jL

    do iuo = 1,nuoL
     do j = 1,numhg(iuo)
      ind = listhptrg(iuo) + j
      juo = listhg(ind)
      do i = 1,nn
       if(juo.gt.nuotot*(i-1).and.juo.le.nuotot*(i-1)+nuoL) then
        indL = listdptrL(iuo) + 1
        juoL = listdL(indL)
        jL = 1
        do while (((juoL-nuoL*(i-1)).ne.(juo-nuotot*(i-1))) .and.(jL .le. numdL(iuo)))
         indL = listdptrL(iuo) + jL
         juoL = listdL(indL)
         jL = jL + 1
        enddo
        if ((juoL-nuoL*(i-1)).eq.(juo-nuotot*(i-1))) then
         NeqL=NeqL+1
         iequivL(NeqL,1)=ind
         iequivL(NeqL,2)=indL
!    write(200,'(7i7)') ind,indL,juoL,juo,i,juoL-nuoL*(i-1),
!          juo-nuotot*(i-1)
        endif
       endif
       if(juo .gt. nuotot*i-nuoR .and. juo .le. nuotot*i.and. Set_RhoBoundaryOverlap_Leads) then
        indL = listdptrL(iuo) + 1
        juoL = listdL(indL)
        jL = 1
        do while ((juoL-nuoL*(2*nn+i-1)).ne.(juo-(nuotot*i-nuoR)) .and. (jL .le. numdL(iuo)))
         indL = listdptrL(iuo) + jL
         juoL = listdL(indL)
         jL = jL + 1
        enddo
        if ((juoL-nuoL*(2*nn+i-1)).eq.(juo-(nuotot*i-nuoR))) then
         NeqOffL=NeqOffL+1
         iequivOffL(NeqOffL,1)=ind
         iequivOffL(NeqOffL,2)=indL
! write(201,'(7i7)') ind,indL,juoL,juo,i,juoL-nuoL*(2*nn+i-1),
!          juo-(nuotot*i-nuoR)
         endif
       endif
      enddo
     enddo
    enddo
end subroutine rho_boundaryL2


subroutine rho_boundaryR2(nn,nuoR,maxnhR,listdptrR,listdR,numdR,nuotot,maxnhg,numhg,listhptrg,listhg,neqR,iequivR,nuoL,Set_RhoBoundaryOverlap_Leads,NeqOffR,iequivOffR)

    implicit none

    integer nuoR,nuotot,maxnhg,maxnhR,numhg(nuotot),listhptrg(nuotot),listhg(maxnhg),listdptrR(nuoR),listdR(maxnhR),neqR,iequivR(maxnhR,2),nuoL,NeqOffR,iequivOffR(maxnhR,2),numdR(nuoR)
    logical Set_RhoBoundaryOverlap_Leads
    integer iuo,i,j,ind,juo,nn,indR,juoR,jR

     do iuo = nuotot-nuoR+1,nuotot
           do j = 1, numhg(iuo)
            ind = listhptrg(iuo) + j
            juo = listhg(ind)
            do i = 1,nn
             if(juo.gt.nuotot*i-nuoR.and.juo.le.nuotot*i) then
              indR = listdptrR(iuo-(nuotot-nuoR)) + 1
              juoR = listdR(indR)
              jR = 1
              do while ((juoR-nuoR*(i-1)).ne.(juo-(nuotot*i-nuoR)) .and.(jR .le. numdR(iuo-(nuotot-nuoR))))
               indR = listdptrR(iuo-(nuotot-nuoR)) + jR
               juoR = listdR(indR)
               jR = jR + 1
              enddo
              if ((juoR-nuoR*(i-1)).eq.(juo-(nuotot*i-nuoR))) then
               NeqR=NeqR+1
               iequivR(NeqR,1)=ind
               iequivR(NeqR,2)=indR
              endif
             endif
             if(juo.gt.nuotot*(i-1).and.juo.le.nuotot*(i-1)+nuoL .and.Set_RhoBoundaryOverlap_Leads) then
              indR = listdptrR(iuo-(nuotot-nuoR)) + 1
              juoR = listdR(indR)
              jR = 1
              do while ((juoR-nuoL*(nn+i-1)).ne.(juo-nuotot*(i-1)) .and.(jR .le. numdR(iuo-(nuotot-nuoR))))
               indR = listdptrR(iuo-(nuotot-nuoR)) + jR
               juoR = listdR(indR)
               jR = jR + 1
              enddo
              if ((juoR-nuoL*(nn+i-1)).eq.(juo-nuotot*(i-1))) then
               NeqOffR=NeqOffR+1
               iequivOffR(NeqOffR,1)=ind
               iequivOffR(NeqOffR,2)=indR
              endif
             endif
            enddo
           enddo
          enddo



end subroutine rho_boundaryR2

  subroutine boundary_elements_mpi2(neql,maxnhl,iequivL,maxnhg, listhg,listig,listhptrg,nuotot,dnew,DMbulkL, NspinRealInputMatrix, maxnh,numh,nuo,listhptr,listh,DMlead,iBias,flag)
 
  integer neql,maxnhl,nuotot,iequivL(maxnhL,2),maxnhg,nuo, listhg(maxnhg),listig(maxnhg),listhptrg(nuotot),NspinRealInputMatrix,maxnh, numh(nuo),listhptr(nuo),listh(maxnh)
  double precision DMbulkL(maxnhL,NspinRealInputMatrix),  dnew(maxnh,NspinRealInputMatrix),dnewele
  integer ind,mpierror,iuo,jo,bnode,i,j,ind2,jo2,iio,  ispin
  integer :: iopt
  double precision, optional :: DMlead(maxnhL,NspinRealInputMatrix), iBias
  character(len=*), optional :: flag

#ifdef MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

  
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
         Dnew(ind2,:)=DMbulkL(iequivL(ind,2),:)+  (iBias/2.d0)*DMlead(iequivL(ind,2),:)
      ELSE if(iopt.eq.3) then
         Dnew(ind2,:)=DMbulkL(iequivL(ind,2),:)- (iBias/2.d0)*DMlead(iequivL(ind,2),:)
      ENDIF
    endif
  enddo

  end subroutine boundary_elements_mpi2



subroutine printhsd(S,maxnh,numh,listhptr,listh,n1loc,N1,indxuo,no,nam)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF
#ifdef MPI
  use parallel
#endif

  implicit none


  integer maxnh,NspinComplexMatrix,NspinRealInputMatrix,n1,n1loc
  integer nspinMin
  double precision S(maxnh)

  integer io,j,no,ind,jo,iuo,juo,indxuo(no),i,numh(n1loc),listhptr(n1loc),listh(maxnh),iio,ii,ispin
  integer ind2,is,maxnelerow,ind4,mynode,nnodes,mpierror,bnode
  CHARACTER(LEN=*) :: nam

#ifdef MPI
  call MPI_Comm_Rank(negfo_comm,myNode,MPIerror)
  call MPI_Comm_Size(negfo_comm,NNodes,MPIerror)
#else
  MyNode = 0
  NNodes = 1
#endif


  do io = 1,n1

#ifdef MPI
    CALL MPI_BARRIER(negfo_comm, MPIerror)
#endif

    call WhichNodeOrb(io,NNodes,BNode)
    if (MyNode.eq.BNode) then
      call GlobalToLocalOrb(io,MyNode,NNodes,iio)
      do j = 1,numh(iio)
        ind = listhptr(iio) + j
!        juo = indxuo(listh(ind))
        write(12347,*)nam,io,listh(ind),s(ind)
      enddo
    endif

  enddo

end subroutine printhsd


subroutine Xa2Xo(xa,na_s,xo,no_u,no_s,iaorb)

integer, intent(in)   :: na_s, no_u,no_s
real(kdp), intent(in) :: xa(3,na_s)
integer, intent(in)   :: iaorb(no_s)

real(kdp), intent(out) :: xo(3,no_u)

integer i

do i=1,no_u
!  write(12347,*)"xo=",i,iaorb(i),xa(1,iaorb(i)),xa(2,iaorb(i)),xa(3,iaorb(i))
  xo(:,i)=xa(:,iaorb(i))
enddo

end subroutine Xa2Xo


END MODULE mNegf_Interface

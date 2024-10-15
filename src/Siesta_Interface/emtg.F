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
      subroutine emtg(NspinRealInputMatrix, nuotot, nuo, no,
     .    xa, na_u,na_s,
     .    maxnh, numh,
     .    listhptr, listh, indxuo,iaorb, H, S, temp, nsc, xij,slabel,
     .    NEnergR, NEnergIC, NEnergIL, NPoles,
     .    Delta, EnergLB, Iter, last, istep, inicoor, ivv, VBias,
     .    VInitial, VFinal, SpinCL, NSlices, TrCoeff, CalcIETS,
     .    Dnew,Enew,ldos,StructureChanged, idyn, tmdskip, tmdsampling)
C *********************************************************************
C Subroutine to calculate the transport in the extended molecule.
C Gamma-point version for non-polarized or spin-polarized Hamiltonians.
C Written by V. M. Garcia-Suarez, September 2003.
C Departamento de Fisica
C Universidad de Oviedo
C e-mail: victor@condmat.uniovi.es
C ***************************** HISTORY *******************************
C Original version:	June 2003
C 
C Altered to include MPI parallelization	February 2004
C by Alexandre Reily Rocha, Trinity College Dublin
C e-mail: rochaa@tcd.ie
! ------------------------------------
! Add argument inicoor, instead of 0, idyn, tmdskip, and tmdsampling for MD
! Meilin Bai, Dec 2012
C **************************** INPUT **********************************
C integer NspinRealInputMatrix          : Number of spin components (1 or 2)
C integer nuo            : Number of basis orbitals local to node
C integer maxnh          : Maximum number of orbitals interacting  
C integer numh(nuo)      : Number of nonzero elements of each row 
C                          of Hamiltonian matrix
C integer listhptr(nuo)  : Pointer to each row (-1) of the
C                          Hamiltonian matrix
C integer listh(maxnh)   : Nonzero Hamiltonian-matrix element  
C                          column indexes for each matrix row
C real*8  H(maxnh,NspinRealInputMatrix) : Hamiltonian in sparse form
C real*8  S(maxnh)       : Overlap in sparse form
C real*8  temp           : Electronic temperature
C integer nsc(2)         : Number of unit cells along parallel
C                          directions
C character*20 slabel    : System label
C integer NEnergR        : Number of points to be used in the integral
C                          along the real energy axis
C integer NEnergIC       : Number of points to be used in the integral
C                          along the complex energy circle
C integer NEnergIL       : Number of points to be used in the integral
C                          along the complex energy line
C integer NPoles         : Number of poles in the Fermi function
C real*8  Delta          : Small imaginary part
C real*8  EnergLB        : Lower bound for the contour integral
C integer Iter           : SCF iteration
C logical last           : If last SCF iteration, calculate the current
C integer istep          : Molecular dynamics step 
C integer ivv            : Bias potential step 
C real*8  VBias          : Bias potential
C real*8  VInitial       : Initial value of the bias potential
C real*8  VFinal         : Final value of the bias potential
C integer SpinCL         : True = Leads with different spin
C integer NSlices        : Number of bulk slices to include in the
C                          extended molecule
C logical TrCoeff        :Calculate the transmission coefficients?
c logical CalcIETS       :Calculate the IETS signal
C **************************** OUTPUT *********************************
C real*8  Dnew(maxnh,NspinRealInputMatrix) : Output density matrix
C *********************************************************************
      
      use precision
      use sys
      use negfmod
      use mTypes
      use mMatrixUtil
      use mMPI_NEGF
      use set_rhobd
      use ScissorOperator, only : sco,scoSCF,SCOSetHamiltonianBlock,
     .    SCOApplyK
      use mCurrDistTotal,  only : CurrentDistributionMatrix_Vectors,
     .     CurrentDistributionMatrix
      use mNegf_Interface, only : Xa2Xo

      implicit none

C  Parallel Variables
      integer
     .  node, Nodes, BNode, ie, ierror, io, iio, jo, nuog,ispin,
     .  no,indxuo(no)
      integer, intent(in) :: iaorb(no)
#ifdef MPI
      integer
     .  MPIerror,is,istatus(MPI_STATUS_SIZE)
#endif
      integer, save :: maxnhg
      integer, dimension(:), allocatable, save :: numhg,listhptrg,listhg
      integer, dimension(:), allocatable, save :: listig
      double precision, dimension(:,:), allocatable, save ::
     .  Hnew, Dnewloc
      integer, intent(in)  :: na_u,na_s
      double precision,intent(in) :: xa(3,na_s)

      integer
     .  maxnh, nuo, NspinRealInputMatrix, listh(maxnh), numh(nuo),
     .  listhptr(nuo),
     .  nsc(2), NEnergR, NEnergIC, NEnergIL, Npoles,
     .  Iter, istep, inicoor, ivv, SpinCL, NSlices, nuotot,
     .  idyn, tmdskip, tmdsampling
      integer NspinComplexMatrix,NspinBlocks
      integer, save ::
     .  NeqL, NeqR, NeqOffL, NeqOffR
      double precision
     .  H(maxnh,NspinRealInputMatrix), S(maxnh), temp, Delta, 
     .  EnergLB, VBias,
     .  VInitial, VFinal, Dnew(maxnh,NspinRealInputMatrix),kpoint(3),
     .  xij(3,maxnh)
      double precision Enew(maxnh,NspinRealInputMatrix)
      character*20
     .  slabel
      logical
     .  last, TrCoeff,ldos,StructureChanged, CalcIETS
      type(ioType) :: iout

C  Internal variables 
      integer
     .  iuo, j, juo, ind, iu1, iu2, nspinL, nspinR,nn,
     .  juoL, jL, indL, juoR, jR, indR,ind2
      integer, save::
     .   nuoL, nuoR, maxnhL, maxnhR
      integer, allocatable, save::
     .  numdL(:), listdptrL(:), listdL(:), numdR(:), listdptrR(:),
     .  listdR(:),iequivL(:,:),iequivR(:,:),iequivOffL(:,:),
     .  iequivOffR(:,:)
      double complex, allocatable::
     .  rhoG(:,:,:)
      double precision, allocatable, save::
     .  DMbulkL(:,:), DMbulkR(:,:)
      logical, save::
     .  periodic
      character*20
     .  slabelL, slabelR
      external
     .  negfk, memory, io_assign, io_close, dmbk
      type(matrixTypeGeneral), allocatable :: rhogeneralp(:)
      type(matrixTypeGeneral), allocatable :: ematgeneralp(:)
      type(matrixTypeGeneral), allocatable :: hgeneralp(:)
      type(matrixTypeGeneral) :: sgeneralp
      integer nelerow,maxnelerow
      integer nmat
      integer*4:: sc_0,sc_1,sc_r,sc_m
      integer i

      double precision ,allocatable :: DMImag(:,:),OmegaImag(:,:)
      double precision, allocatable :: xo(:,:)
!zrx,
      double precision, allocatable, save::
     .  EDMbulkL(:,:), EDMbulkR(:,:)
!zrx,

#ifdef MPI
      call MPI_Comm_Rank(negf_comm,Node,MPIerror)
      call MPI_Comm_Size(negf_comm,Nodes,MPIerror)
#else
      Node = 0
      Nodes = 1
#endif

      em_Last_SCF_Step=last

      iout%isDebug=.false.
      NspinComplexMatrix=NspinRealInputMatrix
      NspinBlocks=NspinRealInputMatrix
      nmat=nuotot
      if(NspinRealInputMatrix >= 4)then
        NspinBlocks=4
        if(.not.negfon)then
          NspinComplexMatrix=4
        else
          NspinComplexMatrix=1
          nmat=2 * nuotot
        endif
      endif
      EM_NSPINBlocks=NspinBlocks



      allocate(rhogeneralp(NspinComplexMatrix))
      allocate(ematgeneralp(NspinComplexMatrix))
      allocate(hgeneralp(NspinComplexMatrix))

!  write(12347,*)"allocating dmimag"
      if(curr_dist)then
        allocate(DMImag(maxnh,NspinRealInputMatrix))
        allocate(OmegaImag(maxnh,NspinRealInputMatrix))
        DMImag = 0.d0
        OmegaImag = 0.d0
    
        allocate(xo(3,nuotot))
        call Xa2Xo(xa,na_s,xo,nuotot,no,iaorb)
      !    do i=1,nuotot
      !      write(12347,*)"xorb=",i,xo(1,i),xo(2,i),xo(3,i)
      !    enddo
      endif
 
 


      Dnew=0.D0
      Enew=0.D0

      if(.true..or..not.allocated(numhg))then
C Find Hamiltonian and overlap 
C Allocate local memory for global list arrays
      allocate(numhg(nuotot))
      call memory('A','I',nuotot,'emtg')
      allocate(listhptrg(nuotot))
      call memory('A','I',nuotot,'emtg')

C Globalise numh
      do io = 1,nuotot
       call WhichNodeOrb(io,Nodes,BNode)
       if (Node.eq.BNode) then
         call GlobalToLocalOrb(io,Node,Nodes,iio)
C         print*, iio,io, nuo, nuotot, numh(iio), numhg(io)
         numhg(io) = numh(iio)
       endif
#ifdef MPI
       call MPI_Bcast(numhg(io),1,MPI_integer,BNode,
     .   negf_comm,MPIerror)
#endif
      enddo

      listhptrg(1) = 0
      do io = 2,nuotot
        listhptrg(io) = listhptrg(io-1) + numhg(io-1)
      enddo
C Globalse listh
      maxnhg = listhptrg(nuotot)+numhg(nuotot)
      allocate(listhg(maxnhg))
      allocate(listig(maxnhg))
C      call memory('A','I',maxnhg,'diagkp')
      do io = 1,nuotot
        listig(listhptrg(io)+1:listhptrg(io)+numhg(io)) = io
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do jo = 1,numhg(io)
            listhg(listhptrg(io)+1:listhptrg(io)+numhg(io)) = 
     .        listh(listhptr(iio)+1:listhptr(iio)+numh(iio))
          enddo
        endif
#ifdef MPI
        call MPI_Bcast(listhg(listhptrg(io)+1),numhg(io),MPI_integer,
     .    BNode,negf_comm,MPIerror)
#endif
      enddo
      endif

C Create new distribution of H and S
      nuog = nuotot
      call memory('A','D',maxnhg,'diagkp')

      kpoint=0.d0

      if(emtimings)then
       CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif

      call convertmatrixsiestatosmeagolgeneral(H,S,maxnh,numh,
     .    listhptr,listh,nuo,xij,NspinRealInputMatrix,
     .    NspinComplexMatrix,nuotot,nmat,
     .    indxuo,no,kpoint(:),hgeneralp,sgeneralp,rhogeneralp,
     .    ematgeneralp,emforces,
     .    maxnelerow,nnodes_inverse,mynode_inverse,2)

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_convertmatrix_serial',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


      if(nprocs_hs.ne.1)then
        call convertmatrixsiestatosmeagolgeneral(H,S,maxnh,numh,
     .      listhptr,listh,nuo,xij,NspinRealInputMatrix,
     .      NspinComplexMatrix,nuotot,nmat,
     .      indxuo,no,kpoint(:),hgeneralp,sgeneralp,rhogeneralp,
     .      ematgeneralp,emforces,
     .      maxnelerow,nnodes_inverse,mynode_inverse,3)
        do ispin=1,NspinComplexMatrix
          hgeneralp(ispin)%mattype=2
        enddo
 
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_convertmatrix_parallel',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif

       sgeneralp%mattype=2

!        call PrintMatrixGeneral(sgeneralp,"sgeneralp=",iout)
        call MaxDifferenceMatrixCRS_CRSP(sgeneralp,sgeneralp,
     .      mynode_inverse,"(sgeneralp) ",iout)

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_differences_s',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif


        do ispin=1,NspinComplexMatrix
          call MaxDifferenceMatrixCRS_CRSP(hgeneralp(ispin),
     .        hgeneralp(ispin),mynode_inverse,"(hgeneraludp) ",iout)

          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           't_differences_h',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif




        enddo
      endif



      kpointmod=0D0
      ikpmod=1
      ikpmodK=1
      wkmod=1.0D0


      if (sco.and.(last.or.scoSCF)) then
        SCOSetHamiltonianBlock=.true.
        call SCOApplyK(NspinRealInputMatrix,maxnh,H,S,xij,
     .      kpoint(:),mynode_negfo,nnodes_negfo,nuo,
     .      no,numh,listhptr,listh,indxuo,hgeneralp,sgeneralp,
     .      NspinComplexMatrix)
      else
        SCOSetHamiltonianBlock=.false.
      endif



      call NEGFK(nuotot,NspinBlocks,NspinComplexMatrix,nsc,NEnergR,
     .    NEnergIL, NEnergIC,
     .    Npoles,istep,inicoor,
     .    ivv,Iter,1,SpinCL,VInitial,VFinal,temp,Delta,
     .    VBias,EnergLB,last,NSlices,TrCoeff,CalcIETS,
     .    1,kpoint(1),slabel,1.0D0,
     .    hgeneralp,sgeneralp,rhogeneralp,ematgeneralp,ldos,
     &    idyn, tmdskip, tmdsampling)

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_negfk',(sc_1-sc_0)*1.0d0/sc_r
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

C New density matrix      

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_before_rhonegf_to_rhosiesta',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif

      call rhonegf_to_rhosiesta(maxnelerow,NspinRealInputMatrix,
     .    NspinComplexMatrix,rhogeneralp,Dnew,
     .    maxnh,no,nuo,nuotot,nmat, listh, numh, listhptr,indxuo,
     .    kpoint,xij,1.0D0,timereversal)
      if(emforces)
     .    call rhonegf_to_rhosiesta(maxnelerow,NspinRealInputMatrix,
     .    NspinComplexMatrix,ematgeneralp,
     .    Enew, maxnh,no,nuo,nuotot,nmat, listh, numh, listhptr,indxuo,
     .    kpoint,xij,1.0D0,timereversal)

      if(curr_dist)then
        call rhonegf_to_rhosiestaImaginary(maxnelerow,
     .      NspinRealInputMatrix,  NspinComplexMatrix,rhogeneralp,
     .      DMImag, maxnh,no,nuo,nuotot,nmat, listh, numh,
     .      listhptr,indxuo, kpoint,xij,1.0D0,
     .      timereversal)
        call rhonegf_to_rhosiestaImaginary(maxnelerow,
     .      NspinRealInputMatrix,  NspinComplexMatrix,ematgeneralp,
     .      OmegaImag, maxnh,no,nuo,nuotot,nmat, listh, numh,
     .      listhptr,indxuo, kpoint,xij,1.0D0,
     .      timereversal)
      endif





      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_rhonegf_to_rhosiesta',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


      if(mynode_inverse.eq.0)then
        do ispin=1,NspinComplexMatrix
          call DestroyMatrixGeneral(rhogeneralp(ispin),"emtg",iout)
          call DestroyMatrixGeneral(ematgeneralp(ispin),"emtg",iout)
        enddo
      endif

      if(curr_dist)then
        call CurrentDistributionMatrix(DMImag,OmegaImag,H,S,maxnh,numh,
     .      listhptr,listh,  nuo,xij,xo,NspinRealInputMatrix,
     .      NspinComplexMatrix,nuotot,  indxuo,no,VBias,.true.)
        deallocate(DMImag,OmegaImag,xo)
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtg1',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


      if (Node .EQ. 0) then
!from here everything is executed on node 0
! Add bulk density matrix: this is done on node 0, and then sent to the other nodes

        if (Iter.eq.1 .and. istep.eq. inicoor .and. ivv.eq.0) then
          call io_assign(iu1)
          open(iu1,file='bulklft.DAT',status='old')
          read(iu1,*) slabelL, nuoL, nspinL, maxnhL

          call io_assign(iu2)
          open(iu2,file='bulkrgt.DAT',status='old')
          read(iu2,*) slabelR, nuoR, nspinR, maxnhR

          allocate(DMbulkL(maxnhL,NspinRealInputMatrix))
          call memory('A','Z',maxnhL*NspinRealInputMatrix,'emtk')
          allocate(DMbulkR(maxnhR,NspinRealInputMatrix))
          call memory('A','Z',2*maxnhR*NspinRealInputMatrix,'emtk')

          allocate(numdL(nuoL),listdptrL(nuoL),listdL(maxnhL))
          call memory('A','I',2*nuoL+maxnhL,'emtk')
          allocate(numdR(nuoR),listdptrR(nuoR),listdR(maxnhR))
          call memory('A','I',2*nuoR+maxnhR,'emtk')

          DMbulkL = 0.d0
          call dmbk(nuoL,maxnhL,nspinL,slabelL,
     .     numdL,listdL,listdptrL,DMbulkL,1)
          DMbulkR = 0.d0
          call dmbk(nuoR,maxnhR,nspinR,slabelR,
     .     numdR,listdR,listdptrR,DMbulkR,1)


          if(emforces)then
            allocate(EDMbulkL(maxnhL,NspinRealInputMatrix))
            allocate(EDMbulkR(maxnhR,NspinRealInputMatrix))
            EDMbulkL = 0.d0
            call dmbk(nuoL,maxnhL,nspinL,slabelL,
     .       numdL,listdL,listdptrL,EDMbulkL,2)
            EDMbulkR = 0.d0
            call dmbk(nuoR,maxnhR,nspinR,slabelR,
     .       numdR,listdR,listdptrR,EDMbulkR,2)
          endif

          call io_close(iu1)
          call io_close(iu2)

          periodic = Set_RhoBoundaryOverlap_Leads

          allocate(iequivL(maxnhL,2),iequivR(maxnhR,2))
          if (periodic) then
           allocate(iequivOffL(maxnhL,2),iequivOffR(maxnhR,2))
          else
           allocate(iequivOffL(1,1),iequivOffR(1,1))
          endif

        endif


        if ((Iter.eq.1).and.(StructureChanged)) then


          NeqL=0
          NeqOffL=0
          NeqR=0
          NeqOffR=0

          nn = nsc(1)*nsc(2)

          call rho_boundaryL(nn,nuoL,maxnhL,listdptrL,listdL,numdL,
     .        nuotot,
     .        maxnhg,numhg,listhptrg,listhg,neqL,iequivL,nuoR,
     .        periodic,NeqOffL,iequivOffL)

          call rho_boundaryR(nn,nuoR,maxnhR,listdptrR,listdR,numdR,
     .        nuotot,
     .        maxnhg,numhg,listhptrg,listhg,neqR,iequivR,nuoL,
     .        periodic,NeqOffR,iequivOffR)

        endif

      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtg2',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif



!up to here everything is executed on node 0, only node 0 has the full DM
#ifdef MPI
       if ((Iter.eq.1).and.(StructureChanged)) then
        CALL MPI_BCAST(periodic,1,MPI_LOGICAL,0,
     &      negf_comm,MPIerror)
        call MPI_Bcast(NeqL,1,MPI_integer,0, negf_comm,MPIerror)
        call MPI_Bcast(NeqR,1,MPI_integer,0, negf_comm,MPIerror)
        if(periodic)then
          call MPI_Bcast(NeqOffL,1,MPI_integer,0, negf_comm,MPIerror)
          call MPI_Bcast(NeqOffR,1,MPI_integer,0, negf_comm,MPIerror)
        endif
        call MPI_Bcast(maxnhl,1,MPI_integer,0, negf_comm,MPIerror)
        call MPI_Bcast(maxnhr,1,MPI_integer,0, negf_comm,MPIerror)

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_emtg2b',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif



        if(node.ne.0)then
          if(.not.allocated(iequivL))allocate(iequivL(maxnhL,2))
          if(.not.allocated(iequivR))allocate(iequivR(maxnhR,2))
          if(.not.allocated(DMbulkL))allocate(DMbulkL(maxnhl,
     .        NspinRealInputMatrix))
          if(.not.allocated(DMbulkR))allocate(DMbulkR(maxnhr,
     .        NspinRealInputMatrix))
          if(emforces)then
            if(.not.allocated(EDMbulkL))allocate(EDMbulkL(maxnhl,
     .          NspinRealInputMatrix))
            if(.not.allocated(EDMbulkR))allocate(EDMbulkR(maxnhr,
     .          NspinRealInputMatrix))
          endif
          if(periodic)then
            if(.not.allocated(iequivoffL))allocate(iequivoffL(maxnhL,2))
            if(.not.allocated(iequivoffR))allocate(iequivoffR(maxnhR,2))
          endif
        endif


        call MPI_Bcast(iequivL(1,1),maxnhL*2,MPI_integer,0, 
     .      negf_comm,MPIerror)

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_emtg2d',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif


        call MPI_Bcast(iequivR(1,1),maxnhR*2,MPI_integer,0, 
     .      negf_comm,MPIerror)

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         't_emtg2e',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif



        if(periodic)then
          call MPI_Bcast(iequivoffL(1,1),maxnhL*2,MPI_integer,0, 
     .        negf_comm,MPIerror)
          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           't_emtg2f',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif

          call MPI_Bcast(iequivoffR(1,1),maxnhR*2,MPI_integer,0, 
     .        negf_comm,MPIerror)

          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           't_emtg2g',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif


        endif
      endif
#endif
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtg3',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif



!zrx,because all variables have been declared using 'save'
!DO NOT put this block into IOnode block
      IF((ivv.eq.0).and.(istep.eq. inicoor).and.(iter.eq.1)) THEN
#ifdef MPI
#ifdef NODAT
         call MPI_Bcast(DMbulkL(1,1),maxnhL*NspinRealInputMatrix,
     .       MPI_double_precision,
     .     0,negf_comm,MPIerror)
         call MPI_Bcast(DMbulkR(1,1),maxnhR*NspinRealInputMatrix,
     .       MPI_double_precision,
     .     0,negf_comm,MPIerror)
         if(emforces)then
           call MPI_Bcast(EDMbulkL(1,1),maxnhL*NspinRealInputMatrix,
     .         MPI_double_precision,0,negf_comm,MPIerror)
           call MPI_Bcast(EDMbulkR(1,1),maxnhR*NspinRealInputMatrix,
     .         MPI_double_precision,0,negf_comm,MPIerror)
         endif
#else
         call MPI_Bcast(DMbulkL(1,1),maxnhL*NspinRealInputMatrix,
     .       DAT_double,
     .     0,negf_comm,MPIerror)
         call MPI_Bcast(DMbulkR(1,1),maxnhR*NspinRealInputMatrix,
     .       DAT_double,
     .     0,negf_comm,MPIerror)
         if(emforces)then
           call MPI_Bcast(EDMbulkL(1,1),maxnhL*NspinRealInputMatrix,
     .         DAT_double,
     .         0,negf_comm,MPIerror)
           call MPI_Bcast(EDMbulkR(1,1),maxnhR*NspinRealInputMatrix,
     .         DAT_double, 0,negf_comm,MPIerror)
         endif
#endif
#endif
      ENDIF

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtg4',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif



      if(.not.ldos.and..not.last)then
        call boundary_elements_mpi(neql,maxnhl,iequivL,maxnhg,
     .      listhg,listig,listhptrg,nuotot,dnew,DMbulkL,
     .      NspinRealInputMatrix,
     .      maxnh,
     .      numh,nuo,listhptr,listh)
        call boundary_elements_mpi(neqr,maxnhr,iequivR,maxnhg,
     .      listhg,listig,listhptrg,nuotot,dnew,DMbulkR,
     .      NspinRealInputMatrix,
     .      maxnh,
     .      numh,nuo,listhptr,listh)
        if(periodic)then
          call boundary_elements_mpi(neqoffl,maxnhl,iequivoffL,maxnhg,
     .        listhg,listig,listhptrg,nuotot,dnew,DMbulkL,
     .        NspinRealInputMatrix,
     .        maxnh,numh,nuo,listhptr,listh)
          call boundary_elements_mpi(neqoffr,maxnhr,iequivoffR,maxnhg,
     .        listhg,listig,listhptrg,nuotot,dnew,DMbulkR,
     .        NspinRealInputMatrix,
     .        maxnh,numh,nuo,listhptr,listh)
        endif
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtg5',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


      
!       if(last .and. SetEBD) then
       if(emforces) then
         call boundary_elements_mpi(neql,maxnhl,iequivL,maxnhg,
     .       listhg,listig,listhptrg,nuotot,Enew,EDMbulkL,
     .       NspinRealInputMatrix,maxnh,
     .       numh,nuo,listhptr,listh,DMbulkL,VBias,'L')
         call boundary_elements_mpi(neqr,maxnhr,iequivR,maxnhg,
     .       listhg,listig,listhptrg,nuotot,Enew,EDMbulkR,
     .       NspinRealInputMatrix,maxnh,
     .       numh,nuo,listhptr,listh,DMbulkR,VBias,'R')
        if(periodic)then
          call boundary_elements_mpi(neqoffl,maxnhl,iequivoffL,maxnhg,
     .        listhg,listig,listhptrg,nuotot,Enew,EDMbulkL,
     .        NspinRealInputMatrix,maxnh,
     .        numh,nuo,listhptr,listh,DMbulkL,VBias,'L')
          call boundary_elements_mpi(neqoffr,maxnhr,iequivoffR,maxnhg,
     .        listhg,listig,listhptrg,nuotot,Enew,EDMbulkR,
     .        NspinRealInputMatrix,maxnh,
     .        numh,nuo,listhptr,listh,DMbulkR,VBias,'R')
         endif
       endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtg6',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif



C Free local arrays      
      if(.true..or..false.)then
        call memory('D','I',size(numhg),'emtg')
        deallocate(numhg)
        call memory('D','I',size(listhptrg),'emtg')
        deallocate(listhptrg)
        call memory('D','I',size(listhg),'emtg')
        deallocate(listhg)
        deallocate(listig)
      endif

!      write(12347,*)"emtimings=",emtimings
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtkexit',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


      deallocate(rhogeneralp,ematgeneralp,hgeneralp)
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       't_emtg7',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif




      return
      end


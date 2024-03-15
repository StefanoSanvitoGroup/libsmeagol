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
! UPON COMPLETION OF THE "SMEAGOL ACADEMIC LICENSE".
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
      subroutine hsleads(trnspin, nspin, dnuoL, dnuoR, nsc, iter,
     .                   istep, inicoor, iv, gamma, ik, kpoint, temp,
     .                   H0L, H1L, S0L, S1L, H0R, H1R, S0R, S1R)
C *******************************************************************
C Reads and calculates the Hamiltonians, overlaps and density matrices
C of the leads.
C Written by V. M. Garcia-Suarez.
C Departamento de Fisica
C Universidad de Oviedo
C e-mail: victor@condmat.uniovi.es
C ***************************** HISTORY *******************************
C Original version:	June 2003
! --------------------------------------------
! Add argument inicoor instead of 0 for MD
! Meilin Bai, Dec 2012
!
C *****************************INPUT*********************************
C For simplicity, I list only the left lead variables
C integer trnspin              : True value of the spin
C integer nspin                : Number of spin components
C integer dnuoL                : Number of basis orbitals in the left 
C                                lead, including spin components
C integer dnuoR                : Number of basis orbitals in the right 
C                                lead, including spin components
C integer nsc(2)               : Number of unit cells along parallel
C                                directions
C integer iter                 : Scf iteration in SIESTA
C integer istep                : Molecular dynamics iteration
C integer iv                   : Bias potential iteration
C logical gamma                : Calculation with parallel k-points
C integer ik                   : k-point index
C real*8  kpoint(3)            : Current parallel k point
C real*8  temp                 : Electronic temperature
C *****************************OUTPUT********************************
C complex*8  H0L(dnuoL,dnuoL,nspin) : Hamiltonian in the unit cell of the
C                                left lead
C complex*8  H1L(dnuoL,dnuoL,nspin) : Conexion between unit cells in the
C                                left lead
C complex*8  S0L(dnuoL,dnuoL)    : Overlaps in the unit cell of the
C                                left lead
C complex*8  S1L(dnuoL,dnuoL)    : Overlaps between unit cells in the
C                                left lead
C complex*8  H0R(dnuoR,dnuoR,nspin) : Hamiltonian in the unit cell of the
C                                right lead
C complex*8  H1R(dnuoR,dnuoR,nspin) : Conexion between unit cells in the
C                                right lead
C complex*8  S0R(dnuoR,dnuoR)    : Overlaps in the unit cell of the
C                                right lead
C complex*8  S1R(dnuoR,dnuoR)    : Overlaps between unit cells in the
C                                right lead
C *******************************************************************

      use negfmod, only : ThetaLeadsL,PhiLeadsL,ThetaLeadsR,PhiLeadsR
#ifdef MPI
      use mpi_siesta, only : MPI_COMM_WORLD
#endif
      use ScissorOperator, only: SCOLeads,SCO_Leads


      implicit none
      
      integer
     .  trnspin, nspin, dnuoL, dnuoR, nsc(2), iter, istep, iv, ik,
     &  inicoor
      double precision
     .  kpoint(3), temp
      double complex
     .  H0L(dnuoL,dnuoL,nspin), H1L(dnuoL,dnuoL,nspin),
     .  S0L(dnuoL,dnuoL), S1L(dnuoL,dnuoL), H0R(dnuoR,dnuoR,nspin),
     .  H1R(dnuoR,dnuoR,nspin), S0R(dnuoR,dnuoR), S1R(dnuoR,dnuoR)
      logical
     .  gamma
      character
     .  slabelL*20, slabelR*20

C Internal variables
      integer, save :: nspinL,nspinR
      integer
     .  nscL(2), nscR(2),
     .  iu, iu1, iu2, io, iuo, ind, j, ispin,nspinMin
      integer, save ::
     .  nuoL, nuoR, noL, noR, maxnhL, maxnhR
      integer, allocatable, save::
     .  numhL(:), numhR(:), listhptrL(:), listhptrR(:), indxuoL(:),
     .  indxuoR(:), listhL(:), listhR(:) 
      double precision
     .  efL, efR, tempR, tempL
      double precision, allocatable, save::
     .  xijL(:,:), xijR(:,:), SL(:), SR(:), HL(:,:), HR(:,:)
      character
     .  paste*25
      integer MPIerror
      
      external
     .  io_assign, io_close, hsl, hslk

      if (iter.eq.1 .and. istep.eq. inicoor
     .      .and. iv.eq.0 .and. ik.eq.1) then
C Read data
        call io_assign(iu1)
        open(iu1,file='bulklft.DAT',status='old')
        read(iu1,*) slabelL, nuoL, nspinL, maxnhL, efL, tempL,
     .              nscL(1), nscL(2), noL

        call io_assign(iu2)
        open(iu2,file='bulkrgt.DAT',status='old')
        read(iu2,*) slabelR, nuoR, nspinR, maxnhR, efR, tempR,
     .              nscR(1), nscR(2), noR

c Allocate arrays
        allocate(numhL(nuoL),listhptrL(nuoL))
        call memory('A','I',2*nuoL,'hsleads')
        allocate(indxuoL(noL))
        call memory('A','I',noL,'hsleads')
        allocate(listhL(maxnhL))
        call memory('A','I',maxnhL,'hsleads')
        allocate(xijL(3,maxnhL))
        call memory('A','D',3*maxnhL,'hsleads')

        allocate(numhR(nuoR),listhptrR(nuoR))
        call memory('A','I',2*nuoR,'hsleads')
        allocate(indxuoR(noR))
        call memory('A','I',noR,'hsleads')
        allocate(listhR(maxnhR))
        call memory('A','I',maxnhR,'hsleads')
        allocate(xijR(3,maxnhR))
        call memory('A','D',3*maxnhR,'hsleads')

        allocate(SL(maxnhL))
        call memory('A','D',maxnhL,'hsleads')
        allocate(SR(maxnhR))
        call memory('A','D',maxnhR,'hsleads')
        allocate(HL(maxnhL,nspinl))
        call memory('A','D',maxnhL*nspinl,'hsleads')
        allocate(HR(maxnhR,nspinr))
        call memory('A','D',maxnhR*nspinr,'hsleads')

c Read data of the left lead
        do iuo = 1, nuoL
          read(iu1,*) numhL(iuo), listhptrL(iuo)
        enddo
        do io = 1, noL
          read(iu1,*) indxuoL(io)
        enddo
        do iuo = 1, nuoL
          do j = 1, numhL(iuo)
            ind = listhptrL(iuo) + j
            read(iu1,*) listhL(ind)
            if (.not.gamma) then
              read(iu1,*) xijL(1,ind), xijL(2,ind), xijL(3,ind)
            endif
          enddo
        enddo
      
c Read data of the right lead
        do iuo = 1, nuoR
          read(iu2,*) numhR(iuo), listhptrR(iuo)
        enddo
        do io = 1, noR
          read(iu2,*) indxuoR(io)
        enddo
        do iuo = 1, nuoR
          do j = 1, numhR(iuo)
            ind = listhptrR(iuo) + j
            read(iu2,*) listhR(ind)
            if (.not.gamma) then
              read(iu2,*) xijR(1,ind), xijR(2,ind), xijR(3,ind)
            endif
          enddo
        enddo
      
        call io_close(iu1)
        call io_close(iu2)

c Compare supercell of leads and EM
        if (nsc(1).ne.nscL(1) .or. nsc(2).ne.nscL(2)) then
          write(6,'(a)') 'ERROR: The left supercell along parallel'
          write(6,'(a)') 'directions is different of the supercell'
          write(6,'(a)') 'of the EM. Change the size of it'
          write(6,'(a,2i4)') 'nsc = ', nsc(1), nsc(2)
          write(6,'(a,2i4)') 'nscL = ', nscL(1), nscL(2)
          write(6,'(a)') 'exiting...'
#ifdef MPI
          call MPI_Abort(MPI_Comm_World, 1, MPIerror)
#endif
          stop
        endif
        if (nsc(1).ne.nscR(1) .or. nsc(2).ne.nscR(2)) then
          write(6,'(a)') 'ERROR: The right supercell along parallel'
          write(6,'(a)') 'directions is different of the supercell'
          write(6,'(a)') 'of the EM. Change the size of it'
          write(6,'(a,2i4)') 'nsc = ', nsc(1), nsc(2)
          write(6,'(a,2i4)') 'nscR = ', nscR(1), nscR(2)
          write(6,'(a)') 'exiting...'
#ifdef MPI
          call MPI_Abort(MPI_Comm_World, 1, MPIerror)
#endif
          stop
        endif
          
c Verify if the spin at the leads and the EM is the same
        if(.false.)then
          if (nspinL.ne.trnspin) then
            write(6,'(a)') 'ERROR: The spin at the left lead is not'
            write(6,'(a)') 'the same as in the extended molecule'
            stop
          else if (nspinR.ne.trnspin) then
            write(6,'(a)') 'ERROR: The spin at the right lead is not'
            write(6,'(a)') 'the same as in the extended molecule'
            stop
          endif
        endif

c Check the temperature
        if (tempL-temp.gt.1.d-4) then
          write(6,'(a)')'WARNING: The temperature at the left lead is'
          write(6,'(a)')'not the same as in the extended molecule'
        endif
        if (tempR-temp.gt.1.d-4) then
          write(6,'(a)')'WARNING: The temperature at the right lead is'
          write(6,'(a)')'not the same as in the extended molecule'
        endif

c Initialize overlaps and Hamiltonians
        SL = 0.d0
        SR = 0.d0
        HL = 0.d0
        HR = 0.d0

c Read overlaps and Hamiltonians
        call io_assign(iu)
        open(iu,file=paste(slabelL,'.HST'),status='old')
        do iuo = 1, nuoL
          do j = 1, numhL(iuo) 
            ind = listhptrL(iuo) + j
            read(iu,*) SL(ind)
            read(iu,*) HL(ind,:)
          enddo
        enddo
        call io_close(iu)
        
        call io_assign(iu)
        open(iu,file=paste(slabelR,'.HST'),status='old')
        do iuo = 1, nuoR
          do j = 1, numhR(iuo) 
            ind = listhptrR(iuo) + j
            read(iu,*) SR(ind)
            read(iu,*) HR(ind,:)
          enddo
        enddo
        call io_close(iu)



        if(ThetaLeadsL.ne.0.0D0.or.PhiLeadsL.ne.0.0D0)then
          call rotateDM(HL,maxnhL,nspinl,ThetaLeadsL,PhiLeadsL)
        endif

        if(ThetaLeadsR.ne.0.0D0.or.PhiLeadsR.ne.0.0D0)then
          call rotateDM(HR,maxnhR,nspinr,ThetaLeadsR,PhiLeadsR)
        endif

        if(gamma) then 
c Calculate S0, S1, H0 and H1 (left and right)

          call hsl(nspinl, nspin, nuoL, dnuoL, noL, maxnhL,
     .             numhL, listhptrL, indxuoL, listhL, SL, HL,
     .             S0L, S1L, H0L, H1L)
          call hsl(nspinr, nspin, nuoR, dnuoR, noR, maxnhR,
     .             numhR, listhptrR, indxuoR, listhR, SR, HR,
     .             S0R, S1R, H0R, H1R)
          if(SCO_Leads) then
            call SCOLeads(S0L,H0L, dnuoL,nspin,0,0)
            call SCOLeads(S0R,H0R, dnuoR,nspin,0,0)
          endif
        endif
          
      endif

c Calculate S0, S1, H0 and H1 (left and right) if not gamma point
      if(.not.gamma) then
        call hslk(nspinl, nspin, nuoL, dnuoL, noL, maxnhL, numhL,
     .            listhptrL, indxuoL, listhL, xijL, kpoint, nsc,
     .            SL, HL, S0L, S1L, H0L, H1L)
        call hslk(nspinr, nspin, nuoR, dnuoR, noR, maxnhR, numhR,
     .            listhptrR, indxuoR, listhR, xijR, kpoint, nsc,
     .            SR, HR, S0R, S1R, H0R, H1R)
      endif

      nspinMin=nspin
      if(nspin>2)nspinMin=2

      do ispin=1,nspinMin
        H0L(:,:,ispin)=0.5D0 * (H0L(:,:,ispin)+
     .      DCONJG(TRANSPOSE(H0L(:,:,ispin))))
        H0R(:,:,ispin)=0.5D0 * (H0R(:,:,ispin)+
     .      DCONJG(TRANSPOSE(H0R(:,:,ispin))))
      enddo
      S0L=0.5D0 * (S0L+DCONJG(TRANSPOSE(S0L)))
      S0R=0.5D0 * (S0R+DCONJG(TRANSPOSE(S0R)))

      if(SCO_Leads.and..not.gamma) then
        call SCOLeads(S0L,H0L, dnuoL,nspin,0,0)
        call SCOLeads(S0R,H0R, dnuoR,nspin,0,0)
      endif

     

      end


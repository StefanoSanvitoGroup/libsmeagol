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
!                   READ_EF_BSS  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
      subroutine read_ef_bss(ef_bss,nleads,nspin,ef_lead)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use precision
      use mMPI_NEGF

      implicit none

      integer   iif, iostatv
      integer   nleads,nspin
      double precision ef_bss(nleads,nspin),ef_bssn(nleads,nspin),
     .    ef_bsso(nleads,nspin),ef_lead
      real*8    wmix, wmix_old, wmix_new
#ifdef MPI
      integer  MPIerror
#endif

      iif=12349

      ef_bsso=ef_bss
      if (myhead.eq.0) then
        open(unit=iif,file='DM.EF.dat',iostat=iostatv)
        if (iostatv .eq. 0 ) then
          read(iif,*,iostat=iostatv)ef_bssn
          ef_bssn=(ef_bssn / 13.6057D0 ) + ef_lead
          if (iostatv .eq. 0 ) then
            ef_bss=ef_bssn
          endif
!          write(12347,*) 'initial_ef:  = ',ef_bss
          close(iif,iostat=iostatv) 
          if (iostatv .ne. 0 )
     .        write(6,'("initial ef: error in close iif")')
        endif
      endif

#ifdef MPI
#ifdef NODAT
      call MPI_Bcast(ef_bss(1,1),nleads * nspin,MPI_double_precision,0,
     .    inverseheads_comm,MPIerror)
#else
      call MPI_Bcast(ef_bss(1,1),nleads * nspin,DAT_double,0,
     .    inverseheads_comm,MPIerror)
#endif
#endif

      end subroutine read_ef_bss


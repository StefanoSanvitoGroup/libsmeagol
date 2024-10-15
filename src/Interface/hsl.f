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
      subroutine hsl(trnspin, nspin, nuo, dnuo, no, maxnh,
     .               numh, listhptr, indxuo, listh, S, H,
     .               S0, S1, H0, H1)
C *********************************************************************
C Subroutine which calculates S0, S1, H0 and H1 
C Coded by V. M. Garcia-Suarez
C Departamento de Fisica
C Universidad de Oviedo
C e-mail: victor@condmat.uniovi.es
C ***************************** HISTORY *******************************
C Original version:	June 2003
C **************************** INPUT **********************************
C integer trnspin             : True value of the spin
C integer nspin               : Number of spin components
C integer nuo                 : Number of basis orbitals in unit cell
C integer dnuo                : Number of basis orbitals in unit cell
C                               including spin components
C integer no                  : Number of basis orbitals in supercell
C integer maxnh               : Maximum number of orbitals interacting
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to each row (-1) of the
C                               Hamiltonian matrix
C integer indxuo(no)          : Index of equivalent orbital in unit cell
C integer listh(maxnh)        : Nonzero Hamiltonian-matrix element
C                               column indexes for each matrix row
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  H(maxnh,trnspin)    : Hamiltonian in sparse form
C **************************** OUTPUT *********************************
C complex*8  S0(dnuo,dnuo)    : Overlap in the unit cell
C complex*8  S1(dnuo,dnuo)    : Overlap that connects unit cells along
C                               the transport direction
C complex*8  H0(dnuo,dnuo,nspin) : Hamiltonian in the unit cell
C complex*8  H1(dnuo,dnuo,nspin) : Hamiltonian that connects unit cells
C                               along the transport direction
C *********************************************************************
C
C Modules
C
      use sys 

      implicit none

      integer
     .  trnspin, nspin, nuo, dnuo, no, maxnh, numh(nuo),
     .  listhptr(nuo), indxuo(no), listh(maxnh)
      double precision
     .  S(maxnh), H(maxnh,trnspin)
      double complex
     .  S0(dnuo,dnuo), S1(dnuo,dnuo), H0(dnuo,dnuo,nspin),
     .  H1(dnuo,dnuo,nspin),Sm1(dnuo,dnuo),Hm1(dnuo,dnuo,nspin)

C  Internal variables .............................................
      integer
     .  iuo, juo, j, jo, ind
      double complex
     .  ii
      
      external
     .  memory
C ....................

      ii = (0.d0,1.d0)

c Find S0, S1, H0 and H1
      S0 = 0.0d0
      S1 = 0.0d0
      Sm1 = 0.0d0
      H0 = 0.0d0
      H1 = 0.0d0
      Hm1 = 0.0d0
      do iuo = 1, nuo
        do j = 1, numh(iuo)
          ind = listhptr(iuo) + j
          jo = listh(ind)
          juo = indxuo(jo)

          if(jo.le.nuo) then
            S0(iuo,juo) = S0(iuo,juo) + S(ind)* 0.5D0
            S0(juo,iuo) = S0(juo,iuo) + S(ind)* 0.5D0
            H0(iuo,juo,1) = H0(iuo,juo,1) + H(ind,1)* 0.5D0
            H0(juo,iuo,1) = H0(juo,iuo,1) + H(ind,1)* 0.5D0
            if (trnspin > 1) then
              H0(iuo,juo,2) = H0(iuo,juo,2) + H(ind,2)* 0.5D0
              H0(juo,iuo,2) = H0(juo,iuo,2) + H(ind,2)* 0.5D0
              if (trnspin > 2) then
                H0(iuo,juo,3) = H0(iuo,juo,3) + 
     .              (H(ind,3)+ii * H(ind,4))

                if(trnspin > 4)then
                  H0(iuo,juo,1) = H0(iuo,juo,1) + ii *H(ind,5)*
     .                  0.5D0
                  H0(iuo,juo,2) = H0(iuo,juo,2) + ii *H(ind,6)*
     .                  0.5D0

                  H0(juo,iuo,1) = H0(juo,iuo,1) - ii *H(ind,5)*
     .                  0.5D0
                  H0(juo,iuo,2) = H0(juo,iuo,2) - ii *H(ind,6)*
     .                  0.5D0

                endif

              endif
            endif
          else if(jo.gt.nuo.and.jo.le.nuo*2) then
            S1(iuo,juo) = S1(iuo,juo) + S(ind)
            H1(iuo,juo,1) = H1(iuo,juo,1) + H(ind,1)
            if (trnspin > 1) then
              H1(iuo,juo,2) = H1(iuo,juo,2) + H(ind,2)
              if (trnspin > 2) then
                H1(iuo,juo,3) = H1(iuo,juo,3) +
     .              (H(ind,3)+ii * H(ind,4))
                H1(iuo,juo,4) = H1(iuo,juo,4) +
     .              (H(ind,3)-ii * H(ind,4))
                if(trnspin > 4)then
                  H1(iuo,juo,1) = H1(iuo,juo,1) + ii * H(ind,5)
                  H1(iuo,juo,2) = H1(iuo,juo,2) + ii * H(ind,6)
                endif
              endif
            endif
          else if(jo.gt.nuo*2.and.jo.le.nuo*3) then
            Sm1(iuo,juo) = Sm1(iuo,juo) + S(ind)
            Hm1(iuo,juo,1) = Hm1(iuo,juo,1) + H(ind,1)
            if (trnspin > 1) then
              Hm1(iuo,juo,2) = Hm1(iuo,juo,2) + H(ind,2)
              if (trnspin > 2) then
                Hm1(iuo,juo,3) = Hm1(iuo,juo,3) +
     .              (H(ind,3)+ii * H(ind,4))
                Hm1(iuo,juo,4) = Hm1(iuo,juo,4) +
     .              (H(ind,3)-ii * H(ind,4))
                if(trnspin > 4)then
                  Hm1(iuo,juo,1) = Hm1(iuo,juo,1) + ii * H(ind,5)
                  Hm1(iuo,juo,2) = Hm1(iuo,juo,2) + ii * H(ind,6)
                endif
              endif
            endif

          endif

        enddo
      enddo

      if(trnspin>2) then
        H0(:,:,4) = dconjg(transpose(H0(:,:,3)))
      endif

      s1=0.5D0 * (S1+dconjg(transpose(Sm1)))

      H1(:,:,1)=0.5D0 * (H1(:,:,1)+dconjg(transpose(Hm1(:,:,1))))
      if (trnspin > 1) then
        H1(:,:,2)=0.5D0 * (H1(:,:,2)+dconjg(transpose(Hm1(:,:,2))))
        if (trnspin > 2) then

!          dh34 and dh43 should be 0
!          Sm1=H1(:,:,3)-dconjg(transpose(hm1(:,:,4)))
!          write(*,*)"dh34=",maxval(abs(sm1))
!          Sm1=H1(:,:,4)-dconjg(transpose(hm1(:,:,3)))
!          write(*,*)"dh43=",maxval(abs(sm1))

          H1(:,:,3)=0.5D0 * (H1(:,:,3)+dconjg(transpose(Hm1(:,:,4))))
          H1(:,:,4)=0.5D0 * (H1(:,:,4)+dconjg(transpose(Hm1(:,:,3))))

        endif
      endif

! Hermiticity of H0 (non-collinear spin)
      s0=0.5D0 * (S0+dconjg(transpose(S0)))
      H0(:,:,1)=0.5D0 * (H0(:,:,1)+dconjg(transpose(H0(:,:,1))))
      if (trnspin > 2) then
        H0(:,:,2)=0.5D0 * (H0(:,:,2)+dconjg(transpose(H0(:,:,2))))
      endif



      return
      end



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
!                   SET_GAMMAMP,
!                   SET_GAMMAMP_LR,
!                   SET_VALUES_SIGMAMP,
!                   SET_VALUES_GAMMAMP,
!                   SET_VALUES_SIGMAMPLR,
!                   FINDNNZ_SUBBLOCK,
!                   FINDNNZ_SUBBLOCK_DENSESIGMALR  
! AND
! THE MODULE
!                   MONBOUNDSTATES  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
module mONBoundStates
 use mConstants
 use mTypes
 use mMatrixUtil

 implicit none
  private

  public :: set_gammamp
  public :: set_gammamp_lr


contains

subroutine set_gammamp(n1,nl,nr,nleads,nbss,nebss,deltabss,gammamp,sigmamp,sgeneral)

  use negfmod, only: nebss_bs,deltabss_bs

  character(len=*), parameter :: sMyName="set_gammamp"
  type(matrixTypeGeneral) :: gammamp(nleads)
  type(matrixTypeGeneral) :: sigmamp(nleads)
  type(matrixTypeGeneral) :: sgeneral
  integer, intent(in) :: n1,nl,nr,nleads,nbss
  integer nebss(nleads,2)
  double precision deltabss(nleads)
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)

  type(ioType) :: io
  integer ispin,i,j,i1,i2,i3,i4,ind,ii,jj,nnz,mattype,imin,imax,jmin,jmax

  io%iout=778
  io%isDebug= .false.
  io%memCount=0.0_kdp

  deltabss=0D0
  do i=1,nbss
!    write(12347,*)"i,nebss_bs,deltabss_bs=",i,nebss_bs(i,1),nebss_bs(i,2),deltabss_bs(i)
    nebss(i,:)=nebss_bs(i,:)
    deltabss(i)=deltabss_bs(i)
  enddo

  nebss(nbss+1,1)=1
  nebss(nbss+1,2)=nl
  nebss(nleads,1)=n1-nr+1
  nebss(nleads,2)=n1

  mattype=2

  nnz=sgeneral%matSparse%nnz

  do i1=1,nbss
    call findnnz_subblock(Sgeneral,nebss(i1,1),nebss(i1,1),nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1,nnz)
!    write(12347,*)"nnz=",i1,nnz
    call AllocateMatrixGeneral(nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1,nnz,mattype,gammamp(i1),sMyName,io)
    call AllocateMatrixGeneral(nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1,nnz,mattype,sigmamp(i1),sMyName,io)

    call set_values_gammamp(Sgeneral,nebss(i1,1),nebss(i1,1),gammamp(i1),deltabss(i1))
    call set_values_sigmamp(gammamp(i1),sigmamp(i1))

  enddo

end subroutine set_gammamp


subroutine set_gammamp_lr(n1,nl,nr,nleads,nbss,nebss,gammamp,sigmamp,sigmal,sigmar)

  use negfmod, only: nebss_bs,deltabss_bs

  character(len=*), parameter :: sMyName="set_gammamp_lr"
  type(matrixTypeGeneral) :: gammamp(nleads)
  type(matrixTypeGeneral) :: sigmamp(nleads)
  integer, intent(in) :: n1,nl,nr,nleads,nbss
  integer nebss(nleads,2)
  complex(kdp),intent(in) :: sigmal(nl,nl),sigmar(nr,nr)
  type(matrixTypeGeneral) :: msigmalr(2)

  type(ioType) :: io
  integer ispin,i,j,i1,i2,i3,i4,ind,ii,jj,nnz,mattype,imin,imax,jmin,jmax

  io%iout=778
  io%isDebug= .false.
  io%memCount=0.0_kdp

 
  call AllocateMatrixGeneral(nl,nl,nl*nl,0,msigmalr(1),sMyName,io)
  msigmalr(1)%matdense%a=sigmal
  call AllocateMatrixGeneral(nr,nr,nr*nr,0,msigmalr(2),sMyName,io)
  msigmalr(2)%matdense%a=sigmar

  mattype=2

  do i1=nleads-1,nleads
!    write(12347,*)"nno=",nebss(i1,1),nebss(i1,1),nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1
    call findnnz_subblock_densesigmalr(msigmalr(i1-nleads+2)%matdense,nebss(i1,1),nebss(i1,1),nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1,nnz)
!    write(12347,*)"nnoa=",nebss(i1,1),nebss(i1,1),nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1,nnz
!    write(12347,*)"nnz=",i1,nnz

    call AllocateMatrixGeneral(nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1,nnz,mattype,gammamp(i1),sMyName,io)
    call AllocateMatrixGeneral(nebss(i1,2)-nebss(i1,1)+1,nebss(i1,2)-nebss(i1,1)+1,nnz,mattype,sigmamp(i1),sMyName,io)
!
    call set_values_sigmamplr(msigmalr(i1-nleads+2),nebss(i1,1),nebss(i1,1),gammamp(i1),sigmamp(i1))

  enddo
  call DestroyMatrixGeneral(msigmalr(1),sMyName,io)
  call DestroyMatrixGeneral(msigmalr(2),sMyName,io)

end subroutine set_gammamp_lr


subroutine set_values_sigmamp(gammampi,sigmampi)

 type(matrixTypeGeneral) :: gammampi,sigmampi
 integer ii,jj

 sigmampi%matSparse%q(:)=gammampi%matSparse%q(:)
 sigmampi%matSparse%j(:)=gammampi%matSparse%j(:)
 sigmampi%matSparse%b(:)=(0_kdp,-0.5_kdp) * gammampi%matSparse%b(:)

! do ii=1,sigmampi%matSparse%iRows
!   write(12347,*)"i0=",ii,sigmampi%matSparse%q(ii),sigmampi%matSparse%q(ii+1)-1
!   do jj=sigmampi%matSparse%q(ii),sigmampi%matSparse%q(ii+1)-1
!     write(12347,*)"sigmampi=",ii,sigmampi%matSparse%j(jj),sigmampi%matSparse%b(jj)
!   enddo
! enddo

end subroutine set_values_sigmamp



subroutine set_values_gammamp(Sgeneral,iVert,iHorz,gammampi,delta)

 implicit none
 integer iVert,iHorz
 integer ii,ind,jj,imin,imax,jmin,jmax,ind2
 type(matrixTypeGeneral) :: Sgeneral,gammampi
 double precision delta

 imin=iVert
 imax=iVert+gammampi%matSparse%iRows-1
 jmin=iHorz
 jmax=iHorz+gammampi%matSparse%iCols-1
! write(12347,*)"jmax=",imin,imax,jmin,jmax

 gammampi%matSparse%q(1)=1
 ind2=1
 do ii=imin,imax
   gammampi%matSparse%q(ii-imin+2)=gammampi%matSparse%q(ii-imin+1)
   do ind=Sgeneral%matSparse%q(ii),Sgeneral%matSparse%q(ii+1)-1
     if((Sgeneral%matSparse%j(ind)<=jmax).and.(Sgeneral%matSparse%j(ind)>=jmin))then
       gammampi%matSparse%q(ii-imin+2)=gammampi%matSparse%q(ii-imin+2)+1
       gammampi%matSparse%b(ind2)=2.0_kdp * delta * Sgeneral%matSparse%b(ind)
       gammampi%matSparse%j(ind2)=Sgeneral%matSparse%j(ind)-ihorz+1
! maybe we can add a relation ind2 to ind
! another option is to hermitianize gammap
       ind2=ind2+1
     endif
   enddo
 enddo

! do ii=1,gammampi%matSparse%iRows
!   write(12347,*)"i0=",ii,gammampi%matSparse%q(ii),gammampi%matSparse%q(ii+1)-1
!   do jj=gammampi%matSparse%q(ii),gammampi%matSparse%q(ii+1)-1
!     write(12347,*)"gammampi=",ii,gammampi%matSparse%j(jj),gammampi%matSparse%b(jj)
!   enddo
! enddo


end subroutine set_values_gammamp

subroutine set_values_sigmamplr(msigmalr,iVert,iHorz,gammampi,sigmampi)

 implicit none
 integer iVert,iHorz
 integer ii,ind,jj,imin,imax,jmin,jmax,ind2
 type(matrixTypeGeneral) :: Sgeneral,gammampi,msigmalr,sigmampi
 complex(kdp) :: zi

 zi=(0.0_kdp,1.0_kdp)
 imin=iVert
 imax=iVert+gammampi%matSparse%iRows-1
 jmin=iHorz
 jmax=iHorz+gammampi%matSparse%iCols-1
! write(12347,*)"jmax=",imin,imax,jmin,jmax


 sigmampi%matSparse%q(1)=1
 ind=1
 do ii=1,sigmampi%matSparse%iRows
   sigmampi%matSparse%q(ii+1)=sigmampi%matSparse%q(ii)
   do jj=1,sigmampi%matSparse%iCols
     if(msigmalr%matdense%a(ii,jj).ne.0.0_kdp)then
       sigmampi%matSparse%q(ii+1)=sigmampi%matSparse%q(ii+1)+1
       sigmampi%matSparse%b(ind)=msigmalr%matdense%a(ii,jj)
       gammampi%matSparse%b(ind)=zi * (msigmalr%matdense%a(ii,jj)-DCONJG(msigmalr%matdense%a(jj,ii)))
       sigmampi%matSparse%j(ind)=jj
       ind=ind+1
     endif
   enddo
 enddo

! do ii=1,sigmampi%matSparse%iRows
!   do jj=sigmampi%matSparse%q(ii),sigmampi%matSparse%q(ii+1)-1
!     write(12347,*)"sigmalr=",ii,jj,sigmampi%matSparse%j(jj),sigmampi%matSparse%b(jj)
!   enddo
! enddo


 gammampi%matSparse%q(:)=sigmampi%matSparse%q(:)
 gammampi%matSparse%j(:)=sigmampi%matSparse%j(:)

end subroutine set_values_sigmamplr





subroutine findnnz_subblock(Sgeneral,iVert,iHorz,icols,irows,nnz)

 implicit none
 integer nnz,irows,icols,iVert,iHorz
 integer ii,ind,jj,imin,imax,jmin,jmax
 type(matrixTypeGeneral) :: Sgeneral

 imin=iVert
 imax=iVert+irows-1
 jmin=iHorz
 jmax=iHorz+icols-1

 nnz=0
 do ii=imin,imax
   do ind=Sgeneral%matSparse%q(ii),Sgeneral%matSparse%q(ii+1)-1
     if((Sgeneral%matSparse%j(ind)<=jmax).and.(Sgeneral%matSparse%j(ind)>=jmin))then
       nnz=nnz+1
     endif
   enddo
 enddo

end subroutine findnnz_subblock


subroutine findnnz_subblock_densesigmalr(sigma,iVert,iHorz,icols,irows,nnz)

 implicit none
 integer nnz,irows,icols,iVert,iHorz
 integer ii,ind,jj,imin,imax,jmin,jmax
 type(matrixType) :: sigma

 imin=iVert
 imax=iVert+irows-1
 jmin=iHorz
 jmax=iHorz+icols-1

 nnz=0
 do ii=1,irows
   do jj=1,icols
     if(sigma%a(ii,jj).ne.0.0_kdp)then
       nnz=nnz+1
     endif
   enddo
 enddo

end subroutine findnnz_subblock_densesigmalr





end module mONBoundStates

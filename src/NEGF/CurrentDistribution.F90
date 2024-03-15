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
!                   CURRENTDISTRIBUTION  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

module mCurrentDistribution
 
use mConstants
implicit none

private

public :: CurrentDistributionGeneral
public :: CurrentDistributionMatrix_RhoK

real(kdp), parameter :: rydberg_to_ev=13.6056981_kdp
contains

subroutine CurrentDistributionGeneral(gfgeneral,n1,nl,nr,sigmal,sigmar,hgeneral,sgeneral,NspinBlocks,NspinComplexMatrix,weightc,ispin,ene,ef,v,T,ik,wk)

  use mMPI_NEGF
  use negfmod, only : deauto, curr_fl_L,curr_fr_L, curr_fl_R,curr_fr_R,em_nflux,em_nbfluxStart
  use mTypes
  use mMatrixUtil
  use mONInterface

  implicit none
      
  type(matrixTypeGeneral), intent(in) :: gfgeneral
  integer, intent(in) :: n1,nl,nr,ispin,NspinComplexMatrix,NspinBlocks,ik
  type(matrixTypeGeneral), intent(in) :: hgeneral(NspinComplexMatrix)
  type(matrixTypeGeneral), intent(in) :: sgeneral
  DOUBLE COMPLEX,intent(in) :: sigmal(nl,nl),weightc
  DOUBLE COMPLEX,intent(in) :: sigmar(nr,nr)
  double precision, intent(in) :: ef,v,T,wk
  double complex, intent(in) :: ene
  DOUBLE COMPLEX, allocatable, save :: GF(:,:), sigmal2(:,:),sigmar2(:,:)

  double precision fl,fr,efl,efr
  double precision flr_L,flr_R,fl_L,fr_L,fl_R,fr_R


  integer ii,jj,ind
  integer nbl2_min, nbl2_max,nbr2_min,nbr2_max
  integer, allocatable ::  nb(:,:)

  allocate(nb(4,em_nflux))

  do ii=1,em_nflux
    nb(1,ii)=1
    nb(2,ii)=em_nbfluxStart(ii)-1
    nb(3,ii)=em_nbfluxStart(ii)
    if(NspinComplexMatrix.le.2)then
      nb(4,ii)=n1
    else
      nb(4,ii)=n1/2
    endif
  enddo

  efl=ef+V/2.D0
  efr=ef-V/2.D0
  call fermi_distribution_real(dreal(ene)-efl,T,fl)
  call fermi_distribution_real(dreal(ene)-efr,T,fr)

  if(abs(curr_fl_L) > 1.000000001D0)then
    fl_L=fl
    if(curr_fl_L<0.D0)fl_L=-fl_L
  else
    fl_L=curr_fl_L
  endif
  if(abs(curr_fr_L) > 1.000000001D0)then
    fr_L=fr
    if(curr_fr_L<0.0D0)fr_L=-fr_L
  else
    fr_L=curr_fr_L 
  endif
  flr_L = (fl_L+fr_L)

  if(abs(curr_fl_R) > 1.000000001D0)then
    fl_R=fl
    if(curr_fl_R<0.D0)fl_R=-fl_R
  else
    fl_R=curr_fl_R
  endif
  if(abs(curr_fr_R) > 1.000000001D0)then
    fr_R=fr
    if(curr_fr_R<0.0D0)fr_R=-fr_R
  else
    fr_R=curr_fr_R 
  endif
  flr_R = (fl_R+fr_R)
!xxx: to go to complex energies fl and fr should be used as complex

  write(12347,*)"flr_L,flr_R=",dreal(ene),flr_L,flr_R

   if(gfgeneral%mattype.eq.0)then
      if(NspinBlocks<=2)then
        if(ispin==1)then
          allocate(GF(2*n1,2*n1))
          gf=0.0D0
          gf(1:n1,1:n1)=gfgeneral%matdense%a
          allocate(sigmal2(2*nl,2*nl))
          sigmal2=0.0D0
          sigmal2(1:nl,1:nl)=sigmal
          allocate(sigmar2(2*nr,2*nr))
          sigmar2=0.0D0
          sigmar2(1:nr,1:nr)=sigmar
          return
        else
          gf(n1+1:2*n1,n1+1:2*n1)=gfgeneral%matdense%a

          sigmal2(nl+1:2*nl,nl+1:2*nl)=sigmal
          sigmar2(nr+1:2*nr,nr+1:2*nr)=sigmar
        endif
        call  CurrentDistributionDense(gf,2*n1,2*nl,2*nr,sigmal2,sigmar2,hgeneral,sgeneral,NspinBlocks,NspinComplexMatrix,weightc,flr_L,flr_R,ene,ef,v,nb,em_nflux,ik,wk)
        deallocate(gf,sigmal2,sigmar2)
      else
!        write(12347,*)"calling CurrentDistributionDenseNC"
        call CurrentDistributionDenseNC(gfgeneral%matdense%a,n1,nl,nr,sigmal,sigmar,hgeneral,sgeneral,NspinBlocks,NspinComplexMatrix,weightc,flr_L,flr_R,ene,ef,v,nb,em_nflux,ik,wk)
      endif
!    elseif(gfmattype.eq.2)then
!      call updaterhosparse(rhogeneral(ispin),ematgeneral(ispin),emforces,gf,nl,nr,gfmattype,weightc,cl,cr,weightrho,ene,set_rho_boundary)
    endif
    deallocate (nb)

end subroutine CurrentDistributionGeneral


subroutine CurrentDistributionDense(GF,n1,nl,nr,sigmal,sigmar,hgeneral,sgeneral,NspinBlocks,NspinComplexMatrix,weightc,fl,fr,ene,ef,v,nb,nflux,ik,wk)

  use mMPI_NEGF
  use mTypes
  use mMatrixUtil
  use mONInterface

  implicit none
      
  integer, intent(in) :: n1,nl,nr,NspinComplexMatrix,NspinBlocks,nflux,nb(4,nflux),ik
  type(matrixTypeGeneral), intent(in) :: hgeneral(NspinComplexMatrix)
  type(matrixTypeGeneral), intent(in) :: sgeneral
  DOUBLE COMPLEX,intent(in) :: sigmal(nl,nl),weightc
  DOUBLE COMPLEX,intent(in) :: sigmar(nr,nr)
  DOUBLE COMPLEX, intent(in) :: GF(n1,n1)
  double precision, intent(in) :: fl,fr,ef,v,wk
  double complex, intent(in) :: ene

  DOUBLE COMPLEX, ALLOCATABLE :: gammal(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gammar(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: GF_iter2l(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gf1(:,:),gf2(:,:)
  double complex, allocatable :: b(:)
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  DOUBLE COMPLEX, ALLOCATABLE :: work(:),GF_dag(:,:),GF_less(:,:),Gamma_tot(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: mat1(:,:),mat2(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: Hdense(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: Sdense(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: jmat(:,:)
  double complex j12(4)
  double complex j12b(4)
  DOUBLE COMPLEX, ALLOCATABLE :: rhox(:,:),rhoy(:,:),rhoz(:,:),rho0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: hx(:,:),hy(:,:),hz(:,:),h0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: jx(:,:),jy(:,:),jz(:,:),j0(:,:)
  double precision efl,efr


  integer ii,jj,ind,il
  integer nlh,nrh,n1h

  ALLOCATE(gammal(nl,nl))
  ALLOCATE(gammar(nr,nr))
  gammal=zi* (sigmal-DCONJG(TRANSPOSE(sigmal)))*fl
  gammar=zi* (sigmar-DCONJG(TRANSPOSE(sigmar)))*fr

  nlh=nl/2
  nrh=nr/2
  n1h=n1/2

!  write(12347,*)"n1=",n1,NspinComplexMatrix,NspinBlocks,nl,nr

  ALLOCATE(Gamma_tot(N1,N1),GF_less(N1,N1))
  Gamma_tot=0.0d0

!  write(12347,*)"nl,nlh",nl,nlh


  Gamma_tot(1:nlh,1:nlh)=gammal(1:nlh,1:nlh)
  Gamma_tot(1:nlh,n1h+1:n1h+nlh)=gammal(1:nlh,nlh+1:nl)
  Gamma_tot(n1h+1:n1h+nlh,1:nlh)=gammal(nlh+1:nl,1:nlh)
  Gamma_tot(n1h+1:n1h+nlh,n1h+1:n1h+nlh)=gammal(nlh+1:nl,nlh+1:nl)

  Gamma_tot(n1h-nrh+1:n1h,n1h-nrh+1:n1h)=gammar(1:nrh,1:nrh)
  Gamma_tot(n1h-nrh+1:n1h,n1-nrh+1:n1)=gammar(1:nrh,nrh+1:nr)
  Gamma_tot(n1-nrh+1:n1,n1h-nrh+1:n1h)=gammar(nrh+1:nr,1:nrh)
  Gamma_tot(n1-nrh+1:n1,n1-nrh+1:n1)=gammar(nrh+1:nr,nrh+1:nr)

  ALLOCATE(gf_dag(n1,n1))
  gf_dag=DCONJG(TRANSPOSE(GF))
  CALL ZGEMM('N','N',N1,n1,n1,(1.D0,0.D0),Gamma_tot,N1, gf_dag,n1,(0.D0,0.D0),GF_less,N1)
  CALL ZGEMM('N','N',N1,n1,n1,(1.D0,0.D0),gf,N1, GF_less,n1,(0.D0,0.D0),Gamma_tot,N1)
  GF_less= zi * Gamma_tot !gf_less is actually the density matrix; there would be a 1/2pi, but for the current there woudl be a multiplication by 2pi later, which would cancel this, so we don't multiply it in the first place
!  GF_less=-zi * matmul(gf,MATMUL(Gamma_tot,GF_dag))
  deallocate(gf_dag)

  ALLOCATE(Hdense(n1,n1))

  Hdense=0.0d0
  DO ii=1,n1h
    DO ind=hgeneral(1)%matSparse%q(ii),hgeneral(1)%matSparse%q(ii+1)-1
      Hdense(ii,hgeneral(1)%matSparse%j(ind))=hgeneral(1)%matSparse%b(ind)
      Hdense(n1h+ii,n1h+hgeneral(2)%matSparse%j(ind))=hgeneral(2)%matSparse%b(ind)
      if(NspinComplexMatrix>2)Hdense(ii,n1h+hgeneral(3)%matSparse%j(ind))=hgeneral(3)%matSparse%b(ind)
    ENDDO
  ENDDO
  if(NspinComplexMatrix>2) Hdense(n1h+1:n1,1:n1h)=dconjg(transpose(Hdense(1:n1h,n1h+1:n1)))

  ALLOCATE(sdense(n1,n1))

  Sdense=0.0d0
  DO ii=1,n1h
    DO ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
      sdense(ii,sgeneral%matSparse%j(ind))=sgeneral%matSparse%b(ind)
      sdense(n1h+ii,n1h+sgeneral%matSparse%j(ind))=sgeneral%matSparse%b(ind)
    ENDDO
  ENDDO

  
  allocate(rhoz(n1h,n1h),rho0(n1h,n1h))
  allocate(hz(n1h,n1h),h0(n1h,n1h))
  allocate(jz(n1h,n1h),j0(n1h,n1h))
  if(NspinComplexMatrix>2)then
    allocate(rhox(n1h,n1h),rhoy(n1h,n1h))
    allocate(hx(n1h,n1h),hy(n1h,n1h))
    allocate(jx(n1h,n1h),jy(n1h,n1h))
  endif

  if(NspinComplexMatrix>2)then
    hx=0.5D0*(Hdense(1:n1h,n1h+1:n1)+Hdense(n1h+1:n1,1:n1h))
    hy=(0.0D0,0.5D0)*(Hdense(1:n1h,n1h+1:n1)-Hdense(n1h+1:n1,1:n1h))
    hz=0.5D0*(Hdense(1:n1h,1:n1h)-Hdense(n1h+1:n1,n1h+1:n1))
    h0=0.5D0*(Hdense(1:n1h,1:n1h)+Hdense(n1h+1:n1,n1h+1:n1))

    rhox=0.5D0*(GF_less(1:n1h,n1h+1:n1)+GF_less(n1h+1:n1,1:n1h))
    rhoy=(0.0D0,0.5D0)*(GF_less(1:n1h,n1h+1:n1)-GF_less(n1h+1:n1,1:n1h))
    rhoz=0.5D0*(GF_less(1:n1h,1:n1h)-   GF_less(n1h+1:n1,n1h+1:n1))
    rho0=0.5D0*(GF_less(1:n1h,1:n1h)+   GF_less(n1h+1:n1,n1h+1:n1))
  else
    hz=0.5D0*(Hdense(1:n1h,1:n1h)-Hdense(n1h+1:n1,n1h+1:n1))
    h0=0.5D0*(Hdense(1:n1h,1:n1h)+Hdense(n1h+1:n1,n1h+1:n1))

    rhoz=0.5D0*(GF_less(1:n1h,1:n1h)-   GF_less(n1h+1:n1,n1h+1:n1))
    rho0=0.5D0*(GF_less(1:n1h,1:n1h)+   GF_less(n1h+1:n1,n1h+1:n1))
  endif

  do ii=1,n1h
    do jj=1,n1h

      if(NspinComplexMatrix>2)then
        jx(ii,jj)=2.0D0* (hx(ii,jj)*rho0(jj,ii)+(h0(ii,jj)-ene*sdense(ii,jj))*rhox(jj,ii)&
                       & -hx(jj,ii)*rho0(ii,jj)-(h0(jj,ii)-ene*sdense(jj,ii))*rhox(ii,jj))
        jy(ii,jj)=2.0D0* (hy(ii,jj)*rho0(jj,ii)-(h0(jj,ii)-ene*sdense(jj,ii))*rhoy(ii,jj)&
                       & -hy(jj,ii)*rho0(ii,jj)+(h0(ii,jj)-ene*sdense(ii,jj))*rhoy(jj,ii))
        j0(ii,jj)=2.0D0* ((h0(ii,jj)-ene*sdense(ii,jj))*rho0(jj,ii)-(h0(jj,ii)-ene*sdense(jj,ii))*rho0(ii,jj)&
               &+((hx(ii,jj))*rhox(jj,ii)-(hx(jj,ii))*rhox(ii,jj))&
               &+((hy(ii,jj))*rhoy(jj,ii)-(hy(jj,ii))*rhoy(ii,jj))&
               &+((hz(ii,jj))*rhoz(jj,ii)-(hz(jj,ii))*rhoz(ii,jj)))
      else
        j0(ii,jj)=2.0D0* ((h0(ii,jj)-ene*sdense(ii,jj))*rho0(jj,ii)-(h0(jj,ii)-ene*sdense(jj,ii))*rho0(ii,jj)&
               &+((hz(ii,jj))*rhoz(jj,ii)-(hz(jj,ii))*rhoz(ii,jj)))
      endif
      jz(ii,jj)=2.0D0* (hz(ii,jj)*rho0(jj,ii)-(h0(jj,ii)-ene*sdense(jj,ii))*rhoz(ii,jj)&
                     & -hz(jj,ii)*rho0(ii,jj)+(h0(ii,jj)-ene*sdense(ii,jj))*rhoz(jj,ii))


    enddo
  enddo

  do il=1,nflux
!    write(12347,*)"nboundary=",nb(1,il),nb(2,il),nb(3,il),nb(4,il)
    j12b=0.0D0
    do ii=nb(1,il),nb(2,il)
      do jj=nb(3,il),nb(4,il)
        if(NspinComplexMatrix>2)then
          j12b(1)=j12b(1)+jx(ii,jj)
          j12b(2)=j12b(2)+jy(ii,jj)
        endif
        j12b(3)=j12b(3)+jz(ii,jj)
        j12b(4)=j12b(4)+j0(ii,jj)
      enddo
    enddo
    if(NspinComplexMatrix>2)then
      write(12346,*)"j12x=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(1)),il,ik,v*rydberg_to_ev
      write(12346,*)"j12y=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(2)),il,ik,v*rydberg_to_ev
    endif
    write(12346,*)"j12z=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(3)),il,ik,v*rydberg_to_ev
    write(12346,*)"j120=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(4)),il,ik,v*rydberg_to_ev
  enddo

  deallocate(Hdense,Sdense)


  DEALLOCATE(Gamma_tot,GF_less)

  deallocate(rhoz,rho0)
  deallocate(hz,h0)
  deallocate(jz,j0)
  if(NspinComplexMatrix>2)then
    deallocate(rhox,rhoy)
    deallocate(hx,hy)
    deallocate(jx,jy)
  endif

end subroutine CurrentDistributionDense

subroutine CurrentDistributionDenseNC(GF,n1,nl,nr,sigmal,sigmar,hgeneral,sgeneral,NspinBlocks,NspinComplexMatrix,weightc,fl,fr,ene,ef,v,nb,nflux,ik,wk)

  use mMPI_NEGF
  use mTypes
  use mMatrixUtil
  use mONInterface

  implicit none
      
  integer, intent(in) :: n1,nl,nr,NspinComplexMatrix,NspinBlocks,nflux,nb(4,nflux),ik
  type(matrixTypeGeneral), intent(in) :: hgeneral(NspinComplexMatrix)
  type(matrixTypeGeneral), intent(in) :: sgeneral
  DOUBLE COMPLEX,intent(in) :: sigmal(nl,nl),weightc
  DOUBLE COMPLEX,intent(in) :: sigmar(nr,nr)
  DOUBLE COMPLEX, intent(in) :: GF(n1,n1)
  double precision, intent(in) :: fl,fr,ef,v,wk
  double complex, intent(in) :: ene

  DOUBLE COMPLEX, ALLOCATABLE :: gammal(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gammar(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: GF_iter2l(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gf1(:,:),gf2(:,:)
  double complex, allocatable :: b(:)
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  DOUBLE COMPLEX, ALLOCATABLE :: work(:),GF_dag(:,:),GF_less(:,:),Gamma_tot(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: mat1(:,:),mat2(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: Hdense(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: Sdense(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: jmat(:,:)
  double complex j12(4)
  double complex j12b(4)
  DOUBLE COMPLEX, ALLOCATABLE :: rhox(:,:),rhoy(:,:),rhoz(:,:),rho0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: hx(:,:),hy(:,:),hz(:,:),h0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: jx(:,:),jy(:,:),jz(:,:),j0(:,:)
  double precision efl,efr


  integer ii,jj,ind,il
  integer nlh,nrh,n1h

  nlh=nl/2
  nrh=nr/2
  n1h=n1/2

!  write(12347,*)"n1=",n1,NspinComplexMatrix,NspinBlocks,nl,nr

  ALLOCATE(gammal(nl,nl))
  ALLOCATE(gammar(nr,nr))
!  write(12347,*)"fl,fr",dreal(ene-ef),fl,fr
  gammal=zi* (sigmal-DCONJG(TRANSPOSE(sigmal)))*fl
  gammar=zi* (sigmar-DCONJG(TRANSPOSE(sigmar)))*fr

  ALLOCATE(Gamma_tot(N1,N1),GF_less(N1,N1),gf_dag(n1,n1))
  Gamma_tot=0.0d0

  Gamma_tot(1:nlh,1:nlh)=gammal(1:nlh,1:nlh)
  Gamma_tot(1:nlh,n1h+1:n1h+nlh)=gammal(1:nlh,nlh+1:nl)
!  write(12347,*)"nl,nlh",nl,nlh
  Gamma_tot(n1h+1:n1h+nlh,1:nlh)=gammal(nlh+1:nl,1:nlh)
  Gamma_tot(n1h+1:n1h+nlh,n1h+1:n1h+nlh)=gammal(nlh+1:nl,nlh+1:nl)

  Gamma_tot(n1h-nrh+1:n1h,n1h-nrh+1:n1h)=gammar(1:nrh,1:nrh)
  Gamma_tot(n1h-nrh+1:n1h,n1-nrh+1:n1)=gammar(1:nrh,nrh+1:nr)
  Gamma_tot(n1-nrh+1:n1,n1h-nrh+1:n1h)=gammar(nrh+1:nr,1:nrh)
  Gamma_tot(n1-nrh+1:n1,n1-nrh+1:n1)=gammar(nrh+1:nr,nrh+1:nr)

  gf_dag=DCONJG(TRANSPOSE(GF))

  CALL ZGEMM('N','N',N1,n1,n1,(1.D0,0.D0),Gamma_tot,N1, gf_dag,n1,(0.D0,0.D0),GF_less,N1)
  CALL ZGEMM('N','N',N1,n1,n1,(1.D0,0.D0),gf,N1, GF_less,n1,(0.D0,0.D0),Gamma_tot,N1)
  GF_less= zi * Gamma_tot !gf_less is actually the density matrix; there would be a 1/2pi, but for the current there woudl be a multiplication by 2pi later, which would cancel this, so we don't multiply it in the first place
!  GF_less=-zi * matmul(gf,MATMUL(Gamma_tot,GF_dag))

  ALLOCATE(Hdense(n1,n1))

  Hdense=0.0d0
  DO ii=1,n1h
    DO ind=hgeneral(1)%matSparse%q(ii),hgeneral(1)%matSparse%q(ii+1)-1
      Hdense(ii,hgeneral(1)%matSparse%j(ind))=hgeneral(1)%matSparse%b(ind)
      Hdense(n1h+ii,n1h+hgeneral(2)%matSparse%j(ind))=hgeneral(2)%matSparse%b(ind)
      Hdense(ii,n1h+hgeneral(3)%matSparse%j(ind))=hgeneral(3)%matSparse%b(ind)
    ENDDO
  ENDDO
  Hdense(n1h+1:n1,1:n1h)=dconjg(transpose(Hdense(1:n1h,n1h+1:n1)))

  ALLOCATE(sdense(n1,n1))

  Sdense=0.0d0
  DO ii=1,n1h
    DO ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
      sdense(ii,sgeneral%matSparse%j(ind))=sgeneral%matSparse%b(ind)
      sdense(n1h+ii,n1h+sgeneral%matSparse%j(ind))=sgeneral%matSparse%b(ind)
    ENDDO
  ENDDO

  allocate(rhox(n1h,n1h),rhoy(n1h,n1h),rhoz(n1h,n1h),rho0(n1h,n1h))
  allocate(hx(n1h,n1h),hy(n1h,n1h),hz(n1h,n1h),h0(n1h,n1h))
  allocate(jx(n1h,n1h),jy(n1h,n1h),jz(n1h,n1h),j0(n1h,n1h))

  hx=0.5D0*(Hdense(1:n1h,n1h+1:n1)+Hdense(n1h+1:n1,1:n1h))
  hy=(0.0D0,0.5D0)*(Hdense(1:n1h,n1h+1:n1)-Hdense(n1h+1:n1,1:n1h))
  hz=0.5D0*(Hdense(1:n1h,1:n1h)-Hdense(n1h+1:n1,n1h+1:n1))
  h0=0.5D0*(Hdense(1:n1h,1:n1h)+Hdense(n1h+1:n1,n1h+1:n1))


  rhox=0.5D0*(GF_less(1:n1h,n1h+1:n1)+GF_less(n1h+1:n1,1:n1h))
  rhoy=(0.0D0,0.5D0)*(GF_less(1:n1h,n1h+1:n1)-GF_less(n1h+1:n1,1:n1h))
  rhoz=0.5D0*(GF_less(1:n1h,1:n1h)-   GF_less(n1h+1:n1,n1h+1:n1))
  rho0=0.5D0*(GF_less(1:n1h,1:n1h)+   GF_less(n1h+1:n1,n1h+1:n1))

  do ii=1,n1h
    do jj=1,n1h

      jx(ii,jj)=2.0D0* (hx(ii,jj)*rho0(jj,ii)+(h0(ii,jj)-ene*sdense(ii,jj))*rhox(jj,ii)&
                     & -hx(jj,ii)*rho0(ii,jj)-(h0(jj,ii)-ene*sdense(jj,ii))*rhox(ii,jj))
      jy(ii,jj)=2.0D0* (hy(ii,jj)*rho0(jj,ii)-(h0(jj,ii)-ene*sdense(jj,ii))*rhoy(ii,jj)&
                     & -hy(jj,ii)*rho0(ii,jj)+(h0(ii,jj)-ene*sdense(ii,jj))*rhoy(jj,ii))
      jz(ii,jj)=2.0D0* (hz(ii,jj)*rho0(jj,ii)-(h0(jj,ii)-ene*sdense(jj,ii))*rhoz(ii,jj)&
                     & -hz(jj,ii)*rho0(ii,jj)+(h0(ii,jj)-ene*sdense(ii,jj))*rhoz(jj,ii))

      j0(ii,jj)=2.0D0* ((h0(ii,jj)-ene*sdense(ii,jj))*rho0(jj,ii)-(h0(jj,ii)-ene*sdense(jj,ii))*rho0(ii,jj)&
             &+((hx(ii,jj))*rhox(jj,ii)-(hx(jj,ii))*rhox(ii,jj))&
             &+((hy(ii,jj))*rhoy(jj,ii)-(hy(jj,ii))*rhoy(ii,jj))&
             &+((hz(ii,jj))*rhoz(jj,ii)-(hz(jj,ii))*rhoz(ii,jj)))

    enddo
  enddo

!!  call writemat8(ene,Hdense,n1,n1,0.0D0,"h")
!!  call writemat8(ene,GF_less,n1,n1,0.0D0,"gf_less")

  do il=1,nflux
!    write(12347,*)"nboundary=",nb(1,il),nb(2,il),nb(3,il),nb(4,il)
    j12b=0.0D0
    do ii=nb(1,il),nb(2,il)
      do jj=nb(3,il),nb(4,il)
        j12b(1)=j12b(1)+jx(ii,jj)
        j12b(2)=j12b(2)+jy(ii,jj)
        j12b(3)=j12b(3)+jz(ii,jj)
        j12b(4)=j12b(4)+j0(ii,jj)
      enddo
    enddo
!   if(dimag(j12b(:))>1.0e-9)write(12347,*)"warning: the imaginary part of the current is larger than 1.0e-9"
    write(12346,*)"j12x=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(1)),wk,il,ik,v*rydberg_to_ev
    write(12346,*)"j12y=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(2)),wk,il,ik,v*rydberg_to_ev
    write(12346,*)"j12z=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(3)),wk,il,ik,v*rydberg_to_ev
    write(12346,*)"j120=",rydberg_to_ev*dreal(ene-ef),dreal(j12b(4)),wk,il,ik,v*rydberg_to_ev
  enddo

  deallocate(Hdense,Sdense)


  DEALLOCATE(Gamma_tot,GF_less,gf_dag)

  deallocate(rhox,rhoy,rhoz,rho0)
  deallocate(hx,hy,hz,h0)
  deallocate(jx,jy,jz,j0)


  end subroutine CurrentDistributionDenseNC



subroutine CurrentDistributionMatrix_RhoK(rhogeneral,ematgeneral, hgeneral,sgeneral,NspinComplexMatrix,n1h, v,wk,ik)

  use mMPI_NEGF
  use mTypes
  use mMatrixUtil
  use negfmod, only : em_nflux,em_nbfluxStart

  integer, intent(in)                 :: n1h,NspinComplexMatrix,ik
  type(matrixTypeGeneral), intent(in) :: hgeneral(NspinComplexMatrix),sgeneral
  type(matrixTypeGeneral), intent(in) :: rhogeneral(NspinComplexMatrix)
  type(matrixTypeGeneral), intent(in) :: ematgeneral(NspinComplexMatrix)
  real(kdp), intent(in)               :: v,wk

  complex(kdp), allocatable :: HDense(:,:)
  complex(kdp), allocatable :: OmegaDense(:,:)
  complex(kdp), allocatable :: RhoDense(:,:)
  complex(kdp), allocatable :: SDense(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: rhox(:,:),rhoy(:,:),rhoz(:,:),rho0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: omegax(:,:),omegay(:,:),omegaz(:,:),omega0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: hx(:,:),hy(:,:),hz(:,:),h0(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: jx(:,:),jy(:,:),jz(:,:),j0(:,:)
  integer                     :: n1,i,j,ind,il
  real(kdp) PI
  integer, allocatable        ::  nb(:,:)
  complex(kdp) :: j12b(4)

  PI=2.0_kdp * acos(0.0_kdp)
  n1=2*n1h

  allocate(nb(4,em_nflux))
  do i=1,em_nflux
    nb(1,i)=1
    nb(2,i)=em_nbfluxStart(i)-1
    nb(3,i)=em_nbfluxStart(i)
    nb(4,i)=n1/2
  enddo

  allocate(Hdense(n1,n1))
  Hdense=0.0d0
  do i=1,n1h
    do ind=hgeneral(1)%matSparse%q(i),hgeneral(1)%matSparse%q(i+1)-1
      Hdense(i,hgeneral(1)%matSparse%j(ind))=hgeneral(1)%matSparse%b(ind)
!xxx: make it work also for nonspinpolarized calculations
      Hdense(n1h+i,n1h+hgeneral(2)%matSparse%j(ind))=hgeneral(2)%matSparse%b(ind)
      if(NspinComplexMatrix>2) Hdense(i,n1h+hgeneral(3)%matSparse%j(ind))=hgeneral(3)%matSparse%b(ind)
    enddo
  enddo
  if(NspinComplexMatrix>2) Hdense(n1h+1:n1,1:n1h)=dconjg(transpose(Hdense(1:n1h,n1h+1:n1)))



  ALLOCATE(sdense(n1,n1))
  Sdense=0.0d0
  DO i=1,n1h
    DO ind=sgeneral%matSparse%q(i),sgeneral%matSparse%q(i+1)-1
      sdense(i,sgeneral%matSparse%j(ind))=sgeneral%matSparse%b(ind)
      sdense(n1h+i,n1h+sgeneral%matSparse%j(ind))=sgeneral%matSparse%b(ind)
    ENDDO
  ENDDO

  allocate(Rhodense(n1,n1))
  Rhodense=0.0d0
  do i=1,n1h
    do ind=rhogeneral(1)%matSparse%q(i),rhogeneral(1)%matSparse%q(i+1)-1
      Rhodense(i,rhogeneral(1)%matSparse%j(ind))=rhogeneral(1)%matSparse%b(ind)
      Rhodense(n1h+i,n1h+rhogeneral(2)%matSparse%j(ind))=rhogeneral(2)%matSparse%b(ind)
      if(NspinComplexMatrix>2) Rhodense(i,n1h+rhogeneral(3)%matSparse%j(ind))=rhogeneral(3)%matSparse%b(ind)
    enddo
  enddo
  if(NspinComplexMatrix>2) Rhodense(n1h+1:n1,1:n1h)=dconjg(transpose(Rhodense(1:n1h,n1h+1:n1)))

  allocate(Omegadense(n1,n1))
  Omegadense=0.0d0
  do i=1,n1h
    do ind=ematgeneral(1)%matSparse%q(i),ematgeneral(1)%matSparse%q(i+1)-1
      Omegadense(i,ematgeneral(1)%matSparse%j(ind))=ematgeneral(1)%matSparse%b(ind)
      Omegadense(n1h+i,n1h+ematgeneral(2)%matSparse%j(ind))=ematgeneral(2)%matSparse%b(ind)
      if(NspinComplexMatrix>2) Omegadense(i,n1h+ematgeneral(3)%matSparse%j(ind))=ematgeneral(3)%matSparse%b(ind)
    enddo
  enddo
  if(NspinComplexMatrix>2) Omegadense(n1h+1:n1,1:n1h)=dconjg(transpose(Omegadense(1:n1h,n1h+1:n1)))
  
  allocate(rhoz(n1h,n1h),rho0(n1h,n1h))
  allocate(omegaz(n1h,n1h),omega0(n1h,n1h))
  allocate(hz(n1h,n1h),h0(n1h,n1h))
  allocate(jz(n1h,n1h),j0(n1h,n1h))

  if(NspinComplexMatrix>2) then
    allocate(rhox(n1h,n1h),rhoy(n1h,n1h))
    allocate(omegax(n1h,n1h),omegay(n1h,n1h))
    allocate(hx(n1h,n1h),hy(n1h,n1h))
    allocate(jx(n1h,n1h),jy(n1h,n1h))
  endif

  hz=0.5D0*(Hdense(1:n1h,1:n1h)-Hdense(n1h+1:n1,n1h+1:n1))
  h0=0.5D0*(Hdense(1:n1h,1:n1h)+Hdense(n1h+1:n1,n1h+1:n1))

  rhoz=0.5D0*(Rhodense(1:n1h,1:n1h)- Rhodense(n1h+1:n1,n1h+1:n1))
  rho0=0.5D0*(Rhodense(1:n1h,1:n1h)+ Rhodense(n1h+1:n1,n1h+1:n1))

  Omegaz=0.5D0*(Omegadense(1:n1h,1:n1h)- Omegadense(n1h+1:n1,n1h+1:n1))
  Omega0=0.5D0*(Omegadense(1:n1h,1:n1h)+ Omegadense(n1h+1:n1,n1h+1:n1))


  if(NspinComplexMatrix>2)then
    hx=0.5D0*(Hdense(1:n1h,n1h+1:n1)+Hdense(n1h+1:n1,1:n1h))
    hy=(0.0D0,0.5D0)*(Hdense(1:n1h,n1h+1:n1)-Hdense(n1h+1:n1,1:n1h))

    rhox=0.5D0*(Rhodense(1:n1h,n1h+1:n1)+Rhodense(n1h+1:n1,1:n1h))
    rhoy=(0.0D0,0.5D0)*(Rhodense(1:n1h,n1h+1:n1)-Rhodense(n1h+1:n1,1:n1h))

    omegax=0.5D0*(Omegadense(1:n1h,n1h+1:n1)+Omegadense(n1h+1:n1,1:n1h))
    omegay=(0.0D0,0.5D0)*(Omegadense(1:n1h,n1h+1:n1)-Omegadense(n1h+1:n1,1:n1h))
  endif




!  call writemat8(0.0D0,h0,n1h,n1h,0.0D0,"h0")

  do i=1,n1h
    do j=1,n1h
!!!      write(12346,*)"densek:h0=",i,j,dreal(h0(i,j)),dimag(h0(i,j))
!!!      write(12346,*)"densek:hz=",i,j,dreal(hz(i,j)),dimag(hz(i,j))
!!!      write(12346,*)"densek:s0=",i,j,dreal(sdense(i,j)),dimag(sdense(i,j))
!!!      write(12346,*)"densek:rho0=",i,j,dreal(rho0(i,j)),dimag(rho0(i,j))
!!!      write(12346,*)"densek:rhoz=",i,j,dreal(rhoz(i,j)),dimag(rhoz(i,j))
!!!      write(12346,*)"densek:omega0=",i,j,dreal(omega0(i,j)),dimag(omega0(i,j))
!!!      write(12346,*)"densek:omegaz=",i,j,dreal(omegaz(i,j)),dimag(omegaz(i,j))

!!!!      j0(i,j)=(0.0_kdp,4.0_kdp)* PI * (h0(i,j)*rho0(j,i)-sdense(i,j)*Omega0(j,i)+hz(i,j)*rhoz(j,i)-&
!!!!                                      &(h0(j,i)*rho0(i,j)-sdense(j,i)*Omega0(i,j)+hz(j,i)*rhoz(i,j)))
!!!!      j0(i,j)=(0.0_kdp,4.0_kdp)* PI * (h0(i,j)*rho0(j,i)-sdense(i,j)*Omega0(j,i)+hz(i,j)*rhoz(j,i)-&
!!!!                               &dconjg(h0(i,j)*rho0(j,i)-sdense(i,j)*Omega0(j,i)+hz(i,j)*rhoz(j,i)))


      if(NspinComplexMatrix>2)then
        j0(i,j)=-8.0_kdp* PI * dimag(h0(i,j)*dconjg(rho0(i,j))-sdense(i,j)*dconjg(Omega0(i,j))+hx(i,j)*dconjg(rhox(i,j))+hy(i,j)*dconjg(rhoy(i,j))+hz(i,j)*dconjg(rhoz(i,j)))
        jx(i,j)=-8.0_kdp* PI * dimag(hx(i,j)*dconjg(rho0(i,j))+h0(i,j)*dconjg(rhox(i,j))-sdense(i,j)*dconjg(Omegax(i,j)))
        jy(i,j)=-8.0_kdp* PI * dimag(hy(i,j)*dconjg(rho0(i,j))+h0(i,j)*dconjg(rhoy(i,j))-sdense(i,j)*dconjg(Omegay(i,j)))
      else
        j0(i,j)=-8.0_kdp* PI * dimag(h0(i,j)*dconjg(rho0(i,j))-sdense(i,j)*dconjg(Omega0(i,j))+hz(i,j)*dconjg(rhoz(i,j)))

      endif
      jz(i,j)=-8.0_kdp* PI * dimag(hz(i,j)*dconjg(rho0(i,j))+h0(i,j)*dconjg(rhoz(i,j))-sdense(i,j)*dconjg(Omegaz(i,j)))

    enddo
  enddo


  do il=1,em_nflux
!    write(12347,*)"nboundary=",nb(1,il),nb(2,il),nb(3,il),nb(4,il)
    j12b=0.0D0
    do i=nb(1,il),nb(2,il)
      do j=nb(3,il),nb(4,il)
        if(NspinComplexMatrix>2)then
          j12b(1)=j12b(1)+jx(i,j)
          j12b(2)=j12b(2)+jy(i,j)
        endif
        j12b(3)=j12b(3)+jz(i,j)
        j12b(4)=j12b(4)+j0(i,j)
      enddo
    enddo
    if(NspinComplexMatrix==2)then
      write(12346,*)"JLR_k=",   il,dreal(j12b(4)), dreal(j12b(3)), wk, v*rydberg_to_ev,ik
    elseif(NspinComplexMatrix>2)then
      write(12346,*)"JLR_k=",   il,dreal(j12b(4)), dreal(j12b(1)),dreal(j12b(2)),dreal(j12b(3)), wk, v*rydberg_to_ev,ik
    endif
  enddo


  deallocate(Hdense,Sdense,RhoDense,OmegaDense)

  deallocate(rhoz,rho0)
  deallocate(omegaz,omega0)
  deallocate(hz,h0)
  deallocate(jz,j0)

  end subroutine CurrentDistributionMatrix_RhoK

end module mCurrentDistribution

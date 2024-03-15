module mSigmaFourier

  use mConstants
  use mTypes

  implicit none

  public selfenergy_k
  public ExpandSstates
  private 

  contains 

  SUBROUTINE selfenergy_k(SIDE,ndivxy,n,e,H0,H1,S0,S1,Sigma,nchan,delta,DoFourier)

    use negfmod, only: overwritehs,TransmissionMatrix,em_Last_SCF_Step, EM_NSPINBlocks
    use mSigmaMethod1, only :check_error_sigma, GetSelfEnergy
    use mComputeULR, only : PhiS,SetKPhiS

    CHARACTER(LEN=1),intent(in)  :: SIDE
    integer, intent(in)          :: ndivxy(2)
    integer, intent(in)          :: n
    integer, intent(out)         :: nchan
    complex(kdp),intent(inout)   :: h0(n,n),h1(n,n),s0(n,n),s1(n,n)
    complex(kdp),intent(out)     :: sigma(n,n)
    complex(kdp),intent(in)      :: e
    real(kdp), intent(in)        :: delta
    logical, intent(in)          :: DoFourier

    integer i1,i2,j1,j2,ndiv,ns,istart,jstart,iend,jend,irow(2),ikx,iky,ispinold
    real(kdp), allocatable :: k(:,:,:)
    real(kdp), parameter :: pi=3.14159265358979323846264338327950288_kdp
    complex(kdp), parameter :: ii=(0.0_kdp,1.0_kdp)
    complex(kdp),allocatable,save :: h0s(:,:,:,:),h1s(:,:,:,:),s0s(:,:,:,:),s1s(:,:,:,:)
    complex(kdp),allocatable :: sigmas(:,:,:,:),gfs(:,:,:,:)
    complex(kdp),allocatable :: gf2(:,:),sigma2(:,:),k1(:,:),km1(:,:)
    real(kdp) :: dsigma
    complex(kdp) :: eik
    integer nchank
    integer il
    logical, parameter :: checksigma=.false.
    integer nkx,nky

    if(side.eq.'L')then
      il=1
    else
      il=2
    endif

    ndiv=ndivxy(1)*ndivxy(2)
    ns=n/ndiv
    nkx=1*ndivxy(1) ! here instead of 1 it is possible to use a larger integer number, so that more kpoints along x are used
    nky=1*ndivxy(2) ! here instead of 1 it is possible to use a larger integer number, so that more kpoints along y are used
    allocate(k(2,nkx,nky))
    allocate(sigmas(ns,ns,nkx,nky))
    if(checksigma)allocate(gfs(ns,ns,nkx,nky),gf2(n,n),sigma2(n,n),k1(n,n),km1(n,n))

!    k(1)=kpointmod(1) * lucmod(1) * 0.5D0
!    k(2)=kpointmod(1) * lucmod(1) * 0.5D0+2.0_kdp * pi / 2.0_kdp

    do ikx=1,nkx
      do iky=1,nky
        k(1,ikx,iky)=(2.0_kdp * (ikx-1)) * pi /(1.0_kdp * nkx)
        k(2,ikx,iky)=(2.0_kdp * (iky-1)) * pi /(1.0_kdp * nky)
      enddo
    enddo

!    write(12347,*)"n=",n,ndiv,ndivxy(1),ndivxy(2),nkx,nky
!    write(12347,*)"ik,k,l=",ikpmod,kpointmod,lucmod

    if(DoFourier)then
!      write(12347,*)"nck=",n,ndiv,ndivxy(1),ndivxy(2),nkx,nky
!      write(12347,*)"ck,ik,k,l=",ikpmod,kpointmod,lucmod

      if(allocated(h0s)) deallocate(s1s,s0s,h1s,h0s)
      if(.not.allocated(h0s))allocate(s1s(ns,ns,nkx,nky),s0s(ns,ns,nkx,nky),h1s(ns,ns,nkx,nky),h0s(ns,ns,nkx,nky))
      call HSFourier(n,H0,H1,S0,S1,H0S,H1S,S0S,S1S,ndivxy,ns,nkx,nky,k,overwritehs)
    endif

!xxx: check nkxy
    if(TransmissionMatrix.and.em_Last_SCF_Step) call SetKPhiS(PhiS(il),ndivxy,k)

    nchan=0
    do ikx=1,nkx
      do iky=1,nky
        if(TransmissionMatrix.and.em_Last_SCF_Step) then
          PhiS(il)%ik(1)=ikx
          PhiS(il)%ik(2)=iky
        endif
        CALL GetSelfEnergy(side,ns ,e,h0s(:,:,ikx,iky), h1s(:,:,ikx,iky),s0s(:,:,ikx,iky),s1s(:,:,ikx,iky),sigmas(:,:,ikx,iky) ,nchank,delta)
        nchan=nchan+nchank
        if(checksigma)call surfacegf(ns,gfs(:,:,ikx,iky),e,h0s(:,:,ikx,iky),s0s(:,:,ikx,iky),sigmas(:,:,ikx,iky))
      enddo
    enddo

    if(.false.)then
      if(TransmissionMatrix.and.em_Last_SCF_Step) then
!xxx: update nkxy
        call WriteFourierPhiS(PhiS(il),"nphis: ")
      endif
    endif
!inverse Fourier transform
    sigma=0.0_kdp
    do i1=1,ndivxy(1)
      do i2=1,ndivxy(2)
        do j1=1,ndivxy(1)
          do j2=1,ndivxy(2)

            istart=(i1-1) * ns+(i2-1) * ndivxy(1) * ns+1
            iend=i1 * ns+(i2-1) * ndivxy(1) * ns
            jstart=(j1-1) * ns+(j2-1) * ndivxy(1) * ns+1
            jend=j1 * ns+(j2-1) * ndivxy(1) * ns

            do ikx=1,nkx
              do iky=1,nky
                eik=exp(-ii * (k(1,ikx,iky)* (j1-i1)+k(2,ikx,iky)* (j2-i2)))
                sigma(istart:iend,jstart:jend)=sigma(istart:iend,jstart:jend)+ eik * sigmas(:,:,ikx,iky)
              enddo
            enddo

          enddo
        enddo
      enddo
    enddo
    sigma=sigma/(1.0_kdp * nkx*nky)
    deallocate(sigmas)


    if(checksigma)then

!inverse Fourier transform
      gf2=0.0_kdp
      do i1=1,ndivxy(1)
        do i2=1,ndivxy(2)
          do j1=1,ndivxy(1)
            do j2=1,ndivxy(2)

              istart=(i1-1) * ns+(i2-1) * ndivxy(1) * ns+1
              iend=i1 * ns+(i2-1) * ndivxy(1) * ns
              jstart=(j1-1) * ns+(j2-1) * ndivxy(1) * ns+1
              jend=j1 * ns+(j2-1) * ndivxy(1) * ns

              do ikx=1,nkx
                do iky=1,nky
                  eik=exp(-ii * (k(1,ikx,iky)* (j1-i1)+k(2,ikx,iky)* (j2-i2)))
                  gf2(istart:iend,jstart:jend)=gf2(istart:iend,jstart:jend)+ eik * gfs(:,:,ikx,iky)
                enddo
              enddo

            enddo
          enddo
        enddo
      enddo
      gf2=gf2/(1.0_kdp * nkx*nky)

      call SurfacePDOS(gf2,s0,n,e)

      k1=e * s1 - h1
      km1=e * DCONJG(TRANSPOSE(s1)) - DCONJG(TRANSPOSE(h1)) 

      if(side.eq.'L')then
        sigma2=matmul(km1,matmul(gf2,k1))
      else
        sigma2=matmul(k1,matmul(gf2,km1))
      endif
      write(12347,*)"maxdsigmaft=",maxval(abs(sigma-sigma2))/maxval(abs(sigma)),maxval(abs(sigma-sigma2)),maxval(abs(sigma))
!---
      write(12347,*)"dsigmanew="
      call check_error_sigma(dsigma,.true.,e,side,n,  H0,H1,S0,S1,sigma)

      CALL GetSelfEnergy(side,n ,e,H0, H1,S0,S1,sigma2 ,nchank,delta)
      write(12347,*)"dsigmaoldorig="
      call check_error_sigma(dsigma,.true.,e,side,n,  h0,H1,S0,S1,sigma2)

      write(12347,*)"maxdsigmainout=",maxval(abs(sigma-sigma2))/maxval(abs(sigma)),maxval(abs(sigma-sigma2)),maxval(abs(sigma))

      deallocate(km1,k1,gf2,gfs,sigma2)
    endif
    deallocate(k)

  end SUBROUTINE selfenergy_k

  SUBROUTINE HSFourier(n,H0,H1,S0,S1,H0S,H1S,S0S,S1S,ndivxy,ns,nkx,nky,k,overwritehs)

    integer, intent(in) :: n,ndivxy(2),nkx,nky
    complex(kdp),intent(inout) :: h0(n,n),h1(n,n),s0(n,n),s1(n,n)
    complex(kdp),intent(out) :: h0s(ns,ns,nkx,nky),h1s(ns,ns,nkx,nky),s0s(ns,ns,nkx,nky),s1s(ns,ns,nkx,nky)
    real(kdp),intent(in) :: k(2,nkx,nky)
    logical, intent(in) :: overwritehs

    integer i1,i2,j1,j2,ndiv,ns,istart,jstart,iend,jend,irow(2),ikx,iky,i3,j3
    integer istart2,iend2,jstart2,jend2
    integer iRowMin,iRowMax,jRowMin,jRowMax
    real(kdp), parameter :: pi=3.14159265358979323846264338327950288_kdp
    complex(kdp), parameter :: ii=(0.0_kdp,1.0_kdp)
    complex(kdp),allocatable :: h02(:,:),h12(:,:),s02(:,:),s12(:,:)
    integer nx,ny,nxy,i1set,deltai
    complex(kdp) :: eik

    ndiv=ndivxy(1)*ndivxy(2)
    ns=n/ndiv
    allocate(h02(n,n),s02(n,n),h12(n,n),s12(n,n))

!forward Fourier transform
    h0s=0.0_kdp
    h1s=0.0_kdp
    s0s=0.0_kdp
    s1s=0.0_kdp

!default
    iRowMin=1
    iRowMax=ndivxy(1)
!or else these numbers will be set in the input file
!    iRowMin=5
!    iRowMax=5

    jRowMin=1
    jRowMax=ndivxy(2)

    nx=iRowMax-iRowMin+1
    ny=jRowMax-jRowMin+1
    nxy=nx*ny
    do ikx=1,nkx
      do iky=1,nky

        do i1=1,ndivxy(1)
          do i2=1,ndivxy(2)

            do j1=iRowMin,iRowMax
              do j2=jRowMin,jRowMax

                irow(1)=j1
                irow(2)=j2
                eik=exp(ii * (k(1,ikx,iky)* dble(i1-irow(1))+k(2,ikx,iky)* dble(i2-irow(2))))

                istart=(irow(1)-1) * ns+(irow(2)-1) * ndivxy(1) * ns+1
                iend=irow(1) * ns+(irow(2)-1) * ndivxy(1) * ns
                jstart=(i1-1) * ns+(i2-1) * ndivxy(1) * ns+1
                jend=(i1) * ns+(i2-1) * ndivxy(1) * ns

!!!!-->                i1set=5
!!!!-->                irow(1)=i1set
!!!!-->                irow(2)=j2
!!!!-->                deltai=j1-i1set
!!!!-->
!!!!-->                istart=(irow(1)-1) * ns+(irow(2)-1) * ndivxy(1) * ns+1
!!!!-->                iend=irow(1) * ns+(irow(2)-1) * ndivxy(1) * ns
!!!!-->                jstart=(i1-1-deltai) * ns+(i2-1) * ndivxy(1) * ns+1
!!!!-->                jend=(i1-deltai) * ns+(i2-1) * ndivxy(1) * ns
!!!!-->!                write(12347,*)"jstart=",i1,j1,istart,jstart,istart2,jstart2,i1set,deltai
!!!!-->                if(jstart>n)then
!!!!-->                  jstart=jstart-n
!!!!-->                  jend=jend-n
!!!!-->                elseif(jstart<=0)then
!!!!-->                  jstart=jstart+n
!!!!-->                  jend=jend+n
!!!!-->                endif
!!!!-->
!                do i3=1,ns
!                  do j3=1,ns
!                    write(12347,*)"hcompare=",i3,j3,h0(istart+i3-1,jstart+j3-1),h0(istart2+i3-1,jstart2+j3-1),dreal(h0(istart+i3-1,jstart+j3-1)-h0(istart2+i3-1,jstart2+j3-1)),dimag(h0(istart+i3-1,jstart+j3-1)-h0(istart2+i3-1,jstart2+j3-1))
!                  enddo
!                enddo
!                write(12347,*)"h0compare=",istart,jstart,istart2,jstart2,maxval(abs(h0(istart2:iend2,jstart2:jend2)-h0(istart:iend,jstart:jend)))
!                write(12347,*)"h1compare=",istart,jstart,istart2,jstart2,maxval(abs(h1(istart2:iend2,jstart2:jend2)-h1(istart:iend,jstart:jend)))
!                write(12347,*)"s0ompare=",istart,jstart,istart2,jstart2,maxval(abs(s0(istart2:iend2,jstart2:jend2)-s0(istart:iend,jstart:jend)))
!                write(12347,*)"s1compare=",istart,jstart,istart2,jstart2,maxval(abs(s1(istart2:iend2,jstart2:jend2)-s1(istart:iend,jstart:jend)))
!        
                h0s(:,:,ikx,iky)=h0s(:,:,ikx,iky)+h0(istart:iend,jstart:jend) * eik
!                h0s(:,:,ikx,iky)=h0s(:,:,ikx,iky)+0.5D0*(h0(istart:iend,jstart:jend)* eik +DCONJG(TRANSPOSE(h0(istart:iend,jstart:jend)* eik))) 
                s0s(:,:,ikx,iky)=s0s(:,:,ikx,iky)+s0(istart:iend,jstart:jend) * eik
                h1s(:,:,ikx,iky)=h1s(:,:,ikx,iky)+h1(istart:iend,jstart:jend) * eik
                s1s(:,:,ikx,iky)=s1s(:,:,ikx,iky)+s1(istart:iend,jstart:jend) * eik

              enddo
            enddo

          enddo
        enddo

!        do i1=1,ns
!          do i2=1,ns
!            write(12347,*)"h0sb=",i1,i2,dreal(h0s(i1,i2,ikx,iky)),dimag(h0s(i1,i2,ikx,iky)),dreal(h0s(i1,i2,ikx,iky)-DCONJG(h0s(i2,i1,ikx,iky))),dimag(h0s(i1,i2,ikx,iky)-DCONJG(h0s(i2,i1,ikx,iky))),ikx,iky
!          enddo
!        enddo
!

        if(.true.)then
          h0s(:,:,ikx,iky)=0.5D0 * (h0s(:,:,ikx,iky)+DCONJG(TRANSPOSE(h0s(:,:,ikx,iky))))
          s0s(:,:,ikx,iky)=0.5D0 * (s0s(:,:,ikx,iky)+DCONJG(TRANSPOSE(s0s(:,:,ikx,iky))))
        endif

!        do i1=1,ns
!          do i2=1,ns
!            write(12347,*)"h0s=",i1,i2,h0s(i1,i2,ikx,iky),ikx,iky
!          enddo
!        enddo


      enddo
    enddo

    h0s=h0s/(1.0_kdp * nxy)
    h1s=h1s/(1.0_kdp * nxy)
    s0s=s0s/(1.0_kdp * nxy)
    s1s=s1s/(1.0_kdp * nxy)


!inverse Fourier transform
    h02=0.0_kdp
    s02=0.0_kdp
    h12=0.0_kdp
    s12=0.0_kdp
    do i1=1,ndivxy(1)
      do i2=1,ndivxy(2)
        do j1=1,ndivxy(1)
          do j2=1,ndivxy(2)

            istart=(i1-1) * ns+(i2-1) * ndivxy(1) * ns+1
            iend=i1 * ns+(i2-1) * ndivxy(1) * ns
            jstart=(j1-1) * ns+(j2-1) * ndivxy(1) * ns+1
            jend=j1 * ns+(j2-1) * ndivxy(1) * ns

            do ikx=1,nkx
              do iky=1,nky

                eik=exp(-ii * (k(1,ikx,iky)* (j1-i1)+k(2,ikx,iky)* (j2-i2)))
                h02(istart:iend,jstart:jend)=h02(istart:iend,jstart:jend)+ h0s(:,:,ikx,iky) * eik
                s02(istart:iend,jstart:jend)=s02(istart:iend,jstart:jend)+ s0s(:,:,ikx,iky) * eik
                h12(istart:iend,jstart:jend)=h12(istart:iend,jstart:jend)+ h1s(:,:,ikx,iky) * eik
                s12(istart:iend,jstart:jend)=s12(istart:iend,jstart:jend)+ s1s(:,:,ikx,iky) * eik

!                write(12347,*)"k12=",ikx,iky,i1,i2,j1,j2,k(1,ikx,iky),k(2,ikx,iky),dreal(eik),dimag(eik)
!                do i3=1,ns
!                  do j3=1,ns
!                    write(12347,*)"h0s=",i3,j3,dreal(h0s(i3,j3,ikx,iky)),dimag(h0s(i3,j3,ikx,iky)),h02(istart+i3-1,jstart+j3-1)
!                  enddo
!                enddo

              enddo
            enddo
!            write(12347,*)"end kloop ",istart,iend,jstart,jend
!            do i3=1,ns
!              do j3=1,ns
!                write(12347,*)"h02_out=",i3,j3,h02(istart+i3-1,jstart+j3-1)/(1.0_kdp * nkx*nky),h0(istart+i3-1,jstart+j3-1),(h02(istart+i3-1,jstart+j3-1)/(1.0_kdp * nkx*nky))-h0(istart+i3-1,jstart+j3-1)
!              enddo
!            enddo
!            write(12347,*)



          enddo
        enddo
      enddo
    enddo
!    write(12347,*)"end h0s output"
    h02=h02/(1.0_kdp * nkx*nky)
    h12=h12/(1.0_kdp * nkx*nky)
    s02=s02/(1.0_kdp * nkx*nky)
    s12=s12/(1.0_kdp * nkx*nky)

!    do i3=1,n
!      do j3=1,n
!        write(12347,*)"h02_final=",i3,j3,dreal(h0(i3,j3)),dimag(h0(i3,j3)),dreal(h02(i3,j3)),dimag(h02(i3,j3))
!      enddo
!    enddo
 

    write(12347,*)"maxdhinout=",maxval(abs(h0-h02)),maxval(abs(h1-h12)),maxval(abs(s0-s02)),maxval(abs(s1-s12))

!    do i1=1,n
!      do i2=1,n
!        write(12347,*)"h0b=",i1,i2,abs(h02(i1,i2)-h0(i1,i2)),h02(i1,i2),h0(i1,i2)
!      enddo
!    enddo

    if(overwritehs)then
      h0=h02
      h1=h12
      s0=s02
      s1=s12
    endif

    deallocate(h02,s02,h12,s12)

  end SUBROUTINE HSFourier

  subroutine surfacegf(n,gf,e,h0,s0,sigma)

    integer, intent(in)     :: n
    complex(kdp),intent(in) :: sigma(n,n),h0(n,n),s0(n,n)
    complex(kdp),intent(out) :: gf(n,n)
    complex(kdp),intent(in) :: e
 
    INTEGER , ALLOCATABLE :: ilu(:)
    complex(kdp), ALLOCATABLE :: work(:)
    integer nwork,info

    gf=e * s0 - h0 - sigma

    nwork=8
    allocate(work(nwork*n))
    allocate(ilu(n))
    call ZGETRF( n, n, gf, n, ilu, INFO )
    call ZGETRI( n, gf,n, ilu, WORK, nwork * n, INFO )
    deallocate(ilu)
    deallocate(work)

  end subroutine surfacegf

  SUBROUTINE WriteFourierPhiS(Phi,nam)

  use mTypes

  CHARACTER(LEN=*), intent(in) :: nam
  type(FourierScatteringStates), intent(in) :: Phi
  integer ikx,iky

  write(12347,*)nam,"ndivxy,n,ns=",Phi%ndivxy(1),Phi%ndivxy(2),Phi%n,Phi%ns,DREAL(Phi%e),DIMAG(Phi%e),Phi%side
  do ikx=1,Phi%ndivxy(1)
    do iky=1,Phi%ndivxy(2)
      write(12347,*)nam,"ikx,iky,kx,ky=",ikx,iky,Phi%kxy(1,ikx,iky),Phi%kxy(2,ikx,iky)
      call WritePhiS(Phi%FourierSstates(ikx,iky),nam)
    enddo
  enddo

  end SUBROUTINE WriteFourierPhiS


  SUBROUTINE WritePhis(Phi,nam)

  use negfmod, only: EM_NSPINBlocks

  type(ScatteringStates), intent(in) :: Phi
  CHARACTER(LEN=*), intent(in) :: nam
 
  integer i1,i2,n_in,n_out

  write(12347,*)nam,"nPhi=",Phi%n,Phi%nTstates(1),Phi%nTstates(2)

  do i1=1,Phi%nTstates(1)

    write(12347,*)nam,"k_in=",Phi%k_in(i1)
    write(12347,*)nam,"v_in=",Phi%v_in(i1)
    if(EM_NSPINBlocks==4)then
      write(12347,*)nam,"mu_in=",Phi%mu_in(1,i1),Phi%mu_in(2,i1),Phi%mu_in(3,i1),Phi%mu_in(4,i1)
    endif
    write(12347,*)nam,"phi_in=",Phi%phi_in(:,i1)

  enddo

  do i1=1,Phi%nTstates(2)
 
    write(12347,*)nam,"k_out=",Phi%k_out(i1)
    write(12347,*)nam,"v_out=",Phi%v_out(i1)
    if(EM_NSPINBlocks==4)then
      write(12347,*)nam,"mu_out=",Phi%mu_out(1,i1),Phi%mu_out(2,i1),Phi%mu_out(3,i1),Phi%mu_out(4,i1)
    endif
    write(12347,*)nam,"phi_out=",Phi%phi_out(:,i1)

  enddo


  end SUBROUTINE WritePhiS

  SUBROUTINE ExpandSstates(FourierPhiS,side,ndivxy,n,ns,e,k)

  use mTypes
  use mTransmissionDecompositionRemove, only : ConvertExpandFourierSstates

  type(FourierScatteringStates), intent(inout) :: FourierPhiS
  integer, intent(in)         :: n,ns
  CHARACTER(LEN=1),intent(in) :: SIDE
  complex(kdp),intent(in)     :: e
  integer, intent(in)         :: ndivxy(2)
  real(kdp),intent(in)        :: k(2,ndivxy(1),ndivxy(2))

  integer i1,i2,j1,j2,istart,iend,jstart,jend,ikx,iky
  integer il

  complex(kdp) :: eik
  real(kdp), parameter :: pi=3.14159265358979323846264338327950288_kdp
  complex(kdp), parameter :: ii=(0.0_kdp,1.0_kdp)
 
  
  allocate(FourierPhiS%Sstates(ndivxy(1),ndivxy(2)))

!xxx nkxy
  do ikx=1,ndivxy(1)
    do iky=1,ndivxy(2)
      call ExpandSstatesGet(FourierPhiS%FourierSstates(ikx,iky),FourierPhiS%sstates(ikx,iky),ndivxy,ikx,iky,k)
!      write(12347,*)"phi3=",FourierPhiS%Sstates(ikx,iky)%phi_in
    enddo
  enddo

  if(side.eq.'L')then
    il=1
  else
    il=2
  endif

  call ConvertExpandFourierSstates(il,FourierPhiS)

  end SUBROUTINE ExpandSstates





  SUBROUTINE ExpandSstatesGet(FourierSstates,Sstates,ndivxy,ikx,iky,k)

  use negfmod, only: EM_NSPINBlocks

  type(ScatteringStates), intent(in) :: FourierSstates
  type(ScatteringStates), intent(inout) :: Sstates
  integer, intent(in)         :: ndivxy(2)
  integer, intent(in)         :: ikx,iky
  real(kdp),intent(in)        :: k(2,ndivxy(1),ndivxy(2))

  integer ns,n,n_in,n_out
  integer i1,i2,j1,j2,istart,iend,jstart,jend
  complex(kdp) :: eik
  complex(kdp), parameter :: ii=(0.0_kdp,1.0_kdp)

  ns=FourierSstates%n
  n=ns*ndivxy(1)*ndivxy(2)

  Sstates%n=n
  Sstates%nTstates(1)=FourierSstates%nTstates(1)
  Sstates%nTstates(2)=FourierSstates%nTstates(2)

  n_in=Sstates%nTstates(1)
  n_out=Sstates%nTstates(2)

  if(n_in>0)then
    allocate(Sstates%phi_in(n,n_in))
    allocate(Sstates%phit_in(n,n_in))
    allocate(Sstates%k_in(n_in))
    allocate(Sstates%v_in(n_in))
    allocate(Sstates%mu_in(4,n_in))

    Sstates%k_in=FourierSstates%k_in
    Sstates%v_in=FourierSstates%v_in
    if(EM_NSPINBlocks==4)then
      Sstates%mu_in=FourierSstates%mu_in
    endif


  endif

  if(n_out>0)then
    allocate(Sstates%phi_out(n,n_out))
    allocate(Sstates%phit_out(n,n_out))
    allocate(Sstates%k_out(n_out))
    allocate(Sstates%v_out(n_out))
    allocate(Sstates%mu_out(4,n_out))

    Sstates%k_out=FourierSstates%k_out
    Sstates%v_out=FourierSstates%v_out
    if(EM_NSPINBlocks==4)then
      Sstates%mu_out=FourierSstates%mu_out
    endif

  endif


  do i1=1,ndivxy(1)
    do i2=1,ndivxy(2)

      istart=(i1-1) * ns+(i2-1) * ndivxy(1) * ns+1
      iend=i1 * ns+(i2-1) * ndivxy(1) * ns

      eik=exp(-ii * (k(1,ikx,iky)* (1-i1)+k(2,ikx,iky)* (1-i2)))
!      write(12347,*)"ikxy=",ikx,iky,i1,i2,n_in,dreal(eik),dimag(eik),k(1,ikx,iky),k(2,ikx,iky)
      if(n_in>0)then
        Sstates%phi_in(istart:iend,:)=eik * FourierSstates%phi_in
        Sstates%phit_in(istart:iend,:)=FourierSstates%phit_in/DCONJG(eik) 
      endif
      if(n_out>0)then
        Sstates%phi_out(istart:iend,:)=eik * FourierSstates%phi_out
        Sstates%phit_out(istart:iend,:)=FourierSstates%phit_out/DCONJG(eik) 
      endif

    enddo
  enddo
   
  if(n_in>0)then
!    Sstates%phit_in=Sstates%phit_in/(1.0D0 * ndivxy(1) * ndivxy(2))

    Sstates%phi_in=Sstates%phi_in/(1.0D0 * ndivxy(1) * ndivxy(2))
!    Sstates%phit_in=Sstates%phit_in*(1.0D0 * ndivxy(1) * ndivxy(2))
  endif
  if(n_out>0)then
!    Sstates%phit_out=Sstates%phit_out/(1.0D0 * ndivxy(1) * ndivxy(2))

    Sstates%phi_out=Sstates%phi_out/(1.0D0 * ndivxy(1) * ndivxy(2))
!    Sstates%phit_out=Sstates%phit_out*(1.0D0 * ndivxy(1) * ndivxy(2))
  endif


  end SUBROUTINE ExpandSstatesGet


  subroutine SurfacePDOS(gf2,s0,n,e)

  integer, intent(in)      :: n
  complex(kdp),intent(in)  :: gf2(n,n)
  complex(kdp),intent(in)  :: s0(n,n)
  complex(kdp),intent(in)  :: e

  complex(kdp),allocatable :: gf2s0(:,:)
  real(kdp)                :: dosgf
  integer                  :: i1
  real(kdp), parameter     :: pi=3.14159265358979323846264338327950288_kdp

  allocate(gf2s0(n,n))

  gf2s0=matmul(gf2,s0)
  
  dosgf=0.0_kdp
  do i1=1,n
    dosgf=dosgf-dimag(gf2s0(i1,i1))/(pi)
  enddo
  write(12347,*)"dos=",dreal(e),dosgf

  deallocate(gf2s0)


  end subroutine SurfacePDOS

end module mSigmaFourier

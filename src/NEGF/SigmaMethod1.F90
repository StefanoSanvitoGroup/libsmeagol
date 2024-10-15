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
!                   OPTIONS_KOFE,
!                   ALPHAOFZ,
!                   SVDNOISEK1K0,
!                   SOLVEQEP,
!                   SPLITLGRG,
!                   PRINTCOMPLEXBANDS,
!                   PRINTCOMPLEXBANDSVEC,
!                   SIGMAVINV,
!                   CHECKDSIGMA,
!                   KOFEC,
!                   EVSIGMA,
!                   DOSTRC,
!                   VG_DOS_DOSVV,
!                   SDSIGMADE,
!                   KOFEWRAP,
!                   KOFE_SVDLR2,
!                   CHECK_ERROR_SIGMA  
! AND
! THE MODULE
!                   MSIGMAMETHOD1  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

module mSigmaMethod1

  use mConstants

  implicit none
  private
  
  public :: GetSelfEnergy
  public :: kofewrap
  public :: printcomplexbands
  public :: printcomplexbandsvec
  public :: dostrc
  public :: check_error_sigma
  public :: evsigma
  
  contains

      
  subroutine options_kofe(usehinv,usevinv,tolki,svdtolmax, svdtolmin, dsigmamax,callsvd,rnoise,skipsvd,complexbands,eimag,svdtolzi,dosleads, dosk_local, doskvv_local, SigmaWideBand_local,NRunSigmaMax_local,dsigmade_local,TransmissionMatrix_local)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod, only: m_usehinv,m_usevinv,m_tolki,m_svdtolmax, m_svdtolmin, m_dsigmamax,m_callsvd,m_rnoise,&
                      &m_skipsvd,m_complexbands,deltaimag,m_svdtolzi,m_dosleads,&
                      &dosk, doskvv, SigmaWideBand,NRunSigmaMax,dsigmade,TransmissionMatrix

    implicit none
    LOGICAL usehinv,usevinv,callsvd,complexbands, dosleads
    DOUBLE PRECISION tolki,svdtolmax,svdtolmin,dsigmamax,rnoise, skipsvd,eimag,svdtolzi
    logical  dsigmade_local,TransmissionMatrix_local
    DOUBLE PRECISION :: dosk_local,doskvv_local,SigmaWideBand_local
    integer :: NRunSigmaMax_local



    usehinv=      m_usehinv      
    usevinv=      m_usevinv
    tolki=        m_tolki
    svdtolmax=    m_svdtolmax
    svdtolmin=    m_svdtolmin
    dsigmamax=    m_dsigmamax
    callsvd=      m_callsvd
    rnoise=       m_rnoise
    skipsvd=      m_skipsvd
    complexbands= m_complexbands
    eimag       = deltaimag 
    svdtolzi    = m_svdtolzi
    dosleads    = m_dosleads
    dosk_local              = dosk
    doskvv_local            = doskvv
    SigmaWideBand_local     = SigmaWideBand
    NRunSigmaMax_local      = NRunSigmaMax
    dsigmade_local          = dsigmade
    TransmissionMatrix_local= TransmissionMatrix

  end subroutine options_kofe


  SUBROUTINE GetSelfEnergy(SIDE,N2,Ei,H0,H1,S0,S1,Sigma,nrchan,delta)

! **********************************************************************
! Calculates the self-energies, based on the singularity-free scheme,
! Ref: I. Rungger and S. Sanvito, Phys Rev B 78, 035407 (2008)
!
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

  use negfmod, only: pdosgs, skipright
 
  IMPLICIT NONE

  INTEGER, intent(in)          :: N2
  CHARACTER(LEN=1), intent(in) :: SIDE
  DOUBLE COMPLEX, intent(in)   :: Ei
  DOUBLE PRECISION, intent(in) :: delta
  DOUBLE COMPLEX, DIMENSION (N2,N2), intent(in) ::  S0,S1,H0,H1
  INTEGER, intent(out)         :: nrchan
  DOUBLE COMPLEX, DIMENSION (N2,N2), intent(out) :: Sigma

  DOUBLE COMPLEX, ALLOCATABLE :: k0(:,:),k1(:,:), km1(:,:)
  double complex              :: Ei0
  DOUBLE PRECISION            :: dsigma
  DOUBLE COMPLEX, PARAMETER   :: zi=(0.D0,1.D0)

  allocate(k0(N2,N2),k1(N2,N2),km1(N2,N2))

  Ei0=Ei+ zi * delta

  k0=H0 - Ei0*S0
  k1=H1 - Ei0*S1
  km1=DCONJG(TRANSPOSE(H1)) - Ei0*DCONJG(TRANSPOSE(S1))

  if(skipright.and.side .eq. 'R' )return

!on output of kofewrap "k0" contains the values of the self-energy
  call kofewrap(side,'S',k0,k1,km1,N2,Ei0,dsigma,nrchan,S0,S1)

!the selfenergy is found in the k0 matrix and returned in sigma
  sigma=k0
  deallocate(k0,k1,km1)

  if(pdosgs)then
    if(side.eq.'L')then
      call PDOSLeads(n2,H0,H1,S0,S1,Ei0,sigma)
    endif
  endif

  END SUBROUTINE GetSelfEnergy


  subroutine PDOSLeads(n,H0,H1,S0,S1,ene,sigma)

  INTEGER, intent(in)          :: n
  DOUBLE COMPLEX, intent(in)   :: ene
  DOUBLE COMPLEX, DIMENSION (n,n), intent(in) ::  S0,S1,H0,H1,sigma
  
  INTEGER        :: j
  DOUBLE COMPLEX, ALLOCATABLE :: k0(:,:),k1(:,:), km1(:,:)

  INTEGER   info, IPIV(n)
  DOUBLE COMPLEX tracegs0, tracegsm1
  DOUBLE PRECISION dsigma,dsigmar, maxsigma,dgammam
  DOUBLE COMPLEX, allocatable :: Gr_2(:,:),gr_2t(:,:),rhom1(:,:),rho0(:,:),k1ob(:,:)
  DOUBLE COMPLEX, allocatable :: work(:)
 
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793D0

  allocate(k0(n,n))
  k0=H0 - ene*S0
  k1=H1 - ene*S1
  km1=DCONJG(TRANSPOSE(H1)) - ene*DCONJG(TRANSPOSE(S1))

  dgammam=0D0
  do j=1,n
    dgammam=dgammam+DREAL(sigma(j,j))
  enddo

  k0=-k0-sigma
  CALL ZGETRF(n,n,k0,n,IPIV,INFO)
  allocate(work(n**2))
  CALL ZGETRI(n,k0,n,IPIV,work,n**2,INFO)
  deallocate(work)
  allocate(gr_2(n,n))
  gr_2=k0

  allocate(k1ob(n,n))
  k1ob=matmul(km1,gr_2)
  k0=matmul(k1ob,k1)
  k1ob=sigma-k0
  deallocate(k0)
  dsigma=MAXVAL(ABS(k1ob))
  deallocate(k1ob)
  maxsigma=MAXVAL(ABS(sigma))
  dsigmar=dsigma/maxsigma
  
  allocate(gr_2t(n,n))
  allocate(rhom1(n,n))
  allocate(rho0(n,n))
  gr_2t=DCONJG(TRANSPOSE(gr_2))
  rho0=(0.5D0 * zi /pi) * (gr_2-gr_2t)
  rhom1=(0.5D0 * zi /pi) * (matmul(gr_2t,matmul(DCONJG(TRANSPOSE(k1)),gr_2))-matmul(gr_2,matmul(km1,gr_2)))
  deallocate(gr_2,gr_2t)

  rho0=matmul(rho0,S0)
  rhom1=matmul(rhom1,S1)

  tracegs0=0D0
  do j=1,n
    tracegs0=tracegs0+rho0(j,j)
  enddo
  tracegsm1=0D0
  do j=1,n
    tracegsm1=tracegsm1+rhom1(j,j)
  enddo
  deallocate(rhom1,rho0)
  write(12347,'(a," ",e12.5," ",9(e12.5," "))') "trrss=",DREAL(ene), DREAL(tracegs0),DREAL(tracegs0+tracegsm1),  DREAL(tracegsm1),DIMAG(tracegsm1),dsigma,dsigmar, maxsigma,DIMAG(ene),dgammam

  end subroutine PDOSLeads


  subroutine kofewrap(side,gf,k0,k1in,km1in,n,ene, dsigma,nrchan,S0,S1)


! **********************************************************************
! Calculates the self-energies, based on the singularity-free scheme,
! Ref: I. Rungger and S. Sanvito, Phys Rev B 78, 035407 (2008)
!
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod, only : em_Last_SCF_Step,ikpmod
 
    use mComputeULR, only :ComputeVgMuPhi
    IMPLICIT NONE

    INTEGER n,nrchan,i
    CHARACTER(LEN=1) :: side,gf
    DOUBLE COMPLEX, DIMENSION (n,n) :: k0,k1in,k0in,km1in,vrl,S0,S1, vrlt
    DOUBLE COMPLEX, DIMENSION(n) :: zvo
    DOUBLE COMPLEX  ene
    DOUBLE PRECISION  dsigma,svdtolzi2,tolki
    LOGICAL, ALLOCATABLE, SAVE :: firstcall(:)
    LOGICAL, SAVE :: usehinv, usevinv ,complexbands,callsvd,dosleads
    DOUBLE PRECISION, SAVE :: tolkisave, svdtolmax,svdtolmin,dsigmamax,rnoise,skipsvd,svdtolzi,eimag
    logical, save :: dsigmade,TransmissionMatrix
    DOUBLE PRECISION, SAVE :: dosk,doskvv,SigmaWideBand
    DOUBLE COMPLEX, ALLOCATABLE :: vlevout(:,:)
    integer, save :: NRunSigmaMax


    if(.not.allocated(firstcall))then
      allocate(firstcall(1))

      call options_kofe(usehinv,usevInv,tolkisave,svdtolmax, svdtolmin, dsigmamax,callsvd,rnoise,skipsvd,complexbands,eimag,svdtolzi,dosleads,dosk, doskvv, SigmaWideBand,NRunSigmaMax,dsigmade,TransmissionMatrix)
      dosk=0D0
      doskvv=0D0

      firstcall(1)=.false.
    endif
    svdtolzi2=svdtolzi
    tolki=tolkisave

    k0in=k0
    allocate(vlevout(1,1))
    if(svdtolzi2.eq.0D0)then
      if(SigmaWideBand==0.0D0)then
        if(dsigmade) then
          deallocate(vlevout)
          allocate(vlevout(n,n))
        endif
        call kofec(side,gf,k0,k1in,km1in,n,ene, dsigma,nrchan,vrl,vrlt,zvo,usehinv,usevinv,tolki,svdtolmax,svdtolmin, dsigmamax,callsvd,rnoise,skipsvd,complexbands,dsigmade,em_Last_SCF_Step,ikpmod,vlevout)
      else
        write(12347,*)"Wide band approximation for self-energies",SigmaWideBand
        k0=0.0D0
        do i=1,n
          k0(i,i)=(0.0D0,-1.0D0)*SigmaWideBand
        enddo
      endif
    else
        call kofe_svdlr2(side,k0,k1in,km1in,n,ene,svdtolzi2, dsigma,nrchan,SigmaWideBand,usehinv,usevinv,tolki,svdtolmax,svdtolmin, dsigmamax,callsvd,rnoise,skipsvd,complexbands,dsigmade,em_Last_SCF_Step,ikpmod,vlevout)
        if(dsigma.gt.dsigmamax.and.dsigma.gt.1D1 * eimag)then
          write(12347,*)"svdtr1= 0,",svdtolzi2,dsigma,1D1 * eimag
          do i=1,NRunSigmaMax
            svdtolzi2=svdtolzi2 / 1D2
            write(*,*)"warning: changing SVD tolerance for selfenergy to",svdtolzi2
            write(12347,*)"warning: changing SVD tolerance for selfenergy to",svdtolzi2
            if(svdtolzi2.lt.1d-16)exit
            k0=k0in
            call kofe_svdlr2(side,k0,k1in,km1in,n,ene, svdtolzi2, dsigma,nrchan,SigmaWideBand,usehinv,usevinv,tolki,svdtolmax,svdtolmin, dsigmamax,callsvd,rnoise,skipsvd,complexbands,dsigmade,em_Last_SCF_Step,ikpmod,vlevout)
            if(dsigma.lt.dsigmamax.or.dsigma.lt.1D1 * eimag) then
              write(*,*)"Selfenergy calculated to the required accuracy with the updated SVD tolerance"
              write(12347,*)"Selfenergy calculated to the required accuracy with the updated SVD tolerance"
              exit
            endif
            write(12347,*)"svdtr1=",i,svdtolzi2,dsigma,1D1 * eimag
            if(i==3)then
             write(*,*)"warning: selfenergy was not calculated to the required accuracy"
             write(12347,*)"warning: selfenergy was not calculated to the required accuracy"
            endif
          enddo
          write(12347,*)"svdtexit=",svdtolzi2,dsigma
        endif
    endif

    if(dsigmade) call sdsigmade(n,ene,zvo,vrl,vrlt,S0,S1,k0in,k1in,km1in,k0,side,vlevout)
    if((side.eq.'L').and.dosleads)call vg_dos_dosvv(n,ene,zvo,vrl,S0,S1,km1in,k1in,tolki)

    if(em_Last_SCF_Step.and.TransmissionMatrix)then
      call ComputeVgMuPhi(side,n,ene,S0,S1,km1in,k1in)
    endif
    deallocate(vlevout)

  end subroutine kofewrap


  subroutine kofec(side,gf,k0,k1in,km1in,n,ene, dsigmam,nrchanm,vrgout, vrgbout,zvrout,usehinv,usevinv,tolki,svdtolmax,svdtolmin, dsigmamax,callsvdin,rnoise,skipsvd,complexbands,dsigmade,em_Last_SCF_Step,ik,vlevout)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod, only : TransmissionMatrix
    use mComputeULR, only : ComputeULR

    IMPLICIT NONE

    LOGICAL, intent(in) :: usehinv, usevinv ,callsvdin,complexbands,dsigmade,em_Last_SCF_Step
    DOUBLE PRECISION, intent(in) :: svdtolmax,svdtolmin,dsigmamax,tolki, rnoise,skipsvd
    integer, intent(in) :: ik
    INTEGER n,nrchan,nrchanm,i,j
    INTEGER   i1,i2,vml(n),vmr(n),il,iloc,iroc, nls,nrs
    INTEGER   INFO,IPIV(n)


    DOUBLE PRECISION :: svdtol,dsigma,signs,rworkc(2 * n), dsigmam,svdtolm,anorml,anormr,rcondl,rcondr, anormgf,rcondgf,dsigmar,dsigmarm,pi

    DOUBLE COMPLEX, DIMENSION(2 * n,2 * n) :: m2s,m2sb
    DOUBLE COMPLEX, DIMENSION(n,2 * n) :: vr,vl
    DOUBLE COMPLEX, DIMENSION (n,n) :: k0,k1,grfunct,k0input, matout,k1in,km1in,km1,al,ar,vlg,vrg,vlgb,vrgb,vrgout, vrgbout,sigma,k1o,k1ob
    double complex, allocatable :: vlev(:,:)
    DOUBLE COMPLEX, DIMENSION (n,n) :: vlevout
    DOUBLE COMPLEX, DIMENSION(n) :: vri,vritd,vli,vlitd,zvl,zvr,zvrout
    DOUBLE COMPLEX worki(n**2),walpha(2 * n),zv(2 * n)

    DOUBLE COMPLEX  ene,zvu

    CHARACTER(LEN=1) :: side,gf
    LOGICAL  callsvd,usehinv1,calcsucceed,wrotebands,computedulr
    DOUBLE PRECISION ZLANGE

!    EXTERNAL  zlange


    pi=3.141592653589793D0
    nrchanm=0
    usehinv1=.false.
    calcsucceed=.false.
    wrotebands=.false.
    computedulr=.false.
!    write(*,*)"dsigmam=",dsigmamax,rnoise,skipsvd
    callsvd=callsvdin


    If (side .eq. 'R') then
      signs=-1D0
    else
      signs=1D0
    endif

    If (side .eq. 'R') then
      k0input=km1in
      km1in=k1in
      k1in=k0input
    endif
    k0input=k0

    svdtol=0D0
    svdtolm=svdtol
    dsigmam=0D0
    dsigmarm=0D0

    1020  k0=k0input
    k1=k1in
    km1=km1in

    call svdnoisek1k0(k0,k1,km1, callsvd,svdtolm,svdtol,svdtolmin,skipsvd,usehinv,n, rnoise,k1o,ipiv, worki)


    1100 call solveqep(usehinv,usehinv1,dsigmade,m2s,m2sb,k1o, k0,k1,km1,zv,n,info,vl,vr)


    if(info.ne.0.and.svdtol.lt.svdtolmax)then 
      write(12347,*)"warning: increasing svdtol info", DREAL(ene),svdtol,info
      goto 1020
    elseif(info.ne.0)then 
      write(12347,*)"warning: last not increasing svdtol info",  DREAL(ene),svdtol,info
      goto 1030
    endif


    call splitlgrg(n,nrchan,zv,k1,km1,k1o,signs,vml,vmr, ene,nrs,nls,side,svdtol,vr,iroc,iloc,tolki)

    if((nrs.ne.nls).or.(nrs+nls.ne.2 * n))then
      write(12347,*)"warning:nrs,nls=", DREAL(ene),nrs,nls,iroc,iloc,side,svdtol

      if(svdtol.lt.svdtolmax)then 
        write(12347,*)"warning:increasing svdtol",  DREAL(ene),svdtol
        goto 1020
      else 
        write(12347,*)"warning:problem in last",DREAL(ene)
        goto 1030
      endif
    endif

    if(allocated(vlev))deallocate(vlev)
    if(dsigmade)then
      allocate(vlev(n,n))
    else
      allocate(vlev(1,1))
    endif
    do il=1,n
      vrg(:,il)=vr(:,vmr(il))
      vlg(:,il)=vr(:,vml(il))
      zvr(il)=zv(vmr(il))
      zvl(il)=zv(vml(il))
      if(dsigmade)then
        if(side.eq.'L')then
          vlev(:,il)=vl(:,vml(il))
        else
          vlev(:,il)=vl(:,vmr(il))
        endif
      endif
    enddo



    vrgb=vrg
    anormr=ZLANGE('1',n,n,vrgb,n,rworkc)
    CALL ZGETRF(n,n,vrgb,n,IPIV,INFO)
    CALL ZGECON('1',n,vrgb,n,anormr,rcondr,walpha,rworkc,INFO)
    CALL ZGETRI(n,vrgb,n,IPIV,worki,n**2,INFO)

    vlgb=vlg
    anorml=ZLANGE('1',n,n,vlgb,n,rworkc)
    CALL ZGETRF(n,n,vlgb,n,IPIV,INFO)
    CALL ZGECON('1',n,vlgb,n,anorml,rcondl,walpha,rworkc,INFO)
    CALL ZGETRI(n,vlgb,n,IPIV,worki,n**2,INFO)


    If (side .eq. 'L') Then

      al=0D0
      do il=1,n
        vli(:)=vlg(:,il)
        vlitd(:)=vlgb(il,:)
        do i1=1,n
          do i2=1,n
            al(i1,i2)=al(i1,i2)+ vli(i1) * vlitd(i2) / zvl(il)
          enddo
        enddo
      enddo

    else

      ar=0D0
      do il=1,n
        zvu=1D0/zvr(il)
        vri(:)=vrg(:,il)
        vritd(:)=vrgb(il,:)
        do i1=1,n
          do i2=1,n
            ar(i1,i2)=ar(i1,i2)+ vri(i1) * vritd(i2) * zvu
          enddo
        enddo
      enddo

    endif



    if(usevinv)then

      call sigmavinv(side,n,ar,al,vrg,vrgb,vlg,vlgb,zvr,zvl, grfunct,sigma,k1,km1,k1o,k1ob,ipiv,worki,signs)

    else
      if(side.eq.'L')then
        sigma=matmul(km1,al)
!        sigma=matmul(DCONJG(TRANSPOSE(k1in)),al)
      else
        sigma=matmul(km1,ar)
!        sigma=matmul(DCONJG(TRANSPOSE(k1in)),ar)
      endif

    endif


!    If (side .eq. 'L') Then
      if(complexbands.and.(.not.wrotebands))then
        call printcomplexbands(side,n,zvr,zvl,ene,svdtol,ik)
        wrotebands=.true.
      endif
!    endif

    if(em_Last_SCF_Step.and.TransmissionMatrix.and.(.not.computedulr))then
      call ComputeULR(side,n,zvr,zvl,ene,vrg,vrgb,vlg,vlgb,sigma,tolki)
!      computedulr=.true.
    endif


    grfunct=-k0input-sigma

    anormgf=ZLANGE('1',n,n,grfunct,n,rworkc)
    CALL ZGETRF(n,n,grfunct,n,IPIV,INFO)
    CALL ZGECON('1',n,grfunct,n,anormgf,rcondgf,walpha,rworkc,INFO)
    CALL ZGETRI(n,grfunct,n,IPIV,worki,n**2,INFO)

    k1ob=matmul(km1in,grfunct)
    k1o=matmul(k1ob,k1in)
    k1ob=sigma-k1o
    dsigma=MAXVAL(ABS(k1ob))
    dsigmar=dsigma/MAXVAL(ABS(sigma))

    call checkdsigma(rcondr,rcondl,rcondgf,svdtol, dsigma,dsigmar,dsigmamax,dsigmarm,dsigmam,svdtolm, n,k1o,sigma,anorml,anormr,ene,svdtolmax, anormgf,calcsucceed,side,gf,nrchanm,vrgout,vrgbout, zvrout, nrchan,tolki,vlg,vlgb,zvl,vlev, vrg,vrgb,zvr, matout,grfunct, info,ik,dsigmade,vlevout)
    deallocate(vlev)

    if(info.ne.0)goto 1020


    1030  If (side .eq. 'R') then
      k0input=km1in
      km1in=k1in
      k1in=k0input
    endif

    if(.not.calcsucceed)then
      write(*,*)"warning: calculation of selfenergies failed for the used SVD tolerance"
      write(*,*)"Side ",side," energy= ",dreal(ene),dimag(ene)
      write(*,*)"warning: using wide-band limit for the selfenergy for the used SVD tolerance"
      write(12347,*)"warning: calculation of selfenergies failed for the used SVD tolerance"
      write(12347,*)"Side ",side," energy= ",dreal(ene),dimag(ene)
      write(12347,*)"warning: using wide-band limit for the selfenergy for the used SVD tolerance"

      dsigmam=10D0
      dsigmarm=10D0

      k0=0.0D0
      do i=1,n
        k0(i,i)=(0.0D0,-1.0D0)*1.0e-2
      enddo

      return
    endif

    if(dsigmarm.gt.1d3 * dsigmamax.or.dsigmam.gt.dsigmamax)then
      write(12347,*)"f,sv=",side,dreal(ene),svdtolm,dsigmam,dsigmarm, MAXVAL(ABS(matout)),ik
    endif


    k0=matout


  end subroutine kofec


  subroutine kofe_svdlr2(side,k0,k1,km1,n,ene, svdtol, dsigma,nrchan,SigmaWideBand,usehinv,usevinv,tolki,svdtolmax,svdtolmin, dsigmamax,callsvd,rnoise,skipsvd,complexbands,dsigmade,em_Last_SCF_Step,ik,vlevout)

! **********************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
 
    use negfmod, only: TransmissionMatrix
    use mComputeULR, only : PhiS,ComputePhiSVD

    IMPLICIT NONE
    INTEGER n,m,m2,nrchan
    integer, intent(in) :: ik
    CHARACTER(LEN=1) :: side
    DOUBLE COMPLEX, DIMENSION (n,n) :: k0,k1,grf,km1
    double precision, intent(in):: SigmaWideBand
    LOGICAL usehinv,usevinv,callsvd,complexbands,dsigmade,em_Last_SCF_Step
    DOUBLE PRECISION tolki,svdtolmax,svdtolmin,dsigmamax,rnoise, skipsvd
    double complex vlevout(n,n)

    INTEGER   i,i1
    double complex maxv
    DOUBLE PRECISION dsigma
    INTEGER   INFO,nwork
    DOUBLE COMPLEX, ALLOCATABLE :: k1o(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: mata(:,:),matb(:,:),matc(:,:),matd(:,:),matdi(:,:), v1c(:,:),v1n(:,:),vm1c(:,:),vm1n(:,:), v1e(:,:),vm1e(:,:),v0e(:,:)
    INTEGER , ALLOCATABLE :: ilu(:)
    DOUBLE COMPLEX, ALLOCATABLE :: k0t(:,:),k1t(:,:),km1t(:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: k0r2(:,:),vrlt(:,:),zvot(:), vrlt2(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
    DOUBLE COMPLEX, ALLOCATABLE ::u(:,:),vt(:,:),u2(:,:),vt2(:,:), work(:)
    DOUBLE PRECISION, ALLOCATABLE :: sv(:),sv2(:)
    DOUBLE PRECISION :: svmax,svdtol
    DOUBLE COMPLEX  ene

    nwork=8
    allocate(k0t(n,n),k1t(n,n),km1t(n,n),u(n,n),vt(n,n),u2(n,n), vt2(n,n),k1o(n,n),sv(n),sv2(n),work(nwork*n),rwork(5 * n))

    k1o=k1
    if(side.eq.'L')then
      CALL ZGESVD( 'N', 'A', n, n, k1o, n,sv, u, n, vt,n, work, nwork * n, rwork, INFO )
    else
      CALL ZGESVD( 'A', 'N', n, n, k1o, n,sv, u, n, vt,n, work, nwork * n, rwork, INFO )
    endif

    svmax=sv(1)

    m=n
    do i=1,n
      if(sv(i)/svmax.lt.svdtol)then
        m=i-1
        exit
      endif
    enddo

    k1o=DCONJG(TRANSPOSE(km1))
    if(side.eq.'L')then
      CALL ZGESVD( 'N', 'A', n, n, k1o, n,sv2, u2, n, vt2,n, work, nwork * n, rwork, INFO )
    else
      CALL ZGESVD( 'A', 'N', n, n, k1o, n,sv2, u2, n, vt2,n, work, nwork * n, rwork, INFO )
    endif

    deallocate(rwork)

    svmax=sv2(1)

    m2=n
    do i=1,n
      if(sv2(i)/svmax.lt.svdtol)then
        m2=i-1
        exit
      endif
    enddo

    write(12347,*)"msvd2=",DREAL(ene),m,m2,n,svdtol

!    rmmax=0D0
!    do i=1,n
!      if(abs(sv(i)-sv2(i)).gt.rmmax)rmmax=abs(sv(i)-sv2(i))
!    enddo
!    write(12347,*)"sdiffmax=",rmmax

    if(m2.lt.m)m=m2
    
    if(m.eq.n)then
!      write(*,*)"no need to decimate, k1 is invertable"
!      call f77flush
      deallocate(k0t,k1t,km1t)
      allocate(v0e(n,n),v1e(n,n),vm1e(n,n))
      v0e=k0
      v1e=k1
      vm1e=km1
      goto 1000
    endif


    if(side.eq.'R')then
      allocate(v1c(m,m),v1n(m,n-m),vm1c(m,m),vm1n(n-m,m))
      k1t=matmul(DCONJG(TRANSPOSE(u)),matmul(k1,u2))
      km1t=matmul(DCONJG(TRANSPOSE(u)),matmul(km1,u2))
      k0t=matmul(DCONJG(TRANSPOSE(u)),matmul(k0,u2))
 
      v1c=k1t(1:m,1:m)
      v1n=k1t(1:m,m+1:n)
      vm1c=km1t(1:m,1:m)
      vm1n=km1t(m+1:n,1:m)
    else
      allocate(v1c(m,m),v1n(n-m,m),vm1c(m,m),vm1n(m,n-m))
      k1t=matmul(vt2,matmul(k1,DCONJG(TRANSPOSE(vt))))
      km1t=matmul(vt2,matmul(km1,DCONJG(TRANSPOSE(vt))))
      k0t=matmul(vt2,matmul(k0,DCONJG(TRANSPOSE(vt))))
 
      v1c=k1t(1:m,1:m)
      v1n=k1t(m+1:n,1:m)
      vm1c=km1t(1:m,1:m)
      vm1n=km1t(1:m,m+1:n)
    endif


    allocate(mata(m,m),matb(m,n-m),matc(n-m,m),matd(n-m,n-m), matdi(n-m,n-m),v1e(m,m),vm1e(m,m),v0e(m,m))
    mata=k0t(1:m,1:m)
    matb=k0t(1:m,m+1:n)
    matc=k0t(m+1:n,1:m)
    matd=k0t(m+1:n,m+1:n)
    matdi=matd
    allocate(ilu(n-m))
    call ZGETRF( n-m, n-m, matdi, n-m, ilu, INFO )
    call ZGETRI( n-m, matdi,n-m, ilu, WORK, nwork * n, INFO )
    deallocate(ilu)
    deallocate(k0t,k1t,km1t)

    if(side.eq.'R')then
      v1e=v1c-matmul(v1n,matmul(matdi,matc))
      vm1e=vm1c-matmul(matb,matmul(matdi,vm1n))
      v0e=mata-matmul(matb,matmul(matdi,matc)) -matmul(v1n,matmul(matdi,vm1n))
    else
      v1e=v1c-matmul(matb,matmul(matdi,v1n))
      vm1e=vm1c-matmul(vm1n,matmul(matdi,matc))
      v0e=mata-matmul(matb,matmul(matdi,matc)) -matmul(vm1n,matmul(matdi,v1n))
    endif

    1000  allocate(k0r2(m,m))
    allocate(vrlt(m,m),vrlt2(m,m),zvot(m))

    k0r2=v0e
    call addnoise(v1e,m,1D-5 * svdtol * svmax)
    call addnoise(vm1e,m,1D-5 * svdtol * svmax)
    call addnoise(k0r2,m,1D-5 * svdtol * svmax)


    if(SigmaWideBand==0.0D0)then
      call kofec(side,'S',k0r2,v1e,vm1e,m,ene, dsigma,nrchan,vrlt,vrlt2,zvot,usehinv,usevinv,tolki,svdtolmax,svdtolmin, dsigmamax,callsvd,rnoise,skipsvd,complexbands,dsigmade,em_Last_SCF_Step,ik,vlevout)
    else
      write(12347,*)"Wide band approximation for self-energies",SigmaWideBand
      k0r2=0.0D0
      do i=1,m
        k0r2(i,i)=(0.0D0,-1.0D0)*SigmaWideBand
      enddo
    endif

    if(m.eq.n)then
      k0=k0r2
      deallocate(u,vt,u2,vt2,k1o,sv,sv2,work)
      deallocate(vrlt,vrlt2,zvot)
      deallocate(v0e,v1e,vm1e,k0r2)
      return
    endif

    if(em_Last_SCF_Step.and.TransmissionMatrix)then
      if(side.eq.'L')then
        call ComputePhiSVD(side,n,m,matdi,matc,v1n,vt,PhiS(1)%FourierSstates(PhiS(1)%ik(1),PhiS(1)%ik(2)),1)
        call ComputePhiSVD(side,n,m,matdi,matc,v1n,vt,PhiS(1)%FourierSstates(PhiS(1)%ik(1),PhiS(1)%ik(2)),2)
      else
        call ComputePhiSVD(side,n,m,matdi,matc,vm1n,u,PhiS(2)%FourierSstates(PhiS(2)%ik(1),PhiS(2)%ik(2)),1)
        call ComputePhiSVD(side,n,m,matdi,matc,vm1n,u,PhiS(2)%FourierSstates(PhiS(2)%ik(1),PhiS(2)%ik(2)),2)
      endif
    endif
   
    if(side.eq.'R')then
      v0e=-matmul(v1n,matmul(matdi,vm1n))
    else
      v0e=-matmul(vm1n,matmul(matdi,v1n))
    endif

    k0r2=k0r2+v0e

    k1o=0D0
    k1o(1:m,1:m)=k0r2
    if(side.eq.'R')then
      k1o=matmul(matmul(u,k1o),DCONJG(TRANSPOSE(u2)))
    else
      k1o=matmul(matmul(DCONJG(TRANSPOSE(vt2)),k1o),vt)
    endif

    grf=-k0-k1o
    allocate(ilu(n))
    call ZGETRF( n, n, grf, n, ilu, INFO )
    call ZGETRI( n, grf,n, ilu, WORK, nwork * n, INFO )
    deallocate(ilu)
    if(side.eq.'L')then
      grf=matmul(km1,matmul(grf,k1))
    else
      grf=matmul(k1,matmul(grf,km1))
    endif
    grf=k1o-grf
    dsigma=MAXVAL(ABS(grf))


    k0=k1o

    deallocate(v1c,vm1c)
    deallocate(mata,matb,matc,matd)


    deallocate(u,vt,u2,vt2,k1o,sv,sv2,work,v1n,vm1n, matdi,v1e,vm1e,v0e,k0r2)
    deallocate(vrlt,vrlt2,zvot)


  end subroutine kofe_svdlr2


  subroutine alphaofz(alpha,kappa,zev)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************


    DOUBLE PRECISION   alpha,kappa,lz,zr,zi,alpha0,pi
    DOUBLE COMPLEX     zev

    pi=3.141592653589793D0

    lz=CDABS(zev)
    zr=DREAL(zev)
    zi=DIMAG(zev)
    if(zr.ge.0D0.and.zi.ge.0D0) then
      alpha0=0D0
    elseif(zr.ge.0D0.and.zi.lt.0D0) then
      zr=-DIMAG(zev)
      zi=DREAL(zev)
      alpha0=-pi/2D0
    elseif(zr.lt.0D0.and.zi.ge.0D0) then
      zr=DIMAG(zev)
      zi=-DREAL(zev)
      alpha0=pi/2D0
    elseif(zr.lt.0D0.and.zi.lt.0D0) then
      zr=-DREAL(zev)
      zi=-DIMAG(zev)
      alpha0=-pi
    endif

    alpha=ATAN(zi/zr)+alpha0
!    alpha=ATAN((zi/lz)/(zr/lz))+alpha0
    kappa=DLOG(lz)

  end subroutine alphaofz



  subroutine svdnoisek1k0(k0,k1,km1, callsvd,svdtolm,svdtol,svdtolmin,skipsvd,usehinv,n, rnoise,k1o,ipiv, worki)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    IMPLICIT NONE

    integer     n,info
    INTEGER   IPIV(n)
    DOUBLE COMPLEX   worki(n**2)
    double complex k1(n,n),km1(n,n),k0(n,n),k1o(n,n)
    logical callsvd,usehinv
    double precision svdtol,svdtolmin,skipsvd,rnoise,svdtolm


    if(callsvd)then
      if(svdtol.eq.0D0)svdtol=svdtolmin/skipsvd
      svdtol=svdtol * skipsvd
!      write(*,*)"svdtol=",svdtol
      if(.not.usehinv)then
        call svdm(k1,k1o,n,svdtol,rnoise * svdtol,.false.,.true.)
!        call addnoise(k1,n,rnoise * svdtol)
      else
        call svdm(k1,k1o,n,svdtol,rnoise * svdtol,.true.,.false.)

!        call addnoise(k0,n,rnoise * svdtol)
!        call addnoise(k1,n,rnoise * svdtol)
!        call addnoise(km1,n,rnoise * svdtol)
!        k1o=k1
!        CALL ZGETRF(n,n,k1o,n,IPIV,INFO)
!        CALL ZGETRI(n,k1o,n,IPIV,worki,n**2,INFO)

      endif
!      k0=0.5D0 * (k0+DCONJG(TRANSPOSE(k0)))
      if(rnoise.ne.0D0)then
        call addnoise(k0,n,rnoise * svdtol)
        call addnoise(km1,n,rnoise * svdtol)
      endif
    else
      callsvd=.true.
      svdtol=0D0
      svdtolm=svdtol

      if(usehinv)then
!        k0=0.5D0 * (k0+DCONJG(TRANSPOSE(k0)))
        k1o=k1
        CALL ZGETRF(n,n,k1o,n,IPIV,INFO)
        CALL ZGETRI(n,k1o,n,IPIV,worki,n**2,INFO)
      endif
      if(rnoise.ne.0D0)then
        call addnoise(k0,n,rnoise * 1d1 * svdtolmin)
!        If (side .eq. 'R') call addnoise(k1,n,rnoise * 1d-17)
!        If (side .eq. 'L') call addnoise(km1,n,rnoise * 1d-17)
        call addnoise(k1,n,rnoise * 1d1 * svdtolmin)
        call addnoise(km1,n,rnoise * 1d1 * svdtolmin)
      endif
    endif

!    km1=DCONJG(transpose(k1))



  end subroutine svdnoisek1k0


  subroutine solveqep(usehinv,usehinv1,dsigmade,m2s,m2sb,k1o, k0,k1,km1,zv,n,info,vl,vr)
    
! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    IMPLICIT NONE

    integer i,i1,n,info
    double complex k1o(n,n),k0(n,n),k1(n,n),km1(n,n),m2s(2*n,2*n), zv(2*n),vr2(2 * n, 2 * n),vl2(2 * n, 2 * n),m2sb(2*n,2*n), vr(n,2*n),vl(n,2*n)
    logical usehinv,usehinv1,dsigmade
    DOUBLE COMPLEX   walpha(2 * n),wbeta(2 * n)
    DOUBLE PRECISION  normvr,normvr1,normvr2
    double complex maxv


    if(usehinv.or.usehinv1)then
      usehinv1=.false.

      m2s=0D0
      m2s(1:n,1:n)=-matmul(k1o,k0)
      m2s(1:n,1+n:2*n)=-matmul(k1o,km1)
!      m2s(1+n:2 * n,1:2 * n)=0D0
      do i=1,n
        m2s(n+i,i)=1D0
      enddo
      call geigenvalues2(m2s,zv,2 * n,vr2,info)

    else
      m2s=0D0
      m2s(1:n,1:n)=k0(1:n,1:n)
      m2s(1:n,1+n:2*n)=km1(1:n,1:n)
      do i=1,n
        m2s(n+i,i)=1D0
      enddo
      m2sb=0D0
      m2sb(1:n,1:n)=-k1
      do i=1,n
        m2sb(n+i,n+i)=1D0
      enddo

      if(.not.dsigmade)then
        call geigenvalues( m2s,m2sb, walpha,wbeta, 2 * n,vr2,info)
        vl2=1D0
      else
        call geigenvalueslr( m2s,m2sb, walpha,wbeta, 2 * n,vr2,vl2, info)
      endif


!+++ the following lines should be commented out if H1^-1 is to be calculated
!+++        if(info.ne.0.and.svdtol.ne.0D0)then
!+++          call svdm3(k1,k1o,n,svdtol)
!+++          km1=DCONJG(transpose(k1))
!+++          usehinv1=.true.
!+++          goto 1100
!+++        endif
!+++ end H1^-1

      do i1=1,2 * n
        if(ABS(wbeta(i1)).eq.0D0)then
          zv(i1)=1D20
        elseif(ABS(walpha(i1)).eq.0D0)then
          zv(i1)=1D-20
        else
          zv(i1)=walpha(i1)/wbeta(i1)
        endif
      enddo

    endif



    do i1=1,2 * n

      normvr1=SQRT(DOT_PRODUCT(vr2(1:n,i1),vr2(1:n,i1)))
      normvr2=SQRT(DOT_PRODUCT(vr2(n+1:2*n,i1),vr2(n+1:2*n,i1)))

      if(normvr1.ge.normvr2)then
        vr(:,i1)=vr2(1:n,i1)
        normvr=normvr1
      else
        vr(:,i1)=vr2(n+1:2 * n,i1)
        normvr=normvr2
      endif
      vr(:,i1)=vr(:,i1)/normvr

      maxv=0.0D0
      do i=1,n
        if(ABS(vr(i,i1)).gt.abs(maxv))maxv=vr(i,i1)
      enddo
      vr(:,i1)=vr(:,i1)/maxv

      vl(:,i1)=vl2(1:n,i1)
      normvr=SQRT(DOT_PRODUCT(vl(:,i1),vl(:,i1)))
      vl(:,i1)=vl(:,i1)/normvr
    enddo


  end subroutine solveqep


  subroutine splitlgrg(n,nrchan,zv,k1,km1,k1o,signs,vml,vmr, ene,nrs,nls,side,svdtol,vr,iroc,iloc,tolki)


! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      
    IMPLICIT NONE

    integer n,ir,il,iroc,iloc,nrchan,i1,i2,nrs,nls,itol,vmr(n),vml(n)
    double precision tolkil(8),tolki,alphak,lz,svdtol,signs
    double complex vg,ene,zv(2 * n),k1(n,n),km1(n,n),k1o(n,n),vri(n), vrbuf(n), vr(n,2*n)
    CHARACTER(LEN=1) :: side
     

    tolkil(1)=1d-6
    tolkil(2)=1d-5
    tolkil(3)=2d-4
    tolkil(4)=1d-7
    tolkil(5)=1d-8
    tolkil(6)=1d-9
    tolkil(7)=1d-10
    tolkil(8)=1d-3

    do itol=1,8

      ir=1
      il=1
      iroc=0
      iloc=0
      nrchan=0
      tolki=tolkil(itol)

      do i1=1,2 * n
        call alphaofz(alphak,lz,zv(i1))
        if(ABS(lz).lt.tolki)then
!          write(*,*)"open channel:",zv(i1)
          k1o=signs * (k1 * zv(i1) - km1 / zv(i1))
          vri=vr(:,i1)
          vrbuf=matmul(k1o,vri)
          vg=0D0

          do i2=1,n
            vg=vg+(0D0,1D0) *  DCONJG(vri(i2)) * vrbuf(i2)
          enddo

          if(DREAL(vg).gt.0D0)then
!            write(12347,*)"right going vg=",DREAL(vg),dimag(vg)
            if(ir.gt.n)exit
            vmr(ir)=i1
            ir=ir+1
            iroc=iroc+1
            nrchan=nrchan+1
          else
!            write(12347,*)"left going vg=",DREAL(vg),dimag(vg)
            if(il.gt.n)exit
            vml(il)=i1
            il=il+1
            iloc=iloc+1
          endif
        else

!!!          k1o=signs * (k1 * zv(i1) - km1 / zv(i1))
!!!          vri=vr(:,i1)
!!!          vrbuf=matmul(k1o,vri)
!!!          vg=0D0
!!!
!!!          do i2=1,n
!!!            vg=vg+(0D0,1D0) *  DCONJG(vri(i2)) * vrbuf(i2)
!!!          enddo

!          write(*,*)"closed channel:",zv(i1)
          if(signs * lz.lt.0D0)then
!!!            write(12347,*)"right decaying vg=",DREAL(vg),dimag(vg)
!            write(*,*)"right decaying state"
            if(ir.gt.n)exit
            vmr(ir)=i1
            ir=ir+1
          else
!!!            write(12347,*)"left decaying vg=",DREAL(vg),dimag(vg)
!            write(*,*)"left decaying state"
            if(il.gt.n)exit
            vml(il)=i1
            il=il+1
          endif
        endif

      enddo

      nrs=ir-1
      nls=il-1

      if(iloc.eq.iroc.and.nrs.eq.nls)then
        exit
!      elseif(iloc.ne.iroc)then
!        write(12347,*)"warning:changing tolki; npropagating-left not equal npropagating-right=", DREAL(ene),nrs,nls,iroc,iloc,side,svdtol,tolki
!      elseif(nrs.ne.nls)then
!        write(12347,*)"warning:changing tolki; nrs,nls=", DREAL(ene),nrs,nls,iroc,iloc,side,svdtol,tolki
      endif

    enddo

    if(iloc.ne.iroc)then
      write(12347,*)"warning: npropagating-left not equal npropagating-right=", DREAL(ene),nrs,nls,iroc,iloc,side,svdtol,tolki
    elseif(nrs.ne.nls)then
      write(12347,*)"warning: nchannels-left not equal nchannels-right=", DREAL(ene),nrs,nls,iroc,iloc,side,svdtol,tolki
    endif

  end subroutine splitlgrg


  subroutine checkdsigma(rcondr,rcondl,rcondgf,svdtol, dsigma,dsigmar,dsigmamax,dsigmarm,dsigmam,svdtolm, n,k1o,sigma,anorml,anormr,ene,svdtolmax, anormgf,calcsucceed,side,gf,nrchanm,vrgout,vrgbout, zvrout, nrchan,tolki,vlg,vlgb,zvl,vlev, vrg,vrgb,zvr,  matout,grfunct, info,ik,dsigmade,vlevout)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    implicit none
    
    integer, intent(in) :: ik
    integer info,n,nrchanm,nrchan
    double precision relcond,rcondr,rcondl,rcondgf,dsigmartrue,dsigmar,svdtol,dsigma,dsigmamax,anorml,anormr,svdtolmax, anormgf,tolki,dsigmarm,dsigmam,svdtolm
    double complex ene,k1o(n,n),sigma(n,n),matout(n,n),grfunct(n,n), vrgout(n,n),vlg(n,n),vrgbout(n,n),vlgb(n,n),zvrout(n), zvl(n),vlev(n,n), zvr(n),vrgb(n,n),vrg(n,n)
    double complex vlevout(n,n)
    logical calcsucceed
    logical, intent(in) :: dsigmade
    CHARACTER(LEN=1) :: gf,side


    if(side.eq.'L')then
      relcond=rcondl/rcondgf
    else
      relcond=rcondr/rcondgf
    endif
    if(relcond.gt.1D0)then
      dsigmartrue=dsigmar / relcond
    else
      dsigmartrue=dsigmar
    endif

    if(.false.)then
      write(12347,*) "dsi=",side,dreal(ene),svdtol,dsigma,dsigmar,dsigmamax, dsigmartrue,relcond, MAXVAL(ABS(k1o)),MAXVAL(ABS(sigma)),anorml,anormr,rcondl, rcondr,anormgf,rcondgf, ik
   
      if(relcond.gt.1D0)then
        write(12347,*) "dsiss=",side,dreal(ene),svdtol,dsigma,dsigmar,dsigmamax, dsigmartrue,relcond, MAXVAL(ABS(k1o)),MAXVAL(ABS(sigma)),anorml,anormr,rcondl, rcondr,anormgf,rcondgf, ik
      endif
    endif

    if(dsigmar.lt.dsigmamax.or.dsigmar.lt.dsigmarm.or.dsigmarm.eq.0D0) then 
      calcsucceed=.true.
      nrchanm=nrchan
      if( gf .eq. 'G')then
        matout=grfunct
      else
        matout=sigma
      endif
      if(side.eq.'L')then
        vrgout=vlg
        vrgbout=vlgb
        zvrout=zvl
        if(dsigmade)vlevout=vlev
      else
        vrgout=vrg
        vrgbout=vrgb
        zvrout=zvr
        if(dsigmade)vlevout=vlev
      endif
      if(dsigmar.lt.dsigmamax)then
        dsigmam=dsigma
        dsigmarm=dsigmar
        svdtolm=svdtol
      endif
    endif

    info=0
    if(dsigmar.gt.dsigmamax.and.svdtol.lt.svdtolmax)then 
      write(12347,*)"dsir=",side,dreal(ene),svdtol,dsigma,dsigmar, MAXVAL(ABS(k1o)),MAXVAL(ABS(sigma)), ik
      if(dsigmar.lt.dsigmarm.or.dsigmarm.eq.0D0)then
        dsigmam=dsigma
        dsigmarm=dsigmar
        svdtolm=svdtol
      endif
      info=1
    elseif(svdtol.ne.0D0)then
      write(12347,*) "dsirl=",side,dreal(ene),svdtol,dsigma,dsigmar, MAXVAL(ABS(k1o)),MAXVAL(ABS(sigma)), ik
      if(dsigmar.lt.dsigmarm.or.dsigmarm.eq.0D0)then
        dsigmam=dsigma
        dsigmarm=dsigmar
        svdtolm=svdtol
      endif
    endif

     

  end subroutine checkdsigma


  subroutine check_error_sigma(dsigma,write_error,energ_sigma,side,n,H0_L,H1_L,S0_L,S1_L,sigma)

    use mConstants

    CHARACTER(LEN=1),intent(in) :: side
    logical,intent(in) :: write_error
    integer, intent(in)     :: n
    complex(kdp),intent(in) :: sigma(n,n),h0_L(n,n),h1_L(n,n),s0_L(n,n),s1_L(n,n)
    complex(kdp),intent(in) :: energ_sigma
    real(kdp),intent(out) :: dsigma
 
    complex(kdp) :: grf(n,n),k1(n,n),k0(n,n),km1(n,n)
    INTEGER , ALLOCATABLE :: ilu(:)
    complex(kdp), ALLOCATABLE :: work(:)
    integer nwork,info

    k0=energ_sigma * s0_L - h0_L
    k1=energ_sigma * s1_L - h1_L
    km1=energ_sigma * DCONJG(TRANSPOSE(s1_L)) - DCONJG(TRANSPOSE(h1_L))

    grf=k0-sigma

    nwork=8
    allocate(work(nwork*n))
    allocate(ilu(n))
    call ZGETRF( n, n, grf, n, ilu, INFO )
    call ZGETRI( n, grf,n, ilu, WORK, nwork * n, INFO )
    deallocate(ilu,work)
    if(side.eq.'L')then
      grf=matmul(km1,matmul(grf,k1))
    else
      grf=matmul(k1,matmul(grf,km1))
    endif

    grf=sigma-grf
    dsigma=MAXVAL(ABS(grf))

    if(write_error) write(12347,*)"dsigma=",DREAL(energ_sigma),DIMAG(energ_sigma),side,dsigma

  end subroutine check_error_sigma


  subroutine printcomplexbands(side,n,zvr,zvl,ene,svdtol,ik)


! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************


    use negfmod, only: ef_l,kpointmod

    IMPLICIT none

    integer n,il
    integer, intent(in) :: ik
    double complex zvr(n),zvl(n),zvu,ene,ener
    double precision alphak,lz,svdtol
    CHARACTER(LEN=1) :: side

    ener=(ene-ef_l) * 13.6057D0

    do il=1,n
      If (side .eq. 'R') then
        zvu=1D0/zvr(il)
      else
        zvu=zvr(il)
      endif
      call alphaofz(alphak,lz,zvu)
      If (side .eq. 'R') then
        write(12346,*)"k_r_rg=",DREAL(ener),DIMAG(ener),-lz,alphak, svdtol,ik,kpointmod(1),kpointmod(2)

      else
        write(12346,*)"k_l_rg=",DREAL(ener),DIMAG(ener),-lz,alphak, svdtol,ik,kpointmod(1),kpointmod(2)
      endif
      If (side .eq. 'R') then
        zvu=1D0/zvl(il)
      else
        zvu=zvl(il)
      endif
      call alphaofz(alphak,lz,zvu)
      If (side .eq. 'L') then
        write(12346,*)"k_l_lg=",DREAL(ener),DIMAG(ener),-lz,alphak, svdtol,ik,kpointmod(1),kpointmod(2)
      else
        write(12346,*)"k_r_lg=",DREAL(ener),DIMAG(ener),-lz,alphak, svdtol,ik,kpointmod(1),kpointmod(2)
      endif
    enddo

  end subroutine printcomplexbands


  subroutine printcomplexbandsvec(side,n,zvr,zvl,ene,svdtol,vrg,ik)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************


    use negfmod, only: ef_l,kpointmod

    IMPLICIT none

    integer n,il,i1
    integer, intent(in) :: ik
    double complex zvr(n),zvl(n),zvu,ene,ener,vrg(2*n,2*n)
    double precision alphak,lz,svdtol,tol,maxkappa
    CHARACTER(LEN=1) :: side

    maxkappa=8d0

    tol=1d-1
    ener=(ene-ef_l) * 13.6057D0

    do il=1,n
      If (side .eq. 'R') then
        zvu=1D0/zvr(il)
      else
        zvu=zvr(il)
      endif
      call alphaofz(alphak,lz,zvu)
      If (side .eq. 'R') then
        do i1=1,2 * n
          if((abs(vrg(i1,il)).ge.tol).and.(abs(lz).le.maxkappa))then
          write(12347,*)"k_r_rg=",DREAL(ener),DIMAG(ener),-lz,alphak, i1,il,abs(vrg(i1,il)), svdtol,ik,kpointmod(1),kpointmod(2)
          endif
        enddo
      else
        do i1=1,2 * n
          if((abs(vrg(i1,il)).ge.tol).and.(abs(lz).le.maxkappa))then
          write(12347,*)"k_l_rg=",DREAL(ener),DIMAG(ener),-lz,alphak, i1,il,abs(vrg(i1,il)), svdtol,ik,kpointmod(1),kpointmod(2)
          endif
        enddo
      endif
      If (side .eq. 'R') then
        zvu=1D0/zvl(il)
      else
        zvu=zvl(il)
      endif
      call alphaofz(alphak,lz,zvu)
      If (side .eq. 'L') then
        do i1=1,2 * n
          if((abs(vrg(i1,il)).ge.tol).and.(abs(lz).le.maxkappa))then
          write(12347,*)"k_l_lg=",DREAL(ener),DIMAG(ener),-lz,alphak,i1,il,abs(vrg(i1,il)),svdtol,ik,kpointmod(1),kpointmod(2)
          endif
        enddo
      else
        do i1=1,2 * n
          if((abs(vrg(i1,il)).ge.tol).and.(abs(lz).le.maxkappa))then
          write(12347,*)"k_r_lg=",DREAL(ener),DIMAG(ener),-lz,alphak,i1,il,abs(vrg(i1,il)),svdtol,ik,kpointmod(1),kpointmod(2)
          endif
        enddo
      endif
    enddo

  end subroutine printcomplexbandsvec


  subroutine sigmavinv(side,n,ar,al,vrg,vrgb,vlg,vlgb,zvr,zvl, grfunct,sigma,k1,km1,k1o,k1ob,ipiv,worki,signs)


! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    IMPLICIT none

    integer i1,n,il,i2,i,info,ipiv(n)
    double complex zvr(n),zvl(n),zvu
    double complex k1(n,n),km1(n,n),k1o(n,n),k1ob(n,n),vri(n),vli(n), vritd(n),vlitd(n),vlg(n,n),vrg(n,n),vlgb(n,n),vrgb(n,n), grfunct(n,n),sigma(n,n), ar(n,n),al(n,n),worki(n**2)
    double precision signs
    CHARACTER(LEN=1) :: side



    if(side.eq.'L')then
 
        ar=0D0
        do il=1,n
          vri(:)=vrg(:,il)
          vritd(:)=vrgb(il,:)
          do i1=1,n
            do i2=1,n
              ar(i1,i2)=ar(i1,i2)+ vri(i1) * vritd(i2) * zvr(il)
            enddo
          enddo
        enddo
 
        k1ob=-matmul(al,ar)
        do i=1,n
          k1ob(i,i)=1D0+k1ob(i,i)
        enddo

      else

        al=0D0
        do il=1,n
          zvu=1D0/zvl(il)
          vli(:)=vlg(:,il)
          vlitd(:)=vlgb(il,:)
          do i1=1,n
            do i2=1,n
              al(i1,i2)=al(i1,i2)+ vli(i1) * vlitd(i2) / zvu
            enddo
          enddo
        enddo

        k1ob=-matmul(ar,al)
        do i=1,n
          k1ob(i,i)=1D0+k1ob(i,i)
        enddo

      endif


      if(side.eq.'L')then

        ar=0D0
        do il=1,n
          vri(:)=vrg(:,il)
          vritd(:)=vrgb(il,:)
          do i1=1,n
            do i2=1,n
              ar(i1,i2)=ar(i1,i2)+ vri(i1) * vritd(i2) / zvr(il)
            enddo
          enddo
        enddo

      else

        al=0D0
        do il=1,n
          zvu=1D0/zvl(il)
          vli(:)=vlg(:,il)
          vlitd(:)=vlgb(il,:)
          do i1=1,n
            do i2=1,n
              al(i1,i2)=al(i1,i2)+ vli(i1) * vlitd(i2) * zvu
            enddo
          enddo
        enddo

      endif

      k1o=matmul(km1,signs * (ar-al))
      CALL ZGETRF(n,n,k1o,n,IPIV,INFO)
      CALL ZGETRI(n,k1o,n,IPIV,worki,n**2,INFO)

      grfunct=matmul(k1ob,k1o)

      k1o=matmul(grfunct,k1)
      sigma=matmul(km1,k1o)



  end subroutine sigmavinv


  subroutine evsigma(nl,sigmal,nr,sigmar,ispin,ik,EiCompL)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
 
    implicit none
    
    integer i1,nl,nr,ispin,ik
    DOUBLE COMPLEX :: zvl(NL),zvr(NR)
    double complex sigmal(nl,nl),sigmar(nr,nr),EiCompL,EiCompR
    double precision  dsigma,dsigmamin

    call geigenvalues1( SigmaL,zvl, NL)

    do i1=1, NL
      if(ISPIN.eq.1)then
        write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sl1=",ik,DREAL(EiCompL),  DREAL(zvl(i1)),DIMAG(zvl(i1))
      else
        write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sl2=",ik,DREAL(EiCompL),  DREAL(zvl(i1)),DIMAG(zvl(i1))
      endif
    enddo
    dsigma=MAXVAL(ABS(zvl))
    dsigmamin=MINVAL(ABS(zvl))
    if(ISPIN.eq.1)then
      write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sl1M=",ik,DREAL(EiCompL),dsigma, dsigmamin
    else
      write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sl2M=",ik,DREAL(EiCompL),dsigma, dsigmamin
    endif



    call geigenvalues1( SigmaR, zvr, NR)

    do i1=1, NR
      if(ISPIN.eq.1)then
        write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sr1=",ik,DREAL(EiCompR), DREAL(zvr(i1)),DIMAG(zvr(i1))
      else
        write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sr2=",ik,DREAL(EiCompL), DREAL(zvr(i1)),DIMAG(zvr(i1))
      endif
    enddo
    dsigma=MAXVAL(ABS(zvr))
    dsigmamin=MINVAL(ABS(zvr))
    if(ISPIN.eq.1)then
      write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sr2=",ik,DREAL(EiCompL),"sr1M=",ik,DREAL(EiCompL),dsigma,dsigmamin
    else
      write(12347,'(a4," " ,i9," " ,e12.5," " ,e12.5," " ,e12.5)')"sr2M=",ik,DREAL(EiCompL),dsigma, dsigmamin
    endif


  end subroutine evsigma


  subroutine dostrc(n,sigmal,sigmar,h0,h1,s0,s1,ei, ispin,dos,ef,leadspdos)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
 
    implicit none
    
    integer n,info,ispin,ii
    INTEGER, DIMENSION (n) :: IPIV
    double complex sigmal(n,n),sigmar(n,n),h0(n,n),h1(n,n), s0(n,n),s1(n,n)
    double complex k0(n,n),k1(n,n),km1(n,n),k1o(n,n),k1ob(n,n)
    double complex grf(n,n),g01(n,n),g10(n,n)
    DOUBLE COMPLEX, DIMENSION (n,n) :: pdos2,pdos3
    double complex work(n,n)
    double complex ei
    double precision dosk2,dosk3, dos,signspin,totpdos,ef,pi
    logical  leadspdos

    pi=3.141592653589793D0

    k0=h0 - s0 * ei
    k1=h1 - s1 * ei
    km1=DCONJG(TRANSPOSE(h1)) - DCONJG(TRANSPOSE(s1)) * ei


    k1o=-k0-sigmal-sigmar
    CALL ZGETRF(n,n,k1o,n,IPIV,info)
    CALL ZGETRI(n,k1o,n,IPIV,work,n**2,info)
    grf=k1o

    k1o=-k0-sigmal
    CALL ZGETRF(n,n,k1o,n,IPIV,info)
    CALL ZGETRI(n,k1o,n,IPIV,work,n**2,info)
    g01=matmul(matmul(k1o,k1),grf)
 
    k1o=-k0-sigmar
    CALL ZGETRF(n,n,k1o,n,IPIV,info)
    CALL ZGETRI(n,k1o,n,IPIV,work,n**2,info)
    g10=matmul(matmul(k1o,km1),grf)

    k1o=((0D0,0.5D0) /pi) * (grf-DCONJG(TRANSPOSE(grf)))
    k1ob=matmul(k1o,s0)
    dosk2=0D0
    if(leadspdos)pdos2=k1ob
    do ii=1,n
      dosk2=dosk2+k1ob(ii,ii)
    enddo

    k1o=((0D0,0.5D0) /pi) * (g10-DCONJG(TRANSPOSE(g01)))
    k1ob=matmul(k1o,s1)

    dosk3=0D0
    if(leadspdos)pdos3=k1ob
    do ii=1,n
      dosk3=dosk3+DREAL(k1ob(ii,ii))
    enddo

    k1o=((0D0,0.5D0) /pi) * (g01-DCONJG(TRANSPOSE(g10)))
    k1ob=matmul(k1o,DCONJG(TRANSPOSE(s1)))

    if(leadspdos)pdos3=pdos3+k1ob
    do ii=1,n
      dosk3=dosk3+DREAL(k1ob(ii,ii))
    enddo

    dos = (dosk2+dosk3) / 13.6057D0
    if(leadspdos)pdos2=pdos2+pdos3

    if(leadspdos)then
      if(ISPIN.eq.1)then
        signspin=1D0
      else
        signspin=-1D0
      endif

      totpdos=0D0
      do ii=1,n
        write(12346,*)"bulkdos_sr=",ii,ispin, (ei-ef)* 13.6057D0, signspin * DREAL(pdos2(ii,ii))/13.657D0
        totpdos=totpdos+DREAL(pdos2(ii,ii))
      enddo
      write(12346,*)"totbulkdos_sr=",ispin,(ei-ef)* 13.6057D0, signspin * totpdos/13.6057D0
      
    endif

  end subroutine dostrc


  subroutine vg_dos_dosvv(n,ene,zvr,vrg,S0,S1,km1,k1,tolki)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
 
    use negfmod

    IMPLICIT NONE

    INTEGER n
    DOUBLE COMPLEX, DIMENSION (n,n) :: k1,km1,S0,S1,k1t,vrg,s01k
    DOUBLE COMPLEX zvr(n),vri(n),vrbuf(n),vg,ene
    DOUBLE PRECISION dos,dosvv,alphak,lz,pi
    double precision, intent(in) :: tolki

    INTEGER   i,j

    pi=3.141592653589793D0

    dos=0D0
    dosvv=0D0
    do i=1,n
      call alphaofz(alphak,lz,zvr(i))
      if(ABS(lz).lt.tolki)then
        k1t=(k1 * zvr(i) - km1 / zvr(i))
        vri=vrg(:,i)
        
        s01k=s0+s1 * zvr(i) + DCONJG(TRANSPOSE(s1))/ zvr(i)
        vrbuf=matmul(s01k,vri)
        vg=0D0
        do j=1,n
          vg=vg+ DCONJG(vri(j)) * vrbuf(j)
        enddo
        vg=sqrt(vg)
        vri=vri/vg

        vrbuf=matmul(k1t,vri)
        vg=0D0
        do j=1,n
          vg=vg+(0D0,1D0) *  DCONJG(vri(j)) * vrbuf(j)
        enddo
        vg=vg* 13.6057d0 ! this transforms the units of the group velocity from Ry to eV

        dos=dos+1D0/(abs(vg) * pi)
        dosvv=dosvv+abs(vg)/ pi
      endif

    enddo

    dosk=dos
    doskvv=dosvv

  end subroutine vg_dos_dosvv


  subroutine sdsigmade(n,ene,zvr,vrg,vrgt,S0,S1,k0,k1,km1,sigma, side,vlevout)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
 

    use negfmod

    IMPLICIT NONE

    INTEGER n
    DOUBLE COMPLEX, DIMENSION (n,n) :: k0,k1,km1,S0,S1,k1t,vrg,s01k, dmnde,psipsit,mn,xn,dppde,sigma,vrgt,dppdes,dtde,tbar,  mdsigmade
    DOUBLE COMPLEX zvr(n),vri(n),vritd(n),vrbuf(n),vg(n),ene,  work(4 * n),trxn,vrild(n),ascale,vgs,vgsev
    DOUBLE PRECISION dos,dosvv,alphak,lz,pi,psimax
    double complex, intent(in) :: vlevout(n,n)

    INTEGER   i,j,i1,i2,npsimax(n),ipiv(n),INFO
    DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
    CHARACTER(LEN=1) :: side


    pi=3.141592653589793D0

    dos=0D0
    dosvv=0D0
    dppdes=0D0
    do i=1,n
        if(side.eq.'L')then
          k1t=(-k0 - 2D0 * km1 / zvr(i))
        else
          k1t=(k0 + 2D0 * k1 / zvr(i))
        endif
        vrild=DCONJG(vlevout(:,i))
        vri=vrg(:,i)
        vritd=vrgt(i,:)

        if(side.eq.'L')then
          s01k=s0+s1 * zvr(i) + DCONJG(TRANSPOSE(s1))/ zvr(i)
        else
          s01k=s0+s1 / zvr(i) + DCONJG(TRANSPOSE(s1))* zvr(i)
        endif
        vrbuf=matmul(s01k,vri)
        vgs=0D0
        do j=1,n
          vgs=vgs+ vrild(j) * vrbuf(j)
        enddo
        vgs=sqrt(vgs)
        vri=vri/vgs
        vrild=vrild/vgs
!        write(12347,*)"vrvla1=",DREAL(ene),DREAL(vri(1)),
!   .        DREAL(vrild(1)), DIMAG(vri(1)),DIMAG(vrild(1))
!        write(12347,*)"vrvla2=",DREAL(ene),DREAL(vri(2)),
!   .        DREAL(vrild(2)), DIMAG(vri(2)),DIMAG(vrild(2))
        vritd=vritd*vgs

        vrbuf=matmul(k1t,vri)
        vgs=0D0
        do j=1,n
          vgs=vgs+(0D0,1D0) *  vrild(j) * vrbuf(j)
        enddo
        vgsev=vgs* 13.6057d0 ! this transforms the units of the group velocity from Ry to eV
        write(12347,*)"vgds=",DREAL(ene),DREAL(vgsev),DIMAG(vgsev), alphak,lz,vgsev**2

        vg(i)=vgs

        dmnde=-s01k+zi * k1t / vg(i)
!        dmnde=dmnde / zvr(i)
        if(side.eq.'L')then
          mn=k0+k1 * zvr(i) + km1 / zvr(i)
        else
          mn=k0+k1 / zvr(i) + km1 * zvr(i)
        endif

        write(12347,*)"dmn=",DREAL(ene),DREAL(dmnde(1,1)), DIMAG(dmnde(1,1)),DREAL(mn(1,1)), DIMAG(mn(1,1))


        psimax=0D0
        npsimax(i)=1
        do i1=1,n
          if(ABS(vri(i1)).gt.psimax)then
            psimax=ABS(vri(i1))
            npsimax(i)=i1
          endif
        enddo
!        npsimax=1
        write(12347,*)"npsimax=",npsimax(i)





        ascale=1D0
        vritd=vritd*(ascale * vri(npsimax(i)))
        vri=vri/(ascale * vri(npsimax(i)))

        do i1=1,n
          do i2=1,n
            psipsit(i1,i2)=vri(i1) * vritd(i2)
          enddo
        enddo

        trxn=0D0
        do i1=1,n
          trxn=trxn+psipsit(i1,i1)
        enddo
        write(12347,*)"trpsipsit=",trxn

!        dmnde=-MATMUL(dmnde,psipsit)

        mn(npsimax(i),:)=0D0
        mn(:,npsimax(i))=0D0
        mn(npsimax(i),npsimax(i))=1D0
        dmnde(npsimax(i),:)=0D0

        call ZGETRF( n, n, mn, n, ipiv, INFO )
        call ZGETRI( n, mn,n, ipiv, WORK, 4 * n, INFO )

        xn(:,i)=-matmul(matmul(mn,dmnde),vri)

        write(12347,*)"dx=",DREAL(ene),DREAL(xn(1,i)), DIMAG(xn(1,i)),DREAL(xn(2,i)),DIMAG(xn(2,i)), DREAL(vri(1)),DIMAG(vri(1)), DREAL(vri(2)),DIMAG(vri(2))

        do i1=1,n
          do i2=1,n
            dppdes(i1,i2)=dppdes(i1,i2)+xn(i1,i) * vritd(i2)
          enddo
        enddo


!        deikde=zi * zvr(i) / vgs

!        write(12347,*)"deikde=",DREAL(ene),DREAL(deikde),
!   .        DIMAG(deikde),DREAL(zvr(i)),DIMAG(zvr(i)),i

    enddo

    dtde=0D0
    tbar=0D0
    do i=1,n
      vri=vrg(:,i)
      vritd=vrgt(i,:)

      vritd=vritd*(ascale * vri(npsimax(i)))
      vri=vri/(ascale * vri(npsimax(i)))

      do i1=1,n
        do i2=1,n
          psipsit(i1,i2)=vri(i1) * vritd(i2)
        enddo
      enddo

      dppde=0D0
      do i1=1,n
        do i2=1,n
          dppde(i1,i2)=xn(i1,i) * vritd(i2)
        enddo
      enddo

      dppde=dppde-matmul(psipsit,dppdes)

      write(12347,*)"dppde=",DREAL(ene),DREAL(dppde(1,1)), DIMAG(dppde(1,1)),DREAL(psipsit(1,1)), DIMAG(psipsit(1,1))


      if(side.eq.'L')then

        dtde=dtde+(-zi /(zvr(i) * vg(i))) * psipsit + (1D0/zvr(i)) * dppde

        tbar=tbar+ (1D0/zvr(i)) * psipsit 


      else
        dtde=dtde+(zi /(zvr(i) * vg(i))) * psipsit + (1D0/zvr(i)) * dppde

        tbar=tbar+ (1D0/zvr(i)) * psipsit 
      endif


    enddo
    
    if(side.eq.'L')then
      mdsigmade= - matmul(DCONJG(TRANSPOSE(s1)),tbar) + matmul(km1,dtde)
      write(12347,*)"dsigmade=",DREAL(ene), DREAL(mdsigmade(n,n)),DIMAG(mdsigmade(n,n)), DREAL(sigma(n,n)), DIMAG(sigma(n,n)),1
    else
      mdsigmade= - matmul(s1,tbar) + matmul(k1,dtde)
      write(12347,*)"dsigmade=",DREAL(ene), DREAL(mdsigmade(n,n)),DIMAG(mdsigmade(n,n)), DREAL(sigma(n,n)), DIMAG(sigma(n,n)),2
    endif

    


  end subroutine sdsigmade

 

end module mSigmaMethod1

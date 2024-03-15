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
!                   ROTATEPARTIALDM2,
!                   ROTATEPARTIALDM3,
!                   ROTATEPARTIALDM,
!                   ROTATEDM  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine rotatePartialDM2(dm,maxnh,numh,listhptr,listh,NspinRealInputMatrix,no,N1,n1loc,indxuo,thetadeg,phideg,istart,iend,mynode,nnodes)

  use mConstants
  use mMPI_NEGF
#ifdef MPI
  use parallel
#endif


  implicit none

  integer, intent(in) :: istart,iend
  integer,intent(in) :: maxnh,NspinRealInputMatrix,n1,n1loc
  real(kdp), intent(inout) :: dm(maxnh,NspinRealInputMatrix)
  real(kdp), intent(in) :: thetadeg,phideg
  integer,intent(in):: no,indxuo(no),numh(n1loc),listhptr(n1loc),listh(maxnh),mynode,nnodes
!  real(kdp) :: dmorig(maxnh,NspinRealInputMatrix)

  real(kdp) :: phi,theta,mtot(3),dm3new,dm4new
  real(kdp) :: PI
  integer :: is
  real(kdp) :: dm0xyz(maxnh,4)
  real(kdp) :: rotMat(3,3),buf1,buf2,buf3,buf4


  integer iluo,ind,jguo,j,iguo

  PI=2.0_kdp * acos(0.0_kdp)

  phi=phideg * PI/180.0_kdp
  theta=thetadeg * PI/180.0_kdp

  rotmat=0.0_kdp
  rotmat(1,1)=cos(theta)
  rotmat(1,3)=sin(theta)
  rotmat(2,2)=1.0_kdp
  rotmat(3,1)=-sin(theta)
  rotmat(3,3)=cos(theta)

  dm0xyz(:,1)=0.5_kdp * (dm(:,1)+dm(:,2))
  if(NspinRealInputMatrix.gt.2)then
    dm0xyz(:,2)=dm(:,3)
    dm0xyz(:,3)=-dm(:,4)
  else
    dm0xyz(:,2)=0.0D0
    dm0xyz(:,3)=0.0D0
  endif
  dm0xyz(:,4)=0.5_kdp * (dm(:,1)-dm(:,2))

!  dmorig=dm

  do iluo = 1,n1loc
#ifdef MPI
    call LocalToGlobalOrb(iluo,MyNode,NNodes,iguo)
#else
    iguo=iluo
#endif
    do j = 1,numh(iluo)
      ind = listhptr(iluo) + j
      jguo = indxuo(listh(ind))
!      write(*,*)"indi",iguo,jguo,istart,iend,mynode
      if((iguo>=istart .and. iguo <= iend).and.&
&         (jguo>=istart .and. jguo <= iend))then
        if(NspinRealInputMatrix.gt.2)then
          buf1=rotmat(1,1) * dm0xyz(ind,2)+rotmat(1,2) * dm0xyz(ind,3)+rotmat(1,3) * dm0xyz(ind,4)
          buf2=rotmat(2,1) * dm0xyz(ind,2)+rotmat(2,2) * dm0xyz(ind,3)+rotmat(2,3) * dm0xyz(ind,4)
        endif
        buf3=rotmat(3,1) * dm0xyz(ind,2)+rotmat(3,2) * dm0xyz(ind,3)+rotmat(3,3) * dm0xyz(ind,4)

!        buf1=dm0xyz(ind,2)
!        buf2=dm0xyz(ind,3)
!        buf3=dm0xyz(ind,4)

        dm(ind,1)=dm0xyz(ind,1)+buf3
        dm(ind,2)=dm0xyz(ind,1)-buf3
        if(NspinRealInputMatrix.gt.2)then
          dm(ind,3)=buf1
          dm(ind,4)=-buf2
        endif


!        write(12345,*)"dmmdmorig=",ind,iguo,jguo,dm(ind,1)-dmorig(ind,1),dm(ind,2)-dmorig(ind,2),dm(ind,3)-dmorig(ind,3),dm(ind,4)-dmorig(ind,4)
!        write(12345,*)"dm1=",ind,iguo,jguo,dm(ind,1),dmorig(ind,1),dm(ind,1)-dmorig(ind,1)
!        write(12345,*)"dm2=",ind,iguo,jguo,dm(ind,2),dmorig(ind,2),dm(ind,2)-dmorig(ind,2)
!        write(12345,*)"dm3=",ind,iguo,jguo,dm(ind,3),dmorig(ind,3),dm(ind,3)-dmorig(ind,3)
!        write(12345,*)"dm4=",ind,iguo,jguo,dm(ind,4),dmorig(ind,4),dm(ind,4)-dmorig(ind,4)


      endif
       
    enddo
  enddo

end subroutine rotatePartialDM2


subroutine rotatePartialDM3(dm,maxnh,numh,listhptr,listh,NspinRealInputMatrix,no,N1,n1loc,indxuo,thetadeg,phideg,istart,iend,mynode,nnodes)

  use mConstants
  use mMPI_NEGF
  use negfmod, only :outinfo
#ifdef MPI
  use parallel
#endif


  implicit none

  integer, intent(in) :: istart,iend
  integer,intent(in) :: maxnh,NspinRealInputMatrix,n1,n1loc
  real(kdp), intent(inout) :: dm(maxnh,NspinRealInputMatrix)
  real(kdp), intent(in) :: thetadeg,phideg
  integer,intent(in):: no,indxuo(no),numh(n1loc),listhptr(n1loc),listh(maxnh),mynode,nnodes


  real(kdp) :: dmorig(maxnh,NspinRealInputMatrix)
  real(kdp) :: phi,theta,mtot(3),dm3new,dm4new
  real(kdp) :: PI
  integer :: is


  integer iluo,ind,jguo,j,iguo

  PI=2.0_kdp * acos(0.0_kdp)

  phi=phideg * PI/180.0_kdp
  theta=thetadeg * PI/180.0_kdp
!  theta=0.0_kdp * PI

  mtot(3)=1.0_kdp
  mtot(1)=sin(theta) * cos(phi) * mtot(3)
  mtot(2)=sin(theta) * sin(phi) * mtot(3)
  mtot(3)=cos(theta) * mtot(3)
!  write(*,*)"angles=",sin(theta),cos(theta),sin(phi),cos(phi)
!  write(*,*)"mxyz=",mtot(1),mtot(2),mtot(3)
  dmorig=dm

  do iluo = 1,n1loc
#ifdef MPI
    call LocalToGlobalOrb(iluo,MyNode,NNodes,iguo)
#else
    iguo=iluo
#endif

    do j = 1,numh(iluo)
      ind = listhptr(iluo) + j
      jguo = indxuo(listh(ind))
!      write(*,*)"indi",iguo,jguo,istart,iend,mynode
      if((iguo>=istart .and. iguo <= iend).and.&
&         (jguo>=istart .and. jguo <= iend))then
        dm(ind,3)=0.5_kdp * (dm(ind,1)-dm(ind,2))
        dm(ind,4)=0.5_kdp * (dm(ind,1)+dm(ind,2))
        dm(ind,1)=dm(ind,4)+dm(ind,3)
        dm(ind,2)=dm(ind,4)-dm(ind,3)
        dm(ind,4)=0.0_kdp
        dm(ind,3)=0.0_kdp
!      elseif((iguo>=istart .and. iguo <= iend).or.&
!&         (jguo>=istart .and. jguo <= iend))then
!        dm3new=0.5_kdp * (dm(ind,1)-dm(ind,2))
!        dm4new=0.5_kdp * (dm(ind,1)+dm(ind,2))
!        dm(ind,1)=0.5_kdp * (dm(ind,1)+dm4new+mtot(3)*dm3new)
!        dm(ind,2)=0.5_kdp * (dm(ind,2)+dm4new-mtot(3)*dm3new)
!        dm(ind,4)=0.5_kdp * (dm(ind,4)-mtot(2)*dm3new)
!        dm(ind,3)=0.5_kdp * (dm(ind,3)+mtot(1)*dm3new)
      if (outinfo) then
        write(12345,*)"dmmdmorig=",ind,iguo,jguo,dm(ind,1)-dmorig(ind,1),dm(ind,2)-dmorig(ind,2),dm(ind,3)-dmorig(ind,3),dm(ind,4)-dmorig(ind,4)
        write(12345,*)"dm1=",ind,iguo,jguo,dm(ind,1),dmorig(ind,1),dm(ind,1)-dmorig(ind,1)
        write(12345,*)"dm2=",ind,iguo,jguo,dm(ind,2),dmorig(ind,2),dm(ind,2)-dmorig(ind,2)
        write(12345,*)"dm3=",ind,iguo,jguo,dm(ind,3),dmorig(ind,3),dm(ind,3)-dmorig(ind,3)
        write(12345,*)"dm4=",ind,iguo,jguo,dm(ind,4),dmorig(ind,4),dm(ind,4)-dmorig(ind,4)
    endif
      endif
       
    enddo
  enddo

end subroutine rotatePartialDM3



subroutine rotatePartialDM(dm,maxnh,numh,listhptr,listh,NspinRealInputMatrix,no,N1,n1loc,indxuo,thetadeg,phideg,istart,iend,mynode,nnodes)

  use mConstants
  use mMPI_NEGF
#ifdef MPI
  use parallel
#endif


  implicit none

  integer, intent(in) :: istart,iend
  integer,intent(in) :: maxnh,NspinRealInputMatrix,n1,n1loc
  real(kdp), intent(inout) :: dm(maxnh,NspinRealInputMatrix)
  real(kdp), intent(in) :: thetadeg,phideg
  integer,intent(in):: no,indxuo(no),numh(n1loc),listhptr(n1loc),listh(maxnh),mynode,nnodes


  real(kdp) :: dmorig(maxnh,NspinRealInputMatrix)
  real(kdp) :: phi,theta,mtot(3),dm3new,dm4new
  real(kdp) :: PI
  integer :: is


  integer iluo,ind,jguo,j,iguo

  PI=2.0_kdp * acos(0.0_kdp)

  phi=phideg * PI/180.0_kdp
  theta=thetadeg * PI/180.0_kdp
!  theta=0.0_kdp * PI

  mtot(3)=1.0_kdp
  mtot(1)=sin(theta) * cos(phi) * mtot(3)
  mtot(2)=sin(theta) * sin(phi) * mtot(3)
  mtot(3)=cos(theta) * mtot(3)
!  write(*,*)"angles=",sin(theta),cos(theta),sin(phi),cos(phi)
!  write(*,*)"mxyz=",mtot(1),mtot(2),mtot(3)
  dmorig=dm

  do iluo = 1,n1loc
#ifdef MPI
    call LocalToGlobalOrb(iluo,MyNode,NNodes,iguo)
#else
    iguo=iluo
#endif
    do j = 1,numh(iluo)
      ind = listhptr(iluo) + j
      jguo = indxuo(listh(ind))
!      write(*,*)"indi",iguo,jguo,istart,iend,mynode
      if((iguo>=istart .and. iguo <= iend).and.&
&         (jguo>=istart .and. jguo <= iend))then
        dm(ind,3)=0.5_kdp * (dm(ind,1)-dm(ind,2))
        dm(ind,4)=0.5_kdp * (dm(ind,1)+dm(ind,2))
        dm(ind,1)=dm(ind,4)+mtot(3)*dm(ind,3)
        dm(ind,2)=dm(ind,4)-mtot(3)*dm(ind,3)
        dm(ind,4)=-mtot(2)*dm(ind,3)
        dm(ind,3)=mtot(1)*dm(ind,3)
!      elseif((iguo>=istart .and. iguo <= iend).or.&
!&         (jguo>=istart .and. jguo <= iend))then
!        dm3new=0.5_kdp * (dm(ind,1)-dm(ind,2))
!        dm4new=0.5_kdp * (dm(ind,1)+dm(ind,2))
!        dm(ind,1)=0.5_kdp * (dm(ind,1)+dm4new+mtot(3)*dm3new)
!        dm(ind,2)=0.5_kdp * (dm(ind,2)+dm4new-mtot(3)*dm3new)
!        dm(ind,4)=0.5_kdp * (dm(ind,4)-mtot(2)*dm3new)
!        dm(ind,3)=0.5_kdp * (dm(ind,3)+mtot(1)*dm3new)
        write(*,*)"dmmdmorig=",ind,iguo,jguo,dm(ind,1)-dmorig(ind,1),dm(ind,2)-dmorig(ind,2),dm(ind,3)-dmorig(ind,3),dm(ind,4)-dmorig(ind,4)
      endif
       
    enddo
  enddo

end subroutine rotatePartialDM

subroutine rotateDM(dm,maxnh,nspin,thetadeg,phideg)

  use mConstants

  implicit none

  integer, intent(in) :: nspin,maxnh
  real(kdp), intent(in) :: thetadeg,phideg
  real(kdp), intent(inout) :: dm(maxnh,nspin)

  real(kdp) :: phi,theta,mtot(3)
  real(kdp) :: PI
  integer :: is
  real(kdp), allocatable :: dmbuf(:,:)

  PI=2.0_kdp * acos(0.0_kdp)

  phi=phideg * PI/180.0_kdp
  theta=thetadeg * PI/180.0_kdp
!  theta=0.0_kdp * PI

  mtot(3)=1.0_kdp
  mtot(1)=sin(theta) * cos(phi) * mtot(3)
  mtot(2)=sin(theta) * sin(phi) * mtot(3)
  mtot(3)=cos(theta) * mtot(3)
!  write(*,*)"angles=",sin(theta),cos(theta),sin(phi),cos(phi)
!  write(*,*)"mxyz=",mtot(1),mtot(2),mtot(3)

  if(nspin.gt.2)then
    dm(:,3)=0.5_kdp * (dm(:,1)-dm(:,2))
    dm(:,4)=0.5_kdp * (dm(:,1)+dm(:,2))
    dm(:,1)=dm(:,4)+mtot(3)*dm(:,3)
    dm(:,2)=dm(:,4)-mtot(3)*dm(:,3)
    dm(:,4)=-mtot(2)*dm(:,3)
    dm(:,3)=mtot(1)*dm(:,3)
  else
    allocate(dmbuf(maxnh,2))
    dmbuf(:,1)=0.5_kdp * (dm(:,1)-dm(:,2))
    dmbuf(:,2)=0.5_kdp * (dm(:,1)+dm(:,2))
    dm(:,1)=dmbuf(:,2)+mtot(3)*dmbuf(:,1)
    dm(:,2)=dmbuf(:,2)-mtot(3)*dmbuf(:,1)
    deallocate(dmbuf)
  endif

end subroutine rotateDM

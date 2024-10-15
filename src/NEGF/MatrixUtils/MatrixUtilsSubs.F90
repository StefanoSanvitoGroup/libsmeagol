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
!                   INDDENSETOINDSPARSEGENERAL,
!                   INDDENSETOINDSPARSEGENERAL_LEAD,
!                   FINDIJ,
!                   FINDNNZGF_BS,
!                   FINDNNZGF2P,
!                   FINDNNZGF2,
!                   FINDNNZGF,
!                   MATHERMITIAN,
!                   REORDER_NC,
!                   SPARSETODENSE,
!                   SPARSETODENSE_NC,
!                   PRINTSPARSE2DENSEREORDEREDNC,
!                   MATHERMITIANCRS  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine inddensetoindsparsegeneral(iuo,juo,ind2,mat)

 use mTypes

 implicit none
 integer iuo,juo,ind2
 integer ii
 type(matrixTypeGeneral) ::  mat

 do ind2=mat%matSparse%q(iuo),mat%matSparse%q(iuo+1)-1
   if(mat%matSparse%j(ind2).eq.juo)return
 enddo
 ind2=0

end subroutine inddensetoindsparsegeneral


! binary search
! Meilin Bai, 08 June 2013
function bsearch(a, n, x)
implicit none
integer, intent(in) :: a(n), n, x
integer bsearch, i, p, q 
p=1
q=n
do while (p .le. q)
    i =(p+q)/2
    if (x .eq. a(i)) then
        bsearch =i
        return
    elseif (x .lt. a(i)) then
        q =i-1
    else
        p =i+1
    endif
end do
bsearch =0
return
end function bsearch
    
logical function issorted(a, n, opt)
    implicit none
    integer, intent(in) :: a(n), n, opt
    integer i
    if (opt .eq. 1) then
       do i=2,n
           if (a(i-1) .gt. a(i)) then
               issorted =.false.
               return
           endif
       enddo
     elseif (opt .eq. -1) then
       do i=2,n
           if (a(i-1) .lt. a(i)) then
               issorted =.false.
               return
           endif
       enddo
     endif
     issorted =.true.
     return
 end function issorted



subroutine inddensetoindsparsegeneral_lead(iuo,juo,ind2,mat)

 use mTypes

 implicit none
 integer iuo,juo,ind2
 integer ii
 type(matrixTypeGeneral) ::  mat

 if(iuo<1.or.mat%matSparse%iRows<iuo)then
   ind2=0
   return
 endif
 do ii=mat%matSparse%q(iuo),mat%matSparse%q(iuo+1)-1
   ind2=ii
   if(mat%matSparse%j(ind2).eq.juo)return
 enddo
 ind2=0

end subroutine inddensetoindsparsegeneral_lead


subroutine findij(ind,nuotot,maxnhg,listhg,listhptrg,i,j)

 implicit none
   integer ind,nuotot,maxnhg,listhptrg(nuotot),i,j,listhg(maxnhg)
   integer ii,jj

  do ii=1,nuotot
    if(listhptrg(ii).ge.ind)exit
  enddo
  i=ii-1
  j=listhg(ind)

end subroutine findij


subroutine findnnzgf_bs(nnz,nnzrow,irows,icols,nleads,sigmamp,nebss,hgenerals)

 use mTypes
 use negfmod, only: outinfo

 implicit none
 integer nnz,irows,icols,nleads
 integer nnzrow(irows),ii,ind,jj,i1,i2,i3
 type(matrixTypeGeneral) :: sigmamp(nleads)
 type(matrixTypeGeneral) :: hgenerals
 integer nebss(nleads,2)

 nnzrow=0
 nnz=0
 do ii=1,irows
   columnsloop: do jj=1,icols
     call inddensetoindsparsegeneral(ii,jj,ind,hgenerals)
     if(ind.ne.0)then
       if(abs(hgenerals%matSparse%b(ind)).ne.0D0)then
         nnzrow(ii)=nnzrow(ii)+1
         cycle columnsloop
       endif
     endif
     
     do i1=1,nleads

       call inddensetoindsparsegeneral_lead(ii-nebss(i1,1)+1,jj-nebss(i1,1)+1,ind,sigmamp(i1))
       if (outinfo) write(12347,*)"findsparsesetind=",ii,jj,i1,ind
       if(ind.ne.0)then
         if(sigmamp(i1)%matSparse%b(ind).ne.0D0)then
           nnzrow(ii)=nnzrow(ii)+1
           cycle columnsloop
         endif
       endif
     enddo

   enddo columnsloop
   nnz=nnz+nnzrow(ii)
 enddo

end subroutine findnnzgf_bs



subroutine findnnzgf2P(nnz,nnzrow,irows,icols,nl,nr,sigmal,sigmar,hgenerals)

 use mTypes

 implicit none
 integer nnz,irows,icols,nl,nr
 integer nnzrow(irows),ii,ind,jj,nj,iiglobal,irowsglobal
 double complex sigmal(nl,nl),sigmar(nr,nr)
 type(matrixTypeGeneral) :: hgenerals

 double complex w(icols)
 integer idxrow(icols)

 integer  iend,istart


 istart=hgenerals%matSparseP%matSparse%iVert
 iend=istart+hgenerals%matSparseP%matSparse%iRows-1
 irowsglobal=hgenerals%matSparseP%iRowsGlobal

 nnzrow=0
 nnz=0
 w=0D0
 do ii=1,irows

   iiglobal=ii+istart-1
   nj=0

   do ind=hgenerals%matSparseP%matSparse%q(ii),hgenerals%matSparseP%matSparse%q(ii+1)-1
     if(hgenerals%matSparseP%matSparse%b(ind).ne.0D0)then
       jj=hgenerals%matSparseP%matSparse%j(ind)
       w(jj)=hgenerals%matSparseP%matSparse%b(ind)
       nj=nj+1
       idxrow(nj)=jj
     endif
   enddo

   if( iiglobal <= nl)then
     do jj=1,nl
       if(sigmal(iiglobal,jj).ne.0D0)then
         if(w(jj).ne.0D0)then
           w(jj)=w(jj)+sigmal(iiglobal,jj)
         else
           w(jj)=sigmal(iiglobal,jj)
           nj=nj+1
           idxrow(nj)=jj
         endif
       endif
     enddo
   endif

   if(iiglobal > iRowsGlobal-nr)then
     do jj=icols-nr+1,icols
       if(sigmar(iiglobal-iRowsGlobal+nr,jj-icols+nr).ne.0D0)then
         if(w(jj).ne.0D0)then
           w(jj)=w(jj)+sigmar(iiglobal-iRowsGlobal+nr,jj-icols+nr)
         else
           w(jj)=sigmar(iiglobal-iRowsGlobal+nr,jj-icols+nr)
           nj=nj+1
           idxrow(nj)=jj
         endif
       endif
     enddo
   endif

   nnzrow(ii)=nj
   nnz=nnz+nj

   do jj=1,nj
     w(idxrow(jj))=0D0
   enddo

 enddo

end subroutine findnnzgf2P


subroutine findnnzgf2(nnz,nnzrow,irows,icols,nl,nr,sigmal,sigmar,hgenerals)

 use mTypes

 implicit none
 integer nnz,irows,icols,nl,nr
 integer nnzrow(irows),ii,ind,jj,nj
 double complex sigmal(nl,nl),sigmar(nr,nr)
 type(matrixTypeGeneral) :: hgenerals

 double complex w(icols)
 integer idxrow(icols)

 nnzrow=0
 nnz=0
 w=0D0
 do ii=1,irows

   nj=0

   do ind=hgenerals%matSparse%q(ii),hgenerals%matSparse%q(ii+1)-1
!     if(.true..or.hgenerals%matSparse%b(ind).ne.0D0)then
     if(hgenerals%matSparse%b(ind).ne.0D0)then
       jj=hgenerals%matSparse%j(ind)
       w(jj)=hgenerals%matSparse%b(ind)
       nj=nj+1
       idxrow(nj)=jj
     endif
   enddo

   if(ii <= nl)then
     do jj=1,nl
       if(sigmal(ii,jj).ne.0D0)then
         if(w(jj).ne.0D0)then
           w(jj)=w(jj)+sigmal(ii,jj)
         else
           w(jj)=sigmal(ii,jj)
           nj=nj+1
           idxrow(nj)=jj
         endif
       endif
     enddo
   endif

   if(ii > irows-nr)then
     do jj=icols-nr+1,icols
       if(sigmar(ii-irows+nr,jj-icols+nr).ne.0D0)then
         if(w(jj).ne.0D0)then
           w(jj)=w(jj)+sigmar(ii-irows+nr,jj-icols+nr)
         else
           w(jj)=sigmar(ii-irows+nr,jj-icols+nr)
           nj=nj+1
           idxrow(nj)=jj
         endif
       endif
     enddo
   endif

   nnzrow(ii)=nj
   nnz=nnz+nj

   do jj=1,nj
     w(idxrow(jj))=0D0
   enddo

 enddo

end subroutine findnnzgf2


subroutine findnnzgf(nnz,nnzrow,irows,icols,nl,nr,sigmal,sigmar,hgenerals)

 use mTypes

 implicit none
 integer nnz,irows,icols,nl,nr
 integer nnzrow(irows),ii,ind,jj
 double complex sigmal(nl,nl),sigmar(nr,nr)
 type(matrixTypeGeneral) :: hgenerals

 nnzrow=0
 nnz=0
 do ii=1,irows
   do jj=1,icols
       call inddensetoindsparsegeneral(ii,jj,ind,hgenerals)
       if(ind.ne.0)then
         if(abs(hgenerals%matSparse%b(ind)).ne.0D0)then
           nnzrow(ii)=nnzrow(ii)+1
           cycle
         endif
       endif
       if((ii.le.nl).and.(jj.le.nl))then
         if(sigmal(ii,jj).ne.0D0)then
           nnzrow(ii)=nnzrow(ii)+1
           cycle
         endif
       endif
       if((ii.gt.irows-nr).and.(jj.gt.icols-nr))then
         if(sigmar(ii-irows+nr,jj-icols+nr).ne.0D0)then
           nnzrow(ii)=nnzrow(ii)+1
           cycle
         endif
       endif
   enddo
   nnz=nnz+nnzrow(ii)
 enddo

end subroutine findnnzgf


  SUBROUTINE mathermitian(mat,matdagger)

  use mTypes

  implicit none

  type(matrixTypeGeneral) :: mat,matdagger

  integer i,j,ind,irows,icols,nnzcol(mat%icols),inddagger

  irows=mat%irows
  icols=mat%icols

  nnzcol=0
  do i=1,irows
    do ind=mat%matSparse%q(i),mat%matSparse%q(i+1)-1
      j=mat%matSparse%j(ind)
      nnzcol(j)=nnzcol(j)+1
    enddo
  enddo

  matdagger%matSparse%q(1)=1
  do i=1,icols
    matdagger%matSparse%q(i+1)=matdagger%matSparse%q(i)+nnzcol(i)
  enddo

  nnzcol=0
  do i=1,irows
    do ind=mat%matSparse%q(i),mat%matSparse%q(i+1)-1
      j=mat%matSparse%j(ind)
      nnzcol(j)=nnzcol(j)+1
      inddagger=matdagger%matSparse%q(j)+nnzcol(j)-1
      matdagger%matSparse%j(inddagger)=i
      matdagger%matSparse%b(inddagger)=DCONJG(mat%matSparse%b(ind))
    enddo
  enddo

  end SUBROUTINE mathermitian

  subroutine reorder_nc(mat,n)

    use mTypes

    implicit none
    integer, intent(in) :: n
    double complex, intent(inout) :: mat(n,n)

    double complex matb(n,n)

    integer i1,j1,nh,ip1,ip2,jp1,jp2

    nh=n/2

    matb=mat
!    mat=1.0d+10
    do i1=1,n

      if(mod(i1,2)==1)then
        ip1=(i1+1)/2
      else
        ip1=i1/2+nh
      endif
      do j1=1,n
  
        if(mod(j1,2)==1)then
          jp1=(j1+1)/2
        else
          jp1=j1/2+nh
        endif

!        if(mat(ip1,jp1).ne.1.0d+10)write(*,*)"matdouble",ip1,jp1,j1,j1
!        write(12347,*)"m1=",ip1,jp1,j1,j1
        mat(ip1,jp1)=matb(i1,j1)

      enddo
      
    enddo

  end subroutine reorder_nc


  subroutine sparsetodense(mats,mat,n)

    use mTypes

    implicit none
    integer, intent(in) :: n
    type(matrixSparseType), intent(in) :: mats
    double complex, intent(out) :: mat(n,n)

    integer i1,ind

    mat=0.0D0
    do i1=1,n
      do ind=mats%q(i1),mats%q(i1+1)-1
        mat(i1,mats%j(ind))=mats%b(ind)
      enddo
    enddo

  end subroutine sparsetodense


  subroutine sparsetodense_nc(matg,mat,n,nspin)

    use mTypes

    implicit none
    integer, intent(in) :: n,nspin
    type(matrixTypeGeneral), intent(in) :: matg(nspin)
    double complex, intent(out) :: mat(n,n)

    double complex, allocatable :: mats(:,:,:)
    integer i1,i2,n2,ispin,ind

    n2=n/2
    allocate(mats(n2,n2,nspin))

    mats=0.0D0
!    mats=54321000000.D0

    do i1=1,n2
      do ind=matg(1)%matSparse%q(i1),matg(1)%matSparse%q(i1+1)-1
        do ispin=1,nspin
          mats(i1,matg(1)%matSparse%j(ind),ispin)=matg(ispin)%matSparse%b(ind)
        enddo
      enddo
    enddo

    mat=0.0D0
    mat(1:n2,1:n2)=mats(:,:,1)
    mat(n2+1:2* n2,n2+1:2*n2)=mats(:,:,2)
    mat(1:n2,n2+1:2*n2)=mats(:,:,3)
    mat(n2+1:2 * n2,1:n2)=mats(:,:,4)

!    put option to set this true or false in subroutine
!    mat(n2+1:2 * n2,1:n2)=dconjg(transpose(mat(1:n2,n2+1:2*n2)))

    deallocate(mats)

  end subroutine sparsetodense_nc




subroutine PrintSparse2DenseReorderedNC(mat,n,nspinmat,nspin,NspinComplexMatrix,operation,nam)

  use mTypes

  integer, intent(in) :: n,nspin,NspinComplexMatrix,operation,nspinmat
  type(matrixTypeGeneral), intent(in) :: mat(nspinmat)
  CHARACTER(LEN=*), intent(in) :: nam
  double complex, allocatable ::  matbuf(:,:)

!  return !comment out  if output is wanted
  if(nspin<4)return

  allocate(matbuf(2*n,2*n))
!  write(12347,*)"start output of ",nam,"N=",n
  matbuf=0.0D0
  if(NspinComplexMatrix==4)then
    if(operation==1)then
!      write(12346,*)"opi1",nam,nspinmat,nspin,NspinComplexMatrix
      call sparsetodense(mat(1)%matSparse,matbuf(1:n,1:n),n)
      call sparsetodense(mat(2)%matSparse,matbuf(n+1:2*n,n+1:2*n),n)
      call sparsetodense(mat(3)%matSparse,matbuf(1:n,n+1:2*n),n)
!      if(.true.)then
      if(.false.)then
        matbuf(n+1:2*n,1:n)=DCONJG(transpose(matbuf(1:n,n+1:2*n)))
      else
        call sparsetodense(mat(4)%matSparse,matbuf(n+1:2*n,1:n),n)
      endif
    elseif(operation==2)then
!      write(12346,*)"opi2",nam,nspinmat,nspin,NspinComplexMatrix
      call sparsetodense(mat(1)%matSparse,matbuf(1:n,1:n),n)
      matbuf(n+1:2*n,n+1:2*n)=matbuf(1:n,1:n)
    elseif(operation==3)then
!      write(12346,*)"opi3",nam,nspinmat,nspin,NspinComplexMatrix
      matbuf=mat(1)%matdense%a
    endif
  else
    call sparsetodense(mat(1)%matSparse,matbuf,2*n)
    call reorder_nc(matbuf,2*n)
  endif

  call writemat9(1,0.0D0,matbuf,2*n,2*n,0.0d0,nam)
  deallocate(matbuf)


end subroutine PrintSparse2DenseReorderedNC

  SUBROUTINE mathermitianCRS(mat,matdagger)

  use mTypes

  implicit none

  type(MatrixSparseType), intent(in) :: mat
  type(MatrixSparseType), intent(inout) :: matdagger

  integer i,j,ind,irows,icols,nnzcol(mat%icols),inddagger

  irows=mat%irows
  icols=mat%icols

  nnzcol=0
  do i=1,irows
    do ind=mat%q(i),mat%q(i+1)-1
      j=mat%j(ind)
      nnzcol(j)=nnzcol(j)+1
    enddo
  enddo

  matdagger%q(1)=1
  do i=1,icols
    matdagger%q(i+1)=matdagger%q(i)+nnzcol(i)
  enddo

  nnzcol=0
  do i=1,irows
    do ind=mat%q(i),mat%q(i+1)-1
      j=mat%j(ind)
      nnzcol(j)=nnzcol(j)+1
      inddagger=matdagger%q(j)+nnzcol(j)-1
      matdagger%j(inddagger)=i
      matdagger%b(inddagger)=DCONJG(mat%b(ind))
    enddo
  enddo

  end SUBROUTINE mathermitianCRS


  subroutine SetRhoToZero(NspinComplexMatrix,matgeneral)
  
  use mTypes

  integer NspinComplexMatrix
  type(matrixTypeGeneral) :: matgeneral(NspinComplexMatrix)

  integer ispin

  do ispin=1,NspinComplexMatrix
    matgeneral(ispin)%matSparse%b(:)=0.0D0
  enddo

  end subroutine SetRhoToZero





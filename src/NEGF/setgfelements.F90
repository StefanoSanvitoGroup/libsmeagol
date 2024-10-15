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
!                   SETGFELEMENTSSPARSE3P,
!                   SETGFELEMENTSSPARSE3,
!                   SETGFELEMENTSGENERAL_NC,
!                   SETGFELEMENTSGENERAL_NC_BS,
!                   SETGFELEMENTSDENSE_BS,
!                   SETGFELEMENTSSPARSE2_BS,
!                   SETGFELEMENTSSPARSE2,
!                   SETGFELEMENTSDENSE_BS_NC,
!                   SETGFELEMENTSDENSE_NC,
!                   SETGFELEMENTSDENSE  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

subroutine setgfelementssparse3P(ei,gfsparse,nl,nr,sigmal,sigmar,hgenerals,sgeneral)

 use mTypes

 implicit none

 type(matrixTypeGeneral) :: gfsparse
 type(matrixTypeGeneral) :: hgenerals,sgeneral
 integer i0,irows,icols,nl,nr
 double complex ei
 integer ii,ind,jj,nj
 double complex sigmal(nl,nl),sigmar(nr,nr)

 double complex w(hgenerals%icols)
 integer idxrow(hgenerals%icols)

 integer  iend,istart,iRowsGlobal,iiglobal


 istart=hgenerals%matSparseP%matSparse%iVert
 iend=istart+hgenerals%matSparseP%matSparse%iRows-1
 irowsglobal=hgenerals%matSparseP%iRowsGlobal


 icols=gfsparse%iCols
 irows=gfsparse%matSparseP%matSparse%iRows

 w=0D0

 gfsparse%matSparseP%matSparse%q(1)=1
 i0=1
 do ii=1,irows

   iiglobal=ii+istart-1
   nj=0

   do ind=hgenerals%matSparseP%matSparse%q(ii),hgenerals%matSparseP%matSparse%q(ii+1)-1
     if(hgenerals%matSparseP%matSparse%b(ind).ne.0D0)then
       jj=hgenerals%matSparseP%matSparse%j(ind)
       w(jj)=-hgenerals%matSparseP%matSparse%b(ind)+ei * sgeneral%matSparseP%matSparse%b(ind)
       nj=nj+1
       idxrow(nj)=jj
     endif
   enddo

   if(iiglobal <= nl)then
     do jj=1,nl
       if(sigmal(iiglobal,jj).ne.0D0)then
         if(w(jj).ne.0D0)then
           w(jj)=w(jj)-sigmal(iiglobal,jj)
         else
           w(jj)=-sigmal(iiglobal,jj)
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
           w(jj)=w(jj)-sigmar(iiglobal-iRowsGlobal+nr,jj-icols+nr)
         else
           w(jj)=-sigmar(iiglobal-iRowsGlobal+nr,jj-icols+nr)
           nj=nj+1
           idxrow(nj)=jj
         endif
       endif
     enddo
   endif

   gfsparse%matSparseP%matSparse%q(ii+1)=gfsparse%matSparseP%matSparse%q(ii)+nj

   i0=gfsparse%matSparseP%matSparse%q(ii)
   do ind=gfsparse%matSparseP%matSparse%q(ii),gfsparse%matSparseP%matSparse%q(ii+1)-1
     gfsparse%matSparseP%matSparse%j(ind)=idxrow(ind-i0+1)
     gfsparse%matSparseP%matSparse%b(ind)=w(idxrow(ind-i0+1))
   enddo


   do jj=1,nj
     w(idxrow(jj))=0D0
   enddo

 enddo

end subroutine setgfelementssparse3P


subroutine setgfelementssparse3(ei,gfsparse,nl,nr,sigmal,sigmar,hgenerals,sgeneral)

 use mTypes

 implicit none

 type(matrixTypeGeneral) :: gfsparse
 type(matrixTypeGeneral) :: hgenerals,sgeneral
 integer i0,irows,icols,nl,nr
 double complex ei
 integer ii,ind,jj,nj
 double complex sigmal(nl,nl),sigmar(nr,nr)

 double complex w(hgenerals%icols)
 integer idxrow(hgenerals%icols)


 icols=gfsparse%iCols
 irows=gfsparse%iRows

 w=0D0

 gfsparse%matSparse%q(1)=1
 i0=1
 do ii=1,irows

   nj=0

   do ind=hgenerals%matSparse%q(ii),hgenerals%matSparse%q(ii+1)-1
!     if(.true..or.hgenerals%matSparse%b(ind).ne.0D0)then
     if(hgenerals%matSparse%b(ind).ne.0D0)then
       jj=hgenerals%matSparse%j(ind)
       w(jj)=-hgenerals%matSparse%b(ind)+ei * sgeneral%matSparse%b(ind)
       nj=nj+1
       idxrow(nj)=jj
     endif
   enddo

   if(ii <= nl)then
     do jj=1,nl
       if(sigmal(ii,jj).ne.0D0)then
         if(w(jj).ne.0D0)then
           w(jj)=w(jj)-sigmal(ii,jj)
         else
           w(jj)=-sigmal(ii,jj)
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
           w(jj)=w(jj)-sigmar(ii-irows+nr,jj-icols+nr)
         else
           w(jj)=-sigmar(ii-irows+nr,jj-icols+nr)
           nj=nj+1
           idxrow(nj)=jj
         endif
       endif
     enddo
   endif

   gfsparse%matSparse%q(ii+1)=gfsparse%matSparse%q(ii)+nj

   i0=gfsparse%matSparse%q(ii)
   do ind=gfsparse%matSparse%q(ii),gfsparse%matSparse%q(ii+1)-1
     gfsparse%matSparse%j(ind)=idxrow(ind-i0+1)
     gfsparse%matSparse%b(ind)=w(idxrow(ind-i0+1))
   enddo


   do jj=1,nj
     w(idxrow(jj))=0D0
   enddo

 enddo

end subroutine setgfelementssparse3


subroutine setgfelementsgeneral_nc(ei,nspin,ispin,mat,nnz,irows,nl,nr,sigmal,sigmar,hgeneral,sgeneral)

 use mTypes

 implicit none
 integer, intent(in) :: nnz,irows,nl,nr,ispin,nspin
 type(matrixTypeGeneral) :: mat
 type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
 double complex ei
 double complex sigmal(nl,nl),sigmar(nr,nr)

 if(mat%mattype == 0)then
   if(nspin <= 2)then
     call setgfelementsdense(ei,mat%matdense,nnz,irows,nl,nr,sigmal,sigmar,hgeneral(ispin),sgeneral,ispin)
   else
     call setgfelementsdense_nc(ei,nspin,mat%matdense,nnz,irows,nl,nr,sigmal,sigmar,hgeneral,sgeneral)
   endif
 elseif(mat%mattype == 2)then
   call setgfelementssparse3(ei,mat,nl,nr,sigmal,sigmar,hgeneral(ispin),sgeneral)
 elseif(mat%mattype == 3)then
   call setgfelementssparse3P(ei,mat,nl,nr,sigmal,sigmar,hgeneral(ispin),sgeneral)
 endif

end subroutine setgfelementsgeneral_nc


subroutine setgfelementsgeneral_nc_bs(ei,nspin,ispin,mat,nnz,irows,nl,nr,nbss,nleads,sigmamp,nebss,hgeneral,sgeneral)

 use mTypes

 implicit none
 integer, intent(in) :: nnz,irows,nl,nr,ispin,nspin,nleads,nbss
 type(matrixTypeGeneral) :: mat
 type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
 type(matrixTypeGeneral), intent(in) :: sigmamp(nleads)
 integer, intent(in) :: nebss(nleads,2)
 double complex ei

 if(mat%mattype == 0)then
   if(nspin <= 2)then
     call setgfelementsdense_bs(ei,ispin,mat,nleads,sigmamp,nebss,hgeneral(ispin),sgeneral)
   else
!     call setgfelementsdense_bs_nc(ei,nspin,mat%matdense,nnz,irows,nl,nr,sigmal,sigmar,hgeneral,sgeneral)
!     call setgfelementsdense_bs_nc(ei,ispin,mat%matdense,nnz,irows,nl,nr,nleads,sigmamp,nebss,hgeneral(ispin),sgeneral)
     call setgfelementsdense_bs_nc(ei,nspin,mat%matdense,nbss,nleads,sigmamp,nebss,hgeneral,sgeneral)
   endif
 endif

end subroutine setgfelementsgeneral_nc_bs


subroutine setgfelementsdense_bs(ei,ispin,gf,nleads,sigmamp,nebss,hgenerals,sgeneral)

 use mTypes

 implicit none
 integer irows,icols,ispin,nleads,nebss(nleads,2)
 type(matrixTypeGeneral) :: gf,sigmamp(nleads)
 type(matrixTypeGeneral) :: hgenerals,sgeneral
 double complex ei
 integer ii,ind,jj,indadd,ind2,i1,i2
 logical elementset



 icols=gf%iCols
 irows=gf%iRows
 gf%matdense%a=0D0
 do ii=1,irows

   columnsloop: do jj=1,icols

     call inddensetoindsparsegeneral(ii,jj,ind,sgeneral)
     if(ind.ne.0)then
       if(abs(hgenerals%matSparse%b(ind)).ne.0D0)then
         gf%matdense%a(ii,jj)=ei * sgeneral%matSparse%b(ind)-hgenerals%matSparse%b(ind)
       endif
     endif
     
     do i1=1,nleads

       call inddensetoindsparsegeneral_lead(ii-nebss(i1,1)+1,jj-nebss(i1,1)+1,ind,sigmamp(i1))
       if(ind.ne.0)then
         if(sigmamp(i1)%matSparse%b(ind).ne.0D0)then
           gf%matdense%a(ii,jj)=gf%matdense%a(ii,jj)-sigmamp(i1)%matSparse%b(ind)
         endif
       endif
     enddo

   enddo columnsloop
 enddo

end subroutine setgfelementsdense_bs


subroutine setgfelementssparse2_bs(ei,gfsparse,nnz,irows,nnzrow,nleads,sigmamp,nebss,hgenerals,sgeneral)

 use mTypes
 use negfmod, only: outinfo

 implicit none
 integer nnz,irows,icols,nnzrow(irows),nleads,nebss(nleads,2)
 type(matrixTypeGeneral) :: gfsparse,sigmamp(nleads)
 double complex ei
 integer ii,ind,jj,indadd,ind2,i1,i2
 logical elementset
 type(matrixTypeGeneral) :: hgenerals,sgeneral

 icols=gfsparse%iCols
 irows=gfsparse%iRows

 gfsparse%matSparse%q(1)=1
 ind2=1
 do ii=1,irows

   columnsloop: do jj=1,icols

     elementset=.false.

     call inddensetoindsparsegeneral(ii,jj,ind,sgeneral)
     if(ind.ne.0)then
       if(abs(hgenerals%matSparse%b(ind)).ne.0D0)then
         gfsparse%matSparse%b(ind2)=ei * sgeneral%matSparse%b(ind)-hgenerals%matSparse%b(ind)
         if (outinfo) write(12347,*)"gfsparseset1=",ii,jj,gfsparse%matSparse%b(ind2)
         gfsparse%matSparse%j(ind2)=sgeneral%matSparse%j(ind)
         elementset=.true.
       endif
     endif
     
     do i1=1,nleads

       call inddensetoindsparsegeneral_lead(ii-nebss(i1,1)+1,jj-nebss(i1,1)+1,ind,sigmamp(i1))
       if (outinfo) write(12347,*)"gfsparsesetind=",ii,jj,i1,ind
       if(ind.ne.0)then
         if(sigmamp(i1)%matSparse%b(ind).ne.0D0)then
           gfsparse%matSparse%b(ind2)=gfsparse%matSparse%b(ind2)-sigmamp(i1)%matSparse%b(ind)
           if (outinfo) write(12347,*)"gfsparseset=",ii,jj,gfsparse%matSparse%b(ind2)
           gfsparse%matSparse%j(ind2)=sigmamp(i1)%matSparse%j(ind)+nebss(i1,1)-1
           elementset=.true.
         endif
       endif
     enddo

     if(elementset)then
       gfsparse%matSparse%q(ii+1)=ind2+1
       ind2=ind2+1
     endif

   enddo columnsloop
   gfsparse%matSparse%q(ii+1)=ind2
 enddo



end subroutine setgfelementssparse2_bs


subroutine setgfelementssparse2(ei,gfsparse,nnz,irows,nnzrow,nl,nr,sigmal,sigmar,hgenerals,sgeneral)

 use mTypes

 implicit none
 type(matrixTypeGeneral) :: gfsparse
 integer nnz,irows,icols,nl,nr,nnzrow(irows)
 double complex ei
 integer ii,ind,jj,indadd,ind2
 double complex sigmal(nl,nl),sigmar(nr,nr)
 type(matrixTypeGeneral) :: hgenerals,sgeneral

 icols=gfsparse%iCols
 irows=gfsparse%iRows

 nnz=1
 do ii=1,irows
   gfsparse%matSparse%q(ii)=nnz
   nnz=nnz+nnzrow(ii)
 enddo
 gfsparse%matSparse%q(irows+1)=nnz

 ind=1
 gfsparse%matSparse%b=0D0

 do ii=1,irows
   do jj=1,icols
     indadd=0

     call inddensetoindsparsegeneral(ii,jj,ind2,sgeneral)
     if(ind2.ne.0)then
       if(abs(hgenerals%matSparse%b(ind2)).ne.0D0)then
         gfsparse%matSparse%j(ind)=jj

         gfsparse%matSparse%b(ind)= ei*sgeneral%matSparse%b(ind2)-hgenerals%matSparse%b(ind2)
         indadd=1
       endif
     endif

     if((ii.le.nl).and.(jj.le.nl))then
       if(sigmal(ii,jj).ne.0D0)then
         gfsparse%matSparse%j(ind)=jj
         gfsparse%matSparse%b(ind)= gfsparse%matSparse%b(ind)-sigmal(ii,jj)
         indadd=1
       endif
     endif
     if((ii.gt.irows-nr).and.(jj.gt.icols-nr))then
         if(sigmar(ii-irows+nr,jj-icols+nr).ne.0D0)then
         gfsparse%matSparse%j(ind)=jj
         gfsparse%matSparse%b(ind)= gfsparse%matSparse%b(ind)-sigmar(ii-irows+nr,jj-icols+nr)
         indadd=1
       endif
     endif
     ind=ind+indadd
   enddo
 enddo

end subroutine setgfelementssparse2

subroutine setgfelementsdense_bs_nc(ei,nspin,gfdense,nbss,nleads,sigmamp,nebss,hgeneral,sgeneral)

 use mTypes

 implicit none
 integer, intent (in) :: nspin,nleads,nbss
 type(matrixType) :: gfdense
 double complex ei
 integer ii,ind,jj,indadd
 type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
 integer irows,icols,nl2,nr2,nl,nr,ni,ni2
 type(matrixTypeGeneral), intent(in) :: sigmamp(nleads)
 integer, intent(in) :: nebss(nleads,2)
 double complex, allocatable ::  sigmal(:,:),sigmar(:,:),sigmai(:,:)


 icols=gfdense%iCols/2
 irows=gfdense%iRows/2

 gfdense%a=0D0
 do ii=1,iRows
   do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
     gfdense%a(ii,sgeneral%matSparse%j(ind))= ei*sgeneral%matSparse%b(ind)-hgeneral(1)%matSparse%b(ind)
     gfdense%a(irows+ii,icols+sgeneral%matSparse%j(ind))= ei*sgeneral%matSparse%b(ind)-hgeneral(2)%matSparse%b(ind)
     gfdense%a(ii,icols+sgeneral%matSparse%j(ind))= -hgeneral(3)%matSparse%b(ind)
   enddo
 enddo
 gfdense%a(irows+1:2 * irows,1:icols)=dconjg(transpose(gfdense%a(1:irows,icols+1:2 * icols))) 
 
 nl=sigmamp(nleads-1)%iRows
 nr=sigmamp(nleads)%iRows
 nl2=nl/2
 nr2=nr/2

!write(12347,*)"nlr=",nl,nr,nleads
 allocate(sigmal(nl,nl),sigmar(nr,nr))

 call sparsetodense(sigmamp(nleads-1)%MatSparse,sigmal,nl)
 call sparsetodense(sigmamp(nleads)%MatSparse,sigmar,nr)


 gfdense%a(1:NL2,1:NL2)=gfdense%a(1:NL2,1:NL2) - sigmal(1:nl2,1:nl2)
 gfdense%a(1:NL2,iCols+1:iCols+NL2)=gfdense%a(1:NL2,iCols+1:iCols+NL2) - sigmal(1:nl2,nl2+1:2 * nl2)
 gfdense%a(iRows+1:iRows+NL2,1:NL2)=gfdense%a(iRows+1:iRows+NL2,1:NL2) - sigmal(nl2+1:2 * nl2,1:nl2)
 gfdense%a(iRows+1:iRows+NL2,iCols+1:iCols+NL2)=gfdense%a(iRows+1:iRows+NL2,iCols+1:iCols+NL2) - sigmal(nl2+1:2 * nl2,nl2+1:2 *nl2)
 
 
 gfdense%a(iRows-NR2+1:iRows,iCols-NR2+1:iCols)= gfdense%a(iRows-NR2+1:iRows,iCols-NR2+1:iCols)-sigmar(1:nr2,1:nr2)
 gfdense%a(iRows-NR2+1:iRows,iCols+iCols-NR2+1:iCols+iCols)= gfdense%a(iRows-NR2+1:iRows,iCols+iCols-NR2+1:iCols+iCols)-sigmar(1:nr2,nr2+1:nr2+nr2)
 gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols-NR2+1:iCols)= gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols-NR2+1:iCols)-sigmar(nr2+1:nr2+nr2,1:nr2)
 gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols+iCols-NR2+1:iCols+iCols)= gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols+iCols-NR2+1:iCols+iCols)-sigmar(nr2+1:nr2+nr2,nr2+1:nr2+nr2)


 deallocate(sigmal,sigmar)

 do ii=1,nbss
!   write(12347,*)"ibs=",ii,nebss(ii,1),nebss(ii,2)
   ni=nebss(ii,2)-nebss(ii,1)+1
   allocate(sigmai(ni,ni))
   call sparsetodense(sigmamp(ii)%MatSparse,sigmai,ni)
   gfdense%a(nebss(ii,1):nebss(ii,2),nebss(ii,1):nebss(ii,2))=gfdense%a(nebss(ii,1):nebss(ii,2),nebss(ii,1):nebss(ii,2)) - sigmai(1:ni,1:ni)
   gfdense%a(irows+nebss(ii,1):irows+nebss(ii,2),irows+nebss(ii,1):irows+nebss(ii,2))=gfdense%a(irows+nebss(ii,1):irows+nebss(ii,2),irows+nebss(ii,1):irows+nebss(ii,2)) - sigmai(1:ni,1:ni)
   deallocate(sigmai)
 enddo

end subroutine setgfelementsdense_bs_nc


subroutine setgfelementsdense_nc(ei,nspin,gfdense,nnz,irowstotal,nl,nr,sigmal,sigmar,hgeneral,sgeneral)

 use mTypes
 use ScissorOperator, only : SCO_istart, SCO_nob, SCO_Hblock, SCOSetHamiltonianBlock
 use negfmod, only: outinfo

 implicit none
 integer, intent (in) :: nnz,nl,nr,nspin,irowstotal
 type(matrixType) :: gfdense
 double complex ei
 integer ii,ind,jj,indadd
 double complex sigmal(nl,nl),sigmar(nr,nr)
 type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral
 integer irows,icols,nl2,nr2

 icols=gfdense%iCols/2
 irows=gfdense%iRows/2
 nl2=nl/2
 nr2=nr/2

 if(.not.SCOSetHamiltonianBlock)then
   gfdense%a=0D0
   do ii=1,iRows
     do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
       gfdense%a(ii,sgeneral%matSparse%j(ind))= ei*sgeneral%matSparse%b(ind)-hgeneral(1)%matSparse%b(ind)
       gfdense%a(irows+ii,icols+sgeneral%matSparse%j(ind))= ei*sgeneral%matSparse%b(ind)-hgeneral(2)%matSparse%b(ind)
       gfdense%a(ii,icols+sgeneral%matSparse%j(ind))= -hgeneral(3)%matSparse%b(ind)
     enddo
   enddo
  !yyy: check for complex energies, we have to take transpose only of h, not e s-h!!
   gfdense%a(irows+1:2 * irows,1:icols)=dconjg(transpose(gfdense%a(1:irows,icols+1:2 * icols))) 
 else
   if (outinfo) write(12347,*)"Applying scissor operator to Hamiltonian", dreal(ei),dimag(ei)

   gfdense%a=0D0
   do ii=1,irows
     do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
       gfdense%a(ii,sgeneral%matSparse%j(ind))= ei*sgeneral%matSparse%b(ind)
       gfdense%a(irows+ii,icols+sgeneral%matSparse%j(ind))= ei*sgeneral%matSparse%b(ind)
     enddo
   enddo
   
   gfdense%a(1:irows,1:icols)=gfdense%a(1:irows,1:icols)-SCO_Hblock(:,:,1)
   gfdense%a(irows+1:2*irows,icols+1:2*icols)=gfdense%a(irows+1:2*irows,icols+1:2*icols)-SCO_Hblock(:,:,2)
   gfdense%a(1:irows,icols+1:2*icols)=gfdense%a(1:irows,icols+1:2*icols)-SCO_Hblock(:,:,3)
   gfdense%a(irows+1:2*irows,1:icols)=gfdense%a(irows+1:2*irows,1:icols)-SCO_Hblock(:,:,4)

 endif




 gfdense%a(1:NL2,1:NL2)=gfdense%a(1:NL2,1:NL2) - sigmal(1:nl2,1:nl2)
 gfdense%a(1:NL2,iCols+1:iCols+NL2)=gfdense%a(1:NL2,iCols+1:iCols+NL2) - sigmal(1:nl2,nl2+1:2 * nl2)
 gfdense%a(iRows+1:iRows+NL2,1:NL2)=gfdense%a(iRows+1:iRows+NL2,1:NL2) - sigmal(nl2+1:2 * nl2,1:nl2)
 gfdense%a(iRows+1:iRows+NL2,iCols+1:iCols+NL2)=gfdense%a(iRows+1:iRows+NL2,iCols+1:iCols+NL2) - sigmal(nl2+1:2 * nl2,nl2+1:2 *nl2)
 
 
 gfdense%a(iRows-NR2+1:iRows,iCols-NR2+1:iCols)= gfdense%a(iRows-NR2+1:iRows,iCols-NR2+1:iCols)-sigmar(1:nr2,1:nr2)
 gfdense%a(iRows-NR2+1:iRows,iCols+iCols-NR2+1:iCols+iCols)= gfdense%a(iRows-NR2+1:iRows,iCols+iCols-NR2+1:iCols+iCols)-sigmar(1:nr2,nr2+1:nr2+nr2)
 gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols-NR2+1:iCols)= gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols-NR2+1:iCols)-sigmar(nr2+1:nr2+nr2,1:nr2)
 gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols+iCols-NR2+1:iCols+iCols)= gfdense%a(iRows+iRows-NR2+1:iRows+iRows,iCols+iCols-NR2+1:iCols+iCols)-sigmar(nr2+1:nr2+nr2,nr2+1:nr2+nr2)


end subroutine setgfelementsdense_nc



subroutine setgfelementsdense(ei,gfdense,nnz,irows,nl,nr,sigmal,sigmar,hgenerals,sgeneral,ispin)

 use mTypes
 use ScissorOperator, only : SCO_istart, SCO_nob, SCO_Hblock, SCOSetHamiltonianBlock
 use negfmod, only: outinfo
! use negfmod, only: ikpmod

 implicit none
 type(matrixType) :: gfdense
 integer, intent(in) :: ispin
 integer nnz,irows,icols,nl,nr
 double complex ei
 integer j,ii,ind,jj,indadd,iend
 double complex sigmal(nl,nl),sigmar(nr,nr)
 type(matrixTypeGeneral) :: hgenerals,sgeneral
! double complex mat(irows,irows),matblock(sco_nob,sco_nob) !only for debugging purposes

 icols=gfdense%iCols
 irows=gfdense%iRows

 gfdense%a=0D0
 do ii=1,iRows
   do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1
     gfdense%a(ii,sgeneral%matSparse%j(ind))= ei*sgeneral%matSparse%b(ind)-hgenerals%matSparse%b(ind)
   enddo
 enddo

!****************************begin SCO code block
 if(SCOSetHamiltonianBlock)then
   if (outinfo) write(12347,*)"Applying scissor operator to Hamiltonian",  SCO_istart,SCO_istart+SCO_nob-1,dreal(ei),dimag(ei)

   iend=SCO_istart+SCO_nob-1

!******************debugging************************************************************
!    call sparsetodense(hgenerals%matSparse,mat,iRows)
!    matblock=mat(SCO_istart:iend,SCO_istart:iend)
!***************************************************************************************

   do ii=SCO_istart,iend

     do ind=sgeneral%matSparse%q(ii),sgeneral%matSparse%q(ii+1)-1

       j=sgeneral%matSparse%j(ind)
       if(j>=SCO_istart .and. j <= iend)then
         gfdense%a(ii,j)= ei*sgeneral%matSparse%b(ind)
       endif

     enddo

!******************debugging************************************************************
!    do j=SCO_istart,iend
!      gfdense%a(ii,j)= gfdense%a(ii,j) -  matblock(ii-SCO_istart+1,j-SCO_istart+1)
!    enddo
!!!***************************************************************************************
     do j=SCO_istart,iend
         gfdense%a(ii,j)= gfdense%a(ii,j)-SCO_Hblock(ii-SCO_istart+1,j-SCO_istart+1,ispin)
     enddo
   
   enddo
!!!!******************debugging************************************************************
!!!  do ii=1,SCO_nob
!!!    do j=1,SCO_nob
!!!      if(abs(dreal(matblock(ii,j)-SCO_Hblock(ii,j,ispin)))>1e-8)then
!!!        write(12347,*)"deltah=",ii,j,ispin,ikpmod,dreal(matblock(ii,j)-SCO_Hblock(ii,j,ispin)),dimag(matblock(ii,j)-SCO_Hblock(ii,j,ispin)),matblock(ii,j),SCO_Hblock(ii,j,ispin)
!!!      endif
!!!    enddo
!!!  enddo

!   call stopnegf
!!!!***************************************************************************************
 endif

 gfdense%a(1:NL,1:NL)=gfdense%a(1:NL,1:NL)- sigmal(:,:)
 gfdense%a(iRows-NR+1:iRows,iCols-NR+1:iCols)= gfdense%a(iRows-NR+1:iRows,iCols-NR+1:iCols)-sigmar(:,:)

end subroutine setgfelementsdense




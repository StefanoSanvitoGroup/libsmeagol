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
!                   RHO_BOUNDARYL,
!                   RHO_BOUNDARYR  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine rho_boundaryL(nn,nuoL,maxnhL,listdptrL,listdL,numdL,nuotot,maxnhg,numhg,listhptrg,listhg,neqL,iequivL,nuoR,periodic,NeqOffL,iequivOffL)

    implicit none

    integer nuoL,nuotot,maxnhg,maxnhL,numhg(nuotot),listhptrg(nuotot),listhg(maxnhg),listdptrL(nuoL),listdL(maxnhL),neqL,iequivL(maxnhL,2),nuoR,NeqOffL,iequivOffL(maxnhL,2),numdL(nuoL)
    logical periodic
    integer iuo,i,j,ind,juo,nn,indL,juoL,jL

    do iuo = 1,nuoL
     do j = 1,numhg(iuo)
      ind = listhptrg(iuo) + j
      juo = listhg(ind)
      do i = 1,nn
       if(juo.gt.nuotot*(i-1).and.juo.le.nuotot*(i-1)+nuoL) then
        indL = listdptrL(iuo) + 1
        juoL = listdL(indL)
        jL = 1
        do while (((juoL-nuoL*(i-1)).ne.(juo-nuotot*(i-1))) .and.(jL .le. numdL(iuo)))
         indL = listdptrL(iuo) + jL
         juoL = listdL(indL)
         jL = jL + 1
        enddo
        if ((juoL-nuoL*(i-1)).eq.(juo-nuotot*(i-1))) then
         NeqL=NeqL+1
         iequivL(NeqL,1)=ind
         iequivL(NeqL,2)=indL
!    write(200,'(7i7)') ind,indL,juoL,juo,i,juoL-nuoL*(i-1),
!          juo-nuotot*(i-1)
        endif
       endif
       if(juo .gt. nuotot*i-nuoR .and. juo .le. nuotot*i.and. periodic) then
        indL = listdptrL(iuo) + 1
        juoL = listdL(indL)
        jL = 1
        do while ((juoL-nuoL*(2*nn+i-1)).ne.(juo-(nuotot*i-nuoR)) .and. (jL .le. numdL(iuo)))
         indL = listdptrL(iuo) + jL
         juoL = listdL(indL)
         jL = jL + 1
        enddo
        if ((juoL-nuoL*(2*nn+i-1)).eq.(juo-(nuotot*i-nuoR))) then
         NeqOffL=NeqOffL+1
         iequivOffL(NeqOffL,1)=ind
         iequivOffL(NeqOffL,2)=indL
! write(201,'(7i7)') ind,indL,juoL,juo,i,juoL-nuoL*(2*nn+i-1),
!          juo-(nuotot*i-nuoR)
         endif
       endif
      enddo
     enddo
    enddo
end subroutine rho_boundaryL


subroutine rho_boundaryR(nn,nuoR,maxnhR,listdptrR,listdR,numdR,nuotot,maxnhg,numhg,listhptrg,listhg,neqR,iequivR,nuoL,periodic,NeqOffR,iequivOffR)

    implicit none

    integer nuoR,nuotot,maxnhg,maxnhR,numhg(nuotot),listhptrg(nuotot),listhg(maxnhg),listdptrR(nuoR),listdR(maxnhR),neqR,iequivR(maxnhR,2),nuoL,NeqOffR,iequivOffR(maxnhR,2),numdR(nuoR)
    logical periodic
    integer iuo,i,j,ind,juo,nn,indR,juoR,jR

     do iuo = nuotot-nuoR+1,nuotot
           do j = 1, numhg(iuo)
            ind = listhptrg(iuo) + j
            juo = listhg(ind)
            do i = 1,nn
             if(juo.gt.nuotot*i-nuoR.and.juo.le.nuotot*i) then
              indR = listdptrR(iuo-(nuotot-nuoR)) + 1
              juoR = listdR(indR)
              jR = 1
              do while ((juoR-nuoR*(i-1)).ne.(juo-(nuotot*i-nuoR)) .and.(jR .le. numdR(iuo-(nuotot-nuoR))))
               indR = listdptrR(iuo-(nuotot-nuoR)) + jR
               juoR = listdR(indR)
               jR = jR + 1
              enddo
              if ((juoR-nuoR*(i-1)).eq.(juo-(nuotot*i-nuoR))) then
               NeqR=NeqR+1
               iequivR(NeqR,1)=ind
               iequivR(NeqR,2)=indR
              endif
             endif
             if(juo.gt.nuotot*(i-1).and.juo.le.nuotot*(i-1)+nuoL .and.periodic) then
              indR = listdptrR(iuo-(nuotot-nuoR)) + 1
              juoR = listdR(indR)
              jR = 1
              do while ((juoR-nuoL*(nn+i-1)).ne.(juo-nuotot*(i-1)) .and.(jR .le. numdR(iuo-(nuotot-nuoR))))
               indR = listdptrR(iuo-(nuotot-nuoR)) + jR
               juoR = listdR(indR)
               jR = jR + 1
              enddo
              if ((juoR-nuoL*(nn+i-1)).eq.(juo-nuotot*(i-1))) then
               NeqOffR=NeqOffR+1
               iequivOffR(NeqOffR,1)=ind
               iequivOffR(NeqOffR,2)=indR
              endif
             endif
            enddo
           enddo
          enddo









end subroutine rho_boundaryR



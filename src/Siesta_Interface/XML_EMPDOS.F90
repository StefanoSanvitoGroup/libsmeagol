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
!                   XML_EMPDOS  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

  SUBROUTINE xml_empdos(iufileempdos,i)

  use atmfuncs, only: symfio,labelfis,cnfigfio,lofio,mofio,zetafio
  use negfmod, only: em_nuo,em_nau,em_isa,em_iaorb,em_iphorb, gamma_negf

  implicit none

  integer iufileempdos,i,atm_spec(em_nuo)
  character atm_label(em_nuo)*20

  atm_label(i)=labelfis(em_isa(em_iaorb(i)))
  atm_spec(i)=em_isa(em_iaorb(i))

  WRITE(iufileempdos,'(a8)') "<orbital"
  WRITE(iufileempdos,'(a8,i25,a1)') " index=""",i,""""
  WRITE(iufileempdos,'(a13,i25,a1)') " atom_index=""",em_iaorb(i), """"
  WRITE(iufileempdos,'(a10,a,a1)') " species=""",trim(atm_label(i)), """"
  WRITE(iufileempdos,'(a)')" position=""   0.000000   0.000000   0.000000"""
  WRITE(iufileempdos,'(a4,i25,a1)') " n=""",cnfigfio(atm_spec(i), em_iphorb(i)), """"
  WRITE(iufileempdos,'(a4,i25,a1)') " l=""",lofio(atm_spec(i), em_iphorb(i)), """"
  WRITE(iufileempdos,'(a4,i25,a1)') " m=""",mofio(atm_spec(i), em_iphorb(i)), """"
  WRITE(iufileempdos,'(a4,i25,a1)') " z=""",zetafio(atm_spec(i), em_iphorb(i)), """"
  WRITE(iufileempdos,'(a1)') ">"
  WRITE(iufileempdos,'(a6)') "<data>"

  end SUBROUTINE xml_empdos





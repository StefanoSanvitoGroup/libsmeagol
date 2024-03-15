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
!                   ADDMULTDARRAYOMP,
!                   SETCONSTDARRAYOMP,
!                   ADDARRAYOMP,
!                   TRANSPOSEMATRIXCOMP,
!                   CONJUGATETRANSPOSEMATRIXCOMP,
!                   COPYMATRIXCOMP,
!                   SUBSTRACTMULTIPLYCMATRIXCOMP  
! AND
! THE MODULE
!                   MMATRIXUTILOMP  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

module mMatrixUtilOMP

  implicit none
  private
  
  public ::  AddMultDArrayOMP
  public ::  SetConstDArrayOMP
  public ::  AddArrayOMP
  public ::  TransposeMatrixCOMP
  public ::  ConjugateTransposeMatrixCOMP
  public ::  CopyMatrixCOMP
  public ::  SubstractMultiplyCMatrixCOMP

  contains

  subroutine AddMultDArrayOMP(v1,v2,a,n)

  integer, intent(in) :: n
  double precision, intent(in) :: a
  double complex, intent(inout) :: v1(n)
  double complex, intent(in) :: v2(n)

  integer i

!$omp parallel do default(shared) private(i) schedule(dynamic)
  do i=1,n 
    v1(i)=v1(i)+a*v2(i)
  enddo
!$omp end parallel do



  end subroutine AddMultDArrayOMP

  subroutine SetConstDArrayOMP(v1,a,n)

  integer, intent(in) :: n
  double precision, intent(in) :: a
  double complex, intent(out) :: v1(n)

  integer i

!$omp parallel do default(shared) private(i) schedule(dynamic)
  do i=1,n 
    v1(i)=a
  enddo
!$omp end parallel do



  end subroutine SetConstDArrayOMP


  subroutine AddArrayOMP(v1,v2,n)

  integer, intent(in) :: n
  double complex, intent(inout) :: v1(n)
  double complex, intent(in) :: v2(n)

  integer i

!$omp parallel do default(shared) private(i) schedule(dynamic)
  do i=1,n 
    v1(i)=v1(i)+v2(i)
  enddo
!$omp end parallel do



  end subroutine AddArrayOMP

  subroutine TransposeMatrixCOMP(m1,m2,n1,n2)

  integer, intent(in) :: n1,n2
  double complex, intent(out) :: m1(n1,n2)
  double complex, intent(in) :: m2(n2,n1)

  integer i1,i2

!$omp parallel do default(shared) private(i1,i2) schedule(dynamic)
  do i2=1,n2 
    do i1=1,n1 
      m1(i1,i2)=m2(i2,i1)
    enddo
  enddo
!$omp end parallel do



  end subroutine TransposeMatrixCOMP



  subroutine ConjugateTransposeMatrixCOMP(m1,m2,n1,n2)

  integer, intent(in) :: n1,n2
  double complex, intent(out) :: m1(n1,n2)
  double complex, intent(in) :: m2(n2,n1)

  integer i1,i2

!$omp parallel do default(shared) private(i1,i2) schedule(dynamic)
  do i2=1,n2 
    do i1=1,n1 
      m1(i1,i2)=DCONJG(m2(i2,i1))
    enddo
  enddo
!$omp end parallel do



  end subroutine ConjugateTransposeMatrixCOMP

  subroutine CopyMatrixCOMP(m1,m2,n1,n2)

  integer, intent(in) :: n1,n2
  double complex, intent(out) :: m1(n1,n2)
  double complex, intent(in) :: m2(n1,n2)

  integer i1,i2

!$omp parallel do default(shared) private(i1,i2) schedule(dynamic)
  do i2=1,n2 
    do i1=1,n1 
      m1(i1,i2)=m2(i1,i2)
    enddo
  enddo
!$omp end parallel do

  end subroutine CopyMatrixCOMP

  subroutine SubstractMultiplyCMatrixCOMP(m0,m1,m2,n1,n2,a)

  integer, intent(in) :: n1,n2
  double complex, intent(in) :: a
  double complex, intent(out) :: m0(n1,n2)
  double complex, intent(in) :: m1(n1,n2)
  double complex, intent(in) :: m2(n1,n2)

  integer i1,i2

!$omp parallel do default(shared) private(i1,i2) schedule(dynamic)
  do i2=1,n2 
    do i1=1,n1 
      m0(i1,i2)=a*(m1(i1,i2)-m2(i1,i2))
    enddo
  enddo
!$omp end parallel do


  end subroutine SubstractMultiplyCMatrixCOMP





end module mMatrixUtilOMP

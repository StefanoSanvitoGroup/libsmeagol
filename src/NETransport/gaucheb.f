! 
! Copyright (c) Smeagol Authors:
! A. R. Rocha, V. Garcia-Suarez, S. Bailey, C. J. Lambert, J. Ferrer and
! S. Sanvito 2003-2005
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
! SMEAGOL IS DISTRIBUTED ONLY THROUGH THE OFICIAL WEBSITE (www.smeagol.tcd.ie)
! UPON COMPLETION OF THE "SMEAGOL ACADEMIC LICENSE".
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
	subroutine GauCheb(X0,Xf,X,W,n)

C *****************************************************************
C Calculates the points and weights for numerical integration
C using Gauss-Chebychev polynomials
C
C Written by Alexandre Reily Rocha, June 2003 
C Computational Spintronics Group
C Trinity College Dublin
C e-mail: rochaa@tcd.ie
C ********************* HISTORY ***********************************
C Original version:	June 2003
C ******************* INPUTS **************************************
C X0		: Start of integration range
C Xf		: End of integration range
C n		: number of points
C ******************* OUTPUT **************************************
C X(n)		: Vector containing the points for integration
C W(n)		: Weights for the points of integration
C *****************************************************************

	double precision, parameter :: Pi = 3.141592654d0
	integer :: n
	double precision :: x0, xf
	double precision :: X(n), W(n)

	if (xf .lt. x0) then
	 minus=-1.d0
	else
	 minus = 1.d0
	endif

	if (n.gt.0) then
	 do i=1,n
	  yj=dcos(Pi*i/(n+1))
	  W(i) = minus*(Xf-X0)*Pi*SQRT(1-yj**2)/
     &     (2.D0*(n+1))
	  X(i) = (yj*(Xf-X0) + Xf +X0)/2.D0
	 enddo
	endif

	end

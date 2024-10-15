module negfcoop

implicit none

type coopcohp
   integer :: nbond
   integer :: nom
   logical :: ccoop
   integer, allocatable :: ia1(:), ia2(:)
   integer, allocatable :: no1(:), no2(:) 
   integer, allocatable :: io1(:,:), io2(:,:) 
end type coopcohp

type(coopcohp), save :: coopinfo

contains

!!! Get COOP/COHP option
subroutine read_coop_option(iu)

    use atomlist
    use fdf

#ifdef MPI
    use mpi_siesta
#endif
    implicit none

    integer iu, i, j
    integer Node, Nodes


#ifdef MPI
    integer    MPIerror, INFO
#endif

! Get Node number
#ifdef MPI
    call MPI_Comm_Rank(MPI_Comm_World,Node,MPIerror)
    call MPI_Comm_Size(MPI_Comm_World,Nodes,MPIerror)
#else
    Node = 0
#endif

    coopinfo%nom =0 ! max number of orbital
! Bonds
    allocate(coopinfo%ia1(coopinfo%nbond), coopinfo%ia2(coopinfo%nbond))

    if (Node ==0 ) then
       if (fdf_block('EM.COOPBonds', iu)) then
           do i =1, coopinfo%nbond
               read(iu, *) coopinfo%ia1(i), coopinfo%ia2(i)
           enddo
       else
           write(6,'(a)') 'No bonds defined. Stop!'
#ifdef MPI
           call MPI_Abort( MPI_Comm_World, INFO, MPIerror )
#else
           stop
#endif
       endif ! fdf_block
    endif ! IOnode

#ifdef MPI
    call MPI_Bcast(coopinfo%ia1(1),coopinfo%nbond, MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(coopinfo%ia2(1),coopinfo%nbond, MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

! orbitals
    ! calculate on each node
    do i=1,coopinfo%nbond
       if (coopinfo%nom .lt. lasto(coopinfo%ia1(i))-lasto(coopinfo%ia1(i)-1)) &
           coopinfo%nom =lasto(coopinfo%ia1(i))-lasto(coopinfo%ia1(i)-1)
    
       if (coopinfo%nom .lt. lasto(coopinfo%ia2(i))-lasto(coopinfo%ia2(i)-1)) &
           coopinfo%nom =lasto(coopinfo%ia2(i))-lasto(coopinfo%ia2(i)-1)
    enddo

    allocate(coopinfo%io1(coopinfo%nom,coopinfo%nbond),coopinfo%io2(coopinfo%nom,coopinfo%nbond))
    allocate(coopinfo%no1(coopinfo%nbond),coopinfo%no2(coopinfo%nbond))

    coopinfo%io1 =0
    coopinfo%io2 =0

    ! read orbital
    if (Node .eq. 0) then
        if (fdf_block('EM.COOPOrbitalsOfAtoms', iu)) then
! two lines for per bond
           do i=1,coopinfo%nbond
              ! Atom 1
              read(iu,*) coopinfo%no1(i),(coopinfo%io1(j,i),j=1,coopinfo%no1(i))
              if (coopinfo%no1(i) .eq. -1) then
                  coopinfo%no1(i) = lasto(coopinfo%ia1(i))-lasto(coopinfo%ia1(i)-1)
                  do j=1,coopinfo%no1(i)
                     coopinfo%io1(j,i)= 1+lasto(coopinfo%ia1(i)-1)
                  enddo
              else
                  coopinfo%io1(1:coopinfo%no1(i),i) =coopinfo%io1(1:coopinfo%no1(i),i) +lasto(coopinfo%ia1(i)-1)
              endif

              ! Atom 2
              read(iu,*) coopinfo%no2(i),(coopinfo%io2(j,i),j=1,coopinfo%no2(i))
              if (coopinfo%no2(i) .eq. -1) then
                  coopinfo%no2(i) = lasto(coopinfo%ia2(i))-lasto(coopinfo%ia2(i)-1)
                  do j=1,coopinfo%no2(i)
                     coopinfo%io2(j,i)= 1+lasto(coopinfo%ia2(i)-1)
                  enddo
              else
                  coopinfo%io2(1:coopinfo%no2(i),i) =coopinfo%io2(1:coopinfo%no2(i),i) +lasto(coopinfo%ia2(i)-1)
              endif
           enddo
       else
           do i=1,coopinfo%nbond
              coopinfo%no1(i)=lasto(coopinfo%ia1(i))-lasto(coopinfo%ia1(i)-1)
              coopinfo%no2(i)=lasto(coopinfo%ia2(i))-lasto(coopinfo%ia2(i)-1)

              do j=1,coopinfo%no1(i)
                 coopinfo%io1(j,i)=lasto(coopinfo%ia1(i)-1) +j
              enddo
              do j=1,coopinfo%no2(i)
                 coopinfo%io2(j,i)=lasto(coopinfo%ia2(i)-1) +j
              enddo
           enddo
       endif ! fdf_block
    endif ! Node
#ifdef MPI
    call MPI_Bcast(coopinfo%io1(1,1), coopinfo%nbond*coopinfo%nom, MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(coopinfo%io2(1,1), coopinfo%nbond*coopinfo%nom, MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(coopinfo%no1(1),coopinfo%nbond, MPI_integer,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(coopinfo%no2(1),coopinfo%nbond, MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
    return

end subroutine read_coop_option

!coop
subroutine calc_coop(N1, NL, NR, N1Half, NLHalf, NRHalf, nspinblocks, nspinzmatrix, wk, GammaL, GammaR, gfgeneral, hgeneral,sgeneral, nbond, co1, co2, ch1, ch2)
  use mTypes

  implicit none
  integer, intent(in) :: N1, NL, NR, N1Half, NLHalf, NRHalf, nspinblocks, nbond, nspinzmatrix
  type(matrixTypeGeneral), intent(in) :: gfgeneral, hgeneral(nspinzmatrix), sgeneral
  double complex, intent(in) :: GammaL(NL,NL), GammaR(NR,NR)
  double precision ::  co1(nbond), co2(nbond), ch1(nbond), ch2(nbond), wk

  integer j, k
  double complex :: GammaLFull(N1,N1), GammaRFull(N1,N1)

  GammaLFull=(0.d0, 0.d0)
  GammaRFull=(0.d0, 0.d0)
  if (NspinBlocks .le. 2) then
      GammaLFull(1:NL,1:NL) =GammaL
      GammaRFull(N1-NR+1:N1,N1-NR+1:N1) =GammaR
  else
      do j=1,2
         do k=1,2
            GammaLFull((j-1)*N1Half+1:(j-1)*N1Half+NLHalf,(k-1)*N1Half+1:(k-1)*N1Half+NLHalf) =GammaL((j-1)*NLHalf+1:j*NLHalf,(k-1)*NLHalf+1:k*NLHalf)
            GammaRFull(j*N1Half-NRHalf+1:j*N1Half,k*N1Half-NRHalf+1:k*N1Half) =GammaR((j-1)*NRHalf+1:j*NRHalf,(k-1)*NRHalf+1:k*NRHalf)
          enddo
      enddo
  endif
  call gfcoop(N1,NspinBlocks, Nspinzmatrix, wk, GammaLFull, GammaRFull, gfgeneral, hgeneral, sgeneral, nbond,co1,co2,ch1,ch2)
  
  return
end subroutine calc_coop

! calculate COOP/COHP
subroutine gfcoop(N,nspinblocks, nspinzmatrix, wk, GammaL, GammaR, gf, hgeneral, sgeneral, nbond, co1, co2, ch1,ch2)

  use mTypes

  implicit none
  integer N, nspinblocks, nbond, nspinzmatrix
  type(matrixTypeGeneral) :: gf, hgeneral(nspinzmatrix), sgeneral
  double precision  co1(nbond), co2(nbond), ch1(nbond), ch2(nbond), wk

  integer i,j, ii,jj,ind, ibond, nloop, ij, i0, j0, ish
  double precision, parameter ::  PI=3.141592653589793D0, RYDBERG=13.6056981D0
  double complex :: GammaL(N,N), GammaR(N,N)
  double complex, parameter :: ZZ=(0.d0, 0.d0), ZO=(1.d0,0.d0)

  if (nspinblocks .le. 2) then
      nloop =1
  else
      nloop =2
  endif


  ! get AL, AR
  CALL ZHEMM('R','U',N,N,ZO,GammaL+ZZ,N,gf%matdense%a,N, ZZ,GammaL,N)
  CALL ZGEMM('N','C',N,N,N,ZO,GammaL+ZZ,N,gf%matdense%a,N, ZZ,GammaL,N)
!
  CALL ZHEMM('R','U',N,N,ZO,GammaR+ZZ,N,gf%matdense%a,N, ZZ,GammaR,N)
  CALL ZGEMM('N','C',N,N,N,ZO,GammaR+ZZ,N,gf%matdense%a,N, ZZ,GammaR,N)

!
  do ibond =1, nbond
    do ij=1,nloop
       do i =1,coopinfo%no1(ibond)
          i0 =coopinfo%io1(i, ibond)
          ii =i0+(ij-1)*N/2
       ! For noncollinear spin, ignore cross elements
      ! do il=1,nloop
       !   ii =i0 +(il-1)*N/2
          do j=1,coopinfo%no2(ibond)
             j0 =coopinfo%io2(j,ibond)
             jj =j0 +(ij-1)*N/2
             call inddensetoindsparsegeneral(j0, i0, ind, sgeneral) 
             if (ind .ne. 0) then
                co1(ibond) =co1(ibond) +dble(GammaL(ii,jj)*sgeneral%matSparse%b(ind))*wk/(2*PI)/RYDBERG
                co2(ibond) =co2(ibond) +dble(GammaR(ii,jj)*sgeneral%matSparse%b(ind))*wk/(2*PI)/RYDBERG
             endif
      
             ish = ij*ij
             call inddensetoindsparsegeneral(j0, i0, ind, hgeneral(ish)) 
             if (ind .ne. 0) then
                ch1(ibond) =ch1(ibond) +dble(GammaL(ii,jj)*hgeneral(ish)%matSparse%b(ind))*wk/(2*PI)
                ch2(ibond) =ch2(ibond) +dble(GammaR(ii,jj)*hgeneral(ish)%matSparse%b(ind))*wk/(2*PI)
             endif
          enddo
       enddo
    enddo
  enddo

  return
end subroutine gfcoop

! output
subroutine output_coop(slabel, nheads, nep, nt, iv, nk, nspin, nbond, V, ep, coop1, coop2, cohp1, cohp2)
   integer, intent(in) :: iv, nheads, nep, nk, nspin, nbond, nt
   double precision, intent(in) :: ep(nep), V, coop1(nep, nspin, nbond), coop2(nep, nspin,nbond),cohp1(nep, nspin,nbond),cohp2(nep, nspin, nbond)
   character, intent(in) :: slabel*20

   character costr*8, coopfile*35, paste*35
   integer ibond, fid_coop, ispin, j, i, indt
   external paste

   do ibond=1,nbond
      write( costr, '(i7 )' ) iv
      costr = paste(costr, '.')
      coopfile = paste(costr,slabel)
      coopfile = paste(coopfile,'.b')
      write( costr, '(i7 )' ) ibond
      coopfile = paste(coopfile,TRIM(ADJUSTL(costr)))
      coopfile = paste(coopfile,'_')
      write( costr, '(i7 )' ) coopinfo%ia1(ibond)
      coopfile = paste(coopfile,TRIM(ADJUSTL(costr)))
      coopfile = paste(coopfile,'-')
      write( costr, '(i7 )' ) coopinfo%ia2(ibond)
      coopfile = paste(coopfile,TRIM(ADJUSTL(costr)))
      coopfile = paste(coopfile,'.COOP')

      call io_assign(fid_coop)

      OPEN(UNIT=fid_coop,FILE=coopfile,status='unknown',RECL=100000)
      WRITE(fid_coop,'(a6,f11.4,a14,i4)') '# V = ',V, '    k-points: ',nk
      WRITE(fid_coop,'(a15,3x)',ADVANCE='NO') '#    Energy    '
      DO ISPIN=1,NSPIN
        WRITE(fid_coop,'(2(a11, i1,4x))',ADVANCE='NO') '  COOP_A1_s',ispin, '  COOP_A2_s',ispin
      ENDDO
      WRITE(fid_coop,'(3x)',ADVANCE='NO') 
      DO ISPIN=1,NSPIN
        WRITE(fid_coop,'(2(a11, i1,4x))',ADVANCE='NO') '  COHP_A1_s',ispin, '  COHP_A2_s',ispin
      ENDDO

      WRITE(fid_coop,*) 
      WRITE(fid_coop,'(a, i5)') '# Bond : ', ibond 
      WRITE(fid_coop,'(a, i5, a)', advance='no') '# Atom 1 : ', coopinfo%ia1(ibond), ', '
      WRITE(fid_coop,'(3x,a)', advance='no') 'Orbital : '
      do j=1,coopinfo%no1(ibond)
         WRITE(fid_coop,'(x,i5)', advance='no') coopinfo%io1(j,ibond)
      enddo
      WRITE(fid_coop,*) 

      WRITE(fid_coop,'(a, i5,a)', advance='no') '# Atom 2 : ', coopinfo%ia2(ibond), ', '
      WRITE(fid_coop,'(3x,a)', advance='no') 'Orbital : '
      do j=1,coopinfo%no2(ibond)
         WRITE(fid_coop,'(x,i5)', advance='no') coopinfo%io2(j,ibond)
      enddo
      WRITE(fid_coop,*) 

      DO I=1,nep
      ! parallel segment for energy is one point 
        indt=MOD(i-1,nheads) * nt+ (i-1)/nheads + 1
        write(fid_coop, '(e15.5,3x)',ADVANCE='NO') ep(i)* 13.6056981D0
        DO ISPIN=1,NSPIN
          WRITE(fid_coop,'(2(e15.5,x))',ADVANCE='NO') coop1(indt,ispin,ibond), coop2(indt,ispin,ibond)
        ENDDO
        WRITE(fid_coop,'(3x)',ADVANCE='NO') 
        DO ISPIN=1,NSPIN
          WRITE(fid_coop,'(2(e15.5,x))',ADVANCE='NO') cohp1(indt,ispin,ibond), cohp2(indt,ispin,ibond)
        ENDDO
        write(fid_coop,*)
      ENDDO !energy
          
      call io_close(fid_coop)
   enddo ! ibond

   return
end subroutine output_coop

! free memory
subroutine free_coop
  deallocate(coopinfo%ia1, coopinfo%ia2)
  deallocate(coopinfo%no1, coopinfo%no2)
  deallocate(coopinfo%io1, coopinfo%io2)
  return
end subroutine free_coop

end module negfcoop

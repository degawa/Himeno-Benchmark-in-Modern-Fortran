!This program is Himeno benchmark problem written in Modern Fortran style.
!In this program, global variables are eliminated 
!and some variable names and subroutine names are refined.
!The execution performance is almost the same as the original version, 
!but it is about 2% slower in my environment (computer and compiler).
!
!For the original version of the Himeno benchmark, 
!please refer to the URLs below:
!http://accc.riken.jp/supercom/himenobmt/
!http://accc.riken.jp/en/supercom/himenobmt/
!
!This program is free open-source software distributed under LGPL version 2
!or any later version, inheriting the license of the original version of the Himeno benchmark.

program HimenoBMTxp_F90
    implicit none

    real(4),dimension(:,:,:)  ,allocatable :: p   !pressure
    real(4),dimension(:,:,:,:),allocatable :: a   !coefficient matrix for p(i+1), p(j+1), p(k+1), p(ijk)
    real(4),dimension(:,:,:,:),allocatable :: b   !coefficient matrix for cross derivative terms
    real(4),dimension(:,:,:,:),allocatable :: c   !coefficient matrix for p(i-1), p(j-1), p(k-1)
    real(4),dimension(:,:,:)  ,allocatable :: bnd !control variable for boundaries and objects
    real(4),dimension(:,:,:)  ,allocatable :: src !source term of Poisson equation
    real(4),dimension(:,:,:)  ,allocatable :: wrk !working area

    integer :: mimax,mjmax,mkmax
    integer ::  imax, jmax, kmax

    integer :: numItr
    integer :: count,countRate,countMax
    real(4) :: flop,mflops,score,error
    real(8) :: time_begin_s,time_end_s,time_elapsed_s,dt
    
    !Parameters related to performance measurments
    real(4),parameter :: FlopToMFlop = 1e-6
    real(4),parameter :: numFlopPerPoint = 34.0  ![operations]
    real(4),parameter :: MFlopsPenIII600 = 82.84 ![MFLOPS]

    call readGridParameter(mimax, mjmax, mkmax, imax, jmax, kmax)

    ! Initializing matrixes
    Initialize : block
        integer :: k, numPoints
        allocate(p  (mimax,mjmax,mkmax)  ,source = 0.0)
        allocate(a  (mimax,mjmax,mkmax,4),source = 0.0)
        allocate(b  (mimax,mjmax,mkmax,3),source = 0.0)
        allocate(c  (mimax,mjmax,mkmax,3),source = 0.0)
        allocate(bnd(mimax,mjmax,mkmax)  ,source = 0.0)
        allocate(src(mimax,mjmax,mkmax)  ,source = 0.0)
        allocate(wrk(mimax,mjmax,mkmax)  ,source = 0.0)

        a  (1:imax,1:jmax,1:kmax,1:3)=1.0
        a  (1:imax,1:jmax,1:kmax,4  )=1.0/6.0
        b  (1:imax,1:jmax,1:kmax, : )=0.0
        c  (1:imax,1:jmax,1:kmax, : )=1.0
        bnd(1:imax,1:jmax,1:kmax)    =1.0
        do k=1,kmax
            p(:,:,k)=sngl((k-1)**2)/sngl((kmax-1)**2)
        end do
        numPoints = (kmax-2)*(jmax-2)*(imax-2)
        flop=sngl(numPoints)*numFlopPerPoint
    end block Initialize

    print *," mimax=",mimax," mjmax=",mjmax," mkmax=",mkmax
    print *,"  imax=", imax,"  jmax=", jmax,"  kmax=",kmax

    call system_clock(count,countRate,countMax)
    dt= 1.0/dble(countRate)
    print "(a,e12.5)","Time measurement accuracy : ",dt

    !Rehearsal measurment to estimate the number of iterations
    Rehearsal : block
        numItr=3
        print *," Start rehearsal measurement process."
        print *," Measure the performance in 3 times."

        !Jacobi iteration
        time_begin_s = getCurrentTime()
        call jacobi(p,error, a,b,c,bnd,src,wrk,numItr)
        time_end_s   = getCurrentTime()
        
        time_elapsed_s = time_end_s - time_begin_s
        if(time_elapsed_s < dt) stop "error : execution time is not correct. The grid size may be too small."

        mflops = flop*FlopToMFlop / (time_elapsed_s/dble(numItr))
        print *,"  MFLOPS:",mflops,"  time(s):",time_elapsed_s,error
    end block Rehearsal
    !end Rehearsal measurment

    !Acatual measurment
    Actual : block
        !ExecTime specifys the measuring period in sec
        real(4),parameter :: ExecTime=60.0 !sec

        !set the number of Iterations so that the execution time is roughly ExecTime sec
        numItr=int( ExecTime/(time_elapsed_s/dble(numItr)) )
        print *,"Now, start the actual measurement process."
        print *,"The loop will be excuted in",numItr," times."
        print *,"This will take about one minute."
        print *,"Wait for a while."
        
        !Jacobi iteration
        time_begin_s = getCurrentTime()
        call jacobi(p,error, a,b,c,bnd,src,wrk,numItr)
        time_end_s   = getCurrentTime()

        time_elapsed_s = time_end_s - time_begin_s
        mflops = flop*FlopToMFlop / (time_elapsed_s/dble(numItr))
        score = mflops/MFlopsPenIII600
        
        print *," Loop executed for ",numItr," times"
        print *," Error :",error
        print *," MFLOPS:",mflops, "  time(s):",time_elapsed_s
        print *," Score based on Pentium III 600MHz :",score
    end block Actual

    deallocate(p  )
    deallocate(a  )
    deallocate(b  )
    deallocate(c  )
    deallocate(bnd)
    deallocate(src)
    deallocate(wrk)

    contains

    function getCurrentTime() result(currentTime)
        implicit none
        integer :: count,countRate,countMax
        real(8) :: currentTime
        
        call system_clock(count,countRate,countMax)
        currentTime= dble(count)/dble(countRate)
    end function getCurrentTime

    subroutine readGridParameter(mimax, mjmax, mkmax, imax, jmax, kmax)
        implicit none
        integer,intent(inout) :: mimax
        integer,intent(inout) :: mjmax
        integer,intent(inout) :: mkmax
        integer,intent(inout) ::  imax
        integer,intent(inout) ::  jmax
        integer,intent(inout) ::  kmax
        
        character(10) :: size
        
        print *,"Select Grid-size:"
        print *,"           XS (64x32x32)"
        print *,"           S  (128x64x64)"
        print *,"           M  (256x128x128)"
        print *,"           L  (512x256x256)"
        print *,"           XL (1024x512x512)"
        print "(A,$)"," Grid-size = "
        read(*,*) size
        
        call setGridSize(mimax,mjmax,mkmax,size)
        
        imax = mimax - 1
        jmax = mjmax - 1
        kmax = mkmax - 1
        
    end subroutine readGridParameter

    subroutine setGridSize(mimax,mjmax,mkmax,size)
        implicit none
        integer,intent(inout) :: mimax
        integer,intent(inout) :: mjmax
        integer,intent(inout) :: mkmax
        character(*),intent(in) :: size

        select case(size)
            case("XS","xs")
                mimax = 64 + 1
                mjmax = 32 + 1
                mkmax = 32 + 1
            case("S", "s")
                mimax = 128 + 1
                mjmax =  64 + 1
                mkmax =  64 + 1
            case("M", "m")
                mimax = 256 + 1
                mjmax = 128 + 1
                mkmax = 128 + 1
            case("L", "l")
                mimax = 512 + 1
                mjmax = 256 + 1
                mkmax = 256 + 1
            case("XL","xl")
                mimax = 1024 + 1
                mjmax =  512 + 1
                mkmax =  512 + 1
            case default
                stop "Unexpected GridSize"
        end select
    end subroutine setGridSize

    subroutine jacobi(p,error, a,b,c,bnd,src,wrk,numItr)
        implicit none
        real(4),intent(inout) :: error

        real(4),intent(inout) :: p  (:,:,:)    !pressure
        real(4),intent(in   ) :: a  (:,:,:,:)  !coefficient matrix for p(i+1), p(j+1), p(k+1), p(ijk)
        real(4),intent(in   ) :: b  (:,:,:,:)  !coefficient matrix for cross derivative term
        real(4),intent(in   ) :: c  (:,:,:,:)  !coefficient matrix for p(i+1), p(j+1), p(k+1)
        real(4),intent(in   ) :: bnd(:,:,:)    !control variable for boundaries and objects
        real(4),intent(in   ) :: src(:,:,:)    !source term of Poisson equation
        real(4),intent(inout) :: wrk(:,:,:)    !working area
        integer,intent(in   ) :: numItr
        
        integer :: loop
        integer :: i,j,k
        real(4) :: pnew,dp
        real(4),parameter :: rlx=0.8  !relaxation parameter
        integer,parameter :: x   = 1
        integer,parameter :: y   = 2
        integer,parameter :: z   = 3
        integer,parameter :: ctr = 4
        integer,parameter :: xy = 1
        integer,parameter :: yz = 2
        integer,parameter :: zx = 3
        
        !These variables are not necessary because the variables defined in the main routine can be referred.
        !But, I recommend to declare from the viewpoint of subroutine completeness.
        integer :: imax
        integer :: jmax
        integer :: kmax
        imax = ubound(p,x)-1
        jmax = ubound(p,y)-1
        kmax = ubound(p,z)-1

        do loop=1,numItr
            error= 0.0
            do k=2,kmax-1
                do j=2,jmax-1
                    do i=2,imax-1
                        pnew =  a(i,j,k,x )*  p(i+1,j  ,k  ) &
                               +a(i,j,k,y )*  p(i  ,j+1,k  ) &
                               +a(i,j,k,z )*  p(i  ,j  ,k+1) &
                               +b(i,j,k,xy)*( p(i+1,j+1,k  )-p(i+1,j-1,k  )  &
                                             -p(i-1,j+1,k  )+p(i-1,j-1,k  )) &
                               +b(i,j,k,yz)*( p(i  ,j+1,k+1)-p(i  ,j-1,k+1)  &
                                             -p(i  ,j+1,k-1)+p(i  ,j-1,k-1)) &
                               +b(i,j,k,zx)*( p(i+1,j  ,k+1)-p(i-1,j  ,k+1)  &
                                             -p(i+1,j  ,k-1)+p(i-1,j  ,k-1)) &
                               +c(i,j,k,x )*  p(i-1,j  ,k  ) &
                               +c(i,j,k,y )*  p(i  ,j-1,k  ) &
                               +c(i,j,k,z )*  p(i  ,j  ,k-1) &
                               +src(i,j,k)
                        
                        dp = (pnew*a(i,j,k,ctr) - p(i,j,k))*bnd(i,j,k)
                        error = error + dp*dp
                        wrk(i,j,k) = p(i,j,k) + rlx*dp
                    end do
                end do
            end do
            p(2:imax-1,2:jmax-1,2:kmax-1) = wrk(2:imax-1,2:jmax-1,2:kmax-1)
        end do
        !End of iteration
    end subroutine jacobi
end program HimenoBMTxp_F90

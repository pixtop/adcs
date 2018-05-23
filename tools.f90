  subroutine parse_args(n, m, version, imat, percentage, disp)

    implicit none
    integer          :: m, n, version, imat, disp, len
    real(kind(1.d0)) :: percentage

    character(len=20)  :: str
    integer :: i, cln

    ! default values
    version    = 0
    imat       = 1
    n          = -1
    m          = -1
    disp       = 0
    percentage = 0.7


    cln = command_argument_count()

    if (cln/2 .lt. 3) then
       call print_use()
       stop
    else
       do i=1, cln, 2
          call get_command_argument(i,value=str,length=len)

          select case(str(1:len))
          case('-n')
             call get_command_argument(i+1,value=str,length=len)
             read(str(1:len),*)n
          case('-m')
             call get_command_argument(i+1,value=str,length=len)
             read(str(1:len),*)m
          case('-imat')
             call get_command_argument(i+1,value=str,length=len)
             read(str(1:len),*)imat
          case('-v')
             call get_command_argument(i+1,value=str,length=len)
             read(str(1:len),*)version
          case('-per')
             call get_command_argument(i+1,value=str,length=len)
             read(str(1:len),*)percentage
          case('-disp')
             call get_command_argument(i+1,value=str,length=len)
             read(str(1:len),*)disp
          case default
             write(*,'("Wrong argument!")')
             write(*,'(" ")')
             call print_use()
             stop
          end select
       end do
    end if

    if(n .le. 0) then
       write(*,'(" ")')
       write(*,'("================================================='//&
         & '==================================================")')
       write(*,'(" ")')
       write(*,'("n is mandatory arguments and should be positive!")')
       write(*,'(" ")')
       call print_use()
       stop
    end if

    if((version .ge. 1).and.(m.lt.0)) then
       write(*,'(" ")')
       write(*,'("================================================='//&
         & '==================================================")')
       write(*,'(" ")')
       write(*,'("m is mandatory for solver version greater than 1")')
       write(*,'(" ")')
       call print_use()
       stop
    end if
    
    return
  end subroutine parse_args

  subroutine print_use()
    write(*,'("================================================='//&
         & '==================================================")')
    write(*,'(" ")')
    write(*,'("This program should be invoked like this:        '//&
         &'                                                                                 ")')
    write(*,'("./main args                                      '//&
         &'                                                                                 ")')
    write(*,'("where args can be                                '//&
         &'                                                                                 ")')
    write(*,'("-n n     : the size of the matrix A'//&
         &'                                                                                 ")')
    write(*,'("-m m     : the maximum number of eigenvalues we c'//&
         &'ompute (v=1) or the size of the subspace (v >= 2)                                ")')
    write(*,'("-imat i  : the matrix type 1, 2, 3 or 4 (default=1)'//&
         &'                                                                                 ")')
    write(*,'("-v v     : which solver version to use (default=0'//&
         &', i.e., LAPACK DGESVD)                                                           ")')
    write(*,'("           v = 0 : the LAPACK DGESVD")')
    write(*,'("           v = 1 : the power method with deflation")')
    write(*,'("           v = 2 : the iteration subspace method with a fixed number of eigenpairs to be computed")')
    write(*,'("           v = 3 : the iteration subspace method with a percentage of the trace")')
    write(*,*) 
    write(*,'("-per x.y : the amount of information to retain (d'//&
         &'efault=0.7, i.e., 70% of the trace)                                              ")')
    write(*,'("-disp d  : verbosity (default=0)")')
    write(*,*) 
    write(*,'("           d = 0, few information")')
    write(*,'("           d = 1, print eigenvalues")')
    write(*,'("           d = 2, print eigenvalues and accuracy")')
    write(*,'("                                                 '//&
         &'                                                                                 ")')
    write(*,'("Note that -n is mandatory and -m when (v>=1)            '//&
         &'                                                                                 ")')
    write(*,'("Example:                                         '//&
         &'                                                                                 ")')
    write(*,'(" ./main -n 50 -m 20 -per 0.1 -imat 2 -v 0 '//&
         &'-disp 1                                                                         ")')
    write(*,'("================================================='//&
         & '==================================================")')
  end subroutine print_use

  !!========================================================================
  !! orthogonalization using gram schmidt procedure
  !!========================================================================
  subroutine gram_schmidt(u_in, n, m, u_out)
    implicit none
    integer,          intent(in)                 :: n,m
    double precision, dimension(n,m),intent(in)  :: u_in
    double precision, dimension(n,m),intent(out) :: u_out    

    !! local variables
    integer :: i,j

    do i = 1, m
       u_out(:,i) = u_in(:,i)
    end do

    do j = 1, m
       if(j .gt. 1)then
          do i = 1, j-1
             u_out(:,j) = u_out(:,j) - dot_product(u_out(:,j), u_out(:,i))*u_out(:,i)
          end do
       end if
       u_out(:,j) = u_out(:,j)/sqrt(dot_product(u_out(:,j),u_out(:,j)))
    end do

  end subroutine gram_schmidt

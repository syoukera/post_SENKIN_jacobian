program calc_jacobian
      implicit none

      ! unit number
      integer, parameter :: unit_stdout = 6
      integer, parameter :: unit_cklink = 10
      integer, parameter :: unit_skdata = 11
      integer, parameter :: unit_jac    = 12

      ! chemkin work array
      integer len_int_cklen, len_real_cklen, len_char_cklen
      integer,      allocatable :: int_ckwk(:)
      real(8),      allocatable :: real_ckwk(:)
      character(6), allocatable :: char_ckwk(:)*16

      ! chemkin index
      integer mm, kk, ii, nfit

      ! simulation result
      real(8) time
      real(8) temperature_K
      real(8) pressure_atm
      real(8), allocatable :: x(:)
      real(8), allocatable :: y(:)

      integer i
      
      open(unit_cklink, form='unformatted', file='input/cklink')
      open(unit_skdata, form='formatted', file='input/skout_datasheet')
      open(unit_jac,    form='formatted', file='output/jacobian')

      !   ------- initialize chemkin work array ---------

      call cklen(unit_cklink, unit_stdout, len_int_cklen, len_real_cklen, len_char_cklen)

      allocate(int_ckwk(len_int_cklen))
      allocate(real_ckwk(len_real_cklen))
      allocate(char_ckwk(len_char_cklen))

      call ckinit(len_int_cklen, len_real_cklen, len_char_cklen, unit_cklink, &
                  unit_stdout, int_ckwk, real_ckwk, char_ckwk)

      call ckindx(int_ckwk, real_ckwk, mm, kk, ii, nfit)

      !   ------- load simulation result ---------

      allocate(x(kk), y(kk))

      ! read each line
      do 
            read(unit_skdata, *, end=999) time, pressure_atm, temperature_K,  &
                                        (x(i), i = 1, kk)
! C
! C     CONVERT X TO Y
! C
!           CALL CKXTY  (X, IWORK, RWORK, Y)
! C
! C     CALCULATE ENTHALPY
! C
!           CALL CKHBMS (T, Y, IWORK, RWORK, HBMS)   ! mass units
!           CALL CKHBML (T, X, IWORK, RWORK, HBML) ! molar units
! C
! C     CALCULATE ROP
! C
!           CALL CKQYP  (P, T, Y, IWORK, RWORK, Q)
! C
! C     PRINT OUT ENTHALPY
! C
            write(unit_jac, *) time
      enddo

999   continue

end program calc_jacobian
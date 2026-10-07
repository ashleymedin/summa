! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! **************************************************************************************************
! MODFLOW 6 parameters for calibration, as multipliers on the model's own input files.
!
! A parameter scales numbers in plain-text input files of the MODFLOW model directory, in one of
! two forms:
!   array  every value of an OPEN/CLOSE array file (K, K33, SY, ...), or only those whose cell holds
!          a given zone in a zone file of the same shape
!   list   one column of every data row in a list file, or in one block of it (UZF PACKAGEDATA's
!          VKS, a DRN conductance, ...)
! Each trial rewrites the files it names in its own copy of the model directory, always from the
! master copy, so no trial inherits another's values.  Several parameters may scale one file; their
! multipliers combine.  Nothing here calls libmf6, so it builds without MODFLOW.
! **************************************************************************************************
module mf6_parameters

  USE nr_type,    only: i4b, rkind, lgt
  USE data_types, only: mf6_param_info

  implicit none
  private

  public :: check_mf6_parameters
  public :: write_mf6_parameters
  public :: get_line, split_line, upper

  integer(i4b), parameter :: maxLine = 65536   ! longest input line read

contains

  ! **************************************************************************************************
  ! Refuse, before any trial runs, a MODFLOW parameter whose files, zone or column the model lacks.
  ! **************************************************************************************************
  subroutine check_mf6_parameters(params, model_dir, err, message)
    type(mf6_param_info), intent(in)  :: params(:)
    character(*),         intent(in)  :: model_dir    ! master MODFLOW model directory
    integer(i4b),         intent(out) :: err
    character(*),         intent(out) :: message
    real(rkind),  allocatable :: values(:)
    integer(i4b), allocatable :: zones(:)
    integer(i4b) :: i, j, k
    logical(lgt) :: exists
    character(len=1024) :: cmessage

    err=0
    message='check_mf6_parameters/'

    do i=1,size(params)
      if(.not.allocated(params(i)%files))then
        message=trim(message)//'MODFLOW parameter "'//trim(params(i)%name)//'" names no files'; err=20; return
      endif
      if(size(params(i)%files) == 0)then
        message=trim(message)//'MODFLOW parameter "'//trim(params(i)%name)//'" names no files'; err=20; return
      endif
    enddo

    do i=1,size(params)
      associate(p => params(i))
      if(p%column > 0 .and. len_trim(p%zone_file) > 0)then
        message=trim(message)//'MODFLOW parameter "'//trim(p%name)//'" scales a list column, which takes no zone_file'
        err=20; return
      endif
      if(len_trim(p%zone_file) > 0 .neqv. p%zone >= 0)then
        message=trim(message)//'MODFLOW parameter "'//trim(p%name)//'" needs both zone_file and zone, or neither'
        err=20; return
      endif
      if(len_trim(p%block) > 0 .and. p%column <= 0)then
        message=trim(message)//'MODFLOW parameter "'//trim(p%name)//'" names a block but no column'
        err=20; return
      endif

      do j=1,size(p%files)
        inquire(file=model_dir//'/'//trim(p%files(j)), exist=exists)
        if(.not.exists)then
          message=trim(message)//'MODFLOW parameter "'//trim(p%name)//'": no file '//trim(p%files(j))// &
                  ' in '//model_dir; err=20; return
        endif

        ! a file is scaled as an array or as a list, not both
        do k=1,size(params)
          if(any(params(k)%files == p%files(j)) .and. (params(k)%column > 0 .neqv. p%column > 0))then
            message=trim(message)//'MODFLOW parameters "'//trim(p%name)//'" and "'//trim(params(k)%name)// &
                    '" scale '//trim(p%files(j))//' as an array and as a list'; err=20; return
          endif
        enddo

        if(p%column > 0)then
          ! every data row of the block reaches the column with a number
          call scale_list(model_dir//'/'//trim(p%files(j)), '', [p], [1._rkind], err, cmessage)
          if(err/=0)then; message=trim(message)//'MODFLOW parameter "'//trim(p%name)//'": '//trim(cmessage); return; endif
        else
          call read_array(model_dir//'/'//trim(p%files(j)), values, err, cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          if(len_trim(p%zone_file) > 0)then
            call read_zones(model_dir//'/'//trim(p%zone_file), zones, err, cmessage)
            if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
            if(size(zones) /= size(values))then
              write(message,'(a,i0,a,i0,a)') trim(message)//'MODFLOW parameter "'//trim(p%name)//'": '// &
                    trim(p%zone_file)//' holds ',size(zones),' zones and '//trim(p%files(j))//' ',size(values),' values'
              err=20; return
            endif
            if(.not.any(zones == p%zone))then
              write(message,'(a,i0)') trim(message)//'MODFLOW parameter "'//trim(p%name)//'": '// &
                    trim(p%zone_file)//' has no cell in zone ',p%zone
              err=20; return
            endif
          endif
        endif
      enddo
      end associate
    enddo

  end subroutine check_mf6_parameters

  ! **************************************************************************************************
  ! Write one trial's MODFLOW parameters into run_dir, scaling each named file from model_dir.
  ! multiplier(i) belongs to params(i); a multiplier of one leaves the master's values.
  ! **************************************************************************************************
  subroutine write_mf6_parameters(params, multiplier, model_dir, run_dir, err, message)
    type(mf6_param_info), intent(in)  :: params(:)
    real(rkind),          intent(in)  :: multiplier(:)
    character(*),         intent(in)  :: model_dir    ! master MODFLOW model directory
    character(*),         intent(in)  :: run_dir      ! this instance's copy, rewritten
    integer(i4b),         intent(out) :: err
    character(*),         intent(out) :: message
    character(len=256), allocatable :: files(:)
    real(rkind),        allocatable :: values(:), scale(:)
    integer(i4b),       allocatable :: zones(:)
    logical(lgt),       allocatable :: onFile(:)
    integer(i4b) :: i, j, iFile
    character(len=1024) :: cmessage

    err=0
    message='write_mf6_parameters/'

    ! each file once, however many parameters scale it
    allocate(files(0))
    do i=1,size(params)
      do j=1,size(params(i)%files)
        if(.not.any(files == params(i)%files(j))) files=[files, params(i)%files(j)]
      enddo
    enddo

    do iFile=1,size(files)
      onFile=[(any(params(i)%files == files(iFile)), i=1,size(params))]

      if(any(onFile .and. params(:)%column > 0))then
        call scale_list(model_dir//'/'//trim(files(iFile)), run_dir//'/'//trim(files(iFile)), &
                        pack(params, onFile), pack(multiplier, onFile), err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        cycle
      endif

      call read_array(model_dir//'/'//trim(files(iFile)), values, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      allocate(scale(size(values))); scale=1._rkind
      do i=1,size(params)
        if(.not.onFile(i)) cycle
        if(len_trim(params(i)%zone_file) == 0)then
          scale=scale*multiplier(i)
        else
          call read_zones(model_dir//'/'//trim(params(i)%zone_file), zones, err, cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          where(zones == params(i)%zone) scale=scale*multiplier(i)
        endif
      enddo
      call write_array(run_dir//'/'//trim(files(iFile)), values*scale, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      deallocate(scale)
    enddo

  end subroutine write_mf6_parameters

  ! **************************************************************************************************
  ! Scale the columns params name in the data rows of a list file, writing the result to out_path;
  ! a blank out_path only checks that every row the parameters touch has a number in their column.
  ! **************************************************************************************************
  subroutine scale_list(in_path, out_path, params, multiplier, err, message)
    character(*),         intent(in)  :: in_path
    character(*),         intent(in)  :: out_path
    type(mf6_param_info), intent(in)  :: params(:)
    real(rkind),          intent(in)  :: multiplier(:)
    integer(i4b),         intent(out) :: err
    character(*),         intent(out) :: message
    character(len=:), allocatable :: line, word, block
    character(len=64), allocatable :: tokens(:)
    character(len=32) :: cval
    real(rkind)  :: val
    integer(i4b) :: iu, ou, ios, i, nRow, iLine
    integer(i4b), allocatable :: nHit(:)
    logical(lgt) :: writing

    err=0; message=''
    writing=len_trim(out_path) > 0
    allocate(nHit(size(params))); nHit=0

    open(newunit=iu, file=in_path, status='old', action='read', iostat=ios)
    if(ios/=0)then; message='cannot open '//in_path; err=20; return; endif
    if(writing)then
      open(newunit=ou, file=out_path, status='replace', action='write', iostat=ios)
      if(ios/=0)then; message='cannot write '//out_path; err=20; close(iu); return; endif
    endif

    block=''; iLine=0
    do
      call get_line(iu, line, ios)
      if(ios/=0) exit
      iLine=iLine+1
      call split_line(line, tokens)

      ! comments, blank lines and block delimiters pass through
      if(size(tokens) == 0)then
        if(writing) write(ou,'(a)') line
        cycle
      endif
      word=upper(tokens(1))
      if(word(1:1) == '#' .or. word(1:1) == '!')then
        if(writing) write(ou,'(a)') line
        cycle
      endif
      if(word == 'BEGIN' .and. size(tokens) > 1)then
        block=upper(tokens(2))
        if(writing) write(ou,'(a)') line
        cycle
      endif
      if(word == 'END')then
        block=''
        if(writing) write(ou,'(a)') line
        cycle
      endif

      nRow=0
      do i=1,size(params)
        if(params(i)%column <= 0) cycle
        if(len_trim(params(i)%block) > 0 .and. upper(params(i)%block) /= block) cycle
        if(len_trim(params(i)%block) == 0 .and. len(block) > 0) cycle
        if(size(tokens) < params(i)%column)then
          write(message,'(a,i0,a,i0)') trim(in_path)//' line ',iLine,' has no column ',params(i)%column
          err=20; exit
        endif
        read(tokens(params(i)%column),*,iostat=ios) val
        if(ios/=0)then
          write(message,'(a,i0,a,i0,a)') trim(in_path)//' line ',iLine,' column ',params(i)%column, &
                ' is not a number: '//trim(tokens(params(i)%column))
          err=20; exit
        endif
        write(cval,'(es24.16e3)') val*multiplier(i)
        tokens(params(i)%column)=adjustl(cval)
        nRow=nRow+1; nHit(i)=nHit(i)+1
      enddo
      if(err/=0) exit

      if(writing)then
        if(nRow == 0)then
          write(ou,'(a)') line
        else
          write(ou,'(*(1x,a))') (trim(tokens(i)), i=1,size(tokens))
        endif
      endif
    enddo
    close(iu)
    if(writing) close(ou)
    if(err/=0) return

    ! a block or column that matched nothing is a misnamed one
    do i=1,size(params)
      if(params(i)%column > 0 .and. nHit(i) == 0)then
        message='no data rows'
        if(len_trim(params(i)%block) > 0) message=trim(message)//' in block '//trim(params(i)%block)
        message=trim(message)//' of '//in_path
        err=20; return
      endif
    enddo

  end subroutine scale_list

  ! **************************************************************************************************
  ! Every number in a free-format array file, expanding list-directed repeat counts (n*value).
  ! **************************************************************************************************
  subroutine read_array(path, values, err, message)
    character(*),             intent(in)  :: path
    real(rkind), allocatable, intent(out) :: values(:)
    integer(i4b),             intent(out) :: err
    character(*),             intent(out) :: message
    character(len=64), allocatable :: tokens(:)
    character(len=:),  allocatable :: line
    real(rkind)  :: val
    integer(i4b) :: iu, ios, i, ix, nRep, n

    err=0; message=''
    allocate(values(1024)); n=0
    open(newunit=iu, file=path, status='old', action='read', iostat=ios)
    if(ios/=0)then; message='cannot open '//path; err=20; return; endif
    do
      call get_line(iu, line, ios)
      if(ios/=0) exit
      call split_line(line, tokens)
      do i=1,size(tokens)
        nRep=1
        ix=index(tokens(i),'*')
        if(ix > 0)then
          read(tokens(i)(1:ix-1),*,iostat=ios) nRep
          if(ios==0) read(tokens(i)(ix+1:),*,iostat=ios) val
        else
          read(tokens(i),*,iostat=ios) val
        endif
        if(ios/=0)then
          message=path//' holds something that is not a number: '//trim(tokens(i))//'; only text arrays can be scaled'
          err=20; close(iu); return
        endif
        do while(n+nRep > size(values))
          values=[values, values]
        enddo
        values(n+1:n+nRep)=val
        n=n+nRep
      enddo
    enddo
    close(iu)
    values=values(1:n)

  end subroutine read_array

  ! **************************************************************************************************
  ! The zone of every cell in a zone file of the same shape as the arrays it selects from.
  ! **************************************************************************************************
  subroutine read_zones(path, zones, err, message)
    character(*),              intent(in)  :: path
    integer(i4b), allocatable, intent(out) :: zones(:)
    integer(i4b),              intent(out) :: err
    character(*),              intent(out) :: message
    real(rkind), allocatable :: values(:)

    call read_array(path, values, err, message)
    if(err/=0) return
    zones=nint(values)

  end subroutine read_zones

  ! **************************************************************************************************
  ! Write a free-format array file, at full double precision so a multiplier of one leaves it unchanged.
  ! **************************************************************************************************
  subroutine write_array(path, values, err, message)
    character(*), intent(in)  :: path
    real(rkind),  intent(in)  :: values(:)
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: ou, ios

    err=0; message=''
    open(newunit=ou, file=path, status='replace', action='write', iostat=ios)
    if(ios/=0)then; message='cannot write '//path; err=20; return; endif
    write(ou,'(10(1x,es24.16e3))') values
    close(ou)

  end subroutine write_array

  ! **************************************************************************************************
  ! One line of a file, of any length.
  ! **************************************************************************************************
  subroutine get_line(iu, line, ios)
    integer(i4b),                  intent(in)  :: iu
    character(len=:), allocatable, intent(out) :: line
    integer(i4b),                  intent(out) :: ios
    character(len=1024) :: buf
    integer(i4b) :: nRead

    line=''
    do
      read(iu,'(a)',advance='no',size=nRead,iostat=ios) buf
      line=line//buf(1:nRead)
      if(is_iostat_eor(ios))then; ios=0; return; endif
      if(ios/=0) return
      if(len(line) > maxLine)then; ios=-1; return; endif
    enddo

  end subroutine get_line

  ! **************************************************************************************************
  ! The words of a line, separated by blanks, tabs or commas.
  ! **************************************************************************************************
  subroutine split_line(line, tokens)
    character(*),                   intent(in)  :: line
    character(len=64), allocatable, intent(out) :: tokens(:)
    integer(i4b) :: i, i0
    logical(lgt) :: sep

    allocate(tokens(0))
    i0=0
    do i=1,len(line)+1
      sep=.true.
      if(i <= len(line)) sep=index(' ,'//achar(9), line(i:i)) > 0
      if(sep)then
        if(i0 > 0) tokens=[tokens, line(i0:i-1)]
        i0=0
      else if(i0 == 0)then
        i0=i
      endif
    enddo

  end subroutine split_line

  ! **************************************************************************************************
  ! Upper case, for MODFLOW's case-insensitive keywords.
  ! **************************************************************************************************
  pure function upper(s) result(u)
    character(*), intent(in) :: s
    character(len=len_trim(adjustl(s))) :: u
    integer(i4b) :: i, ic
    u=trim(adjustl(s))
    do i=1,len(u)
      ic=iachar(u(i:i))
      if(ic >= iachar('a') .and. ic <= iachar('z')) u(i:i)=achar(ic-32)
    enddo
  end function upper

end module mf6_parameters

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
! MODFLOW 6 observations as calibration targets.
!
! A target with variable = "modflow_obs" scores one HEAD observation of a GWF model's OBS6 package,
! named by obs_name.  find_mf6_obs finds it in the model directory's input files, and read_mf6_obs
! reads its values back from the CONTINUOUS FILEOUT CSV MODFLOW writes, one row per time step, time
! in the model's time unit.  Nothing here calls libmf6, so it builds without MODFLOW.
! **************************************************************************************************
module mf6_observations

  USE nr_type,        only: i4b, rkind, lgt
  USE mf6_parameters, only: get_line, split_line, upper

  implicit none
  private

  public :: is_mf6_obs
  public :: find_mf6_obs
  public :: read_mf6_obs

contains

  ! **************************************************************************************************
  ! Report whether a target's variable names a MODFLOW 6 observation.
  ! **************************************************************************************************
  pure function is_mf6_obs(name) result(isObs)
    character(*), intent(in) :: name
    logical(lgt)             :: isObs
    isObs = trim(name) == 'modflow_obs'
  end function is_mf6_obs

  ! **************************************************************************************************
  ! Find HEAD observation obs_name in the OBS6 files of the GWF models of mfsim.nam in model_dir.
  ! csv_file is the CONTINUOUS FILEOUT holding it, relative to the model directory.
  ! **************************************************************************************************
  subroutine find_mf6_obs(model_dir, obs_name, csv_file, err, message)
    character(*),                  intent(in)  :: model_dir
    character(*),                  intent(in)  :: obs_name
    character(len=:), allocatable, intent(out) :: csv_file
    integer(i4b),                  intent(out) :: err
    character(*),                  intent(out) :: message
    character(len=256), allocatable :: gwfFiles(:), obsFiles(:), words(:)
    character(len=64),  allocatable :: tokens(:)
    character(len=:),   allocatable :: block, fileout
    integer(i4b) :: iGwf, iObs, iu, iTok
    logical(lgt) :: isHead
    character(len=1024) :: cmessage

    err=0
    message='find_mf6_obs/'

    call block_entries(model_dir//'/mfsim.nam', 'MODELS', 'GWF6', gwfFiles, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    allocate(obsFiles(0))
    do iGwf=1,size(gwfFiles)
      call block_entries(model_dir//'/'//trim(gwfFiles(iGwf)), 'PACKAGES', 'OBS6', words, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      obsFiles=[obsFiles, words]
    enddo
    if(size(obsFiles) == 0)then
      message=trim(message)//'no GWF model in '//trim(model_dir)//'/mfsim.nam lists an OBS6 file'
      err=20; return
    endif

    do iObs=1,size(obsFiles)
      call open_input(model_dir//'/'//trim(obsFiles(iObs)), iu, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      block=''; fileout=''
      do
        call next_tokens(iu, tokens, err)
        if(err/=0)then; err=0; exit; endif
        if(upper(tokens(1)) == 'BEGIN' .and. size(tokens) >= 2)then
          block=upper(tokens(2))
          fileout=''
          if(block == 'CONTINUOUS' .and. size(tokens) >= 4) fileout=unquote(tokens(4))
          do iTok=5,size(tokens)
            if(upper(tokens(iTok)) == 'BINARY') fileout=''
          enddo
        else if(upper(tokens(1)) == 'END')then
          block=''
        else if(block == 'CONTINUOUS' .and. upper(tokens(1)) == upper(obs_name))then
          close(iu)
          if(size(tokens) < 2)then
            isHead=.false.
          else
            isHead = upper(tokens(2)) == 'HEAD'
          endif
          if(.not.isHead)then
            message=trim(message)//'MODFLOW observation "'//trim(obs_name)//'" in '//trim(obsFiles(iObs))// &
                    ' is not a HEAD observation'
            err=20; return
          endif
          if(len(fileout) == 0)then
            message=trim(message)//'MODFLOW observation "'//trim(obs_name)//'" in '//trim(obsFiles(iObs))// &
                    ' is not in a CONTINUOUS FILEOUT text file'
            err=20; return
          endif
          csv_file=fileout
          return
        endif
      enddo
      close(iu)
    enddo

    message=trim(message)//'no OBS6 file of the model in '//trim(model_dir)//' names observation "'//trim(obs_name)//'"'
    err=20

  end subroutine find_mf6_obs

  ! **************************************************************************************************
  ! Read the column of observation obs_name from a MODFLOW 6 observation CSV: times and values, row by row.
  ! **************************************************************************************************
  subroutine read_mf6_obs(csv_path, obs_name, times, values, err, message)
    character(*),             intent(in)  :: csv_path
    character(*),             intent(in)  :: obs_name
    real(rkind), allocatable, intent(out) :: times(:), values(:)
    integer(i4b),             intent(out) :: err
    character(*),             intent(out) :: message
    character(len=:),  allocatable :: line
    character(len=64), allocatable :: tokens(:)
    integer(i4b) :: iu, ios, iCol, nRow, iRow
    character(len=1024) :: cmessage

    err=0
    message='read_mf6_obs/'

    call open_input(csv_path, iu, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call get_line(iu, line, ios)
    if(ios/=0)then
      message=trim(message)//trim(csv_path)//' has no header'; err=20; close(iu); return
    endif
    call split_line(line, tokens)
    iCol=0
    do iRow=2,size(tokens)
      if(upper(tokens(iRow)) == upper(obs_name)) iCol=iRow
    enddo
    if(iCol == 0)then
      message=trim(message)//trim(csv_path)//' has no column "'//trim(obs_name)//'"'; err=20; close(iu); return
    endif

    nRow=0
    do
      call get_line(iu, line, ios)
      if(ios/=0) exit
      if(len_trim(line) > 0) nRow=nRow+1
    enddo
    rewind(iu)
    call get_line(iu, line, ios)

    allocate(times(nRow), values(nRow))
    iRow=0
    do while(iRow < nRow)
      call get_line(iu, line, ios)
      if(ios/=0) exit
      if(len_trim(line) == 0) cycle
      iRow=iRow+1
      call split_line(line, tokens)
      if(size(tokens) < iCol) exit
      read(tokens(1),*,iostat=ios) times(iRow)
      if(ios==0) read(tokens(iCol),*,iostat=ios) values(iRow)
      if(ios/=0) exit
    enddo
    close(iu)
    if(iRow < nRow .or. ios/=0)then
      write(message,'(a,i0)') trim(message)//trim(csv_path)//' cannot be read at row ', iRow+1
      err=20; return
    endif

  end subroutine read_mf6_obs

  ! **************************************************************************************************
  ! The second word of every line starting with key in block blockName of a MODFLOW 6 input file.
  ! **************************************************************************************************
  subroutine block_entries(path, blockName, key, entries, err, message)
    character(*),                    intent(in)  :: path
    character(*),                    intent(in)  :: blockName
    character(*),                    intent(in)  :: key
    character(len=256), allocatable, intent(out) :: entries(:)
    integer(i4b),                    intent(out) :: err
    character(*),                    intent(out) :: message
    character(len=64), allocatable :: tokens(:)
    logical(lgt) :: inBlock
    integer(i4b) :: iu

    allocate(entries(0))
    call open_input(path, iu, err, message)
    if(err/=0) return
    inBlock=.false.
    do
      call next_tokens(iu, tokens, err)
      if(err/=0)then; err=0; exit; endif
      if(size(tokens) >= 2 .and. upper(tokens(1)) == 'BEGIN')then
        inBlock = upper(tokens(2)) == upper(blockName)
      else if(upper(tokens(1)) == 'END')then
        inBlock=.false.
      else if(inBlock .and. size(tokens) >= 2 .and. upper(tokens(1)) == upper(key))then
        entries=[character(len=256) :: entries, unquote(tokens(2))]
      endif
    enddo
    close(iu)

  end subroutine block_entries

  ! **************************************************************************************************
  ! The words of the next line that is neither blank nor a comment; err /= 0 at the end of the file.
  ! **************************************************************************************************
  subroutine next_tokens(iu, tokens, err)
    integer(i4b),                   intent(in)  :: iu
    character(len=64), allocatable, intent(out) :: tokens(:)
    integer(i4b),                   intent(out) :: err
    character(len=:), allocatable :: line

    do
      call get_line(iu, line, err)
      if(err/=0) return
      line=adjustl(line)
      if(len_trim(line) == 0) cycle
      if(line(1:1) == '#' .or. line(1:1) == '!') cycle
      call split_line(trim(line), tokens)
      if(size(tokens) > 0) return
    enddo

  end subroutine next_tokens

  ! **************************************************************************************************
  ! Open an existing text file for reading.
  ! **************************************************************************************************
  subroutine open_input(path, iu, err, message)
    character(*), intent(in)  :: path
    integer(i4b), intent(out) :: iu
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    message=''
    open(newunit=iu, file=trim(path), status='old', action='read', iostat=err)
    if(err/=0)then; message='cannot open '//trim(path); err=20; endif

  end subroutine open_input

  ! **************************************************************************************************
  ! A file name without the quotes MODFLOW allows around it.
  ! **************************************************************************************************
  pure function unquote(s) result(u)
    character(*), intent(in)      :: s
    character(len=:), allocatable :: u
    u=trim(adjustl(s))
    if(len(u) >= 2)then
      if((u(1:1) == "'" .or. u(1:1) == '"') .and. u(len(u):len(u)) == u(1:1)) u=u(2:len(u)-1)
    endif
  end function unquote

end module mf6_observations

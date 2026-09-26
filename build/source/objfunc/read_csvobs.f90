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
! Observations held in a dated CSV file
!
! Satellite products are distributed as dated columns far more often than as NetCDF: the GRACE
! terrestrial-water-storage anomalies this was written for arrive as one row per month, with a
! column per processing centre, and an empty field wherever a month has no solution.
!
! The file is read by column name, so one file can serve several targets - the JPL, CSR and GSFC
! columns of the same GRACE record, say - and the date column is whichever column comes first,
! named or not, as those products are written with an unnamed index column.
! **************************************************************************************************
module read_csvobs_module

  USE nr_type, only: i4b,rkind,lgt

  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

  implicit none
  private

  public :: read_csv_observations

  integer(i4b), parameter :: maxLine = 4096   ! longest CSV line this reads

contains

  ! **************************************************************************************************
  ! Read one named column of a dated CSV file.
  !
  ! The first column holds the date, as YYYY-MM-DD; the requested column holds the values.  Times are
  ! returned as days since 1900-01-01, which is what the time alignment understands, and an empty
  ! field becomes a quiet NaN, which is what the metrics understand as missing.
  ! **************************************************************************************************
  subroutine read_csv_observations(obs_path,obs_file,vname_obs,timeObs,valObs,timeUnits,valUnits,err,message)
    character(*), intent(in)                    :: obs_path      ! path to the observation file
    character(*), intent(in)                    :: obs_file      ! observation file
    character(*), intent(in)                    :: vname_obs     ! name of the column to read
    real(rkind), allocatable, intent(out)       :: timeObs(:)    ! observation time coordinate
    real(rkind), allocatable, intent(out)       :: valObs(:)     ! observed values
    character(len=:), allocatable, intent(out)  :: timeUnits     ! units and reference time
    character(len=:), allocatable, intent(out)  :: valUnits      ! units of the observed values
    integer(i4b), intent(out)                   :: err           ! error code
    character(*), intent(out)                   :: message       ! error message
    character(len=maxLine) :: line
    character(len=len(vname_obs)) :: wanted
    integer(i4b) :: unt
    integer(i4b) :: iostat
    integer(i4b) :: iCol,nCol
    integer(i4b) :: iRow,nRow
    integer(i4b) :: year,month,day
    logical(lgt) :: file_exists
    character(len=256) :: field
    character(len=256) :: cmessage

    err=0
    message='read_csv_observations/'

    ! the units of a CSV column are not carried by the file.  The caller states them on the target,
    ! and an empty string tells the alignment there is nothing to check.
    timeUnits='days since 1900-01-01'
    valUnits=''

    inquire(file=trim(obs_path)//trim(obs_file),exist=file_exists)
    if(.not.file_exists)then
      message=trim(message)//'observation file does not exist: '//trim(obs_path)//trim(obs_file)
      err=20; return
    endif

    open(newunit=unt,file=trim(obs_path)//trim(obs_file),status='old',action='read',iostat=iostat)
    if(iostat/=0)then
      message=trim(message)//'unable to open observation file: '//trim(obs_path)//trim(obs_file)
      err=20; return
    endif

    ! the header names the columns; find the one asked for
    read(unt,'(A)',iostat=iostat) line
    if(iostat/=0)then
      message=trim(message)//'observation file is empty: '//trim(obs_path)//trim(obs_file)
      err=20; close(unt); return
    endif
    wanted=adjustl(vname_obs)
    iCol=0
    nCol=count_fields(line)
    call find_column(line,trim(wanted),iCol)
    if(iCol < 2)then
      message=trim(message)//'column "'//trim(vname_obs)//'" is not in '//trim(obs_file)// &
              ' (the first column holds the date)'
      err=20; close(unt); return
    endif

    ! count the rows, then read them
    nRow=0
    do
      read(unt,'(A)',iostat=iostat) line
      if(iostat/=0) exit
      if(len_trim(line)==0) cycle
      nRow=nRow+1
    enddo
    if(nRow < 1)then
      message=trim(message)//'observation file holds no rows: '//trim(obs_path)//trim(obs_file)
      err=20; close(unt); return
    endif

    allocate(timeObs(nRow),valObs(nRow),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the observed series'
      close(unt); return
    endif

    rewind(unt)
    read(unt,'(A)',iostat=iostat) line     ! the header again
    iRow=0
    do
      read(unt,'(A)',iostat=iostat) line
      if(iostat/=0) exit
      if(len_trim(line)==0) cycle
      iRow=iRow+1

      ! the date, from the first column
      call get_field(line,1,field)
      call parse_date(trim(field),year,month,day,err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)//' in '//trim(obs_file)
        close(unt); return
      endif
      timeObs(iRow)=real(days_since_1900(year,month,day),rkind)

      ! the value, from the column asked for.  A month with no solution is left empty in these
      ! products, and reads as missing rather than as a number.
      call get_field(line,iCol,field)
      if(len_trim(field)==0)then
        valObs(iRow)=ieee_value(0._rkind,ieee_quiet_nan)
      else
        read(field,*,iostat=iostat) valObs(iRow)
        if(iostat/=0) valObs(iRow)=ieee_value(0._rkind,ieee_quiet_nan)
      endif
    enddo

    close(unt)

  end subroutine read_csv_observations

  ! **************************************************************************************************
  ! Count the comma-separated fields on a line.
  ! **************************************************************************************************
  pure function count_fields(line) result(nField)
    character(*), intent(in) :: line
    integer(i4b)             :: nField
    integer(i4b) :: i

    nField=1
    do i=1,len_trim(line)
      if(line(i:i)==',') nField=nField+1
    enddo

  end function count_fields

  ! **************************************************************************************************
  ! Return the n-th comma-separated field of a line, stripped of blanks and quotes.
  ! **************************************************************************************************
  subroutine get_field(line,iField,field)
    character(*), intent(in)  :: line
    integer(i4b), intent(in)  :: iField
    character(*), intent(out) :: field
    integer(i4b) :: i,iStart,iCount

    field=''
    iCount=1
    iStart=1
    do i=1,len_trim(line)+1
      if(i==len_trim(line)+1 .or. line(i:i)==',')then
        if(iCount==iField)then
          if(i > iStart) field=adjustl(line(iStart:i-1))
          call strip_quotes(field)
          return
        endif
        iCount=iCount+1
        iStart=i+1
      endif
    enddo

  end subroutine get_field

  ! **************************************************************************************************
  ! Find the column with a given header name, counting from one.
  ! **************************************************************************************************
  subroutine find_column(header,name,iCol)
    character(*), intent(in)  :: header
    character(*), intent(in)  :: name
    integer(i4b), intent(out) :: iCol
    character(len=256) :: field
    integer(i4b) :: i

    iCol=0
    do i=1,count_fields(header)
      call get_field(header,i,field)
      if(trim(field)==trim(name))then
        iCol=i
        return
      endif
    enddo

  end subroutine find_column

  ! **************************************************************************************************
  ! Remove surrounding quotes from a field.
  ! **************************************************************************************************
  subroutine strip_quotes(field)
    character(*), intent(inout) :: field
    integer(i4b) :: n

    n=len_trim(field)
    if(n < 2) return
    if((field(1:1)=='"' .and. field(n:n)=='"') .or. &
       (field(1:1)=="'" .and. field(n:n)=="'")) field=field(2:n-1)

  end subroutine strip_quotes

  ! **************************************************************************************************
  ! Parse a YYYY-MM-DD date, with or without a trailing time.
  ! **************************************************************************************************
  subroutine parse_date(text,year,month,day,err,message)
    character(*), intent(in)  :: text
    integer(i4b), intent(out) :: year,month,day
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: iostat
    integer(i4b) :: i
    character(len=64) :: work

    err=0
    message='parse_date/'
    year=0; month=0; day=0

    ! the delimiters become blanks, so an unpadded date reads as readily as a padded one and a
    ! trailing time of day is simply ignored
    work=adjustl(text)
    do i=1,len(work)
      select case(work(i:i))
        case ('-',':','T','t','/'); work(i:i)=' '
      end select
    enddo

    read(work,*,iostat=iostat) year,month,day
    if(iostat/=0 .or. month < 1 .or. month > 12 .or. day < 1 .or. day > 31)then
      message=trim(message)//'"'//trim(text)//'" is not a YYYY-MM-DD date'
      err=20; return
    endif

  end subroutine parse_date

  ! **************************************************************************************************
  ! Days from 1900-01-01 to a calendar date, by the usual civil-to-day-number identity.
  ! **************************************************************************************************
  pure function days_since_1900(year,month,day) result(days)
    integer(i4b), intent(in) :: year,month,day
    integer(i4b)             :: days

    days=day_number(year,month,day)-day_number(1900,1,1)

  end function days_since_1900

  ! **************************************************************************************************
  ! Day number of a proleptic Gregorian date, counted from an arbitrary fixed epoch.
  ! **************************************************************************************************
  pure function day_number(year,month,day) result(n)
    integer(i4b), intent(in) :: year,month,day
    integer(i4b)             :: n
    integer(i4b) :: y,m,era,yoe,doy,doe

    y=year
    m=month
    if(m <= 2) y=y-1
    if(y >= 0)then
      era=y/400
    else
      era=(y-399)/400
    endif
    yoe=y-era*400
    if(m > 2)then
      doy=(153*(m-3)+2)/5+day-1
    else
      doy=(153*(m+9)+2)/5+day-1
    endif
    doe=yoe*365+yoe/4-yoe/100+doy
    n=era*146097+doe

  end function day_number

end module read_csvobs_module

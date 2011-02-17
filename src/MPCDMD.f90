module MPCDMD
  implicit none
contains

  subroutine MPCDMD_info(short)
    logical, intent(in), optional :: short
    logical :: short_var
    include 'MPCDMD_version.h'

    if (present(short)) then
       short_var = short
    else
       short_var = .false.
    end if
    
    if (short_var) then
       write(*,*) 'MPCDMD> Version/date : ', trim(adjustl(MPCDMD_version_date))
    else
       write(*,*) 'MPCDMD> MPCDMD package'
       write(*,*) 'MPCDMD> (C) 2011 P. de Buyl'
       write(*,*) 'MPCDMD> Version/date : ', trim(adjustl(MPCDMD_version_date))
       write(*,*) 'MPCDMD> Git commit   : ', trim(MPCDMD_commit)
       write(*,*) 'MPCDMD> Built on     : ', trim(MPCDMD_machine)
       write(*,*) 'MPCDMD> Compiler     : ', trim(MPCDMD_compiler)
    end if


  end subroutine MPCDMD_info

end module MPCDMD

module sys
  
  type sys_t
     character(len=12) :: name
     integer :: N_max
     integer :: N_groups
     integer :: N_species
     integer, allocatable :: N(:)
     double precision, allocatable :: mass(:)
     double precision, allocatable :: oo_mass(:)
  end type sys_t

contains

  subroutine config_sys(sys_var, name, CF)
    use ParseText
    implicit none
    type(sys_t), intent(out) :: sys_var
    character(len=*), intent(in) :: name
    type(PTo), intent(in) :: CF

    character(len=24) :: temp_name
    character(len=9) :: format1
    character(len=15) :: format2
    integer :: len_name, i

    sys_var%name = name

    ! Reading the max number of particles and the number of species
    write(format1,'(a2,i02.2,a2,i02.2,a1)') '(a',len(name),',a',len('N_max'),')'
    write(temp_name,format1) name,'N_max'
    sys_var % N_max = PTread_i(CF,trim(temp_name))
    write(format1,'(a2,i02.2,a2,i02.2,a1)') '(a',len(name),',a',len('N_groups'),')'
    write(temp_name,format1) name,'N_groups'
    sys_var % N_groups = PTread_i(CF,trim(temp_name))
    write(format1,'(a2,i02.2,a2,i02.2,a1)') '(a',len(name),',a',len('N_species'),')'
    write(temp_name,format1) name,'N_species'
    sys_var % N_species = PTread_i(CF,trim(temp_name))


    allocate(sys_var%N(0:sys_var%N_species))
    sys_var%N = 0
    
    ! Reading the number of particles of each species
    do i=1,sys_var%N_species
       write(format2,'(a2,i02.2,a2,i02.2,a7)') '(a',len(name),',a',len('N'),',i02.2)'
       write(temp_name,format2) name,'N',i
       sys_var % N(i)  = PTread_i(CF,trim(temp_name))
    end do

    ! Checking that the max number of particles is not exceeded
    sys_var % N(0) = sum(sys_var%N(1:sys_var%N_species))
    if ( sys_var%N(0) > sys_var%N_max ) then
       write(*,*) 'Total number of particles for different species exceeds N_max'
       write(*,*) 'system ', sys_var%name
       stop
    end if

    allocate(sys_var%mass(sys_var%N_species))
    allocate(sys_var%oo_mass(sys_var%N_species))
    
    ! Reading the masses of the different species
    do i=1,sys_var%N_species
       write(format2,'(a2,i02.2,a2,i02.2,a7)') '(a',len(name),',a',len('mass'),',i02.2)'
       write(temp_name,format2) name,'mass',i
       sys_var % mass(i)  = PTread_d(CF,trim(temp_name))
       sys_var % oo_mass(i) = 1.d0 / sys_var%mass(i)
    end do

  end subroutine config_sys

end module sys

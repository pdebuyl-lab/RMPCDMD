subroutine test_config_group
  use group
  use ParseText
  implicit none

  type(group_t) :: g
  type(PTo) :: CF

  call PTparse(CF,'group_in',7)

  call config_group(g, 1, CF)

  call PTkill(CF)

  open(8,file='out')
  write(8,*) 'group01 = ', g % g_type
  write(8,*) 'group01N = ', g % N
  if (g%g_type .eq. DIMER_G) then
     write(8,*) 'group01length = ', g % dimer_length
  end if
  close(8)

end subroutine test_config_group

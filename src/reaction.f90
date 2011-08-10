
!> This module implements a "reac_t" datatype for solute-catalyzed chemical
!! reactions.
module reaction
  implicit none

  !> reac_t contains information about a given reaction.
  !! The properties are about the number of products (currently 1 or 2) and the
  !! immediate (or not) character of the reaction. The products are also given.
  type reac_t
     !> Toggles the reaction on
     logical(kind=1) :: on
     !> If set to .true. , a selected solvent will wait to be outside of the
     !! LJ cut-off to actually proceed in the reaction.
     logical(kind=1) :: at_exit
     !> If set to .true. , the reaction will produce two products, else only
     !! one product is created.
     logical(kind=1) :: two_products
     !>  If set to .true., the difference in u_int is taken into account when reacting.
     logical(kind=1) :: thermal
     !> Species of the first product.
     integer :: product1
     !> Species of the second product (only used when two_products==.true.).
     integer :: product2
     !> Reaction rate.
     double precision :: rate
  end type reac_t

contains

  !> Sets up a reaction variable.
  !!
  !! @param CF A ParseText configuration file.
  !! @param reac_var The resulting reaction variable.
  !! @param at_i The atom index.
  !! @param so_i The solvent index.
  subroutine config_reaction(CF, reac_var, at_i, so_i)
    use ParseText
    implicit none
    type(PTo), intent(in) :: CF
    type(reac_t), intent(out) :: reac_var
    integer, intent(in) :: at_i
    integer, intent(in) :: so_i
    
    character(len=5) :: at_so_s

    write(at_so_s, '(i2.02,a,i2.02)') at_i, '_', so_i
    
    reac_var % on = PTread_l(CF, 'reac'//at_so_s//'_on')

    if (reac_var % on ) then

       reac_var % at_exit = PTread_l(CF, 'reac'//at_so_s//'_at_exit')
       reac_var % two_products = PTread_l(CF, 'reac'//at_so_s//'_two_products')
       reac_var % thermal = PTread_l(CF, 'reac'//at_so_s//'_thermal')
       
       reac_var % product1 = PTread_i(CF, 'reac'//at_so_s//'_product1')
       if (reac_var % two_products) reac_var % product2 = PTread_i(CF, 'reac'//at_so_s//'_product2')
    
       reac_var % rate = PTread_d(CF, 'reac'//at_so_s//'_rate')

    end if

  end subroutine config_reaction

end module reaction

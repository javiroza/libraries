! Subrutina write --> Escriu dues línies en blanc en un arxiu
subroutine write(arxiu)
    ! arxiu --> número de l'arxiu
    implicit none
    integer arxiu

    write(arxiu,*) ""
    write(arxiu,*) ""

    return
end subroutine
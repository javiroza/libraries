! Aquest arxiu conté 4 subrutines amb mètodes diferents per calcular arrels
! de funcions donades.

program testo
    implicit none
    double precision a,b,eps,xarrel
    integer niter
    external f

    a=-1.2d0
    b=-1.1d0
    eps=0.0000001

    call Secant(a,b,eps,f,niter,xarrel)
    print*,xarrel
end program testo

! Subrutina RegulaFalsi --> retorna la posició d'un zero d'una funció donada
subroutine RegulaFalsi(A,B,eps,function,niter,xarrel)
    ! A,B --> Extrems de l'interval on hi ha l'arrel
    ! eps --> precisió desitjada
    ! function --> funció(subrutina) de la qual s'ha de trobar les arrels
    ! niter --> Nombre d'iteracions (comptador)
    ! xarrel --> coordenada x de l'arrel (a calcular)
    ! Nota: el criteri de convergència ha de comparar dues iteracions successives
    !       i no (B-A) amb eps.
    implicit none
    double precision A,B,eps,xarrel
    integer niter
    double precision C,fu,fa,fb,fc
    integer i

    niter=1
    call function(A,fa)
    call function(B,fb)

    if ((fa*fb).lt.0.d0) then
        11  C=(A*fb-B*fa)/(fb-fa)
        if (min(C-A,B-C).gt.eps) then
            niter=niter+1
            call function(C,fc)
            if ((fa*fc).lt.0.d0) then
                B=C
            else 
                A=C
            endif
            goto 11
        endif
    else
        print*,"No hi ha canvi de signe en l'interval seleccionat."
    endif

    xarrel = C
    return
end subroutine RegulaFalsi

! Subrutina Bisection --> retorna la posició d'un zero d'una funció donada
subroutine Bisection(A,B,eps,function,niter,xarrel) 
    ! A,B --> Extrems de l'interval on hi ha l'arrel
    ! eps --> precisió desitjada
    ! funciton --> funció de la qual s'ha de trobar les arrels
    ! niter --> nombre d'iteracions (comptador)
    ! xarrel --> coordenada x de l'arrel (a calcular)
    implicit none
    double precision A,B,eps,xarrel
    integer niter
    double precision C,fu,fa,fb,fc
    integer maxiter,i

    maxiter = nint(log((B-A)/eps)/log(2.d0)) ! Nombre màxim d'iteracions

    ! Algoritme de la bisecció
    do i=1,maxiter-1    
        call function(A,fa)
        call function(B,fb)
    ! Hi ha canvi de signe? Si és així, seguim
        if ((fa*fb).lt.0.d0) then
            C = (A+B)/2.d0
    ! L'interval és tant petit com es volia? Si és així, parem
            if ((B-A).lt.eps) then
                niter = i
                xarrel = C
                exit
            else
                call function(C,fc)
    ! Hem trobat la solució exacta? Si és així, parem         
                if (fc.eq.0.d0) then 
                    niter = i
                    xarrel = C
                    exit
                else
    ! En cas contrari, continuem iterant                    
                    if ((fa*fc).lt.0.d0) then
                        B = C
                    else if ((fc*fb).lt.0.d0) then
                        A = C
                    endif
                endif
            endif
        else
            print*,"No hi ha canvi de signe en l'interval donat."
            exit
        endif
    enddo

    xarrel = C
    niter = i
    return
end subroutine Bisection

! Subrutina NewtonRap --> retorna la posició d'un zero d'una funció donada
subroutine NewtonRap(x0,eps,function,niter,xarrel)
    ! x0 --> Punt d'inici
    ! eps --> Precisió desitjada
    ! function --> funció(subrutina) de la qual s'ha de trobar les arrels
        ! Nota: cal que la funció retorni la seva derivada a més del valor en el punt
    ! xarrel --> coordenada x de l'arrel (a calcular)
    implicit none
    double precision x0,x1,eps,xarrel,fu,dfu
    integer niter
    niter = 1

 12 call function(x0,fu,dfu)
    x1 = x0-(fu/dfu)
    if (abs(x1-x0).le.eps) then
        xarrel = x1
    else
        x0 = x1
        niter = niter+1
        goto 12
    endif

    return
end subroutine NewtonRap

! Subrutina Secant --> retorna la posició d'un zero d'una funció donada
subroutine Secant(x0,x1,eps,function,niter,xarrel)
    ! x0,x1 --> Punts d'inici
    ! eps --> Precisió desitjada
    ! function --> funció(subrutina) de la qual s'ha de trobar les arrels
    ! niter --> Nombre d'iteracions (comptador)
    ! xarrel --> coordenada x de l'arrel (a calcular)
    implicit none
    double precision x0,x1,eps,xarrel
    integer niter
    double precision x,fu,fx0,fx1
    niter = 1

    ! Inicialitzem el mètode amb dos punts inicials prou propers entre si
    call function(x0,fx0)
    call function(x1,fx1)

 13 x = x1-fx1*((x1-x0)/(fx1-fx0))

    if (abs(x-x1).le.eps) then
        xarrel = x
    else
        x0=x1
        fx0=fx1
        x1=x
        call function(x,fx1)
        niter=niter+1
        goto 13
    endif

    return
end subroutine Secant

subroutine f(x,fx)
    implicit none
    double precision x,fx

    fx=x**3.d0

    return
end subroutine f
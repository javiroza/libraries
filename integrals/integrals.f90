program testo
    implicit none
    double precision a,b,integral,f,M,pi,error,L,p
    integer k,N
    common/cts/pi,L
    external f,p
    call srand(20034276)
    pi=dacos(-1.d0)

    a=0.d0
    b=pi
    k=24
    N=2**k
    M=1.1d0

    call acceptrebuigintegral(a,b,N,M,f,integral,error)
    print*,integral
    print*,error
end program testo

! Subrutina trapezis --> Calcula una integral 1-D per trapezis
subroutine trapezis(x1,x2,k,funci,integral)
    ! x1,x2 --> Extrems de l'interval d'integració
    ! k --> N=2**k intervals
    ! funci --> funció a integrar
    ! integral --> valor de la integral (a calcular)
    ! Nota: s'assumeixen intervals de discretització iguals, h
    implicit none
    double precision x1,x2,funci,integral
    integer k
    double precision h
    integer N,i

    N = 2**k
    h = (x2-x1)/dble(N)
    integral = 0.d0

    do i=1,N-1 ! No es tenen en compte els extrems
        integral = integral + funci(x1+i*h)
    enddo

    ! Afegim els extrems i multipliquem per h
    integral = (integral + funci(x1)/2.d0 + funci(x2)/2.d0)*h

    return
end subroutine trapezis

! Subrutina simpson --> Calcula una integral 1-D per Simpson
subroutine simpson(x1,x2,k,funci,integral)
    ! x1,x2 --> Extrems de l'interval d'integració
    ! k --> N=2**k intervals
    ! funci --> funció a integrar
    ! integral --> valor de la integral (a calcular)
    ! Nota: s'assumeixen intervals de discretització iguals, h
    implicit none
    double precision x1,x2,funci,integral,h
    integer k
    double precision integral1,integral2
    integer N,i

    N = 2**k
    h = (x2-x1)/dble(N)
    integral1=0.d0
    integral2=0.d0

    do i=1,N-1,2
        integral1=integral1+funci(x1+i*h)
    enddo
    integral1=integral1*4.d0
    do i=2,N-2,2
        integral2=integral2+funci(x1+i*h)
    enddo

    integral2=integral2*2.d0
    integral=(integral1+integral2+funci(x1)+funci(x2))*(h/3.d0)

    return
end subroutine simpson

! Subrutina simpson38 --> Calcula una integral 1-D per Simpson
subroutine simpson38(x1,x2,k,funci,integral)
    ! x1,x2 --> Extrems de l'interval d'integració
    ! k --> N=2**k intervals
    ! funci --> funció a integrar
    ! integral --> valor de la integral (a calcular)
    ! Nota: s'assumeixen intervals de discretització iguals, h
    ! Nota2: aquest mètode no val la pena; milor fer servir simpson normal
    implicit none
    double precision x1,x2,funci,integral,h
    integer k
    integer N,i

    N = 2**k
    h = (x2-x1)/dble(N)
    integral=0.d0

    do i=0,N-3,3
        integral=integral+funci(x1+i*h)+3.d0*funci(x1+(i+1)*h)+ &
        3.d0*funci(x1+(i+2)*h)+funci(x1+(i+3)*h)
    enddo

    integral=integral*3.d0*h/8.d0

    return
end subroutine simpson38

! Subrutina boole --> Calcula una integral 1-D per Boole 
subroutine boole(x1,x2,k,funci,integral)
    ! x1,x2 --> Extrems de l'interval d'integració
    ! k --> N=2**k intervals
    ! funci --> funció a integrar
    ! integral --> valor de la integral (a calcular)
    ! Nota: s'assumeixen intervals de discretització iguals, h
    implicit none
    double precision x1,x2,funci,integral
    integer k
    double precision h
    integer N,i

    N=2**k
    h=(x2-x1)/dble(N)
    integral=0.d0

    do i=0,N-4,4
        integral=integral+7.d0*funci(x1+i*h)+32.d0*funci(x1+(i+1)*h)
        integral=integral+7.d0*funci(x1+(i+4)*h)+32.d0*funci(x1+(i+3)*h)
        integral=integral+12.d0*funci(x1+(i+2)*h)
    enddo

    integral=integral*(2.d0*h/45.d0)

    return
end subroutine boole

! Subrutina gl2 --> Calcula una integral 1-D amb Gauss-Legendre 2
subroutine gl2(x1,x2,funci,integral)
    ! x1,x2 --> Extrems de l'interval d'integració
    ! funci --> funció a integrar
    ! integral --> valor de la integral (a calcular)
    ! Nota: per a una major eficàcia, la funció a integrar s'hauria de poder
    !       aproximar bé per un polinomi de grau n<=3
    implicit none
    double precision x1,x2,funci,integral
    double precision xi ! Variable auxiliar

    !u=-1.d0+2.d0*((x-x1)/(x2-x1)) --> Canvi de variable
    xi=1.d0/dsqrt(3.d0)

    integral=funci(x1+((-xi+1.d0)/2.d0)*(x2-x1))+funci(x1+((xi+1.d0)/2.d0)*(x2-x1))
    integral=integral*0.5d0*(x2-x1) ! Jacobià del canvi

    return
end subroutine

! Subrutina acceptrebuigintegral --> Calcula una integral 1-D amb accept-rebuig
subroutine acceptrebuigintegral(a,b,N,M,funci,integral,error)
    ! a,b --> Extrems de l'interval d'integració
    ! N --> nombre de números aleatoris que es vol generar
    ! M --> cota superior
    ! funci --> funció a integrar
    ! integral --> Valor de la integral (a calcular)
    ! error --> Estimació de l'error comès (output)
    ! Nota: aquest mètode només serveix en dominis acotats on la funció a
    !       integrar sigui positiva
    implicit none
    double precision a,b,M,funci,integral,error
    integer N
    double precision x,p
    integer counter_dins,counter_fora
    counter_dins=0
    counter_fora=0

    do while ((counter_dins+counter_fora).lt.N)
        x=(b-a)*rand()+a 
        p=M*rand()
        if (funci(x).ge.p) then 
            counter_dins=counter_dins+1
        else 
            counter_fora=counter_fora+1
        endif
    enddo

    integral=M*(b-a)*(counter_dins/dble(counter_dins+counter_fora))
    error=(M*(b-a)/dsqrt(dble(N)))*dsqrt((counter_dins/dble(N))*(1.d0-(counter_dins/dble(N))))

    return
end subroutine acceptrebuigintegral

! Subrutina montecarlocru --> Calcula una integral 1-D amb Montecarlo "cru"
subroutine montecarlocru(a,b,N,funci,integral,error)
    ! a,b --> Extrems de l'interval d'integració (input)
    ! N --> Nombre de números aleatoris que es volen generar (input)
    ! funci --> Funció a integrar (input)
    ! integral --> Valor de la integral (a calcular) (output)
    ! error --> Error comès en el càlcul de integral (output)
    implicit none
    double precision a,b,funci,integral,error
    integer N
    double precision x
    integer i
    integral=0.d0
    error=0.d0

    do i=1,N
        x=rand()
        integral=integral+funci((b-a)*x+a)
        error=error+(funci((b-a)*x+a))**2.d0
    enddo

    integral=integral*(b-a)/dble(N)
    error=(dble(N))**(-0.5d0)*((error*((b-a)**2.d0)/dble(N))-integral**2.d0)**0.5d0

    return
end subroutine montecarlocru

double precision function f(x)
    implicit none
    double precision x
    double precision pi,L
    common/cts/pi,L

    f=(dsin(x))**2.d0

    return
end function f

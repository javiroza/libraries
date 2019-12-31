program testo
    implicit none
    double precision f,polinterpol,x_p,interpol
    double precision, allocatable :: x(:),funci(:)
    integer N,i
    external f

    N=100
    allocate(x(N))
    allocate(funci(N))

    do i=0,N
        x(i)=i*5.d0/100.d0
        funci(i)=f(x(i))
    enddo

    print*,funci(100)
    x_p=1.256d0
    print*,interpol(N,x,funci,x_p)
end program testo
        

! Subrutina interpol --> Calcula el valor d'una funció interpolada amb Lagrange
double precision function interpol(N,x,funci,x_p) result(polinterpol)
    ! N --> dimensió dels vectors x i funci
    ! x --> vector amb els valors de la variable independent (input)
    ! funci --> vector amb els imatges del vector x (input)
    ! x_p --> punt on es vol avaluar el polinomi interpolador
    ! polinterpol --> polinomi interpolador avaluat a x_p (output)
    implicit none
    double precision x(N),funci(N),x_p
    integer N
    double precision sumand
    integer i,k

    polinterpol=0.d0
    sumand=1.d0
    do i=0,N
        do k=0,N
            if (k.ne.i) then
                sumand=sumand*((x_p-x(k))/(x(i)-x(k)))
            endif
        enddo
        sumand=sumand*funci(i)
        polinterpol=polinterpol+sumand
        sumand=1.d0
    enddo

    return
end function interpol

! Funció f --> funció de prova
double precision function f(x)
    implicit none
    double precision x

    f=x**2.d0

    return
end function f

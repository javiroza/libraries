program test1
    implicit none
    double precision, allocatable :: y0(:),y1(:)
    integer nvar,i
    double precision g,l,h,f
    common/cts/g,l
    external f
    g=9.81d0
    l=1.d0
    nvar=2
    allocate(y0(nvar))
    allocate(y1(nvar))
    h=0.001d0

    y0=(/0.d0,1.d0/)

    call euler(0.d0,h,nvar,y0,y1)
    print*,y1

end program test1

! Subrutina euler --> Calcula un pas del mètode d'euler N-D
subroutine euler(x,dx,nequs,yyin,yyout)
    ! x --> Variable independent del problema
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt anterior
    ! yyout --> Vector amb les dades del punt següent
    ! Nota: aquesta subrutina necessita la subrutina derivades per funcionar
    implicit none
    double precision x,dx,yyin(nequs),yyout(nequs)
    integer nequs
    double precision dyout(nequs)

    call derivades(nequs,x,yyin,dyout)
    yyout=yyin+dx*dyout

    return
end subroutine euler

! Subrutina eulerm --> Calcula un pas del mètode d'euler millorat N-D
subroutine eulerm(x,dx,nequs,yyin,yyout)
    ! x --> Variable independent del problema
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt 2 vegades anterior
    ! yyout --> Vector amb les dades del punt següent
    ! Nota1: aquesta subrutina necessita la subrutina derivades per funcionar
    ! Nota2: el primer pas s'ha de donar amb Euler simple
    implicit none
    double precision x,dx,yyin(nequs),yyout(nequs)
    integer nequs
    double precision dyout(nequs)

    call derivades(nequs,x,yyin,dyout)
    yyout=yyin+2.d0*dx*dyout

    return
end subroutine eulerm

! Subrutina predcorr --> Calcula un pas del mètode de predictor-corrector
subroutine predcorr(x,dx,nequs,yyin,yyout)
    ! x --> Variable independent del problema
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt 2 vegades anterior
    ! yyout --> Vector amb les dades del punt següent
    ! Nota1: aquesta subrutina necessita la subrutina derivades per funcionar
    ! Nota2: el primer pas s'ha de donar amb Euler simple
    implicit none
    double precision x,dx,yyin(nequs),yyout(nequs)
    integer nequs
    double precision dyout(nequs),y_pred(nequs)

    ! Càlcul de y_pred
    call derivades(nequs,x,yyin,dyout)
    y_pred=yyin+dx*dyout
    ! Càlcul de la part de yyout que conté f(x,y_pred)
    call derivades(nequs,x,y_pred,dyout)
    yyout=yyin+0.5*dx*dyout
    ! Càlcul de la part de yyout que conté f(x,yyin)
    call derivades(nequs,x,yyin,dyout)
    yyout=yyout+0.5d0*dx*dyout

    return
end subroutine predcorr

! Subrutina RLSTN3 --> Calcula un pas del mètode de Ralston de tercer ordre per un sistema
! d'n equacions de primer ordre acoblades
subroutine RLSTN3(x,dx,nequs,yyin,yyout)
    ! x --> Variable independent del problema
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt anterior
    ! yyout --> Vector amb les dades del punt següent
    implicit none
    integer nequs,i
    double precision x,dx,yyin(nequs),yyout(nequs)
    double precision k1(nequs),k2(nequs),k3(nequs)
    double precision dyout(nequs) ! Vector mut necessari per cridar la subrutina derivades

    ! Càlcul dels vectors k1,k2,k3 
    call derivades(nequs,x,yyin,dyout)
    k1=dyout
    call derivades(nequs,x+dx/2.d0,yyin+dx/2.d0*k1,dyout)
    k2=dyout
    call derivades(nequs,x+0.75d0*dx,yyin+(3.d0*dx/4.d0)*k2,dyout)
    k3=dyout

    ! Càlcul del vector yyout 
    yyout=yyin+dx/9.d0*(2.d0*k1+3.d0*k2+4.d0*k3)

    return 
end subroutine RLSTN3

! Subrutina RK4 --> Calcula un pas del mètode de Runge-Kutta d'ordre 4 
subroutine RK4(x,dx,nequs,yyin,yyout)
    ! x --> Variable independent del problema
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt anterior
    ! yyout --> Vector amb les dades del punt següent
    implicit none
    integer nequs,i
    double precision x,dx,yyin(nequs),yyout(nequs)
    double precision k1(nequs),k2(nequs),k3(nequs),k4(nequs)
    double precision dyout(nequs) ! Vector mut necessari per cridar la subrutina derivades

    ! Càlcul dels vectors k1,k2,k3 
    call derivades(nequs,x,yyin,dyout)
    k1=dyout
    call derivades(nequs,x+dx/2.d0,yyin+dx/2.d0*k1,dyout)
    k2=dyout
    call derivades(nequs,x+dx/2.d0,yyin+dx/2.d0*k2,dyout)
    k3=dyout
    call derivades(nequs,x+dx,yyin+dx*k3,dyout)
    k4=dyout

    ! Càlcul del vector yyout 
    yyout=yyin+dx/6.d0*(k1+2.d0*k2+2.d0*k3+k4)

    return 
end subroutine RK4

! Subrutina cash_karp --> Calcula un pas del mètode de Cash-Karp (RK6)
subroutine cash_karp(x,dx,nequs,yyin,yyout)
    ! x --> Variable independent del problema
    ! dx --> Pas 
    ! nequs --> Nombre d'equacions 
    ! yyin --> Vector amb les dades del punt anterior
    ! yyout --> Vector amb les dades del punt següent
    implicit none
    integer nequs,i
    double precision x,dx,yyin(nequs),yyout(nequs)
    double precision k1(nequs),k2(nequs),k3(nequs),k4(nequs),k5(nequs),k6(nequs)
    double precision dyout(nequs) ! Vector mut necessari per cridar la subrutina derivades

    ! Càlcul dels vectors k1,k2,k3 
    call derivades(nequs,x,yyin,dyout)
    k1=dx*dyout
    call derivades(nequs,x+dx*0.2d0,yyin+0.2d0*k1,dyout)
    k2=dx*dyout
    call derivades(nequs,x+dx*0.3d0,yyin+(3.d0/40.d0)*k1+(9.d0/40.d0)*k2,dyout)
    k3=dx*dyout
    call derivades(nequs,x+dx*(3.d0/5.d0),yyin+0.3d0*k1-0.9d0*k2+(6.d0/5.d0)*k3,dyout)
    k4=dx*dyout
    call derivades(nequs,x+dx,yyin-(11.d0/54.d0)*k1+2.5d0*k2-(70.d0/27.d0)*k3+(35.d0/27.d0)*k4,dyout)
    k5=dx*dyout
    call derivades(nequs,x+dx*(7.d0/8.d0),yyin+k1*(1631.d0/55296.d0)+k2*(175.d0/512.d0)+k3*(575.d0/13824.d0) & 
        +k4*(44275.d0/110592.d0)+k5*(253.d0/4096.d0),dyout)
    k6=dx*dyout

    ! Càlcul del vector yyout 
    yyout=yyin+(37.d0/378.d0)*k1+(250.d0/621.d0)*k3+(125.d0/594.d0)*k4+(512.d0/1771.d0)*k6

    return 
end subroutine cash_karp

! Subrutina solver --> Calcula una iteració per resoldre l'equació de Poisson
subroutine solver(T_old,T_new,h,funci,icontrol)
    ! T_old --> Matriu amb els valors a recalcular
    ! T_new --> Matriu amb els valors nous calculats
    ! h --> Interval de particionat de la malla
    ! funci --> "fonts" de l'equació de Poisson
    ! icontrol --> Variable de control del mètode de resolució
        ! 0 --> Jacobi
        ! 1 --> Gauss-Seidel
        ! 2 --> Sobrerelaxació
    ! Nota: si en comptes de Poisson es vol resoldre Laplace, només s'ha de treure el terme h² de cada mètode
    implicit none
    integer i,j,Nx,Ny,icontrol
    double precision T_old(Nx,Ny),T_new(Nx,Ny),h,funci,Lx,Ly,w
    common/cts/Lx,Ly,Nx,Ny,w

    if (icontrol.eq.0) then ! Jacobi
        do i=2,Nx-1
            do j=2,Ny-1
                T_new(i,j)=0.25d0*(T_old(i,j+1)+T_old(i,j-1)+T_old(i+1,j)+T_old(i-1,j)+(h**2.d0)*funci(i*h,j*h))
            enddo
        enddo
    else if (icontrol.eq.1) then ! Gauss-Seidel
        do i=2,Nx-1
            do j=2,Ny-1
                T_new(i,j)=0.25d0*(T_old(i,j+1)+T_old(i,j-1)+T_old(i+1,j)+T_old(i-1,j)+(h**2.d0)*funci(i*h,j*h))
                T_old(i,j)=T_new(i,j)
            enddo
        enddo
    else if (icontrol.eq.2) then ! Sobrerelaxació
        do i=2,Nx-1
            do j=2,Ny-1
                T_new(i,j)=T_old(i,j)+0.25d0*w*(T_old(i,j+1)+T_old(i,j-1)+T_old(i+1,j)+T_old(i-1,j)-4.d0*T_old(i,j)) 
                T_new(i,j)=T_new(i,j)+0.25d0*w*(h**2.d0)*funci(i*h,j*h) ! Aquesta línea és la continuació de la de dalt
                T_old(i,j)=T_new(i,j)
            enddo
        enddo
    endif

    return 
end subroutine solver

! Subrutina derivades --> Calcula la derivada de la funció (vectorial) a trobar
subroutine derivades(nequ,x,yin,dyout)
    ! nequ --> nombre d'equacions
    ! x --> Variable independent del problema
    ! yin --> Vector amb les dades del punt on es calcula la derivada
    ! dyout --> Vector amb les derivades
    implicit none
    integer nequ
    double precision x,yin(nequ),dyout(nequ)
    double precision g,l
    common/cts/g,l

    dyout(1)=yin(2)
    dyout(2)=-(g/l)*dsin(yin(1))

    return
end subroutine derivades



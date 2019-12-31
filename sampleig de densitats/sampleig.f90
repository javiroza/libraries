program testo
    implicit none
    double precision, allocatable :: xnums(:)
    integer ndat

    ndat=1000
    allocate(xnums(ndat))

    call boxmuller(ndat,xnums)
    print*,xnums(345)
end program testo

! Subrutina acceptrebuig --> genera un vector amb ndat números distribuïts segons funci
subroutine acceptrebuig(ndat,xnums,a,b,M,funci)
    ! ndat --> nombre de números aleatoris que es vol generar (input)
    ! xnums --> vector que contindrà els ndat números aleatoris (output)
    ! a,b --> extrems de l'interval on està definida la nova variable aleatòria (input)
    ! M --> cota superior (input)
    ! funci --> densitat de probabilitat segons la qual està distribuida la nova var. aleat. (output)
    ! Nota --> la subrutina també calcula la variància i la desviació estàndard de la nova
    ! distribució. Es pot esborrar tranquilament aquesta part sense afectar la generació de nombres. 
    implicit none
    double precision xnums(ndat),a,b,M,funci,x1,x2,p,x
    double precision valmitj,var,desvest
    integer ndat,iseed,counter,i
    counter=0 ! Variable que porta el compte dels nombres aleatoris generats

    ! Generació dels nombres aleatoris en distribució uniforme [0,1]
 1  x1=rand()
    x2=rand()
    ! Canvi de variable x1-->x, x2-->p 
    x=(b-a)*x1+a 
    p=M*x2 
    ! Ara, x és U(a,b) i p és U(0,M). Iniciem la comprovació
    if (funci(x).ge.p) then 
        xnums(counter)=x
        counter=counter+1
        if (counter.le.ndat) then 
            goto 1
        endif
    else 
        goto 1
    endif

    ! Càlcul del valor mitjà
    valmitj=0.d0
    do i=0,ndat
        valmitj=valmitj+xnums(i)
    enddo
    valmitj=valmitj/dble(ndat)

    ! Càlcul de la variància
    var=0.d0
    do i=0,ndat
        var=var+(xnums(i)-valmitj)**2.d0
    enddo
    var=var/dble(ndat)

    ! "Càlcul" de la desviació estàndard
    desvest=var**0.5d0

    print*,""
    print*,"Resultats exercici 1:"
    print*,valmitj,var,desvest
    return
end subroutine acceptrebuig

! Subrutina boxmuller --> genera un vector amb ndat números amb distribució normal centrada a 0
subroutine boxmuller(ndat,xnormal)
    ! ndat --> nombre de números aleatoris que es vol generar (input)
    ! xnormal --> vector que contindrà els ndat números aleatoris (output)
    implicit none
    double precision xnormal(ndat)
    integer ndat
    double precision pi,r,phi
    integer i
    pi=dacos(-1.d0)

    do i=1,ndat-1,2
        r=dsqrt(-2.d0*log(rand()))
        phi=2.d0*pi*rand()
        xnormal(i)=r*dcos(phi)
        xnormal(i+1)=r*dsin(phi)
    enddo

    return
end subroutine boxmuller

! Subrutina HISTOGRAMA --> agafada del campus virtual, genera un histograma
SUBROUTINE HISTOGRAMA(NDAT,XDATA,XA,XB,NBOX,XHIS,VHIS,ERRHIS,BOXSIZE,IERR)
    ! ndat --> longitud del vector XDATA (input)
    ! xdata --> vector amb els nombres aleatoris (input)
    ! xa,xb --> extrems de l'interval en què està definida la variable aleatòria (input)
    ! nbox --> nombre de caixes de l'histograma (input)
    ! xhis,vhis,errhis --> vectors generats automàticament per la subrutina, són el
    !                      propi histograma (output)
    ! boxsize --> mida de cada caixa (automàtic)
    ! ierr --> paràmetre de control de la validesa dels inputs
    ! Exemple d'utilització de la subrutina: pre-pràctica 5, Exercici 1 + histotest.gnu
    IMPLICIT NONE
    ! INPUT/OUTPUT VARIABLES
    INTEGER NDAT,NBOX
    DOUBLE PRECISION XDATA(NDAT),XA,XB
    DOUBLE PRECISION XHIS(NBOX),VHIS(NBOX),ERRHIS(NBOX)
    INTEGER IERR

    INTEGER I,IBOX,ICOUNT
    DOUBLE PRECISION BOXSIZE

    IF (XA.GE.XB) THEN 
        IERR=1
        RETURN
    ENDIF
    ! BOX SIZE
    BOXSIZE=(XB-XA)/NBOX

    ! COUNTS NUMBER OF POINTS WITHIN THE INTERVAL XA,XB
    ICOUNT=0

    ! SETS ALL TO ZERO
    DO I=1,NBOX
        VHIS(I)=0
        ERRHIS(I)=0
    ENDDO

    ! WE RUN THROUGH THE DATASET
    DO I=1,NDAT
    ! CHECKS IF DATA LIES WITHIN XA,XB
    IF (XDATA(I).GE.XA.AND.XDATA(I).LE.XB) THEN 
        IBOX=INT((XDATA(I)-XA)/BOXSIZE)+1
    ! PUTS XB INTO THE LAST BOX, IF NEEDED
        IF (IBOX.EQ.NBOX+1) IBOX=NBOX 

            VHIS(IBOX)=VHIS(IBOX)+1
            ICOUNT=ICOUNT+1
    ENDIF
    ENDDO

    IF (ICOUNT.EQ.0) THEN 
       IERR=2
       RETURN
    ENDIF

    IERR=0
    PRINT*,"ACCEPTED:",ICOUNT," OUT OF:",NDAT

    DO I=1,NBOX
    ! CENTRAL VALUE OF THE BAR
        XHIS(I)=XA+BOXSIZE/2.D0+(I-1)*BOXSIZE
    !  ERROBAR, STANDARD DEVIATION OF CORRESPONDING BINOMIAL
        ERRHIS(I)=SQRT(VHIS(I)/ICOUNT*(1.D0-VHIS(I)/ICOUNT))/BOXSIZE / SQRT(DBLE(ICOUNT))
    ! NORMALIZED VALUE OF THE BAR
    VHIS(I)=VHIS(I)/ICOUNT/BOXSIZE
    ENDDO
END SUBROUTINE HISTOGRAMA

! Subrutina write --> Escriu dues línies en blanc en un arxiu
subroutine write(arxiu)
    ! arxiu --> número de l'arxiu
    implicit none
    integer arxiu

    write(arxiu,*) ""
    write(arxiu,*) ""

    return
end subroutine

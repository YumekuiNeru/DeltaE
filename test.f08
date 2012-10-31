INTEGER FUNCTION prime(n) result(res)
    IMPLICIT none
    INTEGER, intent(in) :: n
    INTEGER, DIMENSION(n):: primes
    INTEGER :: i,j,exists=1
    INTEGER(kind=16):: x,y,xy
    res = 1
    primes=0
    DO y=1,int(float(n)**(1.0/4.0)),1
why:    DO x=1,int(sqrt(float(n))),1
            xy = x**2 + y**4
            IF (xy > n) THEN
                CYCLE why
            ELSEIF (xy == 2) THEN
13              exists=1
                DO i=1,res,1 !check if xy is already recorded
                    IF (primes(i) == 0) THEN
                        exists = 0
                        !print *,x,y,xy
                        EXIT
                    ELSEIF (primes(i) == xy) THEN
                        exists=1
                        EXIT
                    ENDIF
                ENDDO
                IF (exists /= 1) THEN
                    primes(res) = xy
                    res = res + 1
                ENDIF
                CYCLE why
            ENDIF
            DO j=2,int(sqrt(float(xy)))+1,1 !check if xy is prime
                IF (mod(xy,j) == 0) THEN
                    CYCLE why
                ENDIF
            ENDDO
            GOTO 13 !do the array checking again (goto avoids repeated code)
        ENDDO why
    ENDDO
    res = res-1 !(because the current position contains 0)
ENDFUNCTION prime

PROGRAM nprime
    IMPLICIT none
    INTEGER :: prime,res,nn,k,n
    !provide amount of test cases
    READ*,nn
    DO k=1,nn,1
        !PRINT "(A10)","Provide n:"
        READ *,n
        res = prime(n)
        PRINT "(i0)",res
        !PRINT "(A10i0A23i0A20)","There are ",res," possible primes below ",n," on the form x^2+y^4"
    ENDDO
ENDPROGRAM nprime
PROGRAM DeltaE
    IMPLICIT none
    REAL(KIND=4), DIMENSION(3) :: sRGB1, sRGB2, XYZ1, XYZ2, Lab1, Lab2
    INTEGER, DIMENSION(3) :: RGB1,RGB2
    INTEGER :: testing = 0
    REAL(KIND=4) :: diff
    !RGB1 = (/136, 220, 73/) !reference colour
    !RGB2 = (/88, 99, 111/) !colour to compare to
    
    !sRGB1 = RGB2sRGB(RGB1)
    !sRGB2 = RGB2sRGB(RGB2)
    
    !XYZ1 = sRGB2XYZ(sRGB1)
    !XYZ2 = sRGB2XYZ(sRGB2)
    
    !Lab1 = XYZ2LAB(XYZ1)
    !Lab2 = XYZ2LAB(XYZ2)
    call test()
    !diff = dE_single(Lab1, Lab2)
    !PRINT*,""
    !PRINT*,"Colour 1:",RGB1
    !PRINT*,"Colour 1:",Lab1
    
    !PRINT*,"Colour 2:",RGB2
    !PRINT*,"Colour 2:",Lab2
    
    !PRINT*,"Colourdifference (lower=less diff): ",diff
    
CONTAINS
    SUBROUTINE test()
        IMPLICIT none
        REAL(KIND=4), DIMENSION(34,39) :: cases
        REAL(KIND=4), DIMENSION(3) :: temp1, temp2
        REAL(KIND=4) :: deltares, diff
        INTEGER :: i
        
        ! 1   2   3   4   5   6   7   8  ...   39 
        !L1, a1, b1, L2, a2, b2, C1, C2, ..., dE2000    
        OPEN (13,file="testcases.csv")
        READ (13,*)
        DO i=1,34,1
            READ (13,*) cases(i,:)
        ENDDO
        PRINT "(A4A12A10A16)","test","Expected","Result","Difference"
        DO i=1,34,1
            temp1(1:3) = (/cases(i,1:3)/)
            temp2(1:3) = (/cases(i,4:6)/)
            deltares = dE_single(temp1, temp2)
            diff = abs(cases(i,39) - deltares)
            PRINT "(i3,f12.4,f12.4,f12.4)",i,cases(i,39), deltares, diff
        ENDDO
    ENDSUBROUTINE test
!-------------------------------------------------------------------------------
    FUNCTION RGB2sRGB(RGB) result(sRGB)
        IMPLICIT none
        INTEGER, DIMENSION(3), INTENT(in) :: RGB
        REAL(KIND=4), DIMENSION(3) :: sRGB
        sRGB = RGB/255.0
    ENDFUNCTION RGB2sRGB
!-------------------------------------------------------------------------------
    FUNCTION sRGB2XYZ(sRGB) RESULT(XYZ)
        IMPLICIT none
        REAL(KIND=4), DIMENSION(3), intent(in) :: sRGB
        REAL(KIND=4), DIMENSION(3) :: XYZ, CLIN
        REAL(KIND=4), DIMENSION(3,3) :: component
        REAL(KIND=4) :: a=0.055, C
        INTEGER :: i

        !sRGB D65
        DATA component /0.4124564,0.2126729,0.0193339,&
                        0.3575761,0.7151522,0.1191920,& !<-- do not read this as a matrix
                        0.1804375,0.0721750,0.9503041/
        !component(1,:) = (/0.4124,0.3576,0.1805/)
        !component(2,:) = (/0.2126,0.7152,0.0722/) ! <-- Matrix will look like this and stuffs =w=
        !component(3,:) = (/0.0193,0.1192,0.9505/)
        
        CLIN = sRGB
        DO i=1,3,1
            C = CLIN(i)
            IF (C <= 0.04045) THEN
                CLIN(i) = C/12.92
            ELSE
                CLIN(i) = ((C+a)/(1.+a))**2.4
            ENDIF
        ENDDO
        XYZ = matmul(component,CLIN)
    ENDFUNCTION sRGB2XYZ        
!-------------------------------------------------------------------------------
        REAL(KIND=4) FUNCTION ft(t) result(res)
            IMPLICIT none
            REAL(KIND=4), intent(in) :: t
            REAL(KIND=4) :: one_third=1.0/3.0 
            
            IF (t > (6./29.)**3.0) THEN
                res = t**one_third
            ELSE
                res = one_third*((29./6.)**(2.))*t + 4./29.
            ENDIF
        ENDFUNCTION ft
        
    FUNCTION XYZ2Lab(XYZ) RESULT(Lab)
        IMPLICIT none
        REAL(KIND=4), DIMENSION(3), INTENT(in) :: XYZ
        REAL(KIND=4), DIMENSION(3) :: Lab
        !#CIE XYZ tristimulus values of the reference white point
        REAL(KIND=4) :: Xn=0.95047,Yn=1.00000,Zn=1.08883
        
        !L = 116.0*ft(Y/Yn) - 16
        !a = 500.0*(ft(X/Xn) - ft(Y/Yn))
        !b = 200.0*(ft(Y/Yn) - ft(Z/Zn))
        !XYZ   = [1X, 2Y, 3Z]
        Lab(1) = 116.0*ft(XYZ(2)/Yn) - 16
        Lab(2) = 500.0*(ft(XYZ(1)/Xn) - ft(XYZ(2)/Yn))
        Lab(3) = 200.0*(ft(XYZ(2)/Yn) - ft(XYZ(3)/Zn))
    ENDFUNCTION XYZ2Lab
!-------------------------------------------------------------------------------    
    !Given two colours, ref and poll, returns their deltaE (difference)
    REAL(KIND=4) FUNCTION dE_single(rL, pL) result(dE00)
        IMPLICIT none                  
        REAL(KIND=4), DIMENSION(3), INTENT(in) :: rL, pL !ref_Lab, poll_Lab
        REAL(KIND=4), PARAMETER :: PI = 3.1415926535898
        REAL(KIND=4), PARAMETER :: deg2rad = PI/180.0
        REAL(KIND=4), PARAMETER :: rad2deg = 180.0/PI
        REAL(KIND=4), PARAMETER :: PRECISE = 0.0001
        REAL(KIND=4) :: c1s_ab, c2s_ab, Csbar_ab, G, a1p, a2p, C1p, C2p, h1p, h2p, h2p_h1p, &
                        dLp, dCp, dhp, dHHp, Lbarp, Cbarp, h1ph2p, h1p_h2p, hbarp, T, d0, &
                        R_C, S_L, S_C, S_H, R_T, K_L, K_C, K_H
        
        !i2 = poll, i1 = ref
        !Lab = [1L, 2a, 3b]
        C1s_ab = sqrt((rL(2)**2) + (rL(3)**2))
        C2s_ab = sqrt((pL(2)**2) + (pL(3)**2))
        
        Csbar_ab = (C1s_ab + C2s_ab)/2.0
        
        G = 0.5*(1.0-sqrt((Csbar_ab**7)/(Csbar_ab**7+25.0**7)))
        
        a1p = (1+G)*rL(2)
        a2p = (1+G)*pL(2)
        
        C1p = sqrt(a1p**2 + rL(3)**2)
        C2p = sqrt(a2p**2 + pL(3)**2)
        
        IF (rL(3) == 0 .AND. a1p == 0) THEN
            h1p = 0
        ELSE
            h1p = atan2(rL(3),a1p)
            IF (h1p < 0) THEN
                h1p = h1p + 2*pi
            ENDIF
            h1p = rad2deg*h1p
        ENDIF
        IF (pL(3) == 0 .AND. a2p == 0) THEN
            h2p = 0
        ELSE 
            h2p = atan2(pL(3),a2p)
            IF (h2p < 0) THEN
                h2p = h2p + 2*pi
            ENDIF
            h2p = rad2deg*h2p
        ENDIF
        
        dLp = pL(1) - rL(1)
        dCp = C2p - C1p
        
        h2p_h1p = h2p - h1p

        IF (C1p*C2p == 0) THEN
            dhp = 0
        ELSEIF (abs(h2p_h1p) - 180 <= PRECISE) THEN
            dhp = h2p_h1p
        ELSEIF (h2p_h1p - 180 > PRECISE) THEN
            dhp = h2p_h1p - 360
        ELSEIF (h2p_h1p + 180 < PRECISE) THEN
            dhp = h2p_h1p + 360
        ENDIF
        
        dHHp = 2.0*sqrt(C1p*C2p)*sin(deg2rad*dhp/2.0)
        
        Lbarp = (rL(1) + pL(1))/2.0
        Cbarp = (C1p + C2p)/2.0
        
        h1ph2p = h1p + h2p
        h1p_h2p = h1p - h2p

        IF (C1p == 0 .OR. C2p == 0) THEN
            hbarp = h1ph2p
        ELSEIF (abs(h1p_h2p) - 180 <= PRECISE) THEN
            hbarp = (h1ph2p)/2.0
        ELSEIF (abs(h1p_h2p) - 180 > PRECISE .AND. (h1ph2p - 360 < PRECISE)) THEN
            hbarp = (h1ph2p + 360.0)/2.0
        ELSEIF (abs(h1p_h2p) - 180 > PRECISE .AND. (h1ph2p - 360 >= PRECISE)) THEN
            hbarp = (h1ph2p - 360.0)/2.0
        ENDIF
        
        T = 1 - 0.17*cos(deg2rad*(hbarp-30.0)) + 0.24*cos(deg2rad*(2*hbarp)) + &
            0.32*cos(deg2rad*(3*hbarp + 6.0)) - 0.20*cos(deg2rad*(4*hbarp-63.0))
        
        d0 = 30.0*exp(-((hbarp-275.0)/25.0)**2.0)
        
        R_C = 2*sqrt((Cbarp**7.0)/((Cbarp**7)+25.0**7))
        R_T = -sin(deg2rad*(2.0*d0))*R_C
        
        S_L = 1.0 + ((0.015*((Lbarp-50)**2))/(sqrt(20.0+((Lbarp-50.0)**2))))
        S_C = 1.0 + 0.045*Cbarp
        S_H = 1.0 + 0.015*Cbarp*T
        
        K_L = 1
        K_C = 1
        K_H = 1
        
        dE00 = sqrt((dLp/(K_L*S_L))**2 + (dCp/(K_C*S_C))**2 + (dHHp/(K_H*S_H))**2 + R_T*(dCp/(K_C*S_C))*(dHHp/(K_H*S_H)))
    ENDFUNCTION dE_single
    
!-------------------------------------------------------------------------------    
    !Given one colour (ref) and an array of colours (polls), return the 
    !REAL(KIND=4) FUNCTION dE_multi(rL, pL) result(dE00)
   
    
!-------------------------------------------------------------------------------    
ENDPROGRAM DeltaE
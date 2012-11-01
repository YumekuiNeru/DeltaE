PROGRAM DeltaE
    IMPLICIT none
    REAL(kind=8), DIMENSION(3) :: sRGB1, sRGB2, XYZ1, XYZ2
    REAL(kind=8), DIMENSION(5) :: LABCH1, LABCH2
    INTEGER, DIMENSION(3) :: RGB1,RGB2
    REAL(kind=8) :: diff
    
    RGB1 = (/136, 220, 73/) !reference colour
    RGB2 = (/88, 99, 111/) !colour to compare to
    
    sRGB1 = RGB2sRGB(RGB1)
    sRGB2 = RGB2sRGB(RGB2)
    
    XYZ1 = sRGB2XYZ(sRGB1)
    XYZ2 = sRGB2XYZ(sRGB2)
    
    LABCH1 = XYZ2LABCH(XYZ1)
    LABCH2 = XYZ2LABCH(XYZ2)
    
    diff = dE_single(LABCH1, LABCH2)
    PRINT*,""
    PRINT*,"Colour 1:",RGB1
    PRINT*,"Colour 1:",LABCH1
    
    PRINT*,"Colour 2:",RGB2
    PRINT*,"Colour 2:",LABCH2
    
    PRINT*,"Colourdifference (lower=less diff): ",diff
    
CONTAINS
!-------------------------------------------------------------------------------
        FUNCTION RGB2sRGB(RGB) result(sRGB)
        IMPLICIT none
        INTEGER, DIMENSION(3), INTENT(in) :: RGB
        REAL(KIND=8), DIMENSION(3) :: sRGB
        sRGB = RGB/255.0
    ENDFUNCTION RGB2sRGB
!-------------------------------------------------------------------------------
    FUNCTION sRGB2XYZ(sRGB) RESULT(XYZ)
        IMPLICIT none
        REAL(kind=8), DIMENSION(3), intent(in) :: sRGB
        REAL(kind=8), DIMENSION(3) :: XYZ, CLIN
        REAL(kind=8), DIMENSION(3,3) :: component
        REAL(kind=8) :: a=0.055, C
        INTEGER :: i
        
        
        !sRGB D65
        data component /0.4124564,0.2126729,0.0193339,&
                        0.3575761,0.7151522,0.1191920,& !<-- do not read this as a matrix
                        0.1804375,0.0721750,0.9503041/
        !component(1,:) = (/0.4124,0.3576,0.1805/)
        !component(2,:) = (/0.2126,0.7152,0.0722/) ! <-- Matrix will look like this and stuffs =w=
        !component(3,:) = (/0.0193,0.1192,0.9505/)
        !print*,component(1,:)
        !print*,component(2,:)
        !print*,component(3,:)
        
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
        REAL(kind=8) FUNCTION ft(t) result(res)
            IMPLICIT none
            REAL(kind=8), intent(in) :: t
            REAL(kind=8) :: A=(6./29.)**3.0,& 
                            B=(29./6.)**(2.),&
                            C=4./29.,& 
                            one_third=1.0/3.0 !to be honest, 1./3. is used twice.
            !IF (t > ((6./29.)**(3.))) THEN
            !    return (t**(1./3.))
            !ELSE
            !    return ((1./3.)*((29./6.)**(2.))*t+(4./29.))
            !ENDIF
            IF (t > A) THEN
                res = t**one_third
            ELSE
                res = one_third*B*t+C
            ENDIF
        ENDFUNCTION ft
    FUNCTION XYZ2LABCH(XYZ) RESULT(LABCH)
        REAL(kind=8), DIMENSION(3), INTENT(in) :: XYZ
        REAL(kind=8), DIMENSION(5) :: LABCH
        !#CIE XYZ tristimulus values of the reference white point
        REAL(kind=8) :: Xn=0.95047,Yn=1.00000,Zn=1.08883
        !L = 116.0*ft(Y/Yn) - 16
        !A = 500.0*(ft(X/Xn) - ft(Y/Yn))
        !B = 200.0*(ft(Y/Yn) - ft(Z/Zn))
        !C = sqrt((A**2)+(B**2))
        !H = atan2(B,A) ! *!*
        
        !XYZ   = [1X, 2Y, 3Z]
        !LABCH = [1L, 2A, 3B, 4C, 5H]
        
        LABCH(1) = 116.0*ft(XYZ(2)/Yn) - 16
        LABCH(2) = 500.0*(ft(XYZ(1)/Xn) - ft(XYZ(2)/Yn))
        LABCH(3) = 200.0*(ft(XYZ(2)/Yn) - ft(XYZ(3)/Zn))
        LABCH(4) = sqrt((LABCH(2)**2)+(LABCH(3)**2))
        LABCH(5) = atan2(LABCH(3),LABCH(2))
    ENDFUNCTION XYZ2LABCH
!-------------------------------------------------------------------------------    
        REAL(kind=8) FUNCTION rad(deg) result(res)
            IMPLICIT None
            REAL(kind=8), PARAMETER :: PI = 3.1415926535898
            REAL(kind=8), intent(in) :: deg
            res = (deg/180.0)*pi
        ENDFUNCTION rad
    !Given two colours, ref and poll, returns their deltaE (difference)
    !http://en.wikipedia.org/wiki/ΔE_(color_space)#CIEDE2000
    !REAL(kind=4) FUNCTION dE_single(ref_LABCH, poll_LABCH)
    REAL(kind=8) FUNCTION dE_single(rL, pL)
        IMPLICIT NONE
        !REAL(kind=4), DIMENSION(5), INTENT(in) :: ref_LABCH, poll_LABCH
        REAL(kind=8), DIMENSION(5), INTENT(in) :: rL, pL
        REAL(kind=8), PARAMETER :: PI = 3.1415926535898
        REAL(kind=8) :: dLp, Lbar, Cbar, a1p, a2p, C1p, C2p, dCp, Cbarp, h1_h2, &
                        h1p, h2p, dhp, dHHp, HHbarp, T, S_L, S_C, S_H, R_t, d0, R_C, dE00
        !2 = poll, 1 = ref
        !LABCH = [1L, 2a, 3b, 4C, 5h]
        
        !dLp = L2 - L1
        dLp = pL(1) - rL(1)
        !Lbar = (L1 + L2)/2
        Lbar = (rL(1) + pL(1))/2.0
        Cbar = (rL(4) + pL(4))/2.0
        
        a1p = rL(2) + (rL(2)/2.0)*(1.0-sqrt((Cbar**7)/((Cbar**7)+(25.0**7))))
        a2p = pL(2) + (pL(2)/2.0)*(1.0-sqrt((Cbar**7)/((Cbar**7)+(25.0**7))))
        C1p = sqrt((a1p**2) + (rL(3)**2))
        C2p = sqrt((a2p**2) + (pL(3)**2))
        dCp = C2p - C1p
        Cbarp = (C1p + C2p)/2.0
        
        IF (a1p == 0 .AND. rL(3) == 0) THEN
            h1p = 0
        ELSE
            h1p = atan2(rL(3),a1p)
            !print*,h1p,"h1p"
            IF (h1p < 0) THEN
                h1p = h1p + 2*pi
            ENDIF
            h1p = mod((180.0/pi)*h1p, 360.0)
            !print*,h1p,"h1p"
        ENDIF
        
        IF (a2p == 0 .AND. pL(3) == 0) THEN
            h2p = 0
        ELSE
            h2p = atan2(pL(3),a2p)
            !print*,h2p,"h2p"
            IF (h2p < 0) THEN
                h2p = h2p + 2*pi
            ENDIF
            h2p = mod((180.0/pi)*h2p, 360.0)
            !print*,h2p,"h2p"
        ENDIF
        
        IF (C1p == 0 .OR. C2p == 0) THEN
            dhp = 0
        ELSE
            h1_h2 = abs(h1p-h2p)
            IF (h1_h2 <= 180.0) THEN
                dhp = h2p - h1p
            ELSEIF ((h2p - h1p) > 180.0) THEN
                dhp = (h2p-h1p)-360.0
            ELSEIF ((h2p-h1p) < -180.0) THEN
                dhp = (h2p-h1p)+360.0
            ENDIF
        ENDIF
        !dH' *and not* dh'
        dHHp = 2.0*sqrt(C1p*C2p)*sin(rad(dhp/2.0))
        
        
        IF (C1p == 0 .OR. C2p == 0) THEN
            HHbarp = h1p + h2p
        ELSEIF (abs(h1p-h2p) <= 180.0) THEN
            HHbarp = (h1p + h2p)/2.0
        ELSEIF (abs(h1p - h2p) > 180.0 .AND. ((h1p+h2p) < 360.0)) THEN
            HHbarp = (h1p + h2p + 360.0)/2.0
        ELSEIF (abs(h1p - h2p) > 180.0 .AND. ((h1p + h2p) >= 360.0)) THEN
            HHbarp = (h1p + h2p - 360.0)/2.0
        ENDIF
        
        T = 1.0 - 0.17*cos(rad(HHbarp-30.0))+0.24*cos(rad(2*HHbarp))+0.32*cos(rad(3*HHbarp+6.0))-0.20*cos(rad(4*HHbarp-63.0))
        
        S_L = 1.0 + (0.015*((Lbar-50)**2))/(sqrt(20.0+((Lbar-50.0)**2)))
        S_C = 1.0 + 0.045*Cbarp
        S_H = 1.0 + 0.015*Cbarp*T
        
        R_C = 2*sqrt((Cbarp**7)/((Cbarp**7)+25.0**7))
        print*,R_C
        d0 = 30.0**(-(((HHbarp-275.0)/25.0)**2))
        print*,d0
        R_T = -R_C*sin(rad(2*d0))
        !R_T = -2.0*sqrt((Cbarp**7)/((Cbarp**7)+(25.0**7)))*sin(60.0**(-((HHbarp-275.0)/(25))**2))
        
        dE00 = sqrt(((dLp/(S_L))**2) + ((dCp/(S_C))**2) + ((dHp/S_h)**2) + R_T*((dCp/S_C)*(dHp/S_H)))
    ENDFUNCTION dE_single
    
!-------------------------------------------------------------------------------    
    !Given one colour (ref) and an array of colours (polls), return the 
    !FUNCTION dE_multi(ref_LABCH, [poll_LABCH])
        
    !ENDFUNCTION dE_multi
    
!-------------------------------------------------------------------------------    


ENDPROGRAM DeltaE
    
    
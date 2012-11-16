PROGRAM DeltaE
    IMPLICIT none
    REAL(kind=8), DIMENSION(3) :: sRGB1, sRGB2, XYZ1, XYZ2
    REAL(kind=8), DIMENSION(5) :: LABCH1, LABCH2
    INTEGER, DIMENSION(3) :: RGB1,RGB2
    REAL(kind=8) :: diff
    
    
    !RGB1 = (/136, 220, 73/) !reference colour
    !RGB2 = (/88, 99, 111/) !colour to compare to
    
    !sRGB1 = RGB2sRGB(RGB1)
    !sRGB2 = RGB2sRGB(RGB2)
    
    !XYZ1 = sRGB2XYZ(sRGB1)
    !XYZ2 = sRGB2XYZ(sRGB2)
    
    !LABCH1 = XYZ2LABCH(XYZ1)
    !LABCH2 = XYZ2LABCH(XYZ2)
    REAL(kind=8) :: l1,a1,b1,c1, l2,a2,b2,c2
    l1 = 50.000
    l2 = l1
    a1 = 2.6772
    a2 = 0.0000
    b1 = -79.7751
    b2 = -82.7485
    c1 = sqrt(a1**2 + b1**2)
    c2 = sqrt(a2**2 + b2**2)
    LABCH1 = (/l1,a1,b1,c1,0.0_8/)
    LABCH2 = (/l2,a2,b2,c2,0.0_8/)
    
    diff = dE_single(LABCH1, LABCH2)
    !PRINT*,""
    !PRINT*,"Colour 1:",RGB1
    !PRINT*,"Colour 1:",LABCH1
    
    !PRINT*,"Colour 2:",RGB2
    !PRINT*,"Colour 2:",LABCH2
    
    !PRINT*,"Colourdifference (lower=less diff): ",diff
    
CONTAINS
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
            REAL(kind=8) :: one_third=1.0/3.0 
            
            IF (t > (6./29.)**3.0) THEN
                res = t**one_third
            ELSE
                res = one_third*((29./6.)**(2.))*t + 4./29.
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
    !Given two colours, ref and poll, returns their deltaE (difference)
    !http://en.wikipedia.org/wiki/ΔE_(color_space)#CIEDE2000
    !REAL(kind=4) FUNCTION dE_single(ref_LABCH, poll_LABCH)
    REAL(kind=8) FUNCTION dE_single(rL, pL)
        IMPLICIT NONE
        !REAL(kind=4), DIMENSION(5), INTENT(in) :: ref_LABCH, poll_LABCH
        REAL(kind=8), DIMENSION(5), INTENT(in) :: rL, pL
        REAL(kind=8), PARAMETER :: PI = 3.1415926535898
        REAL(kind=8), PARAMETER :: deg2rad = PI/180.0
        REAL(kind=8), PARAMETER :: rad2deg = 180.0/PI
        
        !REAL(kind=8) :: dLp, Lbar, Cbar, a1p, a2p, C1s, C2s, dCp, Cbarp, h1_h2, &
        !                h1p, h2p, dhp, dHHp, HHbarp, T, S_L, S_C, S_H, R_t, d0, R_C, dE00
        REAL(KIND=8) :: c1s_ab, c2s_ab, Csbar_ab, G, a1p, a2p, C1p, C2p, h1p, h2p, h2p_h1p, &
                        dLp, dCp, dhp, dHHp, Lbarp, Cbarp, h1ph2p, h1p_h2p, hbarp, T, d0, &
                        R_C, S_L, S_C, S_H, R_T, K_L, K_C, K_H, dE00
        !2 = poll, 1 = ref
        !LABCH = [1L, 2a, 3b, 4C, 5h]
        
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
            h1p = h1p*rad2deg
        ENDIF
        IF (pL(3) == 0 .AND. a2p == 0) THEN
            h1p = 0
        ELSE 
            h2p = atan2(pL(3),a2p)
            IF (h2p < 0) THEN
                h2p = h2p + 2*pi
            ENDIF
            h2p = h2p*rad2deg
        ENDIF



        dLp = pL(1) - rL(1)
        dCp = C2p - C1p
        
        h2p_h1p = h2p - h1p
        IF (C1p*C2p == 0) THEN
            dhp = 0
        ELSEIF (abs(h2p_h1p) <= 180) THEN
            dhp = h2p_h1p
        ELSEIF (h2p_h1p > 180) THEN
            dhp = h2p_h1p - 360
        ELSEIF (h2p_h1p < -180) THEN
            dhp = h2p_h1p + 360
        ENDIF
        dHHp = 2.0*sqrt(C1p*C2p)*sin(dhp*deg2rad/2.0)
        
        
        Lbarp = (rL(1) + pL(1))/2.0
        Cbarp = (C1p + C2p)/2.0
        
        h1ph2p = h1p + h2p
        h1p_h2p = h1p - h2p
        IF (C1p*C2p == 0) THEN
            hbarp = h1ph2p
        ELSEIF (abs(h1p_h2p) <= 180.0) THEN
            hbarp = (h1ph2p)/2.0
        ELSEIF (abs(h1p_h2p) > 180.0 .AND. (h1ph2p < 360.0)) THEN
            hbarp = (h1ph2p + 360.0)/2.0
        ELSEIF (abs(h1p_h2p) > 180.0 .AND. (h1ph2p >= 360.0)) THEN
            hbarp = (h1ph2p - 360.0)/2.0
        ENDIF
        
        T = 1 - 0.17*cos(deg2rad*(hbarp-30.0)) + 0.24*cos(deg2rad*(2*hbarp)) + &
            0.32*cos(deg2rad*(3*hbarp + 6.0)) - 0.20*cos(deg2rad*(4*hbarp-63.0))
        
        d0 = 30.0*exp(-((hbarp-275.0)/25.0)**2.0)
        
        R_C = 2*sqrt((Cbarp**7.0)/((Cbarp**7)+25.0**7))
        
        S_L = 1.0 + ((0.015*((Lbarp-50)**2))/(sqrt(20.0+((Lbarp-50.0)**2))))
        S_C = 1.0 + 0.045*Cbarp
        S_H = 1.0 + 0.015*Cbarp*T
        
        R_T = -sin(deg2rad*(2.0*d0))*R_C
        
        K_L = 1
        K_C = 1
        K_H = 1
        
        dE00 = sqrt((dLp/(K_L*S_L))**2 + (dCp/(K_C*S_C))**2 + (dHHp/(K_H*S_H))**2 + R_T*(dCp/(K_C*S_C))*(dHHp/(K_H*S_H)))
        
        
        
        print"(Af0.4)","L1      ",rL(1)
        print"(Af0.4)","L2      ",PL(1)
        print"(Af0.4)","a1      ",rL(2)
        print"(Af0.4)","a2      ",pL(2)
        print"(Af0.4)","b1      ",rL(3)
        print"(Af0.4)","b2      ",pL(3)
        print"(Af0.4)","a1p     ",a1p
        print"(Af0.4)","a2p     ",a2p
        print"(Af0.4)","c1p     ",c1p
        print"(Af0.4)","c2p     ",c2p
        print"(Af0.4)","h1p     ",h1p
        print"(Af0.4)","h2p     ",h2p
        print"(Af0.4)","hbarp   ",hbarp
        print"(Af0.4)","G       ",G
        print"(Af0.4)","T       ",T
        print"(Af0.4)","S_L     ",S_L
        print"(Af0.4)","S_C     ",S_C
        print"(Af0.4)","S_H     ",S_H
        print"(Af0.4)","R_T     ",R_T
        print"(Af0.4)","de00    ",de00
        
    ENDFUNCTION dE_single
    
!-------------------------------------------------------------------------------    
    !Given one colour (ref) and an array of colours (polls), return the 
    !FUNCTION dE_multi(ref_LABCH, [poll_LABCH])
        
    !ENDFUNCTION dE_multi
    
!-------------------------------------------------------------------------------    


ENDPROGRAM DeltaE
    
    
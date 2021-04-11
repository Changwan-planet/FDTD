PROGRAM FDTD_1D_6
IMPLICIT NONE

INTEGER :: a = 0   !!X coordinate
INTEGER :: i = 0
INTEGER :: j = 0   
INTEGER :: time = 0

!Distance=C_S*Maxtrix_SIZE!
INTEGER, PARAMETER :: M_S=2.0d+5 !Maxtrix_SIZE, Place
INTEGER, PARAMETER :: N = 2.0d+5 !Number of time step
REAL*8,  PARAMETER :: f = 5.0d+6 !Frequency


!E-FIELD!
REAL*8 :: E_z(0:M_S) = 0.0d+0
REAL*8 :: J_z(0:M_S) = 0.0d+0!Independent source

!H-FIELD!
REAL*8 :: H_y(0:M_S) = 0.0d+0
REAL*8 :: M_y(0:M_S) = 0.0d+0 !Independent source 
!As for x-directed, z-polarized TEM mode, M_source=0 for all time

REAL*8, PARAMETER :: pi=Acos(-1.0)

!Material Propertiy!
REAL*8, PARAMETER :: Conduct = 0.0d+0  
REAL*8, PARAMETER :: Magloss = 0.0d+0              !! Magnetic loss  
REAL*8, PARAMETER :: Permit = 8.854 * (10.0)**(-12)!! Permitivity  Free space=1/(36*pi))*10**(-9) 
REAL*8, PARAMETER :: Permeat = (4.0*pi)*(10.0)**(-7)     !! Permeativity Free space=(4*pi)*10**(-7) 

!Lunar environment!
!Real, Parameter :: Conduct = 7*10**(-4)!! Conductivity  
!Real, Parameter :: Permit  = 2.87      !! Permitivity 
!Real, Parameter :: Permeat = 1.012     !! Permeativity  

!SI unit!
REAL*8, PARAMETER :: E_v = 1.0 / sqrt(Permit * Permeat)
REAL*8, PARAMETER :: C_S = 3.0        ![m]
REAL*8, PARAMETER :: T_D = 1.0d-9     ![s]

REAL*8  ::  pp1 = (1.0d-5)/T_D
REAL*8  ::  pp2 = (1.0d-4-1.d-5)/T_D
REAL*8  ::  pp3 = (1.0d-4)/T_D


!E-Filed Component!
REAL*8 :: C_a(0:M_S) = (1.0-((Conduct*T_D)/(2.0*Permit))) / (1.0+((Conduct*T_D)/(2.0*Permit)))
REAL*8 :: C_b(0:M_S) = T_D/(Permit*C_S) / (1.0+(Conduct*T_D)/(2.0*Permit))

!H-Field Component!
REAL*8 :: D_a(0:M_S) =(1.0-((Magloss*T_D)/(2.0*Permeat))) / (1.0+((Magloss*T_D)/(2.0*Permeat)))
REAL*8 :: D_b(0:M_S) = (T_D/(Permeat*C_S)) / (1.0+(Magloss*T_D/(2.0*Permeat)))

PRINT*,"Permit",Permit
PRINT*,"Permeat",Permeat
PRINT*,"Conduct",Conduct
PRINT*,"Magloss",Magloss
PRINT*,"================================"
PRINT*,"C_a",C_a(1)
PRINT*,"C_b",C_b(2)
PRINT*,"D_a",D_a(3)
PRINT*,"D_b",D_b(4)
PRINT*,"================================"
PRINT*,"M_S (Matrix Size)=",M_S
PRINT*,"C_S (Cell Size)=",C_S
PRINT*,"E_v (velocity) [m/s]=", E_v

PRINT*,"T_D (Time Difference) [s] ",T_D
PRINT*,"N (Number of time_step)",N
PRINT*,"Frequency [Hz]=", f
PRINT*,"Wavelength (Distance of one period)",E_v/f
!Print*,"Distance of a wavelength",E_v*N*T_D/1    !k=1
PRINT*,"Number of C_S of one period",(E_v/f)/C_S
PRINT*,"time_step of one period", 1.0/f 
PRINT*,"Distance = M_S*C_S [m]", M_S*C_S

!Path!
OPEN (UNIT=20,FILE='FDTD_1D_E.txt', STATUS='replace')
OPEN (UNIT=21,FILE="FDTD_1D_H.txt", STATUS='replace')
OPEN (UNIT=22,FILE="Source_1D.txt", STATUS='replace')

OPEN (UNIT=23,FILE="FDTD_1D_E1.txt",STATUS='replace')
OPEN (UNIT=24,FILE="FDTD_1D_H1.txt",STATUS='replace')

!In the Do-loop
    !Left size is n+1 and Right size is n in the time domain.

!Time!    
 DO time = 0,N    

!Elecric Field!
      DO i=1,M_S        
         E_z(i) = C_a(i)*E_z(i)+C_b(i)*(H_y(i)-H_y(i-1)-J_z(i)*C_S)      
         !n+1!  = !n!                                 !Time index   
      END DO

!First ouput
      IF (time==1.0/(f*T_D)) THEN
             DO i=0,M_S
                Write(23,*) E_z(i)
             END DO 
      END IF

!Final output  
      IF (time==N) THEN
             DO i=0,M_S
                WRITE(20,*) E_z(i)
                !Write(20,100) E_z(i)
                !100 FORMAT(E15.7)
                !Write(20,*) Log(E_z(i)**2)  !Divergence Check
             END DO 
      END IF 

!Source ( Boundary Condition )
!===========================================================================      
       !E_z(0) = 1-cos(2*pi*f*time)
       !Print*,"E_z(0)=", E_z(0)    
!===========================================================================   

!Long-time sources boudary condition(Cell_size=10000)
!===========================================================================
!      If (time<=(1.0d-4/T_D)) then
!          E_z(0) = 1-cos(2*pi*f*time*T_D)
!      Else
!          E_z(0) = 0
!      End if 
!===========================================================================

!One-time sources (Boundary condition)
!===========================================================================
      IF (time<=1.0/(f*T_D)) THEN
         E_z(0) = 1-COS(2*pi*f*time*T_D)
      ELSE 
         E_z(0) = 0
      END IF 
!===========================================================================




!Weigted Long-time Sources (Boundary sources)
!===========================================================================


!       If (time<=pp1) then
!          E_z(0) = Weight_fun1(time)*(1-cos(2*pi*f*time*T_D-pi))
       
!       Else if ((pp1<time).AND.(time<pp2)) then
!          E_z(0) = 1-cos(2*pi*f*time*T_D)

!       Else if ((pp2<=time).AND.(time<=pp3)) then
!          E_z(0) = Weight_fun2(time)*(1-cos(2*pi*f*time*T_D+pi))

!       Else 
!          E_z(0) = 0
     
!       End If

!===========================================================================

!One-time sources (Boundary condition)
!===========================================================================
IF (time<=1.0/(f*T_D)) THEN
    E_z(0) = 1-COS(2*pi*f*time*T_D)
END IF 

!ELSE 
!         E_z(0) = 0
!      End if 
!===========================================================================


WRITE(22,*) E_z(0)

!Magnetic Field 
      DO i=0,M_S-1
         H_y(i) = D_a(i)*H_y(i)+D_b(i)*(E_z(i+1)-E_z(i)-M_y(i)*C_S)   
         !n+2!  = !n+1!                                 !Time index
      END DO   

!First output
      IF(time==1.0/(f*T_D)) THEN
             DO i=0,M_S
                WRITE(24,*) H_y(i)
             END DO 
      END IF

!Final output
      IF(time==N) THEN
             DO i=0,M_S   
                WRITE(21,*) H_y(i)
                !Write(21,100) H_y(i)
                !Write(21,*) Log(H_y(j)**2) !Divergence Chek
             END DO
      END IF 

 END DO

CONTAINS

  FUNCTION Weight_Fun1 (time) 

           REAL*8 :: Weight_fun1
           INTEGER, INTENT (IN) :: time

            Weight_fun1=1.0/2.0*COS(time/1.0d+4*pi+pi)+1.0/2.0 
            !Weight_fun1=time/10000.0

  END FUNCTION Weight_Fun1

  FUNCTION Weight_Fun2 (time)
           
           REAL*8 :: Weight_fun2
           INTEGER, INTENT (IN) :: time

           Weight_fun2=1.0/2.0*COS(((2.0d+5-time)/1.0d+4)*pi+pi)+1.0/2.0
           !Weight_fun2=1/(200000-time)
            
  END FUNCTION Weight_Fun2

END PROGRAM FDTD_1D_6

PROGRAM FDTD_3D
IMPLICIT NONE

!!!!
INTEGER :: a = 0   !!X coordinate
INTEGER :: i = 0
INTEGER :: j = 0   
INTEGER :: k = 0
INTEGER :: time = 0


!     +++++++++++++++++++++++++++
!++++++Distance=C_S*Maxtrix_SIZE++++++
!     +++++++++++++++++++++++++++
INTEGER, PARAMETER :: M_S = 12.0d+1 !Maxtrix_SIZE, Place
                                   !This value is time-consuming.
                                   !If this is too small, it's unrealistic 
                                   ! but fast. 


INTEGER, PARAMETER :: N = 1.0d+3   !Number of time step
REAL*8,  PARAMETER :: f = 5.0d+6   !Frequency
!REAL*8,  PARAMETER :: f = 1.0d+0   !Frequency


!     +++++++++
!++++++E-FIELD++++++
!     +++++++++ 
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: E_x = 0.0d+0
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: E_y = 0.0d+0
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: E_z = 0.0d+0


REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: J_x = 0.0d+0    !Independent source
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: J_y = 0.0d+0    !Independent source
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: J_z = 0.0d+0    !Independent source


!     +++++++++
!++++++H-FIELD++++++
!     +++++++++
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: H_x = 0.0d+0
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: H_y = 0.0d+0
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: H_z = 0.0d+0


REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: M_x = 0.0d+0    !Independent source 
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: M_y = 0.0d+0    !Independent source 
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: M_z = 0.0d+0    !Independent source 


!As for x-directed, z-polarized TEM mode, M_source=0 for all time
REAL*8, PARAMETER :: pi=Acos(-1.0)

!     ++++++++++++++++++++
!++++++Material Propertiy++++++
!     ++++++++++++++++++++
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2), PARAMETER :: Conduct = 0.0d+0

REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2), PARAMETER :: Magloss = 0.0d+0        
   
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2), PARAMETER :: Permit = 8.854*(10.0)**(-12) 
!Permitivity  Free space=1/(36*pi))*10**(-9) 
  
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2), PARAMETER :: Permeat = (4.0*pi)*(10.0)**(-7) 
!Permeativity Free space=(4*pi)*10**(-7) 


!Lunar environment!
!Real, Parameter :: Conduct = 7*10**(-4)!! Conductivity  
!Real, Parameter :: Permit  = 2.87      !! Permitivity 
!Real, Parameter :: Permeat = 1.012     !! Permeativity  

!     +++++++++
!++++++SI unit++++++
!     +++++++++
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2)  :: E_v 
REAL*8, PARAMETER :: C_S = 1.0                              ![m]
                           !If this is too small, it's realistic but slow.
                           !If this is too big, it's fast but errorneous.

REAL*8, PARAMETER :: T_S = 1.0d-9                           ![s]


REAL*8  ::  pp1 = (1.0d-5)/T_S
REAL*8  ::  pp2 = (1.0d-4-1.d-5)/T_S
REAL*8  ::  pp3 = (1.0d-4)/T_S


!     +++++++++++++++++++
!++++++E-Filed Component++++++
!     +++++++++++++++++++    
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: C_a 
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: C_b 


!     +++++++++++++++++++
!++++++H-Field Component++++++
!     +++++++++++++++++++
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: D_a 
REAL*8, DIMENSION(0:M_S+2, 0:M_S+2, 0:M_S+2) :: D_b 


E_v = 1.0 / sqrt(Permit* Permeat)

C_a = (1.0-((Conduct*T_S)/(2.0*Permit))) / (1.0+((Conduct*T_S)/(2.0*Permit)))
C_b = (T_S/(Permit*C_S)) / (1.0+(Conduct*T_S/(2.0*Permit)))

D_a = (1.0-((Magloss*T_S)/(2.0*Permeat))) / (1.0+((Magloss*T_S)/(2.0*Permeat)))
D_b = (T_S/(Permeat*C_S)) / (1.0+(Magloss*T_S/(2.0*Permeat)))
 
    
PRINT*,"Permit(1,1,1)",Permit(1,1,1)
PRINT*,"Permeat(1,1,1)",Permeat(1,1,1)
PRINT*,"Conduct(1,1,1)",Conduct(1,1,1)
PRINT*,"Magloss(1,1,1)",Magloss(1,1,1)
PRINT*,"================================"
PRINT*,"C_a(0,0,0)",C_a(0,0,0)
PRINT*,"C_b(0,0,0)",C_b(0,0,0)
PRINT*,"D_a(1,1,1)",D_a(3,3,3)
PRINT*,"D_b(1,1,1)",D_b(4,4,4)
PRINT*,"================================"
PRINT*,"M_S (Matrix Size)=",M_S
PRINT*,"C_S (Cell Size)=",C_S
PRINT*,"E_v(1,1,1) (velocity) [m/s]=", E_v(1,1,1)

!Print*,"Wavelength=",E_v(1,1,1)*N*T_D/(C_S*k) 
!                         C(velocity of light) x N (number of time step) x T_D 
!         Wavelength = ----------------------------------------------------------
!                             C_S (Cell size) x k (number of wavelength)
!
!         d = C X N X T_D (Total wave from distance)
!         d/k (distance of a wavelength)
!         N x T_D    (travleing time)

PRINT*,"T_S (Time Difference) [s] ",T_S
PRINT*,"N (Number of time_step)",N
PRINT*,"Frequency [Hz]=", f
PRINT*,"Wavelength (Distance of one period)",E_v(1,1,1)/f
Print*,"Distance of a wavelength",E_v(1,1,1)*N*T_S/1                         !k=1
PRINT*,"Number of C_S of one period",(E_v(1,1,1)/f)/C_S
PRINT*,"time_step of one period", 1.0/f 
PRINT*,"Distance = M_S*C_S [m]", M_S*C_S


!     ++++++,
!++++++Path++++++
!     ++++++
OPEN (UNIT=20,FILE='FDTD_3D_E.txt', STATUS='replace')
OPEN (UNIT=21,FILE="FDTD_3D_H.txt", STATUS='replace')
OPEN (UNIT=22,FILE="Source_3D.txt", STATUS='replace')
OPEN (UNIT=23,FILE="FDTD_3D_E1.txt", STATUS='replace')
OPEN (UNIT=24,FILE="FDTD_3D_H1.txt", STATUS='replace')



!In the Do-loop
!Left size is n+1 and Right size is n in the time domain.

!     ++++++
!++++++Time++++++    
!     ++++++

DO time = 0,N    

!     +++++++++++++++
!++++++Elecric Field++++++
!     +++++++++++++++

      DO i = 0, M_S-2
         DO j = 0, M_S-2
            DO k = 0, M_S-2



                 !n+1!  = !n!                                 !Time index 
                   E_x (i+2,j+3,k+3) = C_a(i+2,j+3,k+3) * E_x(i+2,j+3,k+3) + &
                                      &C_b(i+2,j+3,k+3) * (H_z(i+2,j+4,k+3) - H_z(i+2,j+2,k+3)+ &
                                                          &H_y(i+2,j+3,k+2) - H_y(i+2,j+3,k+4)- &
                                                          &J_z(i+2,j+3,k+3) * C_S)





                 !n+1!  = !n!                                 !Time index 
                   E_y (i+1,j+4,k+3) = C_a(i+1,j+4,k+3) * E_y (i+1,j+4,k+3) + &
                                      &C_b(i+1,j+4,k+3) * (H_x(i+1,j+4,k+4) - H_x(i+1,j+4,k+2)+ &
                                                          &H_z(i,j+4,k+3) - H_z(i+2,j+4,k+3) - &
                                                          &J_z(i+1,j+4,k+3) * C_S)






                 !n+1!  = !n!                                 !Time index  
                   E_z (i+1,j+3,k+4) = C_a(i+1,j+3,k+4) * E_z(i+1,j+3,k+4) + &
                                      &C_b(i+1,j+3,k+4) * (H_y(i+2,j+3,k+4) - H_y(i,j+3,k+4)+ &
                                                         &H_x(i+1,j+2,k+4) - H_x(i+1,j+4,k+4)- &
                                                         &J_z(i+1,j+3,k+4) * C_S)                    


            END DO
         END DO
      END DO   

!     +++++++++++++      
!++++++First ouput++++++
!     +++++++++++++
       IF (time==1.5/(f*T_S)) THEN
!      IF (time==3*1.0/(f*T_S)) THEN
!      IF (time==1) THEN
         ! 1.0/ (f*T_S)= 200 time step = 1 Wavelenth 
         
             DO i=0,M_S+2
                DO j=0,M_S+2
                   DO k=0,M_S+2
            
                      WRITE(20,*) E_x(i,j,k), E_y(i,j,k), E_z(i,j,k)
               
                   END DO  
                END DO 
             END DO 

      END IF

!     ++++++++++++++      
!++++++Final output++++++
!     ++++++++++++++
      IF (time==N) THEN
             i=0
             j=0
             k=0
             DO i=0,M_S+2
                DO j=0,M_S+2
                   DO k=0,M_S+2

                      WRITE(23,*) E_x(i,j,k), E_y(i,j,k), E_z(i,j,k)
                !Write(20,100) E_z(i)
                !100 FORMAT(E15.7)
                !Write(20,*) Log(E_z(i)**2)  !Divergence Check
                   END DO 
                END DO
             END DO   
      END IF 

!Source ( Boundary Condition )
!===========================================================================      
!       E_z(0,0,0) = 1-cos(2*pi*f*time)
!       Print*,"E_z(0,0,0)=", E_z(0,0,0)    
!===========================================================================   

!Long-time sources boudary condition(Cell_size=10000)
!===========================================================================
!     If (time<=(1.0d-4/T_D)) then
!          E_z(0) = 1-cos(2*pi*f*time*T_D)
!     Else
!          E_z(0) = 0
!     End if 
!===========================================================================


!One-time sources (Boundary condition)
!===========================================================================
      IF (time<=1.0/(f*T_S)) THEN
         E_x(2,3,3) = 1-COS(2*pi*f*time*T_S)
 !        E_y(1,4,3) = 1-COS(2*pi*f*time*T_S)
 !        J_z(10,10,10) = 1-COS(2*pi*f*time*T_S)
      ELSE 
         E_x(2,3,3) = 0
 !        E_y(1,4,3) = 0
 !        J_z(10,10,10) = 0
      END IF  
!===========================================================================


!Weigted Long-time Sources (Boundary sources)
!===========================================================================
!     IF (time<=pp1) THEN
!          E_z(1,3,4) = Weight_fun1(time)*(1-cos(2*pi*f*time*T_S-pi))
!       
!     ELSE IF ((pp1<time).AND.(time<pp2)) THEN
!          E_z(1,3,4) = 1-cos(2*pi*f*time*T_S)
!
!     ELSE IF ((pp2<=time).AND.(time<=pp3)) THEN
!          E_z(1,3,4) = Weight_fun2(time)*(1-cos(2*pi*f*time*T_S+pi))
!
!     ELSE
!          E_z(1,3,4) = 0
!     
!     END IF
!===========================================================================

!Hard source 1, 

!Hard source 2, a lowpass Gaussian pulse with finite dc content
!         IF(time
!     E_z(1,3,4) = E_0 * EXP(-(time-time_0)/time_decay)**2

!PRINT *, E_x(2,3,3)
WRITE(22,*) E_x(2,3,3), E_y(1,4,3), E_z(1,3,4)

!WRITE(22,*) E_x(2,3,3), E_y(0,4,3), E_z(1,3,4)
!WRITE(22,*) H_x(1,4,4), H_y(2,3,4), H_z(2,4,3)
 
!     ++++++++++++++++    
!++++++Magnetic Field++++++
!     ++++++++++++++++
     DO i = 0, M_S-3
         DO j = 0, M_S-3
            DO k = 0, M_S-3
  



                             !n+2!  = !n+1!              !Time index
                   H_x(i+1,j+4,k+4) = D_a(i+1,j+4,k+4) * H_x(i+1,j+4,k+4) + &
                                    & D_b(i+1,j+4,k+4) * (E_y(i+1,j+4,k+5) - E_y(i+1,j+4,k+3) + &
                                                        &E_z(i+1,j+3,k+4) - E_z(i+1,j+5,k+4) - &
                                                        &M_x(i+1,j+4,k+4) * C_S)   





                   H_y(i+2,j+3,k+4) = D_a(i+2,j+3,k+4) * H_y(i+2,j+3,k+4)+ &
                                    & D_b(i+2,j+3,k+4) * (E_z(i+3,j+3,k+4) - E_z(i+1,j+3,k+4) + &
                                                         &E_x(i+2,j+3,k+3) - E_x(i+2,j+3,k+5) - &
                                                         &M_y(i+2,j+3,k+4)*C_S)   






                   H_z(i+2,j+4,k+3) = D_a(i+2,j+4,k+3) * H_z(i+2,j+4,k+3) + &
                                    & D_b(i+2,j+4,k+3) * (E_x(i+2,j+5,k+3) - E_x(i+2,j+3,k+3) + & 
                                                         &E_y(i+1,j+4,k+3) - E_y(i+3,j+4,k+3)- &
                                                         &M_z(i+2,j+4,k+3) * C_S)   
 


            END DO   
         END DO
      END DO   

!     ++++++++++++++ 
!++++++First output++++++
!     ++++++++++++++      
      IF(time==1.0/(f*T_S)) THEN
!       IF(time==400) THEN
            DO i=0,M_S+2
               DO j=0,M_S+2
                 DO k=0,M_S+2

                    WRITE(24,*) H_x(i,j,k), H_y(i,j,k), H_z(i,j,k)
             
                 END DO
               END DO
             END DO
!                     WRITE(24,*) H_x(1,4,4), H_y(2,3,4), H_z(2,4,3)
      END IF

!     ++++++++++++++
!++++++Final output++++++
!     ++++++++++++++
      IF(time==N) THEN
             i=0
             j=0
             k=0
             DO i=0,M_S
                DO j=0,M_S
                   DO k=0,M_S
                      WRITE(21,*) H_x(i,j,k), H_y(i,j,k), H_z(i,j,k)
                      !Write(21,100) H_y(i)
                      !Write(21,*) Log(H_y(j)**2) !Divergence Chek
             
                   END DO
                END DO
             END DO   
      END IF 

 END DO

CONTAINS

  FUNCTION Weight_Fun1 (time) 

           REAL*8 :: Weight_fun1
           INTEGER, INTENT (IN) :: time

            Weight_fun1=1.0/2.0*cos(time/1.0d+4*pi+pi)+1.0/2.0 
            !Weight_fun1=time/10000.0

  END FUNCTION Weight_Fun1

  FUNCTION Weight_Fun2 (time)
           
           REAL*8 :: Weight_fun2
           INTEGER, INTENT (IN) :: time

           Weight_fun2=1.0/2.0*cos(((2.0d+5-time)/1.0d+4)*pi+pi)+1.0/2.0
           !Weight_fun2=1/(200000-time)
            
  END FUNCTION Weight_Fun2

END PROGRAM FDTD_3D

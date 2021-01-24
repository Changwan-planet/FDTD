Program FDTD_1D_3
Implicit none

Integer :: a = 0   !!X coordinate
Integer :: i = 0
Integer :: j = 0   
Integer :: time = 0

!Distance=C_S*Maxtrix_SIZE!
Integer, Parameter :: M_S=2.0d+4 !Maxtrix_SIZE, Place
Integer, Parameter :: N = 2.0d+5 !Number of time step
Real*8,  Parameter :: f = 5.0d+6 !Frequency


!E-FIELD!
Real*8 :: E_z(0:M_S) = 0.0d+0
Real*8 :: J_z(0:M_S) = 0.0d+0!Independent source

!H-FIELD!
Real*8 :: H_y(0:M_S) = 0.0d+0
Real*8 :: M_y(0:M_S) = 0.0d+0 !Independent source 
!As for x-directed, z-polarized TEM mode, M_source=0 for all time

Real*8, Parameter :: pi=Acos(-1.0)

!Material Propertiy!
Real*8, Parameter :: Conduct = 0.0d+0  
Real*8, Parameter :: Magloss = 0.0d+0              !! Magnetic loss  
Real*8, Parameter :: Permit = 8.854 * (10.0)**(-12)!! Permitivity  Free space=1/(36*pi))*10**(-9) 
Real*8, Parameter :: Permeat = (4.0*pi)*(10.0)**(-7)     !! Permeativity Free space=(4*pi)*10**(-7) 

!Lunar environment!
!Real, Parameter :: Conduct = 7*10**(-4)!! Conductivity  
!Real, Parameter :: Permit  = 2.87      !! Permitivity 
!Real, Parameter :: Permeat = 1.012     !! Permeativity  

!SI unit!
Real*8, Parameter :: E_v = 1.0 / sqrt(Permit * Permeat)
Real*8, Parameter :: C_S = 3.0        ![m]
Real*8, Parameter :: T_D = 1.0d-9     ![s]

Real*8  ::  pp1 = (1.0d-5)/T_D
Real*8  ::  pp2 = (1.0d-4-1.d-5)/T_D
Real*8  ::  pp3 = (1.0d-4)/T_D


!E-Filed Component!
Real*8 :: C_a(0:M_S) = (1.0-((Conduct*T_D)/(2.0*Permit))) / (1.0+((Conduct*T_D)/(2.0*Permit)))
Real*8 :: C_b(0:M_S) = T_D/(Permit*C_S) / (1.0+(Conduct*T_D)/(2.0*Permit))

!H-Field Component!
Real*8 :: D_a(0:M_S) =(1.0-((Magloss*T_D)/(2.0*Permeat))) / (1.0+((Magloss*T_D)/(2.0*Permeat)))
Real*8 :: D_b(0:M_S) = (T_D/(Permeat*C_S)) / (1.0+(Magloss*T_D/(2.0*Permeat)))

Print*,"Permit",Permit
Print*,"Permeat",Permeat
Print*,"Conduct",Conduct
Print*,"Magloss",Magloss
Print*,"================================"
Print*,"C_a",C_a(1)
Print*,"C_b",C_b(2)
Print*,"D_a",D_a(3)
Print*,"D_b",D_b(4)
Print*,"================================"
Print*,"M_S (Matrix Size)=",M_S
Print*,"C_S (Cell Size)=",C_S
Print*,"E_v (velocity) [m/s]=", E_v

!Print*,"Wavelength=",E_v*N*T_D/(C_S*k) 
!                         C(velocity of light) x N (number of time step) x T_D 
!         Wavelength = ----------------------------------------------------------
!                             C_S (Cell size) x k (number of wavelength)
!
!         d = C X N X T_D (Total wave from distance)
!         d/k (distance of a wavelength)
!         N x T_D    (travleing time)

Print*,"T_D (Time Difference) [s] ",T_D
Print*,"N (Number of time_step)",N
Print*,"Frequency [Hz]=", f
Print*,"Wavelength (Distance of one period)",E_v/f
!Print*,"Distance of a wavelength",E_v*N*T_D/1    !k=1
Print*,"Number of C_S of one period",(E_v/f)/C_S
Print*,"time_step of one period", 1.0/f 
Print*,"Distance = M_S*C_S [m]", M_S*C_S

!Path!
Open (unit=20,file='FDTD_1D_E.txt', status='replace')
Open (unit=21,file="FDTD_1D_H.txt", status='replace')
Open (unit=22,file="Source.txt", status='replace')

Open (unit=23,file="FDTD_1D_E1.txt",status='replace')
Open (unit=24,file="FDTD_1D_H1.txt",status='replace')

!In the Do-loop
    !Left size is n+1 and Right size is n in the time domain.

!Time!    
 Do time = 0,N    

!Elecric Field!
      Do i=1,M_S        
         E_z(i) = C_a(i)*E_z(i)+C_b(i)*(H_y(i)-H_y(i-1)-J_z(i)*C_S)      
         !n+1!  = !n!                                 !Time index   
      End do

!First ouput
      If (time==1.0/(f*T_D)) then
             Do i=0,M_S
                Write(23,*) E_z(i)
             End do 
      End if

!Final output  
      If (time==N) then
             Do i=0,M_S
                Write(20,*) E_z(i)
                !Write(20,100) E_z(i)
                !100 FORMAT(E15.7)
                !Write(20,*) Log(E_z(i)**2)  !Divergence Check
             End do 
      End if 

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
!      If (time<=1.0/(f*T_D)) then
!         E_z(0) = 1-cos(2*pi*f*time*T_D)
!      Else 
!         E_z(0) = 0
!      End if 
!===========================================================================

!Weigted Long-time Sources (Boundary sources)
!===========================================================================


       If (time<=pp1) then
          E_z(0) = Weight_fun1(time)*(1-cos(2*pi*f*time*T_D-pi))
       
       Else if ((pp1<time).AND.(time<pp2)) then
          E_z(0) = 1-cos(2*pi*f*time*T_D)

       Else if ((pp2<=time).AND.(time<=pp3)) then
          E_z(0) = Weight_fun2(time)*(1-cos(2*pi*f*time*T_D+pi))

       Else 
          E_z(0) = 0
     
       End If

!===========================================================================


Write(22,*) E_z(0)

!Magnetic Field 
      Do i=0,M_S-1
         H_y(i) = D_a(i)*H_y(i)+D_b(i)*(E_z(i+1)-E_z(i)-M_y(i)*C_S)   
         !n+2!  = !n+1!                                 !Time index
      End do   

!First output
      If(time==1.0/(f*T_D)) then
             Do i=0,M_S
                Write(24,*) H_y(i)
             End do 
      End if

!Final output
      If(time==N) then
             Do i=0,M_S   
                Write(21,*) H_y(i)
                !Write(21,100) H_y(i)
                !Write(21,*) Log(H_y(j)**2) !Divergence Chek
             End do
      End if 

 End do

CONTAINS

  Function Weight_Fun1 (time) 

           Real*8 :: Weight_fun1
           INTEGER, INTENT (IN) :: time

            Weight_fun1=cos(time*pi-pi) 
            !Weight_fun1=time/10000.0

  End Function Weight_Fun1

  Function Weight_Fun2 (time)
           
           Real*8 :: Weight_fun2
           INTEGER, INTENT (IN) :: time

           Weight_fun2=cos((200000-time)*pi+pi)
           !Weight_fun2=1/(200000-time)
            
  End Function Weight_Fun2

End program FDTD_1D_3

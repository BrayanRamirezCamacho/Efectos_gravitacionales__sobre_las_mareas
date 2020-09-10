!!El objetivo de este programa es el de aproximar la solución de una
!!ecuación diferencial (E.D) utilizado el método de Runge-Kutta de
!!orden 4 (RK4).
!!-------------------------------------------------------------------

!Establecemos un módulo con los objetos a utilizar en el programa.
MODULE objetos

  IMPLICIT NONE

  SAVE

  !!Kind con una precisión de 16 para mayor resolución en el sistema.
  INTEGER,PARAMETER::d = SELECTED_REAL_KIND( 16 )

  !Diccionario de parámetros
  REAL(d)::G = 6.674E-11                  !!Constante gravitacional de la Tierra ((N*m^2)/kg^2).


  !!Tipo público para simular objetos celestes (estáticos).
  TYPE,PUBLIC::celeste

     REAL(d)::m                 !!Masa del cuerpo celeste (kg).
     REAL(d)::r                 !!Radio del cuerpo celeste (m).

  END TYPE celeste

  !!Tipo público para almacenar los datos de la simulación.
   TYPE,PUBLIC::simulacion

     INTEGER::n                  !!Número de datos a simular.
     INTEGER::ttot               !!Tiempo total de la simulación.
     REAL(d)::dt                 !!Ancho de paso.

  END TYPE simulacion
  
CONTAINS

  !!Establecemos las E.D a resolver que dictan la posición de la marea en cada eje (x,y,z).

  !!Para el eje x:
  FUNCTION ax(x,r1,m1,r2,m2)
    REAL(d),INTENT(IN)::x,r1,m1,r2,m2      !!Argumentos de la ecuación.
    REAL(d)::ax                              !!E.D a resolver.
    
    ax= -((G*m1)/r1**3)*x -G*m2*x*(r2**2+r1**2-2*x*r2)**(-3/2) +G*m2*r2*(r2**2+r1**2-2*x*r2)**(-3/2)  - ((G*m2)/r2**2)

  END FUNCTION ax

  !!Para el eje y:
  FUNCTION ay(y,x,r1,m1,r2,m2)
    REAL(d),INTENT(IN)::y,x,r1,m1,r2,m2    !!Argumentos de la ecuación.
    REAL(d)::ay                              !!E.D a resolver.

    ay = -((G*m1)/r1**3)*y -G*m2*y*(r2**2+r1**2-2*x*r2)**(-3/2)
  
  END FUNCTION ay

  FUNCTION az(z,x,r1,m1,r2,m2)
    REAL(d),INTENT(IN)::z,x,r1,m1,r2,m2    !!Argumentos de la ecuación.
    REAL(d)::az                              !!E.D a resolver.

    az = -((G*m1)/r1**3)*z -G*m2*z*(r2**2+r1**2-2*x*r2)**(-3/2)

  END FUNCTION az

  !Establecemos una función para calcular la magnitud de un vector.
  FUNCTION mag(x,y,z)
    REAL(d),INTENT(IN)::x,y,z          !!Argumentos de la función.
    REAL(d)::mag                 !!Función que calcula la magnitud de la función.

    mag = (x**2 + y**2 + z**2)**(1./2. )

  END FUNCTION mag

  

END MODULE objetos

PROGRAM mareas

  !!Llamamos al módulo objetos
  USE objetos

  !!Diccionario de variables
  IMPLICIT NONE

  REAL(d),ALLOCATABLE::kr1(:),kr2(:),kr3(:),kr4(:),kr5(:)    !!Vectores k para posiciones.
  REAL(d),ALLOCATABLE::kv1(:),kv2(:),kv3(:),kv4(:),kv5(:)    !!Vectores k para velocidades.
  REAL(d),ALLOCATABLE::r(:,:)                           !!Vector de posición.
  REAL(d),ALLOCATABLE::rt(:)                            !!Vector de posiión en un instante t.
  REAL(d),ALLOCATABLE::ru(:)                            !!Vector unitario de posición en un instante t.
  REAL(d),ALLOCATABLE::v(:,:)                           !!Vector de velocidad.
  REAL(d),ALLOCATABLE::fm(:)                            !!Fuerza gravitacional sobre la marea en un instante t.
  REAL(d),ALLOCATABLE::hm(:)                            !!Altura de la marea en un instante t.
  REAL(d)::l                                            !!Longitud del "resorte".
  REAL(d)::s                                            !!Diferencia de posición con respecto a la posición inicial de la marea.
  INTEGER::i                                             !!Contador.
  TYPE(celeste)::tr                                      !!Objeto del tipo celeste al cual simularemos sus mareas.
  TYPE(celeste)::lun                                     !!Objeto del tipo celeste cuyas fuerzas gravitacionales ejercen "fuerza de marea" sobre el otro.
  TYPE(simulacion)::sim                                 !!Objeto del tipo simulación.
  REAL(d)::a

  !!Establecemos los valores de cada cuerpo celeste.
  tr%m = 5.972E24            !!m1
  tr%r = 6371E3              !!r1
  lun%m   = 7.349E22         !!m2
  lun%r   = 384.4E6          !!r2
  
  !!Establecemos los parámetros de nuestra aproximación.
  PRINT*,"Introduce el número de valores a simular:"
  READ*,sim%n
  PRINT*,"Introduce el tiempo de duración de la simulación:"
  READ*,sim%ttot

  !!Establecemos un ancho de paso a partir de nuestros parámetros.
  sim%dt = REAL(sim%ttot)/REAL(sim%n)

  !!Alocamos los arreglos según el número de datos:
  !Posición
  ALLOCATE(r(0:2,0:sim%n))
  ALLOCATE(rt(0:2))
  ALLOCATE(ru(0:2))
  ALLOCATE(kr1(0:2))
  ALLOCATE(kr2(0:2))
  ALLOCATE(kr3(0:2))
  ALLOCATE(kr4(0:2))
  ALLOCATE(kr5(0:2))

  !Velocidad
  ALLOCATE(v(0:2,0:sim%n))
  ALLOCATE(kv1(0:2))
  ALLOCATE(kv2(0:2))
  ALLOCATE(kv3(0:2))
  ALLOCATE(kv4(0:2))
  ALLOCATE(kv5(0:2))

  !Fuerza
  ALLOCATE(fm(0:2))

  !Altura de mareas
  ALLOCATE(hm(0:sim%n))

  
  !!Fijamos las condiciones iniciales:
  !Posición inicial
  r(0,0) = tr%r/(3**(1./2.))
  r(1,0) = tr%r/(3**(1./2.))
  r(2,0) = tr%r/(3**(1./2.))

  !Velocidad inicial
  v(0,0) = 0.0
  v(1,0) = 0.0
  v(2,0) = 0.0

  !Altura inicial de la marea
  hm(0) = mag(r(0,0),r(1,0),r(2,0)) - tr%r

  !!Para resolver la ecuación, podemos visualizar el movimiento de mareas como
  !!un oscilador, cuya longitud original l es el radio de la tierra.
  DO i=1,sim%n
     rt = r(:,i-1)
     l = mag(rt(0),rt(1),rt(2))
     s = l - tr%r
     ru = rt/l
     fm(0) = ax(ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(1) = ay(ru(1)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(2) = az(ru(2)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m

     kr1 = sim%dt*v(:,i-1)
     kv1 = sim%dt*fm(:)!/tr%m

     rt = r(:,i-1) + (1./2.)*kr1
     l = mag(rt(0),rt(1),rt(2))
     s = l - tr%r
     ru = rt/l
     fm(0) = ax(ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(1) = ay(ru(1)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(2) = az(ru(2)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m

     kr2 = sim%dt*(v(:,i-1) + (1./2.)*kv1)
     kv2 = sim%dt*fm(:)!/tr%m

     rt = r(:,i-1) + (1./2.)*kr2
     l = mag(rt(0),rt(1),rt(2))
     s = l - tr%r
     ru = rt/l
     fm(0) = ax(ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(1) = ay(ru(1)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(2) = az(ru(2)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m

     kr3 = sim%dt*(v(:,i-1) + (1./2.)*kv2)
     kv3 = sim%dt*fm(:)!/tr%m

     rt = r(:,i-1) + kr3
     l = mag(rt(0),rt(1),rt(2))
     s = l - tr%r
     ru = rt/l
     fm(0) = ax(ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(1) = ay(ru(1)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m
     fm(2) = az(ru(2)*s,ru(0)*s,l,tr%m,lun%r,lun%m)!*tr%m

     kr4 = sim%dt*(v(:,i-1) + kv3)
     kv4 = sim%dt*fm(:)!/tr%m

     kr5 = kr1 + 2*kr2 + 2*kr3 + kr4
     kv5 = kv1 + 2*kv2 + 2*kv3 + kv4

    !Establecemos el nuevo valor de posición y velocidad en x utilizando el método de RK4:
     r(:,i) = r(:,i-1) + (1/6.)*kr5
     v(:,i) = v(:,i-1) + (1/6.)*kv5
     l = mag(rt(0),rt(1),rt(2))
     s = l - tr%r
     hm(i) = s
  !   PRINT*,i,r(1,i),hm(i),r(1,i)-hm(i),v(1,i)
  END DO

   !Escribimos los valores en un archivo:
     OPEN(1,FILE="res",STATUS="UNKNOWN",ACTION="WRITE")
     DO i=0,sim%n
        !   WRITE(1,*)sim%dt*i,r(:,i),v(:,i),hm(i)
        WRITE(1,*)sim%dt*i,hm(i)
     END DO
     CLOSE(1)
          
  
  
END PROGRAM mareas

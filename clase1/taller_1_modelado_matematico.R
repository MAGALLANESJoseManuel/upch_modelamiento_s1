# instalar las librerias necesarias 
install.packages('deSolve')
# cargar libreria 
library(deSolve)

# crear una funcion que resuelva el sistema de ecuacines 

sir.model.closed <- function (t, x, params) {
  S <- x[1]
  I <- x[2]
  R <- x[3]
  
  with( #we can simplify code using "with"
    as.list(params), #this argument to "with" lets us use the variable names
     { #the system of rate equations
      dS <- -beta*S*I
      dI <- beta*S*I-gamma*I
      dR <- gamma*I
      dx <- c(dS,dI,dR) #combine results into a single vector dx
      list(dx) #return result as a list
       }
    )
}

# ceramos lo neceario para correr mi modelo y resolverlo 
times <- seq(0,120,by=2) #function seq returns a sequence
params <- c(beta=0.3,gamma=1/7) #function c "c"ombines values into a vector
xstart <- c(S=9999/10000,I=1/10000,R=0) #initial conditions

out <- as.data.frame(ode(xstart,times,sir.model.closed,params))


op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1)) #set graphical parameters
plot(I~time,data=out,type='b') #plot the I variable against time
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T) #re-set graphical parameters
plot(I~S,data=out,type='b',yaxt='n',xlab='S') #plot phase portrait
par(op)
sol <- as.data.frame(deSolve::ode(xstart,times,sir.model.closed,params))

with(sol, {
    # plotting the time series of susceptibles:
    plot(time, S, type = "l", col = "blue",
         xlab = "time (days)", ylab = "number of people")
    # adding the time series of infectious:
    lines(time, I, col = "red")
    # adding the time series of recovered:
    lines(time, R, col = "green")
})

# adding a legend:
legend("right", c("susceptibles", "infectious", "recovered"),
       col = c("blue", "red", "green"), lty = 1, bty = "n")
####ejercicio 1


 params2 <- c(beta=0.5,gamma=1/7)
 params3 <- c(beta=0.3,gamma=0.0857)
 params4 <- c(beta=0.6,gamma=2/7)




# Resolviendo y trazando vemos que:
#  1. Aumentar β (curva roja) hace que el brote ocurra más rápido 
 # (el crecimiento inicial en el aumento de individuos infectados 
 # es más rápido y el pico es más alto); 
 # además, el tamaño total de la epidemia es mayor (es decir
 #                                                                                                                                                                                   el área bajo la curva).
# 2. Disminuir γ (a un nivel que logre el mismo R0 utilizado en el caso anterior; curva azul) de manera similar
# conduce a una epidemia más grande, sin embargo, el despegue inicial no es tan rápido y el pico se retrasa en comparación
# con el aumento de β. El retrato de fases es el mismo que en el caso anterior, lo que significa que la
# epidemia sigue el mismo camino, pero lo hace más lentamente.
# 3. Combinando estos dos resultados, anticipamos el resultado que obtenemos al aumentar β y γ por un
# factor común de dos (curva verde): Como el R0 no se ve afectado, el retrato de fases no cambia,
# pero la epidemia ocurre más rápido (ya que el valor de β ha aumentado).



out2 <- as.data.frame(lsoda(xstart, times, sir.model.closed, params2))
out3 <- as.data.frame(lsoda(xstart, times, sir.model.closed, params3))
out4 <- as.data.frame(lsoda(xstart, times, sir.model.closed, params4))
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
plot(I ~ time, data=out, type='l', lwd=2, ylim=c(0, 0.45))
lines(I ~ time, data=out2, type='l', lwd=2, col='red')
lines(I ~ time, data=out3, type='l', lwd=2, col='blue')
lines(I ~ time, data=out4, type='l', lwd=2, col='green')
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
plot(I ~ S, data=out, type='l', lwd=2, log='x', yaxt='n', xlab='S', ylim=c(0.01, 0.45), xlim=c(0.01, 1))
lines(I ~ S, data=out2, type='l', lwd=2, col='red')
lines(I ~ S, data=out3, type='l', lwd=2, col='blue')
lines(I ~ S, data=out4, type='l', lwd=2, col='green')
legend('topleft', col=c('red','blue','green'), legend=c('Incrementar beta', 'Disminuir gamma', 'Incrementar beta y gamma por un factor común de dos'), lwd=2)
par(op)

#EJercicio 2

# Aquí volvemos a nuestro conjunto inicial de valores de parámetros (β = 0.3, γ = 1/7).
# Lo que cambiamos son las condiciones iniciales. Consideraremos un caso con un número 
# cien veces mayor de infectados iniciales, pero sin inmunidad en la población (S = 9900/10000, I = 100/10000, R = 0),
# un caso con el número original de infectados iniciales y algunos individuos inmunes
# (S = 9899/10000, I = 1/10000, R = 100), y un caso con un gran número de infectados 
# iniciales y un gran número de individuos inmunes (S = 6400/10000, I = 100/10000, R = 3500).

 xstart5 <- c(S=9900/10000,I=100/10000,R=0)
 xstart6 <- c(S=9899/10000,I=1/10000,R=100)
 xstart7 <- c(S=6400/10000,I=100/10000,R=3500)
 out5 <- as.data.frame(lsoda(xstart5,times,sir.model.closed,params))
 out6 <- as.data.frame(lsoda(xstart6,times,sir.model.closed,params))
 out7 <- as.data.frame(lsoda(xstart7,times,sir.model.closed,params))
 op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
 plot(I~time,data=out,type='l', lwd=2, ylim=c(0, 0.3))
 lines(I~time,data=out5,type='l', lwd=2, col='red')
 lines(I~time,data=out6,type='l', lwd=2, col='blue')
 lines(I~time,data=out7,type='l', lwd=2, col='green')
 par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
 plot(I~S,data=out,type='l',lwd=2, log='x',yaxt='n',xlab='S', ylim=c(0.01, 0.3), xlim=c(0.01, 1))
 lines(I~S,data=out5,type='l', lwd=2, col='red')
 lines(I~S,data=out6,type='l', lwd=2, col='blue')
 lines(I~S,data=out7,type='l', lwd=2, col='green')
 legend('topleft', col=c('red','blue','green'), 
        legend=c('I(0) large', 'R(0) large', 'R(0) very large'))
        par(op)




# Lo que vemos es que aumentar el número inicial de infectados incrementa la 
# tasa de despegue (curva roja comparada con la curva negra), que aumentar el 
# número de inmunes en una cantidad similar tiene muy poco efecto (curva azul 
# comparada con la curva negra), y que aumentar el número de inmunes en gran medida reduce 
# significativamente el brote, es decir, es más pequeño y más lento que en una 
# población totalmente susceptible (curva verde comparada con la curva negra). 
# En este caso, hemos comenzado con una fracción de inmunes (35%) 
# aproximadamente la mitad de la que se requeriría para lograr la inmunidad de 
# rebaño (1 - 1/R0 = 1 - γ/β ≈ 0.524).
# Lo que vemos es que aumentar el número inicial de infectados incrementa la tasa 
# de despegue (curva roja comparada con la curva negra), que aumentar el número 
# de inmunes en una cantidad similar tiene muy poco efecto (curva azul comparada 
# con la curva negra), y que aumentar el número de inmunes en gran medida reduce 
# significativamente el brote, es decir, es más pequeño y más lento que en una 
# población totalmente susceptible (curva verde comparada con la curva negra). 
# En este caso, hemos comenzado con una fracción de inmunes (35%) aproximadamente 
# la mitad de la que se requeriría para lograr la inmunidad de rebaño (1 - 1/R0 = 1 - γ/β ≈ 0.524).


# Ejercicio 3   

# Modifique los códigos dados para estudiar la dinámica de un modelo SIR demográficamente abierto.
# Para estudiar una población demográficamente abierta, añadimos una tasa de natalidad per cápita (µ). 
# necesitamos modificar nuestra función de derivadas para incluir esta tasa en cada ecuación y 
# modificar nuestro vector de parámetros para reflejar los diferentes valores de µ que deseamos estudiar.
# Primero, escribimos una nueva función para las derivadas.Modifique los códigos dados para estudiar
# la dinámica de un modelo SIR demográficamente abierto.
# Para estudiar una población demográficamente abierta, añadimos una tasa de natalidad per cápita (µ). 
#Necesitamos modificar nuestra función de derivadas para incluir esta tasa en cada ecuación y modificar
#nuestro vector de parámetros para reflejar los diferentes valores de µ que deseamos estudiar.
# Primero, escribimos una nueva función para las derivadas.


 sir.model.open <- function (t, x, params) { #here we begin a function with three arguments
   S <- x[1] #create local variable S, the first element of x
   I <- x[2] #create local variable I
   R <- x[3] #create local variable R
   with( #we can simplify code using "with"
     as.list(params), #this argument to "with" lets us use the variable names
     { #the system of rate equations
       dS <- mu*(S+I+R) -beta*S*I - mu*S
       dI <- beta*S*I-gamma*I - mu*I
       dR <- gamma*I - mu*R
       dx <- c(dS,dI,dR) #combine results into a single vector dx
       list(dx) #return result as a list
       }
     )
 }
 
 
 # Ahora modificamos el vector de parámetros de entrada. Dado que γ = 0.7, podemos
 # imaginar que estamos modelando una infección aguda donde los parámetros tienen 
 # unidades de días, es decir, un período infeccioso de 7 días. Estudiaremos tres niveles 
 # de µ: µ = 1/365, µ = 1/3650, µ = 1/36500. Estos corresponderían a duraciones 
 # de vida que para los vertebrados son (1 año), intermedia (10 años) o de larga duración (100 años).
 
 
  params8 <- c(beta=0.3,gamma=1/7, mu=1/52)
  params9 <- c(beta=0.3,gamma=1/7, mu=1/520)
  params10 <- c(beta=0.3,gamma=1/7, mu=1/5200)
  
  
   out8 <- as.data.frame(lsoda(xstart,times,sir.model.open,params8))
   out9 <- as.data.frame(lsoda(xstart,times,sir.model.open,params9))
   out10 <- as.data.frame(lsoda(xstart,times,sir.model.open,params10))
   op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))
  plot(I~time,data=out,type='l', lwd=2, ylim=c(0, 0.3))
   lines(I~time,data=out8,type='l', lwd=2, col='red')
   lines(I~time,data=out9,type='l', lwd=2, col='blue')
   lines(I~time,data=out10,type='l', lwd=2, col='green')
   par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)
   plot(I~S,data=out,type='l',lwd=2, log='x',yaxt='n',xlab='S', ylim=c(0.01, 0.3), xlim=c(0.01, 1))
   lines(I~S,data=out8,type='l', lwd=2, col='red')
   lines(I~S,data=out9,type='l', lwd=2, col='blue')
   lines(I~S,data=out10,type='l', lwd=2, col='green')
  legend('topright', col=c('red','blue','green'), 
         legend=c('short lifespan', 'intermediate lifespan', 'lo'),
         lty = 1, bty = "n")
  par(op)
  
  out8 <- as.data.frame(lsoda(xstart,times,sir.model.open,params8))

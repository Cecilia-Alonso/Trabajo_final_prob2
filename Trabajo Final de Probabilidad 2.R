# Trabajo de Probabilidad II: Simulación del Teorema 9 de (Bollobas, 2013)------

library(igraph)
library(tidyverse)
library(flow)

options(scipen=999)# Para elimiar notación científica

# Inicialmente analizamos una realización, para luego simular ----

# En línea con la demostración del teorema, supondremos que ω(n) no es demasiado grande, digamos ω(n) ≤ log n
# y n es lo suficientemente grande para garantizar que ω(n) ≥ 10
n <- 25000
wn <- log(n)-1/n
pln <- (log(n)-wn)/n
pun <- (log(n)+wn)/n

#Comprobaremos que Gpℓ esta desconectado y Gpu esta conectado

set.seed(42)
G_pl <- erdos.renyi.game(n, pln,type ="gnp")
is.connected(G_pl) # Grafo de pln esta desconectado, por lo que en el modelo Erdos Renyi es una funcion de umbral inferior
G_pu <- erdos.renyi.game(n, pun,type ="gnp") 
is.connected(G_pu) # Grafo de pun esta conectado, por lo que en el modelo Erdos Renyi es una funcion de umbral superior

#Analizamos qué ocurre con una probabilidad un poco mayores/menores a pln y pun
G_pl_menos <- erdos.renyi.game(n, 0,type ="gnp")
is.connected(G_pl_menos) #desconectado
G_pl_mas <- erdos.renyi.game(n, pln+0.0001,type ="gnp")
is.connected(G_pl_mas) #desconectado

G_pu_menos <- erdos.renyi.game(n, pun-0.0001,type ="gnp")
is.connected(G_pl_menos) #desconectado
G_pu_mas <- erdos.renyi.game(n, pun+0.0001,type ="gnp")
is.connected(G_pu_mas) #conectado
G_pu_mass <- erdos.renyi.game(n, pun+0.01,type ="gnp")
is.connected(G_pu_mass) #conectado

# Simulación-----

# Debemos elegir la cantidad de simulaciones a realizar ("s" que implican cambios de semilla), 
# el "n" del teorema, y la cantidad de probabilidades en el recorrido de 0 a (pun*1.1) que queremos evaluar
# para posteriormente hacer un gráfico de las probabilidades contra la proporción de grafos conectados. 

n_sim_seed_p <- function(s, n, cantp){ 
  conectado <- matrix(NA,cantp+1,2)
  conect_seed <- matrix(NA, 1,s)
  wn <- log(n)-1/n
  pln <- (log(n)-wn)/n
  pun <- (log(n)+wn)/n
  i=0
  for (p in seq(from=0, to=(pun*1.1), by=(pun*1.1)/cantp)) { 
    i=i+1
    for (j in 1:s) { 
      set.seed(j)
      G_p <- erdos.renyi.game(n, p,type ="gnp")
      conect_seed[1,j]  <-  as.numeric(is.connected(G_p))
    }
    conectado[i,1] <- p
    conectado[i,2] <- mean(conect_seed)
  }
  conectado2 <- as.data.frame(conectado)
  names(conectado2) <- c("p", "prom G_p conectado")
  return(conectado2)
}

#Esquema de la función anterior
flow_view(n_sim_seed_p)

#Aplico la función generada para determinado valor de "s", "n" y "cantp"
n <- 40000
s <- 10
cantp <- 100
conectado <- n_sim_seed_p(s,n,cantp) #s, n, cantp
pln <- (log(n)-(log(n)-1/n))/n
pun <- (log(n)+(log(n)-1/n))/n
G_pl <- erdos.renyi.game(n, pln,type ="gnp")
G_pu <- erdos.renyi.game(n, pun,type ="gnp")
con_G_pl <- as.numeric(is.connected(G_pl))
con_G_pu <- as.numeric(is.connected(G_pu))

# Grafico los resultados
conectado %>%
  ggplot(aes(x = p, y = `prom G_p conectado`)) +
  geom_point() +
  geom_point(aes(x = pln, y = con_G_pl), shape = 23, size = 3, fill="red",col="red")+
  geom_point(aes(x = pun, y = con_G_pu), shape = 23, size = 3, fill="blue", col="blue")+
  theme(axis.ticks = element_line(size = 2), legend.title = element_blank(),legend.position="bottom",legend.text=element_text(size=12,face="bold"), axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=12)) +
  ggtitle(paste("Resultados de la simulación n =",n, ", s=",s, ", cantp=", cantp)) + xlab("Probabilidad") + ylab("Proporción de G_p conectados")+
  annotate("text", x = pln, y = 0.05, label = "pln")+
  annotate("text", x = pun, y = 0.95, label = "pun")

# Extra: Otras funciones vinculadas -----

# Debemos elegir el "n" inicial y el "n" final con el que queremos calcular "wn", "pln", "pun"
# y evaluar si Gpℓ esta desconectado y Gpu esta conectado
n_sim <- function(n1, n){
  conectado <- matrix(NA, n-n1+1,6)
  for (i in n1:n) {
    conectado[i-n1+1,1] <-i
    wn <- log(i-n1+1)-1/i
    conectado[i-n1+1,2] <-wn
    pln <- (log(i)-wn)/i
    conectado[i-n1+1,3] <-pln
    pun <- (log(i)+wn)/i
    conectado[i-n1+1,4] <-pun
    G_pl <- erdos.renyi.game(i, pln,type ="gnp")
    conectado[i-n1+1,5] <- is.connected(G_pl)
    G_pu <- erdos.renyi.game(i, pun,type ="gnp")
    conectado[i-n1+1,6] <- is.connected(G_pu)
    }
  conectado <- as.data.frame(conectado)
  names(conectado) <- c("n", "wn", "pln", "pun", "G_pl conectado", "G_pu conectado")
  return(conectado)
}

extra1 <- n_sim(100,200)

# A la función anterior se le introduce la posibilidad de generar varias realizaciones con distinta semilla y analizar el porcentaje en las que resultó conectado el grafo con pln y pun. 
n_sim_seed <- function(s, n1, n){ # s es la cantidad de cambios de semilla, y n va de n1 a n, p es la cantidad de probabilidades 
  conectado <- matrix(NA, n-n1+1,6)
  conect_seed_G_pl <- matrix(NA, 1,s)
  conect_seed_G_pu <- matrix(NA, 1,s)
  for (i in n1:n) { #este loop va a ir cambiando el n
    conectado[i-n1+1,1] <-i
    wn <- log(i-n1+1)-1/i
    conectado[i-n1+1,2] <-wn
    pln <- (log(i)-wn)/i
    pun <- (log(i)+wn)/i
    conectado[i-n1+1,3] <-pln
    conectado[i-n1+1,4] <-pun
    
    for (j in 1:s) { #esto va a ir cambiando las semillas
      set.seed(s+j)
      G_pl <- erdos.renyi.game(i, pln,type ="gnp")
      conect_seed_G_pl[1,j]  <-  as.numeric(is.connected(G_pl))
      G_pu <- erdos.renyi.game(i, pun,type ="gnp")
      conect_seed_G_pu[1,j]  <-  as.numeric(is.connected(G_pu))
    }
    
    conectado[i-n1+1,5] <- mean(conect_seed_G_pl)
    conectado[i-n1+1,6] <- mean(conect_seed_G_pu)
  }
  conectado2 <- as.data.frame(conectado)
  names(conectado2) <- c("n", "wn", "pln", "pun", "prom G_pl conectado", "prom G_pu conectado")#,"G_pl trans", "G_pu trans"
  return(conectado2)
}

extra2 <- n_sim_seed(10,25000,25010)


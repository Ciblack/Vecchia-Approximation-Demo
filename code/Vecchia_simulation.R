library(ggplot2)
library(gridExtra)
library(akima)
library(plotly)
library(metR)
library(GpGp)
library(Rfast)

set.seed(42)
ny=90
nx=90
grille <- data.frame(
  X = rep(seq(from = 0, to = 100, length.out = nx), ny),
  Y = rep(seq(from = 0, to = 100, length.out = ny), each = nx)
)


Covariance_matrix=GpGp::matern45_isotropic(covparms =c(50,5,0.05) ,locs =as.matrix(grille) )  # definition de la matrice de variance covariance,  le vecteur covparams contient les parametres de la fonction : seuil, portée, pepite
chol=Rfast::cholesky(Covariance_matrix,parallel = TRUE) # decomposition de cholesky
Y=rnorm(dim(grille)[1]) # simulation de gausiennes independantes
grille$Z=as.vector(Y%*%chol) # integration de la variable simulée Z, dans le jeu de données "grille". 


#plot 3D 
interp_data <- interp(x = grille$X, y = grille$Y, z = grille$Z, duplicate = "mean") # Créer la surface 3D
plot_ly(x = interp_data$x, y = interp_data$y, z = interp_data$z, type = "surface")


ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = Z), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7)))  +
  labs(title = "Champ gaussien simulé")+
  theme(plot.title = element_text(hjust = 0.5))

rm(nx,ny,Y,chol,Covariance_matrix,interp_data)




## Données observées
N=sample(dim(grille)[1],dim(grille)[1]*0.1)
p1=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = Z), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7))) +
  geom_point(
    data = grille[N,],
    aes(x = X, y = Y),
    color = "black", # Border color
    fill = NA,       # Transparent fill
    shape = 21       # Shape with both border and fill
  ) +
  labs(title = "Données observées")+
  theme(plot.title = element_text(hjust = 0.5))
p1

## krigeage exact
result_exact_krig=exact_krig(param = c(50,5,0.05),locs_obs =grille[N,-3] ,Y =grille[N,3] ,locs_pred = grille[,-3])

## krigeage vecchia methode 1, C=5
result_vecchia_5=krig_vecchia_method1(param = c(50,5,0.05),locs_obs =grille[N,-3] ,Y =grille[N,3] ,locs_pred = grille[,-3],C= 5)

## krigeage vecchia methode 1, C=50
result_vecchia_50=krig_vecchia_method1(param = c(50,5,0.05),locs_obs =grille[N,-3] ,Y =grille[N,3] ,locs_pred = grille[,-3],C = 50)

## krigeage vecchia methode 1, C=100
result_vecchia_100=krig_vecchia_method1(param = c(50,5,0.05),locs_obs =grille[N,-3] ,Y =grille[N,3] ,locs_pred = grille[,-3],C = 100)


# krigeage vecchia methode 2, C=5 
result_vecchia2_5=krig_vecchia_method2(param = c(50,5,0.05),locs_obs =grille[N,-3] ,Y =grille[N,3] ,locs_pred = grille[,-3],C= 5)

# krigeage vecchia methode 2, C=50
result_vecchia2_50=krig_vecchia_method2(param = c(50,5,0.05),locs_obs =grille[N,-3] ,Y =grille[N,3] ,locs_pred = grille[,-3],C = 50)

# krigeage vecchia methode 2, C=100
result_vecchia2_100=krig_vecchia_method2(param = c(50,5,0.05),locs_obs =grille[N,-3] ,Y =grille[N,3] ,locs_pred = grille[,-3],C = 100)



### plots 
p2=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = result_exact_krig), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7)))  +
  labs(title = "Prediction méthode exacte")+
  theme(plot.title = element_text(hjust = 0.5))


p3=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = as.numeric(result_vecchia_5)), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7)))  +
  labs(title = "Prediction Vecchia méthode 1, C = 5") +
  theme(plot.title = element_text(hjust = 0.5))

p4=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = as.numeric(result_vecchia_50)), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7))) +
  labs(title = "Prediction Vecchia méthode 1, C = 50")+
  theme(plot.title = element_text(hjust = 0.5))


p5=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = as.numeric(result_vecchia_100)), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7)))  +
  labs(title = "Prediction Vecchia méthode 1, C = 100")+
  theme(plot.title = element_text(hjust = 0.5))




p6=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = as.numeric(result_vecchia2_5)), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7)))  +
  labs(title = "Prediction Vecchia méthode 2, C = 5") +
  theme(plot.title = element_text(hjust = 0.5))

p7=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = as.numeric(result_vecchia2_50)), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7))) +
  labs(title = "Prediction Vecchia méthode 2, C = 50")+
  theme(plot.title = element_text(hjust = 0.5))


p8=ggplot() +
  geom_contour_fill(data = grille, aes(x = X, y = Y, z = as.numeric(result_vecchia2_100)), alpha = 0.5) + 
  scale_fill_gradientn(colours = rev(rainbow(7)))  +
  labs(title = "Prediction Vecchia méthode 2, C = 100")+
  theme(plot.title = element_text(hjust = 0.5))


grid.arrange(p1, p2,p3,p4, ncol = 2)

grid.arrange(p1, p2,p3,p5, ncol = 2)

grid.arrange(p1, p2,p6,p8, ncol = 2)

grid.arrange(p1, p2,p6,p3, ncol = 2)





grille$Obersvation=0
grille$Obersvation[N]=1

grille$exact_krig=result_exact_krig
grille$Vecchia1_5=as.numeric(result_vecchia_5)
grille$Vecchia1_50=as.numeric(result_vecchia_50)
grille$Vecchia1_100=as.numeric(result_vecchia_100)

grille$Vecchia2_5=as.numeric(result_vecchia2_5)
grille$Vecchia2_50=as.numeric(result_vecchia2_50)
grille$Vecchia2_100=as.numeric(result_vecchia2_100)


save(grille,file="Full_data")

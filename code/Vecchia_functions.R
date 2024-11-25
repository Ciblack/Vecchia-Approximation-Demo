library(Matrix)
library(GpGp)
library(Rfast)

Get_lower=function(param,C,coords){
  
  V=GpGp::matern45_isotropic(param,locs = as.matrix(coords)) # gives the covariance matrix for specific parameters
  
  
  NN=find_ordered_nn(coords,C)
  result= lapply(1:dim(coords)[1], FUN = function(i) get_B_F(i, V[rev(NN[i,][!is.na(NN[i,])]),rev(NN[i,][!is.na(NN[i,])])],NN=NN))  # no Parallel
  #result= future.apply::future_lapply(1:dim(coords)[1], FUN = function(i) get_B_F(i, V[rev(NN[i,][!is.na(NN[i,])]),rev(NN[i,][!is.na(NN[i,])])],NN=NN))  #  Parallel
  
  
  B_matrix = do.call(rbind, lapply(result, function(x) x$Bsub))
  
  
  B_matrix <-  Matrix(B_matrix , sparse = T) # Construction de Bs
  
  F_values <- sapply(result, function(x) x$Fsub)                          # Construction de Fs       
  inv_F <- Matrix( diag(1 / (F_values)) , sparse = T)          
  
  return(list(B_matrix=B_matrix,inv_F=inv_F))
}

get_B_F=function(position,sub_V,NN){
  n=dim(NN)[1]
  if(position==1){
    Bsub <- rep(0,n)  
    Bsub[1]=1
    Fsub=sub_V
  }else{
    q=dim(sub_V)[1]
    chol=(Rfast::cholesky(sub_V[1:(q-1),1:(q-1)])) 
    
    ## Construction de Bsi
    X=t(sub_V[1:(q-1),q])%*%forwardsolve(chol,diag(q-1),upper.tri = TRUE)         
    Bsi=(X%*%forwardsolve( t(chol),x = diag(q-1)))
    
    Bsub <- rep(0,n)                   
    Bsub[position]=1
    Bsub[ rev( NN[position,][!is.na(NN[position,])][-1] )  ]=-Bsi
    
    ## Construction de Fsi
    X=forwardsolve(t(chol),sub_V[1:(q-1),q],upper.tri = FALSE)       
    N=sub_V[1:(q-1),q]%*%forwardsolve((chol),X,upper.tri = TRUE)
    
    Fsub=sub_V[q,q]- N
  }
  return(list(Bsub=Bsub,Fsub=Fsub))
}


exact_krig=function(param,locs_obs,Y,locs_pred){
  
  all_locs=rbind(locs_obs,locs_pred)
  Covariance_matrix=GpGp::matern45_isotropic(covparms =param ,locs =as.matrix(all_locs) )  # definition de la matrice de variance covariance,  le vecteur covparams contient les parametres de la fonction : seuil, portée, pepite
  
  npred=dim(locs_pred)[1]
  nobs=dim(locs_obs)[1]
  sigma12=Covariance_matrix[(nobs+1):(nobs+npred),1:nobs]
  sigma=Covariance_matrix[1:nobs,1:nobs]
  inv_sigma=solve(sigma)
  
  result=sigma12%*% inv_sigma%*%Y
  return(result)
}


krig_vecchia_method1=function(param,locs_obs,Y,locs_pred,C){
  
  lower=Get_lower(param = param,C = C,coords = locs_obs)
  
  Q=t(lower$B_matrix)%*%lower$inv_F%*%lower$B_matrix
  all_locs=rbind(locs_obs,locs_pred)
  Covariance_matrix=GpGp::matern45_isotropic(covparms =param ,locs =as.matrix(all_locs) )  # definition de la matrice de variance covariance,  le vecteur covparams contient les parametres de la fonction : seuil, portée, pepite
  
  npred=dim(locs_pred)[1]
  nobs=dim(locs_obs)[1]
  sigma12=Covariance_matrix[(nobs+1):(nobs+npred),1:nobs]
  
  result=sigma12%*% Q%*%Y
  return(result)
}


krig_vecchia_method2=function(param,locs_obs,Y,locs_pred,C){

  
  coords_fulldata=rbind(locs_pred,locs_obs)
  lower=Get_lower(param = param,C = C,coords = coords_fulldata)
  Q=t(lower$B_matrix)%*%lower$inv_F%*%(lower$B_matrix)
  U=t(sqrt(lower$inv_F)%*%lower$B_matrix)
  W=Q[1: dim(locs_pred)[1],1:dim(locs_pred)[1]]
  U_r.=U[(dim(locs_pred)[1]+1 ):dim(coords_fulldata)[1],] 
  U_l.=U[1: dim(locs_pred)[1],]  
  
  return(   -solve(W)%*%U_l.%*%t(U_r.)%*%( Y   ))
  
}


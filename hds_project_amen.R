###########################
#Funzioni utili
###########################


threshold <- function(Y,thresh){  
   
   #34*(rmean[i]+rmean[j])/(ncol(Y))
   
   #rmean <- rowMeans(Y)
   
   for(i in 1:nrow(Y)){
      
      for(j in 1:ncol(Y)){
         
         if(Y[i,j]>thresh ){
            
            Y[i,j] = 1;
            
            
         }
         else{
            Y[i,j] = 0
         }
         
      }
      
   }
   
   
   # for(j in 1:ncol(Y)){
   # 
   #   Y[j,j] = NA;
   # 
   # 
   # }
   
   return(Y);
   
}


threshold_soft <- function(Y,thresh){  
   
   for(i in 1:nrow(Y)){
      
      for(j in 1:ncol(Y)){
         
         if(Y[i,j]<thresh ){
            
            Y[i,j] = 0;
            
            
         }
         
      }
      
   }
   
   
   return(Y);
   
}



simulate <- function(ame_model,psi,n_samples,...){  # simulate statistics
   
   statistics = list();
   
   
   for(i in 1:n_samples){
      
      Y_sample <- simY_bin(EZ = ame_model$EZ,rho = 1);
      statistics[[i]] <- psi(Y_sample,...);
      
   }
   
   return(statistics);
   
};



simulate_communities <- function(ame_model,community_detection_algo,nsim,...){
   
   
   Y_sample <- simY_bin(EZ = ame_model$EZ,rho = 1);
   
   c2 <- community_detection_algo(graph_from_adjacency_matrix(as.matrix(
      Y_sample),mode = 'undirected'),...);
   
   m <- length(c2$names);
   
   I <- diag(m);
   
   comm_matrix <- NULL;
   
   for(j in 1:length(c2$membership)){
      
      
      comm_matrix <- rbind(comm_matrix,I[c2$membership[j],]);
      
   }
   
   
   
   
   comms <- list(comm_matrix);
   
   mu <- comm_matrix;
   
   
   for(i in 2:nsim){
      
      Y_sample <- simY_bin(EZ = ame_model$EZ,rho = 1);
      
      c2 <- community_detection_algo(graph_from_adjacency_matrix(as.matrix(
         Y_sample),mode = 'undirected'),...);
      
      
      # if(m < max(c2$membership)){
      #    
      #    m <- max(c2$membership);
      #    
      # }
      
      comm_matrix <- NULL;
      # I <- diag(m);
      
      for(j in 1:length(c2$membership)){
         
         
         comm_matrix <- rbind(comm_matrix,I[c2$membership[j],]);
         
      }
      
      
      
      
      
      comms[[i]] <- comm_matrix ;
      
      mu <- mu+comm_matrix;
      
      
   }
   
   mu <- mu/nsim;
   
   k = 1;
   
   while(k<m & mean(mu[,k] != 0) != 0){
      
      k <- k+1
      
   }
   
   #return(comms);
   
   rownames(mu) <- c2$names;
   
   return(list(comms,mu[,1:k-1]));
   
   
}


###########################
# Giallo zafferano
###########################


#install.packages('ergm')
library(ergm)
#install.packages('statnet')
library(statnet)
Y = read.csv('sociogiallo.csv')
# rownames(Y) <- Y[,1];
# Y = Y[,-1];
rownames(Y) <- colnames(Y);
#colnames(Y) <- rownames(Y)
#View(Y)



View(Y)

plot.sociomatrix(as.matrix(Y),diaglab = FALSE,cex.lab = 1,srt=90,
                 drawlab = F);

Y_net <- as.network(as.matrix(Y));

p_Y <- ergm(Y_net ~ edges+triangles ,reference = ~ Bernoulli,symmetric = TRUE);

library(ergm.count)

p_Y <- ergm(Y_net ~ edges,reference = ~ Poisson,symmetric = TRUE);

#install.packages('amen')
library(amen)

p.ame.Y <- ame(as.matrix(Y,dimnames = list( colnames(Y),colnames(Y)  )),
               family = 'ord',symmetric = TRUE);
summary(p.ame.Y)

p.ame.Y$EZ    # estimated expectation of Z matrix, so:

lambda <- exp(p.ame.Y$EZ)   # is the matrix containing the estimates of the lambda poiss parameters


#plot( apply(Y,1,mean,na.rm=TRUE) , type = 'l' ) 




#################################################
# Choosing the threshold for binarization
#################################################

Y_base = read.csv('sociogiallo.csv');
rownames(Y) <- colnames(Y);
edge_density_giallo = NULL;

for(thresh in 1:max(Y_base)){

   Y = Y_base;
   # rownames(Y) <- Y[,1];
   # Y = Y[,-1];
   #rownames(Y) <- colnames(Y);
   #colnames(Y) <- rownames(Y)
   #View(Y)
   
   Y = threshold(Y,thresh)
   
   
   edge_density_giallo = c(edge_density_giallo,mean(Y==1));
   
   
}


plot(edge_density_giallo)

edge_density_giallo[12]

Y_giallo = threshold(Y_base,12);

rownames(Y_giallo) <- colnames(Y_giallo);






Y_base = read.csv('sociotasty.csv');
rownames(Y) <- colnames(Y);
edge_density_tasty = NULL;

for(thresh in 1:max(Y_base)){
   
   Y = Y_base;
   # rownames(Y) <- Y[,1];
   # Y = Y[,-1];
   #rownames(Y) <- colnames(Y);
   #colnames(Y) <- rownames(Y)
   #View(Y)
   
   Y = threshold(Y,thresh)
   
   
   edge_density_tasty = c(edge_density_tasty,mean(Y==1));
   
   
}


plot(edge_density_tasty)

edge_density_tasty[12]

Y_tasty = threshold(Y_base,12);

rownames(Y_tasty) <- colnames(Y_tasty);



#########################################################
# Fitting AME models 
#########################################################

library(amen)

Y_base = read.csv('sociogiallo.csv');
rownames(Y_base) <- colnames(Y_base);


Y_giallo = Y_base;

for(j in 1:ncol(Y_giallo)){
   
   Y_giallo[j,j] = NA;
}


p.ame.Y_giallo <- ame(as.matrix(Y_giallo,dimnames = list( colnames(Y_giallo),colnames(Y_giallo)  )),
               family ='bin',symmetric = TRUE,R = 2);


dev.new()
par(mfrow= c(1,1))

plot(p.ame.Y_giallo$U,xlab = 'U1',ylab = 'U2')
text(x = p.ame.Y_giallo$U[,1],y = p.ame.Y_giallo$U[,2] ,labels=colnames(Y_giallo), cex=0.9, font=2)


dev.new()
par(mfrow= c(1,1))

netplot(Y = as.matrix(Y_giallo),plot.iso=FALSE,plotnames=TRUE)


dev.new()
par(mfrow= c(1,1))

circplot(Y = as.matrix(p.ame.Y_giallo$YPM>0.4),plotnames=TRUE)



#sim_Y_giallo <- simY_bin(EZ = p.ame.Y_giallo$EZ,rho = 1)   # rho = 1 for symmetric simulation




library(igraph)
#install.packages('lsa')
library(lsa)



c2 <- cluster_walktrap(graph_from_adjacency_matrix(as.matrix(Y_giallo),
                                                       mode = 'undirected')
                           );   # comm det for the real graph


#simulation of the communities
comm_giallo <- simulate_communities(p.ame.Y_giallo,cluster_leading_eigen,100);


print(comm_giallo[[2]])   # mean community matrix

interaction_comm_matrix <- 0

for(sample in comm_giallo[[1]]){   # computing community interaction matrix
   
   for(r in 1:nrow(sample)){
      
      interaction_comm_matrix <- interaction_comm_matrix+
         outer(sample[,r],sample[,r])/length(comm_giallo[[1]]);
   }
   
   
}

colnames(interaction_comm_matrix) <- colnames(Y_giallo);
rownames(interaction_comm_matrix) <- colnames(Y_giallo);


plot.sociomatrix(threshold_soft(as.matrix(interaction_comm_matrix),0.3),diaglab = FALSE,cex.lab = 1,srt=90,
                 drawlab = F);



View(as.matrix(interaction_comm_matrix))    # element (i,j) estimates the probability of i and j of
                                 # being in the same community according to ame


################################
# Tasty
################################



Y_base = read.csv('sociotasty.csv');
rownames(Y_base) <- colnames(Y_base);

Y_tasty = Y_base/sum(Y_base);

for(j in 1:ncol(Y_tasty)){
   
   Y_tasty[j,j] = NA;
}


p.ame.Y_tasty <- ame(Y = as.matrix(Y_tasty,dimnames = list( colnames(Y_tasty),colnames(Y_tasty)  )),
                     family = 'bin',symmetric = T,dcor = FALSE,R = 2);


dev.new()
par(mfrow= c(1,1))

plot(p.ame.Y_tasty$U,xlab = 'U1',ylab = 'U2')
text(x = p.ame.Y_tasty$U[,1],y = p.ame.Y_tasty$U[,2] ,labels=colnames(Y_tasty), cex=0.9, font=2)


dev.new()
par(mfrow= c(1,1))

#circplot(Y = as.matrix(Y_giallo),plotnames=TRUE)
netplot(Y = as.matrix(Y_tasty),plot.iso=FALSE,plotnames=TRUE)


dev.new()
par(mfrow= c(1,1))

netplot(Y = as.matrix(Y_tasty),plot.iso=FALSE,plotnames=TRUE)





comm_tasty <- simulate_communities(p.ame.Y_tasty,cluster_leading_eigen,100);

# mean_comm_giallo <- 0;
# 
# for(j in 1:length(comms_giallo)){
#    
#    mean_comm_giallo = mean_comm_giallo+comms_giallo[[j]]/length(comms_giallo)
#    
# }

print(comm_tasty[[2]])   #mean community matrix

interaction_comm_matrix_tasty <- 0

for(sample in comm_tasty[[1]]){   # community interaction matrix
   
   for(r in 1:nrow(sample)){
      
      interaction_comm_matrix_tasty <- interaction_comm_matrix_tasty+
         outer(sample[,r],sample[,r])/length(comm_tasty[[1]]);
   }
   
   
}

colnames(interaction_comm_matrix_tasty) <- colnames(Y_tasty);
rownames(interaction_comm_matrix_tasty) <- colnames(Y_tasty);


plot.sociomatrix(threshold_soft(as.matrix(interaction_comm_matrix_tasty),0.3),diaglab = FALSE,cex.lab = 1,srt=90,
                 drawlab = F);



View(as.matrix(interaction_comm_matrix_tasty))    # element (i,j) estimates the probability of i and j of
# being in the same community according to ame



# Y_tasty = read.csv('tasty_matrix.csv')
# rownames(Y_tasty) <- Y_tasty[,1];
# Y_tasty = Y_tasty[,-1];
# colnames(Y_tasty) <- rownames(Y_tasty)
# View(Y_tasty)
# 
# 
# 
# 
# p.ame.Y_tasty <- ame(as.matrix(Y_tasty,dimnames = list( colnames(Y_tasty),colnames(Y_tasty)  )),
#                family = 'bin',symmetric = TRUE);
# 
# summary.ame(p.ame.Y_tasty)
# 
# p.ame.Y_tasty$EZ    # estimated expectation of Z matrix, so:
# 
# lambda <- exp(p.ame.Y_tasty$EZ)   # is the matrix containing the estimates of the lambda poiss parameters


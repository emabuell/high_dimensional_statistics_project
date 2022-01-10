###########################
#Funzioni utili
###########################


threshold <- function(Y,thresh){  
   
   #34*(rmean[i]+rmean[j])/(ncol(Y))
   
   #rmean <- rowMeans(Y)
   
   for(i in 1:nrow(Y)){
      
      for(j in 1:ncol(Y)){
         
         if(i!=j & Y[i,j]>thresh){
            
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
         
         if(i!=j & Y[i,j]<thresh ){
            
            Y[i,j] = 0;
            
            
         }
         
      }
      
   }
   
   
   return(Y);
   
}



simulate_stat <- function(ame_model,psi,n_samples,...){  # simulate statistics
   
   statistics = list();
   
   
   for(i in 1:n_samples){
      
      Y_sample <- simY_bin(EZ = ame_model$EZ,rho = 1);
      statistics[[i]] <- psi(Y_sample,...);
      
   }
   
   return(statistics);
   
};


simulate_stat_ergm <- function(ergm_model,psi,n_samples,...){  # simulate statistics
   
   statistics = list();
   
   
   for(i in 1:n_samples){
      
      Y_sample <- graph_from_adjacency_matrix(
         as.matrix(simulate(ergm_model,nsim = 1,output = 'network')),
         mode = 'undirected'
         );
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



simulate_communities_ergm <- function(ergm_model,community_detection_algo,nsim,...){
   
   
   Y_sample <- as.matrix(simulate(ergm_model,nsim = 1,output = 'network'));
   
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
      
      Y_sample <- as.matrix(simulate(ergm_model,nsim = 1,output = 'network'));
      
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
# ERGM Giallo zafferano
###########################


#install.packages('ergm')
library(ergm)
library(amen)
#install.packages('statnet')
library(statnet)
Y_giallo = read.csv('sociogiallo.csv')
# rownames(Y) <- Y[,1];
# Y = Y[,-1];
rownames(Y_giallo) <- colnames(Y_giallo);
#colnames(Y) <- rownames(Y)
#View(Y)
sum_giallo <- sum(Y_giallo);
Y_giallo <- Y_giallo/sum_giallo;



View(Y_giallo)

plot.sociomatrix(as.matrix(Y_giallo),diaglab = FALSE,cex.lab = 1,srt=90,
                 drawlab = F);

Y_net <- as.network(as.matrix(Y_giallo),directed = T);

p_Y <- ergm(Y_net ~ edges+triangles ,reference = ~ Bernoulli,symmetric = TRUE);

library(ergm.count)

#p_Y <- ergm(Y_net ~ edges,reference = ~ Poisson,symmetric = TRUE);

p_Y_giallo_ergm <- ergm(Y_net ~ edges+triangles+ sender + mutual ,symmetric = TRUE,
            estimate = 'CD',reference = ~ Bernoulli,
            control = control.ergm(CD.maxit = 100)
            ); # p1 model


comm_giallo_ergm <- simulate_communities_ergm(p_Y_giallo_ergm,cluster_fast_greedy,200);


interaction_comm_matrix_giallo_ergm <- 0;

for(sample in comm_giallo_ergm[[1]]){   # computing community interaction matrix
   
   for(r in 1:nrow(sample)){
      
      interaction_comm_matrix_giallo_ergm <- interaction_comm_matrix_giallo_ergm+
         outer(sample[,r],sample[,r])/length(comm_giallo_ergm[[1]]);
   }
   
   
}

colnames(interaction_comm_matrix_giallo_ergm) <- colnames(Y_giallo);
rownames(interaction_comm_matrix_giallo_ergm) <- colnames(Y_giallo);


dev.new()
par(mfrow = c(1,1));
netplot(Y = threshold(interaction_comm_matrix_giallo_ergm,0.5),plot.iso=FALSE,plotnames=TRUE,
        ncol="Blue",main = 'Giallo community interaction graph ergm')



# giallo stats gof

d <- simulate_stat_ergm(p_Y_giallo_ergm,degree,50);
mdeg <- 0

for(el in d){
   
   mdeg <- mdeg+el/length(d);
   
}



#############################################
# ERGM tasty
#############################################


Y_tasty = read.csv('sociotasty.csv')
# rownames(Y) <- Y[,1];
# Y = Y[,-1];
rownames(Y_tasty) <- colnames(Y_tasty);
#colnames(Y) <- rownames(Y)
#View(Y)
Y_tasty <- Y_tasty/sum(Y_tasty);

Y_net <- as.network(as.matrix(Y_tasty),directed = T);


p_Y_tasty_ergm <- ergm(Y_net ~ edges+triangles+ sender + mutual ,symmetric = TRUE,
                        estimate = 'CD',reference = ~ Bernoulli,
                        control = control.ergm(CD.maxit = 100)
); # p1 model


comm_tasty_ergm <- simulate_communities_ergm(p_Y_tasty_ergm,cluster_fast_greedy,200);


interaction_comm_matrix_tasty_ergm <- 0;

for(sample in comm_tasty_ergm[[1]]){   # computing community interaction matrix
   
   for(r in 1:nrow(sample)){
      
      interaction_comm_matrix_tasty_ergm <- interaction_comm_matrix_tasty_ergm+
         outer(sample[,r],sample[,r])/length(comm_tasty_ergm[[1]]);
   }
   
   
}

colnames(interaction_comm_matrix_tasty_ergm) <- colnames(Y_tasty);
rownames(interaction_comm_matrix_tasty_ergm) <- colnames(Y_tasty);


dev.new()
par(mfrow = c(1,1));
netplot(Y = threshold(interaction_comm_matrix_tasty_ergm,0.5),plot.iso=FALSE,plotnames=TRUE,
        ncol="Blue",main = 'Tasty community interaction graph ergm')


v <- as.matrix(simulate(p_Y_tasty_ergm,nsim = 1,type = 'network'))/100;
for(i in 2:100){
   
   v <- v+ as.matrix(simulate(p_Y_tasty_ergm,nsim = 1,type = 'network'))/100;
}

netplot(Y = v,plot.iso=FALSE,plotnames=TRUE,
        ncol="Blue",main = 'Tasty community interaction graph ergm')



# #install.packages('amen')
# library(amen)
# 
# p.ame.Y <- ame(as.matrix(Y,dimnames = list( colnames(Y),colnames(Y)  )),
#                family = 'ord',symmetric = TRUE);
# summary(p.ame.Y)
# 
# p.ame.Y$EZ    # estimated expectation of Z matrix, so:
# 
# lambda <- exp(p.ame.Y$EZ)   # is the matrix containing the estimates of the lambda poiss parameters


#plot( apply(Y,1,mean,na.rm=TRUE) , type = 'l' ) 




#################################################
# Choosing the threshold for binarization
#################################################
# 
# Y_base = read.csv('sociogiallo.csv');
# rownames(Y) <- colnames(Y);
# edge_density_giallo = NULL;
# 
# for(thresh in 1:max(Y_base)){
# 
#    Y = Y_base;
#    # rownames(Y) <- Y[,1];
#    # Y = Y[,-1];
#    #rownames(Y) <- colnames(Y);
#    #colnames(Y) <- rownames(Y)
#    #View(Y)
#    
#    Y = threshold(Y,thresh)
#    
#    
#    edge_density_giallo = c(edge_density_giallo,mean(Y==1));
#    
#    
# }
# 
# 
# plot(edge_density_giallo)
# 
# edge_density_giallo[12]
# 
# Y_giallo = threshold(Y_base,12);
# 
# rownames(Y_giallo) <- colnames(Y_giallo);
# 
# 
# 
# 
# 
# 
# Y_base = read.csv('sociotasty.csv');
# rownames(Y) <- colnames(Y);
# edge_density_tasty = NULL;
# 
# for(thresh in 1:max(Y_base)){
#    
#    Y = Y_base;
#    # rownames(Y) <- Y[,1];
#    # Y = Y[,-1];
#    #rownames(Y) <- colnames(Y);
#    #colnames(Y) <- rownames(Y)
#    #View(Y)
#    
#    Y = threshold(Y,thresh)
#    
#    
#    edge_density_tasty = c(edge_density_tasty,mean(Y==1));
#    
#    
# }
# 
# 
# plot(edge_density_tasty)
# 
# edge_density_tasty[12]
# 
# Y_tasty = threshold(Y_base,12);
# 
# rownames(Y_tasty) <- colnames(Y_tasty);



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
comm_giallo <- simulate_communities(p.ame.Y_giallo,cluster_louvain,100);


print(comm_giallo[[2]])   # mean community matrix

interaction_comm_matrix_giallo <- 0

for(sample in comm_giallo[[1]]){   # computing community interaction matrix
   
   for(r in 1:nrow(sample)){
      
      interaction_comm_matrix_giallo <- interaction_comm_matrix_giallo+
         outer(sample[,r],sample[,r])/length(comm_giallo[[1]]);
   }
   
   
}

colnames(interaction_comm_matrix_giallo) <- colnames(Y_giallo);
rownames(interaction_comm_matrix_giallo) <- colnames(Y_giallo);


plot.sociomatrix(threshold_soft(as.matrix(interaction_comm_matrix_giallo),0.3),diaglab = FALSE,cex.lab = 1,srt=90,
                 drawlab = F);



View(as.matrix(interaction_comm_matrix_giallo))    # element (i,j) estimates the probability of i and j of
                                 # being in the same community according to ame



dev.new()
par(mfrow = c(1,1));
netplot(Y = threshold(interaction_comm_matrix_giallo,0.5),plot.iso=FALSE,plotnames=TRUE,
        ncol="Blue",main = 'Giallo community interaction graph',seed  = 19)


# gof stats amen

d <- simulate_stat(p.ame.Y_giallo,rowSums,n_samples = 200,na.rm = TRUE)
mdeg <- 0

for(el in d){
   
   mdeg <- mdeg+el/length(d);
   
}

hist(mdeg,xlab = 'degree',freq = F,
     main = 'simulated mean degree distribution giallo\n amen model')
abline(v = mean(rowSums(threshold(Y_giallo,0),na.rm = T)),col = 'red')


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


################################
# Latent factors interpretation
################################


dev.new()
par(mfrow= c(1,1))

plot(p.ame.Y_tasty$U,xlab = 'U1',ylab = 'U2')
text(x = p.ame.Y_tasty$U[,1],y = p.ame.Y_tasty$U[,2] ,labels=colnames(Y_tasty), cex=0.9, font=2)

plot(p.ame.Y_giallo$U,xlab = 'U1',ylab = 'U2',col = 'red')
text(x = p.ame.Y_giallo$U[,1],y = p.ame.Y_giallo$U[,2] ,labels=colnames(Y_giallo), cex=0.9, font=2)

U_tasty <- p.ame.Y_tasty$U;
U_giallo <- p.ame.Y_giallo$U;

delta <- 0.01;
delta_theta <- delta*pi;
rtheta <- matrix(data = c(cos(delta_theta),sin(delta_theta),
                          -sin(delta_theta),cos(delta_theta)),
                 nrow = 2
                 );



errors <- c();

U_rot_giallo <- U_giallo

for(i in 1:(2/delta)){
   
   U_rot_giallo <- U_rot_giallo%*%rtheta;
   
   errors <- c(errors,mean((U_rot_giallo-U_tasty)^2));
};

plot(x = delta*1:(2/delta),y = errors,xlab = 'theta',type = 'l',main = 'least square angle');
abline(v = delta*which.min(errors),col = 'red')
text(x = delta*which.min(errors),y = 0.03,labels = delta*which.min(errors))


delta_theta <- delta*which.min(errors);
U_rot_giallo <- U_giallo%*%matrix(data = c(cos(delta_theta),sin(delta_theta),
                                           -sin(delta_theta),cos(delta_theta)),
                                  nrow = 2);

# rand_pert <- matrix(rnorm(length(U_tasty),sd = 0.01),ncol = 2);
# U_tasty <- U_tasty + rand_pert
# U_giallo <- U_giallo + rand_pert;

#par(mar = c(5, 4, 1.4, 10))
#par(mar = 4.1*c(1,1,1,1))
#par(mar = c(0, 0, 0, 0))
plot(U_rot_giallo,xlab = 'U_1',ylab = 'U_2',pch = '',main = 'latent effects for amen model')
text(x = U_rot_giallo[,1],y = U_rot_giallo[,2] ,labels=colnames(Y_giallo), cex=0.8, font=1)
points(U_tasty,xlab = 'U_1',ylab = 'U_2',col = 'red',pch = '')
text(x = U_tasty[,1],y = U_tasty[,2] ,labels=colnames(Y_tasty), cex=0.8, font=1,col = 'red')
legend("bottomright", legend=c("giallo", "tasty"),
       col=c("black", "red"),lty = 1, cex=0.6)


###############################

dev.new()
par(mfrow= c(1,1))

#circplot(Y = as.matrix(Y_giallo),plotnames=TRUE)
netplot(Y = as.matrix(Y_tasty),plot.iso=FALSE,plotnames=TRUE)


dev.new()
par(mfrow= c(1,1))

netplot(Y = as.matrix(Y_tasty),plot.iso=FALSE,plotnames=TRUE)





comm_tasty <- simulate_communities(p.ame.Y_tasty,cluster_louvain,100);

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


dev.new()
par(mfrow = c(1,1));
netplot(Y = threshold(interaction_comm_matrix_tasty,0.5),plot.iso=FALSE,plotnames=TRUE,
        ncol="Blue",main = 'Tasty community interaction graph')


# boxplot ingredients degree simulated

library(data.table)

d_tasty <- simulate_stat(p.ame.Y_tasty,rowSums,n_samples = 300,na.rm = TRUE)
d_giallo <- simulate_stat(p.ame.Y_giallo,rowSums,n_samples = 300,na.rm = TRUE)

d_data_frame_tasty <- list();
d_data_frame_giallo <- list();

for( i in 1:length(d_giallo)){
   
   d_data_frame_tasty[[i]] <- as.data.frame(t(d_tasty[[i]]));
   d_data_frame_giallo[[i]] <- as.data.frame(t(d_giallo[[i]]));
   
}

d_data_frame_giallo <- rbindlist(d_data_frame_giallo);
d_data_frame_tasty <- rbindlist(d_data_frame_tasty);

d_data_frame_giallo$label <- 'giallo';
d_data_frame_tasty$label <- 'tasty';

d <- rbind.data.frame(d_data_frame_giallo,d_data_frame_tasty);
library(reshape2)
d.m  <- melt(d, id.var = "label");


require(ggplot2)
ggplot(data = d.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=label))+
   ggtitle("Ame simulated degrees")+ 
   #theme(plot.title = element_text(hjust = 1))+
   theme_bw(base_size = 12) +
   theme(axis.text.x=element_text(angle=45,hjust=1))+
   labs(y="degree", x = "ingredient")
   




# ggplot(data = d.m, aes(x=variable, y=value)) + geom_point(aes(colour=label),alpha = 0.1)+
#    ggtitle("Ame simulated degrees")+ 
#    theme(plot.title = element_text(hjust = 1))+
#    theme_bw(base_size = 12) +
#    theme(axis.text.x=element_text(angle=45,hjust=1))+
#    labs(y="degree", x = "ingredient")

# fused lasso on ingredients


#delta <- as.data.frame(-d_data_frame_giallo[,-c('label')]+
                                             #d_data_frame_tasty[,-c('label')]);
delta <- as.data.frame(t(colSums(Y_tasty>0,na.rm = T)-colSums(Y_giallo>0,na.rm = T)));
#delta <- t(delta);
# plot(colMeans(delta),xaxt = 'n' ,xlab = 'ingredient',
#      ylab = 'degree tasty-degree giallo');

plot(x = 1:ncol(delta),y = delta,xaxt = 'n' ,xlab = '',
           ylab = 'degree tasty - degree giallo',
            main = 'difference in degree');

#lines(fused_delta$mu,col  ='red')
points(1:length(fit.sol1),fit.sol1,col = 'red')
axis(1, at=1:length(colnames(delta)), labels=colnames(delta), cex.axis=0.8,srt = 45,
     las = 2)


library(flsa);  

fused_delta <- flsa(colMeans(delta));
fit.sol1 <-flsaGetSolution(fused_delta,lambda1 = 8,lambda2 = 0)


#install.packages('mvmeta')
#library(mvmeta)


#boxplot(rbindlist(d_data_frame))
#abline(as.data.frame(t(rowSums(threshold(Y_tasty,0),na.rm = T))),col = 'red')


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

mu <- 0;

for(i in 1:200){
   
   Y_sample <- simY_bin(EZ = p.ame.Y_giallo$EZ,rho = 1);
   mu <- mu + Y_sample/200;
   
};

#mu <- mu/sum(mu,na.rm = T)

hist(as.matrix ( (Y_giallo-mu)/sqrt(mu*(1-mu)) ))



########################
#simulation for the gof
########################

####################
#sd rowmean
####################

sd_ame_giallo <- simulate_stat(p.ame.Y_giallo,rowMeans,200,na.rm  = T);
sd_ame_tasty <- simulate_stat(p.ame.Y_tasty,rowMeans,200,na.rm  = T);

sd_ame_giallo <- Map(sd,sd_ame_giallo);
sd_ame_tasty <- Map(sd,sd_ame_tasty);

sd_ame_giallo_vec <- c()
sd_ame_tasty_vec <- c()
for(i in 1:length(sd_ame_giallo)){
   sd_ame_giallo_vec <- c(sd_ame_giallo_vec,sd_ame_giallo[[i]]);
   sd_ame_tasty_vec <- c(sd_ame_tasty_vec,sd_ame_tasty[[i]]);
}


hist(sd_ame_giallo_vec,xlab = 'sd degrees',freq = F,
     main = 'simulated sd of degrees giallo\n amen model')
abline(v = sd(rowMeans(Y_giallo>0,na.rm = T)),col = 'red')


hist(sd_ame_tasty_vec,xlab = 'sd degrees',freq = F,
     main = 'simulated sd of degrees tasty\n amen model')
abline(v = sd(rowMeans(Y_tasty>0,na.rm = T)),col = 'red')



####################
# degree 
####################

nsim <- 200;
d_giallo <- simulate_stat(p.ame.Y_giallo,rowSums,n_samples = nsim,na.rm = TRUE);
d_tasty <- simulate_stat(p.ame.Y_tasty,rowSums,n_samples = nsim,na.rm = TRUE);
mdeg_giallo <- 0;
mdeg_tasty <- 0;

for(i in 1:nsim){
   
   mdeg_giallo <- mdeg_giallo+d_giallo[[i]]/length(d_giallo);
   mdeg_tasty <- mdeg_tasty+d_tasty[[i]]/length(d_tasty);
}

hist(mdeg_giallo,xlab = 'degree',freq = F,
     main = 'simulated mean degree distribution giallo\n amen model')
abline(v = mean(rowSums(Y_giallo>0,na.rm = T)),col = 'red')


hist(mdeg_tasty,xlab = 'degree',freq = F,
     main = 'simulated mean degree distribution tasty\n amen model')
abline(v = mean(rowSums(Y_tasty>0,na.rm = T)),col = 'red')

################
# Grand mean
################

nsim <- 500;
gm_giallo <- simulate_stat(p.ame.Y_giallo,mean,n_samples = nsim,na.rm = TRUE);
gm_tasty <- simulate_stat(p.ame.Y_tasty,mean,n_samples = nsim,na.rm = TRUE);



gm_giallo_vec <- c();
gm_tasty_vec <- c();

for(i in 1:nsim){
   
   gm_giallo_vec <- c(gm_giallo_vec,gm_giallo[[i]]);
   gm_tasty_vec <- c(gm_tasty_vec,gm_tasty[[i]]);
}

hist(gm_giallo_vec,xlab = 'grand mean',freq = F,
     main = 'simulated grand mean distribution giallo\n amen model')
abline(v = mean(Y_giallo>0,na.rm = T),col = 'red')


hist(gm_tasty_vec,xlab = 'grand mean',freq = F,
     main = 'simulated grand mean distribution tasty\n amen model')
abline(v = mean(Y_tasty>0,na.rm = T),col = 'red')







netplot(Y = as.matrix(U_giallo%*%t(U_giallo) > 100*mean(U_giallo%*%t(U_giallo))),
        plot.iso=FALSE,plotnames=TRUE,
        main = 'positive latent effect interaction giallo',ncol = 'blue')



netplot(Y = as.matrix(U_tasty%*%t(U_tasty) > 80*mean(U_tasty%*%t(U_tasty))),
        plot.iso=FALSE,plotnames=TRUE,
        main = 'positive latent effect interaction tasty',ncol = 'blue')


delta_U_dot <- U_tasty%*%t(U_tasty)-U_giallo%*%t(U_giallo);

plot.sociomatrix(abs(as.matrix(delta_U_dot))>2*mean(abs(delta_U_dot)),diaglab = FALSE,cex.lab = 0.5,srt=90,
                 drawlab = T);



#####################
# Fused for latent vars
#####################


library(flsa);  

fused_delta <- flsa(delta_U_dot);
fit.sol1 <-flsaGetSolution(fused_delta,lambda1 = 0.08,lambda2 = 0.001)


plot(x = 1:ncol(delta_U_dot)^2,y = delta_U_dot,xaxt = 'n' ,xlab = '',
     ylab = 'latent dots tasty - giallo',
     main = 'difference in latent interaction');

#lines(fused_delta$mu,col  ='red')
points(1:length(fit.sol1),fit.sol1,col = 'red',cex = 0.6)
axis(1, at=1:length(colnames(delta_U_dot))^2, labels=outer(colnames(Y_giallo),colnames(Y_giallo),FUN = paste), cex.axis=0.6,srt = 45,
     las = 2);


outer(colnames(Y_giallo),colnames(Y_giallo),FUN = paste)[fit.sol1>0]

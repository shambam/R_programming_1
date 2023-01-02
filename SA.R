install.packages("progress")
library(progress)

ycc <- read.delim("Spellman_Yeast_Cell_Cycle.tsv",row.names=1,header=T,sep="\t")
ycc <- as.matrix(ycc)

scale.01 <- function(v){
    sc.01 <- (v-min(v))/(max(v)-min(v))
    sc.01
}

ycc.01 <- t(apply(ycc,1,scale.01))


calc.V.tot <- function(m,clus){

    clus.V <- NULL

    for(i in 1:max(clus)){
        #print(i)
        clus.d <- m[which(clus==i),]
        #clus.V <- c(clus.V,sqrt(sum(dist(clus.d)^2)))
        clus.V <- c(clus.V,sum(dist(clus.d)))
    }
    
    sum(clus.V)/max(clus)
}


calc.V.tot(ycc.01,kmeans(ycc.01,12)$cluster)

calc.exp <- function(V.new,V.old,T){
    exp(-((V.new-V.old)/T))
}

Temp <-0.3
cool <- 0.9999
K <- 9
clusters <- sample(1:K,nrow(ycc.01),replace = T)

clusters.o <- clusters

V.old <- calc.V.tot(ycc.01,clusters)
V.old
V.start <- V.old
choice <- 1:K
Iter <- 500000
#for(i in 1:nrow(flips)){

pb <- progress_bar$new(total = Iter)

for(i in 1:Iter){   

    clusters.new <- clusters

    row.id <- sample(1:nrow(ycc),1)
    
    clusters.new[row.id] <- sample(choice,1)

    V.new <- calc.V.tot(ycc.01,clusters.new)

    if(V.new==V.old){
        #Temp <- Temp*cool
        next
    }

    if(V.new < V.old){
        clusters <-  clusters.new
        V.old <- V.new
        
    }else{
       
       if(calc.exp(V.new,V.old,Temp) > runif(1)){
        clusters <- clusters.new
        V.old <- V.new
       }else{
            next
       }
    }
    
    #pb$tick()
    Sys.sleep(1 / Iter)
    #print(V.old)
    {cat("\r",V.old)}

    Temp <- Temp*cool
}
calc.V.tot(ycc.01,kmeans(ycc.01,9)$cluster)

Temp
V.old

par(mfrow=c(3,3))

for(i in 1:9){
  
  ycc.c <- ycc.01[which(clusters==i),]
  plot(ycc.c[1,],ty="l",ylim=range(ycc.c))
  apply(ycc.c,1,lines)
  
}




for(j in 1:100000){

    Temp <- Temp*cool
    print(Temp)

}

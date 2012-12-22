\name{data_simu_TcGSA}
\docType{data}
\alias{data_simu_TcGSA}
\alias{expr_1grp}
\alias{expr_2grp}
\alias{design}
\alias{gmt_sim}

\title{
Simulated Data for TcGSA
}

\description{
Simulated data for 5 genesets of 50 genes each.  Gene expression is simulated at 5 time points for 10 patients.
}

\usage{
data(data_simu_TcGSA)
}


\value{
  \item{expr_1grp}{
See Details.
}
  \item{expr_2grp}{
See Details.
}
  \item{design}{
a data frame with 5 variables:\itemize{
  \item \code{Patient_ID}: a factor that contains the patient ID.

  \item \code{TimePoint}: a numeric vector or a factor that contains the time points at which gene expression was measured.

\item \code{sample_name}: a character vector with the names of the sample (corresponding to the names of the columns of \code{expr_1grp} and of \code{expr_2grp}).

\item \code{group.var}: a factor that indicates to which of the 2 treatment groups each sample belongs to.

\item \code{Group_paired_ID}
a random paired identifier for paired couples (one in each of the 2 treatment groups) of patients.
}
}
  \item{gmt_sim}{
a \bold{gmt} object containing the gene sets definition.  See \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{GMT definition on www.broadinstitute.org}.
}
}

\details{

In \code{expr_1grp} all patients belong to the same unique treatment group.  The first 2 gene sets are simulated under the null hypothesis.  The gene sets 3 and 4 are simulated under the alternative hypothesis that there is a significant homogeneous time trend within the gene set.  The gene set 5 is simulated under the alternative hypothesis that there are significant heterogeneous time trends within the gene set.
\cr

In \code{expr_2grp} all patients belong to 2 treatment groups.  The 5 first patients belong to the treatment group '\code{T}', The 5 other patients belong to the treatment group '\code{C}'.  The first 2 gene sets are simulated under the null hypothesis that there is no difference in the time trend between the 2 treatment goups.  The gene sets 3 and 4 are simulated under the alternative hypothesis that there are significantly different homogeneous time trends within the gene set between the 2 treatment groups.  The gene set 5 are simulated under the alternative hypothesis that there are significantly different heterogeneous time trends between the 2 treatment groups within the gene set.

}

\source{
This is simulated data.
}

\author{
Boris P. Hejblum
}

\seealso{
\code{\link{TcGSA.LR}}
}

\examples{

data(data_simu_TcGSA)
summary(expr_1grp)
summary(design)
gmt_sim

\dontrun{
simu_data<-function(H=c("H0","H1_homogeneous","H1_heterogeneous"), G=50,
P=20, I=8, N.grp=1){
  N=P*G*I
  gene <- NULL
  for (g in 1:G){
    gene <- c(gene,rep(as.character(g),P*I))
  }
  gene <- as.factor(gene)
  
  patient <- NULL
  for (p in 1:P){
    patient <- c(patient,rep(as.character(p),I))
  }
  patient <- rep(patient,G)
  patient <- as.factor(patient)
  
  t1 <- NULL
  t1 <- rep(0:(I-1),G*P)
  t1 <- t1/10
  t2 <- t1^2
  t3 <- t1^3
  
  group <- (which(rmultinom(G, size=1, prob=c(1,1,1))==1, arr.ind=T)[,1]-2)
  
  mu <- runif(n=1, min=0, max=7)
  beta <- runif(n=G, min=0, max=3)
  random_G <- rnorm(n=G,mean=0,sd=0.05)
  random_P <- rnorm(n=P,mean=0,sd=0.1)
  sigma <- runif(n=1, min=0, max=0.5)
  
  fixed <- numeric(3)
  fixed[1] <- rnorm(n=1, mean=0, sd=2)
  fixed[2] <- rnorm(n=1, mean=0, sd=2)
  fixed[3] <- rnorm(n=1, mean=1, sd=2)
  random_T1 <- group*fixed[1] + rnorm(n=G, mean=0, sd=0.2)
  random_T2 <- group*fixed[2] + rnorm(n=G, mean=0, sd=0.2)
  random_T3 <- group*fixed[3] + rnorm(n=G, mean=0, sd=0.2)
  
  y <- numeric(N)
  
  
  
  
  
  if(N.grp==1){
    if(H=="H0"){
      for (g in 1:G){
        for (p in 1:P){
          for (i in 1:I){
            obs <- (g-1)*P*I+(p-1)*I+i
            y[obs] <- (mu + beta[g] + random_G[g] + random_P[p]
                       + rnorm(n=1,mean=0,sd=sigma))
          }
        }
      }
    }else if(H=="H1_homogeneous"){
      for (g in 1:G){
        for (p in 1:P){
          for (i in 1:I){
            obs <- (g-1)*P*I+(p-1)*I+i
            y[obs] <- (mu + beta[g] + random_G[g] + random_P[p] 
                       + rnorm(n=1,mean=0,sd=sigma) 
                       + fixed[1]*t1[obs] + fixed[2]*t2[obs]+ fixed[3]*t3[obs])
          }
        }
      }
    }else if(H=="H1_heterogeneous"){
      beta <- beta + (-group+1)
      for (g in 1:G){
        for (p in 1:P){
          for (i in 1:I){
            obs <- (g-1)*P*I+(p-1)*I+i
            y[obs] <- (mu + beta[g] + random_G[g] + random_P[p] 
                       + rnorm(n=1,mean=0,sd=sigma)
                       + (fixed[1]+random_T1[g])*t1[obs]
                       + (fixed[2]+random_T2[g])*t2[obs]
                       + (fixed[3]+random_T3[g])*t3[obs])
          }
        }
      }
    }
  }else if(N.grp>1){
    fixed_grp <- list()
    random_grp <- list()
    grp_tt <- as.factor(patient)
    levels(grp_tt) <- as.character(rep(1:N.grp, length(levels(grp_tt))/2))
    grp_tt <- as.numeric(as.character(grp_tt))
    for(i in 1:N.grp){
      fixed_grp[[i]] <- numeric(3)
      fixed_grp[[i]][1] <- rnorm(n=1, mean=0, sd=2)
      fixed_grp[[i]][2] <- rnorm(n=1, mean=0, sd=2)
      fixed_grp[[i]][3] <- rnorm(n=1, mean=0, sd=2)
      random_grp[[i]] <- list()
      random_grp[[i]][[1]] <- group*fixed[1] + rnorm(n=G, mean=0, sd=0.2)
      random_grp[[i]][[2]] <- group*fixed[2] + rnorm(n=G, mean=0, sd=0.2)
      random_grp[[i]][[3]] <- group*fixed[3] + rnorm(n=G, mean=0, sd=0.2)
    } 
    if(H=="H0"){
      for (g in 1:G){
        for (p in 1:P){
          for (i in 1:I){
            obs <- (g-1)*P*I+(p-1)*I+i
            y[obs] <- (mu + beta[g] + random_G[g] + random_P[p]
                       + rnorm(n=1,mean=0,sd=sigma)
                       + (fixed[1]+random_T1[g])*t1[obs]
                       + (fixed[2]+random_T2[g])*t2[obs]
                       + (fixed[3]+random_T3[g])*t3[obs])
          }
        }
      }
    }else if(H=="H1_homogeneous"){
      for (g in 1:G){
        for (p in 1:P){
          for (i in 1:I){
            obs <- (g-1)*P*I+(p-1)*I+i
            y[obs] <- (mu + beta[g] + random_G[g] + random_P[p] 
                       + rnorm(n=1,mean=0,sd=sigma)
                       + (fixed_grp[[grp_tt[obs]]][1])*t1[obs]
                       + (fixed_grp[[grp_tt[obs]]][2])*t2[obs]
                       + (fixed_grp[[grp_tt[obs]]][3])*t3[obs])
          }
        }
      }
    }else if(H=="H1_heterogeneous"){
      beta <- beta + (-group+1)
      for (g in 1:G){
        for (p in 1:P){
          for (i in 1:I){
            obs <- (g-1)*P*I+(p-1)*I+i
            y[obs] <- (mu + beta[g] + random_G[g] + random_P[p] 
                       + rnorm(n=1,mean=0,sd=sigma)
                       + (fixed_grp[[grp_tt[obs]]][1] 
                          + random_grp[[grp_tt[obs]]][[1]][g])*t1[obs]
                       + (fixed_grp[[grp_tt[obs]]][2] 
                          + random_grp[[grp_tt[obs]]][[2]][g])*t2[obs]
                       + (fixed_grp[[grp_tt[obs]]][3] 
                          + random_grp[[grp_tt[obs]]][[3]][g])*t3[obs])
          }
        }
      }
    }
  }else{
    cat("ERROR: wrong value for the 'N.grp' argument\n")
  }
  return(cbind.data.frame(y,gene,patient,t1,t2,t3))
}


expr_1grp <- NULL
for( i in 1:2){
  data_tomerge <- simu_data("H0", G=50, P=10, I=5, N.grp=1)
  wide <- dcast(data_tomerge, "gene ~ patient + t1", value.var="y")
  expr_1grp <- rbind.data.frame(expr_1grp, wide[,-1])
  print(dim(expr_1grp))
}
for( i in 1:2){
  data_tomerge <- simu_data("H1_homogeneous", G=50, P=10, I=5, N.grp=1)
  wide <- dcast(data_tomerge, "gene ~ patient + t1", value.var="y")
  rownames(wide) <- (nrow(wide)*(i-1)+1):(nrow(wide)*i)
  expr_1grp <- rbind.data.frame(expr_1grp, wide[,-1])
  print(dim(expr_1grp))
}
for( i in 1:1){
  data_tomerge <- simu_data("H1_heterogeneous", G=50, P=10, I=5, N.grp=1)
  wide <- dcast(data_tomerge, "gene ~ patient + t1", value.var="y")
  rownames(wide) <- (nrow(wide)*(i-1)+1):(nrow(wide)*i)
  expr_1grp <- rbind.data.frame(expr_1grp, wide[,-1])
  print(dim(expr_1grp))
}
rownames(expr_1grp) <- 1:250
colnames(expr_1grp) <- 1:50

expr_2grp <- NULL
for( i in 1:2){
  data_tomerge <- simu_data("H0", G=50, P=10, I=5, N.grp=2)
  wide <- dcast(data_tomerge, "gene ~ patient + t1", value.var="y")
  expr_2grp <- rbind.data.frame(expr_2grp, wide[,-1])
  print(dim(expr_2grp))
}
for( i in 1:2){
  data_tomerge <- simu_data("H1_homogeneous", G=50, P=10, I=5, N.grp=2)
  wide <- dcast(data_tomerge, "gene ~ patient + t1", value.var="y")
  rownames(wide) <- (nrow(wide)*(i-1)+1):(nrow(wide)*i)
  expr_2grp <- rbind.data.frame(expr_2grp, wide[,-1])
  print(dim(expr_2grp))
}
for( i in 1:1){
  data_tomerge <- simu_data("H1_heterogeneous", G=50, P=10, I=5, N.grp=2)
  wide <- dcast(data_tomerge, "gene ~ patient + t1", value.var="y")
  rownames(wide) <- (nrow(wide)*(i-1)+1):(nrow(wide)*i)
  expr_2grp <- rbind.data.frame(expr_2grp, wide[,-1])
  print(dim(expr_2grp))
}
rownames(expr_2grp) <- 1:250
colnames(expr_2grp) <- 1:50


Patient_ID <- rep(1:10, 5)[order(rep(1:10, 5))]
TimePoint <- rep(1:5, 10)
genesets <- rep(1:5, 50)[order(rep(1:5, 50))]
group.var <- as.factor(Patient_ID)
levels(group.var) <- as.character(rep(c("T","C"), length(levels(group.var))/2))
group.var <- as.factor(as.character(group.var))
Group_ID_paired <- c(rep(1:5, 5), rep(1:5, 5))[order(rep(1:5, 5*2))]
design <- cbind.data.frame("Patient_ID" = as.factor(paste("P", 
                                            Patient_ID, sep = "")), 
                           "TimePoint" = TimePoint, 
                           "Sample_name"=colnames(expr_1grp),
                           "group.var" = group.var, 
                           "Group_ID_paired" = as.factor(paste("cpl", 
                                                 Group_ID_paired, sep = ""))
                           )
                           
gmt_sim <- list()
for(gs in 1:5){
  gmt_sim[["genesets"]][[gs]] <- ((gs-1)*50+1):(gs*50)
}
gmt_sim[["geneset.names"]] <- paste("Gene set", 1:5)
gmt_sim[["geneset.descriptions"]] <- c(
    rep("This is a simulated gene set under H0", 3),
    rep("This is a simulated gene set under H1 (homogeneous)", 1), 
    rep("This is a simulated gene set under H1 (heterogeneous)", 1)
  )
class(gmt_sim) <- "GSA.genesets"                           


#save(expr_1grp, expr_2grp, design, gmt_sim, 
      file="data_simu_TcGSA.RData", compress="xz")
}
}

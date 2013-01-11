estim.TcGSA <-
function(tcgsa, expr, Patient_ID, TimePoint){
# require(splines)
  
  func <- tcgsa$func_form
  group.var <- tcgsa$group.var
  gmt <- tcgsa$GeneSets_gmt
  estim_expr <- list()
  ngs <- length(gmt$genesets)
  pat <- as.character(unique(Patient_ID))
  npat <- length(unique(Patient_ID))
  tim <- unique(TimePoint)
  ntim <- length(unique(TimePoint))
  measuredProbes <- rownames(expr)
  
  t1 <- tim/10
  
  if(is.null(group.var)){
    if(func=="linear"){
      for (gs in 1:ngs){
        if(!is.na(tcgsa$Estimations[["RanEf"]][[gs]])[1]){
          if(!is.null(tcgsa$Estimations[["RanEf"]][[gs]]$probe)){
            prob <- rownames(tcgsa$Estimations[["RanEf"]][[gs]]$probe)
          }else{
            prob <- intersect(gmt$genesets[[gs]], measuredProbes)
          }
          nprob <- length(prob)
          estim_expr[[gs]] <- array(dim=c(nprob, npat, ntim), dimnames=list(prob, pat, tim))
          
          if(nprob<2){
            i=1
            for (j in 1:npat){
              for (k in 1:ntim){
                estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                            + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                            + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID"][[1]][pat[j],]
                )
              }
            }
          }else{
            i=1
            for (j in 1:npat){
              for (k in 1:ntim){
                estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                              + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t1"] 
                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                )
              }
            }
            for (i in 2:nprob){
              for (j in 1:npat){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")] 
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t1"] 
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }
            }
          }
          cat(paste(gs,"/", ngs," gene set dynamics estimated ", "\n", sep=""))
        }
      }
    }
    else if(func=="cubic"){
      t2 <- t1^2
      t3 <- t1^3
      for (gs in 1:length(gmt$genesets)){
        if(!is.na(tcgsa$Estimations[["RanEf"]][[gs]])[1]){
          if(!is.null(tcgsa$Estimations[["RanEf"]][[gs]]$probe)){
            prob <- rownames(tcgsa$Estimations[["RanEf"]][[gs]]$probe)
          }else{
            prob <- intersect(gmt$genesets[[gs]], measuredProbes)
          }
          nprob <- length(prob)
          estim_expr[[gs]] <- array(dim=c(nprob, npat, ntim), dimnames=list(prob, pat, tim))
          
          if(nprob<2){
            i=1
            for (j in 1:npat){
              for (k in 1:ntim){
                estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                            + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                            + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID"][[1]][pat[j],]
                )
              }
            }
          }else{
            i=1
            for (j in 1:npat){
              for (k in 1:ntim){
                estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                              + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t1"] + t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t2"] + t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t3"] 
                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                )
              }
            }
            for (i in 2:nprob){
              for (j in 1:npat){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")]
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t1"] + t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t2"] + t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"t3"] 
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }
            }
          }
          cat(paste(gs,"/", ngs," gene set dynamics estimated ", "\n", sep=""))
        }
      }
    }
    else if(func=="splines"){
      # library(splines)
      nk = ceiling(length(unique(t1))/4)
      noeuds = quantile(t1, probs=(c(0:(nk+1))/(nk+1))[-c(1,(nk+1+1))])
      bsplines <- as.data.frame(bs(t1, knots = noeuds, degree = 3, Boundary.knots = range(t1)), intercept = FALSE)
      colnames(bsplines) <- paste("spline_t", colnames(bsplines) , sep="")
      bsplines <- bsplines*10
      for (gs in 1:length(gmt$genesets)){
        if(!is.na(tcgsa$Estimations[["RanEf"]][[gs]])[1]){
          if(!is.null(tcgsa$Estimations[["RanEf"]][[gs]]$probe)){
            prob <- rownames(tcgsa$Estimations[["RanEf"]][[gs]]$probe)
          }else{
            prob <- intersect(gmt$genesets[[gs]], measuredProbes)
          }
          nprob <- length(prob)
          estim_expr[[gs]] <- array(dim=c(nprob, npat, ntim), dimnames=list(prob, pat, tim))
          
          if(nprob<2){
            i=1
            for (j in 1:npat){
              for (k in 1:ntim){
                estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                            + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                            + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID"][[1]][pat[j],]
                )
              }
            }
          }
          else{
            i=1
            for (j in 1:npat){
              for (k in 1:ntim){
                estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                              + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                              + bsplines$spline_t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t5"]
                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                )
              }
            }
            for (i in 2:nprob){
              for (j in 1:npat){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")]
                                                                + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                + bsplines$spline_t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"spline_t5"]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }
            }
          }
          cat(paste(gs,"/", ngs," gene set dynamics estimated ", "\n", sep=""))
        }
      }
    }
  }
  else{
    
    patients_group <- acast(cbind.data.frame(Patient_ID, "Group"=group.var), margins="Group", formula=Patient_ID~., fun.aggregate=function(x){sum(as.numeric(x))/ntim})[,1]
    patients_group <- as.factor(patients_group[pat])
    levels(patients_group) <- levels(group.var)
    
    if(func=="linear"){
      for (gs in 1:ngs){
          if(!is.na(tcgsa$Estimations[["RanEf"]][[gs]])[1]){
            if(!is.null(tcgsa$Estimations[["RanEf"]][[gs]]$probe)){
              prob <- rownames(tcgsa$Estimations[["RanEf"]][[gs]]$probe)
            }else{
              prob <- intersect(gmt$genesets[[gs]], measuredProbes)
            }
          nprob <- length(prob)
          estim_expr[[gs]] <- array(dim=c(nprob, npat, ntim), dimnames=list(prob, pat, tim))
          
          if(nprob<2){
            i=1
            for (j in 1:npat){
              if (as.character(patients_group[j])==levels(group.var)[1]){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["probe"][[1]][pat[j],]
                  )
                }
              }else{
                for (k in 1:ntim){
                  thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                              + tcgsa$Estimations[["FixEf"]][[gs]][thisgr]
                                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t1", sep=":")]
                                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][pat[j],]
                  )
                }
              }
            }
          }else{
            i=1
            for (j in 1:npat){
              if (as.character(patients_group[j])==levels(group.var)[1]){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"] 
                                                                + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }else{
                for (k in 1:ntim){
                  thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                + tcgsa$Estimations[["FixEf"]][[gs]][thisgr]
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t1", sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"] 
                                                                + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }
            }
            for (i in 2:nprob){
              for (j in 1:npat){
                if (as.character(patients_group[j])==levels(group.var)[1]){
                  for (k in 1:ntim){
                    estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")] 
                                                                  + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"] 
                                                                  + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                    )
                  }
                }else{
                  for (k in 1:ntim){
                    thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                    estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")]
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][thisgr]
                                                                  + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"]
                                                                  + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t1", sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"] 
                                                                  + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                    )
                  }
                }
              }
            }
          }
          cat(paste(gs,"/", ngs," gene set dynamics estimated ", "\n", sep=""))
        }
      }
    }
    else if(func=="cubic"){
      t2 <- t1^2
      t3 <- t1^3
      for (gs in 1:length(gmt$genesets)){
          if(!is.na(tcgsa$Estimations[["RanEf"]][[gs]])[1]){
            if(!is.null(tcgsa$Estimations[["RanEf"]][[gs]]$probe)){
              prob <- rownames(tcgsa$Estimations[["RanEf"]][[gs]]$probe)
            }else{
              prob <- intersect(gmt$genesets[[gs]], measuredProbes)
            }
          nprob <- length(prob)
          estim_expr[[gs]] <- array(dim=c(nprob, npat, ntim), dimnames=list(prob, pat, tim))
          
          if(nprob<2){
            i=1
            for (j in 1:npat){
              if (as.character(patients_group[j])==levels(group.var)[1]){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID"][[1]][pat[j],]
                  )
                }
              }else{
                for (k in 1:ntim){
                  thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                              + tcgsa$Estimations[["FixEf"]][[gs]][thisgr] 
                                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                              + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t1", sep=":")] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t2", sep=":")] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t3", sep=":")]
                                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID"][[1]][pat[j],]
                  )
                }
              }
            }
          }else{
            i=1
            for (j in 1:npat){
              if (as.character(patients_group[j])==levels(group.var)[1]){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")] + t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t2", thisgr, sep=":")] + t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t3", thisgr, sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }else{
                for (k in 1:ntim){
                  thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                + tcgsa$Estimations[["FixEf"]][[gs]][thisgr] 
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t1", sep=":")] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t2", sep=":")] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t3", sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")] + t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t2", thisgr, sep=":")] + t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t3", thisgr, sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }
            }
            for (i in 2:nprob){
              for (j in 1:npat){
                if (as.character(patients_group[j])==levels(group.var)[1]){
                  for (k in 1:ntim){
                    estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]  
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")]
                                                                  + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                  + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")] + t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t2", thisgr, sep=":")] + t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t3", thisgr, sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                    )
                  }
                }else{
                  for (k in 1:ntim){
                    thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                    estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]  
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")]
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][thisgr] 
                                                                  + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t1"] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t2"] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["t3"]
                                                                  + t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t1", sep=":")] + t2[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t2", sep=":")] + t3[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "t3", sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                  + t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t1", thisgr, sep=":")] + t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t2", thisgr, sep=":")] + t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("t3", thisgr, sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                    )
                  }
                }
              }
            }
          }
          cat(paste(gs,"/", ngs," gene set dynamics estimated ", "\n", sep=""))
        }
      }
    }
    else if(func=="splines"){
      #  library(splines)
      nk = ceiling(length(unique(t1))/4)
      noeuds = quantile(t1, probs=(c(0:(nk+1))/(nk+1))[-c(1,(nk+1+1))])
      bsplines <- as.data.frame(bs(t1, knots = noeuds, degree = 3, Boundary.knots = range(t1)), intercept = FALSE)
      colnames(bsplines) <- paste("spline_t",colnames(bsplines) , sep="")
      bsplines <- bsplines*10
      for (gs in 1:length(gmt$genesets)){
          if(!is.na(tcgsa$Estimations[["RanEf"]][[gs]])[1]){
            if(!is.null(tcgsa$Estimations[["RanEf"]][[gs]]$probe)){
              prob <- rownames(tcgsa$Estimations[["RanEf"]][[gs]]$probe)
            }else{
              prob <- intersect(gmt$genesets[[gs]], measuredProbes)
            }
          nprob <- length(prob)
          estim_expr[[gs]] <- array(dim=c(nprob, npat, ntim), dimnames=list(prob, pat, tim))
            
          if(nprob<2){
            i=1
            for (j in 1:npat){
              if (as.character(patients_group[j])==levels(group.var)[1]){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                              + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                              + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                              + bsplines$spline_t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t1", thisgr, sep=":")] + bsplines$spline_t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t2", thisgr, sep=":")] + bsplines$spline_t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t3", thisgr, sep=":")] + bsplines$spline_t4[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t4", thisgr, sep=":")] + bsplines$spline_t5[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t5", thisgr, sep=":")]
                                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }else{
                for (k in 1:ntim){
                  thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                              + tcgsa$Estimations[["FixEf"]][[gs]][thisgr]
                                                                              + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                              + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t1", sep="")] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t2", sep="")] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t3", sep="")] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t4", sep="")] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t5", sep="")]
                                                                              + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID"][[1]][pat[j],]
                  )
                }
              }
            }
          }else{
            i=1
            for (j in 1:npat){
              if (as.character(patients_group[j])==levels(group.var)[1]){
                for (k in 1:ntim){
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                + bsplines$spline_t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t1", thisgr, sep=":")] + bsplines$spline_t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t2", thisgr, sep=":")] + bsplines$spline_t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t3", thisgr, sep=":")] + bsplines$spline_t4[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t4", thisgr, sep=":")] + bsplines$spline_t5[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t5", thisgr, sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }else{
                for (k in 1:ntim){
                  thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                  estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                + tcgsa$Estimations[["FixEf"]][[gs]][thisgr]
                                                                + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t1", sep="")] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t2", sep="")] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t3", sep="")] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t4", sep="")] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t5", sep="")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                + bsplines$spline_t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t1", thisgr, sep=":")] + bsplines$spline_t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t2", thisgr, sep=":")] + bsplines$spline_t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t3", thisgr, sep=":")] + bsplines$spline_t4[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t4", thisgr, sep=":")] + bsplines$spline_t5[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t5", thisgr, sep=":")]
                                                                + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                  )
                }
              }
            }
            for (i in 2:nprob){
              for (j in 1:npat){
                if (as.character(patients_group[j])==levels(group.var)[1]){
                  for (k in 1:ntim){
                    estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"] 
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")]
                                                                  + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"] 
                                                                  + bsplines$spline_t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t1", thisgr, sep=":")] + bsplines$spline_t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t2", thisgr, sep=":")] + bsplines$spline_t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t3", thisgr, sep=":")] + bsplines$spline_t4[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t4", thisgr, sep=":")] + bsplines$spline_t5[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t5", thisgr, sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                    )
                  }
                }else{
                  for (k in 1:ntim){
                    thisgr <- paste("Group", as.character(patients_group[j]), sep="")
                    estim_expr[[gs]][prob[i], pat[j], as.character(tim[k])] <- (tcgsa$Estimations[["FixEf"]][[gs]]["(Intercept)"]
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][paste("probe", prob[i], sep="")]
                                                                  + tcgsa$Estimations[["FixEf"]][[gs]][thisgr]
                                                                  + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t1"] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t2"] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t3"] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t4"] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]]["spline_t5"]
                                                                  + bsplines$spline_t1[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t1", sep=":")] + bsplines$spline_t2[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t2", sep=":")] + bsplines$spline_t3[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t3", sep=":")] + bsplines$spline_t4[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t4", sep=":")] + bsplines$spline_t5[k]*tcgsa$Estimations[["FixEf"]][[gs]][paste(thisgr, "spline_t5", sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,"(Intercept)"]
                                                                  + bsplines$spline_t1[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t1", thisgr, sep=":")] + bsplines$spline_t2[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t2", thisgr, sep=":")] + bsplines$spline_t3[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t3", thisgr, sep=":")] + bsplines$spline_t4[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t4", thisgr, sep=":")] + bsplines$spline_t5[k]*tcgsa$Estimations[["RanEf"]][[gs]]$probe[i,paste("spline_t5", thisgr, sep=":")]
                                                                  + tcgsa$Estimations[["RanEf"]][[gs]]["Patient_ID:probe"][[1]][paste(pat[j], prob[i], sep=":"),]
                    )
                  }
                }
              }
            }
          }
          cat(paste(gs,"/", ngs," gene set dynamics estimated ", "\n", sep=""))
        }
      }
    }
  }
  names(estim_expr) <- gmt$geneset.names
  return(estim_expr)
}

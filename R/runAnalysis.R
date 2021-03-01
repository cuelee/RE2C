### runAnalysis.R
### created by Cue Hyunkyu Lee

# necessary inputs <-
# 1. work_path(pwd)
# 2. inputfile(full path)
# 3. outputfile as path (out.txt as the default)
# 4. correlationfile as path (out.txt as the default)

# # 5. parallel computing (default = 1)
#
# ## runAnalysis is designed to use multiple cores
# parallelCapability <- function(no_cores){
#    ## library parallel should be existed in your R library folder
#    library(parallel)
#    print(paste("maximum number of cores: ",detectCores(),sep=""),quote = FALSE)
#    print(paste("no_cores: ",no_cores,sep=""),quote = FALSE)
#    if(no_cores > detectCores()) {stop("use small no_cores")} else{return(TRUE)}
#    stop("unknown error")
# }

parse_aLine = function(aLine){
  tryCatch({
    data <- strsplit(aLine,split="\\s+")[[1]]
    rsid <- data[1]
    
    ## accept argument NA
    readL <- data[-1][which_betas]
    nonNA <- !(is.na(readL) <- readL == "NA")
    
    if (all(nonNA)) {
      beta <- as.numeric(data[-1][which_betas])
      stders <- as.numeric(data[-1][which_stders])
      sigma <- diag(stders) %*% cor %*% diag(stders)
      nstudy = length(beta)
      
      out <- RE2C(beta,stders, sigma=sigma, stat2_cdf=as.numeric(RE2C_table[nstudy-1,])*RE2C_corr, FEt.lows=as.numeric(FEtlow_table[nstudy-1,]), NR = NR,tau2.zero.prob=tau2.zero.prob)
      LSout <- LS(beta,stders,cor)
      
      LSs<- LSout[1]
      LSse<- LSout[2]
      LSp<- LSout[3]
      
      RE2Cout<-out[[1]]
      stat1 <- RE2Cout[1]
      stat2 <- RE2Cout[2]
      RE2Cp_cond <- RE2Cout[3]
      RE2Cp <- RE2Cout[4]
      
      if (isdiag){
        Qout<-out[[2]]
        Q <- Qout[1]
        Qp <- pchisq(q=Q,df=(nstudy-1),ncp=0,lower.tail=F)
        
        Isq <- floor(max(100*(Q-(nstudy-1))/Q,0))
        Q = round(Q,2)
        return(list(out=T, text = paste(rsid, nstudy, LSs, LSse, LSp, Isq, Q, Qp, stat1, stat2, RE2Cp_cond, RE2Cp,sep = ' ')))
      } else {
        Q <- NA;  Qp <- NA;  Isq <- NA
        return(list(out=T, text = paste(rsid, nstudy, LSs, LSse, LSp, stat1, stat2, RE2Cp_cond, RE2Cp,sep = ' ')))
      }
      
    } else { 
      beta <- as.numeric(data[-1][which_betas][nonNA])
      stders <- as.numeric(data[-1][which_stders][nonNA])
      sigma <- diag(stders) %*% cor[nonNA,nonNA] %*% diag(stders)
      
      temp_ns = sum(nonNA)
      if(temp_ns < 2) { cat("Input should contain more than two studies\n : ", data)}
      
      temp_NR <- test_n(temp_ns)
      
      ## correction factors
      temp_RE2Cor <- RE2Cor.list[[min(temp_ns-1,6)]]
      temp_tau2prob_cor <- tau2prob_corlist[[min(temp_ns-1,6)]]
      
      ## run correction.function(cor)
      temp_correction.list <- correction.function(cor[nonNA,nonNA],temp_RE2Cor,temp_tau2prob_cor)
      temp_RE2C_corr <- temp_correction.list[[1]];   temp_tau2.zero.prob <- temp_correction.list[[2]]
      
      
      out <- RE2C(beta,stders, sigma=sigma, stat2_cdf=as.numeric(RE2C_table[temp_ns-1,])*temp_RE2C_corr, FEt.lows=as.numeric(FEtlow_table[temp_ns-1,]), NR = temp_NR,tau2.zero.prob=temp_tau2.zero.prob)
      LSout <- LS(beta,stders,cor[nonNA,nonNA])
      
      LSs<- LSout[1]
      LSse<- LSout[2]
      LSp<- LSout[3]
      
      RE2Cout<-out[[1]]
      stat1 <- RE2Cout[1]
      stat2 <- RE2Cout[2]
      RE2Cp_cond <- RE2Cout[3]
      RE2Cp <- RE2Cout[4]
      
      if (isdiag){
        Qout<-out[[2]]
        Q <- Qout[1]
        Qp <- pchisq(q=Q,df=(temp_ns-1),ncp=0,lower.tail=F)
        
        Isq <- floor(max(100*(Q-(temp_ns-1))/Q,0))
        Q = round(Q,2)
        return(list(out=T, text = paste0(rsid, temp_ns, LSs, LSse, LSp, Isq, Q, Qp, stat1, stat2, RE2Cp_cond, RE2Cp)))
      } else {
        Q <- NA;  Qp <- NA;  Isq <- NA
        return(list(out=T, text = paste0(rsid, temp_ns, LSs, LSse, LSp, stat1, stat2, RE2Cp_cond, RE2Cp)))
      }
      
    }
    return(list(out=F,m = 'ERROR'))
  }
  , error = function(e) {print(e); return(list(out = F, m = paste0(rsid,': ', e)))}
  , warning = function(w) {print(w); return(list(out = F, m = paste0(rsid,': ', w)))}
  )
}

## load compiler
require(compiler,warn.conflicts=TRUE,quietly=TRUE)
require(parallel,warn.conflicts=TRUE,quietly=TRUE)

## getting RE2C.sh input arguments
args <- commandArgs(trailingOnly = T)
if(args[1]!="NA") {work_path <- args[1]} else { stop ("RE2C.sh provided no work_path")}
if(args[2]!="NA") {input_file <- args[2]} else { stop ("RE2C.sh provided no inputFile_path")}
if(args[3]!="NA") {output_file <- paste(args[3],".txt",sep="");log_file <- paste(args[3],".log",sep="")
} else {stop("RE2C.sh provided no outputFile_path")}
if(args[4]!="NA") {cor_file<- args[4]} else {cor_file <- NA }
if(args[5]!="NA") {log_file<- args[5]} else {stop("Unknown problem occurs")}
if(args[6]!=1) {ncores <- min(detectCores(), as.numeric(args[6]))} else {ncores <- 1}

# work_path='/Users/cuelee/Dropbox/github/RE2C'
# input_file='/Users/cuelee/Dropbox/github/RE2C/example/example_input.txt'
# output_file='/Users/cuelee/Dropbox/github/RE2C/out'
# cor_file='/Users/cuelee/Dropbox/github/RE2C/example/example_cor.txt'
# log_file='/Users/cuelee/Dropbox/github/RE2C/out.log'

## print log on bash and log_file
sink( log_file , append = TRUE, split = TRUE)

## print number of cores used in this analysis
cat('\n Number of CPU cores used in this analysis:', ncores, 'cores')

## read RData
Rdata.path <- file.path(work_path,"data/RE2C.RData")
if (file.exists(Rdata.path)){ 
  ## RE2C uses tabulated tables
  load(Rdata.path)
  cat("\n RE2C.RData file has been loaded successfully")
} else{ stop("\n Reading .RData is unsuccessful")}


## load RE2C
source(paste(work_path,"/R/RE2C.R",sep=""),echo=F,local=TRUE)

## load LS 
source(paste(work_path,"/R/LS.R",sep=""),echo=F,local=TRUE)

## check the first line of the input_file
titleconn <- file(input_file,"r")
title.line <- readLines(titleconn, n=1); close.connection(titleconn)
title.dat <- strsplit(title.line,split="\\s+")[[1]]
if(length(title.dat[-1])%%2==0){n_input<-length(title.dat[-1])/2
  } else{stop("\n The input file format should be RSID beta1 stders1 beta2 stders2... ")}

which_betas <- (1:n_input)*2-1;   which_stders <- (1:n_input)*2
remove(title.line); remove(title.dat)

## load NR 
NR <- test_n(n_input)

## check correlation matrix file name
if(is.na(cor_file)) {cor <- diag(n_input)} else {cor <- as.matrix(read.table(cor_file,header=F))}

## check diagonality of the cor matrix 
isdiag <- all(cor[!diag(nrow(cor))] == FALSE)

## RE2C uses tabulated matrices to correct stat2_cdf and FEt.lows
correction.function <- function(correlationMatrix,RE2Cor,tau2prob_cor){
  mean_cor <- mean(correlationMatrix[lower.tri(correlationMatrix)])
  int_corr <- floor(mean_cor*10)
  float_corr <- (mean_cor*10) - int_corr
  corr_coeff <- int_corr+1
  aset <- c(corr_coeff,corr_coeff+1)
  #diff(matrix(x)) If x is a matrix then the difference operations are carried out on each column separately.
  RE2C_corr<- RE2Cor[corr_coeff,] + diff(RE2Cor[aset,])*float_corr
  tau2.zero.prob_corr <- tau2prob_cor[corr_coeff]+ diff(tau2prob_cor[aset])*float_corr
  tau2.zero.prob <- raw.tau2.zero.prob*tau2.zero.prob_corr
  return(list(RE2C_corr,tau2.zero.prob))
}

## correction factors
RE2Cor <- RE2Cor.list[[min(n_input-1,6)]]
tau2prob_cor <- tau2prob_corlist[[min(n_input-1,6)]]

## run correction.function(cor)
correction.list <- correction.function(cor,RE2Cor,tau2prob_cor)
RE2C_corr <- correction.list[[1]];   tau2.zero.prob <- correction.list[[2]]

## read and run RE2C
cat("\n Run RE2C... \n")
inst <- proc.time()
fileconn <- file(input_file,"r")

## create the column label
if(isdiag){
  write(c("rsID", "Nstudy", "LSs", "LSse", "LSp", "I2", "Q", "Qp", "RE2Cs1", "RE2Cs2", "RE2Cp*", "RE2Cp"),file=output_file,ncolumns = 12,append=T)
  } else {
  write(c("rsID", "Nstudy", "LSs", "LSse", "LSp", "RE2Cs1", "RE2Cs2", "RE2Cp*", "RE2Cp"),file=output_file,ncolumns = 9,append=T)}

res = mclapply(readLines(fileconn, warn = FALSE), FUN = parse_aLine, mc.silent = F, mc.cores = ncores)

for (i in 1:length(res)){
  if(res[[i]]$out){
    write(res[[i]]$text, file = output_file, append = T)
  }
}

close.connection(fileconn)
cat("\n Analysis Completed.\n")
cat(" Computation time: ",paste((proc.time()-inst)[3],"sec  \n\n",collapse=" "))


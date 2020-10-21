##functions
getparam<-function(ofile){  
  Res = read.table(ofile)
  nu  = Res[1,]
  Res = Res[-1,]
  M   = ncol(Res)
  ss  = seq(1,M*K,by=K)
  ee  = c((ss-1)[2:M],M*K)
  a   = array(data=NA,dim=c(M,M,K))
  for (m in 1:M){
    for (ell in 1:M){
      a[m,ell,] = Res[ss[ell]:ee[ell],m]
    }
  }
  return(a)  
}

##merge forward and bacward influence
merge_fb <-function(af,ab){  
  K    = dim(af)[3]
  M    = dim(af)[1]
  aa   = array(data=NA,dim=c(M,M,2*K))
  for (m in 1:M){
    for (ell in 1:M){
      aa[m,ell,] = c(ab[m,ell,K:1],af[m,ell,1:K])
      }
  }
  return(aa)
}


##launch hawkes exec
hawkes<-function(dim ,bed_files,outputName){
    ##define dimension hawkes parameters
    if(dim == "f"){
        arg="-f"
        dimension="forward"
    } else if (dim == "b"){
        arg="-b"
        dimension="backward"
    }
    
    ##write hw option file
    hk_file_option = paste(rep(arg,length(my_bed_files)),
                           my_bed_files,collapse=" ")
    
    ##build output file name
    file_out<-  str_c(outdir,outputName,sep ="/")%>%
                  str_c(dimension,sep="_")%>%
                    str_c('txt',sep=".")
    
    ##command assembly
    command = str_c(exec, 
                    hk_file_option,
                    "-histogram",
                    K, delta,
                    "-kernel",kernel,
                    "-lambda",lambda,
                    ">",file_out,
                    sep=" ") 
    ##run cmd
    system(command)
    
    return(file_out)
}



##plot resultats
plot_res<-function(bed_files,a){
  
  x    = c(-seq(delta,K*delta,by=delta)[K:1],seq(0,(K-1)*delta,by=delta)[1:K]) 
  M    = length(bed_files)
  
  ##bed name for plot title
  pattern <- ".+/(.+)\\.\\w+$"
  bed_files.names<-sub(pattern, replacement = "\\1", my_bed_files)
  
  ##plot res
  par(mfrow=c(M,M))
  for (m in 1:M){
    for (ell in 1:M){    
      ymin = min(a[m,ell,])
      ymax = max(a[m,ell,])
      plot(x,a[m,ell,],type="s",ylim=c(min(a) , max(a)),xlab="",ylab="",main=paste("influence of ",bed_files.names[ell]," on ",bed_files.names[m],sep=""))
      #plot(x,a[m,ell,],type="s",ylim=c(ymin,ymax),xlab="",ylab="",main=paste("influence of ",bed_files.names[ell]," on ",bed_files.names[m],sep=""))       
      abline(h=0,col="gray")
      
      lines(x,a[m,ell,],type="s",ylim=c(ymin,ymax),xlab="",ylab="","s")        
    
  }
  }
}


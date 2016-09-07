artificial.levels <-function (y, rem, time, tail = TRUE, type = 1) 
{

    l <- length(rem)
    J <- log2(l+2)	# to get max no. scales.  Assumes (n-2) lifting steps.

    remove<-NULL

    p<-list()

    if(type == 1){ 	
	p <- artlev(y,rem)
    }
    else{
	if(type==2){
		    a0 <- 0.5
	}
	else{
		a0 <- min(diff(time))
	}
	
	boundaries<-a0 * 2^(0:J) #vector to define scale intervals

	boundaries[1]<-0

	for(j in 1:J){
        	r <- which( (y < boundaries[j+1]) & (y>= boundaries[j]) )
        	p[[j]]<-rem[sort(r)]
	}

	l<-sapply(p,length)

	l10<-which(l<10)

	if(tail){
		if(length(l10)>0){
                	lm<-min(l10)
                	a<-unlist(p[lm:length(p)])
                	p<-p[-(lm:length(p))]
                	p[[lm]]<-a
        	}
	}
    }
    
    p
}




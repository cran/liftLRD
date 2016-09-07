liftHurst <-function(x, grid = 1:length(x), model = "FGN", ntraj = 50, tradonly = FALSE, cutoffs=0, cut.fine = TRUE, 
efun = meanmo, afun = idj, altype = 1, tail = TRUE, normalise = TRUE, level=0.05, bc = TRUE, vc = TRUE, jsc = TRUE, BHonly = TRUE, verbose = FALSE,...){

levvar<-function(nj,max=5000){

z2njh<-sum(1/((0:max)+nj/2)^2)

v<-z2njh/(log(2))^2

v

}


gj<-function(nj){

bias<-digamma(nj/2)/log(2) - log2(nj/2)

bias
}

if(!is.matrix(x)){x<-matrix(x)}

if(ncol(x)==2){
	grid<-x[,1]
	x<-x[,2]
}

n<-length(x)

trajmat<-matrix(0,ntraj,n-2)

weights <-NULL

beta <-NULL

betamat <-NULL

emat<-scmat<- NULL

jstarc<-1

if(bc){
	if(verbose){
		cat("doing bias correction...\n")
	}
}

if(vc){
	if(verbose){
	        cat("doing weighted regression...\n")
	}
}

if(jsc){
	if(verbose){
		cat("doing j* coefficient computation...\n")
	}
}

# generate trajectories

if(!tradonly){
	for(i in 1:ntraj){
        	trajmat[i,] <-sample(1:n, n-2, replace = FALSE)
	}
	if(verbose){
		cat("generated trajectories.\n")
	}

	for(k in 1:ntraj){
		if(verbose){
			if((k%%5)==0){
				cat(k,"...")
			}
		}
		xlift <-fwtnpperm(grid, x, mod = trajmat[k,], do.W = TRUE, ...)
		v <-sqrt(diag(tcrossprod(xlift$W)))
		rem <-xlift$removelist
		coeff <-xlift$coeff
		if(normalise){
			coeff <-coeff/v
		}
		scales <-xlift$lengthsremove

		scalesa <-rep(NA,n)
		scalesa[rem] <-scales

		al <-artificial.levels(scales, rem, grid, tail = tail, type = altype)
		levno <-length(al)

		energies<-energiesu<-NULL
		scalesx<-NULL
		
		# for jsc 
                l2vec<-NULL
                levvec<-NULL

		nj <- sapply(al,length)

		for(j in 1:levno){
			# for efun, use: mean2, mad2 ...
			energies[j] <-efun(coeff[al[[j]]])

			# for afun, use: meanj, medj, idj
			scalesx[j] <-afun(scalesa[al[[j]]], j = j)

			# now account for approximate linear relationship between log2(I) and (art.) level j*
                	if(jsc){
                        	l2vec<-c(l2vec,log2(scalesa[al[[j]]]))
                        	levvec<-c(levvec,rep(scalesx[j],times=nj[j]))
                	}
		}


		if(jsc){
        		#### do j* coefficient estimation
        		jstarc<-lm(l2vec~log2(levvec))$coef[2]
		}

		# get bias correction.
                if(bc){
                        gg <- sapply(nj,gj)
                }
                else{
                     	gg<-0
                }

                # get variance weights for regression.
                if(vc){
                         levw <- 1/sapply(nj,levvar)
                }
                else{
                     	levw<-rep(1,times=levno)
                }


		for(i in 1:length(cutoffs)){
			if(cut.fine){
				index<-(1+cutoffs[i]):levno
			}
			else{
				index<- 1:(levno-cutoffs[i])
			}

			beta[i]<-lm( (log2(energies)-gg)[index]~log2(scalesx[index]),weights=levw[index])$coef[2]
			if(jsc){
#		        	cat("\nmultiplying by...",jstarc,"\n")
				beta[i]<-beta[i]*jstarc
			}
		}
		betamat<-rbind(betamat,beta)
	}
	cat("\n")

	beta <-apply(betamat, 2, mean)

	Hs <-Hfrombeta(betamat, model = model)
	sds <-apply(Hs, 2, sd)

	H <-apply(Hs, 2, mean)

	ci <-t(apply(Hs,2,bootci,level=level))
        bandH <- cbind(beta, H, sds, ci)

}
else{
	xnlift <-fwtnp(grid, x, do.W = TRUE, ...)

	v <-sqrt(diag(tcrossprod(xnlift$W)))
	rem <-xnlift$removelist
	coeff <-xnlift$coeff
	if(normalise){
		coeff <-coeff/v
	}

	scales <-xnlift$lengthsremove

	scalesa <-rep(NA, n)
	scalesa[rem] <-scales

	al <-artificial.levels(scales, rem, grid, tail = tail, type = altype)
	levno <-length(al)

	energies <-energiesu<-NULL
	scalesx <-NULL

	# for jsc 
        l2vec<-NULL
        levvec<-NULL

	nj <- sapply(al,length)

	for(j in 1:levno){
		# for efun, use: mean2, mad2 ...
		energies[j] <-efun(coeff[al[[j]]])

		# for afun, use: meanj, medj, idj
		scalesx[j] <-afun(scalesa[al[[j]]], j = j)

		# now account for approximate linear relationship between log2(I) and (art.) level j*
               	if(jsc){
                       	l2vec<-c(l2vec,log2(scalesa[al[[j]]]))
                       	levvec<-c(levvec,rep(scalesx[j],times=nj[j]))
               	}
	}


	if(jsc){
        	#### do j* coefficient estimation
        	jstarc<-lm(l2vec~log2(levvec))$coef[2]
	}

	# get bias correction.
        if(bc){
                gg <- sapply(nj,gj)
        }
        else{
        	gg<-0
        }

        # get variance weights for regression.
        if(vc){
                levw <- 1/sapply(nj,levvar)
        }
        else{
        	levw<-rep(1,times=levno)
        }


	for(i in 1:length(cutoffs)){
		if(cut.fine){
			index<-(1+cutoffs[i]):levno
		}
		else{
			index<- 1:(levno-cutoffs[i])
		}

		beta[i]<-lm( (log2(energies)-gg)[index]~log2(scalesx[index]),weights=levw[index])$coef[2]
		if(jsc){
#	        	cat("\nmultiplying by...",jstarc,"\n")
			beta[i]<-beta[i]*jstarc
		}
	}

	H <-Hfrombeta(beta, model = model)

	bandH<-cbind(beta, H)
}

	if(BHonly){
		return(bandH)
	}
	else{
		return(list(bandH=bandH,energies=(log2(energies)-gg)[index], scales = log2(scalesx[index]), weights = levw[index], 
		jstarc = jstarc  ))
	}
	
}

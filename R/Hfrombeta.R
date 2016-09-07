Hfrombeta <-function(beta, model = c("FBM", "FGN", "ID")){

if(model=="FBM"){
	H <- abs((beta-1)/2)
}
else{
	H <- (beta+1)/2
}

H
}

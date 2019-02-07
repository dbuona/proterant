
extract_coefs<-function(x){
  rownames_to_column(as.data.frame(x$coefficients),"trait")
  
}

extract_CIs<-function(x){
  filter(rownames_to_column(as.data.frame(t(as.data.frame(x$bootconfint95))),"trait"),trait!="alpha")
  
}  
   

as.data.frame(t(as.data.frame(z.funct.drought.noint$bootconfint95))),"trait"),trait!="alpha")



  z.funct.drought.noint$bootconfint95

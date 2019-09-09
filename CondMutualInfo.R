computeCmi <- function(v1,v2,vcs)
{
  if(nargs() == 2)
  {
    c1 <-  var(v1)
    c2 <-  var(v2)
    v <-cbind(v1,v2)
    c3 <- det(cov(v))
    cmiv <- 0.5*log(c1*c2/c3) 
  }
  else if(nargs() == 3)
  {   
    v <-cbind(v1,vcs)
    c1=det(cov(v))
    v <-cbind(v2,vcs)
    c2=det(cov(v))
    if(ncol(vcs >1))
       {
         c3=det(cov(vcs))
    }
    else
    {
      c3 =var(vcs)
    }
    v <-cbind(v1,v2,vcs)
    c4=det(cov(v))
    cmiv=0.5*log((c1*c2)/(c3*c4))       	
  }
  return (cmiv)
}
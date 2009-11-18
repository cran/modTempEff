seas<-function(x,ndx=stop("please, provide `ndx' in seas()")){
        r<-x
        attr(r,"ndx")<-ndx
        return(r)
        }

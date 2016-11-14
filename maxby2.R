genelist<-"genelist400k.txt"
paths<-c("/pbtech_mounts/mezeylab_store/abs2017/tcgaov","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/brca","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/thca","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/lgg","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/luad","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/kirc","/pbtech_mounts/mezeylab_store/abs2017/tcgapros","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/coad","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/lusc","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/ucec","/pbtech_mounts/mezeylab_store/abs2017/tcgacancers/hnsc")

for( ii in 1:11) {
eval(parse(text=paste('bclist<-read.table("',paths[ii],'/',genelist,'",stringsAsFactors=F)',sep='')),envir=topenv())

for ( jj in (ii+1):11) {
   eval(parse(text=paste('colist<-read.table("',paths[jj],'/',genelist,'",stringsAsFactors=F)',sep='')),envir=topenv())

     common<-numeric(0)
     s2<-unique(unlist(colist))
     s1<-unique(unlist(bclist))
     nambchord<-intersect(s1,s2)
     nambchord<-as.data.frame(nambchord,stringsAsFactors=F)
   
for( kl in 1:nrow(nambchord)) {
   if( length(which(bclist == nambchord[kl,1])) !=0 &&  length(which(colist == nambchord[kl,1])) !=0 && nambchord[kl,1] != "?") {

   bcc<-bclist[which(bclist == nambchord[kl,1],arr.ind=T)[,1],]   
   coc<-colist[which(colist == nambchord[kl,1],arr.ind=T)[,1],] 
                                                           
   for ( i in 1:nrow(bcc)) {
     if( bcc[i,1] !=nambchord[kl,1]) {
     temp<-as.character(bcc[i,1])
     bcc[i,1]<-as.character(bcc[i,2])
     bcc[i,2]<-as.character(temp)
    }
   }

   for ( i in 1:nrow(coc)) {
     if( coc[i,1] !=nambchord[kl,1]) {
     temp<-as.character(coc[i,1])
     coc[i,1]<-as.character(coc[i,2])
     coc[i,2]<-as.character(temp)
     }
   }
   common<-rbind(common,length(intersect(bcc[,2],coc[,2])))
   }
   else common<-rbind(common,0)
   
   if( kl%% 100 == 0) print(kl)
  }


  ids<-which(common == max(common))
  for( counter in 1:length(ids)) {
  kl = ids[counter]
  bcc<-bclist[sort(which(bclist == nambchord[kl,1],arr.ind=T)[,1]),]   
  coc<-colist[sort(which(colist == nambchord[kl,1],arr.ind=T)[,1]),] 
                                                           
  for ( i in 1:nrow(bcc)) {
    if( bcc[i,1] !=nambchord[kl,1]) {
    temp<-as.character(bcc[i,1])
    bcc[i,1]<-as.character(bcc[i,2])
    bcc[i,2]<-as.character(temp)
    }
  }
  
  for ( i in 1:nrow(coc)) {
     if( coc[i,1] !=nambchord[kl,1]) {
     temp<-as.character(coc[i,1])
     coc[i,1]<-as.character(coc[i,2])
     coc[i,2]<-as.character(temp)
     }
  }
   
   tl<-character()
   tl<-cbind(tl,intersect(bcc[,2],coc[,2])) 
   tl<-rbind(tl,"tfactor is")
   tl<-rbind(tl,nambchord[kl,1])
  
   filename<-paste(gsub('.*\\/','',paths[ii]),gsub('.*\\/','',paths[jj]),counter)
   eval(parse(text=paste('write.table(tl,"/pbtech_mounts/mezeylab_store/abs2017/coexpcommoncancers/dc/400k/by2/',filename,'",sep="\t",col.names=F,row.names=F,quote=F)',sep='')),envir=topenv())
  }
  }#jj loop
  }#ii loop
     
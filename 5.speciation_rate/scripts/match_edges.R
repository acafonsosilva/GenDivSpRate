map_rates = function (sub_tree, sub_rates, whole_tree, rates = rep(NaN, nrow(whole_tree$edge))){
  sub_n = sub_tree$Nnode + 1
  
  matched_nodes = c(sapply(sub_tree$tip.label, function(x){which(whole_tree$tip.label == x)}),rep(-1,sub_tree$Nnode))
  matched_edges = rep(-1, nrow(sub_tree$edge))
  
  for(i in nrow(sub_tree$edge):1){
    sub_end_node = sub_tree$edge[i,2]
    sub_base_node = sub_tree$edge[i,1]
    
    whole_end_node = matched_nodes[sub_end_node]
    whole_edge = which(whole_tree$edge[,2] == whole_end_node)
    
    matched_edges[i] = whole_edge
    if(matched_nodes[sub_base_node] == -1){
      matched_nodes[sub_base_node] = whole_tree$edge[whole_edge,1]
    }
  }

  rates[matched_edges] = sub_rates
  return(rates)
}

map_rates_tipNroot = function (sub_tree, sub_rates, whole_tree, rates = rep(NaN, nrow(whole_tree$edge))){
  sub_n = sub_tree$Nnode + 1
  whole_n = whole_tree$Nnode + 1
  
  matched_nodes = c(sapply(sub_tree$tip.label, function(x){which(whole_tree$tip.label == x)}),rep(-1,sub_tree$Nnode))
  matched_edges = rep(-1, nrow(sub_tree$edge))
  
  for(i in nrow(sub_tree$edge):1){
    sub_end_node = sub_tree$edge[i,2]
    sub_base_node = sub_tree$edge[i,1]
    
    whole_end_node = matched_nodes[sub_end_node]
    whole_edge = which(whole_tree$edge[,2] == whole_end_node)
    
    matched_edges[i] = whole_edge
    if(matched_nodes[sub_base_node] == -1){
      matched_nodes[sub_base_node] = whole_tree$edge[whole_edge,1]
    }
  }
  
  rates[matched_edges] = sub_rates
  tips = matched_nodes[1:sub_n]
  root = matched_nodes[sub_n + 1]
  
  return(list(rates=rates,tips = tips, root = root))
}

plot.with.rate.withNaNs=function(phylo,rate1,NaN_color="black",rate2=NULL,same.scale=T,main=NULL,lwd=1,log=F,scale=NULL,leg=F,show.tip.label = F, ...){
  Colors = scico(100, palette = 'batlow')
 # Colors = colorRampPalette(rev(c('darkred',brewer.pal(n = 8, name = "Spectral"),'darkblue')))(100)
 # Colors = colorRampPalette(c("blue4","steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 )
  if(is.null(rate2)){
    if(log) {
      rate1=log(rate1)
      if(! is.null(scale)) scale=log(scale)
    }
    if(isTRUE(all.equal(rep(as.numeric(rate1[1]),length(rate1)),as.numeric(rate1)))){
      col=rep(1,length(rate1))
      plot.phylo(phylo, type="fan", direction = 'upwards', edge.color = Colors[col], show.tip.label = F,main=main,edge.width =lwd)
      if(log){
        if(leg) image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T, legend.shrink	= 0.5)
      }else{
        if(leg) image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T, legend.shrink	= 0.5)
      }
    }else{
      col = round( (rate1 - min(c(scale,rate1),na.rm=T)) / diff(range(c(scale,rate1),na.rm=T))*99   )+1
      col=Colors[col]
      col[is.na(col)]=NaN_color
      plot.phylo(phylo, type="fan", direction = 'upwards', edge.color = col, show.tip.label = show.tip.label,main=main,edge.width =lwd,...)
      if(log){
        min=min(c(scale,rate1),na.rm = T)
        max=max(c(scale,rate1),na.rm = T)
        m10=floor(min/log(10))
        M10=ceiling(max/log(10))
        if((M10-m10)<4){
          ticks=c(1,2,5)
        }else{
          ticks=1
        }
        ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
        lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
        if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
        if(leg) {
          image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks), legend.shrink	= 0.5)
        }else{
          return(list(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks), legend.shrink	= 0.5))
        }
      }else{
        if(leg) image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T, legend.shrink	= 0.5)
      }
    }
  }else{
    if(log){
      rate1=log(rate1)
      rate2=log(rate2)
    }
    if(same.scale){
      min=min(min(rate1),min(rate2))
      max=max(max(rate1),max(rate2))
      par(mfrow=c(1,2))
      col = round(( (rate1 - min) / (max-min))*99   )+1
      plot.phylo(phylo, type="fan", direction = 'upwards', edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
      col = round(( (rate2 - min) / (max-min))*99   )+1
      plot.phylo(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
      par(mfrow=c(1,1))
      if(log){
        m10=floor(min/log(10))
        M10=ceiling(max/log(10))
        if((M10-m10)<4){
          ticks=c(1,2,5)
        }else{
          ticks=1
        }
        ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
        lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
        if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
        # ticks=seq(min,max,length.out = 5)
        # image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
      }else{
        # image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T)
      }
    }else{
      par(mfrow=c(1,2))
      if(isTRUE(all.equal(rep(rate1[1],length(rate1)),rate1))){
        col=rep(1,length(rate1))
        plot.phylo(phylo, type="fan", direction = 'upwards', edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          
          # image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
        }else{
          # image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
        }
      }else{
        col = round(( (rate1 - min(rate1)) / (max(rate1)-min(rate1)))*99   )+1
        plot.phylo(phylo, type="fan", direction = 'upwards', edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          min=min(rate1)
          max=max(rate1)
          m10=floor(min/log(10))
          M10=ceiling(max/log(10))
          if((M10-m10)<4){
            ticks=c(1,2,5)
          }else{
            ticks=1
          }
          ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
          lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
          if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
          # image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
        }else{
          # image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
        }
      }
      if(isTRUE(all.equal(rep(rate2[1],length(rate2)),rate2))){
        col=rep(1,length(rate2))
        plot.phylo(phylo, type="fan", direction = 'upwards', edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          # image.plot(z = c(exp(rate2[1]),2*exp(rate2[1])),col = Colors, horizontal=T,legend.only = T)
        }else{
          # image.plot(z = c(rate2[1],2*rate2[1]),col = Colors, horizontal=T,legend.only = T)
        }
      }else{
        col = round(( (rate2 - min(rate2)) / (max(rate2)-min(rate2)))*99   )+1
        plot.phylo(phylo, type="fan", direction = 'upwards', edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          min=min(rate2)
          max=max(rate2)
          m10=floor(min/log(10))
          M10=ceiling(max/log(10))
          if((M10-m10)<4){
            ticks=c(1,2,5)
          }else{
            ticks=1
          }
          ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
          lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
          if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
          # image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
        }else{
          # image.plot(z = as.matrix(rate2),col = Colors, horizontal=T,legend.only = T)
        }
      }
    }
    par(mfrow=c(1,1))
  }
  return(Colors[col])
}


# if(F){
#   # an example
#   whole_tree = read.tree("~/ownCloud/Lab Folder/Odile/Afrotheria.tree")
#   
#   sub_trees=subtrees(whole_tree)
#   
#   # sub_tree need to be a sub_clade of whole_tree, if you have some groups that are not monophyletic tell me,
#   # the function will need modifications.
#   sub_tree1=sub_trees[[10]]
#   sub_rates1=rlnorm(nrow(sub_tree1$edge), 2,0.4)
#   new_rates=map_rates(sub_tree1, sub_rates=sub_rates1, whole_tree, rates=rep(NaN,nrow(whole_tree$edge)))
#   
#   sub_tree2=sub_trees[[50]]
#   sub_rates2=rlnorm(nrow(sub_tree2$edge), 2,0.4)
#   new_rates=map_rates(sub_tree2, sub_rates=sub_rates2, whole_tree, rates=new_rates)
#   
#   par(mfrow=c(1,3))
#   plot.with.rate.withNaNs(sub_tree1, sub_rates1, leg=T, lwd = 2, log=T)
#   plot.with.rate.withNaNs(sub_tree2, sub_rates2, leg=T, lwd = 2, log=T)
#   plot.with.rate.withNaNs(whole_tree, new_rates, NaN_color = "gray90", leg=T, lwd = 2, log=T)
# }

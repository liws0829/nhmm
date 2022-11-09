library(RColorBrewer)
library(reshape2)
Colorbar = c("#CE3C35","#1E90FF","#006400","#FF1493","#9932CC","#FF8C00")
Colorbar = c("#CE3C35","#006400","#9932CC","#FF8C00")
Colorbar = c("#cc7892","#7d9c41","#00a0d8","#ab8f28","#564147","#bea5ab")
Colorbar = c("#564147","#cc7892","#bea5ab")
Colorbar = rev(c("#F7FBFF","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B","darkred"))
rColorbar = rev(Colorbar)
#plot graph
draw_graph_topo<-function(G,seed=1)
{
  #no corelation
  nw=graph_from_adjacency_matrix(G,mode='directed',diag=F)
  set.seed(seed)
  plot(nw,edge.arrow.size=1,layout=layout.gem(nw))
}


##plot graph of 2 time with state transition

# V(nw)$color = factor(apply(data$lt[,,2],1,function(x) which.is.max(x)))
# plot(nw)
# set.seed(seed)
# plot(nw,edge.arrow.size=0.5)


# #plot W graph in heatmap
# gr = ggplot(melt(THETA$W),aes(x=as.factor(Var1),y=as.factor(Var2),fill=value))+
#   geom_tile()+
#   coord_fixed(expand=FALSE)+
#   scale_fill_viridis_c(option="A",begin=0)+
#   xlab("Node")+ylab("Node")+theme(legend.title = element_text("value of interaction"))
# show(gr)

#plot data
##plot density
draw_obs_dens<-function(obs,ltn)
{
  q_ = dim(obs)[2]
  n_ = dim(obs)[1]
  t_ = dim(obs)[3]
  df_data = data.frame(as.vector(aperm(obs,c(1,3,2))))
  df_data$obs = factor(rep(res_name,each=n_*t_))
  df_data$ltn = factor(rep(ltn,q_))

  names(df_data)=c('Value','Observation','State')
  gr = ggplot(data=df_data,aes(x=Value,fill=State))+
    geom_density(alpha=0.5,position = "identity")+
    #scale_fill_manual(values = c("#Ce3c35","#4258a1"))+
    facet_grid(State~Observation)+
    #coord_cartesian(xlim=c(min(df_data$Value),max(df_data$Value)),expand=F)+
    theme(text = element_text(size=14),
        panel.background = element_rect(fill="white",color="black"),
        strip.background = element_rect(color="black")
        )+xlim(c(-5,5))
  show(gr)
}
draw_obs_time<-function(obs,q,stepsize=18,nodelist)
{
  #obs:T*N
  obs = obs[nodelist,q,]
  obs = as.data.frame(t(obs))
  stepsize=max(obs,na.rm = T)*2
  t<-seq(1,nrow(obs))
  colnames(obs) = as.character(nodelist)
  mydata<-melt(as.data.frame(cbind(t,obs)),id="t")
  mydata$variable<-factor(mydata$variable)
  mydata$offset = -as.numeric(mydata$variable)*stepsize
  mydata$density_offset = mydata$value+mydata$offset

  gr = ggplot(mydata,aes(t,density_offset,color=variable))+
    #geom_ribbon(aes(t,ymin=offset,ymax=density_offset,color=variable),colour=NA)+
    geom_line(aes(group=variable,color=variable),lwd=0.8)+
    scale_y_continuous(breaks=seq(-stepsize,-stepsize*ncol(obs),-stepsize),labels=colnames(obs))+
    scale_color_brewer(palette = "Paired")+
    xlab("Time")+
    ylab("Node")+
    theme_classic()+
    theme(
      panel.background = element_rect(fill="white",color=NA),
      panel.grid.major.x = element_line(color = "grey80",size=0.25),
      text=element_text(size=15),
      legend.position = "none"
    )+
    annotate("text",x=0.8*nrow(obs),y=0.5,label=paste("Observation",as.character(q)))
  show(gr)
}


draw_obs_dens_st<-function(obs)
{
  q_ = dim(obs)[2]
  n_ = dim(obs)[1]
  t_ = dim(obs)[3]
  df_data = data.frame(as.vector(aperm(obs,c(1,3,2))))
  #df_data$obs = factor(paste("Response",rep(1:q_,each=n_*t_)))
  df_data$obs = factor(rep(res_name,each=n_*t_))
  #df_data$state = factor(paste("State",rep(ltn,q_)))
  names(df_data)=c('Value','Observation')

  gr = ggplot(df_data,aes(x=Value,fill=Observation)) +
    #geom_density(position="identity",color="black") +
    geom_histogram(color="black",bins = 50)+
    facet_grid(.~ Observation) +
    scale_fill_manual(values = c("#Ce3c35","#4258a1","#57b058","#7b4b99","#46a1cd")) +
    panel_border() +
    ylab("Frequency") +
    theme(panel.background = element_rect(fill="white",color="black"),
          text = element_text(size=14),panel.border = element_rect(color="black"),
          strip.background = element_rect(colour = "black"))+
    guides(fill = "none")
  show(gr)
}
##plot time series

draw_series<-function(obs,nodelist)
{
  n_ = length(nodelist)
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  obs = as.data.frame(as.vector(obs[nodelist,,]))
  obs$Node=factor(paste("Node",rep(nodelist,q_*t_)))
  #obs$Response=factor(paste("Response",rep(rep(1:q_,each=n_),t_)))
  obs$Response=factor(rep(rep(names(yt)[-1],each=n_),t_))
  obs$Time = rep(1:t_,each=n_*q_)
  names(obs)[1]="Value"
  gr = ggplot(obs,aes(x=Time,y=Value,color=Node)) +
    geom_line(lwd=1) +
    facet_grid(Node ~ Response) +
    panel_border() +
    scale_color_manual(values = c("#Ce3c35","#4258a1","#57b058","#7b4b99","#46a1cd")) +
    theme(panel.background = element_rect(fill="white",color="black"),
          text = element_text(size=14),legend.position = "top",
          panel.border = element_rect(color="black"),
          strip.background = element_rect(color="black"))
  show(gr)
}

draw_rmse<-function(df,md)
{
  # q_ = ncol(df)/length(md)
  # se = paste("(",as.character(round(as.vector(apply(df,2,sd)),3)),")",sep="")
  # df = colMeans(df)
  # df = round(df,3)
  # df = as.data.frame(as.vector(df))
  # df$Model= factor(rep(md,each=q_),levels = md)
  # df$Response=factor(paste("Response",rep(1:q_,length(md))))
  # df$SE = se
  # names(df)[1]="RMSE"
  
  gr = ggplot(df,aes(x=Model,y=RMSE,fill=Response))+
    geom_bar(stat="identity",position = position_dodge(0.7),width = 0.4,color="black")+
    geom_text(aes(label=RMSE,y=RMSE+0.5),position = position_dodge(0.7),vjust=0,angle=90,cex=4)+
    #geom_text(aes(label=SE,y=RMSE+0.05),position = position_dodge(1),vjust=0)+
    panel_border() +
    scale_fill_manual(values = c("#Ce3c35","#4258a1","#57b058","#7b4b99","#46a1cd")) +
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(ylim=c(0,max(df$RMSE)+1)) +
    theme(panel.background = element_rect(fill="white",color="black"),
          legend.key.size = unit(6,"pt"), legend.text = element_text(size=12),
          axis.title = element_text(size=14), axis.text = element_text(size=12),
          legend.title = element_text("Response"), legend.position = c(0.15,0.8), 
          legend.direction = "vertical", legend.background = element_rect(linetype = "solid",color="black"))+
    xlab("Model")+ylab("RMSE")
  show(gr)
}


draw_comp<-function(df,md,n_,t_,q_)
{
  
  df = df[,which(names(df)%in%md)]
  df$Node = factor(paste("Node",rep(1:n_,q_*t_)))
  #df$Response = factor(paste("Response",rep(rep(1:q_,each=n_),t_)))
  df$Response = factor(rep(rep(res_name,each=n_),t_))
  df$Time = rep(1:t_,each=n_*q_)
  
  df = reshape2::melt(df,id.vars=c("Node","Response","Time"),measure.vars=md,
                      variable.name="Model", value.name = "Value")
  
  gr = ggplot(df[df$Model%in%c("Real","NHMM"),],aes(x=Value,color=Model,lty=Model))+
    geom_density(position="identity",lwd=1)+
    facet_grid(.~Response)+
    scale_color_manual(values = c("black","#CE3C35"))+
    scale_linetype_manual(values=c("solid","dashed"))+
    theme(text = element_text(size=12),legend.key.size = unit(5,"pt"),
          panel.background = element_rect(fill="white",color="black"),
          panel.border = element_blank(), strip.background = element_rect(color = "black"))+
    ylab("Density")+xlim(c(-5,5))
  show(gr)
  
  num = paste("Node",19)
  df = df[(df$Node==num)&(df$Time>=(t_-100+1)),]
  
  gr = ggplot(df[df$Model%in%c("Real","NHMM"),], aes(x=Time,y=Value,color=Model))+
    geom_line(aes(lty=Model),lwd=1)+
    geom_point(aes(shape=Model),size=2)+
    scale_shape_manual(values=c(16,17))+
    facet_grid(Response~.)+
    scale_color_manual(values = c("black",Colorbar))+
    scale_linetype_manual(values = c("solid","dashed"))+
    theme(panel.background = element_rect(fill="white",color="black"),
          text = element_text(size=12),panel.border = element_blank(),
          strip.background  = element_rect(color="black"))+
    labs(title=num)
  show(gr)
  

  gr = ggplot(df[df$Time>=(t_-100+1),],aes(x=Time,y=Value,color=Model))+
    geom_line(aes(lty=Model),lwd=1)+
    geom_point(aes(size=Model),shape=16)+
    #scale_shape_manual(values=c(16,4))+
    scale_size_manual(values = c(1.5,rep(0,6)))+
    facet_grid(Response~.)+
    scale_color_manual(values = c("black",Colorbar))+
    scale_linetype_manual(values = c("blank","solid",rep("dashed",6)))+
    theme(panel.background = element_rect(fill="white",color="black"),
          text = element_text(size=12),
          panel.border = element_blank(),strip.background  = element_rect(color="black"))+
    labs(title=num)
  show(gr)
  
  
  dfv = df[df$Model=="NHMM",]
  dfv$Value = dfv$Value - df[df$Model=="Real",'Value']
  gr = ggplot(dfv,aes(x=Time,y=Value))+
    geom_segment(aes(x=Time,xend=Time,y=0,yend=Value),size=0.5,color="black",lty="solid")+
    geom_point(aes(fill=Response),size=1.2,color="black",shape=21,stroke=1)+
    facet_grid(Response~.)+
    scale_fill_manual(values = c("#Ce3c35","#4258a1","#57b058","#7b4b99","#46a1cd"))+
    theme(panel.background = element_rect(fill="white",color="black"),
          text = element_text(size=12),panel.grid.major.y=element_line(color="#D3D3D3"),
          panel.border = element_blank(),strip.background = element_rect(color="black"))+
    labs(title = num)+ ylab("Residue")+
    ylim(c(min(dfv$Value)-0.1,max(dfv$Value)+0.1))
    
  show(gr)

}

draw_qq<-function(df,md)
{
  
  gr = list()
  for(i in 1:length(md))
  {
    gr[[i]] = ggplot(data=df,aes(sample=df[,i]))+
      geom_qq(alpha=0.7,color=Colorbar[i])+
      geom_qq_line(color="grey")+
      #geom_point(alpha=0.7,color=Colorbar[i])+
    #geom_abline(color="#D62728FF")+
    xlab(md[i])+
    #ylab("Observed -log10(P-value)")+
    #scale_x_continuous(limits = c(0,7))+
    #scale_y_continuous(limits = c(0,7))+
    theme(
      axis.title = element_text(size=12),
      axis.text = element_text(size=8,color = "black"),
      #axis.line = element_line(size=0.8,color="black"),
      axis.ticks= element_line(size=0.8,colour = "black"),
      panel.grid =element_blank(),
      panel.border = element_rect(fill=NA,size = 0.8),
      panel.background = element_blank())
  }
  plot_grid(plotlist = gr,labels="",label_size = 8)
}

#compare two timeseries
draw_pred<-function(obs,pred,stepsize=18,nodelist)
{
  #obs:T*N
  rmse = round(sqrt(mean((obs-pred)^2)),4)
  rmse_text=paste("RMSE =",as.character(rmse),sep=' ')
  obs = as.data.frame(t(obs))
  pred = as.vector(t(pred))
  stepsize=max(obs)*2
  t<-seq(1,nrow(obs))
  colnames(obs) = nodelist
  mydata<-melt(as.data.frame(cbind(t,obs)),id="t")

  mydata$variable<-factor(mydata$variable)
  mydata$offset = -as.numeric(mydata$variable)*stepsize
  mydata$density_offset = mydata$value+mydata$offset
  mydata$pred = pred+mydata$offset

  gr = ggplot(mydata)+
    #geom_ribbon(aes(t,ymin=offset,ymax=density_offset,color=variable),colour=NA)+
    geom_line(aes(x=`t`,y=`density_offset`,group=`variable`,color=`variable`),lwd=0.8,alpha=0.8)+
    geom_line(aes(x=`t`,y=`pred`,group=`variable`),color="black",lwd=0.8,lty=4)+
    scale_y_continuous(breaks=seq(-stepsize,-stepsize*ncol(obs),-stepsize),labels=seq(1,ncol(obs),1))+
    xlab("Time")+
    ylab("Node")+
    scale_color_brewer(palette = "Paired")+
    theme_classic()+
    theme(
      panel.background = element_rect(fill="white",color=NA),
      panel.grid.major.x = element_line(color = "grey80",size=0.25),
      text=element_text(size=15),
      legend.position = "none"
    )+
    annotate("text",x=0.8*nrow(obs),y=0.5,label=rmse_text)
  show(gr)
}
#plot VEM
##convergence
draw_VEM_convergence<-function(vemloglik)
{
  g_=length(vemloglik)
  df_loglik = as.data.frame(vemloglik)
  gr = ggplot(df_loglik) + theme_bw() +
    geom_line(aes(x=(1:g_),y=df_loglik[,1]),size=1) +
    scale_x_continuous(name="Iteration",breaks=seq(1,g_,as.integer(0.2*g_)),labels=as.character(seq(1,g_,as.integer(0.2*g_)))) +
    ylab("log-Likelihood")
  show(gr)
}

#draw the inference compare
draw_inf<-function(rlt)
{
  n_ = nrow(rlt)
  t_ = ncol(rlt)
  df_rlt = as.data.frame(as.factor(as.vector(rlt)))   #N*T
  df_rlt$n = rep(1:n_,t_)
  df_rlt$t = rep(1:t_,each=n_)
  names(df_rlt)=c("State","Cell","Time")
  gr = ggplot(data = df_rlt,aes(x=`Time`, y=`Cell`,fill=`State`)) +
    geom_raster() +  
    scale_fill_manual(values = Colorbar)+
    #scale_fill_brewer(palette = "BuGn")+
    #scale_fill_gradient(low="#564147",high = "#bea5ab",breaks=1:K) +
    xlab("Time")+
    scale_y_continuous(name="Cell",breaks = seq(1,n_,as.integer(0.2*n_)),labels = as.character(seq(1,n_,as.integer(0.2*n_))),expand = c(0,0)) +
    scale_x_continuous(expand=c(0,0))+
    coord_fixed(ratio=t_/n_*0.4)+
    theme(text = element_text(size=15),
          legend.position = "bottom")

  show(gr)
}

draw_cor<-function(dat)
{
  if(dim(dat)[2]>1)
  {
    dfdat=data.frame()
    for(q in 1:dim(dat)[2])
    {
      dfdat = rbind(dfdat,dat[,q,])
    }
  }
  cordf = cor(t(dfdat))
  corrplot(cordf,method = "color",tl.pos = "n")
}

draw_boxplot<-function(obs,ltn)
{
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  mydata=data.frame()
  for(q in 1:q_)
  {
    response = rep(q,t_)
    for(n in 1:n_)
    {
      middf = cbind(response,obs[n,q,],ltn[n,])
    }
    mydata = rbind(mydata,as.data.frame(middf))
  }
  colnames(mydata) = c("Response","Observation","State")
  mydata$State = paste("State",mydata$State,sep="_")
  mydata$State = as.factor(mydata$State)
  mydata$Response = as.factor(mydata$Response)

  gr = ggplot(data=mydata,aes(x=`State`,y=`Observation`,fill=`Response`)) +
    geom_boxplot()+
    facet_wrap(~`State`,scales="free")

  show(gr)


}

draw_compare_time<-function(obs,nodenum,t,q,colname_str)
{
  #obs: TQ*7

  n = length(colname_str)
  gr = list()
  for(i in 1:q)
  {
    pobs = data.frame(as.vector(as.matrix(obs[((i-1)*t+1):(i*t),],ncol=1)))
    pobs$Time = rep(1:t,n)
    pobs$Model = factor(rep(colname_str,each=t),levels = c("Real",rev(colname_str[-1])))
    names(pobs)[1] = "Observation"

    gr[[i]] = ggplot(pobs,aes(x=`Time`,y=`Observation`,colour=`Model`,linetype=`Model`)) +
      #geom_point(aes(shape=`Model`)) +
      geom_line(aes(size=`Model`)) +
      scale_size_manual(values=c(1,rep(0.6,(n-2)),0.7)) +
      scale_linetype_manual(values = c("solid",rep("dashed",n-2),"twodash")) +
      scale_colour_manual(values = rColorbar) +
      # scale_linetype_manual(values = c("solid","twodash")) +
      #scale_colour_manual(values = c("grey","red2")) +
      theme(text = element_text(size=10),
            panel.background = element_rect(fill="white",color="black")) +
      ylab(paste("Observation",as.character(i),sep=" "))
  }
  lab = sapply(1:q,function(q) paste("Response",as.character(q),sep=" "))
  plot_grid(plotlist = gr,nrow=q,labels =NULL)
}


draw_compare_bar<-function(obs,nodenum,t,q,colname_str)
{
  #obs: TQ*7

  n = length(colname_str)
  gr = list()
  for(i in 1:q)
  {
    pobs = data.frame(as.vector(as.matrix(obs[((i-1)*t+1):(i*t),],ncol=1)))
    pobs$Time = rep(1:t,n)
    pobs$Model = factor(rep(colname_str,each=t),levels = rev(colname_str))
    names(pobs)[1] = "Residue"

    gr[[i]] = ggplot(pobs,aes(x=`Time`,y=`Residue`,group=`Model`)) +
      #geom_point(aes(shape=`Model`)) +
      geom_col(aes(fill=`Model`)) +
      scale_fill_manual(values = rColorbar[2:7]) +
      theme(text = element_text(size=10),
            panel.background = element_rect(fill="white",color=NULL)
      )+
      ylab(paste("Residue",as.character(i),sep=" "))
  }
  lab = sapply(1:q,function(q) paste("Response",as.character(q),sep=" "))
  plot_grid(plotlist = gr,nrow=q,labels = NULL)
}

draw_compare_dens<-function(obs,nodenum,t,q,colname_str)
{
  #obs: TQ*?

  n = length(colname_str)
  gr = list()

  for(i in 1:q)
  {
    bt = (0:(nodenum-1))*q*t + (i-1)*t + 1
    pobs = c()
    for(ni in bt)
    {
      pobs = c(pobs,as.vector(as.matrix(obs[ni:(ni+t-1),],ncol=1)))
    }
    pobs = as.data.frame(pobs)
    # pobs$Time = rep(1:t,n)
    pobs$Model = factor(rep(colname_str,each=t*nodenum))
    names(pobs)[1] = "Density"

    gr[[i]] = ggplot(pobs,aes(x=`Density`,color=`Model`)) +
      #geom_point(aes(shape=`Model`)) +
      geom_density(position = "identity",lwd=0.8) +
      scale_linetype_manual(values = c("dotted","dotdash","longdash","dashed","dotted","solid","twodash")) +
      scale_color_manual(values = c("blue", "#33a02c", "darkgoldenrod1", "skyblue2", "purple", "red2", "black")) +
      theme(text = element_text(size=10),
            panel.background = element_rect(fill="white",color="black")
      )
  }
  lab = sapply(1:q,function(q) paste("Response",as.character(q),sep=" "))
  plot_grid(plotlist = gr,ncol=q,labels = lab,label_size = 14)
}

dftolist<-function(x,theta)
{
  newtheta = theta
  x = as.numeric(x)
  l = length(x)
  i=1
  j=1
  while(i<l)
  {
    n = length(theta[[j]])

    if(!is.null(dim(theta[[j]])))
    {
      #print(dim(theta[[j]]))
      newtheta[[j]] = array(x[i:(i+n-1)],dim(theta[[j]]))
    }

    else
      newtheta[[j]] = x[i:(i+n-1)]
    names(newtheta[[j]]) = NULL
    i = i + n
    j = j + 1
  }
  return(newtheta)
}

#table saving
table_est<-function(dfmae)
{
  dfmae[nrow(dfmae)+1,] = round(colMeans(dfmae),3)
  dfmae[nrow(dfmae)+1,] = round(apply(dfmae,2,sd),3)
  row1 = nrow(dfmae)-1
  row2 = nrow(dfmae)
  mmean = dftolist(round(dfmae[row1,],3),THETA)
  msd = dftolist(round(dfmae[row2,],3),THETA)
  esdf = data.frame()
  for(i in 1:length(THETA))
  {
    if(!is.null(dim(THETA[[i]])))
    {
      tdf = c()
      for(j1 in 1:nrow(THETA[[i]]))
      {
        for(j2 in 1:ncol(THETA[[i]]))
        {
          if(mmean[[i]][j1,j2]!=0)
          {
            temp = paste(names(THETA)[i],"_{",as.character(j1),",",as.character(j2),"}",sep="")
            mmsd = paste(as.character(mmean[[i]][j1,j2]),"(",as.character(msd[[i]][j1,j2]),")",sep="")
            tdf = rbind(tdf,c(temp,mmsd))
          }
        }
      }
      esdf = rbind(esdf,tdf)
    }
    else
    {
      tdf = c()
      for(j in 1:length(THETA[[i]]))
      {
        temp = paste(names(THETA)[i],as.character(j),sep="_")
        mmsd = paste(as.character(mmean[[i]][j]),"(",as.character(msd[[i]][j]),")",sep="")
        tdf = rbind(tdf,c(temp,mmsd))
      }
      esdf = rbind(esdf,tdf)
    }
  }
  filename = paste("estimate_",as.character(N),"_",as.character(K),"_",as.character(T_),".csv",sep="")
  write.table(esdf,file =filename,sep=",")
}

table_fit<-function(dffit,rowname_str)
{
  n = length(rowname_str)
  q = ncol(dffit)/n
  mmean = matrix(round(apply(dffit,2,mean),3),n,q,byrow=T)
  msd = matrix(round(apply(dffit[-nrow(dffit),],2,sd),3),n,q,byrow=T)
  tdf = data.frame()
  for(i in 1:n)
  {
    temp = c()
    for(j in 1:q)
    {
      temp = cbind(temp,paste(as.character(mmean[i,j]),"(",as.character(msd[i,j]),")",sep=""))
    }
    tdf = rbind(tdf,temp)
  }
  tdf = cbind(rowname_str,tdf)
  filename = paste("fitrmse_",as.character(N),"_",as.character(K),"_",as.character(T_),".csv",sep="")
  write.table(tdf,file=filename,row.names = F,sep = ",")
  return(tdf)
}

boxplot_ari<-function(dfari,model_name,choise)
{
  n = nrow(dfari)
  c = ncol(dfari)
  df = data.frame()
  for(i in 1:c)
  {
    temp = cbind(dfari[,i],model_name[i])
    df = rbind(df,temp)
  }
  names(df) = c("ARI","Model")
  df$Model = factor(df$Model,levels=model_name)
  df$ARI = as.numeric(df$ARI)
  gr = ggplot(df,aes(x=`Model`,y=`ARI`,fill=`Model`)) + geom_boxplot(alpha=0.6) +
    scale_fill_manual(values = Colorbar) +
    theme(text = element_text(size=14),
          panel.background = element_rect(fill="white",color="black"))+
    ylab(choise)
  show(gr)
}

df = data.frame(round(rowMeans(rsq[-1,]),4))
df = cbind(df,factor(rmse_model_name,levels=rmse_model_name[order(df,decreasing = T)]))
names(df) = c("R2","Model")
gr = ggplot(df,aes(y=`R2`,x=`Model`,fill=`Model`)) + geom_bar(stat="identity",alpha=0.6) +
  scale_fill_manual(values = Colorbar) +
  theme(text = element_text(size=14),
        panel.background = element_rect(fill="white",color="black")) +
  geom_text(aes(label=R2,y=R2+0.02),position = position_dodge(0.9),vjust=0)
show(gr)


G = graph_from_adjacency_matrix(Z,mode='directed')
# corrplot(Z,is.corr = F,tl.col = "black",tl.pos = "d",tl.cex=0.15,method="color",cl.pos = "n",
#          col = c("white","orange"),order="AOE")
plot(G,layout=layout.auto,edge.arrow.size=0.5,vertex.color="#6E9ECE",
     vertex.frame.color="white",vertex.label.color="black",vertex.label.cex=1,
     vertex.size=15,edge.width=0.3,edge.color="black",layout=layout_with_fr)
dev.off()

draw_sam<-function(obs,re_name,shade)
{
  #T*Q
  df = data.frame(Value=as.vector(obs),
                  Response = factor(rep(re_name,each=nrow(obs)),levels = re_name),
                  Time = rep(1:nrow(obs),ncol(obs)))
  shade$Response = factor(shade$Response,levels = re_name)
  gr = ggplot(df,aes(x=Time,y=Value)) +
    geom_line(aes(color=Response))+
    scale_color_brewer(palette = "Dark2",guide="none") +
    facet_grid(Response~.,scales = "free")+
    geom_rect(data = shade,aes(x=NULL,y=NULL,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="blue",alpha=0.2)+
    #ylim(apply(obs,2,min),apply(obs,2,max))+
    theme(panel.background = element_rect(fill="white",color="black"),
          text = element_text(size=14),panel.border = element_blank())+
    guides(fill="none")
    show(gr)
}

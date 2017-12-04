library(mgcv)
library(ggplot2)
library(GGally)
library(MuMIn)

#各シミュレーションの試行回数
n=5

# Load the data:
file_1 <- read.csv("/Users/takimototomohiro/Desktop/analsys/170123/mean/p=0.1.csv",header=F,col.names=c("salience","sec","1","2","3","4","5","SD_1"))#red
file_2 <- read.csv("/Users/takimototomohiro/Desktop/analsys/170123/mean/p=0.01.csv",header=F,col.names=c("salience","sec","1","2","3","4","5","SD_2"))#blue
file_3 <- read.csv("/Users/takimototomohiro/Desktop/analsys/170123/mean/p=0.001.csv",header=F,col.names=c("salience","sec","1","2","3","4","5","SD_3"))#green
file_4 <- read.csv("/Users/takimototomohiro/Desktop/analsys/170123/mean/p=0.0001.csv",header=F,col.names=c("salience","sec","1","2","3","4","5","SD_4"))#black


#標準偏差表示のための配列整理
x <- 1
y <- 1
SD1 <- c(1:9991)
while (x <1000 ) {         # while ( 条件式 )
SD1[y:(y+9)] <- file_1$SD_1[x]
x <- x + 1             # 条件式が TRUE である限り式が繰り返される
y <- y + 10
}                        # 最初に条件式が FALSE ならば式は 1 回も実行されない
SD1[9991] <- file_1$SD_1[1000]

x <- 1
y <- 1
SD2 <- c(1:9991)
while (x <1000 ) {         # while ( 条件式 )
SD2[y:(y+9)] <- file_2$SD_2[x]
x <- x + 1             # 条件式が TRUE である限り式が繰り返される
y <- y + 10
}                        # 最初に条件式が FALSE ならば式は 1 回も実行されない
SD2[9991] <- file_2$SD_2[1000]

x <- 1
y <- 1
SD3 <- c(1:9991)
while (x <1000 ) {         # while ( 条件式 )
SD3[y:(y+9)] <- file_3$SD_3[x]
x <- x + 1             # 条件式が TRUE である限り式が繰り返される
y <- y + 10
}                        # 最初に条件式が FALSE ならば式は 1 回も実行されない
SD3[9991] <- file_3$SD_3[1000]

x <- 1
y <- 1
SD4 <- c(1:9991)
while (x <1000 ) {         # while ( 条件式 )
SD4[y:(y+9)] <- file_4$SD_4[x]
x <- x + 1             # 条件式が TRUE である限り式が繰り返される
y <- y + 10
}                        # 最初に条件式が FALSE ならば式は 1 回も実行されない
SD4[9991] <- file_4$SD_4[1000]


#一般化加法モデル

gam.model_1<- gam(salience ~ s(sec),sp=1,data=file_1,family=Gamma("log"))
gam.model_2<- gam(salience ~ s(sec),sp=1,data=file_2,family=Gamma("log"))
gam.model_3<- gam(salience ~ s(sec),sp=1,data=file_3,family=Gamma("log"))
gam.model_4<- gam(salience ~ s(sec),sp=1,data=file_4,family=Gamma("log"))



#予測用データの作成 新たな説明変数を作成
new_1<-data.frame(sec=seq(min(file_1$sec),max(file_1$sec),0.1),group="p=0.1")
new_2<-data.frame(sec=seq(min(file_2$sec),max(file_2$sec),0.1),group="p=0.01")
new_3<-data.frame(sec=seq(min(file_3$sec),max(file_3$sec),0.1),group="p=0.001")
new_4<-data.frame(sec=seq(min(file_4$sec),max(file_4$sec),0.1),group="p=0.0001")


#about 1
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_1<-predict(gam.model_1,new_1,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
#conf.lwr_1 <- gam.pred_1$fit-critval*gam.pred_1$se.fit
#conf.upr_1 <- gam.pred_1$fit + critval* gam.pred_1$se.fit
conf.lwr_1 <- gam.pred_1$fit-SD1
conf.upr_1 <- gam.pred_1$fit +SD1

#データフレームに格納
gam.pred_1 <- data.frame(new_1,as.data.frame(gam.pred_1),conf.lwr_1 = conf.lwr_1,conf.upr_1 = conf.upr_1,group="p=0.1" )

#about 2
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_2<-predict(gam.model_2,new_2,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
#conf.lwr_2 <- gam.pred_2$fit-critval*gam.pred_2$se.fit
#conf.upr_2 <- gam.pred_2$fit + critval* gam.pred_2$se.fit
conf.lwr_2 <- gam.pred_2$fit-SD2
conf.upr_2 <- gam.pred_2$fit + SD2
#データフレームに格納
gam.pred_2 <- data.frame(new_2,as.data.frame(gam.pred_2),conf.lwr_2 = conf.lwr_2,conf.upr_2 = conf.upr_2,group="p=0.01" )

#about 3
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_3<-predict(gam.model_3,new_3,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
#conf.lwr_3 <- gam.pred_3$fit-critval*gam.pred_3$se.fit
#conf.upr_3 <- gam.pred_3$fit + critval* gam.pred_3$se.fit
conf.lwr_3 <- gam.pred_3$fit-SD3
conf.upr_3 <- gam.pred_3$fit +SD3
#データフレームに格納
gam.pred_3 <- data.frame(new_3,as.data.frame(gam.pred_3),conf.lwr_3 = conf.lwr_3,conf.upr_3 = conf.upr_3,group="p=0.001" )

#about 4
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_4<-predict(gam.model_4,new_4,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
#conf.lwr_4 <- gam.pred_4$fit-critval*gam.pred_4$se.fit
#conf.upr_4 <- gam.pred_4$fit + critval* gam.pred_4$se.fit
conf.lwr_4 <- gam.pred_4$fit-SD4
conf.upr_4 <- gam.pred_4$fit +SD4
#データフレームに格納
gam.pred_4 <- data.frame(new_4,as.data.frame(gam.pred_4),conf.lwr_4 = conf.lwr_4,conf.upr_4 = conf.upr_4,group="p=0.0001" )


SD1t   <- data.frame(sec=gam.pred_1$sec,lwr=gam.pred_1$conf.lwr_1,upr=gam.pred_1$conf.upr_1)
SD2t   <- data.frame(sec=gam.pred_1$sec,lwr=gam.pred_2$conf.lwr_2,upr=gam.pred_2$conf.upr_2)
SD3t   <- data.frame(sec=gam.pred_1$sec,lwr=gam.pred_3$conf.lwr_3,upr=gam.pred_3$conf.upr_3)
SD4t   <- data.frame(sec=gam.pred_1$sec,lwr=gam.pred_4$conf.lwr_4,upr=gam.pred_4$conf.upr_4)

#予測結果の図示
g <- ggplot() +

  theme_set(theme_minimal(base_size=18,base_family="Helvetica"))+

#  geom_point(data=file_1,aes(x=sec,y=salience), size=0.5,colour=group) +
  geom_line(data=gam.pred_1,aes(x=sec, y=fit,colour=group), size=1.2)+
  geom_ribbon(data=SD1t,aes(x=sec, ymin=lwr, ymax=upr),alpha=0.2, fill="black")+
#  geom_point(data=file_1,aes(x=sec,y=SD_1), size=0.5,colour="black") +

#  geom_point(data=file_2,aes(x=sec,y=salience), size=0.5,colour=group) +
  geom_line(data=gam.pred_2,aes(x=sec, y=fit,colour=group), size=1.2)+
  geom_ribbon(data=SD2t,aes(x=sec, ymin=lwr, ymax=upr),alpha=0.2, fill="black")+
#  geom_ribbon(data=gam.pred_2,aes(x=file_1$sec, ymin=conf.lwr_2, ymax=conf.upr_2),alpha=0.2, fill="black")+

#  geom_point(data=file_3,aes(x=sec,y=salience), size=0.5,colour=group) +
  geom_line(data=gam.pred_3,aes(x=sec, y=fit,colour=group), size=1.2)+
  geom_ribbon(data=SD3t,aes(x=sec, ymin=lwr, ymax=upr),alpha=0.2, fill="black")+
#  geom_ribbon(data=gam.pred_3,aes(x=sec, ymin=conf.lwr_3, ymax=conf.upr_3),alpha=0.2, fill="black")+

#  geom_point(data=file_4,aes(x=sec,y=salience), size=0.5,colour=group) +
  geom_line(data=gam.pred_4,aes(x=sec, y=fit,colour=group), size=1.2)+
  geom_ribbon(data=SD4t,aes(x=sec, ymin=lwr, ymax=upr),alpha=0.2, fill="black")+
#  geom_ribbon(data=gam.pred_4,aes(x=sec, ymin=conf.lwr_4, ymax=conf.upr_4),alpha=0.2, fill="black")+


#setting of graf
  labs(x = "sec", y = "salience", size=2) +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  scale_color_manual(name="",values=c("p=0.1"="blue","p=0.01"="black","p=0.001"="green","p=0.0001"="red"))+
 # theme(legend.position = "top")
  guides(fill=FALSE)

g=g+xlim(0,1000)+ylim(0,15)

#plot
plot(g)

#Save as pdf
quartz.save("mean_sd.pdf",type="pdf")
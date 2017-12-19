library(mgcv)
library(ggplot2)
library(GGally)
library(MuMIn)
library(openxlsx)

#----------------------------------------------------------#
#  rbindCOorder(df1,df2,df3…)
# 列名無視して順序で結合。三つ以上のデータフレームにも対応
# df1の列名を採用する。
# a b c
# a c d
#　↓
# a b c
# a c d
#----------------------------------------------------------#
rbindCOrder <- function(...)  {
        n <- length(list(...))
        temp <- list(...)[[1]]
            names(temp)<-NA
        for (i in 2:n) {
            tmp<-list(...)[[i]]
            names(tmp)<-NA
            temp <- rbind(temp, tmp)
        }
            names(temp)<-names(list(...)[[1]])
        return(temp)
    }
#----------------------------------------------------------#

#title
title <- "170222"

#load the directory 
dir <- "~/babbling/created_data/"
filename <- "170222/"
csvfile <- "csv/"

#filename
id <- title
iterate <- "4000"
ploton <- "0"
IP <- "IP"
phase <-"notseparate"
network <- "lattice"
reward <- "negativereward"



#for GAM
# Load the data:
file_1 <- read.csv(paste(dir, filename,csvfile,"mean_p=0.05_NY_NSTDP.csv", sep = ""),header=F,col.names=c("salience","sec"))#red No_NSTD
file_2 <- read.csv(paste(dir, filename,csvfile,"mean_p=0.05_NY_STDP.csv", sep = ""),header=F,col.names=c("salience","sec"))#blue No_STDP
file_3 <- read.csv(paste(dir, filename,csvfile,"mean_p=0.05_Yoked_NSTDP.csv", sep = ""),header=F,col.names=c("salience","sec"))#green Sc_NSTD
file_4 <- read.csv(paste(dir, filename,csvfile,"mean_p=0.05_Yoked_STDP.csv", sep = ""),header=F,col.names=c("salience","sec"))#black Sc_STDP


gam.model_1<- gam(salience ~ s(sec),sp=1,data=file_1,family=Gamma("log"))
gam.model_2<- gam(salience ~ s(sec),sp=1,data=file_2,family=Gamma("log"))
gam.model_3<- gam(salience ~ s(sec),sp=1,data=file_3,family=Gamma("log"))
gam.model_4<- gam(salience ~ s(sec),sp=1,data=file_4,family=Gamma("log"))
#gam.model_5<- gam(salience ~ s(sec),sp=1,data=file_5,family=Gamma("log"))
#gam.model_6<- gam(salience ~ s(sec),sp=1,data=file_6,family=Gamma("log"))
#gam.model_7<- gam(salience ~ s(sec),sp=1,data=file_7,family=Gamma("log"))
#gam.model_8<- gam(salience ~ s(sec),sp=1,data=file_8,family=Gamma("log"))



#予測用データの作成 新たな説明変数を作成
new_1<-data.frame(sec=seq(min(file_1$sec),max(file_1$sec),0.1),group="Auditory feedback with STDP")
new_2<-data.frame(sec=seq(min(file_2$sec),max(file_2$sec),0.1),group="Auditory feedback without STDP")
new_3<-data.frame(sec=seq(min(file_3$sec),max(file_3$sec),0.1),group="Scramble feedback with STDP")
new_4<-data.frame(sec=seq(min(file_4$sec),max(file_4$sec),0.1),group="Scramble feedback without STDP")
#new_5<-data.frame(sec=seq(min(file_5$sec),max(file_5$sec),0.1),group="p=0.01")
#new_6<-data.frame(sec=seq(min(file_6$sec),max(file_6$sec),0.1),group="p=0.025")
#new_7<-data.frame(sec=seq(min(file_7$sec),max(file_7$sec),0.1),group="p=0.04")
#new_8<-data.frame(sec=seq(min(file_8$sec),max(file_8$sec),0.1),group="p=0.06")

#about 1
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_1<-predict(gam.model_1,new_1,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_1 <- gam.pred_1$fit-critval*gam.pred_1$se.fit
conf.upr_1 <- gam.pred_1$fit + critval* gam.pred_1$se.fit
#データフレームに格納
gam.pred_1 <- data.frame(new_1,as.data.frame(gam.pred_1),conf.lwr_1 = conf.lwr_1,conf.upr_1 = conf.upr_1,group="Auditory feedback with STDP" )

#about 2
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_2<-predict(gam.model_2,new_2,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_2 <- gam.pred_2$fit-critval*gam.pred_2$se.fit
conf.upr_2 <- gam.pred_2$fit + critval* gam.pred_2$se.fit
#データフレームに格納
gam.pred_2 <- data.frame(new_2,as.data.frame(gam.pred_2),conf.lwr_2 = conf.lwr_2,conf.upr_2 = conf.upr_2,group="Auditory feedback without STDP" )

#about 3
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_3<-predict(gam.model_3,new_3,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_3 <- gam.pred_3$fit-critval*gam.pred_3$se.fit
conf.upr_3 <- gam.pred_3$fit + critval* gam.pred_3$se.fit
#データフレームに格納
gam.pred_3 <- data.frame(new_3,as.data.frame(gam.pred_3),conf.lwr_3 = conf.lwr_3,conf.upr_3 = conf.upr_3,group="Scramble feedback with STDP" )

#about 4
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_4<-predict(gam.model_4,new_4,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_4 <- gam.pred_4$fit-critval*gam.pred_4$se.fit
conf.upr_4 <- gam.pred_4$fit + critval* gam.pred_4$se.fit
#データフレームに格納
gam.pred_4 <- data.frame(new_4,as.data.frame(gam.pred_4),conf.lwr_4 = conf.lwr_4,conf.upr_4 = conf.upr_4,group="Scramble feedback without STDP" )


#データフレームの間引き
step <- 1000  #1000毎
gam1 <- gam.pred_1[0:(length(gam.pred_1[,1])/step) *step +1, ]
gam2 <- gam.pred_2[0:(length(gam.pred_2[,1])/step) *step +1, ]
gam3 <- gam.pred_3[0:(length(gam.pred_3[,1])/step) *step +1, ]
gam4 <- gam.pred_4[0:(length(gam.pred_4[,1])/step) *step +1, ]


gamlast1 <- gam.pred_1[length(gam.pred_1[,1]),]
gamlast2 <- gam.pred_2[length(gam.pred_2[,1]),]
gamlast3 <- gam.pred_3[length(gam.pred_3[,1]),]
gamlast4 <- gam.pred_4[length(gam.pred_4[,1]),]

#データフレームの結合
gam_all=rbindCOrder(gam1,gam2,gam3,gam4,gamlast1,gamlast2,gamlast3,gamlast4)



#予測結果の図示


#g <- ggplot() +

 # theme_set(theme_minimal(base_size=18,base_family="Helvetica"))+

#  geom_point(data=file_1,aes(x=sec,y=salience), size=0.5,colour=group) +
#  geom_point(data=gam.pred_1,aes(x=sec, y=fit,group=group), size=1.2)+
#  geom_ribbon(data=gam.pred_1,aes(x=sec, ymin=conf.lwr_1, ymax=conf.upr_1),alpha=0.2, fill="black")+

#  geom_point(data=file_2,aes(x=sec,y=salience), size=0.5,colour=group) +
#  geom_line(data=gam.pred_2,aes(x=sec, y=fit,colour=group), size=1.2)+
#  geom_ribbon(data=gam.pred_2,aes(x=sec, ymin=conf.lwr_2, ymax=conf.upr_2),alpha=0.2, fill="black")+

#  geom_point(data=file_3,aes(x=sec,y=salience), size=0.5,colour=group) +
#  geom_line(data=gam.pred_3,aes(x=sec, y=fit,colour=group), size=1.2)+
#  geom_ribbon(data=gam.pred_3,aes(x=sec, ymin=conf.lwr_3, ymax=conf.upr_3),alpha=0.2, fill="black")+

#  geom_point(data=file_4,aes(x=sec,y=salience), size=0.5,colour=group) +
#  geom_line(data=gam.pred_4,aes(x=sec, y=fit,colour=group), size=1.2)+
#  geom_ribbon(data=gam.pred_4,aes(x=sec, ymin=conf.lwr_4, ymax=conf.upr_4),alpha=0.2, fill="black")+



#予測結果の図示(白黒)

g <- ggplot(gam_all,aes(x=sec,y=fit,group=group)) +

  geom_point(aes(shape = group), size=4)+
  geom_line(aes(linetype = group))+
  geom_errorbar(aes(ymax = conf.upr_1, ymin = conf.lwr_1), width = 30)+








#setting of graf
  theme_bw(base_size=30)+
  labs(x = "sec", y = "salience", size=2) +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  scale_color_manual(name="",values=c("Auditory feedback with STDP"="blue","Auditory feedback without STDP"="green","Scramble feedback with STDP"="red","Scramble feedback without STDP"="yellow"))+
 # theme(legend.position = "top")
  guides(fill=FALSE)

g=g+xlim(0,4000)+ylim(5,11)  +theme(legend.position = "right")

#plot
plot(g)

#Save as pdf

quartz.save(paste(dir, filename,title,".pdf", sep = ""),type="pdf")

#wirte csv

wb <- createWorkbook()
addWorksheet(wb, 'Normal+STDP')
addWorksheet(wb, 'Normal+NSTDP')
addWorksheet(wb, 'Scramble+STDP')
addWorksheet(wb, 'Scramble+NSTDP')




# withFilter = T の時は各列に Excel のフィルタが付きます。
writeData(wb, sheet = 'Normal+STDP', x = gam.pred_1, withFilter=F)
writeData(wb, sheet = 'Normal+NSTDP', x = gam.pred_2, withFilter=F)
writeData(wb, sheet = 'Scramble+STDP', x = gam.pred_3, withFilter=F)
writeData(wb, sheet = 'Scramble+NSTDP', x = gam.pred_4, withFilter=F)




# フォント指定をしない場合、Calibri が適用されます。
# "MS Pゴシック"や"ＭＳ　Ｐゴシック"を指定すると、xlsx ファイルが壊れました。
modifyBaseFont(wb, fontSize = 11, fontColour = "#000000", fontName = "MS PGothic")
# 上書きを許可しない場合、overwrite = F とする。
saveWorkbook(wb, "salhist_se.xlsx", overwrite = T)
openXL("~/Dropbox/osakauniversity/program/analsys/salhist_GAM/Yoked_vs_NY/2000/analsys/salhist_se.xlsx")

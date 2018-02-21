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


moving_average <- function(x, n){
  filter(x, rep(1,n)) / n
}


#title
title <- "180120_Sctime"

#load the directory 
dir <- "~/babbling/created_data/"
filename <- "180120_Sctime/"
csvfile <- "csv/"

#filename
id <- title
iterate <- "2000"
ploton <- "0"
IP <- "LiIP"    #LiIP or threIP or NoIP or tonic
phase <-"randSc"   #nseparate or separate or randSc
network <- "random"   #random or lattice
reward <- "normal"    #nega or normal
feedback <- "fft" #consonant or fft or no
iteratenum <- "1"
p <- "0.03"

#for GAM
# Load the data:
file_1 <- read.csv(paste(dir, filename,csvfile,id,"_",iterate,"_reinforce_100_4_No_",ploton,"_1_",p,"_0.3_NSTD_",IP,"_",phase,"_",network,"_",reward,"_",feedback,".csv", sep = ""),header=F,col.names=c("salience","sec","n_reward","DA_hist"))#red No_NSTD
file_2 <- read.csv(paste(dir, filename,csvfile,id,"_",iterate,"_reinforce_100_4_No_",ploton,"_1_",p,"_0.3_STDP_",IP,"_",phase,"_",network,"_",reward,"_",feedback,".csv", sep = ""),header=F,col.names=c("salience","sec","n_reward","DA_hist"))#blue No_STDP
file_3 <- read.csv(paste(dir, filename,csvfile,id,"_",iterate,"_reinforce_100_4_Sc_",ploton,"_1_",p,"_0.3_NSTD_",IP,"_",phase,"_",network,"_",reward,"_",feedback,".csv", sep = ""),header=F,col.names=c("salience","sec","n_reward","DA_hist"))#green Sc_NSTD
file_4 <- read.csv(paste(dir, filename,csvfile,id,"_",iterate,"_reinforce_100_4_Sc_",ploton,"_1_",p,"_0.3_STDP_",IP,"_",phase,"_",network,"_",reward,"_",feedback,".csv", sep = ""),header=F,col.names=c("salience","sec","n_reward","DA_hist"))#black Sc_STDP
file_5 <- read.csv(paste("/Users/takimoto/babbling/created_data/180120_SctimeNone/csv/180120_Sctime_2000_reinforce_100_4_No_0_1_0.03_0.3_STDP_LiIP_randSc_random_normal_none.csv", sep = ""),header=F,col.names=c("salience","sec","n_reward","DA_hist"))#slategray2 Nofeedback_STDP
file_6 <- read.csv(paste("/Users/takimoto/babbling/created_data/180120_SctimeNone/csv/180120_Sctime_2000_reinforce_100_4_No_0_1_0.03_0.3_NSTD_LiIP_randSc_random_normal_none.csv", sep = ""),header=F,col.names=c("salience","sec","n_reward","DA_hist"))#purple Nofeedback_NSTD






gam.model_1<- gam(salience ~ s(sec),sp=1,data=file_1,family=Gamma("log"))
gam.model_2<- gam(salience ~ s(sec),sp=1,data=file_2,family=Gamma("log"))
gam.model_3<- gam(salience ~ s(sec),sp=1,data=file_3,family=Gamma("log"))
gam.model_4<- gam(salience ~ s(sec),sp=1,data=file_4,family=Gamma("log"))
gam.model_5<- gam(salience ~ s(sec),sp=1,data=file_5,family=Gamma("log"))
gam.model_6<- gam(salience ~ s(sec),sp=1,data=file_6,family=Gamma("log"))



#予測用データの作成 新たな説明変数を作成
new_1<-data.frame(sec=seq(min(file_1$sec),max(file_1$sec),0.1),group="Auditory feedback without STDP")
new_2<-data.frame(sec=seq(min(file_2$sec),max(file_2$sec),0.1),group="Auditory feedback with STDP")
new_3<-data.frame(sec=seq(min(file_3$sec),max(file_3$sec),0.1),group="Scramble feedback without STDP")
new_4<-data.frame(sec=seq(min(file_4$sec),max(file_4$sec),0.1),group="Scramble feedback with STDP")
new_5<-data.frame(sec=seq(min(file_5$sec),max(file_5$sec),0.1),group="Nofeedback with STDP")
new_6<-data.frame(sec=seq(min(file_6$sec),max(file_6$sec),0.1),group="Nofeedback without STDP")


#about 1
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_1<-predict(gam.model_1,new_1,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_1 <- gam.pred_1$fit-critval*gam.pred_1$se.fit
conf.upr_1 <- gam.pred_1$fit + critval* gam.pred_1$se.fit
#データフレームに格納
gam.pred_1 <- data.frame(new_1,as.data.frame(gam.pred_1),conf.lwr_1 = conf.lwr_1,conf.upr_1 = conf.upr_1,group="Auditory feedback without STDP" )

#about 2
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_2<-predict(gam.model_2,new_2,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_2 <- gam.pred_2$fit-critval*gam.pred_2$se.fit
conf.upr_2 <- gam.pred_2$fit + critval* gam.pred_2$se.fit
#データフレームに格納
gam.pred_2 <- data.frame(new_2,as.data.frame(gam.pred_2),conf.lwr_2 = conf.lwr_2,conf.upr_2 = conf.upr_2,group="Auditory feedback with STDP" )

#about 3
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_3<-predict(gam.model_3,new_3,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_3 <- gam.pred_3$fit-critval*gam.pred_3$se.fit
conf.upr_3 <- gam.pred_3$fit + critval* gam.pred_3$se.fit
#データフレームに格納
gam.pred_3 <- data.frame(new_3,as.data.frame(gam.pred_3),conf.lwr_3 = conf.lwr_3,conf.upr_3 = conf.upr_3,group="Scramble feedback without STDP" )

#about 4
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_4<-predict(gam.model_4,new_4,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_4 <- gam.pred_4$fit-critval*gam.pred_4$se.fit
conf.upr_4 <- gam.pred_4$fit + critval* gam.pred_4$se.fit
#データフレームに格納
gam.pred_4 <- data.frame(new_4,as.data.frame(gam.pred_4),conf.lwr_4 = conf.lwr_4,conf.upr_4 = conf.upr_4,group="Scramble feedback with STDP" )


#about 5
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_5<-predict(gam.model_5,new_5,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_5 <- gam.pred_5$fit-critval*gam.pred_5$se.fit
conf.upr_5 <- gam.pred_5$fit + critval* gam.pred_5$se.fit
#データフレームに格納
gam.pred_5 <- data.frame(new_5,as.data.frame(gam.pred_5),conf.lwr_5 = conf.lwr_5,conf.upr_5 = conf.upr_5,group="Nofeedback with STDP" )

#about 6
#予測　目標変数の予測type="response",予測子空間"link", se.fit=True 標準誤差出力、予測値fit,標準誤差se.fit
gam.pred_6<-predict(gam.model_6,new_6,se.fit=TRUE,type="response")
#95%信頼度区間を得る
critval=qnorm(0.975,0,1)
conf.lwr_6 <- gam.pred_6$fit-critval*gam.pred_6$se.fit
conf.upr_6 <- gam.pred_6$fit + critval* gam.pred_6$se.fit
#データフレームに格納
gam.pred_6 <- data.frame(new_6,as.data.frame(gam.pred_6),conf.lwr_6 = conf.lwr_6,conf.upr_6 = conf.upr_6,group="Nofeedback without STDP" )





#データフレームの間引き
step <- 1000  #1000毎
gam1 <- gam.pred_1[0:(length(gam.pred_1[,1])/step) *step +1, ]
gam2 <- gam.pred_2[0:(length(gam.pred_2[,1])/step) *step +1, ]
gam3 <- gam.pred_3[0:(length(gam.pred_3[,1])/step) *step +1, ]
gam4 <- gam.pred_4[0:(length(gam.pred_4[,1])/step) *step +1, ]
gam5 <- gam.pred_5[0:(length(gam.pred_5[,1])/step) *step +1, ]
gam6 <- gam.pred_6[0:(length(gam.pred_6[,1])/step) *step +1, ]

gamlast1 <- gam.pred_1[length(gam.pred_1[,1]),]
gamlast2 <- gam.pred_2[length(gam.pred_2[,1]),]
gamlast3 <- gam.pred_3[length(gam.pred_3[,1]),]
gamlast4 <- gam.pred_4[length(gam.pred_4[,1]),]
gamlast5 <- gam.pred_5[length(gam.pred_5[,1]),]
gamlast6 <- gam.pred_6[length(gam.pred_6[,1]),]



#データフレームの結合(GAM)
gam_all=rbindCOrder(gam1,gam2,gam3,gam4,gam5,gam6,gamlast1,gamlast2,gamlast3,gamlast4,gamlast5,gamlast6)



#予測結果の図示(白黒)

#GAM
g <- ggplot(gam_all,aes(x=sec,y=fit,group=group,colour = group)) +
#ggplot(gam_all,aes(x=sec,y=fit,group=group)) +
  geom_point(aes(shape = group), size=4)+
  geom_line(aes(linetype = group))+
  geom_errorbar(aes(ymax = conf.upr_1, ymin = conf.lwr_1), width = 30)+








#setting of graf
  theme_bw(base_size=20)+
  labs(x = "sec", y = "salience", size=20) +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  scale_color_manual(name="",values=c("Auditory feedback with STDP"="red","Auditory feedback without STDP"="green","Scramble feedback with STDP"="blue","Scramble feedback without STDP"="yellow","Nofeedback with STDP"="slategray","Nofeedback without STDP"="purple"))+
  theme(legend.position = "top")
  guides(fill=FALSE)

g=g+
#xlim(0,iterate)+
ylim(4,12)  +theme(legend.position = "right")



#plot
plot(g)

#Save as pdf
quartz.save(paste(dir, filename,title,iteratenum,"formaster.pdf", sep = ""),type="pdf")





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
saveWorkbook(wb,paste(dir, filename,title,iteratenum,".xlsx", sep = "") ,overwrite = T)

#openXL("~/Dropbox/osakauniversity/program/analsys/salhist_GAM/Yoked_vs_NY/2000/analsys/salhist_se.xlsx")



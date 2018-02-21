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



#for negative reward
sali1 <-data.frame(sec = file_1$sec,salience = file_1$salience,group="Auditory feedback without STDP")
sali2 <-data.frame(sec = file_2$sec,reward = file_2$salience,group="Auditory feedback with STDP")
sali3 <-data.frame(sec = file_3$sec,reward = file_3$salience,group="Scramble feedback without STDP")
sali4 <-data.frame(sec = file_4$sec,reward = file_4$salience,group="Scramble feedback with STDP")



#データフレームの結合(GAM)
sali_all <- rbindCOrder(sali1,sali2,sali3,sali4)

#moving average
sali_all$salience <- moving_average(sali_all$salience,20)
#DA_all$DA_hist <- moving_average(DA_all$DA_hist,50)



#予測結果の図示(白黒)

#GAM
g <- ggplot(sali_all,aes(x=sec,y=salience,group=group,colour = group)) +
#ggplot(gam_all,aes(x=sec,y=fit,group=group)) +
  geom_line(aes(linetype = group))+





#setting of graf
  theme_bw(base_size=20)+
  labs(x = "sec", y = "salience", size=20) +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20))+
  scale_color_manual(name="",values=c("Auditory feedback with STDP"="red","Auditory feedback without STDP"="green","Scramble feedback with STDP"="blue","Scramble feedback without STDP"="yellow"))+
  theme(legend.position = "top")
  guides(fill=FALSE)

g=g+
xlim(20,1980)+
ylim(4,12)  +theme(legend.position = "right")


#plot
plot(g)

#Save as pdf
quartz.save(paste(dir, filename,title,iteratenum,"movingaverage.pdf", sep = ""),type="pdf")





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



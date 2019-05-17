## Aim: Simulate the measurement process
## Status/belongs to: Isotope lab/ PaleoEis Group



library(tidyverse)

#Define the machine characteristics (to be determined from experiments)
SD.MACHINE.d18O     <- 0.025 #standard deviation of random measurement error assumed to be iid and gaussian [Value from PICARRO]
DRIFT.MACHINE.d18O  <- 0.2 / 24 / 30 #linear Drift in permille per hour, estimate by hand from the experiment [Value from PICARRO would be 0.2 permille per 24h = 0.2/24]
MEMORY.MACHINE.d18O <- 0.1 #Memory coefficient of the machine

SD.MACHINE.dD     <- 0.025*8 #standard deviation of random measurement error assumed to be iid and gaussian [Value from PICARRO]
DRIFT.MACHINE.dD  <- 0.2 / 24 / 30 *8  #linear Drift in permille per hour, estimate by hand from the experiment [Value from PICARRO would be 0.2 permille d18O per 24h = 0.2/24]
MEMORY.MACHINE.dD <- 0.1  #Memory coefficient of the machine

DRIFT.VIAL.d18O     <- 0.02 / 100 #Linear Drift in permille per hour and per hole in the Septum
DRIFT.VIAL.dD     <- 0.02 / 100 * 8 #Linear Drift in permille per hour and per hole in the Septum

#Holes in the vial are modeled that an untouched vial will first get one hole... and than just an incremental increase of the hole size
DHOLE.VIAL     <- 0. #Not every injection will create a new hole; this coefficient determines the increase of the Septum holes per injection
DWATER.INJECTION <- 1.9/1000 #Injection volume in ml

TIME.INJECTION <- 9 #Time per injection in minutes



calib.dD <- function(x)
{
    OFFSET <- 3
    SLOPE <- 1.1
    return(x*SLOPE + OFFSET)
}

calib.d18O <- function(x)
{
    OFFSET <- 3
    SLOPE <- 1.1
    return(x*SLOPE + OFFSET)
}





##' @title Simple simulation of assumed PICARRO measurement process
##' @param vial dataframe with columns d18O (in permille),dD (in permille), water.volume (in ml),holes (#holes in the septum), Identifier 
##' @param sequence vector giving the sequence of the measurements
##' @param d18O.machine Initial d18O (permille) in the machine... needed for the drift
##' @return data.frame  with columns "Line","Analysis","Port","d.18_16.Mean"
##' @author Thomas Laepple 
simulate <- function(vial, sequence,  d18O.machine = -37.5, dD.machine = -300)
{

    
    ##Initial values
    runtime.minutes <- 0
    i.injection <- 1
    i.Analysis <- 0
    i.vial <- 99e99
    
    ##Predefine an empty dataframe as output...
    output <- matrix(NA,length(sequence),6)
    colnames(output) <- c("Line","Analysis","Port","d.18_16.Mean","d.D_H.Mean","Identifier.1")
    output <- as.data.frame(output)

    ##Run the simulation
    ##i.step = injection step
    for (i.step in 1:length(sequence))
    {
        runtime.minutes = runtime.minutes + TIME.INJECTION


        ##Every timestep, all vials drift proportional to the number of holes and the time passed (in hours)
        vial$d18O<-vial$d18O + vial$holes*DRIFT.VIAL.d18O*TIME.INJECTION/60
        vial$dD<-vial$dD + vial$holes*DRIFT.VIAL.dD*TIME.INJECTION/60
        

        if (i.vial != sequence[i.step])
        {
            i.injection <- 1
            i.Analysis <- i.Analysis+1;
            } else i.injection <- i.injection + 1;

        ##Choose the vial for sampling
        i.vial <- sequence[i.step]
        
        ##remove water from the vial 
        vial$water.volume[i.vial] <- vial$water.volume[i.vial]-DWATER.INJECTION


        ##Put / increase hole in the septum
        if (vial$holes[i.vial]==0)
        {
            vial$holes[i.vial]=1
        } else vial$holes[i.vial]=vial$holes[i.vial]+DHOLE.VIAL

       
        
        ##perform the measurement of d18O
        measured.d18O <- vial$d18O[i.vial] + rnorm(1,sd=SD.MACHINE.d18O) + DRIFT.MACHINE.d18O*(runtime.minutes/60) - (vial$d18O[i.vial]-d18O.machine) * MEMORY.MACHINE.d18O
        measured.dD <- vial$dD[i.vial] + rnorm(1,sd=SD.MACHINE.dD) + DRIFT.MACHINE.dD*(runtime.minutes/60) - (vial$dD[i.vial]-dD.machine) * MEMORY.MACHINE.dD
        

        d18O.machine <- measured.d18O
        dD.machine <- measured.dD

        ##Save the output
        output$Line[i.step]<-i.step
        output$Inj.Nr[i.step]<-i.injection
        output$Identifier.1[i.step]<-levels(vial$identifier[i.vial])[vial$identifier[i.vial]]
        
        output$Analysis[i.step]<-i.Analysis
        output$Port[i.step]<-paste(i.vial)
        output$d.18_16.Mean[i.step]<-calib.d18O(measured.d18O)
        output$d.D_H.Mean[i.step]<-calib.dD(measured.dD)
    }

return(output)
}


#### Simulating 

#Define what is initially in the vials; here we use 6 vials with standards
vial<-data.frame(d18O=c(-38.8,-38.8,-38.8,-38.8,-38.8,-38.8),
                 dD=c(-38.8,-38.8,-38.8,-38.8,-38.8,-38.8)*8,
                 water.volume=c(1.6,1.6,1.6,1.6,1.6,1.6),holes=c(0,0,0,0,0,0),identifier=c(rep("TD1",6)))

##Define the sequence of the measurements... blocks are just used to simplify the construction; block1 is the repeating 75,75,10 injections; block2 three untouched vials at the end to separate the instrument from the vial drift

block1 <-  c(rep(1,75),rep(2,75),rep(3,10))
block2 <-  c(rep(4,75),rep(5,75),rep(6,75))
sequence <- c(block1,block1,block1,block1,block1,block2)

#Arbitary naming of blocks that is later used to average the blocks
names <- (c(block1,block1+3,block1+6,block1+9,block1+12,block2+12))

experiments <- data.frame(Replicate = 1:4)
runs <- experiments %>% group_by(Replicate) %>% do({simulate(vial,sequence)})


##Basic plots and analyses using tidyR / ggplot

ggplot(runs, aes(x=Line,y=d.18_16.Mean,colour=Port)) +  geom_point()  + ggtitle("all measurements") + facet_wrap(~Replicate)
ggplot(runs, aes(x=Line,y=d.D_H.Mean,colour=Port)) +  geom_point()  + ggtitle("all measurements") + facet_wrap(~Replicate)


meanvals <- runs %>% group_by(Analysis,Replicate) %>% summarize(m=mean(d.18_16.Mean),n=n(),se=sd(d.18_16.Mean)/sqrt(n),Port=Port[1])

meanvals %>% ggplot( aes(x=Analysis,y=m,colour=Port)) +  geom_point() + geom_errorbar(aes(ymin=m-se,ymax=m+se)) +  ggtitle("mean values and standard errors")  + facet_wrap(~Replicate)


######

#Define what is initially in the vials; here we use 6 vials with standards
vial<-data.frame(d18O=c(-0.76,-19.9,-33.83,-50.22,-42.44,-30,-40),
                 dD=c(-5.6,-153.2,-266.9,-392.5,-341.3,-200,-300),
                 water.volume=c(1.6,1.6,1.6,1.6,1.6,1.6,1.6),holes=c(0,0,0,0,0,0,0),identifier=c("HGL-1","NZE","TD1","JASE","DML","PROBE1","PROBE2"))

##Define the sequence of the measurements... blocks are just used to simplify the construction; block1 is the repeating 75,75,10 injections; block2 three untouched vials at the end to separate the instrument from the vial drift

block1 <-   c(1,3,2,4)
block2 <-   c(2,5,4)
block3 <-   c(3,2,4)
sampleblock <- c(6,7)

sequence <- c(rep(block1,each=12), rep(rep(sampleblock,each=3),4),rep(block2,each=3),rep(rep(sampleblock,each=3),4),rep(block3,each=3))


experiments <- data.frame(Replicate = 1:4)
runs <- experiments %>% group_by(Replicate) %>% do({simulate(vial,sequence)})

temp<-  simulate(vial,sequence)
colnames(temp)<-c("Line","Analysis","Port","d(18_16)Mean","d(D_H)Mean"  ,"Identifier 1" ,"Inj Nr")
write.csv(file="SIMULATE_IsoWater_20190228_144452.csv",temp,row.names=FALSE)




##Basic plots and analyses using tidyR / ggplot

ggplot(runs, aes(x=Line,y=d.18_16.Mean,colour=Port)) +  geom_point()  + ggtitle("all measurements") + facet_wrap(~Replicate)

ggplot(runs, aes(x=Line,y=d.D_H.Mean,colour=Port)) +  geom_point()  + ggtitle("all measurements") + facet_wrap(~Replicate)


meanvals <- runs %>% group_by(Analysis,Replicate) %>% summarize(m=mean(d.18_16.Mean),n=n(),se=sd(d.18_16.Mean)/sqrt(n),Port=Port[1])

meanvals %>% ggplot( aes(x=Analysis,y=m,colour=Port)) +  geom_point() + geom_errorbar(aes(ymin=m-se,ymax=m+se)) +  ggtitle("mean values and standard errors")  + facet_wrap(~Replicate)






#### Experiment 2:


#Define what is initially in the vials; here we use 6 vials with standards
vial<-data.frame(d18O=c(-40,-30,-20,-35,-32),water.volume=c(1.6,1.6,1.6,1.6,1.6),holes=c(0,0,0,0,0))

##Define the sequence of the measurements... blocks are just used to simplify the construction; block1 is the repeating 75,75,10 injections; block2 three untouched vials at the end to separate the instrument from the vial drift

calibBlock <-  c(rep(1,10),rep(2,10),rep(3,10))
middleBlock <-     c(rep(4,6),rep(5,6))
sequence <- c(calibBlock,rep(middleBlock,20),calibBlock)

#Arbitary naming of blocks that is later used to average the blocks
names.vial <- seq(sequence)

experiments <- data.frame(Replicate = 1:4)

runs <- experiments %>% group_by(Replicate) %>% do({simulate(vial,sequence,names)})


##Basic plots and analyses using tidyR / ggplot

ggplot(runs, aes(x=Line,y=d.18_16.Mean,colour=Port)) +  geom_point()  + ggtitle("all measurements") + facet_wrap(~Replicate)
meanvals <- runs %>% group_by(Analysis,Replicate) %>% summarize(m=mean(d.18_16.Mean),n=n(),se=sd(d.18_16.Mean)/sqrt(n),Port=Port[1])
meanvals %>% ggplot( aes(x=Analysis,y=m,colour=Port)) +  geom_point() + geom_errorbar(aes(ymin=m-se,ymax=m+se)) +  ggtitle("mean values and standard errors")  + facet_wrap(~Replicate)










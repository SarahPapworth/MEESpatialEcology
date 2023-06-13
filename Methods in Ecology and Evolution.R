#R code to model NSD, BRB and RUF. 

#Start with a tab delimited text file. Coordinates are measured in metres, and the file should have columns the with following data for each location, in chronological order within an individual and trip. 
#Date - date on which location recorded dd/mm/yyyy
#Time - time at which location recorded hh:mm:ss. 
#IND: - unique reference code for each individual
#Latitude - Latitudinal reference in metres
#Longitude - Longitudinal reference in metres
#Trip - unique reference code for each trip

#Load the required packages. 
#Before you use a package for the first time you need to download it:
#install.package("")
#The exception to this is the package ruf
#before you use it the first time, remove the # from the following line
#install.packages("ruf",repos="http://www.csde.washington.edu/~handcock")
#then send to R to load the package
library(trip)
library(stringr)
library(adehabitatHR)
library(nlme)
library(lattice)
library(gmodels)
library(spatstat)
library(ruf)
library(maptools)

#example dataset is in this github repository
#open the file in R - navigate to where the file is saved.
d1<-read.table(file.choose(),header=T)

#To start, change date and time to the class POSIXct, so that the package adehabitat can create a trajectory of your locations.
datetime <- as.POSIXct(paste(d1$Date,d1$Time),
                       format = "%d/%m/%Y %H:%M:%S",
                       "GMT") #reference time: If your locations are recorded in the file using local
#time, use "GMT". If they are global time, specify the hour band location, e.g.
# "America/Lima"


##make a data.frame of the variables latitude and longitude. Here the raw values are divided
#by 1000 so that trajectories are calculated using km as the unit of measurement
coord<-data.frame((d1$Latitude),(d1$Longitude))
# make ltraj: a trajectory of all the relocations
d2<-as.ltraj(coord,datetime,
             d1$IND,        #separate your data by individual.
             burst=d1$Trip, #burst is used to creat subdivisions within an individual.
             typeII=TRUE)       #typeII can be TRUE: radio-track data, or FALSE: not time
#recorded, such as tracks in the snow
summary(d2)
#you should now be able to see a summary of the trajectory data,
#with one line for each burst

#you can now make your trajectory regular, as radio tracks tend to lose
#a few seconds / minutes with each relocation
#firstly add "NA" for each missing location in your trajectory
d3<-setNA(d2,
          as.POSIXct("2009-10-21 16:30:30"), #any time before earliest timedate in huntGPS
          60,            #stating there should be a location every 60th time unit
          tol=30,        #how many time units to search each side of expected location
          units="sec")   #specifying the time units

#you can now make your trajectory regular
#firstly create a reference start time
refda <- strptime("00:00:30", "%H:%M:%S")   #all relocations should be altered
#to occur at 30 seconds past each minute
d4<-sett0(d3, refda,
          60,                         #stating the interval at which relocations should be
          correction.xy =c("none"),   #if "cs" performs location correction based on the
          #assumption the individual moves at a constant speed
          tol=30,       #how many time units to search either side of an expected location
          units = "sec")  #specifying the time units

#to view your regular trajectory of points with NA's
summary(d4)
#now calculating NSD for each point
datansd<-NULL
for(n in 1:length(summary(d4)[,1])) #stating that NSD should be
  #calculated separately for each burst
{
  nsdall<-d4[[n]][,8]             #extracting the NSD for each location
  nsdtimeall<-d4[[n]][,3]         #extracting the time for each location
  nsdtimestartzero<-d4[[n]][,3]-d4[[n]][1,3]
  #extracting the time since trip start for each location
  nsdid<-rep(as.vector(summary(d4)[n,1]),
             length.out=summary(d4)[n,3])
  #extracting the individual associated with each location
  nsdtrip<-rep(as.vector(summary(d4)[n,2]),length.out=summary(d4)[n,3])
  #extracting the trip associated with each location
  datansd1<-data.frame(nsdall,nsdtimeall,nsdtimestartzero,nsdid,nsdtrip)
  #joining all these variables together in a data frame
  datansd<-rbind(datansd,datansd1)
  #joining all the data frames together
}
datansd$zero1<-as.numeric(unclass(datansd$nsdtimestartzero))
# making seconds since trip start numeric
datansd$zerostart<-datansd$zero1/60
#changing the time since trip start from seconds to minutes
datansd$minslitr2<-as.numeric(strftime(as.POSIXlt(datansd$nsdtimeall),
                                       format="%M"))
#making a vector of the hour of the day a location occured
datansd$hdaylitr2<-as.numeric(strftime(as.POSIXlt(datansd$nsdtimeall),
                                       format="%H"))
#making a vector of the minute in an hour a location occured
datansd$minsday<-((datansd$hdaylitr2*60)+datansd$minslitr2)
#calculating the minute in the day a location occured

summary(datansd)

#To select and name the hunting trips
hunt2<-c(5,6,8,10,11,12,14,15,19,22,24,26,27,28,30,31)
#a list of all numbers of all the hunting trips
datansd$hunt<-match(datansd$nsdtrip,hunt2,nomatch=0)
#assigning the value 0 to all non-hunting trips
datansd$hunt[datansd$hunt > 1] <- 1  #assigning the value 1 to all hunting trips
datansd$HUNT<-as.factor(datansd$hunt) #making it a factor
datansd1<-na.omit(datansd)            #remove NA's
datansd1$coordinates<-coord           #add the coordinates for each point
#you now have the dataframe you need (datansd1) to start analysis
#if the computer is slow, you can remove all the data sets you don't need
#to help it speed up
rm(d1)
rm(d2)
rm(d3)
rm(d4)
rm(coord)
rm(datetime)
rm(nsdid)
rm(nsdtimeall)
rm(nsdtimestartzero)
rm(nsdtrip)
rm(refda)


#NSD

#Now you can start modelling NSD using nlme. The book:
#Pinheiro and Bates (2004) Mixed-effects models in S and S-Plus.
#Springer Science: New York, can help, as will
#Bunnefeld at al (2011) A model-driven approach to quantify migration patterns:
#individual, regional and yearly differences.
#Journal of Animal Ecology 80: 466 - 476

#first model the data without random effects using nls, a least squares method
#this will help identify parameter estimates for use with nlme
m1<-nls(nsdall ~  asym /(1+exp((xmidA-zerostart)/scale1)) +
          (-asym / (1 + exp((xmidB-zerostart)/scale2))), #this part defines eqn 1
        start = c(asym=40000,xmidA=10,xmidB=30,scale1=4,scale2=4)
        #these are the starting values for each parameter of the equation
        ,data=na.omit(datansd1))   #this is the data
summary(m1)        #this will print a summary of the converged model
#graphical exploration of the data will help you find sensible starting values
#for each of the parameters asym, xmidA, xmidB, scale1 and scale2.
#to graph nsd again time, use:
xyplot(nsdall~zerostart|nsdtrip,data=datansd1)
#Alternately, you can start with a single individual such as "SA" below, and
#gradually add one individual at a time.
m2<-nls(nsdall ~  asym /(1+exp((xmidA-zerostart)/scale1)) +
          (-asym / (1 + exp((xmidB-zerostart)/scale2))),
        start = c(asym=40000,xmidA=10,xmidB=30,scale1=4,scale2=4)
        ,data=na.omit(datansd1[datansd1$nsdid=="A",]))   #to specify only one individual
summary(m2)
#try various starting values - the model will only converge if the values are
#sufficiently close to the modelled values
#now try and model the data including random effects
#start with no variation in the explanatory variable
m3<-nlme(nsdall ~  asym /(1+exp((xmidA-zerostart)/scale1)) +
           (-asym /(1 + exp((xmidB-zerostart)/scale2))), #the equation
         fixed = list(asym+xmidA+xmidB+scale1+scale2~1),  #fixed effects
         random= asym ~ 1|nsdid, #random effects: asym varies between individuals
         start = c(asym=40000,xmidA=10,xmidB=30,scale1=4,scale2=4)
         #starting vlaues for the parameters in the equation
         ,data=na.omit(datansd1))      #the data
print(AIC(m3))           #this will print the AIC of the converged model
#you can change the random effect structure
m4<-nlme(nsdall ~  asym /(1+exp((xmidA-zerostart)/scale1)) +
           (-asym /(1 + exp((xmidB-zerostart)/scale2))),
         fixed = list(asym+xmidA+xmidB+scale1+scale2~1),
         random= asym ~ 1|nsdid/nsdtrip,        #random effects: asym varies between
         #individuals, and also between trips within a single individual
         start = c(asym=38000,xmidA=9,xmidB=30,scale1=3,scale2=3)
         ,data=na.omit(datansd1))
print(AIC(m4))

#When you have the best random effects structure, you can model the data with
#differences between your groups
m5<-nlme(nsdall ~  asym /(1+exp((xmidA-zerostart)/scale1)) +
           (-asym / (1 + exp((xmidB-zerostart)/scale2))),
         fixed = list(asym+xmidA+xmidB+scale1+scale2~HUNT), #just change this to say HUNT
         random= asym ~ 1|nsdid/nsdtrip,
         start = c(asym=40000,20000,xmidA=10,0,xmidB=30,0,scale1=3,0, scale2=3,0)
         #and remember to add the extra parameters here: the first value for each
         #parameter is the expected value for "0", and the second value is the
         #difference between group "1" and group "0"
         ,data=na.omit(datansd1))
print(AIC(m5))
#now show a summary of the best model
summary(m5)
#you can show the fitted values
fitted(m5)
#normal probability plots
qqnorm(m5)
#the residuals
plot(m5)
#and estimate the 95% confidence intervals for the parameter estimates
#first make a matrix of all the posible constrasts: the matrix below is
#applicable if you are comparing two groups
matrix.contrasts<- rbind(c(1,0,0,0,0,0,0,0,0,0),
                         c(1,1,0,0,0,0,0,0,0,0),
                         c(0,0,1,0,0,0,0,0,0,0),
                         c(0,0,1,1,0,0,0,0,0,0),
                         c(0,0,0,0,1,0,0,0,0,0),
                         c(0,0,0,0,1,1,0,0,0,0),
                         c(0,0,0,0,0,0,1,0,0,0),
                         c(0,0,0,0,0,0,1,1,0,0),
                         c(0,0,0,0,0,0,0,0,1,0),
                         c(0,0,0,0,0,0,0,0,1,1))
#estimate the 95% confidence intervals
estimable(m5, matrix.contrasts, conf.int=0.95) #you can change this to
#estimate different confidence intervals
#you can also make a graph like figure 3
#create a window that is divided in two (one for hunting and one for
#non-hunting trips
par(mfrow=c(2,1))
#plot the best model
#this uses the parameters from the model to predict the curve for NSD
datansd1$pred<-predict(m5,level=0) #remember to put in the correct model here
myPanel <- function(x,y, ...){
  panel.xyplot(x,y, ...)
  dotArgs <- list(...)
  # select the appropriate rows of data and predict and then order them
  predY <- datansd1$pred[dotArgs$subscripts]
  predX <- datansd1$zerostart[dotArgs$subscripts]
  ord <- order(predX)
  predX <-  predX[ord]
  predY <-  predY[ord]
  # add as a panel line
  panel.lines(predX, predY, col='black', type='l',lwd=2)
}
#now plot the data with the predicted curve
xyplot(nsdall ~ zerostart|HUNT, data=datansd1,
       col="grey",    #color for the observed locations
       type='b',      # 'b' shows the locations as dots, with a line connecting
       #successive locations. Can also be 'p' for just the locations, or 'l' for just
       #the line between locations
       ylab=expression(paste('Net squared displacement ',' ', (km^2))), #y axis label
       xlab="Minutes after trip start",                                 #x axis label
       group=nsdtrip,            #grouping factor  - changed from nsdTRIP: important??
       panel=myPanel,            #predicted values from above
       strip=strip.custom(bg="grey", factor.levels=c('Non-hunting trips (n=17)',
                                                     'Hunting trips (n=19)'  )),  #to create a strip at the top to label each group
       scales=list(x=list(alternating=1,
                          at = c(0,10,20,30,40,50,60)),tck=-1,        #locations of marks on the x axis
                   y=list(alternating=1,
                          at=c(0,20000,40000,60000,80000,100000,120000,140000,160000),tck=-1) #locations
                   #of marks on the y axis
       ))

#to select the relevant data identified using NSD
#Group 1: non hunting trips
nothunt<-datansd1[datansd1$HUNT=="0",] #select the non hunting data
nothunt1<-na.omit(nothunt)   #remove the NA's generated by removing hunting data
nothunt1$include[nothunt1$zerostart > 18] <- 1    #select all locations where
#time after trip start is greater than 60
nothunt2<-na.omit(nothunt1)                       #remove the NA's generated
nothunt2$include1[nothunt2$zerostart < 24] <- 1  #select all the locations
#where time after time start is smaller than 265
nothunt3<-na.omit(nothunt2)                       #remove the NA's generated
#Group 2: hunting trips
huntdata<-datansd1[datansd1$HUNT=="1",]             #select the hunting data
huntdata1<-na.omit(huntdata)                      #remove the NA's generated by
#removing non hunting data
huntdata1$include[huntdata1$zerostart < 24] <- 1 #select all the locations
#where time after trip start is smaller than 273
huntdata2<-na.omit(huntdata1)                     #remove NA's generated
#nothunt3 and huntdata2 have an unequal number of column
#(nothunt3 has an additional column named "include1")
#in order to join the two, we need to add an additional column to huntdata2
huntdata2$include1<-huntdata2$include
#join the two data sets together
d5<-rbind(huntdata2,nothunt3)


#BRB
#useful reading includes:
#Benhamou (2011) Dynamic Approach to Space and Habitat Use Based on Biased
#Random Bridges. PLoS ONE 6: e14592
#Benhamou and Cornelis (2010) Incorporating movement behaviour and barriers to
#improve kernel home range space use estimates. Journal of Wildlife Management
#74: 1353 - 1360
#Calenge (2011) Home range estimation in R: the adehabitatHR package.
#from: cran.r-project.org/web/packages/adehabitatHR/vignettes/adehabitatHR.pdf
#now check to see how many locations you have for each individual
summary(d5$nsdid)
#remove individuals from the data set which have too few locations to estimate
#UD using BRB
notenough<-c("B")   #the names
#of individuals with insufficient data
d5$insufficient<-match(d5$nsdid,notenough,nomatch=0)      #label all
#the individuals with sufficient data with a 0
d5$insufficient[d5$insufficient > 1] <- 1   #label all the individuals
#with insufficient data with a "1"
d5$INSUF<-as.factor(d5$insufficient)      #make insufficient a factor
d6<-d5[d5$INSUF=="0",] #select individuals with sufficient data
d7<-na.omit(d6)            #remove the NA's
d7<-d5
#to show a summary of the points for each trip, for each individual, use
table(d7$nsdtrip,d7$nsdid)
#create a new trajectory with refined data set, divided by individual and trip
d8<-as.ltraj(d7$coordinates,d7$nsdtimeall,d7$nsdid,
             burst=d7$nsdid,typeII=TRUE)
summary(d8)
#make a 10m x 10m grid square of study area
xpoints<-c(346540:346670)#specific the extent of the study area on a global grid
xpoints1<-xpoints*10    #in metres, removing the last 2 digits
ypoints<-c(926080:926190)   #do the same for the y axis
ypoints1<-ypoints*10
pts = expand.grid(x = xpoints1, y = ypoints1)   #make the grid
grd.pts = SpatialPixels(SpatialPoints(pts)) #it has to be SpatialPixels to use
#in BRB
#calculate the diffusion parameter D for the BRB. BRB.lik is also available to
#estimate D
diffusion<-BRB.D(d8, #the new trajectory
                 Tmax = 5*60,              #the maximum time between relocations where
                 #smoothing should occur. Measured in seconds, so 120*60 for 120 mins, or 2 hours
                 Lmin = 0)                   #The smallest distance at which an animal should be
#considered moving, and therefore modelled in the UD. 0 if all data is included.
#make a UD using BRB
#first get a good value for hmin - use:
summary(d8[[1]])
#to find the mean distance travelled between locations
d9<-BRB(d8, #the trajectory
        diffusion,            #the diffusion parameter
        Tmax=5*60,  #maximum time between relocations: should be the same as smootherD
        Lmin=0,       #the same as for "diffusion"
        hmin=35,      #minimum smoothing parameter in units of locations
        #should be > mean interlocation distance/2
        grid = grd.pts,  #the gird in which to estimate UD
        b = FALSE,       # If TRUE, the relocation variance progressively merges with
        #the movement component; if FALSE, the relocation variance has a constant weight
        same4all = FALSE, #has to be FALSE if a grid is specified. If no grid is
        #specified, can be TRUE so the UD is estimated in the same area for each
        #individual
        extent=0.1,  #extent of the grid used for estimation
        tau=20)      #frequency of modelled relocations between known points.
#Measured in seconds
kerneloverlaphr(d9,    #to calculate overlap between trips
                method = c("HR"),      #type of overlap. HR is the proportion of the home range
                #of one individual / trip used by another
                percent = 95)          #Use percentage of home range for calculating overlap
#Extract the UD for each individual
d9a<-getvolumeUD(d9[[1]])   #select the UD for the first individual
#calculating the  are of use
IND1Area<-kernel.area(d9a,percent=seq(50,95,by=5))   #to get the area
#(in hectares) inside each % use between 50 and 95, at 5% intervals
IND1Area         #show the areas calculated
d9a1<-as.data.frame(d9a)   #change into a data frame that can be used by ruf
summary(d9a1)                    #check it looks ok
d9a1x<-data.frame(d9a1$Var2,d9a1$Var1,d9a1$n) #change dataframe variable
#order, as X and Y need to be the first two columns to use as.ppp
names(d9a1x)<-c("X","Y","UD")  #change the names to more sensible ones
d9a1x$include1[d9a1x$UD < 99] <- 1
#d9a1x$include1[d9a1x$UD < 99] <- 1  #assign 1 to any grid square where UD < 99
IND1<-na.omit(d9a1x)           #remove the grid squares where UD > 99
#open table and move X and Y to be the first columns
coord<-data.frame(IND1$X,IND1$Y)  #create new data frame with all grid
#coordinates
xysp<-SpatialPoints(coord)      #make the data frame into class "SpatialPoints",
#so it can be used to make a Minimum Convex Polygon
cp<-mcp(xysp,percent=100)      #create the minimum convex polygon
MCP<-as(cp, "owin")     #turn it into class "owin", so it can be used with the
#function "nncross"
Resource<-as.ppp(IND1,MCP)  #create an object of class ppp which specifies
#research area to use with nncross
Community<-c(3466080,9261222)    #location of the community
C1<-as.ppp(Community,MCP)    #make this class ppp so it can be used with nncross
Community1<-nncross(Resource,C1)#calculate distance between locations and the
#community
IND1$community<-Community1[,1]   #add the measurements to the datafile
#alternately you can load shapefiles
R<-readShapeSpatial ("h:\\River1.shp")   #1. load your shapefile
R1<-as.psp(R)     #2. make an object of class psp with your shapefile
River<-nncross(Resource,R1)  #3. for each point in the ppp "Resource",
#calculate the distance to the nearest point in "River"
IND1$river<-River[,1]       #4. add these distances to your datafile with UD
#repeat 1-4 for each landscape feature
#drawing the heatmap of UD shown in Figure 3
image(d9a,           #specify the data to use
      col=heat.colors         #specify color scheme
      (50))                  #specify how many different colors to use

xyzv<-as.image.SpatialGridDataFrame(d9a)     #create an object with the
#information required to add contours of use
contour(xyzv,          #contour information
        levels=c(50),          #which % use contour to add
        drawlabels=FALSE,      #can be TRUE or FALSE. If TRUE, adds a label of the % use
        # of the contour. If FALSE, no label added
        lwd=2,                 #width of the added contour line
        add=TRUE)              #can be TRUE or FALSE. If TRUE, contour line will be
#added to the existing image. If FALSE, a new image will be drawn
contour(xyzv,levels=c(95),              #to add 95% contour
        drawlabels=FALSE,lwd=2,lty=2,add=TRUE)  #with a dashed line (lty=2)
plot(R,col="black",pch=19,cex=40,add=TRUE) #if you wish to add landscape
#features

#RUF

#Read: Marzluff et al (2004) Relating resources to a probabilistic measure of
#space use: Forest fragments and Streller's Jays. Ecology 85: 1441 - 1427
#before starting


#check if your variables need to be transformed
#if UD distribution is heavily biased to higher percentages, consider 100-UD
hist(IND1$UD)                #to view a histogram of UD
IND1$UD2<-100-IND1$UD    #create a new variable where distribution will be
#biased to lower numbers, and therefore can be normalised using the natural log
#if explanatory variables cannot be normalised, or have another unusual
#distributions (i.e. strongly binomial), change them to categorial variables
#now fit each possible model and calculate the AIC
model1 <- ruf.fit(log(UD2)     #you can log, sqrt, asin your response variable
                  ~sqrt(community)+sqrt(river),  #put your explanatory variables here
                  space= ~ X + Y,               #specify which dataframe variables represent
                  #latitude and longitude
                  data=IND1,                 #which dataset to use
                  theta=c(0.2,2),               #which values to use for the Malvern correlation
                  #function. The first number is the range, which is the starting point from which
                  #ruf.fit will choose the best value for the range. It's a good idea to start
                  #with a low number. It is measured in metres. The second value is a smoothness
                  #parameter. It can be 0+, up to 10. It will not be estimated by ruf.fit - you
                  #need to vary it and choose the value which lowers to Malvern logLikelihood
                  standardized=FALSE)           #can be TRUE or FALSE. If FALSE, estimates for
#different indivuals can be compared to calculate a population estimate.
#If TRUE, all the estimates for all variables are shown on the same scale
#(within a single model), and the relative importance of each variable
#can be estimated.
summary(model1)             #to show results

#create a list of needed libraries
libraries <- c('leaps', 'tidyverse', 'grid')

#install libraries that aren't present
install.packages(setdiff(libraries, rownames(installed.packages())))

#load needed libraries
lapply(libraries, library, character.only = TRUE)

#initialize needed variables
bic_coef <- list()
adjr2_coef <- list()
cp_coef <- list()
label = ""

#set seed to ensure reproducibility
set.seed(50)

#set graphical parameters to desired values
par(cex = 0.0001)
par(cex.axis = 8000)

#make a small change to the regsubsets plot function to modify y axis labels as well as change rounding (disclaimer: this function was not written by me, a small edit was made to the existing plot.regsubsets function to allow for more the editing of the y axis labels to show the data more clearly)
plot.regsubsets<-function(x,labels=obj$xnames,main=NULL,
                          scale=c("bic","Cp","adjr2","r2"),
                          col=gray(seq(0,0.9,length=10)),digits=3,...){
  obj<-x
  lsum<-summary(obj)
  par(mar=c(7,5,6,3)+0.1)
  nmodels<-length(lsum$rsq)
  np<-obj$np
  propscale<-FALSE
  sscale<-pmatch(scale[1],c("bic","Cp","adjr2","r2"),nomatch=0)
  if (sscale==0)
    stop(paste("Unrecognised scale=",scale))
  if (propscale)
    stop(paste("Proportional scaling only for probabilities"))
  
  yscale<-switch(sscale,lsum$bic,lsum$cp,lsum$adjr2,lsum$rsq)
  up<-switch(sscale,-1,-1,1,1)
  
  index<-order(yscale*up)
  
  colorscale<- switch(sscale,
                      yscale,yscale,
                      -log(pmax(yscale,0.0001)),-log(pmax(yscale,0.0001)))
  
  image(z=t(ifelse(lsum$which[index,],
                   colorscale[index],NA+max(colorscale)*1.5)),
        xaxt="n",yaxt="n",x=(1:np),y=1:nmodels,xlab="", ylab = scale[1], col=col)
  
  laspar<-par("las")
  on.exit(par(las=laspar))
  par(las=2)
  axis(1,at=1:np,labels=labels)
  #axis(2,at=1:nmodels,labels= FALSE, tick = TRUE)
  axis(2,at=1:nmodels,labels=signif(yscale[index],digits))
  
  if (!is.null(main))
    title(main=main)
  box()
  invisible(NULL)
}

#generate the predictor and noise vectors
pred <- rnorm(n = 100)
noise <- rnorm(n = 100)
#choose values for the beta coefficients
beta <- c(20,30,40,50)

#create the response vector according to the formula provided in the question statement
resp <- beta[1] + beta[2]*pred + beta[3]*pred^2+beta[4]*pred^3+noise

#create predictors from X^1 all the way to X^10 and save them to a dataframe alongside the response
pred_resp_df <- data.frame(pred_one = pred, pred_two = pred^2, pred_three = pred^3, pred_four = pred^4, pred_five = pred^5, pred_six = pred^6, pred_sev = pred^7, pred_eight = pred^8, pred_nine = pred^9, pred_ten = pred^10, response = resp)

#remove unneeded variables
rm(pred, noise, beta, resp)

#use a for loop to perform the subset in default, forward stepwise selection and backward stepwise selection
for (i in 1:3) {
  
  #for the first instance of the loop, subset the dataframe in default conditions
  if (i == 1) {
    subset <- regsubsets(response ~ ., data = pred_resp_df, nvmax = 10)
  } 
  
  #for the second instance of the loop, subset the dataframe using the forward stepwise selection method
  else if (i == 2) {
    subset <- regsubsets(response ~ ., data = pred_resp_df, nvmax = 10, method = "forward")
    label = "forward"
  } 
  
  #for the third instance of the loop, subset the dataframe using the backward stepwise selection method
  else {
    subset <- regsubsets(response ~ ., data = pred_resp_df, nvmax = 10, method = "backward")
    label = "backward"
  }
  
  #save the summary of the subset
  sub_summary <- summary(subset)
  
  #create a dataframe which will be used to produce the line plot, save the number of predictors in one column and the cp value in the other
  plot <- data.frame(predictor_quantity = 1:10, cp = sub_summary$cp)
  #mutate the minimum cp value such that it will be red in the plot 
  plot <- plot %>%
    mutate(color = (min(cp) == cp))
  
  #set the position of the predictors plot
  par(fig = c(0.55, 0.98, 0.73, 0.95))
  
  #plot the predictors that are present in each cp value obtained
  plot(subset,scale = "Cp")
  
  #plot a line graph showing how many predictors yields the lowest cp value
  plot_line <- ggplot(plot, mapping = aes(x = predictor_quantity, y = cp)) + 
          geom_line() + 
          geom_point(aes(color = color, size = 2.5 )) + 
          ggtitle(paste0(label, " Cp vs number of predictors")) + xlab("number of predictors") + 
          ylab("Cp value (lower is better)") +
          scale_color_manual(values = c(NA, "red")) + 
          theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank()) + 
          scale_x_continuous(breaks = seq(0,10,1))

  #set the position of the line plot to be adjacent to the predictors
  vp <- viewport(height=unit(0.3, "npc"), width=unit(0.5, "npc"), 
                             just=c("right","bottom"), 
                             y=0.7, x=0.5)
  #print the line plot
  print(plot_line, vp = vp)
  
  #save the cp suggested coefficients in a list. This list will be saved in another list for each instance of the loop to yield a list of 3 lists (one for default, one for forward stepwise and one for backward stepwise)
  cp_coef[[i]] <- as.list(coef(subset, which.min(sub_summary$cp)))
  
  #create a dataframe which will be used to produce the line plot, save the number of predictors in one column and the bic value in the other
  plot <- data.frame(predictor_quantity = 1:10, bic = sub_summary$bic)
  #mutate the minimum bic value such that it will be red in the plot
  plot <- plot %>%
    mutate(color = (min(bic) == bic))

  #set the position of the predictors plot
  par(fig = c(0.55, 0.98, 0.405, 0.64), new = TRUE)
 
   #plot the predictors that are present in each bic value obtained
  plot(subset,scale = "bic")
  
  #plot a line graph showing how many predictors yields the lowest bic value
  plot_line <- ggplot(plot, mapping = aes(x = predictor_quantity, y = bic)) + 
          geom_line() + 
          geom_point(aes(color = color, size = 2.5 )) + 
          ggtitle(paste0(label, " bic vs number of predictors")) + xlab("number of predictors") + 
          ylab("bic value (lower is better)") +
          scale_color_manual(values = c(NA, "red")) + 
          theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank()) + 
          scale_x_continuous(breaks = seq(0,10,1))
  
  #set the position of the line plot to be adjacent to the predictors
  vp <- viewport(height=unit(0.3, "npc"), width=unit(0.5, "npc"), 
                             just=c("right","bottom"), 
                             y= 0.39, x=0.5)
  #print the line plot
  print(plot_line, vp = vp)

  #save the bic suggested coefficients in a list. This list will be saved in another list for each instance of the loop to yield a list of 3 lists (one for default, one for forward stepwise and one for backward stepwise)
  bic_coef[[i]] <- as.list(coef(subset, which.min(sub_summary$bic)))
  
  #create a dataframe which will be used to produce the line plot, save the number of predictors in one column and the adjusted R squared value in the other
  plot <- data.frame(predictor_quantity = 1:10, adjr2 = sub_summary$adjr2)
  #mutate the maximum adjusted R squared value such that it will be red in the plot 
   plot <- plot %>%
    mutate(color = (max(adjr2) == adjr2))
  
   #set the size of the x axis labels
   par(cex.axis = 8000)
   #set the position of the predictors plot
   par(fig = c(0.55, 0.98, 0.1, 0.34), new = TRUE)
   
  #plot the predictors that are present in each adjusted R squared value obtained
  plot(subset, scale = "adjr2", main = paste0(label, " predictors in adjusted R squared subset"))
  
  #plot a line graph showing how many predictors yields the highest adjusted R squared value
  plot_line <- ggplot(plot, mapping = aes(x = predictor_quantity, y = adjr2)) + 
          geom_line() + 
          geom_point(aes(color = color, size = 2.5 )) + 
          ggtitle(paste0(label, " adjusted R squared vs number of predictors")) + 
          xlab("number of predictors") + ylab("R squared value (higher is better)") +
          scale_color_manual(values = c(NA, "red")) + theme(legend.position = "none") + 
          scale_x_continuous(breaks = seq(0,10,1))  

  #set the position of the line plot to be adjacent to the predictors 
  vp <- viewport(height=unit(0.35, "npc"), width=unit(0.5, "npc"), 
                             just=c("right","bottom"), 
                             y= 0.03, x=0.5)
  #print the line plot 
  print(plot_line, vp = vp)
  
  #save the adjusted R squared suggested coefficients in a list. This list will be saved in another list for each instance of the loop to yield a list of 3 lists (one for default, one for forward stepwise and one for backward stepwise)
  adjr2_coef[[i]] <- as.list(coef(subset, which.max(sub_summary$adjr2)))
  
  #remove unneeded variables
  rm(plot, sub_summary, subset, i, label, plot_line, vp)
}

# packages
library(pROC)

# import data
fdata = read.csv("malaria.csv")
head(fdata)

# remove Nid column
data = fdata[,-1]
head(data)

# check for null values
sum(is.na(data))

# change categorical variables using factor
data$source = factor(data$source)
data$behavior = factor(data$behavior)
data$nettype = factor(data$nettype)
data$district = factor(data$district)
data$work = factor(data$work)

# characteristics of subjects
summary(data$insecticide)
summary(data$stress)
summary(data$health)


# remove outlier
data = data[(data$insecticide < 424.01),]
data = data[(data$health != 37),]
data = data[(data$district != "9Moon"),]
data$district = droplevels(data$district, exclude = "9Moon")
data = data[(data$stress>0),]

# check for null values after removing outliers
sum(is.na(data))

# histograms
par(mfrow = c(1, 3))
hist(data$stress)
hist(data$insecticide, breaks = seq(0,350, by=50))
hist(data$health)
par(mfrow = c(1, 1))

# Plot of Contingency Table ---------------------------------
par(mfrow = c(2, 3))
mosaicplot(table(data$source, data$malaria),
           color = c("Red", "Blue"),
           xlab="Source",
           ylab="Malaria",
           main="Malaria vs Source",
           cex.axis = 0.8)

mosaicplot(table(data$behavior, data$malaria),
           color = c("Red", "Blue"),
           xlab="Behavior",
           ylab="Malaria",
           main="Malaria vs Behavior",
           cex.axis = 0.8)

mosaicplot(table(data$nettype, data$malaria),
           color = c("Red", "Blue"),
           xlab="Net Type",
           ylab="Malaria",
           main = "Malaria vs Net Type",
           cex.axis = 0.8)

mosaicplot(table(data$district, data$malaria),
           color = c("Red", "Blue"),
           xlab="District",
           ylab="Malaria",
           main="Malaria vs District",
           cex.axis = 0.8)

mosaicplot(table(data$work, data$malaria),
           color = c("Red", "Blue"),
           xlab="Work",
           ylab="Malaria",
           main="Malaria vs Work",
           cex.axis = )
par(mfrow = c(1, 1))

# Chi-sq test for categorical variable ---------------------------------
chisq.test(table(data$malaria, data$source))
chisq.test(table(data$malaria, data$behavior))
chisq.test(table(data$malaria, data$nettype))
chisq.test(table(data$malaria, data$district))
chisq.test(table(data$malaria, data$work))

# LR test for continuous variable ---------------------------------
glm.stress = glm(malaria ~ stress, data=data, family="binomial")
glm.insecticide = glm(malaria ~ insecticide, data=data, family="binomial")
glm.health = glm(malaria ~ health, data=data, family="binomial")

pchisq(anova(glm.stress)[2,2], anova(glm.stress)[2,1], lower.tail = FALSE)
pchisq(anova(glm.insecticide)[2,2], anova(glm.insecticide)[2,1], lower.tail = FALSE)
pchisq(anova(glm.health)[2,2], anova(glm.health)[2,1], lower.tail = FALSE)

# slicing-dicing plot of empirical log-odds ---------------------------------

# Stress
stress.frac = factor(cut(data$stress,breaks=seq(0,20, by=2)))
stressdata <- data.frame(data, stress.frac)
e.probs= tapply(data$malaria, stress.frac, mean)
e.logits.stress <- log(e.probs/(1-e.probs))

# insecticide
insecticide.frac = factor(cut(data$insecticide,breaks=seq(0,350, by=50)))
insecticidedata <- data.frame(data, insecticide.frac)
e.probsins= tapply(data$malaria, insecticide.frac, mean)
e.logit.insecticide <- log(e.probsins/(1-e.probsins))

# health
health.frac = factor(cut(data$health, breaks=seq(5,34, by=3)))
health_data <- data.frame(data, health.frac)
e.probsins = tapply(data$malaria, health.frac, mean)
e.logitsins.health <- log(e.probsins/(1-e.probsins))

par(mfrow = c(1, 3))
plot(seq(1,19, by=2), e.logits.stress, xlab="Stress", col = "red",
     pch = 16, main= "Logit for stress", ylab = "Empirical Log Odds")
plot(seq(25, 325, by=50), e.logit.insecticide, xlab="Insecticide", col = "red",
     pch = 16,main= "Logit for insecticide", ylab = "Empirical Log Odds")
plot(seq(3, 34, length.out=length(e.logitsins.health)),e.logitsins.health,
     xlab="Health", col = "red", pch = 16,
     main= "Logit for Health", ylab = "Empirical Log Odds")
par(mfrow = c(1, 1))

# Model Fitting ---------------------------------

# test transformations of insecticide
# insecticide
summary(glm(malaria ~ insecticide, data = data, family="binomial"))
BIC(glm(malaria ~ insecticide, data = data, family="binomial"))

# sqrt insecticide
summary(glm(malaria ~ sqrt(insecticide), data = data, family="binomial"))
BIC(glm(malaria ~ sqrt(insecticide), data = data, family="binomial"))

# insecticide^2
summary(glm(malaria ~ I(insecticide^2), data = data, family="binomial"))
BIC(glm(malaria ~ I(insecticide^2), data = data, family="binomial"))


# Anova tests for Model Variations ---------------------------------

# all the associated variables
model.associated = glm(malaria ~ nettype + district + work + stress + insecticide,
                 family= "binomial",
                 data = data)
summary(model.associated)
BIC(model.associated) # 868.0192
AIC(model.associated) # 831.166

#all associated variables but with sqrt(insecticide)
model.transform<-glm(malaria ~ nettype + district + work + stress + sqrt(insecticide),
    family= "binomial",
    data = data)
summary(model.transform) #AIC = 831.05
BIC(model.transform) #BIC = 867.9061

# remove net type
model.rem.nettype = glm(malaria ~ district + work + stress + sqrt(insecticide),
                        family= "binomial",
                        data = data)
summary(model.rem.nettype)
AIC(model.rem.nettype) # 840.3344
BIC(model.rem.nettype) # 872.5809

#remove district
model.rem.district = glm(malaria ~ nettype + work + stress + sqrt(insecticide),
                         family= "binomial",
                         data = data)
summary(model.rem.district)
AIC(model.rem.district) # 859.722
BIC(model.rem.district) # 887.3624

# remove work

model.rem.work = glm(malaria ~ nettype + district + stress + sqrt(insecticide),
                    family= "binomial",
                    data = data)
summary(model.rem.work)
AIC(model.rem.work) # 829.6139
BIC(model.rem.work) # 857.2538

#remove stress
model.rem.stress = glm(malaria ~ nettype + district + work + sqrt(insecticide),
                     family= "binomial",
                     data = data)
summary(model.rem.stress)
AIC(model.rem.stress) # 900.2832
BIC(model.rem.stress) # 932.5298

#remove insecticide
model.rem.stress = glm(malaria ~ nettype + district + work + stress,
                       family= "binomial",
                       data = data)
summary(model.rem.stress)
AIC(model.rem.stress) # 834.129
BIC(model.rem.stress) # 866.3756

# From the outs above, the model with the best AIC is the
# model with netype, district, stress, and sqrt(insecticide) as variables.

# Adding interaction to best fit model ---------------------------------

# nettype*sqrt(insecticide)
summary(glm(data$malaria~data$nettype+ data$district + data$stress + sqrt(data$insecticide)+ data$nettype:sqrt(data$insecticide), family= "binomial")) #829.56
BIC(glm(data$malaria~data$nettype+ data$district + data$stress + sqrt(data$insecticide)+ data$nettype:sqrt(data$insecticide), family= "binomial")) #861.8088

# nettype*insecticide
summary(glm(data$malaria~data$nettype+ data$district + data$stress + data$insecticide+ data$nettype:data$insecticide, family= "binomial")) #828.65
BIC(glm(data$malaria~data$nettype+ data$district + data$stress + data$insecticide+ data$nettype:data$insecticide, family= "binomial")) # 860.8999

# stress*work
summary(glm(data$malaria~data$nettype+ data$district + data$stress+ data$work + data$insecticide + data$stress:data$work, family= "binomial")) #834.32
BIC(glm(data$malaria~data$nettype+ data$district + data$stress+ data$insecticide + + data$work + data$stress:data$work, family= "binomial")) # 880.3883

#stress*insecticide
summary(glm(data$malaria~data$nettype+ data$district + data$stress + data$insecticide + data$stress:data$insecticide, family= "binomial")) #831.66
BIC(glm(data$malaria~data$nettype+ data$district + data$stress + data$insecticide + data$stress:data$insecticide, family= "binomial")) #863.9018

# take away district from best AIC model
summary(glm(data$malaria~data$nettype + data$stress + data$insecticide+ data$nettype:data$insecticide, family= "binomial")) #858.3
BIC(glm(data$malaria~data$nettype + data$stress + data$insecticide+ data$nettype:data$insecticide, family= "binomial")) # 881.3299

# Parameter estimates for best AIC model
best.fit<-glm(malaria ~ nettype + district + stress + insecticide + nettype*insecticide,
              family= "binomial",
              data = data)
summary(best.fit)

# Selection Algorithm Check ---------------------------------

# Use stepwise regression algorithm as a sanity check and test whether the
# method used above agrees with the algorithm

# fit model with all parameters
glm.all = glm(malaria ~., data = data, family="binomial")
summary(glm.all)
pchisq(956.98-807.19, 13, lower.tail=FALSE)

# remove variables that are not statistically significant in output
# summary above
glm.reduce = glm(malaria ~. - source - health - work, data=data, family="binomial")
summary(glm.reduce)

# LR test for variables that were not statistically significant. p-value shows
# that source, health, work is not statistically significant
anova(glm.reduce, glm.all)
pchisq(3.754, 4, lower.tail=FALSE) # 0.4403207

# create intercept model for selection algorithm
glm.intercept = glm(malaria ~ 1, data=data, family="binomial")

# forward selection algorithm
glm.forward = step(glm.intercept, direction='forward', scope=formula(glm.all))
formula(glm.forward) # malaria ~ stress + district + nettype + insecticide

# backward selection algorithm
glm.backward = step(glm.all, direction='backward')
formula(glm.backward) # malaria ~ stress + insecticide + nettype + district

# both direction
glm.both = step(glm.intercept, direction='both', scope=formula(glm.all))
formula(glm.both) # malaria ~ stress + district + nettype + insecticide

# check to see if variables omitted from glm.both are not statistically
# significant. p-value shows that variables can be omitted
anova(glm.both, glm.all)
pchisq(1.8935, 4, lower.tail = FALSE) # 0.755339

# check for interaction
glm.interaction = glm(malaria ~ stress + district + nettype + insecticide + stress*district + stress*nettype + stress*insecticide + district*nettype + district*insecticide + nettype*insecticide, data = data, family = "binomial")
summary(glm.interaction)

# forward
glm.interact.forward = step(glm.both, direction='forward', scope=formula(glm.interaction))
formula(glm.interact.forward) # malaria ~ stress + district + nettype + insecticide + nettype:insecticide

# backward
glm.interact.backward = step(glm.interaction, direction='backward')
formula(glm.interact.backward) # malaria ~ stress + district + nettype + insecticide + nettype*insecticide

# both direction
glm.interact.both = step(glm.both, direction='both', scope=formula(glm.interaction))
formula(glm.interact.both) # malaria ~ stress + district + nettype + insecticide + nettype*insecticide

# LR to see if removed interaction terms are statistically significant
anova(glm.interact.both, glm.interaction)
pchisq(6.9708, 8, lower.tail = FALSE) # 0.5397864

# LR to see if interaction term nettype:insecticide cab be omitted
# seems that interaction term is not relevant
anova(glm.both, glm.interact.both)
pchisq(3.057, 1, lower.tail = FALSE) # 0.08038997

# Odds Ratio, CI, p-value ---------------------------------

# Odds Ratio
oddsratios<-exp(best.fit$coefficients); oddsratios

# Confidence Interval for Odds Ratio
exp(confint(best.fit))

# p-value
summary(best.fit)$coefficients[,4]

# Classification Table, Goodness of Fit, ROC Curve -----------------------------

# classification table / confusion matrix
best.fitted = best.fit$fitted.values
best.predict = ifelse(best.fit$fitted.values >= 0.5, 1, 0)
table(data$malaria, best.predict)

# goodness of fit
best.fit$deviance
best.fit$df.residual
1-pchisq(814.6534, 733) #0.0189730

# roc curve
plot.roc(data$malaria, best.fit$fitted.values, print.auc=TRUE, quiet=TRUE,
         main = "ROC Curve for Fitted Model")
auc(data$malaria, best.fit$fitted.values) # 0.756


# Plot of Success Probabilities ---------------------------------

#stress at each district with nettypeA and insecticide=140
p.success = tapply(data$malaria, factor(cut(data$stress, seq(0,20, by=2))), mean)

new1 = data.frame(stress = 0:20, district="1North", nettype = "TypeA", insecticide = 140)
pred.dist1 = predict(best.fit, newdata = new1, type="response")
prob.dist1 = exp(pred.dist1)/(1+exp(pred.dist1))

new2 = data.frame(stress = 0:20, district="2East", nettype = "TypeA", insecticide = 140)
pred.dist2 = predict(best.fit, newdata = new2, type="response")
prob.dist2 = exp(pred.dist2)/(1+exp(pred.dist2))

new3 = data.frame(stress = 0:20, district="3South", nettype = "TypeA", insecticide = 140)
pred.dist3 = predict(best.fit, newdata = new3, type="response")
prob.dist3 = exp(pred.dist3)/(1+exp(pred.dist3))

#insecticide for each net type in 2East at stress=10

p.success.insecticide = tapply(data$malaria, factor(cut(data$insecticide,breaks=seq(0,350, by=50))), mean)
new4 = data.frame(insecticide = 0:350, district="2East", nettype = "TypeA", stress = 10)
pred.netA = predict(best.fit, newdata = new4, type="response")
prob.netA = exp(pred.netA)/(1+exp(pred.netA))

new5 = data.frame(insecticide = 0:350, district="2East", nettype = "TypeB", stress = 10)
pred.netB = predict(best.fit, newdata = new5, type="response")
prob.netB = exp(pred.netB)/(1+exp(pred.netB))

par(mfrow=c(1, 2))
plot(x=seq(1,19, by=2), y=p.success, type="n", xlim=c(0, 20),ylim=c(0.4,0.8),xlab="Stress", ylab="Probability of Malaria", main="Probability of Malaria - s and d")
lines(0:20, prob.dist1, col="red", lwd=1.5)
lines(0:20, prob.dist2, col="blue", lwd=1.5)
lines(0:20, prob.dist3, col="green", lwd=1.5)
legend("topleft", c("1North","2East","3South"), col=c("red", "blue", "green"), lty=1)

plot(x=seq(25,325, by=50), y=p.success.insecticide, type="n", xlim=c(0, 350),ylim=c(0.5,0.65),xlab="Insecticide", ylab="Probability of Malaria", main="Probability of Malaria - i and nt")
lines(0:350, prob.netA, col="red", lwd=1.5)
lines(0:350, prob.netB, col="blue", lwd=1.5)
legend("topright", c("NetA","NetB"), col=c("red", "blue"), lty=1)
par(mfrow=c(1, 1))

##3 way plots
# plot for stress and district
new1 = data.frame(stress = 0:20, district="1North", nettype = "TypeA", insecticide = 141.0)
pred.dist1 = predict(best.fit, newdata = new1, type="response")

new2 = data.frame(stress = 0:20, district="2East", nettype = "TypeA", insecticide = 141.0)
pred.dist2 = predict(best.fit, newdata = new2, type="response")

new3 = data.frame(stress = 0:20, district="3South", nettype = "TypeA", insecticide = 141.0)
pred.dist3 = predict(best.fit, newdata = new3, type="response")

new4 = data.frame(stress = 0:20, district="1North", nettype = "TypeB", insecticide = 141.0)
pred.dist4 = predict(best.fit, newdata = new4, type="response")

new5 = data.frame(stress = 0:20, district="2East", nettype = "TypeB", insecticide = 141.0)
pred.dist5 = predict(best.fit, newdata = new5, type="response")

new6 = data.frame(stress = 0:20, district="3South", nettype = "TypeB", insecticide = 141.0)
pred.dist6 = predict(best.fit, newdata = new6, type="response")

par(mfrow=c(1, 2))
plot(0:20, pred.dist1, type="l", ylab="Probability of Malaria", xlab="Stress", col="red", ylim = c(0, 1), main = "Probability of Malaria with TypeA Net")
lines(0:20, pred.dist2, type="l", col="blue")
lines(0:20, pred.dist3, type="l", col="green")
legend("topleft", c("1North","2East","3South"), col=c("red", "blue", "green"), lty=1)

plot(0:20, pred.dist4, type="l", ylab="Probability of Malaria", xlab="Stress", col="red", ylim = c(0, 1), main = "Probability of Malaria with TypeB Net")
lines(0:20, pred.dist5, type="l", col="blue")
lines(0:20, pred.dist6, type="l", col="green")
legend("topleft", c("1North","2East","3South"), col=c("red", "blue", "green"), lty=1)
par(mfrow=c(1, 2))

# plot for insecticide and nettype
new1 = data.frame(stress = 10.40, district="1North", nettype = "TypeA", insecticide = 0:350)
pred.insecticide.1 = predict(best.fit, newdata = new1, type="response")

new2 = data.frame(stress = 10.40, district="2East", nettype = "TypeA", insecticide = 0:350)
pred.insecticide.2 = predict(best.fit, newdata = new2, type="response")

new3 = data.frame(stress = 10.40, district="3South", nettype = "TypeA", insecticide = 0:350)
pred.insecticide.3 = predict(best.fit, newdata = new3, type="response")

new4 = data.frame(stress = 10.40, district="1North", nettype = "TypeB", insecticide = 0:350)
pred.insecticide.4 = predict(best.fit, newdata = new4, type="response")

new5 = data.frame(stress = 10.40, district="2East", nettype = "TypeB", insecticide = 0:350)
pred.insecticide.5 = predict(best.fit, newdata = new5, type="response")

new6 = data.frame(stress = 10.40, district="3South", nettype = "TypeB", insecticide = 0:350)
pred.insecticide.6 = predict(best.fit, newdata = new6, type="response")

par(mfrow=c(1, 3))
plot(0:350, pred.insecticide.1, type="l", ylab="Probability of malaria", col="red", ylim = c(0.15, 0.45), xlab="Insecticide", main="Probability of Malaria at 1North")
lines(0:350, pred.insecticide.4, type="l", col="blue")
legend("topright", c("TypeA", "TypeB"), col=c("red", "blue"), lty=1)

plot(0:350, pred.insecticide.2, type="l", ylab="Probability of malaria", col="red", ylim = c(0.15, 0.6), xlab="Insecticide", main="Probability of Malaria at 2East")
lines(0:350, pred.insecticide.5, type="l", col="blue")
legend("topright", c("TypeA", "TypeB"), col=c("red", "blue"), lty=1)

plot(0:350, pred.insecticide.3, type="l", ylab="Probability of malaria", col="red", ylim = c(0.25, 0.8), xlab="Insecticide", main="Probability of Malaria at 3South")
lines(0:350, pred.insecticide.6, type="l", col="blue")
legend("topright", c("TypeA", "TypeB"), col=c("red", "blue"), lty=1)
par(mfrow=c(1, 1))


library(survival)
library(survminer)
data = read.csv("C:/Users/joash/OneDrive/Desktop/AssignMent/pancancerInfo.csv")
colnames(data)[10] = 'OS'

#Subset Acute Leukemia Cancer
data = data[data$type == 'LAML',]

#Remove Unknown columns
data = data[data$race != '[Not Available]',] 

#Gender
pdf('Gender.pdf',width = 20,height = 10)
fit=survfit(Surv(OS,Event)~gender,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)
dev.off()

#Age
colnames(data)[3] = 'Age'
data$Age[data$Age>50] = 'Above_50'
data$Age[data$Age<=50] = 'Below_50'
pdf('Age.pdf',width = 20,height = 10)
fit=survfit(Surv(OS,Event)~Age,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)
dev.off()

#Gender and Age
pdf('GenderAge.pdf',width = 20,height = 10)
fit=survfit(Surv(OS,Event)~Age + gender,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)
dev.off()


#Age and Race
pdf("DiffAge.pdf", width=20, height=10,onefile = TRUE)
for (i in unique(data$race))
{
  fit=survfit(Surv(OS,Event)~race + Age,data=data[data$race == i,])
  plt = ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)
  print(plt)
}
dev.off()

#Race
pdf('Race.pdf',width = 20,height = 10,onefile = FALSE)
fit=survfit(Surv(OS,Event)~race ,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)
dev.off()

#Race and Gender
pdf("DiffRace.pdf", width=20, height=10,onefile = TRUE)
for (i in unique(data$race))
{
fit=survfit(Surv(OS,Event)~race + gender,data=data[data$race == i,])
plt = ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)
print(plt)
}
dev.off()

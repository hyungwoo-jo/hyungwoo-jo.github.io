library(survival);library(dplyr);library(ggplot2);library(purrr);library(cmprsk)

# makedata --------------------

set.seed(181)
data <- mgus2
data$etime <- with(data, ifelse(pstat==0, futime, ptime))
data$event <- with(data, ifelse(pstat==0, 2*death, 1))
data$e<-with(data, ifelse(pstat==0, 2*death, 1))
data$event <- factor(data$event, 0:2, labels=c("censor", "event", "competingrisk"))
sample_sizes <- c(censor = 8, event = 4, competingrisk = 3)
data1 <- data %>%
  group_split(event) %>%
  map2(.x = ., .y = sample_sizes, ~ slice_sample(.x, n = .y, replace = TRUE)) %>%
  bind_rows()
data1 <- data1 %>%
  mutate(id = paste0(id, "(", event, ")"))
data1


# survfit -------------------

pp<-survfit(Surv(etime, event)~1, data = data1)
pp
summary(pp)

# finegray function in survival package ------------------
pdata <- finegray(Surv(etime, event) ~ ., data=data1)
pdata




# Build G(t) = the KM of the censoring distribution
# An event at time x is not "at risk" for censoring at time x (Geskus 2011)
tt <- sort(unique(data1$etime))  # all the times
ntime <- length(tt)
nrisk <- nevent <- double(ntime)
for (i in 1:ntime) {
  nrisk[i] <- sum((data1$etime > tt[i] & data1$e >0) | (data1$etime >= tt[i] & data1$e==0))
  nevent[i] <- sum(data1$etime == tt[i] & data1$e==0)
}
G <- cumprod(1- nevent/nrisk)
G

# with cmprsk package -----------------------
data_cmprsk <- data1 %>%
  mutate(
    ftime = etime,            
    fstatus = case_when(     
      event == "censor" ~ 0,  
      event == "event" ~ 1,    
      event == "competingrisk" ~ 2 
    )
  )
dat<-data_cmprsk %>% select(ftime, fstatus, sex)
cif_fit <- cuminc(
  ftime = data_cmprsk$ftime,
  fstatus = data_cmprsk$fstatus
)
cif_fit
specific_times <- c(39,49,77,80,90,94,127, 128,184, 222,226, 259, 280, 394)


cif_at_times <- timepoints(cif_fit, specific_times)

cif_at_times

fit<-survfit(Surv(etime, event)~1, data = data1)
fit<-summary(fit, times = specific_times)
fit
plot(cif_fit,
     main = "Cumulative Incidence Function (CIF)",
     xlab = "Time",
     ylab = "CIF",
     col = c("blue", "red", "green"), # 그룹별 색상 지정
     lty = 1:3)   
fgfit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1,
               weight=fgwt, data=pdata, model = T)



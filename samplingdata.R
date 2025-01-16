library(survival)
library(cmprsk)
library(dplyr)

set.seed(123)

# Define sample size
N <- 10000

# Assign Cancer Status: 10 Cancer (1) and 10 Non-Cancer (0)
Cancer_status <- c(rep(0, N/2), rep(1, N/2))

# Assign Sex randomly: 0 = Female, 1 = Male
Sex <- sample(c(0,1), size = N, replace = TRUE)

# Define event rates
lambda_cvd <- 0.2  # CVD death rate (per year)
lambda_non_cvd_no_cancer <- 0.1  # Non-CVD death rate for non-cancer patients
lambda_non_cvd_cancer <- 0.5      # Non-CVD death rate for cancer patients

# Initialize vectors to store times and events
CVD_time <- numeric(N)
non_CVD_time <- numeric(N)
Time <- numeric(N)
Event <- numeric(N)

# Simulate times
for (i in 1:N) {
  # Simulate time to CVD death (same for all)
  CVD_time[i] <- rexp(1, rate = lambda_cvd)
  
  # Simulate time to non-CVD death based on Cancer status
  if (Cancer_status[i] == 0) {
    non_CVD_time[i] <- rexp(1, rate = lambda_non_cvd_no_cancer)
  } else {
    non_CVD_time[i] <- rexp(1, rate = lambda_non_cvd_cancer)
  }
  
  # Determine which event occurs first
  if (CVD_time[i] < non_CVD_time[i]) {
    Time[i] <- CVD_time[i]
    Event[i] <- 1  # CVD Death
  } else {
    Time[i] <- non_CVD_time[i]
    Event[i] <- 2  # Non-CVD Death
  }
}

# Create the dataset
sample_data <- data.frame(
  ID = 1:N,
  Time = round(Time, 2),
  Event = Event,
  Sex = ifelse(Sex == 1, "Male", "Female"),
  Cancer = ifelse(Cancer_status == 1, "Yes", "No"),
  CVD_time = round(CVD_time, 2)  # Counterfactual CVD death time
)

# View the simulated dataset
print("Simulated Dataset:")
print(sample_data)

sample_data <- sample_data %>%
  mutate(
    Sex_num = ifelse(Sex == "Male", 1, 0),
    Cancer_num = ifelse(Cancer == "Yes", 1, 0)
  )


surv_object_csh <- Surv(time = sample_data$Time, event = sample_data$Event == 1)
cox_csh <- coxph(surv_object_csh ~ Cancer_num, data = sample_data)

# Summary of the Cause-Specific Hazard Model
print("Cause-Specific Cox Model Summary:")
print(summary(cox_csh))


ftime <- sample_data$Time
fstatus <- sample_data$Event
covariates <- sample_data %>% select(Cancer_num)
fg_model <- crr(ftime = ftime, fstatus = fstatus, cov1 = covariates)
print("Fine-Gray Model Summary:")
print(summary(fg_model))
group <- sample_data$Cancer_num
cif_cvd <- cuminc(ftime, fstatus, group = group, cencode = 0)
plot(cif_cvd, 
     xlab = "Time (years)", 
     ylab = "Cumulative Incidence",
     main = "Cumulative Incidence of CVD Death by Cancer Status",
     col = c("blue", "purple", 'red', 'orange'),
     lty = 1:2,
     lwd = 2)


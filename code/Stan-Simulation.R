libraries = c("dplyr","magrittr","tidyr","ggplot2","rstan","readxl","loo")
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE) }

require(zoo)
require(lubridate)

base_sz = 12 # base_size parameter
theme_set(theme_bw())

'%&%' = function(x,y) paste0(x,y)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

packageVersion("rstan")
packageVersion("StanHeaders")
rstan::stan_version()

## read files


datafilename = "data/Figure1.xlsx"
#datafilename
read_excel(datafilename, sheet="sample") -> df
names(df)

mean_SI=mean(df$dist)
sd_SI=sd(df$dist)



#stan simulation
stanmaindir = 'Figure1A'
unlink(stanmaindir, recursive=T)
dir.create(stanmaindir)

filenames=stanmaindir%&%'/WAIC.txt'

#No truncation
print("No truncation\n")

#Lognormal distribution

## main dir for Stan simulations
standirname = stanmaindir%&%"/lognormal-no_truncation"

unlink(standirname, recursive=T)
dir.create(standirname)
#standirname

# Dumping data

data_input<-list(
    N = nrow(df),
    E_L = df$E_L,
    E_R = df$E_R,
    S_L = df$S_L,
    S_R = df$S_R
)



# Dumping initial conditions
logmean_SI = log(mean(df$dist))
logsd_SI = log(sd(df$dist))
s_raw = rep(.5, nrow(df))
e_raw = rep(.5, nrow(df))


lognormal_n_model <- stan(warmup = 100, data = data_input,
                        iter = 200,
                        chains = 0,
                        model_code = "
data {
    int<lower = 0> N; // number of records
    vector<lower = 0>[N] E_L;
    vector<lower = 0>[N] E_R;
    vector<lower = 0>[N] S_L;
    vector<lower = 0>[N] S_R;
}

parameters {
    real logmean_SI;
    real logsd_SI;
    vector<lower = 0, upper = 1>[N] s_raw;
    vector<lower = 0, upper = 1>[N] e_raw;
}

transformed parameters {
    real<lower = 0> param2 = sqrt(log((exp(2*(logsd_SI-logmean_SI))+1.0)));
    real param1 = logmean_SI - param2^2/2.0;

    vector<lower = min(S_L), upper = max(S_R)>[N] s;
    vector<lower = min(E_L), upper = max(E_R)>[N] e;

    s = S_L + (S_R - S_L) .* s_raw;
    for (k in 1:N) 
        if (E_R[k] > s[k]) 
            e[k] = E_L[k] + (s[k] - E_L[k]) * e_raw[k];
        else
            e[k] = E_L[k] + (E_R[k] - E_L[k]) * e_raw[k];
}

model {
    logmean_SI ~ std_normal();
    logsd_SI ~ std_normal();

    e_raw ~ normal(0.5, 1.0);
    s_raw ~ normal(0.5, 1.0);

    for (k in 1:N)
        target += lognormal_lpdf(s[k] - e[k] | param1, param2);
}

generated quantities {
    real<lower = 0> mean_SI = exp(param1 + param2^2/2);
    real<lower = 0> sd_SI = sqrt((exp(param2^2)-1)*exp(2*param1+param2^2));

    vector[N] log_likelihood;
    for (k in 1:N) 
        log_likelihood[k] = lognormal_lpdf(s[k] - e[k] | param1, param2);
}
" 
)


lognormal_n_fit<-stan(fit=lognormal_n_model, data = data_input,
              warmup = 10000,
              iter = 20000,
              chains = 10,
              sample_file=standirname%&%'/sample',
              diagnostic_file=standirname%&%'/diagnostic',
              control=list(adapt_delta=0.98),
              cores=10
              )

print("lognormal_n_fit\n")

log_lik<-extract_log_lik(lognormal_n_fit, parameter_name = "log_likelihood", merge_chains = TRUE)
w<-waic(log_lik)
write.table("lognormal",file = filenames,append=T)
write.table(w$estimates,file = filenames,append=T)



#Gamma distribution

## main dir for Stan simulations
standirname = stanmaindir%&%"/gamma-no_truncation"
unlink(standirname, recursive=T)
dir.create(standirname)


# Dumping initial conditions
e_raw_g_n = rep(.2, nrow(df))
s_raw_g_n = rep(.8, nrow(df))
param1 = (mean(df$dist)/sd(df$dist))^2
param2 = mean(df$dist)/(sd(df$dist)^2)
 

# Stan program

gamma_n_model <- stan(warmup = 100, data = data_input,
                        iter = 200,
                        chains = 0,
                        model_code = "
data {
    int<lower = 0> N; // number of records
    vector<lower = 0>[N] E_L;
    vector<lower = 0>[N] E_R;
    vector<lower = 0>[N] S_L;
    vector<lower = 0>[N] S_R;
}

parameters {
    real<lower=0> mean_SI;
    real<lower=0> sd_SI;

    vector<lower = 0, upper = 1>[N] e_raw_g_n;
    vector<lower = 0, upper = 1>[N] s_raw_g_n;
}

transformed parameters {
    real<lower = 0> param1 = (mean_SI/sd_SI)^2;
    real<lower = 0> param2 = mean_SI/(sd_SI^2);

    vector<lower = min(S_L), upper = max(S_R)>[N] s;
    vector<lower = min(E_L), upper = max(E_R)>[N] e;

    s = S_L + (S_R - S_L) .* s_raw_g_n;
    for (k in 1:N) 
        if (E_R[k] > s[k]) 
            e[k] = E_L[k] + (s[k] - E_L[k]) * e_raw_g_n[k];
        else
            e[k] = E_L[k] + (E_R[k] - E_L[k]) * e_raw_g_n[k];
}

model {
    mean_SI ~ normal(5.0, 10.0);
    sd_SI ~ cauchy(0, 5.0);

    e_raw_g_n ~ normal(0.5, 1.0);
    s_raw_g_n ~ normal(0.5, 1.0);

    for (k in 1:N)
        target += gamma_lpdf(s[k] - e[k] | param1, param2);
}

generated quantities {
    vector[N] log_likelihood;
    for (k in 1:N) 
        log_likelihood[k] = gamma_lpdf(s[k] - e[k] | param1, param2);
}
"
)

gamma_n_fit<-stan(fit=gamma_n_model, data = data_input,
               warmup = 10000,
               iter = 20000,
               chains = 10,
               sample_file=standirname%&%'/sample',
               diagnostic_file=standirname%&%'/diagnostic',
               control=list(adapt_delta=0.98),
               cores=10
               )
print("gamma_n_fit\n")

log_lik<-extract_log_lik(gamma_n_fit, parameter_name = "log_likelihood", merge_chains = TRUE)
w<-waic(log_lik)
write.table("gamma",file = filenames,append=T)
write.table(w$estimates,file = filenames,append=T)



##Weibull distribution
## main dir for Stan simulations
standirname = stanmaindir%&%"/weibull-no_truncation"
unlink(standirname, recursive=T)
dir.create(standirname)


# Dumping initial conditions
e_raw_w_n = rep(.2, nrow(df))
s_raw_w_n = rep(.8, nrow(df))
param1 = 1.75

 

# Stan program
weibull_n_model <- stan(warmup = 100, data = data_input,
                    iter = 200,
                    chains = 0,
                    model_code = "
data {
    int<lower = 0> N; // number of records
    vector<lower = 0>[N] E_L;
    vector<lower = 0>[N] E_R;
    vector<lower = 0>[N] S_L;
    vector<lower = 0>[N] S_R;
}

parameters {
    real<lower = 0> mean_SI;
    real<lower = 0> param1;

    vector<lower = 0, upper = 1>[N] e_raw_w_n;
    vector<lower = 0, upper = 1>[N] s_raw_w_n;
}

transformed parameters {
    real<lower = 0> param2 = mean_SI/tgamma(1.0+1.0/param1);

    vector<lower = min(S_L), upper = max(S_R)>[N] s;
    vector<lower = min(E_L), upper = max(E_R)>[N] e;

    s = S_L + (S_R - S_L) .* s_raw_w_n;
    for (k in 1:N) 
        if (E_R[k] > s[k]) 
            e[k] = E_L[k] + (s[k] - E_L[k]) * e_raw_w_n[k];
        else
            e[k] = E_L[k] + (E_R[k] - E_L[k]) * e_raw_w_n[k];
}

model {
    mean_SI ~ normal(5.0, 10.0);
    param1 ~ exponential(0.0001);

    e_raw_w_n ~ normal(0.5, 1.0);
    s_raw_w_n ~ normal(0.5, 1.0);

    for (k in 1:N)
        target += weibull_lpdf(s[k] - e[k] | param1, param2);
}

generated quantities {
    real sd_SI = param2*sqrt(tgamma(1.0+2.0/param1)-(tgamma(1.0+1.0/param1))^2);

    vector[N] log_likelihood;
    for (k in 1:N) 
        log_likelihood[k] = weibull_lpdf(s[k] - e[k] | param1, param2);
}
"
)

weibull_n_fit<-stan(fit=weibull_n_model, data = data_input,
               warmup = 10000,
               iter = 20000,
               chains = 10,
               sample_file=standirname%&%'/sample',
               diagnostic_file=standirname%&%'/diagnostic',
               control=list(adapt_delta=0.98),
               cores=10
               )
print("weibull_n_fit\n")

log_lik<-extract_log_lik(weibull_n_fit, parameter_name = "log_likelihood", merge_chains = TRUE)
w<-waic(log_lik)
write.table("weibull",file = filenames,append=T)
write.table(w$estimates,file = filenames,append=T)

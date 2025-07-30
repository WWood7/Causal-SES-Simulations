source("data_generating_mechanism.r")


# calculate the true parameter values
data1 <- generate_data(1000000, 1)
es_smallv <- calc_nuisance_params_smallv(data1)
data2 <- generate_data(1000000, 2)
es_mediumv <- calc_nuisance_params_mediumv(data2)
data3 <- generate_data(1000000, 3)
es_largev <- calc_nuisance_params_largev(data3)
print(es_smallv)
print(es_mediumv)
print(es_largev)


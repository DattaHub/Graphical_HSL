#########################################################################
## Effective model size technique  of  Piironen  and  Vehtari  (2017) ###
#########################################################################
set.seed(123456789)

### Taylor series coefficients after expanding exp(-m^2*tan^2\eta) up to order of 16 ##
### in \eta about \eta = 0 ############################################################
coeff_vec = c(1,-1,-1/6,11/90,281/2520,4859/113400,10559/7484400,-7061209/681080400,
              -63274829/742996800)
coeff_vec = -1*coeff_vec
coeff_vec = coeff_vec[-1] ### left over coefficients before integration

#################### Integrating the polynomial in \eta #######################
temp_x = pi/2
int_vec = rep(0, length(coeff_vec))
for(i in 1:length(int_vec))
{
  int_vec[i] = (temp_x^(2*i+1))/(2*i+1)
}
coeff_vec = coeff_vec*int_vec
coeff_vec = coeff_vec/sqrt(pi)
coeff_vec ### non-zero c_r's 

##################### Preparing coefficient vector for solving roots ##########
##################### c_0, c_1, c_2, ... , c_31 ###############################
coef_vec_for_root_solving = rep(0, 32)

#### Fix p 
p = 100
coef_vec_for_root_solving[1] = 2/(p-1) - 1
j = 1
for(i in seq(3,31,4))
{
  coef_vec_for_root_solving[i+1] = coeff_vec[j]
  j = j+1
}
################## Solving the order 31 polynomial in m #######################
z = polyroot(coef_vec_for_root_solving)
as.matrix(z)

################ Extracting the real positive values of m #####################
Im_part_of_z = rep(0,length(z))
for(i in 1:length(z))
{
  Im_part_of_z[i] = as.numeric(Im(z[i]))
}
needed_roots_1 = which(abs(Im_part_of_z) < 1e-5)
needed_roots_2 = which(abs(Im_part_of_z) == 0)
needed_roots = c(needed_roots_1, needed_roots_2)
temp_real_part = rep(0,length(needed_roots))
for(i in 1:length(temp_real_part))
{
  temp_real_part[i] = Re(z[needed_roots[i]])
}
real_positive = which(temp_real_part>0)
m_req = rep(0, length(real_positive))
for(i in 1:length(m_req))
{
  m_req[i] = Re(z[needed_roots[real_positive[i]]])
}
m_req
########## Setting sample size and computing 'a' ##############################
## Assumption ## 
sigma_sq = 1  ##
################
sample_size = 120
a = (2*m_req^2*sigma_sq)/(sample_size)

############### The resulting global scale parameter ##########################
a

###############################################################################


# plot_privacy-rarefaction_curves.R
# R code to plot privacy-rarefaction curves for sex-specificity, i.e. to visualise the quantitative results made by privacy-rarefaction
# Copyright (C) 2017, ETH Zurich, Mathias Scharmann
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 		
# If you use this code please cite:
#
# "Scharmann M, Grafe TU, Metali F, Widmer A. (2017) Sex-determination 
# and sex chromosomes are shared across the radiation of dioecious 
# Nepenthes pitcher plants. XXX"
# 	
# contact: mathias.scharmann[-at-]env.ethz.ch or msph52[-at-]gmail.com



plot_candidate_drop <- function (indata, outname) {

attach(indata)

x=n_samples_per_sex

ym=log(male_specific_mean)
m_lower=log(male_specific_mean-male_specific_std)
m_upper=log(male_specific_mean+male_specific_std)

yf=log(female_specific_mean)
f_lower=log(female_specific_mean-female_specific_std)
f_upper=log(female_specific_mean+female_specific_std)

g = c(ym,yf)
ylim_lower = min( g[g!=-Inf] ) - 0.25 
ylim_upper = max( g[g!=-Inf] ) + 0.25 

pdf(sprintf("candidate_drop.%s.pdf", outname), width = 8, height = 4)

plot(0,type='n',
ylim=c(ylim_lower,ylim_upper),
xlim=c(1.0, indata$n_samples_per_sex[length( indata$n_samples_per_sex)] ),
xlab="number of individuals of each sex (stringency)",
ylab="log(sex-specific candidate count) +- SD",
main=outname
)

# shaded areas in background
# find highest point on X-axis for each M,F that is not zero, if any, then take the lower one of those.
# if no zeros for any sex, then fill entire plot area grey.
max_nonzero_idxes_m <- max(which(male_specific_mean != 0))
max_nonzero_idxes_f <- max(which(female_specific_mean != 0))
zero_idxes_m <- which(male_specific_mean == 0)
zero_idxes_f <- which(female_specific_mean == 0)
if (length(cbind(zero_idxes_m,zero_idxes_f)) > 0){
	lightgrey_xmax <- min(cbind(max_nonzero_idxes_m,max_nonzero_idxes_f)) + 0.5
} else {
# no upper limit
lightgrey_xmax <- indata$n_samples_per_sex[length( indata$n_samples_per_sex)] + 100
}

# find lowest point on X-axis at which M and F standard deviations are non-overlapping, if any
M_lower_greater_F_upper <- which((male_specific_mean-male_specific_std) > (female_specific_mean+female_specific_std))
F_lower_greater_M_upper <- which((female_specific_mean-female_specific_std) > (male_specific_mean+male_specific_std))

if (length(cbind(M_lower_greater_F_upper,F_lower_greater_M_upper)) > 0){
	darkgrey_xmax <- min(cbind(M_lower_greater_F_upper,F_lower_greater_M_upper)) - 0.5
} else if ( length(cbind(zero_idxes_m,zero_idxes_f)) > 0)  {
# overlap in STD ranges stops when mean of one sex drops to zero: darkgrey covers all of lightgrey
darkgrey_xmax <- lightgrey_xmax
} else {
# no upper limit
darkgrey_xmax <- indata$n_samples_per_sex[length( indata$n_samples_per_sex)] + 100
}
# xleft,ylower,xright,yupper
rect(0,ylim_lower-100,lightgrey_xmax,ylim_upper**2,lty=0,col = rgb(0.5,0.5,0.5,1/4))
rect(0,ylim_lower-100,darkgrey_xmax,ylim_upper**2,lty=0,col = rgb(0.5,0.5,0.5,1/4))

# add the data
arrows(x, ym , x, m_lower, length=0.05, angle=90, code=2)
arrows(x, ym , x, m_upper, length=0.05, angle=90, code=2)
lines(x, ym)
points(x, ym, pch=19, bg="black")

arrows(x, yf , x, f_lower, length=0.05, angle=90, code=2)
arrows(x, yf , x, f_upper, length=0.05, angle=90, code=2)
lines(x, yf)
points(x, yf, pch = 21, bg="white" )

legend("topright", legend=c("male-specific (Y-hemizygous)", "female-specific (W-hemizygous)"),  pch = c(19, 21), bg = "white")

#legend("bottomleft", title = "background shading", legend=c("no significant difference","significant difference","only one sex has specific markers"),  bg = "white", fill=c(rgb(0.5,0.5,0.5,2/4),rgb(0.5,0.5,0.5,1/4),"white"))


dev.off()

detach(indata)

}


##
args = commandArgs(trailingOnly=TRUE)

infile <- args[1]

mydata <- read.table( infile, sep = "\t", header = TRUE)

plot_candidate_drop(mydata, infile)

print("Done plotting!")
######

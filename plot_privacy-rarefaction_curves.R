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
ylim_lower = min( g[g!=-Inf] ) 
ylim_upper = max( g[g!=-Inf] ) 

pdf(sprintf("candidate_drop.%s.pdf", outname), width = 8, height = 4)

plot(x, ym, pch=19,
ylim=c(ylim_lower,ylim_upper),
xlim=c(1.0, indata$n_samples_per_sex[length( indata$n_samples_per_sex)] ),
xlab="stringency (number of samples per sex)",
ylab="log(sex-specific candidate count) +- SD",
main=outname
)

arrows(x, ym , x, m_lower, length=0.05, angle=90, code=2)
arrows(x, ym , x, m_upper, length=0.05, angle=90, code=2)
lines(x, ym)

#points(x, yf, pch = 21, bg="white" )
lines(x, yf)
#arrows(x, f_lower, x, f_upper, length=0.05, angle=90, code=3)

arrows(x, yf , x, f_lower, length=0.05, angle=90, code=2)
arrows(x, yf , x, f_upper, length=0.05, angle=90, code=2)
points(x, yf, pch = 21, bg="white" )

legend("topright", legend=c("male-specific (Y-hemizygous)", "female-specific (W-hemizygous)"),  pch = c(19, 21))
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

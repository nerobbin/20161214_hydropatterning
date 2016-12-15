#This script calculates growth-sustained tissue water potentials for a hypocotyl.
#It mainly serves as a proof-of-concept test of the method used for water potential estimates in root tissues,
#since empirical measurements of water potentials in hypocotyls are available in the literature
#(Nonami and Boyer, Plant Physiology 1993).

library(plyr)


#### Inputs and parameter values ####
#Read in growth table. 
#Column 1 is distance from hypocotyl apex in mm. 
#Column 2 is relative elemental growth rate in h^-1.
growth<-read.table("REGR.csv", header=TRUE, sep=",")

#Specify tissue layer radii in mm (measured from the center of the organ to the outer edge of the layer).
#Default dimensions were derived from positions of empirical measurements in Nonami and Boyer, 1993.
r_stem<-1.125 #Radius of entire organ
r_c1<-1.062 #Cortex layer 1
r_c2<-0.985 #Cortex layer 2
r_c3<-0.913 #Cortex layer 3
r_c4<-0.821 #Cortex layer 4
r_c5<-0.71 #Cortex layer 5
r_x<-0.567 #Xylem layer
r_p1<-0.430 #Pith layer 1
r_p2<-0.283 #Pith layer 2
r_p3<-0.110 #Pith layer 3

#Specify water potential of xylem water source in MPa.
W_xs<-(-0.06)

#Specify tissue water content as a fraction of the total root volume.
percent_water_content<-0.9

#Specify hydraulic conductivities in m^3 m^-2 s^-1 MPa^-1
L_stem_radial<-2.5*10^-7 #Tissue conductivity for radial water movement
L_stem_longitudinal<-2.5*10^-7 #Tissue conductivity for longitudinal water movement
L_xs<-1 #Xylem water source


#### Pre-processing of inputs ####
#Convert REGR values from growth curve to s^-1, so that units match with hydraulic conductivity units.
growth[,2]<-growth[,2]/3600

#Set segment height as the distance between two points on the growth curve.
#Note that this is a constant value, so the measurements on the input growth curve should be equally spaced.
h<-growth[2,1]-growth[1,1]

#Convert conductivity values to mm^3 mm^-2 s^-1 MPa^-1. 
#This helps to avoid errors when using the "solve()" function to solve the systems of equations below.
L_stem_radial<-L_stem_radial*10^3
L_stem_longitudinal<-L_stem_longitudinal*10^3
L_xs<-L_xs*10^3

#Calculate surface area of each compartment interface in mm^2.
if_c2c1<-(2*pi*r_c2)*h
if_c3c2<-(2*pi*r_c3)*h
if_c4c3<-(2*pi*r_c4)*h
if_c5c4<-(2*pi*r_c5)*h
if_xc5<-(2*pi*r_x)*h
if_xp1<-(2*pi*r_p1)*h
if_p1p2<-(2*pi*r_p2)*h
if_p2p3<-(2*pi*r_p3)*h
if_xsx<-(2*pi*(r_p1+2*(r_x-r_p1)))*h #Xylem source is a thin shell between the outer and inner edge of the xylem layer.

#Calculate radial-cross-sectional area of each compartment in mm^2.
cross_c1<-(pi*(r_c1^2))-(pi*(r_c2^2))
cross_c2<-(pi*(r_c2^2))-(pi*(r_c3^2))
cross_c3<-(pi*(r_c3^2))-(pi*(r_c4^2))
cross_c4<-(pi*(r_c4^2))-(pi*(r_c5^2))
cross_c5<-(pi*(r_c5^2))-(pi*(r_x^2))
cross_x<-(pi*(r_x^2))-(pi*(r_p1^2))
cross_p1<-(pi*(r_p1^2))-(pi*(r_p2^2))
cross_p2<-(pi*(r_p2^2))-(pi*(r_p3^2))
cross_p3<-(pi*(r_p3^2))

#Calculate interface conductivity (mm^3 s^-1 MPa^-1) as bulk conductivity multiplied by surface area of interface.
#Note that the source-to-hypocotyl value treats the source and tissue conductivities as being in series with each other.
L_c2c1<-L_stem_radial*if_c2c1
L_c3c2<-L_stem_radial*if_c3c2
L_c4c3<-L_stem_radial*if_c4c3
L_c5c4<-L_stem_radial*if_c5c4
L_xc5<-L_stem_radial*if_xc5
L_xp1<-L_stem_radial*if_xp1
L_p1p2<-L_stem_radial*if_p1p2
L_p2p3<-L_stem_radial*if_p2p3
L_xsx<-((L_stem_radial*L_xs)/(L_stem_radial+L_xs))*if_xsx
L_c1c1<-L_stem_longitudinal*cross_c1
L_c2c2<-L_stem_longitudinal*cross_c2
L_c3c3<-L_stem_longitudinal*cross_c3
L_c4c4<-L_stem_longitudinal*cross_c4
L_c5c5<-L_stem_longitudinal*cross_c5
L_xx<-L_stem_longitudinal*cross_x
L_p1p1<-L_stem_longitudinal*cross_p1
L_p2p2<-L_stem_longitudinal*cross_p2
L_p3p3<-L_stem_longitudinal*cross_p3


#### Flow rate calculation ####
###Notes
#The following sections calculate the rate of water movement (mm^3 s^-1) through each interface of the hypocotyl.
#This is done by treating the compartments as nodes in an electric circuit, and water flow as electric current.
#Thus, every current can be determined using Kirchhoff's circuit laws. These rules were used to write a system of
#linear equations that is then solved using matrix algebra to determine the unknown currents.
#The coefficients and solutions matrices for this (massive) system of equations are written in a loop, with
#each iteration adding the information for a single longitudinal segment of the root. The steps for the first and 
#last segments are performed slightly differently from intervening segments and are done outside of the loop.
#Once the two matrices are composed, the solve() function is used to determine the unknown flow rates.

###First longitudinal segment
#Generate coefficient matrix using values for first segment
coefficients_flows<-matrix(nrow=18,ncol=27,byrow=TRUE,c(
  1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  -1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,-1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,-1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,-1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,
  1/L_c2c1,0,0,0,0,0,0,0,0,1/L_c1c1,-1/L_c2c2,0,0,0,0,0,0,0,-1/L_c2c1,0,0,0,0,0,0,0,0,
  0,1/L_c3c2,0,0,0,0,0,0,0,0,1/L_c2c2,-1/L_c3c3,0,0,0,0,0,0,0,-1/L_c3c2,0,0,0,0,0,0,0,
  0,0,1/L_c4c3,0,0,0,0,0,0,0,0,1/L_c3c3,-1/L_c4c4,0,0,0,0,0,0,0,-1/L_c4c3,0,0,0,0,0,0,
  0,0,0,1/L_c5c4,0,0,0,0,0,0,0,0,1/L_c4c4,-1/L_c5c5,0,0,0,0,0,0,0,-1/L_c5c4,0,0,0,0,0,
  0,0,0,0,1/L_xc5,0,0,0,0,0,0,0,0,1/L_c5c5,-1/L_xx,0,0,0,0,0,0,0,-1/L_xc5,0,0,0,0,
  0,0,0,0,0,-1/L_xp1,0,0,0,0,0,0,0,0,1/L_xx,-1/L_p1p1,0,0,0,0,0,0,0,1/L_xp1,0,0,0,
  0,0,0,0,0,0,-1/L_p1p2,0,0,0,0,0,0,0,0,1/L_p1p1,-1/L_p2p2,0,0,0,0,0,0,0,1/L_p1p2,0,0,
  0,0,0,0,0,0,0,-1/L_p2p3,0,0,0,0,0,0,0,0,1/L_p2p2,-1/L_p3p3,0,0,0,0,0,0,0,1/L_p2p3,0,
  0,0,0,0,0,0,0,0,1/L_xsx,0,0,0,0,0,1/L_xx,0,0,0,0,0,0,0,0,0,0,0,-1/L_xsx))

#Assign column names.
colnames(coefficients_flows)<-c("c2c1 1","c3c2 1","c4c3 1","c5c4 1","xc5 1","xp1 1","p1p2 1","p2p3 1","xsx 1",
                                "c1c1 1","c2c2 1","c3c3 1","c4c4 1","c5c5 1","xx 1","p1p1 1","p2p2 1","p3p3 1",
                                "c2c1 2","c3c2 2","c4c3 2","c5c4 2","xc5 2","xp1 2","p1p2 2","p2p3 2","xsx 2")

#Generate solutions matrix using values for first segment
solutions_flows<-matrix(nrow=18,ncol=1,byrow=TRUE,c(
  (cross_c1*growth[1,2]*h*percent_water_content), (cross_c2*growth[1,2]*h*percent_water_content),
  (cross_c3*growth[1,2]*h*percent_water_content), (cross_c4*growth[1,2]*h*percent_water_content),
  (cross_c5*growth[1,2]*h*percent_water_content), (cross_x*growth[1,2]*h*percent_water_content),
  (cross_p1*growth[1,2]*h*percent_water_content), (cross_p2*growth[1,2]*h*percent_water_content),
  (cross_p3*growth[1,2]*h*percent_water_content), 0,0,0,0,0,0,0,0,0))

###Intermediary segments
for(i in 2:(nrow(growth)-1)){
  #Create matrices for coefficients and solutions for the current segment.
  coefficients_flows_slice<-matrix(nrow=18,ncol=36,byrow=TRUE, c(
    1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,1,0,0,0,0,0,0,0,-1,-1,0,0,1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,1/L_c2c1,0,0,0,0,0,0,0,0,1/L_c1c1,-1/L_c2c2,0,0,0,0,0,0,0,-1/L_c2c1,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,1/L_c3c2,0,0,0,0,0,0,0,0,1/L_c2c2,-1/L_c3c3,0,0,0,0,0,0,0,-1/L_c3c2,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,1/L_c4c3,0,0,0,0,0,0,0,0,1/L_c3c3,-1/L_c4c4,0,0,0,0,0,0,0,-1/L_c4c3,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,1/L_c5c4,0,0,0,0,0,0,0,0,1/L_c4c4,-1/L_c5c5,0,0,0,0,0,0,0,-1/L_c5c4,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,1/L_xc5,0,0,0,0,0,0,0,0,1/L_c5c5,-1/L_xx,0,0,0,0,0,0,0,-1/L_xc5,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/L_xp1,0,0,0,0,0,0,0,0,1/L_xx,-1/L_p1p1,0,0,0,0,0,0,0,1/L_xp1,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/L_p1p2,0,0,0,0,0,0,0,0,1/L_p1p1,-1/L_p2p2,0,0,0,0,0,0,0,1/L_p1p2,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/L_p2p3,0,0,0,0,0,0,0,0,1/L_p2p2,-1/L_p3p3,0,0,0,0,0,0,0,1/L_p2p3,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/L_xsx,0,0,0,0,0,1/L_xx,0,0,0,0,0,0,0,0,0,0,0,-1/L_xsx))
  colnames(coefficients_flows_slice)<-c(paste("c1c1",i-1), paste("c2c2",i-1), paste("c3c3",i-1), paste("c4c4",i-1),
    paste("c5c5",i-1), paste("xx",i-1), paste("p1p1",i-1), paste("p2p2",i-1), paste("p3p3",i-1), paste("c2c1",i),
    paste("c3c2",i), paste("c4c3",i), paste("c5c4",i), paste("xc5",i), paste("xp1",i), paste("p1p2",i),
    paste("p2p3",i), paste("xsx",i), paste("c1c1",i), paste("c2c2",i), paste("c3c3",i), paste("c4c4",i),
    paste("c5c5",i), paste("xx",i), paste("p1p1",i), paste("p2p2",i), paste("p3p3",i), paste("c2c1",i+1),
    paste("c3c2",i+1), paste("c4c3",i+1), paste("c5c4",i+1), paste("xc5",i+1), paste("xp1",i+1), paste("p1p2",i+1),
    paste("p2p3",i+1), paste("xsx",i+1))
  solutions_flows_slice<-matrix(nrow=18,ncol=1,byrow=TRUE,c(
    (cross_c1*growth[i,2]*h*percent_water_content), (cross_c2*growth[i,2]*h*percent_water_content),
    (cross_c3*growth[i,2]*h*percent_water_content), (cross_c4*growth[i,2]*h*percent_water_content),
    (cross_c5*growth[i,2]*h*percent_water_content), (cross_x*growth[i,2]*h*percent_water_content),
    (cross_p1*growth[i,2]*h*percent_water_content), (cross_p2*growth[i,2]*h*percent_water_content),
    (cross_p3*growth[i,2]*h*percent_water_content),0,0,0,0,0,0,0,0,0))
  
  #Merge the matrices for all segments with the ones from the current segment.
  coefficients_flows<-rbind.fill.matrix(coefficients_flows,coefficients_flows_slice)
  solutions_flows<-rbind(solutions_flows, solutions_flows_slice) 
}

###Last segment
#Create coefficents and solutions matrices for current segment, and append to matrices for all segments as before.
coefficients_flows_slice<-matrix(nrow=9,ncol=18,byrow=TRUE,c(
  1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
  0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,-1,1,0,0,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,0,-1,-1,0,0,1,
  0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,-1,0,0,
  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,-1,0,
  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0))
j<-nrow(growth)
colnames(coefficients_flows_slice)<-c(paste("c1c1",j-1), paste("c2c2",j-1), paste("c3c3",j-1), paste("c4c4",j-1),
  paste("c5c5",j-1), paste("xx",j-1), paste("p1p1",j-1), paste("p2p2",j-1), paste("p3p3",j-1), paste("c2c1",j),
  paste("c3c2",j), paste("c4c3",j), paste("c5c4",j), paste("xc5",j), paste("xp1",j), paste("p1p2",j),
  paste("p2p3",j), paste("xsx",j))
solutions_flows_slice<-matrix(nrow=9,ncol=1,byrow=TRUE,c(
  (cross_c1*growth[j,2]*h*percent_water_content), (cross_c2*growth[j,2]*h*percent_water_content),
  (cross_c3*growth[j,2]*h*percent_water_content), (cross_c4*growth[j,2]*h*percent_water_content),
  (cross_c5*growth[j,2]*h*percent_water_content), (cross_x*growth[j,2]*h*percent_water_content),
  (cross_p1*growth[j,2]*h*percent_water_content), (cross_p2*growth[j,2]*h*percent_water_content),
  (cross_p3*growth[j,2]*h*percent_water_content)))
coefficients_flows<-rbind.fill.matrix(coefficients_flows,coefficients_flows_slice)
solutions_flows<-rbind(solutions_flows, solutions_flows_slice)

#Convert all NA's to 0's in the coefficients table, and solve.
coefficients_flows[is.na(coefficients_flows)]<-0
variables_flows<-(solve(coefficients_flows,solutions_flows))

#### Water potential calculation####
###Notes
#Once all currents within the circuit are known, calculating the potential difference between each node is done
#by dividing the current by the resistance. This logic was used to write a second system of linear equations
#that is used to calculate unknown potentials of each node. As done previously for flow rates, coefficient
#and solution matrices are made in a loop. Solving the system of equations yields the water potentials in
#each compartment of the hypocotyl.

###First segment
#Generate coefficients matrix, with values for first segment.
coefficients_potentials<-matrix(nrow=10,ncol=10,byrow=TRUE, c(
  0,0,0,0,0,0,0,0,0,1,
  0,0,0,0,0,-1,0,0,0,1,
  0,0,0,0,-1,1,0,0,0,0,
  0,0,0,-1,1,0,0,0,0,0,
  0,0,-1,1,0,0,0,0,0,0,
  0,-1,1,0,0,0,0,0,0,0,
  -1,1,0,0,0,0,0,0,0,0,
  0,0,0,0,0,1,-1,0,0,0,
  0,0,0,0,0,0,1,-1,0,0,
  0,0,0,0,0,0,0,1,-1,0))
colnames(coefficients_potentials)<-c("c1 1","c2 1","c3 1","c4 1","c5 1","x 1","p1 1","p2 1","p3 1","xs")

#Generate solutions matrix with values for first segment.
solutions_potentials<-c(W_xs,variables_flows["xsx 1",]/L_xsx,variables_flows["xc5 1",]/L_xc5,
  variables_flows["c5c4 1",]/L_c5c4,variables_flows["c4c3 1",]/L_c4c3,variables_flows["c3c2 1",]/L_c3c2,
  variables_flows["c2c1 1",]/L_c2c1,variables_flows["xp1 1",]/L_xp1,variables_flows["p1p2 1",]/L_p1p2,
  variables_flows["p2p3 1",]/L_p2p3)

###Intermediary segments
for(i in 2:(nrow(growth))){
  #Create coefficients and solutions matrices for current segment.
  coefficients_potentials_slice<-matrix(nrow=9,ncol=18,byrow=TRUE,c(
    1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,
    0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,
    0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,
    0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,
    0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,
    0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1))
  colnames(coefficients_potentials_slice)<-c(paste("c1",i-1),paste("c2",i-1),paste("c3",i-1),paste("c4",i-1),
    paste("c5",i-1),paste("x",i-1),paste("p1",i-1),paste("p2",i-1),paste("p3",i-1),paste("c1",i),paste("c2",i),
    paste("c3",i),paste("c4",i),paste("c5",i),paste("x",i),paste("p1",i),paste("p2",i),paste("p3",i))
  
  #Append matrices for current segment to matrices for all segments.
  coefficients_potentials<-rbind.fill.matrix(coefficients_potentials,coefficients_potentials_slice)
  solutions_potentials<-c(solutions_potentials,c(variables_flows[paste("c1c1",i-1),]/L_c1c1,
    variables_flows[paste("c2c2",i-1),]/L_c2c2,variables_flows[paste("c3c3",i-1),]/L_c3c3,
    variables_flows[paste("c4c4",i-1),]/L_c4c4,variables_flows[paste("c5c5",i-1),]/L_c5c5,
    variables_flows[paste("xx",i-1),]/L_xx,variables_flows[paste("p1p1",i-1),]/L_p1p1,
    variables_flows[paste("p2p2",i-1),]/L_p2p2,variables_flows[paste("p3p3",i-1),]/L_p3p3))
}

#Convert all NA's to 0's in the coefficients table and solve for potentials.
coefficients_potentials[is.na(coefficients_potentials)]<-0
variables_potentials<-solve(coefficients_potentials,solutions_potentials)

#### Post-processing of output water potentials ####
#Drop the known xylem-source water potential from the outputs.
variables_potentials<-as.data.frame(variables_potentials[-10])

#Create matrix to store potentials organized by compartment.
potentials_output<-matrix(nrow=0,ncol=9)

#Transpose the values in variables_potentials and add them to potentials_output so that each column is one compartment.
for(i in 0:(nrow(growth)-1)){
  potentials_output<-rbind(potentials_output,t(variables_potentials[c((9*i)+1,9*i+2,9*i+3,9*i+4,9*i+5,9*i+6,9*i+7,9*i+8,9*i+9),]))
}

#Add growth curve data to the water potential matrix.
potentials_output<-cbind(growth,potentials_output)

#Convert REGR back to h^-1 in output water potential table.
potentials_output[,2]<-potentials_output[,2]*3600

#Add column names to output.
colnames(potentials_output)<-c("Distance (mm)","REGR (h^-1)","C1","C2","C3","C4","C5","X","P1","P2","P3")

#### Save output to file ####
write.table(potentials_output, "potentials_output.csv", row.names=FALSE, sep=",")

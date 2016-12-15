#This script calculates growth-sustained tissue water potentials for a root grown on agar media.

library(plyr)


#### Inputs and parameter values ####
#Read in growth table. 
#Column 1 is distance from root tip in mm. 
#Column 2 is relative elemental growth rate in h^-1.
growth<-read.table("REGR.csv", header=TRUE, sep=",")

#Specify tissue layer radii in mm (measured from the center of the organ to the outer edge of the layer).
r_root<-0.5 #Radius for entire organ
r_cortex<-(0.9*r_root) #Radius for cortex tissue layer
r_stele<-(0.4*r_root) #Radius for stele tissue layer
r_phloem<-(0.75*r_stele) #Radius for phloem water source

#Specify water potentials of water sources in MPa.
W_m1<-(-0.1) #Media 1 (contacts 25% of the root)
W_m2<-(-0.1) #Media 2 (contacts 50% of the root, e.g. air gap between two agar sheets)
W_m3<-(-0.1) #Media 3 (contacts 25% of the root)
W_p<-(-0.1) #Phloem

#Specify tissue water content as a fraction of the total root volume.
percent_water_content<-0.9

#Specify hydraulic conductivities in m^3 m^-2 s^-1 MPa^-1
L_root_radial<-1.15*10^-7 #Root conductivity for radial water movement
L_root_longitudinal<-1.15*10^-7 #Root conductivity for longitudinal water movement
L_m1<-2.5*10^-6 #Media 1
L_m2<-3.66*10^-9 #Media 2
L_m3<-3.66*10^-9 #Media 3
L_p<-2*10^-8 #Phloem


#### Pre-processing of inputs ####
#Convert REGR values from growth curve to s^-1, so that units match with hydraulic conductivity units.
growth[,2]<-growth[,2]/3600

#Convert conductivity values to mm^3 mm^-2 s^-1 MPa^-1. 
#This helps to avoid errors when using the "solve()" function to solve the systems of equations below.
L_root_radial<-L_root_radial*10^3
L_root_longitudinal<-L_root_longitudinal*10^3
L_m1<-L_m1*10^3
L_m2<-L_m2*10^3
L_m3<-L_m3*10^3
L_p<-L_p*10^3

#Create a matrix to store conductivity values for each compartment-interface in the root.
#These values are the bulk conductivity multiplied by the surface area of the interface.
LPs<-matrix(nrow=nrow(growth),ncol=21,NA)

#Column names are written as L_xy, with x referring to the sending compartment and y referring to the receiving compartment.
#When x and y are the same, this refers to longitudinal water movement.
#1, 2, and 3 denote root regions that are contacting medias 1, 2, and 3, respectively.
colnames(LPs)<-c("L_m1e1","L_m2e2","L_m3e3","L_ps","L_e1e2","L_e2e3","L_c1c2","L_c2c3",
  "L_e1c1","L_e2c2","L_e3c3","L_c1s","L_c2s","L_c3s","L_e1e1","L_e2e2","L_e3e3","L_c1c1","L_c2c2","L_c3c3","L_ss")

#Run a loop to fill in the interface conductivity table for each longitudinal segment of the root.
for(i in 1:nrow(growth)){
  #Specify height of the current root segment in mm.
  #For all segments except the last one, the height is the distance from the current to the next position on the growth curve.
  #Since the last segment has no "next position" to refer to, use the height of the second-to-last segment.
  ifelse(i!=nrow(growth), h<-growth[(i+1),1]-growth[i,1], h<-growth[(i),1]-growth[(i-1),1])
  
  #Calculate radial-cross-sectional area of each compartment type in mm^2.
  cross_e<-(0.25*((pi*r_root^2)-(pi*r_cortex^2))) #Epidermis compartment
  cross_c<-(0.25*((pi*r_cortex^2)-(pi*r_stele^2))) #Cortex compartment
  cross_s<-(pi*r_stele^2) #Stele compartment
  
  #Calculate interface conductivity (mm^3 s^-1 MPa^-1) as bulk conductivity multiplied by surface area of interface.
  #Note that source-to-root values treat the source and root conductivities as being in series with each other.
  LPs[i,"L_m1e1"]<-((L_root_radial*L_m1)/(L_root_radial+L_m1))*(2*pi*r_root*0.25)*h
  LPs[i,"L_m2e2"]<-((L_root_radial*L_m2)/(L_root_radial+L_m2))*(2*pi*r_root*0.5)*h
  LPs[i,"L_m3e3"]<-((L_root_radial*L_m3)/(L_root_radial+L_m3))*(2*pi*r_root*0.25)*h
  LPs[i,"L_ps"]<-((L_root_radial*L_p)/(L_root_radial+L_p))*(2*pi*r_phloem)*h
  LPs[i,"L_e1e2"]<-L_root_radial*2*(r_root-r_cortex)*h
  LPs[i,"L_e2e3"]<-L_root_radial*2*(r_root-r_cortex)*h
  LPs[i,"L_c1c2"]<-L_root_radial*2*(r_cortex-r_stele)*h
  LPs[i,"L_c2c3"]<-L_root_radial*2*(r_cortex-r_stele)*h
  LPs[i,"L_e1c1"]<-L_root_radial*(2*pi*r_cortex*0.25)*h
  LPs[i,"L_e2c2"]<-L_root_radial*(2*pi*r_cortex*0.5)*h
  LPs[i,"L_e3c3"]<-L_root_radial*(2*pi*r_cortex*0.25)*h
  LPs[i,"L_c1s"]<-L_root_radial*(2*pi*r_stele*0.25)*h
  LPs[i,"L_c2s"]<-L_root_radial*(2*pi*r_stele*0.5)*h
  LPs[i,"L_c3s"]<-L_root_radial*(2*pi*r_stele*0.25)*h
  LPs[i,"L_e1e1"]<-L_root_longitudinal*cross_e
  LPs[i,"L_e2e2"]<-L_root_longitudinal*2*cross_e
  LPs[i,"L_e3e3"]<-L_root_longitudinal*cross_e
  LPs[i,"L_c1c1"]<-L_root_longitudinal*cross_c
  LPs[i,"L_c2c2"]<-L_root_longitudinal*2*cross_c
  LPs[i,"L_c3c3"]<-L_root_longitudinal*cross_c
  LPs[i,"L_ss"]<-L_root_longitudinal*cross_s
}

#### Flow rate calculation ####
###Notes
#The following sections calculate the rate of water movement (mm^3 s^-1) through each interface of the root.
#This is done by treating the root compartments as nodes in an electric circuit, and water flow as electric current.
#Thus, every current can be determined using Kirchhoff's circuit laws. These rules were used to write a system of
#linear equations that is then solved using matrix algebra to determine the unknown currents.
#The coefficients and solutions matrices for this (massive) system of equations are written in a loop, with
#each iteration adding the information for a single longitudinal segment of the root. The steps for the first and 
#last segments are performed slightly differently from intervening segments and are done outside of the loop.
#Once the two matrices are composed, the solve() function is used to determine the unknown flow rates.

###First longitudinal segment
#Calculate segment height for the first segment.
h<-growth[2,1]-growth[1,1]

#Generate coefficient matrix using values for first segment
coefficients_flows<-matrix(nrow=21, ncol=35, byrow=TRUE, c(
  1,0,0,0,-1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,1,0,0,1,-1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,1,0,0,1,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,-1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,1,-1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,-1/LPs[1,"L_e1e2"],0,1/LPs[1,"L_c1c2"],0,1/LPs[1,"L_e1c1"],-1/LPs[1,"L_e2c2"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,-1/LPs[1,"L_e2e3"],0,1/LPs[1,"L_c2c3"],0,1/LPs[1,"L_e2c2"],-1/LPs[1,"L_e3c3"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,-1/LPs[1,"L_c1c2"],0,0,0,0,1/LPs[1,"L_c1s"],-1/LPs[1,"L_c2s"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,-1/LPs[1,"L_c2c3"],0,0,0,0,1/LPs[1,"L_c2s"],-1/LPs[1,"L_c3s"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  1/LPs[1,"L_m1e1"],-1/LPs[1,"L_m2e2"],0,0,1/LPs[1,"L_e1e2"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,1/LPs[1,"L_m2e2"],-1/LPs[1,"L_m3e3"],0,0,1/LPs[1,"L_e2e3"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,1/LPs[1,"L_m3e3"],-1/LPs[1,"L_ps"],0,0,0,0,0,0,1/LPs[1,"L_e3c3"],0,0,1/LPs[1,"L_c3s"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  1/LPs[1,"L_m1e1"],0,0,0,0,0,0,0,0,0,0,0,0,0,1/LPs[1,"L_e1e1"],0,0,0,0,0,0,-1/LPs[2,"L_m1e1"],0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,-1/LPs[1,"L_e1e2"],0,0,0,0,0,0,0,0,0,1/LPs[1,"L_e1e1"],-1/LPs[1,"L_e2e2"],0,0,0,0,0,0,0,0,0,1/LPs[2,"L_e1e2"],0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,-1/LPs[1,"L_e2e3"],0,0,0,0,0,0,0,0,0,1/LPs[1,"L_e2e2"],-1/LPs[1,"L_e3e3"],0,0,0,0,0,0,0,0,0,1/LPs[2,"L_e2e3"],0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,-1/LPs[1,"L_e1c1"],0,0,0,0,0,1/LPs[1,"L_e1e1"],0,0,-1/LPs[1,"L_c1c1"],0,0,0,0,0,0,0,0,0,0,0,1/LPs[2,"L_e1c1"],0,0,0,0,0,
  0,0,0,0,0,0,-1/LPs[1,"L_c1c2"],0,0,0,0,0,0,0,0,0,0,1/LPs[1,"L_c1c1"],-1/LPs[1,"L_c2c2"],0,0,0,0,0,0,0,0,1/LPs[2,"L_c1c2"],0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,-1/LPs[1,"L_c2c3"],0,0,0,0,0,0,0,0,0,0,1/LPs[1,"L_c2c2"],-1/LPs[1,"L_c3c3"],0,0,0,0,0,0,0,0,1/LPs[2,"L_c2c3"],0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[1,"L_c3s"],0,0,0,0,0,1/LPs[1,"L_c3c3"],-1/LPs[1,"L_ss"],0,0,0,0,0,0,0,0,0,0,0,0,0,1/LPs[2,"L_c3s"]))

#Assign column names.
colnames(coefficients_flows)<-c("m1e1 1","m2e2 1","m3e3 1","ps 1","e1e2 1","e2e3 1","c1c2 1","c2c3 1","e1c1 1",
  "e2c2 1", "e3c3 1", "c1s 1", "c2s 1", "c3s 1", "e1e1 1", "e2e2 1","e3e3 1","c1c1 1","c2c2 1","c3c3 1","ss 1",
  "m1e1 2","m2e2 2","m3e3 2","ps 2","e1e2 2","e2e3 2","c1c2 2","c2c3 2","e1c1 2","e2c2 2", "e3c3 2","c1s 2","c2s 2","c3s 2")

#Generate solutions matrix using values for first segment
solutions_flows<-matrix(nrow=21,ncol=1,byrow=TRUE,c(
  (cross_e*growth[1,2]*h*percent_water_content), (2*cross_e*growth[1,2]*h*percent_water_content), 
  (cross_e*growth[1,2]*h*percent_water_content), (cross_c*growth[1,2]*h*percent_water_content), 
  (2*cross_c*growth[1,2]*h*percent_water_content), (cross_c*growth[1,2]*h*percent_water_content), 
  (cross_s*growth[1,2]*h*percent_water_content),
  0,0,0,0,W_m1-W_m2,W_m2-W_m3,W_m3-W_p,0,0,0,0,0,0,0))

###Intermediary segments
for(i in 2:(nrow(growth)-1)){
  #Re-calculate slice height.
  h<-growth[(i+1),1]-growth[i,1]
  
  #Create matrices for coefficients and solutions for the current segment.
  coefficients_flows_slice<-matrix(nrow=21,ncol=42, byrow=TRUE, c(
    1,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,1,0,0,1,-1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,0,0,0,1,-1,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_e1e2"],0,1/LPs[i,"L_c1c2"],0,1/LPs[i,"L_e1c1"],-1/LPs[i,"L_e2c2"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_e2e3"],0,1/LPs[i,"L_c2c3"],0,1/LPs[i,"L_e2c2"],-1/LPs[i,"L_e3c3"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_c1c2"],0,0,0,0,1/LPs[i,"L_c1s"],-1/LPs[i,"L_c2s"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_c2c3"],0,0,0,0,1/LPs[i,"L_c2s"],-1/LPs[i,"L_c3s"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1/LPs[i,"L_m1e1"],-1/LPs[i,"L_m2e2"],0,0,1/LPs[i,"L_e1e2"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,1/LPs[i,"L_m2e2"],-1/LPs[i,"L_m3e3"],0,0,1/LPs[i,"L_e2e3"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,1/LPs[i,"L_m3e3"],-1/LPs[i,"L_ps"],0,0,0,0,0,0,1/LPs[i,"L_e3c3"],0,0,1/LPs[i,"L_c3s"],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1/LPs[i,"L_m1e1"],0,0,0,0,0,0,0,0,0,0,0,0,0,1/LPs[i,"L_e1e1"],0,0,0,0,0,0,-1/LPs[i+1,"L_m1e1"],0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_e1e2"],0,0,0,0,0,0,0,0,0,1/LPs[i,"L_e1e1"],-1/LPs[i,"L_e2e2"],0,0,0,0,0,0,0,0,0,1/LPs[i+1,"L_e1e2"],0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_e2e3"],0,0,0,0,0,0,0,0,0,1/LPs[i,"L_e2e2"],-1/LPs[i,"L_e3e3"],0,0,0,0,0,0,0,0,0,1/LPs[i+1,"L_e2e3"],0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_e1c1"],0,0,0,0,0,1/LPs[i,"L_e1e1"],0,0,-1/LPs[i,"L_c1c1"],0,0,0,0,0,0,0,0,0,0,0,1/LPs[i+1,"L_e1c1"],0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_c1c2"],0,0,0,0,0,0,0,0,0,0,1/LPs[i,"L_c1c1"],-1/LPs[i,"L_c2c2"],0,0,0,0,0,0,0,0,1/LPs[i+1,"L_c1c2"],0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_c2c3"],0,0,0,0,0,0,0,0,0,0,1/LPs[i,"L_c2c2"],-1/LPs[i,"L_c3c3"],0,0,0,0,0,0,0,0,1/LPs[i+1,"L_c2c3"],0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[i,"L_c3s"],0,0,0,0,0,1/LPs[i,"L_c3c3"],-1/LPs[i,"L_ss"],0,0,0,0,0,0,0,0,0,0,0,0,0,1/LPs[i+1,"L_c3s"]))
  colnames(coefficients_flows_slice)<-c(paste("e1e1",i-1), paste("e2e2",i-1), paste("e3e3",i-1), paste("c1c1",i-1), paste("c2c2",i-1),
    paste("c3c3", i-1), paste("ss", i-1), paste("m1e1", i), paste("m2e2", i), paste("m3e3", i), paste("ps", i), paste("e1e2", i), 
    paste("e2e3", i), paste("c1c2", i), paste("c2c3", i), paste("e1c1", i), paste("e2c2", i), paste("e3c3", i), paste("c1s", i), 
    paste("c2s", i), paste("c3s", i), paste("e1e1", i), paste("e2e2", i), paste("e3e3", i), paste("c1c1", i), paste("c2c2", i), 
    paste("c3c3", i), paste("ss", i), paste("m1e1", i+1), paste("m2e2", i+1), paste("m3e3", i+1), paste("ps", i+1), 
    paste("e1e2", i+1), paste("e2e3", i+1), paste("c1c2", i+1), paste("c2c3", i+1), paste("e1c1", i+1), paste("e2c2", i+1), 
    paste("e3c3", i+1), paste("c1s", i+1), paste("c2s", i+1), paste("c3s", i+1))
  solutions_flows_slice<-matrix(nrow=21,ncol=1, byrow=TRUE, c(
    (cross_e*growth[i,2]*h*percent_water_content), (2*cross_e*growth[i,2]*h*percent_water_content), 
    (cross_e*growth[i,2]*h*percent_water_content), (cross_c*growth[i,2]*h*percent_water_content), 
    (2*cross_c*growth[i,2]*h*percent_water_content), (cross_c*growth[i,2]*h*percent_water_content),
    (cross_s*growth[i,2]*h*percent_water_content), 0,0,0,0,W_m1-W_m2,W_m2-W_m3,W_m3-W_p,0,0,0,0,0,0,0))
  
  #Merge the matrices for all segments with the ones from the current segment.
  coefficients_flows<-rbind.fill.matrix(coefficients_flows, coefficients_flows_slice)
  solutions_flows<-rbind(solutions_flows, solutions_flows_slice)
}

###Last segment
#Height from the previous slice is re-used.
#Create coefficents and solutions matrices for current segment, and append to matrices for all segments as before.
j<-nrow(growth)
coefficients_flows_slice<-matrix(nrow=14,ncol=21, byrow=TRUE, c(
  1,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,-1,0,0,0,0,0,
  0,1,0,0,0,0,0,0,1,0,0,1,-1,0,0,0,-1,0,0,0,0,
  0,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,-1,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,-1,0,0,
  0,0,0,0,1,0,0,0,0,0,0,0,0,1,-1,0,1,0,0,-1,0,
  0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,-1,
  0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,-1/LPs[j,"L_e1e2"],0,1/LPs[j,"L_c1c2"],0,1/LPs[j,"L_e1c1"],-1/LPs[j,"L_e2c2"],0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[j,"L_e2e3"],0,1/LPs[j,"L_c2c3"],0,1/LPs[j,"L_e2c2"],-1/LPs[j,"L_e3c3"],0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[j,"L_c1c2"],0,0,0,0,1/LPs[j,"L_c1s"],-1/LPs[j,"L_c2s"],0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1/LPs[j,"L_c2c3"],0,0,0,0,1/LPs[j,"L_c2s"],-1/LPs[j,"L_c3s"],
  0,0,0,0,0,0,0,1/LPs[j,"L_m1e1"],-1/LPs[j,"L_m2e2"],0,0,1/LPs[j,"L_e1e2"],0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,1/LPs[j,"L_m2e2"],-1/LPs[j,"L_m3e3"],0,0,1/LPs[j,"L_e2e3"],0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,1/LPs[j,"L_m3e3"],-1/LPs[j,"L_ps"],0,0,0,0,0,0,1/LPs[j,"L_e3c3"],0,0,1/LPs[j,"L_c3s"]))
colnames(coefficients_flows_slice)<-c(paste("e1e1",j-1), paste("e2e2",j-1), paste("e3e3",j-1), paste("c1c1",j-1), paste("c2c2",j-1),
  paste("c3c3", j-1), paste("ss", j-1), paste("m1e1", j), paste("m2e2", j), paste("m3e3", j), paste("ps", j), paste("e1e2", j), 
  paste("e2e3", j), paste("c1c2", j), paste("c2c3", j), paste("e1c1", j), paste("e2c2", j), paste("e3c3", j), paste("c1s", j), 
  paste("c2s", j), paste("c3s", j))
solutions_flows_slice<-matrix(nrow=14,ncol=1, byrow=TRUE, c(
  (cross_e*growth[j,2]*h*percent_water_content), (2*cross_e*growth[j,2]*h*percent_water_content), 
  (cross_e*growth[j,2]*h*percent_water_content), (cross_c*growth[j,2]*h*percent_water_content), 
  (2*cross_c*growth[j,2]*h*percent_water_content), (cross_c*growth[j,2]*h*percent_water_content),
  (cross_s*growth[j,2]*h*percent_water_content), 0,0,0,0, W_m1-W_m2, W_m2-W_m3, W_m3-W_p))
coefficients_flows<-rbind.fill.matrix(coefficients_flows,coefficients_flows_slice)
solutions_flows<-rbind(solutions_flows, solutions_flows_slice)

#Convert all NA's to 0's in the coefficients table, and solve.
coefficients_flows[is.na(coefficients_flows)]<-0
variables_flows<-solve(coefficients_flows,solutions_flows)


#### Water potential calculation####
###Notes
#Once all currents within the circuit are known, calculating the potential difference between each node is done
#by dividing the current by the resistance. This logic was used to write a second system of linear equations
#that is used to calculate unknown potentials of each node. As done previously for flow rates, coefficient
#and solution matrices are made in a loop. Solving the system of equations yields the water potentials in
#each compartment of the root.

###First segment
#Generate coefficients matrix, with values for first segment.
coefficients_potentials<-matrix(nrow=11,ncol=11,byrow=TRUE, c(
  1,0,0,0,0,0,0,0,0,0,0,
  0,1,0,0,0,0,0,0,0,0,0,
  0,0,1,0,0,0,0,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,0,
  1,0,0,0,-1,0,0,0,0,0,0,
  0,1,0,0,0,-1,0,0,0,0,0,
  0,0,1,0,0,0,-1,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,-1,
  0,0,0,0,1,0,0,-1,0,0,0,
  0,0,0,0,0,1,0,0,-1,0,0,
  0,0,0,0,0,0,1,0,0,-1,0))
colnames(coefficients_potentials)<-c("m1","m2","m3","p","e1 1","e2 1","e3 1","c1 1","c2 1","c3 1","s 1")

#Generate solutions matrix with values for first segment.
solutions_potentials<-c(W_m1,W_m2,W_m3,W_p,variables_flows["m1e1 1",]/LPs[1,"L_m1e1"], variables_flows["m2e2 1",]/LPs[1,"L_m2e2"],
  variables_flows["m3e3 1",]/LPs[1,"L_m3e3"], variables_flows["ps 1",]/LPs[1,"L_ps"],variables_flows["e1c1 1",]/LPs[1,"L_e1c1"],
  variables_flows["e2c2 1",]/LPs[1,"L_e2c2"], variables_flows["e3c3 1",]/LPs[1,"L_e3c3"])

###Intermediary segments
for(i in 2:(nrow(growth))){
  #Create coefficients and solutions matrices for current segment.
  coefficients_potentials_slice<-matrix(nrow=7,ncol=14,byrow=TRUE,c(
  1,0,0,0,0,0,0,-1,0,0,0,0,0,0,
  0,1,0,0,0,0,0,0,-1,0,0,0,0,0,
  0,0,1,0,0,0,0,0,0,-1,0,0,0,0,
  0,0,0,1,0,0,0,0,0,0,-1,0,0,0,
  0,0,0,0,1,0,0,0,0,0,0,-1,0,0,
  0,0,0,0,0,1,0,0,0,0,0,0,-1,0,
  0,0,0,0,0,0,1,0,0,0,0,0,0,-1))
  colnames(coefficients_potentials_slice)<-c(paste("e1",i-1),paste("e2",i-1),paste("e3",i-1),paste("c1",i-1),
    paste("c2",i-1),paste("c3",i-1),paste("s",i-1),paste("e1",i),paste("e2",i),paste("e3",i),paste("c1",i),
    paste("c2",i),paste("c3",i),paste("s",i))
  
  #Append matrices for current segment to matrices for all segments.
  coefficients_potentials<-rbind.fill.matrix(coefficients_potentials,coefficients_potentials_slice)
  solutions_potentials<-c(solutions_potentials,c(variables_flows[paste("e1e1",i-1),]/LPs[i-1,"L_e1e1"],
    variables_flows[paste("e2e2",i-1),]/LPs[i-1,"L_e2e2"],variables_flows[paste("e3e3",i-1),]/LPs[i-1,"L_e3e3"],
    variables_flows[paste("c1c1",i-1),]/LPs[i-1,"L_c1c1"],variables_flows[paste("c2c2",i-1),]/LPs[i-1,"L_c2c2"],
    variables_flows[paste("c3c3",i-1),]/LPs[i-1,"L_c3c3"],variables_flows[paste("ss",i-1),]/LPs[i-1,"L_ss"]))
}

#Convert all NA's to 0's in the coefficients table and solve for potentials.
coefficients_potentials[is.na(coefficients_potentials)]<-0
variables_potentials<-solve(coefficients_potentials,solutions_potentials)

#### Post-processing of output water potentials ####
#Drop the first 4 values of the potentials, since they correspond to the water sources.
variables_potentials<-as.data.frame(variables_potentials[c(-1:-4)])

#Create matrix to store potentials organized by compartment.
potentials_output<-matrix(nrow=0,ncol=7)

#Transpose the values in variables_potentials and add them to potentials_output so that each column is one compartment.
for(i in 0:(nrow(growth)-1)){
  potentials_output<-rbind(potentials_output,t(variables_potentials[c((7*i)+1,7*i+2,7*i+3,7*i+4,7*i+5,7*i+6,7*i+7),]))
}

#Add growth curve data to the water potential matrix.
potentials_output<-cbind(growth, potentials_output)

#Convert REGR back to h^-1 in output object.
potentials_output[,2]<-potentials_output[,2]*3600

#Add column names to output.
#E, epidermis; C, cortex; S, stele
colnames(potentials_output)<-c("Distance (mm)","REGR (h^-1)","E1","E2","E3","C1","C2","C3","S")

#### Calculate percent phloem water ####
#The following code calculates the amount of water taken up from the phloem water source as a percentage of
#the total water taken up by the root.

#Create objects to store phloem influx and total influx.
phloem_influx<-0
total_influx<-0

#Calculate phloem and total influxes in each root segment in a loop.
#Any flow moving in the opposite direction (root to source) are ignored.
for(i in 1:nrow(growth)){
  #Add phloem influx for current segment to phloem influx for all segments, and to total influx.
  if(as.numeric(variables_flows[paste("ps",i),])>=0){
    phloem_influx<-phloem_influx+as.numeric(variables_flows[paste("ps",i),])
    total_influx<-total_influx+as.numeric(variables_flows[paste("ps",i),])
  }
  #Repeat for flows from media 1, 2, and 3.
  if(as.numeric(variables_flows[paste("m1e1",i),])>=0){
    total_influx<-total_influx+as.numeric(variables_flows[paste("m1e1",i),])
  }
  if(as.numeric(variables_flows[paste("m2e2",i),])>=0){
    total_influx<-total_influx+as.numeric(variables_flows[paste("m2e2",i),])
  }
  if(as.numeric(variables_flows[paste("m3e3",i),])>=0){
    total_influx<-total_influx+as.numeric(variables_flows[paste("m3e3",i),])
  }
}

#Calculate ratio of phloem influx to total influx.
percent_phloem<-phloem_influx/total_influx


#### Save output to file ####
write.table(potentials_output, "potentials_output.csv", row.names=FALSE, sep=",")
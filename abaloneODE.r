# Abalone ODE examples
# Original MatLab code by John Klinck (klinck@ccpo.odu.edu)
# Translated into R by Emil Aalto (aalto@cs.stanford.edu)

# Required packages
library(deSolve);
library(ggplot2);
library(grid);
library(reshape2);

# GLOBAL INDICES
iS = 1;
iI = 2;
iD = 3;
iP = 4;

# Top level function
# Run this with case=1,2a,2b,3a,3b,4
runModel <- function(case="1") {
  paramL = getAbaloneParams(case); #   define model parameters
  initL  = getInitialConditions(case, paramL);
  yInitV = initL$y;  # note: could be quite complex if vectorized array
  tSpanV = initL$t;
  
  # get appropriate function for the case
  caseFunc = getCaseFunction(case);

  # run the ODE solver
  resM = ode(yInitV, tSpanV, caseFunc, paramL, method="ode45");

  plotResults(case, resM, paramL);
}

# Handles all hard-coded params, and returns in a list
getAbaloneParams <- function(case="1") {
  # Overlap between cases 1, 2a, and 2b
  if ((case=="1")||(!is.na(charmatch("2", case)))) {
    # Parameters common to all cases
    depth = 1;   #%   vertical dimension of the volume (m)

    #% infection rate by contact with infections particles
    #%  [infected produced/infectious particle/day]
    IPinfect = 0.003 ;  

    #% infection rate by contact with infectious individuals 
    #%  [infected produced/susceptible animal/day]
    Iinfect = 0.001 ;  

    #% infection rate by contact with dead infectious individuals 
    #%  [infected produced/susceptible animal/day]
    Dinfect = 0.0008 ;  

    #% mortality rate of infected [1/day]
    Imort = 8*10^-2;  

    #% background mortality rate of susceptibles [1/day]
    Bmort = 0;  

    #% removal rate of dead [1/day]
    DeadDecay = 1.5 ;  

    #% infectious particles released by infected [particles/animal/day]
    Irelease  = .015;  

    #% infectious particles released by dead [particles/animal/day]
    Drelease = 1.0;  

    #% removal rate of IP from the environment [1/day]
    IPremove = 0.001;  

    baseL = list(case=case, depth=depth, IPinfect=IPinfect, Iinfect=Iinfect, Dinfect=Dinfect,
                Imort=Imort, Bmort=Bmort, DeadDecay=DeadDecay, Irelease=Irelease,
                Drelease=Drelease, IPremove=IPremove);
  
    caseL = list();
    if (!is.na(charmatch("2", case))) {
      Carry  = 100;
	  #  Carrying capacity in second population [number]
  	  Carry2 = 10;
	  #  Maximum reproduction rate in second population [1/day]
	  Repro2 = .01;

      if (case=="2a") {
	    #       exchange rate from second population  [1/day]
        Imm    = .0;     #% case 2a
	    #       diffusion rate from second population [1/day]
        Diff   = .0;     #% case 2a
	    #       Maximum reproduction rate for susceptable [number]
        ReproS = .02;   #% case 2a
	    #       Maximum reproduction rate for infected [number]
        ReproI = .01;   #% case 2a
      } else { # default
	    #       exchange rate from second population  [1/day]
        Imm    = .001;  #% case 2b
	    #       diffusion rate from second population [1/day]
        Diff   = .001;  #% case 2b
	    #       Maximum reproduction rate for susceptable [number]
        ReproS = .0;   #% case 2b
	    #       Maximum reproduction rate for infected [number]
        ReproI = .0;   #% case 2b
      }
	
	  caseL = list(Carry=Carry, Carry2=Carry2, Repro2=Repro2, 
	              Imm=Imm, Diff=Diff, ReproS=ReproS, ReproI=ReproI);
    }
	return(c(baseL, caseL));
  # separate params for 3a, 3b, and 4
  } else if (!is.na(charmatch("3", case))) return(getCase3Params(case))
  else if (case=="4") return(getCase4Params(case));
}

# Loads case 3a and 3b params
# Called by getAbaloneParams
getCase3Params <- function(case) {
  if (case=="3a") {
    # NOTE: Not tested with other values of Nx and Ny.
	# At the very least, needs new values for UexM and VexM.
	# Rest of the code may or may not work properly
	Nx = 3;   #   number of cross shore grid boxes (columns)
	Ny = 4;   #   number of along shore grid boxes (rows)
	Nvar = 4;
	
	oneM   = matrix(1, nrow=Ny, ncol=Nx);  
	zeroM  = 0*oneM;
	depthM = oneM #   vertical dimension of the volume (m)

	# infection rate by contact with infections particles
	#  [infected produced/infectious particle/day]
	IPinfectM = 0.025*oneM;  

	# infection rate by contact with infectious individuals 
	# [infected produced/susceptible animal/day]
	IinfectM = zeroM;
	#Iinfect = 0.02*oneM;

	#% infection rate by contact with dead infectious individuals 
	#%  [infected produced/susceptible animal/day]
	DinfectM = zeroM;
	#Dinfect = 0.0008*oneM;

	# mortality rate of infected [1/day]
	ImortM = 8*10^-2 * oneM;  

	# background mortality for susceptable [1/day]
	BmortM = zeroM;

	# removal rate of dead [1/day]
	DeadDecayM = 1.5*oneM;

	# infectious particles released by infected [particles/animal/day]
	IreleaseM  = .015*oneM;

	# infectious particles released by dead [particles/animal/day]
	DreleaseM = oneM;

	# removal rate of IP from the environment [1/day]
	IPremoveM = 0.001*oneM;

	# rate of transfer of IP in the E/W direction [1/day]
	#UexM = 0.05*oneM[,1:(Nx-1)]
	# NOTE: loaded the matrix in the MatLab style (row first) to match
	UexM = matrix(c(0.05, 0.05, 0, 0, 0, 0, -0.05, -0.05), nrow=Ny, ncol=(Nx-1), byrow=TRUE);  

	# rate of transfer of IP in the N/S direction [1/day]
	# VexM = 0.2*oneM[1:(Ny-1),];  
	VexM = matrix(c(-0.1, 0, 0.1, -0.1, 0, 0.1, -0.1, 0, 0.1), nrow=(Ny-1), ncol=Nx, byrow=TRUE);

    caseL = list(Nx=Nx, Ny=Ny, Nvar=Nvar, depthM=depthM, IPinfectM=IPinfectM, IinfectM=IinfectM, DinfectM=DinfectM,
                 ImortM=ImortM, BmortM=BmortM, DeadDecayM=DeadDecayM, IreleaseM=IreleaseM,
                 DreleaseM=DreleaseM, IPremoveM=IPremoveM, UexM=UexM, VexM=VexM);
  } else if (case=="3b") {
	# infection rate by contact with infections particles
	#  [infected produced/infectious particle/day]
	#  the convention is IPinfect(source species, infected species)

    # NOTE: Not tested with other values of Nspecies and Nvar.
	# At the very least, needs new values some of the other params.
	# Rest of the code may or may not work properly
	Nspecies  = 2;
	Nvar      = 4;
	oneM      = matrix(1, nrow=Nspecies, ncol=Nspecies);
	zeroM     = 0*oneM;
	IPinfectM = zeroM;
	#  particles from species 1 infect species 2
	IPinfectM[1,2] = 0.005;
	IPinfectM[2,1] = 0.0002;
	#IPinfectM[1,1] = 0.025;
	#IPinfectM[2,2] = 0.025;
	#IPinfectM = 0.025*oneM;

	# infection rate by contact with infectious individuals 
	#  [infected produced/susceptible animal/day]
	IinfectM = zeroM;
	#IinfectM = 0.02*oneM;

	# infection rate by contact with dead infected individuals 
	#  [infected produced/susceptible animal/day]
	DinfectM = zeroM;
	#DinfectM = 0.0008*oneM;

	singleM = matrix(1, nrow=Nspecies, ncol=1);
	#  reproduction rate for susceptibles
	SreproM = c(0.09, 0.05)*singleM;
	#  reproduction rate for infected
	IreproM = c(0.005, 0.005)*singleM;
	#   carrying capacity for the population
	CarryM  = c(150, 150)*singleM;

	# mortality rate of infected [1/day]
	ImortM  = 0.04*singleM;

	# background mortality rate of infected [1/day]
	#Bmort = .08*singleM;
	BmortM = 0*singleM;

	# removal rate of dead [1/day]
	DeadDecayM = 0.5*singleM;

	# infectious particles released by infected [particles/animal/day]
	#IreleaseM  = .015*singleM;
	IreleaseM  = c(0.001, 0.0001)*singleM;
	# infectious particles released by dead [particles/animal/day]
	DreleaseM = singleM;

	# removal rate of IP from the environment [1/day]
	#IPremove = 0.001*singleM;
	IPremoveM = c(0.1, 0.01)*singleM;

    caseL = list(Nspecies=Nspecies, Nvar=Nvar, IPinfectM=IPinfectM, IinfectM=IinfectM, DinfectM=DinfectM,
				 SreproM=SreproM, IreproM=IreproM, CarryM=CarryM,
				 ImortM=ImortM, BmortM=BmortM, DeadDecayM=DeadDecayM, IreleaseM=IreleaseM,
                 DreleaseM=DreleaseM, IPremoveM=IPremoveM);
  }

  return(caseL);
}

# Loads case 4 params
# Called by getAbaloneParams
getCase4Params <- function(case) {
  #   define model parameters and values
  #    some parameters are independent of space, others are arrays
  
  # NOTE: Not tested with other values of Nx, Ny, Nvar and Nspecies.
  # At the very least, needs new values for multiple params.
  # Rest of the code may or may not work properly
  Nvar     = 4;
  Nspecies = 2;
  Nx       = 3;
  Ny       = 4;

  onesM   = matrix(1, nrow=Nspecies, ncol=Nspecies);
  zeroM   = 0*onesM;
  singleM = matrix(1, nrow=Nspecies, ncol=1);

  depthM    = matrix(1,nrow=Ny, ncol=Nx);   #%   vertical dimension of the volume (m)

  #% infection rate by contact with infections particles
  #%  [infected produced/infectious particle/day]
  #%  the convention is IPinfect(source species, infected species)
  IPinfectM = zeroM;
  #%  particles from species 1 infect species 2
  IPinfectM[1,2] = 0.025;
  IPinfectM[2,1] = 0.001;
  #IPinfectM[1,1] = .003;
  #IPinfectM[2,2] = 0.025;
  #IPinfectM = 0.025*onesM;

  # infection rate by contact with infectious individuals 
  #  [infected produced/susceptible animal/day]
  IinfectM = zeroM;
  #IinfectM = 0.02*onesM;

  # infection rate by contact with dead infected individuals 
  #  [infected produced/susceptible animal/day]
  DinfectM = zeroM;
  #DinfectM = 0.0008*onesM;

  #%  reproduction rate for susceptibles
  SreproM = 0*singleM;
  SreproM[1,1] = .01;
  SreproM[2,1] = .01;
  #  reproduction rate for infected
  IreproM = 0*singleM;
  #IreproM[1,1] = .05;
  #IreproM[2,1] = .05;
  #   carrying capacity for the population
  CarryM = 0*singleM;
  CarryM[1,1] = 150;
  CarryM[2,1] = 100;

  # mortality rate of infected [1/day]
  ImortM = 0.08*singleM;

  # background mortality rate of infected [1/day]
  # BmortM = 0.08*singleM;
  BmortM = 0*singleM;

  # removal rate of dead [1/day]
  DeadDecayM = 1.5*singleM;

  # infectious particles released by infected [particles/animal/day]
  #IreleaseM  = .015*singleM;
  IreleaseM  = 0*singleM;
  IreleaseM[1,1] = .1;
  IreleaseM[2,1] = .01;
  # infectious particles released by dead [particles/animal/day]
  DreleaseM = singleM;

  # removal rate of IP from the environment [1/day]
  # PAR.IPremove = 0.001*ones(Nspecies,1);  
  IPremoveM = 0*singleM;
  IPremoveM[1,1] = 0.001;
  IPremoveM[2,1] = 0.001;

  # rate of transfer of IP in the E/W direction [1/day]
  #UexM = matrix(0, nrow=Ny, ncol=Nx-1);  
  #UexM = 0.05*matrix(1, nrow=Ny, ncol=Nx-1);  
  #UexM = [0.02 0.02 ;0 0; -0.02 -0.02 ];  
  UexM = matrix(c(0.02, 0.02, 0, 0, 0, 0, -0.02, -0.02), nrow=Ny, ncol=(Nx-1), byrow=TRUE);  

  # rate of transfer of IP in the N/S direction [1/day]
  #VexM = matrix(0, nrow=Ny-1, ncol=Nx);  
  #VexM = 0.2*matrix(1, nrow=Ny-1, ncol=Nx);  
  #VexM = [-0.02 0  0.02; -0.02 0  0.02];
  VexM = matrix(c(-0.02, 0, 0.02, -0.02, 0, 0.02, -0.02, 0, 0.02), nrow=(Ny-1), ncol=Nx, byrow=TRUE);

  caseL = list(Nx=Nx, Ny=Ny, Nvar=Nvar, Nspecies=Nspecies, depthM=depthM, 
               IPinfectM=IPinfectM, IinfectM=IinfectM, DinfectM=DinfectM, 
			   CarryM=CarryM, SreproM=SreproM, IreproM=IreproM, 
               ImortM=ImortM, BmortM=BmortM, DeadDecayM=DeadDecayM, IreleaseM=IreleaseM,
               DreleaseM=DreleaseM, IPremoveM=IPremoveM, UexM=UexM, VexM=VexM);
  return(caseL);
}

# Returns initial condition vector and time span vector in a list
# y: initial condition vector (could be vectorized array)
# t: time span vector
getInitialConditions <- function(case="1", paramL) {
  # initial conditions, time span to run
  # defaults from MatLab code
  if (case=="1") {
    yInitV = c(S=100, I=1, D=0, P=0); 
    tSpanV = seq(1, 100);  
  } else if (!is.na(charmatch("2", case))) {
    yInitV = c(S=100, I=1, D=0, P=0, S2=100); 
    tSpanV = seq(1, 100);  
  } else if (case=="3a") {
    # multi-cell array
	Ny   = paramL$Ny;
    Nx   = paramL$Nx;
	Nvar = paramL$Nvar;
	yInitA = array(0, dim=c(Ny, Nx, Nvar));
	yInitA[,,iS]   = 100;  # init S everywhere
	yInitA[1,1,iI] = 1;  # starting I in 1,1
	yInitA[1,1,iP] = 10; # starting P in 1,1
	# convert to vector
	dim(yInitA) <- c(Ny*Nx*Nvar);
	yInitV = yInitA;
    tSpanV = seq(1, 30);  #   time span
  } else if (case=="3b") {
    # two-species array
    Nspecies   = paramL$Nspecies;
	Nvar       = paramL$Nvar;
	yInitA = array(0, dim=c(Nspecies, Nvar));
	yInitA[,iS]  = 100;  # init S everywhere
	yInitA[1,iI] = 1;    # starting I in species 1
	# convert to vector
	dim(yInitA) <- c(Nspecies*Nvar);
	yInitV = yInitA;
    tSpanV = seq(1, 300);  #   time span
  } else if (case=="4") {
	Ny         = paramL$Ny;
    Nx         = paramL$Nx;
    Nspecies   = paramL$Nspecies;
	Nvar       = paramL$Nvar;
	yInitA = array(0, dim=c(Nvar, Nspecies, Ny, Nx));
	yInitA[iS,1,,] = 100;  # init S for species 1 everywhere
	yInitA[iS,2,,] = 50;   # init S for species 2 everywhere
	yInitA[iI,1,1,1] = 1; # init I for species 1 in cell 1,1
	# convert to vector
	dim(yInitA) <- c(Nspecies*Nvar*Ny*Nx);
	yInitV = yInitA;
    tSpanV = seq(1, 150);  #   time span
  }

  return(list(y=yInitV, t=tSpanV));
}

######################################
# ODE FUNCTIONS
######################################
# Returns appropriate ODE function for the case
getCaseFunction <- function(case="1") {
  if (case=="1") caseFunc = RHSabalone1
  else if (!is.na(charmatch("2", case))) caseFunc = RHSabalone2
  else if (case=="3a") caseFunc = RHSabalone3a
  else if (case=="3b") caseFunc = RHSabalone3b
  else if (case=="4") caseFunc = RHSabalone4
}

# ODE functions, maintaining MatLab names
RHSabalone1 <- function(time, stateV, paramL) {
  #  RHSabalone1
  #
  #    define the RHS of the model equations. 
  #        variables are 
  #  y(1): S, susceptible, uninfected individuals [number]
  #  y(2): I, infected individuals [number]
  #  y(3): D, dead animals from the infected population [number]
  #  y(4): P, infectious particles close 
  #            to the susceptible population [number]
 
  with(as.list(c(stateV, paramL)), {
    # stateV has named values here, so we can ignore index
    #   time changes of susceptable animals
    dS = - IPinfect*P*S - Iinfect*I*S - Dinfect*D*S - Bmort*S;

    #   time changes of infected animals
    dI = IPinfect*P*S + Iinfect*I*S + Dinfect*D*S - Imort*I;

    #   time changes of dead infected animals
    dD = Imort*I - DeadDecay*D;

    #   time changes of infectious particles
    dP = Irelease*I + Drelease*D - IPremove*P;

    return(list(c(dS, dI, dD, dP)));
  });
}

RHSabalone2 <- function(time, stateV, paramL) {
  #  RHSabalone2
  #
  #    define the RHS of the model equations. 
  #        variables are 
  #  y(1): S, susceptible, uninfected individuals [number]
  #  y(2): I, infected individuals [number]
  #  y(3): D, dead animals from the infected population [number]
  #  y(4): P, infectious particles close 
  #            to the susceptible population [number]
  #  y(5): S2, uninfected individuals in second population [number]
 
  with(as.list(c(stateV, paramL)), {
    # stateV is still named here (S,I,D,P)
    sDiff = S2 - S;
	if (sDiff<0) sDiff = 0;

	# calculate carrying capacity
    K = 1 - (S+I) / Carry;
    if (K<0) K = 0;
	
	# calculate carrying capacity for S2
	K2 = 1 - S2/Carry2;
	if (K2<0) K2 = 0;
	
    #   time changes of susceptible animals
    dS = - IPinfect*depth*P*S - Iinfect*I*S - Dinfect*D*S - Bmort*S +
	       ReproS*K*S + ReproI*K*I + Imm*S2 + Diff*sDiff;

    #   time changes of susceptable animals in 2nd population
    dS2 = Repro2*K2*S2 - Imm*S2 - Diff*sDiff;

    #   time changes of infected animals
    dI = IPinfect*P*S + Iinfect*I*S + Dinfect*D*S - Imort*I;

    #   time changes of dead infected animals
    dD = Imort*I - DeadDecay*D;

    #   time changes of infectious particles
    dP = Irelease*I + Drelease*D - IPremove*P;

    return(list(c(dS, dI, dD, dP, dS2)));
  });
}

RHSabalone3a <- function(time, stateV, paramL) {
  #  RHSabalone3a
  #
  #    define the RHS of the model equations. 
  #        variables are 
  #  y(1): S, susceptible, uninfected individuals [number]
  #  y(2): I, infected individuals [number]
  #  y(3): D, dead animals from the infected population [number]
  #  y(4): P, infectious particles close 
  #            to the susceptible population [number]
  
  with(as.list(c(stateV, paramL)), {
    # data comes in as a state vector
	# ODE solver doesn't handle arrays
	# convert into array for our convenience
	# preserving variable names from the MatLab code
	Ay    = array(data=stateV, dim=c(Ny,Nx,Nvar));
	# array for state t+1
    Aydot = array(data=0, dim=c(Ny,Nx,Nvar));

    for (ix in 1:Nx) {
      for (iy in 1:Ny) {
		# use "matrix" form of parameters
        #   time changes of susceptable animals
        Aydot[iy,ix,iS] = - IPinfectM[iy,ix] * depthM[iy,ix] * Ay[iy,ix,iP] * Ay[iy,ix,iS] -
						    IinfectM[iy,ix] * Ay[iy,ix,iI] * Ay[iy,ix,iS] -
						    DinfectM[iy,ix] * Ay[iy,ix,iD] * Ay[iy,ix,iS] -
						    BmortM[iy,ix] * Ay[iy,ix,iS];

        #   time changes of infected animals
        Aydot[iy,ix,iI] = IPinfectM[iy,ix] * depthM[iy,ix] * Ay[iy,ix,iP] * Ay[iy,ix,iS] + 
						  IinfectM[iy,ix] * Ay[iy,ix,iI] * Ay[iy,ix,iS] +
						  DinfectM[iy,ix] * Ay[iy,ix,iD] * Ay[iy,ix,iS] -
						  ImortM[iy,ix] * Ay[iy,ix,iI];
 
        #   time changes of dead infected animals
        Aydot[iy,ix,iD] = ImortM[iy,ix] * Ay[iy,ix,iI] -
						  DeadDecayM[iy,ix] * Ay[iy,ix,iD];

        #   time changes of infectious particles
        Aydot[iy,ix,iP] = IreleaseM[iy,ix] * Ay[iy,ix,iI] +
						  DreleaseM[iy,ix] * Ay[iy,ix,iD] -
						  IPremoveM[iy,ix] * Ay[iy,ix,iP];
      }
    }
    #             calculate the exchange of IP between grid cells
    #    Uex and Vex represent exchange rates (1/day)
    #    Uex(i,j) > 0 means particles move from cell(i,j) to cell(i+1,j)
    #    Uex(i,j) < 0 means particles move from cell(i+1,j) to cell(i,j)
    #    and similarly for Vex. 
    #        this will be Uex/Vex if it is positive and zero if negative.
    PUex = .5*(abs(UexM) + UexM); 
    PVex = .5*(abs(VexM) + VexM); 
    #%        this will be Uex/Vex if it is positive and zero if negative.
    NUex = .5*(abs(UexM) - UexM); 
    NVex = .5*(abs(VexM) - VexM); 

    #    E/W transfer
    for (ix in (1:(Nx-1))) {
      for (iy in 1:Ny) {
        Aydot[iy,ix,iP] = Aydot[iy,ix,iP] -
			PUex[iy,ix]*Ay[iy,ix,iP] + NUex[iy,ix]*Ay[iy,ix+1,iP];
        Aydot[iy,ix+1,iP] = Aydot[iy,ix+1,iP] +
			PUex[iy,ix]*Ay[iy,ix,iP] - NUex[iy,ix]*Ay[iy,ix+1,iP];
      }
    } 
    #   N/S transfer
    for (ix in 1:Nx) {
      for (iy in 1:(Ny-1)) {
        Aydot[iy,ix,iP] = Aydot[iy,ix,iP] -
			PVex[iy,ix]*Ay[iy,ix,iP] + NVex[iy,ix]*Ay[iy+1,ix,iP];
        Aydot[iy+1,ix,iP] = Aydot[iy+1,ix,iP] +
			PVex[iy,ix]*Ay[iy,ix,iP] - NVex[iy,ix]*Ay[iy+1,ix,iP];
      }
    } 

    # return array to vector format
    dim(Aydot) <- c(length(stateV));
    return(list(Aydot));
  });
}

RHSabalone3b <- function(time, stateV, paramL) {
  #  RHSabalone3a
  #
  #    define the RHS of the model equations. 
  #        variables are 
  #  y(1): S, susceptible, uninfected individuals [number]
  #  y(2): I, infected individuals [number]
  #  y(3): D, dead animals from the infected population [number]
  #  y(4): P, infectious particles close 
  #            to the susceptible population [number]
  
  with(as.list(c(stateV, paramL)), {
    # data comes in as a state vector
	# ODE solver doesn't handle arrays
	# convert into array for our convenience
	# preserving variable names from the MatLab code
	Ay    = array(data=stateV, dim=c(Nspecies,Nvar));
	# array for state t+1
    Aydot = array(data=0, dim=c(Nspecies,Nvar));

	#   calculate the carrying capacity term for each species
	kV = c();
	for (s in 1:Nspecies) {
	  K = 1 - (Ay[s,iS]+Ay[s,iI])/CarryM[s,1];
	  if (K<0) kV = c(kV, 0)
	  else kV = c(kV, K);
	}

	for (s in 1:Nspecies) {
	  #   time changes of susceptable animals
      Aydot[s,iS] = - BmortM[s,1] * Ay[s,iS] +
					  kV[s] * (SreproM[s,1]*Ay[s,iS] + IreproM[s,1]*Ay[s,iI]);
	  #   time changes of infected animals
	  Aydot[s,iI] = - ImortM[s,1] * Ay[s,iI];
	  #   add infection by different species
	  #   convention is Infect[source, infected]
	  for (i in 1:Nspecies) {
		Aydot[s,iS] = Aydot[s,iS] -
					  IPinfectM[i,s]*Ay[i,iP]*Ay[s,iS] -
					  IinfectM[i,s]*Ay[i,iI]*Ay[s,iS] -
					  DinfectM[i,s]*Ay[i,iD]*Ay[s,iS];
		Aydot[s,iI] = Aydot[s,iI] +
					  IPinfectM[i,s]*Ay[i,iP]*Ay[s,iS] +
					  IinfectM[i,s]*Ay[i,iI]*Ay[s,iS] +
					  DinfectM[i,s]*Ay[i,iD]*Ay[s,iS];
      }

	  #   time changes of dead infected animals
	  Aydot[s,iD] = ImortM[s,1]*Ay[s,iI] - DeadDecayM[s,1] * Ay[s,iD];

	  #   time changes of infectious particles
	  Aydot[s,iP] = IreleaseM[s,1]*Ay[s,iI] +
					DreleaseM[s,1]*Ay[s,iD] -
					IPremoveM[s,1]*Ay[s,iP];
    }

	# return array to vector format
    dim(Aydot) <- c(length(stateV));
    return(list(Aydot));
  });
}

RHSabalone4 <- function(time, stateV, paramL) {
  #  RHSabalone4
  #
  #    define the RHS of the model equations. 
  #        variables are 
  #  y(1): S, susceptible, uninfected individuals [number]
  #  y(2): I, infected individuals [number]
  #  y(3): D, dead animals from the infected population [number]
  #  y(4): P, infectious particles close 
  #            to the susceptible population [number]
  
  with(as.list(c(stateV, paramL)), {
    # data comes in as a state vector
	# ODE solver doesn't handle arrays
	# convert into array for our convenience
	# preserving variable names from the MatLab code
	Ay    = array(data=stateV, dim=c(Nvar, Nspecies, Ny, Nx));
	# array for state t+1
    Aydot = array(data=0, dim=c(Nvar, Nspecies, Ny, Nx));

    KK = array(1, dim=c(Nspecies,Ny,Nx));
    #   calculate the carrying capacity term for each area
    for (is in 1:Nspecies) {
      for (iy in 1:Ny) {
        for (ix in 1:Nx) { 
          KK[is,iy,ix] = 1 - (Ay[iS,is,iy,ix]+Ay[iI,is,iy,ix])/CarryM[is,1];
		  # check for zero here?
		  if (KK[is,iy,ix]<0) KK[is,iy,ix] = 0;
        }
      }
    }

    for (s in 1:Nspecies) {
      for (iy in 1:Ny) {
        for (ix in 1:Nx) {
          #   time changes of susceptable animals
          Aydot[iS,s,iy,ix] = - BmortM[s,1]*Ay[iS,s,iy,ix] +
                                KK[s,iy,ix]*(SreproM[s,1]*Ay[iS,s,iy,ix] +
                                IreproM[s,1]*Ay[iI,s,iy,ix]);
          #   time changes of infected animals
          Aydot[iI,s,iy,ix] = - ImortM[s,1] * Ay[iI,s,iy,ix];
          #   add infection by different species
          #   convention is Infect[source, infected]
          for (i in 1:Nspecies) {
             Aydot[iS,s,iy,ix] = Aydot[iS,s,iy,ix] -
								 IPinfectM[i,s]*depthM[iy,ix]*Ay[iP,i,iy,ix]*Ay[iS,s,iy,ix] -
								 IinfectM[i,s]*Ay[iI,i,iy,ix]*Ay[iS,s,iy,ix] -
								 DinfectM[i,s]*Ay[iD,i,iy,ix]*Ay[iS,s,iy,ix];
			 Aydot[iI,s,iy,ix] = Aydot[iI,s,iy,ix] +
								 IPinfectM[i,s]* depthM[iy,ix]*Ay[iP,i,iy,ix]*Ay[iS,s,iy,ix] +
								 IinfectM[i,s]*Ay[iI,i,iy,ix]*Ay[iS,s,iy,ix] +
								 DinfectM[i,s]*Ay[iD,i,iy,ix]*Ay[iS,s,iy,ix];
          }

		  #   time changes of dead infected animals
          Aydot[iD,s,iy,ix] = ImortM[s,1]*Ay[iI,s,iy,ix] - DeadDecayM[s,1]*Ay[iD,s,iy,ix];

		  #   time changes of infectious particles
		  Aydot[iP,s,iy,ix] = IreleaseM[s,1]*Ay[iI,s,iy,ix] +
							  DreleaseM[s,1]*Ay[iD,s,iy,ix] -
							  IPremoveM[s,1]*Ay[iP,s,iy,ix];
        }
      }
	  #             calculate the exchange of IP between grid cells
	  #    Uex and Vex represent exchange rates (1/day)
	  #    Uex[i,j] > 0 means particles move from cell[i,j] to cell[i+1,j]
	  #    Uex[i,j] < 0 means particles move from cell[i+1,j] to cell[i,j]
	  #    and similarly for Vex. 
	  #        this will be Uex/Vex if it is positive and zero if negative.
      PUexM = .5*(abs(UexM) + UexM); 
      PVexM = .5*(abs(VexM) + VexM); 
	  #        this will be Uex/Vex if it is positive and zero if negative.
      NUexM = .5*(abs(UexM) - UexM); 
      NVexM = .5*(abs(VexM) - VexM); 

      #    E/W transfer
      for (ix in 1:(Nx-1)) {
        for (iy in 1:Ny) {
          Aydot[iP,s,iy,ix] = Aydot[iP,s,iy,ix] - 
		                      PUexM[iy,ix]*Ay[iP,s,iy,ix] + NUexM[iy,ix]*Ay[iP,s,iy,ix+1];
          Aydot[iP,s,iy,ix+1] = Aydot[iP,s,iy,ix+1] +
							  PUexM[iy,ix]*Ay[iP,s,iy,ix] - NUexM[iy,ix]*Ay[iP,s,iy,ix+1];
        }
      }
      #%   N/S transfer
      for (ix in 1:Nx) {
        for (iy in 1:(Ny-1)) {
          Aydot[iP,s,iy,ix] = Aydot[iP,s,iy,ix] -
                              PVexM[iy,ix]*Ay[iP,s,iy,ix] + NVexM[iy,ix]*Ay[iP,s,iy+1,ix];
          Aydot[iP,s,iy+1,ix] = Aydot[iP,s,iy+1,ix] +
                              PVexM[iy,ix]*Ay[iP,s,iy,ix] - NVexM[iy,ix]*Ay[iP,s,iy+1,ix];
        }
      }
    }

	# return array to vector format
    dim(Aydot) <- c(length(stateV));
    return(list(Aydot));
  });
}

##################################################
# PLOT FUNCTIONS
##################################################
# Convenience function
# Plots to a generated window, and allows multiple subplots
myPPlot <- function(thePlot, w=16, h=12, plotPos=c(1,1), newWindow=TRUE, plotDim=c(1,1), invertY=FALSE) {
  if (newWindow) {
    windows(w,h); 
    grid.newpage();
    pushViewport(viewport(layout=grid.layout(plotDim[1],plotDim[2])));
  }
  
  plotRow = plotPos[1];
  if (invertY) {
    # invert row order to match MatLab output (y=1 at the bottom, not the top)
    if (plotDim[1]>1) plotRow = plotDim[1] - plotRow + 1;
  }
  print(thePlot, vp=viewport(layout.pos.row=plotRow, layout.pos.col=plotPos[2]));
}

# Produces matching plots to MatLab output, based on case
# Calls appropriate sub-function
plotResults <- function(case, resM, paramL) {

  if (case=="1") plotCase1(case, resM, paramL)
  else if (!is.na(charmatch("2", case))) plotCase2(case, resM, paramL)
  else if (case=="3a") plotCase3a(case, resM, paramL)
  else if (case=="3b") plotCase3b(case, resM, paramL)
  else if (case=="4") plotCase4(case, resM, paramL);
}
  
plotCase1 <- function(case, resM, paramL) {  
  # reform matrix as melted data
  resM = resM[,-1]; # drop the time column
  resD = melt(resM, varnames=c("t", "category"), value.name="count");

  # plot four states over time
  titleV = c("Susceptible", "Infected", "Dead Infected", "Infect Particles");
  catV   = c("S", "I", "D", "P");
  posV   = list("S"=c(1,1), "I"=c(1,2), "D"=c(2,1), "P"=c(2,2));
  for (cat in c(iS, iI, iD, iP)) {
	# get data for specific category
	catD    = resD[resD$category==catV[cat],];
	if (cat==iP) yStr = "number/m3"
	else yStr = "number/m2";
    thePlot = ggplot(data=catD, aes(x=t, y=count));
    thePlot = thePlot + geom_line() + theme_bw() +
	                    labs(title=titleV[cat], x="", y=yStr) +
	                    theme(legend.position="none");
	  
	# plot graph in grid (opening new window for first)
	plotPos = posV[[catV[cat]]];
    myPlot(thePlot, plotPos=plotPos, newWindow=(cat==iS), plotDim=c(2,2));
  }
	
  # plot infection rate over time
  # easier to get from unmelted data
  rateD = data.frame(t=1:nrow(resM));
  rateD$pRate = paramL$IPinfect * resM[,iS] * resM[,iP];
  rateD$iRate = paramL$Iinfect * resM[,iS] * resM[,iI];
  rateD$dRate = paramL$Dinfect * resM[,iS] * resM[,iD];
  # melt again
  meltRateD = melt(as.matrix(rateD[,-1]), varnames=c("t", "rate"), value.name="daily");
	
  thePlot = ggplot(data=meltRateD, aes(x=t, y=daily));
  thePlot = thePlot + geom_line(aes(group=rate, linetype=rate)) +
					scale_linetype_manual(values=c("pRate"="solid", "iRate"="dashed", "dRate"="dotted"), 
										  breaks=c("pRate", "iRate", "dRate"),
										  labels=c("pRate"="P", "iRate"="I", "dRate"="D"));
	                    
  thePlot = thePlot + labs(title="Infection Rate (num/day)", x="day", y="number/day") +
	                  theme(legend.position=c(0.8, 0.5)) +
	  				  theme_bw();
  myPlot(thePlot);					
}	

plotCase2 <- function(case, resM, paramL) {
  # reform matrix as melted data
  resM = resM[,-1]; # drop the time column
  resD = melt(resM, varnames=c("t", "category"), value.name="count");

  # plot four states over time (S2 added to S)
  titleV = c("Susceptible", paste("Infected", case), "Dead Infected", "Infect Particles");
  catV   = c("S", "I", "D", "P");
  posV   = list("S"=c(1,1), "I"=c(1,2), "D"=c(2,1), "P"=c(2,2));
  for (cat in c(iS, iI, iD, iP)) {
	# get data for specific category
	catD    = resD[resD$category==catV[cat],];
	if (cat==iP) yStr = "number/m3"
	else yStr = "number/m2";
    thePlot = ggplot(data=catD, aes(x=t, y=count));
    thePlot = thePlot + geom_line() + theme_bw() +
	                    labs(title=titleV[cat], x="", y=yStr) +
	                    theme(legend.position="none");
	if (cat==iS) thePlot = thePlot + geom_line(data=resD[resD$category=="S2",], aes(linetype="dashed"));
	  
	# plot graph in grid (opening new window for first)
	plotPos = posV[[catV[cat]]];
    myPlot(thePlot, plotPos=plotPos, newWindow=(cat==iS), plotDim=c(2,2));
  }
	
  # plot infection rate over time
  # easier to get from unmelted data
  rateD = data.frame(t=1:nrow(resM));
  rateD$pRate = paramL$IPinfect * resM[,iS] * resM[,iP];
  rateD$iRate = paramL$Iinfect * resM[,iS] * resM[,iI];
  rateD$dRate = paramL$Dinfect * resM[,iS] * resM[,iD];
  # melt again (easier to plot from melted data)
  meltRateD = melt(as.matrix(rateD[,-1]), varnames=c("t", "rate"), value.name="daily");
	
  thePlot = ggplot(data=meltRateD, aes(x=t, y=daily));
  thePlot = thePlot + geom_line(aes(group=rate, linetype=rate)) +
	                  scale_linetype_manual(values=c("pRate"="solid", "iRate"="dashed", "dRate"="dotted"), 
						 				    breaks=c("pRate", "iRate", "dRate"),
											labels=c("pRate"="P", "iRate"="I", "dRate"="D"));
	                    
  thePlot = thePlot + labs(title="Infection Rate (num/day)", x="day", y="number/day") +
	                  theme(legend.position=c(0.8, 0.5)) +
	 				  theme_bw();
  myPlot(thePlot);				

  # plot reproduction
  iS2 = 5;
  rateD$K = 1 - (resM[,iS] + resM[,iI]) / paramL$Carry;
  if (nrow(rateD[rateD$K<0,])>0) rateD[rateD$K<0,]$K = 0;
  rateD$repS = paramL$ReproS * rateD$K * resM[,iS];
  rateD$repI = paramL$ReproI * rateD$K * resM[,iI];

  rateD$K2 = 1 - resM[,iS2] / paramL$Carry2;
  if (nrow(rateD[rateD$K2<0,])>0) rateD[rateD$K2<0,]$K2 = 0;
  # ERROR in MatLab: should be K2 here, not K
  rateD$rep2 = paramL$Repro2 * rateD$K2 * resM[,iS2];
  # melt again
  meltRepoD = melt(as.matrix(rateD[,c("repS", "repI", "rep2")]), varnames=c("t", "rep"));

  thePlot = ggplot(data=meltRepoD, aes(x=t, y=value));
  thePlot = thePlot + geom_line(aes(group=rep, linetype=rep)) +
					  scale_linetype_manual(values=c("repS"="solid", "repI"="dashed", "rep2"="dotted"), 
						  				    breaks=c("repS", "repI", "rep2"),
											labels=c("repS"="S", "repI"="I", "rep2"="S2"));
	                    
  thePlot = thePlot + labs(title="Reproduction Rate (num/day)", x="day", y="number/day") +
	                  theme(legend.position=c(0.8, 0.5)) +
			  		  theme_bw();
  myPlot(thePlot);		
	
  # plot immigration rate
  rateD$D = resM[,iS2] - resM[,iS];
  if (nrow(rateD[rateD$D<0,])>0) rateD[rateD$D<0,]$D = 0;
  rateD$inpI = paramL$Imm * resM[,iS2];
  rateD$inpD = paramL$Diff * rateD$D;
  # melt again
  meltInpD = melt(as.matrix(rateD[,c("inpI", "inpD")]), varnames=c("t", "inp"));

  thePlot = ggplot(data=meltInpD, aes(x=t, y=value));
  thePlot = thePlot + geom_line(aes(group=inp, linetype=inp)) +
				  	  scale_linetype_manual(values=c("inpI"="solid", "inpD"="dashed"), 
						  				    breaks=c("inpI", "inpD"),
											labels=c("inpI"="Im", "inpD"="Diff"));
	                    
  thePlot = thePlot + labs(title="Immigration Rate (num/day)", x="day", y="number/day") +
	                  theme(legend.position=c(0.8, 0.5)) +
	  				  theme_bw();
  myPlot(thePlot);				
}

plotCase3a <- function(case, resM, paramL) {
  resM = resM[,-1]; # drop the time column
  Nt   = nrow(resM);
  Ny   = paramL$Ny;
  Nx   = paramL$Nx;
  Nvar = paramL$Nvar;
  resA = array(resM, dim=c(Nt, Ny, Nx, Nvar));
  resD = melt(resA, varnames=c("t", "y", "x", "category"), value.name="count");
  resD$cell = paste(resD$y, resD$x, sep="_");
	
  # plot four states over time (multiple lines per graph)
  titleV = c("Susceptible", paste("Infected", case), "Dead Infected", "Inf Particles");
  catV   = c("S", "I", "D", "P");
  posV   = list("S"=c(1,1), "I"=c(1,2), "D"=c(2,1), "P"=c(2,2));
  for (cat in c(iS, iI, iD, iP)) {
	# get data for specific category
	# Categories now saved as 1-4, not S-P
	catD    = resD[resD$category==cat,];
	if (cat==iP) yStr = "number/m3"
	else yStr = "number/m2";
    thePlot = ggplot(data=catD, aes(x=t, y=count));
    thePlot = thePlot + geom_line(aes(group=cell, color=cell)) + theme_bw() +
	                    labs(title=titleV[cat], x="", y=yStr) +
	                    theme(legend.position="none");
	  
	# plot graph in grid (opening new window for first)
	plotPos = posV[[catV[cat]]];
    myPlot(thePlot, plotPos=plotPos, newWindow=(cat==iS), plotDim=c(2,2));
  }
  
  # plot three states over time for each cell in a grid
  titleV = c("Susceptible", paste("Infected", case), "Dead Infect", "Infect Particles");
  catV   = c("S", "I", "P");
  for (cat in c(iS, iI, iP)) {
	newWindow = TRUE;
	catD  = resD[resD$category==cat,];
	yLimV = c(0, max(catD$count));
    for (iy in 1:Ny) {
	  # inverting y later
	  if (iy==Ny) titleStr = titleV[cat]
	  else titleStr = "";
	  if (iy==1) xStr = "Days"
	  else xStr = "";
	  for (ix in 1:Nx) {
   	    # get data for specific cell and category
	    # Categories now saved as 1-4, not S-P
		cellD   = catD[catD$cell==paste(iy, ix, sep="_"),];
        thePlot = ggplot(data=catD, aes(x=t, y=count));
        thePlot = thePlot + geom_line() + theme_bw() +
	                        labs(title=titleStr, x=xStr, y="") +
							scale_y_continuous(limits=yLimV) +
	                        theme(legend.position="none");
	  
	    # plot graph in grid (opening new window for first)
	    plotPos = c(iy, ix);
		# make sure to invert the cell Y value to match MatLab (y=1 at the bottom)
        myPlot(thePlot, plotPos=plotPos, newWindow=newWindow, plotDim=c(Ny,Nx), invertY=TRUE);
		newWindow = FALSE;
	  }
	}
  }
}

plotCase3b <- function(case, resM, paramL) {
  # reform matrix as melted data
  resM = resM[,-1]; # drop the time column
  Nt   = nrow(resM);
  Nspecies = paramL$Nspecies;
  Nvar = paramL$Nvar; 
  resA = array(resM, dim=c(Nt, Nspecies, Nvar));
  resD = melt(resA, varnames=c("t", "species", "category"), value.name="count");
  resD$species = as.character(resD$species);

  # plot four states over time
  titleV = c("Susceptible", "Infected", "Dead Infected", "Infect Particles");
  catV   = c("S", "I", "D", "P");
  posV   = list("S"=c(1,1), "I"=c(1,2), "D"=c(2,1), "P"=c(2,2));
  for (cat in c(iS, iI, iD, iP)) {
    # get data for specific category
	catD    = resD[resD$category==cat,];
	if (cat==iP) yStr = "number/m3"
	else yStr = "number/m2";
    thePlot = ggplot(data=catD, aes(x=t, y=count));
    thePlot = thePlot + geom_line(aes(group=species, color=species)) + theme_bw() +
	                    labs(title=titleV[cat], x="", y=yStr) +
				        scale_color_manual(values=c("1"="black", "2"="blue"));
	                    theme(legend.position="none");
	if (cat==iS) thePlot = thePlot + theme(legend.position=c(0.8, 0.8))
	else thePlot = thePlot + theme(legend.position="none");
	  
	# plot graph in grid (opening new window for first)
	plotPos = posV[[catV[cat]]];
    myPlot(thePlot, plotPos=plotPos, newWindow=(cat==iS), plotDim=c(2,2));
  }
}

plotCase4 <- function(case, resM, paramL) {
  # reform matrix as melted data
  resM     = resM[,-1]; # drop the time column
  Nt       = nrow(resM);
  Nspecies = paramL$Nspecies;
  Nvar     = paramL$Nvar; 
  Ny       = paramL$Ny;
  Nx       = paramL$Nx;
  resA = array(resM, dim=c(Nt, Nvar, Nspecies, Ny, Nx));
  resD = melt(resA, varnames=c("t", "category", "species", "y", "x"), value.name="count");
  resD$species = as.character(resD$species);
  resD$cell    = paste(resD$y, resD$x, sep="_");
  
  # plot four states over time
  titleV = c("Susceptible", "Infected", "Dead Infected", "Inf Particles");
  catV   = c("S", "I", "D", "P");
  posV   = list("S"=c(1,1), "I"=c(1,2), "D"=c(2,1), "P"=c(2,2));
  for (cat in c(iS, iI, iD, iP)) {
    # get data for specific category
	catD    = resD[resD$category==cat,];
	if (cat==iP) yStr = "number/m3"
	else yStr = "number/m2";
	for (s in 1:Nspecies) {
	  speciesD = catD[catD$species==s,];
  	  if (s==1) {
  	    # species 1
	    thePlot = ggplot(data=speciesD, aes(x=t, y=count));
        thePlot = thePlot + geom_line(aes(group=cell, color=cell)) + theme_bw() +
	                        labs(title=titleV[cat], x="", y=yStr) +
	                        theme(legend.position="none");
	  } else {
	    # other species
		thePlot = thePlot + geom_line(data=speciesD, aes(group=cell, color=cell));
      } 
    }	  
	# plot graph in grid (opening new window for first)
	plotPos = posV[[catV[cat]]];
    myPlot(thePlot, plotPos=plotPos, newWindow=(cat==iS), plotDim=c(2,2));
  }
  
  # plot three states over time for each cell in a grid
  titleV = c("Susceptible", paste("Infected", case), "Dead Infect", "Inf Particles");
  catV   = c("S", "I", "P");
  for (cat in c(iS, iI, iP)) {
	newWindow = TRUE;
	titleStr  = titleV[cat];
	catD      = resD[resD$category==cat,];
    yLimV     = c(0, max(catD$count));
    for (iy in 1:Ny) {
	  # we are inverting Y later
	  if (iy==Ny) titleStr = titleV[cat]
	  else titleStr = "";
	  if (iy==1) xStr = "Days"
	  else xStr = "";
	  for (ix in 1:Nx) {
   	    # get data for specific cell and category
	    # Categories now saved as 1-4, not S-P
		cellD   = catD[catD$cell==paste(iy, ix, sep="_"),];
        thePlot = ggplot(data=cellD, aes(x=t, y=count));
        thePlot = thePlot + geom_line(aes(group=species, color=species)) + theme_bw() +
	                        labs(title=titleStr, x=xStr, y="") +
							scale_y_continuous(limits=yLimV) +
	                        theme(legend.position="none");
	  
	    # plot graph in grid (opening new window for first)
	    plotPos = c(iy, ix);
		# make sure to invert the cell Y value to match MatLab (y=1 at the bottom)
        myPlot(thePlot, plotPos=plotPos, newWindow=newWindow, plotDim=c(Ny,Nx), invertY=TRUE);
		newWindow = FALSE;
	  }
	}
  }
}
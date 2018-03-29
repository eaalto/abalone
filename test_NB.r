library(ggplot2);
source("../common/myUtil.r");
source("../common/myGraphics.r");

testNBParams <- function(r=1.25, c=0.1, g=.01, maxT=100, doPlot=TRUE) {
  nV = c(100);
  pV = c(1);
  
  nK = 10^5;
  for (i in 1:maxT) {
    nT  = nV[i];
	pT  = pV[i];
	nT1 = r*(1-nT/nK)*nT*exp(-c*pT);
	pT1 = g*nT*(1-exp(-c*pT));
	nV  = c(nV, nT1);
	pV  = c(pV, pT1);
  }

  if (doPlot) {
    par(mfrow=c(1,2));
    plot(1:(maxT+1), nV);
    plot(1:(maxT+1), pV);
  }
  return(list(nV=nV, pV=pV));
}

tuneNBParams <- function(r=1.25, c=0.9, g=.1, tuneC=TRUE, adjBy=0.75, maxT=100) {
  testC = c;
  testG = g;
  totalRuns = 0;
  while (TRUE) {
    totalRuns = totalRuns + 1;
    resL = testNBParams(r=r, c=testC, g=testG, doPlot=FALSE, maxT=maxT);
	if (resL$pV[maxT+1]>0) {
	  myP("Non-zero predator for c=", testC, ",g=", testG);
	  myWin();
      par(mfrow=c(1,2));
      plot(1:(maxT+1), resL$nV);
      plot(1:(maxT+1), resL$pV);
	  break;
	} else if (totalRuns>50) {
	  myP("Max runs reached.");
	  break;
	} else {
	  if (tuneC) testC = testC*adjBy
	  else testG = testG*adjBy;
	}
  }
}

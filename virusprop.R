#install.packages("igraph")
require(igraph)
require(compiler)
require(doSNOW)
require(foreach)
require(parallel)

threads <- detectCores()
cl <- makeCluster(threads)
registerDoSNOW(cl)

setwd("C:\\Users\\Tapas\\Desktop\\VirusPropagation")
#setwd("/home/tapas/Desktop/VirusPropagation")

filename="static.network";

#read the file and create the graph
vertex_list <- read.table(filename, header=T, quote="\"")
header <- names(vertex_list)
no_vertices <- as.numeric(substring(header[1],2))
no_edges <- as.numeric(substring(header[2],2))
g <- graph.data.frame(vertex_list, directed=F)

######################################################################
##Now calculate the eigen values of the adjacency matrix of this graph
######################################################################
ADJ=get.adjacency(g);
func <- function(x, extra=NULL) { as.vector(ADJ %*% x) };
eigenValues = arpack(func, options=list(n=vcount(g), nev=3, 
                                        ncv=8, which="LM", maxiter=200))$values;
cat(paste("Max Eigen value: ", abs(eigenValues[1])));

######################################################################
##1a Calculate the effective strength(s)
######################################################################
beta=0.2#beta1
delta=0.7#delta1
s=abs(eigenValues[1])*beta/delta;
cat(paste("Effective strength for beta=0.2 and delta=0.7: ", s));

######################################################################
#1b Now vary beta over a range and calculate s
######################################################################
betarange=seq(0,1,0.005);
srange=c();
for(beta in betarange)
{
  srange = c(srange, (abs(eigenValues[1])*beta/delta));
}
#Now calculate beta for which s=1
beta_thres = 1*delta/(abs(eigenValues[1]));
windows()
plot(betarange,srange, type="l",xlab="Transmission probability", ylab="Effective strength")
abline(h=1,col="red",lty=2)
abline(v=beta_thres,col="red",lty=2)
points(beta_thres,1,pch = 23, col="red",cex = 1);
text(beta_thres,1,paste("Beta=",round(beta_thres,digits=5)),col="red",pos=4);
title("Transmission probability Vs Effective strength (delta=0.7)");

######################################################################
#1c Now vary delta with beta being fixed
######################################################################
deltarange=seq(0,1,0.005);
srange=c();
beta=0.2;
for(delta in deltarange)
{
  srange = c(srange, (abs(eigenValues[1])*beta/delta));
}
######################################################################
#Now calculate delta for which s=1
######################################################################
delta_thres = (abs(eigenValues[1]))*beta;
if(delta_thres > 1)
  delta_thres=1;
windows()
plot(deltarange,srange, type="l",xlab="Healing probability", ylab="Effective strength")
abline(h=1,col="red",lty=2)
abline(v=delta_thres,col="red",lty=2)
points(delta_thres,1,pch = 23, col="red",cex = 1);
text(delta_thres,1,paste("Delta=",round(delta_thres,digits=5)),col="red",pos=2);
title("Healing probability vs Effective strength(beta=0.2)");

#############################################################
#######1d repeat above with different beta and delta defaults
#############################################################
##1a Calculate the effective strength(s)
beta=0.01#beta2
delta=0.6#delta2
s=abs(eigenValues[1])*beta/delta;

cat(paste("Effective strength with beta=0.01 and delta=0.6: ", s));

#############################################################
#1b Now vary beta over a range and calculate s
#############################################################
betarange=seq(0,1,0.005);
srange=c();
for(beta in betarange)
{
  srange = c(srange, (abs(eigenValues[1])*beta/delta));
}
#############################################################
#Now calculate beta for which s=1
#############################################################
beta_thres = 1*delta/(abs(eigenValues[1]));
windows()
plot(betarange,srange, type="l",xlab="Transmission probability", ylab="Effective strength")
abline(h=1,col="red",lty=2)
abline(v=beta_thres,col="red",lty=2)
points(beta_thres,1,pch = 23, col="red",cex = 1);
text(beta_thres,1,paste("Beta=",round(beta_thres,digits=5)),col="red",pos=4);
title("Transmission probability Vs Effective strength (delta=0.6)");

#############################################################
#1c Now vary delta with beta being fixed
#############################################################
deltarange=seq(0,1,0.005);
srange=c();
beta=0.01;
for(delta in deltarange)
{
  srange = c(srange, (abs(eigenValues[1])*beta/delta));
}

######################################################################
#Now calculate delta for which s=1
######################################################################
delta_thres = (abs(eigenValues[1]))*beta;
if(delta_thres > 1)
  delta_thres=1;
windows()
plot(deltarange,srange, type="l",xlab="Healing probability", ylab="Effective strength")
abline(h=1,col="red",lty=2)
abline(v=delta_thres,col="red",lty=2)
points(delta_thres,1,pch = 23, col="red",cex = 1);
text(delta_thres,1,paste("Delta=",round(delta_thres,digits=5)),col="red",pos=2);
title("Healing probability vs Effective strength(beta=0.01)");

######################################################################
#Q2
######################################################################
simulateTentimesOver100Steps <- function(g, no_vertices, beta, delta, timesteps)
{
  #this will hold the total fraction of infected nodes across multiple iterations
  # so we can average it out at the end.
  totalfractioninfnodes=rep(0,length(timesteps));
  iterations=10;
  foreach(i=1:10) %do% #Run it ten times-- made it 3 for faster run..hardcode it to 10.
  {
    #First randomly infect n/10 nodes
    c=round(0.1*no_vertices);
    infectedNodes = sample(1:no_vertices,c);
    fractionInfNodes = c();
    originfectedNodes = c();
    
    foreach(j = 1:100) %do% 
    {
      #nodesHealed=0;
      originfectedNodes = c(infectedNodes);
      for(node in originfectedNodes)
      {
        friends = neighbors(g,node);
        for(friend in friends)
        {
          #Check if each friend will be infected based on beta
          #if friend is not already infected then no need to add - check infectedNodes
          if(runif(1) <= beta & !(friend %in% infectedNodes))#it is going to be infected
          {
            infectedNodes = c(infectedNodes,friend);
          }
        }
        #Now see if the current node which is infected will be healed
        if(runif(1) <= delta)
        {
          infectedNodes = infectedNodes [!infectedNodes == node];
        }
      }
      fractionInfNodes = c(fractionInfNodes, length(infectedNodes)/no_vertices);
    }
    #windows()
    #plot(timesteps,fractionInfNodes, type="l",xlab="Time steps", ylab="Fraction of infected nodes");
    #title(paste("Iter:", k,"Infected nodes over time \n(beta=", beta, "delta=", delta,")"));
    totalfractioninfnodes = totalfractioninfnodes + fractionInfNodes;
    #print("totalfractioninfnodes:")
    #print (totalfractioninfnodes[i])
  }
  #windows()
  plot(timesteps,totalfractioninfnodes/iterations, type="l",xlab="Time steps", ylab="Fraction of infected nodes");
  print(timesteps)
  print(totalfractioninfnodes/iterations)
  title(paste("Average (10 runs) Infected nodes over time \n(beta=", beta, "delta=", delta,")"));
}
#first run this with beta=0.2 and delta=0.7
#simulateTentimesOver100Steps(g,no_vertices,0.2,0.7,1:100);
#Next run this with beta=0.01 and delta=0.6
#simulateTentimesOver100Steps(g,no_vertices,0.01,0.6,1:100);

######################################################################
#Q3
######################################################################
#return strength of graph .inout parameters ipgraph, beta and delta
strength <- function(g,beta,delta)
{
  ADJ=get.adjacency(g);
  func <- function(x, extra=NULL) { as.vector(ADJ %*% x) };
  eigenValues = arpack(func, options=list(n=vcount(g), nev=3, 
                                          ncv=8, which="LM", maxiter=200))$values;
  #cat(paste("Max Eigen value: ", abs(eigenValues[1])));
  s=abs(eigenValues[1])*beta/delta;
  return(s)
}

#Policy A
KRandomNodes <- function(g, no_vertices, k)
{
  immunizedNodes = sample(1:no_vertices, k);
  totalNodes  = c(1:no_vertices)
  nonImmunizedNodes = totalNodes [!totalNodes %in% immunizedNodes];
  subGraph = subgraph(g, nonImmunizedNodes)
  return(subGraph)
}
subGraphKRandom = KRandomNodes(g, no_vertices, 200)
strength_immunized = strength(subGraphKRandom,0.2,0.7)
simulateTentimesOver100Steps(subGraphKRandom,vcount(subGraphKRandom),0.2,0.7,1:100);

#Policy B
KHighestDegreeNodes <- function(g, no_vertices, k)
{
  topKHighestDeg = sort(degree(g), decreasing=TRUE)[1:k]
  topKHighestDegNodeIndex = attributes(topKHighestDeg)[[1]] # high degree nodes got immunized
  subGraph = g
  subGraph = delete.vertices(subGraph, topKHighestDegNodeIndex)
  return(subGraph)
}
subGraphKHighest = KHighestDegreeNodes(g, no_vertices, 200)
strength_immunized = strength(subGraphKHighest,0.2,0.7)
#simulateTentimesOver100Steps(subGraphKHighest,vcount(subGraphKHighest),0.2,0.7,1:100);

#Policy C
# A more efficient solution is to maintain a max-heap and extract k max node, calling heapify
# k times. So the overall complexity will be O(n) to build heap and extracting K 
# largest nodes is O(K log(N)) and to remove the node from graph depends on the
# underlying implemented data structure.
# T(n) = O(n) + O(K log(N)) + tx, where tx is the time to remove node from adjacency list or adjacency matrix.
highestDegreeNodeIteratively <- function(g, no_vertices, k)
{
  subgraph = g;
  for(i in 1:k) {
    subgraph = KHighestDegreeNodes(subgraph, no_vertices, 1) # k = 1
  }
  return(subgraph)
}
subGraphHighest = highestDegreeNodeIteratively(g, no_vertices, 200)
strength_immunized = strength(subGraphHighest,0.2,0.7)
#simulateTentimesOver100Steps(subGraphHighest,vcount(subGraphHighest),0.2,0.7,1:100);

#Policy D
eigenHeuristic <- function(g, no_vertices, k)
{
  eigenVector = arpack(func, options=list(n=vcount(g), nev=3, ncv=8, which="LM", maxiter=200))$vectors;
  vecord = order(abs(eigenVector[,1]), decreasing=TRUE)[1:k]
  subGraph = g
  subGraph = delete.vertices(subGraph, vecord)
  return(subGraph)
}
eigenGraph = eigenHeuristic(g, no_vertices, 200)
strength_immunized = strength(eigenGraph,0.2,0.7)
#simulateTentimesOver100Steps(eigenGraph,vcount(eigenGraph),0.2,0.7,1:100);

#3(e)
#minimum vaccines or policy A
minvaccines = seq(200,5715,25)
strength_immunized = c()
min_reqd = 0
min_strength = 0
for(vaccines in minvaccines)
{
  subGraphKRandom = KRandomNodes(g, no_vertices, vaccines)
  strength_immunized = c(strength_immunized ,strength(subGraphKRandom,0.2,0.7))
}

for(i in 1:length(strength_immunized))
{
  if(strength_immunized[i] < 1)
  {
    min_reqd = minvaccines[i]
    min_strength = strength_immunized[i]
    break
  }
}
plot(minvaccines,strength_immunized, type="l",xlab="Vaccines applied", ylab="Effective strength")
abline(h=min_strength,col="red",lty=2)
abline(v=min_reqd,col="red",lty=2)
points(min_reqd,1,pch = 23, col="red",cex = 1)
text(min_reqd,1,paste("Min vaccines=",round(min_reqd,digits=5)),col="red",pos=2)
title("Vaccines applied vs Effective strength(beta=0.2,delta = 0.7) for policy A")

#minimum vaccines or policy B
minvaccines = seq(200,5715,25)
strength_immunized = c()
min_reqd = 0
min_strength = 0

for(vaccines in minvaccines)
{
  subGraphKHighest = KHighestDegreeNodes(g, no_vertices, vaccines)
  strength_immunized = c(strength_immunized ,strength(subGraphKHighest,0.2,0.7))
}

# for(vaccines in minvaccines)
# {
#   subGraphKHighest = KHighestDegreeNodes(g, no_vertices, vaccines)
#   matrix = get.adjacency(subGraphKHighest)
#   lambda <- eigen(matrix) 
#   strength_immunized = c(strength_immunized ,max(abs(lambda$values)*0.2/0.7))
# }

for(i in 1:length(strength_immunized))
{
  if(strength_immunized[i] < 1)
  {
    min_reqd = minvaccines[i]
    min_strength = strength_immunized[i]
    break
  }
}
windows()
plot(minvaccines,strength_immunized, type="l",xlab="Vaccines applied", ylab="Effective strength")
abline(h=min_strength,col="red",lty=2)
abline(v=min_reqd,col="red",lty=2)
points(min_reqd,1,pch = 23, col="red",cex = 1)
text(min_reqd,1,paste("Min vaccines=",round(min_reqd,digits=3)),col="red",pos=4)
title("Vaccines applied vs Effective strength(beta=0.2,delta = 0.7) policy B")
stopCluster(cl)

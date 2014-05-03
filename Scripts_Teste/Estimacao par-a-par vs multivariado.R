library(pair.mglmm)
library(testpack)

# Simulação de dados
phi<-1; p<-1.6 

mydat<-data.sim(10,'CP',exp,xi=p,phi=phi)
dat<-mydat$Data


cl<-makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl,{library(testpack);
                 library(SparseGrid);
                 library(numDeriv);
                 library(cplm)})

mainForm<-value~-1+variable:period+(-1+variable|ID)

# Estimaçao par-a-par
system.time({
  m0.1<-pair.mglmm:::glmmMultiCP(mainForm,dat$variable,dat,TRUE,cl)
  JKJ.th<-rcov2(m0.1,mainForm)
})  # 3 var's: 165.88 seg; 10 var's: 2409.51 seg

# Criação da matriz A
tes<-lapply(m0.1,pair.mglmm:::format0CP)
tes.1<-do.call(rbind,tes)

df.l<-sort(unique(tes.1$Parametro))

mat<-lapply(df.l,'==',tes.1$Parametro)
mat<-do.call(rbind,mat)

mat<-sweep(mat,1,rowSums(mat),"/") # matriz A(!)

est<-data.frame(Parametro=df.l,value=mat%*%tes.1$value) # ok!
EPmat<-mat%*%JKJ.th%*%t(mat)

final<-cbind(est,sqrt(diag(EPmat)))

stopCluster(cl)


# Modelo multiivariado
system.time(m1<-cpglmm(mainForm,data=dat)) # 3 var's: 21.14 seg; 10 var's: 1498.53 seg
res.m1<-summary(m1)[[5]]

final2<-final[match(rownames(res.m1),final$Parametro),]
cbind(final2,res.m1[,-3])





#              Parametro         value sqrt(diag(EPmat))    Estimate Std. Error
# 3   variableC1:period1  0.3934203764         0.1506481  0.39898120  0.2481833
# 30  variableC2:period1 -0.2524721766         0.2542322 -0.23419428  0.2359502
# 42  variableC3:period1 -0.1949775907         0.2439394 -0.12951036  0.2736643
# 53  variableC4:period1  0.5026645194         0.2547277  0.52527965  0.2133813
# 63  variableC5:period1  0.4096029726         0.2438091  0.40605692  0.2215142
# 72  variableC6:period1  0.8859556603         0.1578931  0.85022690  0.2338953
# 80  variableC7:period1  0.3814635034         0.2966184  0.36850502  0.2037005
# 87  variableC8:period1  0.0258297695         0.1677403  0.01903817  0.1967003
# 93  variableC9:period1  0.1988708075         0.2490985  0.20191165  0.1977878
# 17 variableC10:period1 -1.0952539576         0.2346118 -1.13095152  0.2626847
# 4   variableC1:period2  0.8426982209         0.1525283  0.84945554  0.2402759
# 31  variableC2:period2 -0.2804815360         0.2500481 -0.22605579  0.2357208
# 43  variableC3:period2 -0.8887732114         0.2693877 -0.91685263  0.2934131
# 54  variableC4:period2 -0.7952895545         0.3455878 -0.77457324  0.2497311
# 64  variableC5:period2 -0.9703670978         0.1907310 -0.94995692  0.2598008
# 73  variableC6:period2 -0.7630833672         0.1860849 -0.77362402  0.2715616
# 81  variableC7:period2 -0.6004103291         0.3127710 -0.58951401  0.2322421
# 88  variableC8:period2  0.6600139965         0.1652702  0.67310533  0.1778960
# 94  variableC9:period2  1.4340207031         0.2093930  1.43202097  0.1682939
# 18 variableC10:period2  0.5822073351         0.2497006  0.58002301  0.2118470
# 5   variableC1:period3 -0.3853815538         0.1769062 -0.34336245  0.2644539
# 32  variableC2:period3 -1.2839020321         0.3156861 -1.23060382  0.2690783
# 44  variableC3:period3 -1.8969112678         0.2660640 -1.90910722  0.3273985
# 55  variableC4:period3  0.6564075398         0.2712325  0.66522795  0.2103268
# 65  variableC5:period3 -1.1461059436         0.2698033 -1.12478254  0.2660733
# 74  variableC6:period3 -0.8064597282         0.2161979 -0.83646801  0.2734918
# 82  variableC7:period3  0.0969305341         0.3148719  0.10410523  0.2107427
# 89  variableC8:period3  1.6350896304         0.1105463  1.62865091  0.1563478
# 95  variableC9:period3  0.6703313302         0.2092001  0.66318492  0.1853685
# 19 variableC10:period3 -0.1445005886         0.2385027 -0.16446245  0.2305766
# 6   variableC1:period4  1.3589845797         0.1659890  1.36362727  0.2327650
# 33  variableC2:period4 -0.0002315348         0.2645248  0.04902118  0.2283206
# 45  variableC3:period4 -0.5862909347         0.2397462 -0.65803602  0.2862968
# 56  variableC4:period4 -0.0119967477         0.2664196 -0.01221702  0.2265626
# 66  variableC5:period4  0.5310034630         0.2417449  0.54994108  0.2183866
# 75  variableC6:period4  0.9324982152         0.1975476  0.90671266  0.2329432
# 83  variableC7:period4  1.6296418179         0.2709644  1.63561948  0.1773242
# 90  variableC8:period4  0.7342862660         0.1450227  0.74451559  0.1760560
# 96  variableC9:period4  1.4001032943         0.2375067  1.39268641  0.1690676
# 20 variableC10:period4  0.9127560352         0.2213431  0.91767513  0.2047962
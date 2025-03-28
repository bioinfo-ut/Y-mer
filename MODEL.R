# ------------------------------------------------------
# Haplotyypide maeaeramine 
# osa A
# mudeli loomine teadaolevate haplotyypide pealt.
# ------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
algandmed  = args[1]
mudelifail = args[2]
diagnostics = args[3]
diagnosticsPATH = args[4]

# diagnostics=TRUE; diagnosticsPATH="c:/temp/"

if (is.null(diagnostics) | is.na(diagnostics)) diagnostics=TRUE else if (diagnostics!=TRUE & diagnostics!=FALSE) {
	print(paste("Third argument (diagnostics) should either be TRUE or FALSE, but you have supplied the value:", diagnostics)); diagnostics=TRUE}

if (is.null(diagnosticsPATH) | is.na(diagnosticsPATH)) {PLOTdiagnostics=FALSE; print("Diagnostics plots will not be produced")} else {PLOTdiagnostics=TRUE;
	print(paste("Diagnostics plots will be save on directory:", diagnosticsPATH))}

if(!require(stringr)){
    install.packages("stringr")
    library(stringr)
}

# kuhu salvestada loodud mudel:
# mudelifail="C:/Maido/centro/NIPT/HAPLOmudel7.RData"

# Kust leida algandmete fail (k-meeride arvud + inimeste haplotbid)
# algandmed="C:/Maido/centro/NIPT/HalbJaHea/mudel/uus_nimekiri_NIPTIGA.txt"

# Andmete sisselugemine

# Andmefail 1 struktuur:
# 1. reas indiviidi id
# 2. reas haplotyyp
# 3. - n. reas k_meeri nimi (esimene tt m mudeli k-meeride jaoks ja n nipti k-meeride jaoks) 
#						ja mitu korda antud k-meeri on ntud erinevates inimestes:

# 	HG01890	HG02982	NA18971	NA20911	
#	A0	A0	C1	C1
# mGCATTAGGTATATCTCCCAATGCTA	10	1	9	100
# mAACAGACCTGCAGCTGAGGGTCCTG	8	7	27	103
# ...
# nTTTTTTTCATGGCCCTAGAAAAAAA	9	17	13	1





and=read.table(algandmed, skip=2, header=F, stringsAsFactors = FALSE)
ID=read.table(algandmed, skip=0,nrows=1, header=F, stringsAsFactors = FALSE)
haplo=read.table(algandmed, skip=1,nrows=1, header=F, stringsAsFactors = FALSE)

ID2=unlist(ID)
haplo2=unlist(haplo)
# table(haplo2)

kmeer_type=substr(and[,1],1,1)
# table(kmeer_type)

kmeer=substr(and[,1],2,nchar(and[,1]))
and2=and[,-1]


th=names(table(haplo2))
n_haplo =length(th)
n_kmeere=length(and2[,1])
n_kmeere_mudel=sum(kmeer_type=="m")
n_kmeere_nipt=sum(kmeer_type=="n")
n_inimesi=length(and2[1,])
n_per_haplo=table(haplo2)
k_kmeer=nchar(and[1,1])-1


katvus0=colMeans(and2[kmeer_type=="n",])
and_korduseid = matrix(NA, nrow=nrow(and2), ncol=ncol(and2))

NIPT_A=str_count(kmeer,"A")
NIPT_T=str_count(kmeer,"T")
abiX=NIPT_A+NIPT_T

sum_AT=as.vector(by(abiX[kmeer_type=="n"], abiX[kmeer_type=="n"], mean))
sum_n=as.vector(by(rep(1, length(abiX[kmeer_type=="n"])), abiX[kmeer_type=="n"], sum))

n1=dnorm(sum_AT, mean=0, sd=5)
n2=dnorm(sum_AT, mean=5, sd=5)
n3=dnorm(sum_AT, mean=10, sd=5)
n4=dnorm(sum_AT, mean=15, sd=5)
n5=dnorm(sum_AT, mean=20, sd=5)
n6=dnorm(sum_AT, mean=26, sd=5)

xxx=0:k_kmeer
newdata=data.frame(n1=dnorm(xxx, mean=0, sd=5), n2=dnorm(xxx, mean=5, sd=5), n3=dnorm(xxx, mean=10, sd=5),
 	n4=dnorm(xxx, mean=15, sd=5), n5=dnorm(xxx, mean=20, sd=5), n6=dnorm(xxx, mean=26, sd=5), sum_n=1)

ATkorrektsioon=matrix(NA, n_inimesi, k_kmeer+1)
Dli = rep(NA, n_inimesi)
Eci2= rep(NA, n_inimesi)
Dci = rep(NA, n_inimesi)
E1_ci= rep(NA, n_inimesi)

print("Applying CG-related corrections")

for (i_inimene in 1:n_inimesi){
# i_inimene=1
sum_count=as.vector(by(and2[kmeer_type=="n",i_inimene], abiX[kmeer_type=="n"], sum))


	m=try ( suppressWarnings({ATkorrigeerimismudel=glm(sum_count~n1+n2+n3+n4+n5+n6+offset(log(sum_n)), family=poisson())}), silent=TRUE)
      if(!inherits(m, "try-error")) {
		ATkorrektsioon[i_inimene,] =predict(ATkorrigeerimismudel, newdata=newdata, type="resp")
      } else {
		ATkorrektsioon[i_inimene,] = rep(katvus0[i_inimene], nrow(newdata))
		print("AT-correction ignored once during model building - cannot fit the model")
	}

# ATkorrigeerimismudel=glm(sum_count~n1+n2+n3+n4+n5+n6+offset(log(sum_n)), family=poisson())
# ATkorrektsioon[i_inimene,]=predict(ATkorrigeerimismudel, newdata=newdata, type="resp")
# Dli=(var(NIPT[,i_inimene+1])-coverage - Dci)/(Eci2+Dci)

# D(l_i) = (D(Y)-E(Y)- D(c_i)*E(l_i)^2)/(E(c_i^2)+D(c_i))
Eci2[i_inimene]= sum((ATkorrektsioon[i_inimene,][sum_AT+1])**2 *sum_count)/sum(sum_count)
abi=(ATkorrektsioon[i_inimene,][sum_AT+1])
E1_ci[i_inimene]= sum(1/abi[abi!=0] *sum_count[abi!=0])/sum(sum_count[abi!=0])
Dci[i_inimene]=Eci2[i_inimene]-katvus0[i_inimene]^2

# (var(and2[kmeer_type=="n", i_inimene])-katvus0[i_inimene] - Dci[i_inimene])/(Eci2[i_inimene])
# (var(and2[kmeer_type=="n", i_inimene])-katvus0[i_inimene] - Dci[i_inimene])/(Eci2[i_inimene]+Dci[i_inimene])

# Dli[i_inimene]=(var(and2[kmeer_type=="n", i_inimene])-katvus0[i_inimene] - Dci[i_inimene])/(Eci2[i_inimene]+Dci[i_inimene])
Dli[i_inimene]=(var(and2[kmeer_type=="n", i_inimene])-katvus0[i_inimene] - Dci[i_inimene])/(Eci2[i_inimene])

if (Dli[i_inimene]<0) Dli[i_inimene]=0


if (PLOTdiagnostics){ 
  failinimi=paste(diagnosticsPATH, "/diag_", ID2[i_inimene], ".png", sep="")
  png(failinimi, width=800, height=600)

  # D(Y|c) = c_i*E(l_i) + c_i^2*D(l_i)
  DY_c=ATkorrektsioon[i_inimene,][sum_AT+1]+(ATkorrektsioon[i_inimene,][sum_AT+1])^2*Dli[i_inimene]
  s_viga=sqrt(DY_c/sum_n)
  plot(0:k_kmeer, ATkorrektsioon[i_inimene,], type="l", lwd=3, xlab="number of A/T letters in kmer", ylab="Sequencing coverage" , ylim=range(c(ATkorrektsioon[i_inimene,],sum_count/sum_n)), main=ID2[i_inimene])
  abline(h=katvus0[i_inimene], lwd=2, col="blue")
  points(sum_AT, sum_count/sum_n, pch=20, col="red", cex=1.5)
  arrows(sum_AT, sum_count/sum_n-1.96*s_viga, sum_AT, sum_count/sum_n+1.96*s_viga, code=3, length=0.1, angle=90)

  legend("topright", c("coverage: A/T corrected", "coverage: uncorrected"), title="Coverage estimate", lwd=c(3,2), col=c("black", "blue"))
  dev.off()
}

# system(paste('"C:/Program Files/Google/Chrome/Application/chrome.exe" ', failinimi, sep=""))

and_korduseid[, i_inimene]=and2[,i_inimene]/ATkorrektsioon[i_inimene,][abiX]
}


if (PLOTdiagnostics){ 
  failinimi=paste(diagnosticsPATH, "/diagnostics1.png", sep="")
  png(failinimi, width=1024, height=800)

  plot(Dli, Dci, pch=".", xlab="Uncorrectable variation (Dli)", ylab="Correctable variation (CG-bias, Dci)")
  text(Dli,Dci, ID2, xpd=NA)
  dev.off()
}
# system(paste('"C:/Program Files/Google/Chrome/Application/chrome.exe" ', failinimi, sep=""))



kauguspiir=100000

# Leiame haplotyypide kaupa korduste arvu keskmised ja dispersioonid
print("Calculating Haplotype means and variances for k-mers")


keskmMAT=matrix(NA, n_kmeere, n_haplo)
dispMAT0 =matrix(NA, n_kmeere, n_haplo)
dispMAT =matrix(NA, n_kmeere, n_haplo)
keskmine_katvus_per_haplo = rep(NA, n_haplo)
dispMAT =matrix(NA, n_kmeere, n_haplo)

for (i in 1:n_haplo){
  keskmMAT[,i]=rowMeans(and_korduseid[,haplo2==th[i]])
  dispMAT0[,i] =(rowMeans(and_korduseid[,haplo2==th[i]]**2)-keskmMAT[,i]^2)*n_per_haplo[i]/(n_per_haplo[i]-1)
  keskmine_katvus_per_haplo[i]=mean(katvus0[haplo2==th[i]])
  # D(gamma) = (   D(Y/c) - mu* E(E(l|i)/c_i) - mu^2 * E(D(l|i))   )   /  (E(D(l|i)) + 1 )
  dispMAT[,i] = dispMAT0[,i] - keskmMAT[,i]*1*mean(E1_ci[haplo2==th[i]]) - keskmMAT[,i]**2*mean(Dli[haplo2==th[i]])
  dispMAT[dispMAT[,i]<0,i]=0
}


print("Mean coverages:")
names(keskmine_katvus_per_haplo)=th
print(keskmine_katvus_per_haplo)


print("Finding distances from haplogroup means")

# Toorkaugused
s_dist=matrix(NA, n_inimesi, n_haplo)
# Ell=matrix(NA, n_inimesi, n_haplo)

# Probability of seeing k-mer that actually do not exist in genome
p_error =0.01/3*(1-0.01/3)**(k_kmeer-1)


# Tskel e inimeste:
for (i_inimene in 1:n_inimesi){
   oodatavKatvus_1x=(ATkorrektsioon[i_inimene,][abiX[kmeer_type=="m"]+1])

   # Tskel e kgi haplotpide:
   for (i in 1:n_haplo){

	keskvek=keskmMAT[kmeer_type=="m",i]
	dispvek=dispMAT[kmeer_type=="m",i]
	DY_c_long = (keskvek*oodatavKatvus_1x ) + (dispvek+keskvek^2)*(oodatavKatvus_1x)^2*Dli[i_inimene]

#    	sigma2_k= dispvek*(oodatavKatvus_1x^2) + DY_c_long +(keskvek^2/keskmine_katvus_per_haplo[i] +dispvek)*(1/n_per_haplo[i]) *oodatavKatvus_1x^2   + 0.01/3*(1-0.01/3)*oodatavKatvus_1x  
   	sigma2_k= dispvek*(oodatavKatvus_1x^2) + DY_c_long +(keskvek^2/keskmine_katvus_per_haplo[i] +dispvek)*(1/n_per_haplo[i]) *oodatavKatvus_1x^2   + (p_error**2*oodatavKatvus_1x + p_error*(1-p_error)*oodatavKatvus_1x) 

   	dnvek0=dnorm(and2[kmeer_type=="m",i_inimene], mean=keskvek*oodatavKatvus_1x, sd=sqrt(sigma2_k), log=T)

   	# Yletaeituvuste kaitse 
   	dnvek0[dnvek0<(-kauguspiir)] = (-kauguspiir)
   	dnvek0[dnvek0>(kauguspiir)]  = (kauguspiir)

	# Kaugus - oodatav kaugus
      s_dist[i_inimene,i]=mean(-dnvek0) + (mean(-0.5*log(2*pi*sigma2_k))-1/2)

#	ajut=mean(-0.5*log(2*pi*sigma2_k))-1/2

#	Ell[i_inimene,i]=ajut

    }  # Tskel e haplogruppideyoutube


} # Tskel e inimeste

# tulemus=(s_dist+Ell)
# tulemus[1:5,]
# s_dist[1:5,]



print("Typical representatives")

mean_distance=matrix(NA, nrow=n_haplo ,ncol=n_haplo)
var_distance=list(n_haplo)
sd_distance=rep(NA, n_haplo)


for (i in 1:n_haplo){ 
  mean_distance[i, ]= colMeans(s_dist[haplo2==th[i],])
  var_distance[[i]]= var(s_dist[haplo2==th[i],])
  sd_distance[i]=sqrt(var_distance[[i]][i,i])
}

meanvector=diag(mean_distance)
# sd_distance

model_version=7

if (diagnostics){
 save(model_version, n_per_haplo, keskmine_katvus_per_haplo, n_haplo, n_kmeere, n_kmeere_mudel, n_kmeere_nipt, n_inimesi, kmeer_type,
	keskmMAT, dispMAT, th, kauguspiir, diagnostics, abiX, sum_AT, sum_n, p_error, k_kmeer,
	mean_distance, sd_distance, Dli, haplo2, ID2, ATkorrektsioon, katvus0, kmeer, var_distance, file=mudelifail)
} else {
 save(model_version, n_per_haplo, keskmine_katvus_per_haplo, n_haplo, n_kmeere, n_kmeere_mudel, n_kmeere_nipt, n_inimesi, kmeer_type,
	keskmMAT, dispMAT, th, kauguspiir, diagnostics, abiX, sum_AT, sum_n, p_error, k_kmeer, file=mudelifail)
}


# ------------------------------------------------------
# Haplotyypide maeaeramine 
# osa B
# haplotyybi maaramine kasutades varasemalt hinnatud mudelit
# ------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
mudelifail = args[1]
algandmed  = args[2]
vastusfail = args[3]
diagnostics = args[4]
diagnosticsPATH = args[5]

# diagnostics=TRUE; diagnosticsPATH="c:/temp/"

if (is.null(diagnostics) | is.na(diagnostics)) diagnostics=TRUE else if (diagnostics!=TRUE & diagnostics!=FALSE) {
	print(paste("Third argument (diagnostics) should either be TRUE or FALSE, but you have supplied the value:", diagnostics)); diagnostics=TRUE}

if (is.null(diagnosticsPATH) | is.na(diagnosticsPATH)) {PLOTdiagnostics=FALSE; print("Diagnostics plots will not be produced")} else {PLOTdiagnostics=TRUE;
	print(paste("Diagnostics plots will be save on directory:", diagnosticsPATH))}

# algandmed="C:/Maido/centro/NIPT/KatsedVer5/V12939_counts.txt"
# algandmed="C:/Maido/centro/NIPT/KatsedVer5/V03482.txt"
# algandmed="C:/Maido/centro/NIPT/HalbJaHea/sim3.txt"
# algandmed="C:/Maido/centro/NIPT/HalbJaHea/mudel/HG00101_counts.txt"
# algandmed="C:/Maido/centro/NIPT/HalbJaHea/sim.txt"
# algandmed="C:/Maido/centro/NIPT/HalbJaHea/sim2.txt"
# mudelifail="C:/Maido/centro/NIPT/HAPLOmudel7.RData"


# mudelifail="C:/Maido/centro/NIPT/HAPLOmudel7.RData"
# mudelifail="C:/Maido/centro/NIPT/KatsedVer5/mudel7_fem.Rdata"

load(mudelifail)
if (model_version!=7) stop("Requires version 7 model!")

# Andmete sisselugemine
and0=read.table(algandmed, header=F, nrows=1, stringsAsFactors = FALSE)

if (prod(sapply(and0, is.character))==1) { 
	print("Sample ID-s on a first row")
	and00=read.table(algandmed, header=F, skip=1, nrows=1, stringsAsFactors = FALSE)
	if (prod(sapply(and00, is.character))==1) { 
		print("Ignoring 2. row") 
		and=read.table(algandmed, skip=2, header=FALSE) 
	    } else {
			and=read.table(algandmed, skip=1, header=FALSE) 
      }
   } else {
	print("No sample ID-s found")
	and=read.table(algandmed, header=F)
}

if (is.character(and[1,1])) {
   print("Ignoring 1. column - assuming it contains k-mers")
   and=and[,-1, drop=FALSE]
}

if (prod(sapply(and0, is.character))==1) {  colnames(and)=unlist(and0) }


if (dim(and)[1]!=n_kmeere) stop("model k-mers do not match datafile k-mers!")

katvus0_proov=colMeans(and[kmeer_type=="n",, drop=FALSE])
names(katvus0_proov)=colnames(and)

print("Sample coverages:")
print(katvus0_proov)

if (sum(katvus0_proov==0)>0) {stop("Y-chromosome coverage estimated to be exactly 0 for some samples!!!")}

n_proove=dim(and)[2]


n1=dnorm(sum_AT, mean=0, sd=5)
n2=dnorm(sum_AT, mean=5, sd=5)
n3=dnorm(sum_AT, mean=10, sd=5)
n4=dnorm(sum_AT, mean=15, sd=5)
n5=dnorm(sum_AT, mean=20, sd=5)
n6=dnorm(sum_AT, mean=26, sd=5)

xxx=0:k_kmeer
newdata=data.frame(n1=dnorm(xxx, mean=0, sd=5), n2=dnorm(xxx, mean=5, sd=5), n3=dnorm(xxx, mean=10, sd=5),
 	n4=dnorm(xxx, mean=15, sd=5), n5=dnorm(xxx, mean=20, sd=5), n6=dnorm(xxx, mean=26, sd=5), sum_n=1)

ATkorrektsioon_proov=matrix(NA, n_proove, k_kmeer+1)
Dli = rep(NA, n_proove)
Eci2= rep(NA, n_proove)
Dci = rep(NA, n_proove)
E1_ci= rep(NA, n_proove)

print("Applying CG-related corrections")

for (i_inimene in 1:n_proove){
# i_inimene=1
sum_count=as.vector(by(and[kmeer_type=="n",i_inimene], abiX[kmeer_type=="n"], sum))

	m=try ( suppressWarnings({ATkorrigeerimismudel=glm(sum_count~n1+n2+n3+n4+n5+n6+offset(log(sum_n)), family=poisson())}), silent=TRUE)
      if(!inherits(m, "try-error")) {
		ATkorrektsioon_proov[i_inimene,] =predict(ATkorrigeerimismudel, newdata=newdata, type="resp")
      } else {
		ATkorrektsioon_proov[i_inimene,] = rep(katvus0_proov[i_inimene], nrow(newdata))
		print("AT-correction ignored once during model building - cannot fit the model")
	}

# ATkorrigeerimismudel=glm(sum_count~n1+n2+n3+n4+n5+n6+offset(log(sum_n)), family=poisson())
# ATkorrektsioon[i_inimene,]=predict(ATkorrigeerimismudel, newdata=newdata, type="resp")

# D(l_i) = (D(Y)-E(Y)- D(c_i)*E(l_i)^2)/(E(c_i^2)+D(c_i))
Eci2[i_inimene]= sum((ATkorrektsioon_proov[i_inimene,][sum_AT+1])**2 *sum_count)/sum(sum_count)
abi=(ATkorrektsioon_proov[i_inimene,][sum_AT+1])
E1_ci[i_inimene]= sum(1/abi[abi!=0] *sum_count[abi!=0])/sum(sum_count[abi!=0])
Dci[i_inimene]=Eci2[i_inimene]-katvus0_proov[i_inimene]^2
Dli[i_inimene]=(var(and[kmeer_type=="n", i_inimene])-katvus0_proov[i_inimene] - Dci[i_inimene])/(Eci2[i_inimene]+Dci[i_inimene])
if (Dli[i_inimene]<0) Dli[i_inimene]=0

if (PLOTdiagnostics){ 
  failinimi=paste(diagnosticsPATH, "/diag_", colnames(and)[i_inimene], ".png", sep="")
  png(failinimi, width=800, height=600)

  # D(Y|c) = c_i*E(l_i) + c_i^2*D(l_i)
  DY_c=ATkorrektsioon_proov[i_inimene,][sum_AT+1]+(ATkorrektsioon_proov[i_inimene,][sum_AT+1])^2*Dli[i_inimene]
  s_viga=sqrt(DY_c/sum_n)
  plot(0:k_kmeer, ATkorrektsioon_proov[i_inimene,], type="l", lwd=3, xlab="number of A/T letters in kmer", ylab="Sequencing coverage" , ylim=range(c(ATkorrektsioon_proov[i_inimene,],sum_count/sum_n)), main=colnames(and)[i_inimene] )
  abline(h=katvus0_proov[i_inimene], lwd=2, col="blue")
  points(sum_AT, sum_count/sum_n, pch=20, col="red", cex=1.5)
  suppressWarnings(arrows(sum_AT, sum_count/sum_n-1.96*s_viga, sum_AT, sum_count/sum_n+1.96*s_viga, code=3, length=0.1, angle=90))

  legend("topright", c("coverage: A/T corrected", "coverage: uncorrected"), title="Coverage estimate", lwd=c(3,2), col=c("black", "blue"))
  dev.off()
} # diagnostics
} # tsükkel üle inimeste

print("Coverage(s):")
print(katvus0_proov)

print("Sequencing coverage uniformity, Dli (smaller is better):")
print(Dli)


# Toorkaugused
s_dist_proov=matrix(NA, n_proove, n_haplo)


# Tsükkel üle inimeste:
for (i_inimene in 1:n_proove){
   oodatavKatvus_1x=(ATkorrektsioon_proov[i_inimene,][abiX[kmeer_type=="m"]+1])

   # Tsükkel üle kõigi haplotüüpide:
   for (i in 1:n_haplo){

	keskvek=keskmMAT[kmeer_type=="m",i]
	dispvek=dispMAT[kmeer_type=="m",i]

	DY_c_long = (keskvek*oodatavKatvus_1x ) + (dispvek+keskvek^2)*(oodatavKatvus_1x)^2*Dli[i_inimene]

   	sigma2_k= dispvek*(oodatavKatvus_1x^2) + DY_c_long +(keskvek^2/keskmine_katvus_per_haplo[i] +dispvek)*(1/n_per_haplo[i]) *oodatavKatvus_1x^2   + (p_error**2*oodatavKatvus_1x + p_error*(1-p_error)*oodatavKatvus_1x) 

   	dnvek0=dnorm(and[kmeer_type=="m",i_inimene], mean=keskvek*oodatavKatvus_1x, sd=sqrt(sigma2_k), log=T)
# 	dnvek0 = (and[kmeer_type=="m",i_inimene]/oodatavKatvus_1x-keskvek)**2

   	# Yletaeituvuste kaitse 

   	dnvek0[dnvek0<(-kauguspiir)] = (-kauguspiir)
   	dnvek0[dnvek0>(kauguspiir)]  = (kauguspiir)

#      s_dist_proov[i_inimene,i]=mean(dnvek0)

#      s_dist_proov[i_inimene,i]=mean(-dnvek0)
#      s_dist_proov[i_inimene,i]=mean(-dnvek0)+mean(-0.5*log(2*pi*sigma2_k))-1/2
      s_dist_proov[i_inimene,i]=mean(-dnvek0)+mean(-0.5*log(2*pi*sigma2_k))-1/2


#	ajut=mean(-0.5*log(2*pi*sigma2_k))-1/2
#	Ell_proov[i_inimene,i]=ajut
    }  # Tsükkel üle haplogruppideyoutube

} # Tsükkel üle inimeste

print("Raw distances:")
rownames(s_dist_proov)=colnames(and)
colnames(s_dist_proov)=th
print(s_dist_proov)


# vahed2=t((apply(s_dist_proov, 1, function(x){x-min(x)})))
# vahed2

parimHaplo_vek=apply(s_dist_proov, 1, function(x){(1:length(x))[x==min(x)]})

Evahed=matrix(NA, n_proove, n_haplo)
Dvahed=matrix(NA, n_proove, n_haplo)

# Tsükkel üle inimeste:
for (i_inimene in 1:n_proove){
   parimHaplo=parimHaplo_vek[i_inimene]
   oodatavKatvus_1x=(ATkorrektsioon_proov[i_inimene,][abiX[kmeer_type=="m"]+1])

   # Tsükkel üle kõigi haplotüüpide:
   for (i in 1:n_haplo){

	keskvek=keskmMAT[kmeer_type=="m",i]
	dispvek=dispMAT[kmeer_type=="m",i]

	keskvek_parim=keskmMAT[kmeer_type=="m", parimHaplo]
	dispvek_parim=dispMAT[kmeer_type=="m", parimHaplo]

	DY_c_long = (keskvek*oodatavKatvus_1x ) + (dispvek+keskvek^2)*(oodatavKatvus_1x)^2*Dli[i_inimene]
   	sigma2_k= dispvek*(oodatavKatvus_1x^2) + DY_c_long +(keskvek^2/keskmine_katvus_per_haplo[i] +dispvek)*(1/n_per_haplo[i]) *oodatavKatvus_1x^2   + (p_error**2*oodatavKatvus_1x + p_error*(1-p_error)*oodatavKatvus_1x) 
   	dnvek0=dnorm(and[kmeer_type=="m",i_inimene], mean=keskvek*oodatavKatvus_1x, sd=sqrt(sigma2_k), log=T)


	DY_c_long_parim = (keskvek_parim*oodatavKatvus_1x ) + (dispvek_parim+keskvek_parim^2)*(oodatavKatvus_1x)^2*Dli[i_inimene]
   	sigma2_k_parim = dispvek_parim*(oodatavKatvus_1x^2) + DY_c_long_parim +(keskvek_parim^2/keskmine_katvus_per_haplo[i] +dispvek_parim)*(1/n_per_haplo[i]) *oodatavKatvus_1x^2   + (p_error**2*oodatavKatvus_1x + p_error*(1-p_error)*oodatavKatvus_1x) 
   	dnvek0_parim =dnorm(and[kmeer_type=="m",i_inimene], mean=keskvek_parim*oodatavKatvus_1x, sd=sqrt(sigma2_k_parim), log=T)

   	# Yletaeituvuste kaitse 
   	dnvek0[dnvek0<(-kauguspiir)] = (-kauguspiir)
   	dnvek0[dnvek0>(kauguspiir)]  = (kauguspiir)

   	dnvek0_parim[dnvek0_parim<(-kauguspiir)] = (-kauguspiir)
   	dnvek0_parim[dnvek0_parim>(kauguspiir)]  = (kauguspiir)

#      s_dist_proov[i_inimene,i]=mean(-dnvek0) + mean(-0.5*log(2*pi*sigma2_k))-1/2
      ajut = -dnvek0 -0.5*log(2*pi*sigma2_k)   + dnvek0_parim + 0.5*log(2*pi*sigma2_k_parim)  
      Evahed[i_inimene,i]=mean(ajut)
      Dvahed[i_inimene,i]=var(ajut)


#	ajut=mean(-0.5*log(2*pi*sigma2_k))-1/2
#	Ell_proov[i_inimene,i]=ajut
    }  # Tsükkel üle haplogruppideyoutube

} # Tsükkel üle inimeste


abimat= pnorm(-Evahed, 0, sqrt(Dvahed/n_kmeere_mudel))
pvalue=rep(NA, n_proove)
alternatives=list(n_proove)

for (i in 1:n_proove){
  abi=abimat[i,-parimHaplo_vek[i]] 
  ind=order(abi[abi>0.05], decreasing=TRUE)
  alternatives[[i]]=(th[-parimHaplo_vek[i]][abi>0.05][ind])
  pvalue[i]=max(abi)
}

print("Most likely haplogroup:")
results=data.frame(sample=colnames(and), coverage=katvus0_proov, haplogroup=th[parimHaplo_vek], pvalue=pvalue, alternatives=unlist(lapply(alternatives, FUN=paste, collapse=",")) )
print(results)

#save(results, file="results.RData")
save(results, file=vastusfail)
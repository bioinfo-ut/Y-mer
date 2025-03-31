# Nõuded:
# Väljakutsumisel esimne argument on sisendfail, teine väljundfail
# Sisendfaili formaat:
# 1. reas inimeste id-d (ei kasutata)
# 2. reas inimeste sood (F/M, tähistus oluline!), inimeste sood algavad 8. veerust!
# edasi kmeeride korduste arvud: 1. tulbas kmeer, korduste arvud (sagedus/katvus) alates 8. tulbast

args = commandArgs(trailingOnly=TRUE)

# sisend = file("C:/Maido/centro/pooliktabel.txt")
sisend = file(args[1])
open(sisend, open = "r")

# valjund = file("C:/Maido/centro/out2.txt")
valjund=file(args[2])
open(valjund, open = "w")


abi1=readLines(sisend, n=1) # Loe sisse inimese ID
inimene=strsplit(abi1, "\t")[[1]][-(1:7)]

abi2=readLines(sisend, n=1) # Loe sisse sugu näitav rida (F/M)
sugu=strsplit(abi2, "\t")[[1]][-(1:7)]

# Vaatluste arv, meeste arv ja naiste arv:
N=length(sugu)
m=sum(sugu=="M")
n=sum(sugu=="F")

konst1=m*(m+1)/2
konst2=+0.5-0.5*n*m
konst3=sqrt(m*n*(N+1)/12)


abi3=readLines(sisend, n=1) 

while  (length(abi3)>0){
  abi3a = strsplit(abi3, "\t")[[1]]
  kmeer = abi3a[1]
  # Vajalikud numbrid alates 8. tulbast (viska esimesed 7 minema)
  numbrid = as.numeric(abi3a[-(1:7)])

  mitupuudu=sum(is.na(numbrid))
  numbrid[is.na(numbrid)]=2
  numbrid=c(numbrid, rep(2, N-length(numbrid)))

  stat=sum(rank(numbrid)[sugu=="M"])-konst1
  pvaartus= 2*pnorm(-abs(stat+konst2)/konst3)
  AUC=stat/(m*n)

  valjundrida = paste(kmeer, stat, AUC, pvaartus)
  writeLines(valjundrida, valjund)  

  abi3=readLines(sisend, n=1) 
}

close(valjund)
close(sisend)

unlink(valjund)
unlink(sisend)

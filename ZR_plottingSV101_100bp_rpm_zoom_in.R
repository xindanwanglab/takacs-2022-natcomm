

###################################  100bp bin #############################
# after removing 3 missing replicons, Borrelia B31 combined genome=1,453,623 bp

input_bin = 1

strain1 = read.csv ("AntiGFP_CJW_Bb101_rep1_modified_clean.csv", header=T)
#strain2 = read.csv ("AntiGFP_CJW_Bb101_rep1_sv101_clean.csv", header=T)
#strain3 = read.csv ("AntiGFP_CJW_Bb101_rep1_19plusplasmid_clean.csv", header=T)

norm_strain=function(a,b,c) {
  d=(a[2]-a[1]+1)%%b
  if (d==0){
    order=c(c$Value[a[1]:a[2]])
    order_norm = order/((sum(c$Value)/42)  / 1e6)
    order_bin = colSums(matrix(order_norm,nrow=b))/b
    return(order_bin)
  }
  else{
    order=c(c$Value[a[1]:(a[2]-d)])
    order_last = c(c$Value[(a[2]-d+1):a[2]]) #last bin is less than b bp, need to normalize differently
    order_last_norm = order_last/((sum(c$Value)/42)  / 1e6)
    order_norm = order/((sum(c$Value)/42)  / 1e6)
    order_bin = colSums(matrix(order_norm,nrow=b))/b
    order_last_bin = sum(order_last_norm)/d
    bin_whole = c(order_bin, order_last_bin)
    return(bin_whole)
  }
}

#bin size
size=100
#combine 5 copies of SV101
#onecopy=c((length(strain3$Value)-7111+1),length(strain3$Value))
#onlysc101=c(1,7111)
cp1=c((length(strain1$Value)-7111+1),length(strain1$Value))
cp2=cp1-7111
cp3=cp2-7111
cp4=cp3-7111
cp5=cp4-7111
n_cp1=norm_strain(cp1,size,strain1)
n_cp2=norm_strain(cp2,size,strain1)
n_cp3=norm_strain(cp3,size,strain1)
n_cp4=norm_strain(cp4,size,strain1)
n_cp5=norm_strain(cp5,size,strain1)

plasmid_combined=n_cp1+n_cp2+n_cp3+n_cp4+n_cp5
#only=norm_strain(onlysc101, size, strain2)
#old=norm_strain(onecopy, size, strain3)
#reorder_plasmid=c(plasmid_combined[21:72],plasmid_combined[1:20])
#reorder_old=c(old[21:72],old[1:20])
#reorder_only=c(only[21:72],only[1:20])



plot(seq(0, 71, 1), plasmid_combined, type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", ylim=c(0,30000), col="black")

#lines(seq(0, 71, 1), reorder_only, type="l", lwd=1.5, col="red")
#lines(seq(0, 71, 1), reorder_old, type="l", lwd=1.5, col="grey")

lines(c(24.98,32.08), c(15000,15000), type="l", lwd=2, col="black")
lines(c(32.39,39.55), c(16000,16000), type="l", lwd=2, col="green")
lines(c(45.46,50.79), c(17000,17000), type="l", lwd=2, col="black")

axis(1, at=c(0,10,20,30,40,50,60,70), labels=c(0,1,2,3,4,5,6,7))
axis(2)


#plotting beginning and end of BB0432gene
#abline(v=c(450174,450887)/100, lty=2, col="red")

box()

dev.copy2pdf(file="3_antiGFPonSV101_Bb101_rep1_modifiedGenome.pdf", width = 15, height = 5)




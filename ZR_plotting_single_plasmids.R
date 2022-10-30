

#importing data
strain1 = read.csv ("antiGFP_Borrelia_CJW_Bb378_rep2_clean.csv", header=T)
strain2 = read.csv ("antimCherry_CJW_Bb379_rep2_clean.csv", header=T)
strain3 = read.csv ("antiGFP_CJW_Bb488-rep2_clean.csv", header=T)
strain4 = read.csv ("antiGFP_CJW_Bb473_rep2_clean.csv", header=T)

#plasmid boundaries
lp28_3= c(910725, 939325)
lp25=c(939326, 963502)
lp28_2=c(963503, 993268)
lp38=c(993269, 1032097)
lp36=c(1032098, 1068946)
lp28_4=c(1068947, 1096269)
lp54=c(1096270, 1149926)
cp26=c(1149927, 1176424)
lp17=c(1176425, 1193245)
lp28_1=c(1193246, 1221400)
cp32_1=c(1221401, 1252150)
cp32_3=c(1252151, 1282373)
cp32_4=c(1282374, 1312672)
cp32_6=c(1312673, 1342510)
cp32_7=c(1342511, 1373310)
cp32_8=c(1373311, 1404195)
cp32_9=c(1404196, 1434846)
lp21=c(1434847, 1453623)


#function to normalize data using RPM
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

#function to plot
plt_strain=function(e,f) {
  s=(e[2]-e[1])%/%size
  plot(seq(0, s, 1), norm_strain(a=e, b=size, c=strain1), type="l", lwd=1.5, axes=FALSE, ylim=c(0, 350), ylab="RPM", xlab="genome position (kb)", col="dark green")

  lines(seq(0, s, 1), norm_strain(a=e, b=size, c=strain2), type="l", lwd=1.5, col="red")

  lines(seq(0, s, 1), norm_strain(a=e, b=size, c=strain3), type="l", lwd=1.5, col="green")

  lines(seq(0, s, 1), norm_strain(a=e, b=size, c=strain4), type="l", lwd=1.5, col="black")

  axis(1)
  axis(2)
  box()
  w=(s*size)%/%5000
  dev.copy2pdf(file=paste(f,"_Bb378_Bb379_Bb488_Bb473_rep2.pdf",sep=""), width = w, height = 5)
  dev.off()
}
#set bin size
size=100

#plotting
plt_strain(e=lp28_3, f="lp28_3")
plt_strain(e=lp25, f="lp25")
plt_strain(e=lp28_2, f="lp28_2")
plt_strain(e=lp38, f="lp38")
plt_strain(e=lp36, f="lp36")
plt_strain(e=lp28_4, f="lp28_4")
plt_strain(e=lp54, f="lp54")
plt_strain(e=cp26, f="cp26")
plt_strain(e=lp17, f="lp17")
plt_strain(e=lp28_1, f="lp28_1")
plt_strain(e=cp32_1, f="cp32_1")
plt_strain(e=cp32_3, f="cp32_3")
plt_strain(e=cp32_4, f="cp32_4")
plt_strain(e=cp32_6, f="cp32_6")
plt_strain(e=cp32_7, f="cp32_7")
plt_strain(e=cp32_8, f="cp32_8")
plt_strain(e=cp32_9, f="cp32_9")
plt_strain(e=lp21, f="lp21")
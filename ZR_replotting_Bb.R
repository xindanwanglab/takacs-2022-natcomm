

#importing data
Bb101 = read.csv ("AntiGFP_CJW_Bb101_rep1_clean.csv", header=T)
Bb378_rep1 = read.csv ("antiGFP_CJW_Bb378_rep1_clean.csv", header=T)
Bb378_rep2 = read.csv ("antiGFP_CJW_Bb378_rep2_clean.csv", header=T)
Bb379_rep1 = read.csv ("antimCherry_CJW_Bb379_rep1_clean.csv", header=T)
Bb379_rep2 = read.csv ("antimCherry_CJW_Bb379_rep2_clean.csv", header=T)
Bb403_GFP = read.csv ("antiGFP_CJW_Bb403_rep1_clean.csv", header=T)
Bb403_mCherry = read.csv ("antimCherry_CJW_Bb403_rep2_clean.csv", header=T)
Bb473_rep1 = read.csv ("antiGFP_CJW_Bb473_rep1_clean.csv", header=T)
Bb473_rep2 = read.csv ("antiGFP_CJW_Bb473_rep2_clean.csv", header=T)
Bb488_rep1 = read.csv ("antiGFP_CJW_Bb488_rep1_clean.csv", header=T)
Bb488_rep2 = read.csv ("antiGFP_CJW_Bb488_rep2_clean.csv", header=T)
Bb519 = read.csv ("antiGFP_CJW_Bb519_rep1_clean.csv", header=T)
Bb520 = read.csv ("antiGFP_CJW_Bb520_rep1_clean.csv", header=T)
Bb524 = read.csv ("antiGFP_CJW_Bb524_rep1_clean.csv", header=T)
Bb525 = read.csv ("antimCherry_CJW_Bb525_rep1_clean.csv", header=T)
Bb601 = read.csv ("antimCherry_CJW_Bb601_rep1_clean.csv", header=T)
Bb610 = read.csv ("AntiGFP_CJW_Bb610_rep1_clean.csv", header=T)
Bb403_nopld_GFP = read.csv ("antiGFP_CJW_Bb403_rep1_nopld_clean.csv", header=T)
Bb403_nopld_mCherry = read.csv ("antimCherry_CJW_Bb403_rep2_nopld_clean.csv", header=T)
Bb403_parBdel_GFP = read.csv ("antiGFP_CJW_Bb403_rep1_parBdel_clean.csv", header=T)
Bb403_parBdel_mCherry = read.csv ("antimCherry_CJW_Bb403_rep2_parBdel_clean.csv", header=T)




#function to normalize data using RPM
norm_strain=function(a) {
  d=(length(a$Value))%%size
  if (d==0){
    order=c(a$Value)
    order_norm = order/((sum(a$Value)/42)  / 1e6)
    order_bin = colSums(matrix(order_norm,nrow=size))/size
    return(order_bin)
  }
  else{
    order=c(a$Value[1:(length(a$Value)-d)])
    order_last = c(a$Value[(length(a$Value)-d+1):length(a$Value)]) #last bin is less than b bp, need to normalize differently
    order_last_norm = order_last/((sum(a$Value)/42)  / 1e6)
    order_norm = order/((sum(a$Value)/42)  / 1e6)
    order_bin = colSums(matrix(order_norm,nrow=size))/size
    order_last_bin = sum(order_last_norm)/d
    bin_whole = c(order_bin, order_last_bin)
    return(bin_whole)
  }
}

#function to plot
plt_whole_strain=function(a,c) {
  #s=(e[2]-e[1])%/%size
  e=norm_strain(a)
  plot(seq(0, (length(e)-1), 1), e, ylim=c(0, 4000), type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", col="black")

  #lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="red")

  #lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="green")

  #lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="black")

  axis(1, at=c(0, 200, 400, 600, 800, 910, 1100, 1300, 1453), labels=c(0, 200, 400, 600, 800, 910, 1100, 1300, 1453), cex.axis=0.6)
  axis(2)
  box()
  #w=(s*size)%/%5000
  dev.copy2pdf(file=paste(c,"_RPM.pdf",sep=""), width = 5, height = 5)
  dev.off()
}
#set bin size
size=1000

#plotting
plt_whole_strain(a=Bb473_rep1,c="AntiGFP_Bb473_rep1_whole_genome")
plt_whole_strain(a=Bb378_rep1,c="AntiGFP_Bb378_rep1_whole_genome")
plt_whole_strain(a=Bb379_rep1,c="AntiGFP_Bb379_rep1_whole_genome")
plt_whole_strain(a=Bb488_rep1,c="AntiGFP_Bb488_rep1_whole_genome")
plt_whole_strain(a=Bb473_rep2,c="AntiGFP_Bb473_rep2_whole_genome")
plt_whole_strain(a=Bb378_rep2,c="AntiGFP_Bb378_rep2_whole_genome")
plt_whole_strain(a=Bb379_rep2,c="AntiGFP_Bb379_rep2_whole_genome")
plt_whole_strain(a=Bb488_rep2,c="AntiGFP_Bb488_rep2_whole_genome")







#funtion to plotting zoom in
size=100
plt_part_strain=function(a,c,y) {
  #s=(e[2]-e[1])%/%size
  e=norm_strain(a)
  plot(seq(0, (length(e)-1), 1), e, xlim=c(4400, 4600), ylim=c(0,y), type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", col="black")
  
  #lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="red")
  
  #lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="green")
  
  #lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="black")
  
  axis(1, at=c(4400, 4450, 4500,4550, 4600), labels=c(440,445,450,455,460))
  axis(2)
  box()
  #w=(s*size)%/%5000
  dev.copy2pdf(file=paste(c,"_RPM.pdf",sep=""), width = 5, height = 5)
  dev.off()
}

#plotting zoom-in plots
plt_part_strain(a=Bb378_rep1,c="AntiGFP_Bb378_rep1_zoom_in", y=8000)
plt_part_strain(a=Bb379_rep1,c="AntimCherry_Bb379_rep1_zoom_in", y=8000)
plt_part_strain(a=Bb524,c="AntiGFP_Bb524_rep1_zoom_in", y=8000)
plt_part_strain(a=Bb525,c="AntimCherry_Bb525_rep1_zoom_in", y=8000)
plt_part_strain(a=Bb488_rep1,c="AntiGFP_Bb488_rep1_zoom_in", y=1000)
plt_part_strain(a=Bb520,c="AntiGFP_Bb520_rep1_zoom_in", y=1000)
plt_part_strain(a=Bb519,c="AntiGFP_Bb519_rep1_zoom_in", y=1000)
plt_part_strain(a=Bb610,c="AntiGFP_Bb610_rep1_zoom_in", y=1000)
plt_part_strain(a=Bb601,c="AntimCherry_Bb601_rep1_zoom_in", y=8000)




#funtion to plotting zoom in with bars
size=100
plt_part_strain_bars=function(a,c,y,g,h,i,j,k,l,m,n,o,p,q,r) {
  #s=(e[2]-e[1])%/%size
  e=norm_strain(a)
  plot(seq(0, (length(e)-1), 1), e, xlim=c(4400, 4600), ylim=c(0,y), type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", col="black")
  
  lines(c(g,h), c(6800,6800), type="l", lwd=2, col="black")
  
  lines(c(i,j), c(6900,6900), type="l", lwd=2, col="black")
  
  lines(c(k,l), c(7000,7000), type="l", lwd=2, col="black")
  
  lines(c(m,n), c(7100,7100), type="l", lwd=2, col="black")
  
  lines(c(o,p), c(7200,7200), type="l", lwd=2, col="black")
  
  lines(c(q,r), c(7300,7300), type="l", lwd=2, col="black")
  
  axis(1, at=c(4400, 4450, 4500,4550, 4600), labels=c(440,445,450,455,460))
  axis(2)
  box()
  #w=(s*size)%/%5000
  dev.copy2pdf(file=paste(c,"_RPM.pdf",sep=""), width = 5, height = 5)
  dev.off()
}

#plt_part_strain_bars(a=Bb378_rep1,c="AntiGFP_Bb378_rep1_zoom_in_bars", y=8000,4483.17,4489.19,4490.28,4497.80,4497.83,4512.37,4513.18,4521.27,4525.11,4532.93,4534.12,4558.44)

#plt_part_strain_bars(a=Bb379_rep1,c="AntimCherry_Bb379_rep1_zoom_in_bars", y=8000,4483.17,4489.19,4490.28,4497.80,4497.83,4504.96,4505.74,4513.83,4517.67,4532.81,4534.00,4558.32)

plt_part_strain_bars(a=Bb524,c="Fig_5F_AntiGFP_Bb524_rep1_zoom_in_bars", y=8000,4483.17,4489.19,4490.28,4497.80,4497.83,4512.37,4513.18,4521.27,4525.28,4549.60,0,0)

plt_part_strain_bars(a=Bb525,c="Fig_5G_AntimCherry_Bb379_rep1_zoom_in_bars", y=8000,4483.17,4489.19,4489.80,4497.89,4501.73,4516.87,4518.06,4542.38,0,0,0,0)

plt_part_strain_bars(a=Bb601,c="Fig_6C_AntimCherry_Bb601_rep1_zoom_in_bars", y=8000,4505.02,4511.04,4512.13,4519.65,4519.68,4526.81,4529.33,4537.15,4538.34,4562.66,0,0)




plt_part_strain_bars=function(a,c,y,g,h,i,j,k,l,m,n,o,p,q,r) {
  #s=(e[2]-e[1])%/%size
  e=norm_strain(a)
  plot(seq(0, (length(e)-1), 1), e, xlim=c(4400, 4600), ylim=c(0,y), type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", col="black")
  
  lines(c(g,h), c(680,680), type="l", lwd=2, col="black")
  
  lines(c(i,j), c(690,690), type="l", lwd=2, col="black")
  
  lines(c(k,l), c(700,700), type="l", lwd=2, col="black")
  
  lines(c(m,n), c(710,710), type="l", lwd=2, col="black")
  
  lines(c(o,p), c(720,720), type="l", lwd=2, col="black")
  
  lines(c(q,r), c(730,730), type="l", lwd=2, col="black")
  
  axis(1, at=c(4400, 4450, 4500,4550, 4600), labels=c(440,445,450,455,460))
  axis(2)
  box()
  #w=(s*size)%/%5000
  dev.copy2pdf(file=paste(c,"_RPM.pdf",sep=""), width = 5, height = 5)
  dev.off()
}


plt_part_strain_bars(a=Bb488_rep1,c="Fig_5H_AntiGFP_Bb488_rep1_zoom_in_bars", y=1000,4483.17,4489.19,4490.28,4505.27,4505.30,4512.43,4513.21,4521.30,4525.14,4532.96,4534.15,4558.47)
plt_part_strain_bars(a=Bb520,c="Fig_5H_AntiGFP_Bb520_rep1_zoom_in_bars", y=1000,4483.17,4489.19,4490.28,4505.27,4505.30,4512.43,4513.21,4521.30,4525.31,4549.63,0,0)
plt_part_strain_bars(a=Bb519,c="Fig_5H_AntiGFP_Bb519_rep1_zoom_in_bars", y=1000,4483.17,4489.19,4490.28,4505.27,4506.07,4514.16,4518.00,4525.82,4527.01,4551.33,0,0)
plt_part_strain_bars(a=Bb610,c="Fig_5H_AntiGFP_Bb610_rep1_zoom_in_bars", y=1000,4483.17,4489.19,4490.28,4505.27,4505.30,4511.86,4512.64,4520.73,4524.57,4532.39,4533.58,4557.90)




#plotting strain with two lines
plt_part_strain_two=function(a,b,c,y,g,h,i,j,k,l,m,n,o,p,q,r) {
  #s=(e[2]-e[1])%/%size
  e=norm_strain(a)
  plot(seq(0, (length(e)-1), 1), e, xlim=c(4400, 4600), ylim=c(0,y), type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", col="blue")
  
  lines(seq(0, (length(e)-1), 1), norm_strain(b), type="l", lwd=1.5, col="red")
  
  lines(c(g,h), c(6800,6800), type="l", lwd=2, col="black")
  
  lines(c(i,j), c(6900,6900), type="l", lwd=2, col="black")
  
  lines(c(k,l), c(7000,7000), type="l", lwd=2, col="black")
  
  lines(c(m,n), c(7100,7100), type="l", lwd=2, col="black")
  
  lines(c(o,p), c(7200,7200), type="l", lwd=2, col="black")
  
  lines(c(q,r), c(7300,7300), type="l", lwd=2, col="black")
  axis(1, at=c(4400, 4450, 4500,4550, 4600), labels=c(440,445,450,455,460))
  axis(2)
  box()
  #w=(s*size)%/%5000
  dev.copy2pdf(file=paste(c,"_RPM.pdf",sep=""), width = 5, height = 5)
  dev.off()
}

plt_part_strain_two(a=Bb403_GFP, b=Bb403_mCherry, c="AntiGFP_Bb403Blue_mCherry_Bb403red_zoom_in", y=8000,4483.17,4489.19,4490.28,4497.80,4497.83,4512.37,4513.18,4521.27,4525.11,4532.93,4534.12,4558.44)
plt_part_strain_two(a=Bb403_nopld_GFP, b=Bb403_nopld_mCherry, c="AntiGFP_Bb403Blue_mCherry_Bb403red_zoom_in_no_plasmids", y=8000,4483.17,4489.19,4490.28,4497.80,4497.83,4512.37,4513.18,4521.27,4525.11,4532.93,4534.12,4558.44)
plt_part_strain_two(a=Bb403_parBdel_GFP, b=Bb403_parBdel_mCherry, c="AntiGFP_Bb403Blue_mCherry_Bb403red_zoom_in_plasmidparB_deleted", y=8000,4483.17,4489.19,4490.28,4497.80,4497.83,4512.37,4513.18,4521.27,4525.11,4532.93,4534.12,4558.44)

##plotting new Bb101

a=Bb101_unmodi = read.csv ("AntiGFP_CJW_Bb101_rep1_clean.csv", header=T)
b=Bb101_modi = read.csv ("AntiGFP_CJW_Bb101_rep1_modified_clean.csv", header=T)

plt_part_strain_bars(a=Bb101_modi,c="2_AntiGFP_Bb101_rep1_modifiedGenome_zoom_in_bars", y=8000,4483.04,4489.06,4490.15,4497.67,4497.70,4504.83,4507.35,4515.17,4516.36,4540.68,0,0)
plt_part_strain_bars(a=Bb101_unmodi,c="2_AntiGFP_Bb101_rep1_unmodifiedGenome_zoom_in_bars", y=8000,4487.08,4493.10,4494.19,4501.71,4501.74,4508.87,4511.39,4519.21,4520.40,4544.72,0,0)
plt_whole_strain(a=Bb101_unmodi,c="1_AntiGFP_Bb101_rep1_unmodifiedGenome_whole_genome")

cp1=c((length(a$Value)-7111+1),length(a$Value))
cp2=cp1-7111
cp3=cp2-7111
cp4=cp3-7111
cp5=cp4-7111

size=1000
pld=c(a$Value[cp1[1]:cp1[2]]+a$Value[cp2[1]:cp2[2]]+a$Value[cp3[1]:cp3[2]]+a$Value[cp4[1]:cp4[2]]+a$Value[cp5[1]:cp5[2]])

order=c(a$Value[1:(length(a$Value)-5*7111)],pld[1:(length(pld)-83)])
order_last = c(pld[(length(pld)-83+1):length(pld)]) #last bin is less than b bp, need to normalize differently
order_last_norm = order_last/((sum(a$Value)/42)  / 1e6)
order_norm = order/((sum(a$Value)/42)  / 1e6)
order_bin = colSums(matrix(order_norm,nrow=size))/size
order_last_bin = sum(order_last_norm)/83
bin_whole = c(order_bin, order_last_bin)

plot(seq(0, 1106, 1), bin_whole, ylim=c(0, 12000), type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", col="black")

#lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="red")

#lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="green")

#lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="black")

axis(1, at=c(0, 200, 400, 600, 800, 910, 1100, 1300), labels=c(0, 200, 400, 600, 800, 910, 1100, 1300), cex.axis=0.6)
axis(2)
box()
#w=(s*size)%/%5000
dev.copy2pdf(file="1_AntiGFP_Bb101_rep1_unmidifiedGenome_whole.pdf", width = 5, height = 5)
dev.off()




order=c(b$Value[1:(length(b$Value)-5*7111)],pld[1:(length(pld)-467)])
order_last = c(pld[(length(pld)-467+1):length(pld)]) #last bin is less than b bp, need to normalize differently
order_last_norm = order_last/((sum(a$Value)/42)  / 1e6)
order_norm = order/((sum(b$Value)/42)  / 1e6)
order_bin = colSums(matrix(order_norm,nrow=size))/size
order_last_bin = sum(order_last_norm)/467
bin_whole = c(order_bin, order_last_bin)

plot(seq(0, 1105, 1), bin_whole, ylim=c(0, 12000), type="l", lwd=1.5, axes=FALSE, ylab="RPM", xlab="genome position (kb)", col="black")

#lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="red")

#lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="green")

#lines(seq(0, s, 1), norm_strain(a=e, b=size, c), type="l", lwd=1.5, col="black")

axis(1, at=c(0, 200, 400, 600, 800, 910, 1100, 1300), labels=c(0, 200, 400, 600, 800, 910, 1100, 1300), cex.axis=0.6)
axis(2)
box()
#w=(s*size)%/%5000
dev.copy2pdf(file="1_AntiGFP_Bb101_rep1_midifiedGenome_whole.pdf", width = 5, height = 5)
dev.off()

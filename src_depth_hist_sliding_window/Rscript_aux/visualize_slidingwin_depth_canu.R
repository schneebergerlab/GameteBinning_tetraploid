library(scales)
calculate_genome_size <- function(win, haphigh, diphigh, triphigh, tetraphigh, hap_cov)
{
  # win: is from read from result of CNV_HQ_v3, 
  #      which calculates depth in sliding windows from output of samtools depth.
  genome_size       <- 0 # raw assembly size
  genome_size_infer <- 0 # infered genome size
  size_details      <- matrix(0, nrow=5, ncol=2) # hap~tetraplotig size
  for (wi in 1:length(win$V1))
  {
    wsize   <- win[wi, 3] - win[wi, 2] + 1
    if(win[wi, 5]<=haphigh)
    {
      # haplotig contribute one copy
      genome_size       <- genome_size + wsize
      genome_size_infer <- genome_size_infer + wsize
      size_details[1, 1]<- size_details[1, 1] + wsize
      size_details[1, 2]<- size_details[1, 2] + 1
    }else
      if(win[wi, 5]>=haphigh+1 & win[wi, 5]<=diphigh)
      {
        # diplotig contribute two copies
        genome_size       <- genome_size + wsize
        genome_size_infer <- genome_size_infer + wsize*2
        size_details[2, 1]<- size_details[2, 1] + wsize*2
        size_details[2, 2]<- size_details[2, 2] + 1
      }else
        if(win[wi, 5]>=diphigh+1 & win[wi, 5]<=triphigh)
        {
          # triplotig contribute three copies
          genome_size       <- genome_size + wsize
          genome_size_infer <- genome_size_infer + wsize*3
          size_details[3, 1]<- size_details[3, 1] + wsize*3
          size_details[3, 2]<- size_details[3, 2] + 1
        }else
          if(win[wi, 5]>=triphigh+1 & win[wi, 5]<=tetraphigh)
          {
            # tetraplotig
            genome_size       <- genome_size + wsize
            genome_size_infer <- genome_size_infer + wsize*4
            size_details[4, 1]<- size_details[4, 1] + wsize*4
            size_details[4, 2]<- size_details[4, 2] + 1
          }else
            if(win[wi, 5]>=tetraphigh+1)
            {
              # replotig
              genome_size       <- genome_size + wsize
              genome_size_infer <- genome_size_infer  + wsize*(round(win[wi, 5]/hap_cov))
              size_details[5, 1]<- size_details[5, 1] + wsize*(round(win[wi, 5]/hap_cov))
              size_details[5, 2]<- size_details[5, 2] + 1
            }  
  }
  # 
  return(list(genome_size, genome_size_infer, size_details))
}
#
path <- "/netscratch/dep_mercier/grp_schneeberger/projects/Potato_single_cells/reference_hequan_assembled_canu/s0_read_alignment_rate_checking/otava_smsc/"

pdf(paste(path, "/multiple_slidingwin_depth_comp_canu.pdf", sep=""), height=12, width=14.26)
par(mai = c(1.0, 1.0, 0.5, 0.5)); # margin: bottom, left, top, right
par(mfrow=c(3,2))

# do not use scientific annotation
options(scipen = 999)
# collect win30k
win30k <- 0
#
for (winsize in c(seq(10000, 50000, 10000), 100000))
#for (winsize in seq(10000, 30000, 10000) )
#for (winsize in 30000)
{
  cat("win size :", winsize, "\n")
  winstep <- winsize
  win  <- read.table(paste(path, "/cnv_winsize", winsize, "_step", winstep,"_hq_reformat.txt", sep=""))
  if(winsize == 30000) win30k <- win
  # set up break number
  breaknum <- round(max(win$V5, win$V6))
  # get stats without plotting
  hwin_raw <- hist(win$V6, breaks=breaknum, plot=F) # raw depth
  hwin_nor <- hist(win$V5, breaks=breaknum, plot=F) # normalized with genome-wide average 103x: raw depth / 103 *100
  ylim_both<- max(hwin_nor$counts, hwin_raw$counts)*1.3
  if(winsize==10000) ylim_both = 10000
  # plot
  hwin_raw <- hist(win$V6, 
                   breaks = breaknum, 
                   col    = alpha("yellowgreen", 0.8),              
                   xlim   = c(1, 600), 
                   ylim   = c(0, ylim_both),
                   border = NA, 
                   plot   = T,
                   xlab   = "Depth",
                   ylab   = "Window count in assembly",
                   main   = paste("Depth at windows of ", winsize, " bp (step: ", winstep, " bp)", sep=""),
                   cex.main=1
  )
  hwin_nor <- hist(win$V5, 
                   breaks = breaknum, 
                   col    = alpha("purple", 0.7), 
                   xlim   = c(1, 600), 
                   border = NA, 
                   plot   = T, 
                   add    = T)
  # legend
  legend("topright",
         pch    = c(15,15),
         col    = c("yellowgreen", "purple"),
         legend = c("raw depth", "raw depth/g.avg * 100"),
         horiz  = FALSE,
         border = "NA",
         bty    = "n",
         cex    = 0.9)  
  #
  if(winsize==10000)
  {
    ## size checking according to the above distributions: hap:[0, 160], dip:[161, 250], trip[251, 350], tetrap[351, 450], rep[451, end]
    ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4118015/: A haplotype is a subset of all alleles on specific chromosomes in the population. A diplotype is a subset of all genotypes on homologous chromosome pairs in the population. A specific diplotype is one variant of all possible combinations of the haplotypes that exist in the population.
    # We define haplotigs, diplotigs, triplotigs and tetraplotigs merged from haplotypes, diplotypes, triplotypes, tetraplotypes
    #
    # Real assembly size: 2.437906533 Gb
    ctgsizes    <- read.table("/netscratch/dep_mercier/grp_schneeberger/projects/Potato_single_cells/reference_hequan_assembled_canu/otava_canu_asm/otava_canu.contigs.fasta.ctgsizes_v2")
    observed_assembly_size <- sum(ctgsizes$V2)
    # expected real genome size
    hap_cov     <- 100
    haphigh     <- 160
    diphigh     <- 250
    triphigh    <- 350
    tetraphigh  <- 450
    list_g_size <- calculate_genome_size(win, haphigh, diphigh, triphigh, tetraphigh, hap_cov)
    ## utg009144l	27651 and utg014901l	20522 missing from samtools depth output -- why? list_g_size[[1]] + +27651+20522 = observed_assembly_size
    text(x = 200, y=ylim_both-1200, labels = paste("Assembly size: ", observed_assembly_size, sep=""), pos=4)
    text(x = 200, y=ylim_both-1800, labels = paste("Inferred size: ", list_g_size[[2]], sep=""), pos=4)
    for (lb in 1:5)
    {
      tigtype   <- "haplotig"
      if(lb == 2)
      {
        tigtype <- "diplotig"
      }
      if(lb == 3)
      {
        tigtype <- "triplotig"
      } 
      if(lb == 4)
      {
        tigtype <- "tetraplotig"
      }     
      if(lb == 5)
      {
        tigtype <- "replotig"
      }      
      text(x = 200, y=ylim_both-lb*800-1800, labels = paste(tigtype, ": ", list_g_size[[3]][lb, 1], " bp, w-count: ", list_g_size[[3]][lb, 2], sep=""), pos=4)
    }
  }
}

##
dev.off()






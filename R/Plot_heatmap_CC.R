Plot_heatmap_CC <- function(tmp, schemaTmp, ColSideColors, mainLab=schemaTmp, legend = NULL, fill = NULL) {
  
  if(!is.null(schemaTmp)){
    i <- which(schemaList %in% schemaTmp)
    genesTmp <- subtypeGeneList[[i]]
    idxGene <- match(tmp$featInfo,genesTmp$geneSymbol)
    geneInfo <- genesTmp[idxGene,]
    colors = t(as.matrix(geneInfo$Color))
  }else{
    genesTmp <- tmp$genesTemp
    idxGene <- match(tmp$featInfo,genesTmp$geneSymbol)
    geneInfo <- genesTmp[idxGene,]
    colors = t(as.matrix(geneInfo$Color))
  }

  
  if(ncol(as.matrix(ColSideColors))<=10) {
    ColSideColorsSize <- ncol(as.matrix(ColSideColors))
  } else {
    ColSideColorsSize <- 10
  }
  
  heatmap.3(main= mainLab
            ,x = as.matrix(tmp$ex)
            ,Rowv = tmp$Rowv
            ,Colv = tmp$Colv
            ,dendrogram = "none" ,trace="none" ,scale="row"
            ,labRow = tmp$featInfo
            ,labCol = tmp$sampID
            ,col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99)
            ,ColSideColors = as.matrix(ColSideColors)
            ,ColSideColorsSize = ColSideColorsSize
            ,RowSideColors = colors
            ,RowSideColorsSize=1
            ,lwid = c(1,5)
            ,lhei = c(1,3)
            ,margins =c(8,10)
            ,keysize=0
            ,cex.lab=0.7
            ,cexRow = 1
            ,cexCol = 0.2
  )
  
  # if(is.null(legend)){
  #   legend = c("DECODER_weight","High","Low",
  #            "PurIST_score", "High","Low",
  #            "PurIST_graded","Strong Classical","Likely Classical","Lean Classical","Lean Basal-like","Likely Basal-like","Strong Basal-like",
  #            "PurIST","Basal-like","Classical",
  #            "Moffitt","Basal-like","Classical",
  #            "Collisson","Classical","Exocrine","QM",
  #            "Bailey","ADEX","Immunogenic","Pancreatic Progenitor","Squamous" ,
  #            "Chan-Seng-Yue","Basal-like A","Basal-like B","Hybrid","Classical A","Classical B",
  #            "Puleo", "Pure Classical","Immune Classical","Desmoplastic","Stroma Activated","Pure Basal-like",
  #            "MS","Absent","Activated","Normal",
  #            "MSI","Activated","Immune","Normal",
  #            "Elyada","apCAF","iCAF","myCAF",
  #            "Maurer","ECM-rich","Immune-rich")
  # }
  
  # if(is.null(fill)){
  #   fill = c("white","black","snow2",
  #            "white","black","snow2",
  #            "white","#0000FF","#6666FF","#CCCCFF","#FFECCB","#FFC965","#FFA500",
  #            "white","orange","blue",
  #            "white","orange","blue",
  #            "white","steelblue1","forestgreen","magenta1" ,
  #            "white","darkorchid1","red","slategray2","gold1",
  #            "white","brown","orange","grey","blue3","skyblue",
  #            "white","dodgerblue4","green4","mediumpurple","darkorange2","red3",
  #            "white","black","brown","skyblue",
  #            "white","brown","red","skyblue",
  #            "white","darkslateblue","coral","darkgreen",
  #            "white","purple3","forestgreen")
  # }
  if(!is.null(legend)){
    legend(xy.coords(x=.8,y=.98),
           legend= legend,
           fill= fill,
           border=FALSE, bty="n",
           y.intersp = 0.85 , cex = 0.8)
  }

}

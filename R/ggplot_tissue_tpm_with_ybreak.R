library(ggplot2)

library(argparser)

#============================== 设置参数 ==============================

argv <- arg_parser('')
argv <- add_argument(argv,"--inputfile", help="该选项填写gene_expression_tumor_normal.py的结果文件,ENST00000606017_TPM_for_plot.tsv")
argv <- add_argument(argv,"--out_dir", help="图片的保存路径")
argv <- parse_args(argv)


inputfile <- argv$inputfile
out_dir <- argv$out_dir
########################测试数据=================================
rm(list=ls())
inputfile <- 'ERBB2_ENSG00000141736_tissue_tpm.tsv'
out_dir <- './'
########################
dt<-read.table(inputfile,header = TRUE,sep="\t")

dt[c(1:4),c(1:4)]
#sample	tissue	TPM	group
#TCGA-19-1787-01	GBM	77.56755385612963	TCGA_tumor
#TCGA-S9-A7J2-01	LGG	13.000444366174035	TCGA_tumor
#TCGA-G3-A3CH-11	LIHC_normal	19.190558913055046	TCGA_normal

f_name <- gsub(".tsv","", basename(inputfile))
y_title  <- gsub("_tissue_tpm","", f_name)

#对GTEX组织自定义重要性顺序 有brain
dt$tissue <- factor(dt$tissue,levels = c("Heart","Lung","Pituitary","Nerve","Thyroid","Adrenal Gland","Salivary Gland","Liver","Kidney","Pancreas","Spleen","Blood","Stomach","Brain","Skin","Colon","Small Intestine","Esophagus","Blood Vessel","Bladder","Breast","Muscle","Adipose Tissue","Cervix Uteri","Fallopian Tube","Testis","Prostate","Vagina","Uterus","Ovary","ACC","BLCA_normal","BLCA","BRCA_normal","BRCA","TNBC","CESC_normal","CESC","CHOL_normal","CHOL","COAD_normal","COAD","DLBC","ESCA_normal","ESCA","GBM_normal","GBM","HNSC_normal","HNSC","KICH_normal","KICH","KIRC_normal","KIRC","KIRP_normal","KIRP","LAML","LGG","LIHC_normal","LIHC","LUAD_normal","LUAD","LUSC_normal","LUSC","MESO","OV","PAAD_normal","PAAD","PCPG_normal","PCPG","PRAD_normal","PRAD","READ_normal","READ","SARC_normal","SARC","SKCM_normal","SKCM","STAD_normal","STAD","TGCT","THCA_normal","THCA","THYM_normal","THYM","UCEC_normal","UCEC","UCS","UVM"))


max_value <- max(dt$TPM)
if (max_value >100) {
  library(ggbreak)
  
  p <- ggplot(dt, aes(tissue, TPM, color = factor(group))) +
    geom_boxplot(outlier.shape = NA) +  # 隐藏箱线图异常点
    geom_point(
      aes(y = TPM),
      position = position_jitterdodge(  # 智能抖动避免重叠
        jitter.width = 0.3,
        dodge.width = 0.8
      ),
      shape = 16,
      size = 0.8,  # 适当增大点尺寸
      alpha = 0.5
    ) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
    #annotate("text", x = -Inf, y = 100, 
             #label = "100", hjust = -0.5, vjust = -0.5,
    #         color = "grey40", size = 4) +
    scale_color_manual(values = c("GTEX" = "#3f60aa",
                                  "TCGA_normal" = "#9ec417",
                                  "TCGA_tumor" = "#cc340c")) +
    
    # 坐标轴双断点设置
    scale_y_continuous(
      limits = c(0, 12000),
      breaks = c(seq(0, 200, 20), seq(1000, 12000, 2000)),  # 精细设置刻度
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_break(
      breaks = c(200, 1000),  # 在100-1000之间断开
      scales = c(0.3, 1),     # 放大比例：0-100占30%高度，1000-12000占70%
      space = 0.5,            # 断点间距
      expand = FALSE
    ) +
    
    theme_bw() +
    labs(y = paste0(y_title, " mRNA Expression (TPM)")) +
    theme(
      # 隐藏右侧y轴（新增部分）
      axis.line.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      legend.position = "none"
    )
}else{
  p <- ggplot(dt, aes(tissue, TPM, color = factor(group))) +
    #geom_point(shape = 16)+ geom_jitter(width = 0.3)+
    geom_boxplot()+
    geom_point(aes(y = TPM),position = position_jitter(width = 0.3),shape = 16,size = 0.2,alpha = 0.7)+
    #scale_y_continuous(limits = c(0, 2800),breaks = c(seq(0,2800,by=400)))+
    scale_colour_manual(values=c("GTEX"="#3f60aa","TCGA_normal"="#9ec417","TCGA_tumor"="#cc340c"))+
    theme_bw()+  #去除背景色
    labs(y = paste0(y_title," mRNA Expression (TPM)"))+
    #scale_y_continuous(limits = c(0, 2800),breaks = c(seq(0,2800,by=400))) +
    theme(axis.text.y = element_text(size=15,face="bold"),axis.text.x = element_text(angle=90, hjust=1, vjust=.5,size=16,face="bold"),legend.position='none',axis.title.y = element_text(size=15,face="bold"),axis.title.x = element_blank()) +
    #theme(panel.border = element_blank(),axis.line = element_line(colour = "black",size = 0.5))+ # 去除边框
    theme(panel.grid =element_blank()) #去除网格线s
  
}
 p
#ggsave(file = paste0(out_dir,'/',f_name,'1.png'),plot = p, width = 49, height = 14,units = 'in')

ggsave(file = paste0(out_dir,'/',f_name,'.png'),plot = p, width = 80, height = 30,units = 'cm')

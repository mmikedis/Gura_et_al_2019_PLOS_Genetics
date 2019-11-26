### make_venn_diagram_with_lines.R
### 11/29/17
### run as:
### make_venn_diagram_with_lines.R REGION
### where REGION is the region in which the peaks are to be analyzed


 
args <- commandArgs()
print(args)

region <- args[6]  
library(VennDiagram)

file_name = paste(region, "_overlapping_peaks" ,sep = "")
numbers = read.table(file_name, sep ="\t", header=TRUE)


colors<-c("#b3e2cd","#fdcdac","#cbd5e8")
output_file_name = paste(region, "_peaks_venn.pdf" ,sep = "")
pdf(output_file_name)
venn.plot  <-draw.triple.venn(numbers$No_peaks[1], numbers$No_peaks[2],numbers$No_peaks[3],numbers$No_peaks[4],numbers$No_peaks[5],numbers$No_peaks[6],numbers$No_peaks[7], category=c("Rep1","Rep2","Rep3"),
        fill = c("white", "white", "white"),,
        col = c("black", "black", "black"),
        #lty = "blank",
        cex = 3,
        cat.cex = 2,
        fontfamily = rep("sans", 7));
grid.draw(venn.plot);
dev.off()
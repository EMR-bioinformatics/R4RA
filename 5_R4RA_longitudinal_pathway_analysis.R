# Make sure that Cytoscape v3.6+' is started and the Cytoscape Apps 
# 'yFiles Layout Algorithms' and 'ClueGO v2.5.2' are installed 
# before running this script

rm(list = ls())
require(xml2)
require(RJSONIO)
require(httr) 
library(writexl)
library(ggplot2)
library(plyr)
library(svglite)
library(ggsoccer)
# NOTE: httr needs the libssl-dev and libcurl4-openssl-dev packages in linux 
# (ubuntu: sudo apt-get install libssl-dev libcurl4-openssl-dev)

# # prepare data
# ymatrix <-  readRDS("matrix_from_glmmSeq_longitudinal_analysis.RDS")
# genes <- rownames(ymatrix[which(ymatrix$col == "FDR_time < 0.05" &
#                                 ymatrix$x < -0.5  &
#                                 ymatrix$y < -0.5  ),])
filename <- "DEGs.FDR_time.0.05.neg.logFC.less.than.-0.5"
# write.csv(genes, paste0("Cluego Analysis/", filename, ".txt"), 
#           quote = F, row.names = F)

# Helper function to transform output into a data frame
text.to.data.frame <- function(table.text) {
  table <- NULL
  rows <- unlist(strsplit(result.table.text, split="\n"))
  header <- t(unlist(strsplit(rows[1], split="\t")))
  for(i in 2:length(rows)) {
    if(is.null(table)) {
      table <- t(unlist(strsplit(rows[i], split="\t")))
    } else {
      table <- rbind(table,t(unlist(strsplit(rows[i], split="\t"))))
    }
  }
  table <- as.data.frame(table, stringsAsFactors = F)
  names(table) <- header
  return(table)
}

#### Basic settings for cyREST ####
home.folder <- "Cluego Analysis"

cluego.home.folder = paste(home.folder,"ClueGOConfiguration","v2.5.5",sep="/")
port.number = 1234
host.address <- "localhost"

# define base urls
cytoscape.base.url = paste("http://",host.address,":", toString(port.number), "/v1", sep="")
cluego.base.url = paste(cytoscape.base.url,"apps","cluego","cluego-manager", sep="/")

#### 0.0 Start up ClueGO in case it is not running yet
response <- POST(url=paste(cytoscape.base.url,"apps","cluego","start-up-cluego",sep="/"), encode = "json")
# wait 2 seconds to make sure ClueGO is started
if(grepl("HTTP ERROR 500",response)) {
  print("wait 2 secs")
  Sys.sleep(2)
}

#### Select the ClueGO Organism to analyze ####
organism.name = "Homo Sapiens"
response <- PUT(url = paste(cluego.base.url, "organisms","set-organism",
                          URLencode(organism.name), sep="/"), encode = "json")

cluster <- paste0(filename, ".txt")

#### Upload IDs for a specific Cluster ####
# Set the number of Clusters
max.input.panel.number = 1

response <- PUT(url = paste(cluego.base.url, "cluster", "max-input-panel",
                          max.input.panel.number,sep="/"), encode = "json")

# Set analysis properties for a cluster
cluster.num = "1"

input.panel.index = cluster.num # set the cluster number
node.shape = "Ellipse" 
cluster.color = "#ff0000"
min.number.of.genes.per.term = 5
min.percentage.of.genes.mapped = 6
no.restrictions = FALSE # TRUE for no restricions in number and percentage per term
response <- PUT(url=paste(cluego.base.url, "cluster","set-analysis-properties",
                          input.panel.index, node.shape, 
                          URLencode(cluster.color,reserved = TRUE),
                          min.number.of.genes.per.term,
                          min.percentage.of.genes.mapped,
                          no.restrictions, sep="/"),
                encode = "json")


# set gene list
gene.list1 <- toJSON(read.table(cluster,as.is=TRUE,header=T)[[1]])
response <- PUT(url=paste(cluego.base.url, "cluster", "upload-ids-list", 
                          URLencode(cluster.num), sep="/"),
                body=gene.list1, encode = "json", content_type_json())


####  Select Ontologies
selected.ontologies <- toJSON(c("3;Ellipse", "5;Ellipse", "7;Ellipse",
                                "9;Ellipse", "12;Ellipse", "13;Ellipse"))
response <- PUT(url=paste(cluego.base.url, "ontologies", "set-ontologies",
                          sep="/"), body=selected.ontologies, 
                encode = "json", content_type_json())

# Use GO term significance cutoff
p.value.cutoff = 1
use.significance.cutoff = TRUE
response <- PUT(url=paste(cluego.base.url,"ontologies",use.significance.cutoff,
                          p.value.cutoff,sep="/"), encode = "json")

#### Run ClueGO Analysis ####
# Run the analysis
analysis.name <- home.folder 
# what to do if there are more than 1000 terms found
analysis.option <- "Cancel and refine selection" 
response <- GET(paste(cluego.base.url, URLencode(analysis.name),
                      URLencode(analysis.option),sep="/"))
Sys.sleep(3)

# 4.1 Get network id (SUID) (CyRest function from Cytoscape)
response <- GET(paste(cytoscape.base.url,"networks","currentNetwork",sep="/"))
current.network.suid <- content(response, encode = "json")$data$networkSUID

# Get ClueGO result table
response <- GET(paste(cluego.base.url,"analysis-results","get-cluego-table",
                      current.network.suid,sep="/"))
result.table.text <- content(response, encode = "text", encoding = "UTF-8")
result.table <- text.to.data.frame(result.table.text)
result.table <- result.table[, c("GOID", "GOTerm", "Nr. Genes",
                                 "Ontology Source", "GOGroups", 
                                 "% Associated Genes", "Associated Genes Found",
                                 "GOLevels", "Term PValue", 
                                 "Term PValue Corrected with Bonferroni step down",
                                 "Group PValue", 
                                 "Group PValue Corrected with Bonferroni step down")]
result.table$`Nr. Genes` <- as.numeric(result.table$`Nr. Genes`)
result.table$`% Associated Genes` <- as.numeric(result.table$`% Associated Genes`)
result.table <- result.table[!(duplicated(result.table$GOTerm)), ]
result.table <- result.table[order(-result.table$`Nr. Genes`, 
                                   -result.table$`% Associated Genes`),]
Sys.sleep(1)

# Remove ClueGO analysis result to reduce memory usage.
response <- DELETE(paste(cluego.base.url, "remove-cluego-analysis-result",
                         current.network.suid, sep="/"))

# Remove all ClueGO analysis results
response <- DELETE(paste(cluego.base.url, "remove-all-cluego-analysis-results",
                         sep="/"))
Sys.sleep(2)

saveRDS(result.table, paste0(filename, "_result.table.RDS"))

################################## PLOTS ##################################
# Downregulated pathways in non responders rituximab patients over time
df <- readRDS("rtx.FDR_resp.time.0.05 - Non-Responder.neg_result.table.RDS" )
df <- df[!(duplicated(df$GOTerm)), ]
df <- subset(df, df$`Term PValue Corrected with Bonferroni step down` < 0.05)
df <- df[order(df$Adjusted.pValue, decreasing = T), ]
df$y.axis <- log10(df$Adjusted.pValue)
df$Response <- "Non Responders"

# Downregulated pathways in non responders rituximab patients over time
df2 <- readRDS("rtx.FDR_resp.time.0.05 - Responder.neg_result.table.RDS")
df2 <- df2[!(duplicated(df2$GOTerm)), ]
df2 <- subset(df2, df2$`Term PValue Corrected with Bonferroni step down` < 0.05)
df2 <- df2[order(df2$Adjusted.pValue, decreasing = T), ]
df2$y.axis <- -log10(df2$Adjusted.pValue)
df2$Response <- "Responders"

df.all <- rbind(df,df2)

title.plot <- "Rituximab over time"

ggplot(df.all, aes(x=GOTerm, y=y.axis,  fill = Response)) +
  geom_bar(stat="identity") +
  coord_flip(clip = "off") +
  scale_y_continuous(breaks = c(-5, 0, 5),
                     labels = c("5", "0", "5")) +
  # direction_label(x_label = 0, y_label = 0, label_length = 6, colour = "dimgray") +
  geom_hline(yintercept = -log10(0.05), linetype='dashed', col = 'black') +
  geom_hline(yintercept =  log10(0.05), linetype='dashed', col = 'black') +
  geom_hline(yintercept = 0, col = 'black') +
  scale_fill_manual(values = c("darkblue", "deeppink4")) +
  theme_minimal() +
  theme(plot.title  = element_text(colour = 'black', size = 13, hjust = 0.7, 
                                   family = "sans", face = "bold"),
        axis.text.y = element_text(colour = 'black', size = 13, family = "sans"),
        axis.text.x = element_text(colour = 'black', size = 9, family = "sans"),
        axis.title.x= element_text(colour = 'black', size = 12, family = "sans"),
        axis.title.y = element_blank()) +
  guides(size=FALSE) + 
  labs(y = "-Log10(Adj P)",
       color = "Adj P",
       size = "Nr. Genes",
       family = "sans",
       title = title.plot)
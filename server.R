library(shiny)
library(reshape2)
# library(data.table)
library(ggplot2)

##### 
# To do:
# -- data from scer_disaggregation_recovery_long.txt DONE
# -- paper reference on user interface DONE
# -- switch to view by orf
# -- specify orf by url

rs_dt_all <- read.table("scer_disaggregation_recovery_long.txt",
                         stringsAsFactors=FALSE,header=TRUE,sep="\t")


theme_set(theme_minimal(base_size=14))

scale_y_pSup <- function(name=NULL) {
    list(scale_y_continuous(name,expand=c(0.01,0),
                       #breaks=seq(0,1,.25),labels=c("0","","0.5","","1"))
                       breaks=0:4),
         theme(axis.title.y=element_text(angle=0)))
}

timeexps <- c("42C10min.30C0min"=0,
              "42C10min.30C20min"=20,
              "42C10min.30C60min"=60,
              "42C10min.30C180min"=180)


plot_recovery_bygene <- function(mygenes=c("OLA1","FBA1"),rdata=rs_dt,
                                 type="colour",ncol=NULL,linesize=0.8,linksize=0.2,
                                 gnames=TRUE,idType=c("gene","orf")) {
    if(idType=="gene") {
        rd <- subset(rdata,gene %in% mygenes)
    } else if(idType=="orf") {
        rd <- subset(rdata,orf %in% mygenes)
    } else {
        stop("idType must be gene or orf")
    }
    mrd <- melt(rd,
                value.name = "ratio",
                id.vars=c(idType,"time"),
                measure.vars=c("ratio.pre","ratio.post"))
    maxtime <- max(mrd$time)
    mrd$gene <- factor(mrd[[idType]],levels=mygenes)
    if (type=="wrap") {
        plt_mrd <- ggplot(data=mrd,
                          aes(x=time,y=ratio)) +
            facet_wrap(~gene,ncol=ncol) +
            geom_line(aes(linetype=variable,group=variable),size=linesize)
    } else if (type=="colour") {
        plt_mrd <- ggplot(data=mrd,
                          aes(x=time,y=ratio)) +
            geom_line(aes(linetype=variable,colour=gene,group=interaction(variable,gene)),
                      size=linesize) 
        if (gnames) {
            plt_mrd <- plt_mrd + 
                geom_text(size=4,data=subset(mrd,time==maxtime & variable=="ratio.pre"),
                          aes(label=gene,colour=gene,y=ratio),hjust=0,x=maxtime+1) 
        }
    }
    plt_mrd <- plt_mrd + 
        geom_hline(yintercept=1,colour="grey50",size=0.5) + 
        geom_hline(yintercept=0,colour="grey50",size=0.5) +
        scale_linetype_manual("Amino acid label",
                              labels=c("ratio.pre"="pre-shock","ratio.post"="post-shock"),
                              values=c("ratio.pre"="solid","ratio.post"="21")) + 
        scale_x_continuous(breaks=timeexps,labels=timeexps,limits=c(0,maxtime*1.2)) + 
        labs(y="ratio in supernatant",x="time at 30C post-shock (mins)") + 
        guides(colour = "none", linetype="legend")
    
    return(plt_mrd)
}


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    output$plot <- renderPlot({
        # time plot of 46C psup
        ids <- toupper(strsplit(gsub(" ", "", input$ids, fixed = TRUE),",")[[1]])
        if (input$replicate==1) {
            rs_dt <- subset(rs_dt_all, time <= as.integer(input$maxtime) & rep=="rep1")
        } else if (input$replicate==2) {
            rs_dt <- subset(rs_dt_all, time <= as.integer(input$maxtime) & rep=="rep2")
        } 
        return(plot_recovery_bygene(ids,rs_dt,idType=input$idType))
    })
#     output$tempPlot <- renderPlot({
#         # temperature plot of 46C psup
#         genes <- strsplit(input$genes,",")[[1]]
#         plotmygenes(genes,errorbars=input$interval)$plot_temp
#     })
})
options(error=function()traceback(2))
options(shiny.maxRequestSize = 100*1024^2)


if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
requiredpackages <- c("heatmaply","STRINGdb","scales","affy","shinyjs","reshape2","gtools","orca")

 install_load <- function(packages){
   for (p in packages) {
     if (p %in% rownames(installed.packages())) {
       library(p, character.only=TRUE)
     } else {
       BiocManager::install(p)
       library(p,character.only = TRUE)
     }
  }
}
install_load(requiredpackages)
        
install_github("https://github.com/MathiasBryggerHansen/omiics_rnaseq.git",repos = BiocManager::repositories())
        
server <- function(input, output) {

  ##########################################
  ##Load annotation data
  pathway_dic <- reactive(readRDS(file = paste0("data/gene_dic_",input$species,".RDS")))
  cancer_gene_FC <- reactive(readRDS(file = paste0("data/",input$species,"_cancer.RDS"))) #FC values from E. Atlas cancer data only vs control
  neuro_gene_FC <- reactive(readRDS(file = paste0("data/",input$species,"_neuro.RDS")))
  probe_library <- reactive(readRDS(paste0("data/",input$species,"_probe_library.RDS")))
  stringdb_id <- list("human" = 9606, "zebrafish" = 7955, "rat" = 10116, "chicken" = 9031, "mouse" = 10090)
  ebi_id <- list("human" = "homo_sapiens","zebrafish" = "danio_rerio","rat" = "rattus_norvegicus", "chicken" = "gallus_gallus", "mouse" = "mus_musculus")
  string_db <- reactive(STRINGdb$new( version="11",score_threshold=200, input_directory="",species = stringdb_id[[input$species]]))
  #translate from ensembl to other gene ids
  ensembl2id <- reactive(return(readRDS(paste0("data/id_table_",input$species,".RDS"))))


  ###########################################
  ##Global variables + helper functions


  p_max <- 10**-2 #used to avoid overload
  dbs <- c("GO_Molecular_Function_2018","GO_Cellular_Component_2018","KEGG_2019_Human","Reactome_2016")

  #Runs count2deseq_analysis() or limma_analysis() for the data inputs.
  gene_results_de <- reactive({
    req(input_data$inp)
    res <- list()
    files <- input_data$inp
    for (i in 1:input$nfiles){
      counts <- files[[paste0("count",i)]]
      pheno <- files[[paste0("pheno",i)]]
      circ <- files[[paste0("circRNA",i)]]
      pheno[[2]] <- NULL
      ids <- row.names(counts)
      if(!grepl(ids[1],pattern = "ENS")){
        ids <- probe_library()$ensembl_gene_id[match(x = ids, probe_library()$probe)]
      }
      row.names(counts) <- make.names(ids,unique = T)
      if(!is.null(circ)){
        res[[paste0("circRNA",toString(i))]] <- de_circ(data = circ,pheno = pheno,i = i)
      }

      if(input[[paste0("raw_counts",i)]]){
        res[[toString(i)]] <- count2deseq_analysis(input, countdata = counts, pheno = pheno)
      }
      else {
        d0 <- DGEList(assays(temp)$counts)
        if(input$gene_filter){
          d0 <- d0[filterByExpr(d0, group=pheno),, keep.lib.sizes=FALSE]
        }
        d0 <- calcNormFactors(d0, method = "TMM")# #TMM normalization #assays(counts) - only from atlas
        #counts <- cpm(counts,log = T)
        f <- factor(pheno[[1]], levels=unique(pheno[[1]]))#phenotypes[[1]]
        if(input$batch_correction){
          batch <- pheno[[3]]
          mm <- model.matrix(~batch+f)
        }
        else{
          mm <- model.matrix(~0+f)
        }
        temp <- voom(d0, mm)
        res[[toString(i)]] <- limma_analysis(countdata = temp,phenotypes = pheno[[1]],design = mm)
      }
    }
    return(res)
  })

  #Adds data from the secondary datasets to the main.

  # input_d <- function(){
  #   files <- list()
  #   for( d in 1:input$nfiles){
  #       if(input[[paste0("combined",d)]]){#test for combined text file
  #         counts_data <- read.csv(input[[paste0("count",d)]][["datapath"]], sep = input[[paste0("sep",d)]], header = T,comment.char = "!") #comment.char = "!" in CEL files
  #       }
  #       else{
  #         d_path <- input[[paste0("count",d)]][["datapath"]]
  #         temp_dir <- paste0("./temp",d)
  #         unlink(temp_dir, recursive=TRUE,force = T) #remove if exists
  #         dir.create(temp_dir)
  #         unzip(d_path,exdir = temp_dir,overwrite = T)
  #         dir_name <- list.files(temp_dir)
  #         counts_data <- read_files(paste0(temp_dir,"/",dir_name),d)
  #       }
  #       if(input$gene_filter){
  #         counts_data <- counts_data[apply(X = counts_data,1, function(x) var(x)!=0),] #remove zero variance genes, Warning in var(x) : NAs introduced by coercion
  #       }
  #       pheno_data <- read.table(input[[paste0("phenotype",d)]][["datapath"]],header = T)#read_pheno(input[[paste0("phenotype",d)]][["datapath"]]) #phenotype data is always assumed to be tabulated, the function handles some errors in read.csv
  #       if(input$gene_id_col & input[[paste0("combined",d)]]){
  #         row.names(counts_data) <- make.names(counts_data[,1],unique = T)
  #         counts_data[,1] <- NULL
  #       }
  #       pheno_keep <- grepl(colnames(counts_data),pattern = paste(pheno_data[[2]],collapse = "|"))
  #       ids <- row.names(counts_data)
  #       if(!grepl(ids[1],pattern = "ENS")){
  #         ids <- probe_library()$ensembl_gene_id[match(x = ids, probe_library()$probe)]
  #         row.names(counts_data) <- make.names(ids,unique = T)
  #       }
  #       counts_data <- counts_data[,c(pheno_keep)]#colnames(counts_data)%in%pheno_data[[2]]] #remove samples not in pheno
  #       colnames(counts_data) <- pheno_data[[2]]
  #       counts_data <- counts_data[order(row.names(counts_data)),] #a way to secure compatible order????
  #       files[[paste0("count",d)]] <- counts_data
  #       files[[paste0("pheno",d)]] <- pheno_data
  #     }
  #   return(files)
  # }

  ############################################################################################################

  # output$volcano_out <- downloadHandler({
  #     file <- "volcano_plot.pdf"
  #   }, content = function(file){
  #     ggsave(file, plot = plotVolcano(), device = "pdf")
  #   })


  output$significant_genes <- downloadHandler({
    file <- "significant_genes.txt"
  }, content = function(file) {
    temp <- gene_results_sign()
    temp$wiki_link <- NULL
    write.table(temp, file, row.names = FALSE,sep = "\t",quote = F)
  })

  output$enrichment_analysis <- downloadHandler({
    file <- "enrichment_analysis.txt"
  }, content = function(file) {
    write.table(enrich_out(), file, row.names = FALSE,sep = "\t",quote = F)
  })


  # output$string_db <- downloadHandler({
  #     file <- "stringdb_interaction_network.pdf"
  #   }, content = function(file){
  #     ggsave(file, plot = plotStringdb(), device = "pdf")
  #   })

  observeEvent(input$start, {
    showNotification("Analysis will be running for some seconds...",type = "message",duration = 10)
  })

  # output$fileInputs <- renderUI({
  #   html_ui = " "
  #   for (i in 1:input$nfiles){#checkboxInput(inputId = paste0("CEL",i), label="Is CEL", FALSE),
  #     html_ui <- paste0(html_ui, fileInput(paste0("count",i), label=paste0("count data ",i)), fileInput(paste0("phenotype",i),
  #     label=paste0("phenotype data ",i)), checkboxInput(inputId = paste0("raw_counts",i), label="Is raw counts", TRUE),
  #     radioButtons(paste0("sep",i), "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = "\t",inline = T),
  #     checkboxInput(inputId = paste0("combined",i), label="Is combined", TRUE),
  #     checkboxInput(inputId = "gene_id_col", label = "The first column has gene ids", value = F),
  #     textInput(inputId = paste0("phen",i),label = "Phenotype id",value = "Case"))
  #   }
  #   if(input$afiles > 0){
  #     for (i in 1:input$afiles){
  #       html_ui <- paste0(html_ui, fileInput(paste0("afile",i), label=paste0("annotation file ",i)),
  #                         checkboxGroupInput(inputId = paste0("anno_range",i),choices = list("Col 2" = 2, "Col 3" = 3, "Col 4" = 4),selected = 2, label="Columns to include in annotation", TRUE),
  #                         radioButtons(paste0("sep_a",i), "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = "\t"))
  #     }
  #   }
  #
  #   HTML(html_ui)
  # })

  output$fileInputs <- renderUI({
    html_ui = " "
    for (i in 1:input$nfiles){#checkboxInput(inputId = paste0("CEL",i), label="Is CEL", FALSE),
      html_ui <- paste0(html_ui, fileInput(paste0("count",i), label=paste0("count data ",i)), fileInput(paste0("phenotype",i),
                                                                                                        label=paste0("phenotype data ",i)), checkboxInput(inputId = paste0("raw_counts",i), label="Is raw counts", TRUE),
                        radioButtons(paste0("sep",i), "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = "\t",inline = T),
                        checkboxInput(inputId = paste0("combined",i), label="Is combined", TRUE),
                        checkboxInput(inputId = "gene_id_col", label = "The first column has gene ids", value = F),
                        textInput(inputId = paste0("phen",i),label = "Phenotype id",value = "Case"))
      if(toString(i)%in%input$circRNA){
        html_ui <- paste0(html_ui, fileInput(inputId = paste0("circRNA",i), label=paste0("circRNA data ",i)))
      }
    }
    if(input$afiles > 0){
      for (i in 1:input$afiles){
        html_ui <- paste0(html_ui, fileInput(paste0("afile",i), label=paste0("annotation file ",i)),
                          checkboxGroupInput(inputId = paste0("anno_range",i),choices = list("Col 2" = 2, "Col 3" = 3, "Col 4" = 4),selected = 2, label="Columns to include in annotation", TRUE),
                          radioButtons(paste0("sep_a",i), "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = "\t"))
      }
    }

    HTML(html_ui)
  })

  #Creates interface to include circRNA
  output$circRNAfiles <- renderUI({
    HTML(paste0(selectInput(inputId = "circRNA", label = "Select the datasets with paired circRNA data",multiple = T,choices = seq(1,input$nfiles))))
  })

  input_data <- reactiveValues()

  observeEvent(input$start, {
    #tryCatch({
    input_data$inp <- input_d(input)
    #},
    #error=function(x){
    #  showNotification("The data provided does not fit the requirements, please retry",type = "message",duration = 10)
    #})
  })



  # gene_results <- reactive({##gene_symbol
  #   req(gene_results_de())
  #   data <- gene_results_de()[["1"]][["test"]]##c("padj",multiple "log2FoldChange","baseMean","pvalue","stat","lfcSE")
  #   data$ensembl_gene_id <- row.names(data)
  #
  #   data <- merge(data,ensembl2id(),by = "ensembl_gene_id",all.x = T)#ensembl2id() includes : gene_symbol, wiki_id, biotype
  #   if(input$species == "human"|input$species == "mouse"){#could be extended with a filter that takes in a TF id, and filters for genes that it interacts with
  #     tf_data <- read.table(paste0("trrust_rawdata.",input$species,".tsv"),header = T)
  #     tf_data$position <- NULL
  #     data$TF <- data$gene_symbol%in%tf_data$gene_symbol
  #   }
  #   #data$gene_symbol <- gsub(x = ifelse(is.na(data$gene_symbol)|data$gene_symbol == "",data$ensembl_gene_id,data$gene_symbol),pattern = "\\..+$", replacement = "")
  #   data <- append_secondary_data(data)#annotate log2FoldChange_i, padj_i, from the secondary datasets, also Atlas Search
  #   data <- append_annotation(data)
  #   if(input$species == "human"|input$species == "mouse"){
  #     paths <- data.frame(matrix(rep(NA,length(data$gene_symbol)*length(names(pathway_dic()))),ncol = length(names(pathway_dic()))))
  #     colnames(paths) <- names(pathway_dic())
  #     for(db in 1:length(names(pathway_dic()))){
  #         for(s in 1:length(data$gene_symbol)){
  #           paths[s,db] <- paste(pathway_dic()[[names(pathway_dic())[db]]][[data$gene_symbol[s]]],collapse = ", ")
  #       }
  #     }
  #     for(db in names(pathway_dic())){
  #      data[[db]] <- paths[,db]
  #     }
  #   }
  #   else{#mouse is could be made an exception, requires an updated dictionary
  #     for(db in names(pathway_dic())){
  #       data[[db]] <- "-"
  #     }
  #   }
  #   #wikigene_ids$gene_symbol <- NULL#maybe just not include it in file
  #   #data <- merge(data, wikigene_ids, by = "ensembl_gene_id",all.x = T)
  #   wiki_link <- sapply(data$wikigene_id, function(x) ifelse(is.na(x),"-",paste0(c("https://www.wikigenes.org/e/gene/e/",x,".html"),collapse = "")))
  #   refs <- paste0("<a href='",  wiki_link, "' target='_blank'>",wiki_link,"</a>")
  #   data$wiki_link <- refs
  #   data$wikigene_id <- NULL
  #   data <- data[order(data$padj),]
  #   return(data)
  # })
  gene_results <- reactive({##gene_symbol
    req(gene_results_de())
    data <- gene_results_de()[["1"]][["test"]]##c("padj",multiple "log2FoldChange","baseMean","pvalue","stat","lfcSE")
    return(annotate_results(input = input, data = data, ensembl2id = ensembl2id(), pathway_dic = pathway_dic(), circ=F))
  })

  gene_data <- reactiveValues()

  observeEvent(input$start, {
    gene_data$df <- gene_results()
  })

  observeEvent(input$Atlas_run,{
    req(gene_results())
    gene_data$df <- update_atlas(atlas_id = T, gene_data$df, input = input)
  })

  observeEvent(input$Atlas_search, {
    req(gene_results())
    gene_data$df <- update_atlas(gene_data$df, input = input)
  })

  observeEvent(input$k_cluster_run, {
    if(input$k_cluster > 21){
      showNotification(paste("Please choose a value between 2 and 20"),type = "message")
    }
    req(gene_results(),input$k_cluster < 21)
    counts <- gene_results_de()[["1"]][["norm_counts"]]
    counts <- data.frame(counts)
    temp <- counts[row.names(counts)%in%sign_g()$ensembl_gene_id,]
    clusters <- sapply(kmeans(temp,centers = input$k_cluster,iter.max = 100)$cluster,toString)
    res <- data.frame(t(rbind(row.names(temp),clusters)))
    colnames(res) <- c("ensembl_gene_id","cluster")
    gene_data$df <- merge(gene_data$df,res,by = "ensembl_gene_id",all.x = T)

  })

  gene_results_circ <- reactive({#tilf?je BSJ og LIN?
    req(!is.null(gene_results_de()[[paste0("circRNA",1)]][["circ_info"]]))
    d <- gene_results_de()[[paste0("circRNA",1)]][["test"]]
    d$ensembl_gene_id <- gene_results_de()[[paste0("circRNA",1)]][["circ_info"]]$ensembl_gene_id
    d$padj <- runif(n = length(d$padj), max = 0.1)/1000 ###temporary
    d$log2FoldChange <- runif(n = length(d$padj),min = -30, max = 30)
    return(annotate_results(input = input,data = d,circ=T,pathway_dic = pathway_dic(),ensembl2id = ensembl2id()))
  })

  gene_results_2 <- reactive({
    req(gene_data$df)
    gene_results <- gene_filtering()
    return(gene_results)
  })

  gene_results_filtered <- reactive({
    req(gene_data$df,eval(parse(text = input$p)))
    return(filter_results(input, gene_data$df))
  })


  gene_results_filtered_circ <- reactive({
    req(!is.null(gene_results_circ()),eval(parse(text = input$p)))
    return(filter_results(input, gene_results_circ()))
  })


  gene_results_sign <- reactive({
    req(gene_results())
    data_sign <- gene_results()[gene_results()$padj < eval(parse(text = input$p)) & abs(gene_results()$log2FoldChange) > log2(input$fc),]
    data_sign <- data_sign[!is.na(data_sign[[1]]),] #Why do columns with NA occur?
    return(data_sign)
  })

  # gene_results_de <- reactive({
  #   req(input_data$inp)
  #   res <- list()
  #   files <- input_data$inp
  #   for (i in 1:input$nfiles){
  #     counts <- files[[paste0("count",i)]]
  #     pheno <- files[[paste0("pheno",i)]]
  #     circ <- files[[paste0("circRNA",i)]]
  #     pheno[[2]] <- NULL
  #     ids <- row.names(counts)
  #     if(!grepl(ids[1],pattern = "ENS")){
  #       ids <- probe_library()$ensembl_gene_id[match(x = ids, probe_library()$probe)]
  #     }
  #     row.names(counts) <- make.names(ids,unique = T)
  #     if(!is.null(circ)){
  #       res[[paste0("circRNA",toString(i))]] <- de_circ(data = circ,pheno = pheno,i = i)
  #     }
  #
  #     if(input[[paste0("raw_counts",i)]]){
  #       res[[toString(i)]] <- count2deseq_analysis(input, countdata = counts, pheno = pheno)
  #     }
  #     else {
  #       d0 <- DGEList(assays(temp)$counts)
  #       if(input$gene_filter){
  #         d0 <- d0[filterByExpr(d0, group=pheno),, keep.lib.sizes=FALSE]
  #       }
  #       d0 <- calcNormFactors(d0, method = "TMM")# #TMM normalization #assays(counts) - only from atlas
  #       #counts <- cpm(counts,log = T)
  #       f <- factor(pheno[[1]], levels=unique(pheno[[1]]))#phenotypes[[1]]
  #       if(input$batch_correction){
  #         batch <- pheno[[3]]
  #         mm <- model.matrix(~batch+f)
  #       }
  #       else{
  #         mm <- model.matrix(~0+f)
  #       }
  #       temp <- voom(d0, mm)
  #       res[[toString(i)]] <- limma_analysis(countdata = temp,phenotypes = pheno[[1]],design = mm)
  #     }
  #   }
  #   return(res)
  # })

  output$FC_data_title_cancer <- renderUI({
    req(gene_results(),input$use_cancer)
    HTML(paste(" ",h3("logFC data in cancer"),sep = '<br/>'))
  })

  output$FC_data_text_cancer <- renderUI({
    req(gene_results(),input$use_cancer)
    s1 <- "This table shows log2(fold change) from a selection of datasets from cancer types."
    s2 <- "You can filter on the values at the top, both by values and on gene id."
    s3 <- "If you select and highlight a number of genes, an enrichment test will run on the subset."
    HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
  })

  observe({
    req(gene_results(),input$use_cancer)#input$gene_results_table_rows_all
    gene_data$fc_cancer <- merge(gene_results()[,c("ensembl_gene_id","gene_symbol","log2FoldChange","padj")],cancer_gene_FC(),by = "ensembl_gene_id",all.x = T)
  })

  output$FC_d_cancer <- renderDataTable({#show FC data for cancer dataset
    req(input$use_cancer,gene_data$fc_cancer)
    datatable(gene_data$fc_cancer, options = list(order = list(2, 'asc'),autoWidth = T,scrollX = TRUE),filter = list(position = 'top')) %>% formatRound(columns = c(3:ncol(gene_data$fc_cancer)), digits = 3)
  })
  output$FC_data_cancer <- renderUI(dataTableOutput("FC_d_cancer"))

  output$FC_data_p_cancer <- renderDataTable({
    req(input$FC_d_cancer_rows_selected,length(input$FC_d_cancer_rows_selected)>1)
    Sys.sleep(2)
    datatable(do.call("rbind", enrichr(gene_results()[input$FC_d_cancer_rows_selected,"gene_symbol"], databases = dbs)))})  ##Find a way to order so it matches with selected

  output$FC_data_pathways_cancer <- renderUI({
    req(input$FC_d_cancer_rows_selected)
    dataTableOutput('FC_data_p_cancer')
  })

  ##NEURO
  output$FC_data_title_neuro <- renderUI({
    req(gene_results(),input$use_neuro)
    HTML(paste(" ",h3("Fold change data of neuronal diseases"),sep = '<br/>'))
  })

  output$FC_data_text_neuro <- renderUI({
    req(gene_results(),input$use_neuro)
    s1 <- "This table shows log2(fold change) from a selection of datasets in neuronal diseases."
    s2 <- "You can filter on the values at the top, both by values and on gene id."
    s3 <- "If you select and highlight a number of genes, an enrichment test will run on the subset."
    HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
  })

  observe({
    req(gene_results(),input$use_neuro)#input$gene_results_table_rows_current
    gene_data$fc_neuro <- merge(gene_results()[,c("ensembl_gene_id","log2FoldChange","padj")],neuro_gene_FC(),by = "ensembl_gene_id",all.x = T)#,"gene_symbol"
  })

  output$FC_d_neuro <- renderDataTable({#show FC data for cancer dataset
    req(input$use_neuro,gene_data$fc_neuro)
    datatable(gene_data$fc_neuro, options = list(order = list(2, 'asc'),autoWidth = T,scrollX = TRUE),filter = list(position = 'top')) %>% formatRound(columns = c(3:ncol(gene_data$fc_neuro)), digits = 3)
  })
  output$FC_data_neuro <- renderUI(dataTableOutput("FC_d_neuro"))

  output$FC_data_p_neuro <- renderDataTable({
    req(input$FC_d_neuro_rows_selected,length(input$FC_d_neuro_rows_selected)>1)
    Sys.sleep(2)
    datatable(do.call("rbind", enrichr(gene_results()[input$FC_d_neuro_rows_selected,"gene_symbol"], databases = dbs)))})  ##Find a way to order so it matches with selected

  output$FC_data_pathways_neuro <- renderUI({
    req(input$FC_d_neuro_rows_selected)
    dataTableOutput('FC_data_p_neuro')
  })

  #################

  stringdb_all <- reactive({
    req(gene_results(),eval(parse(text = input$p))<p_max)#Avoid overload
    gene_filtered <- gene_results()[gene_results()$padj < eval(parse(text = input$p)) & !is.na(gene_results()$padj) & abs(gene_results()$log2FoldChange) > log2(input$fc),]
    gene_filtered <- gene_filtered[!is.na(gene_filtered$gene_symbol),]
    mapped <- string_db()$map( gene_filtered, "gene_symbol", removeUnmappedRows = TRUE ) #maybe problem to only use a subset of genes, maybe should be ordered to be sure
    pal <- scales::gradient_n_pal(colours = c("blue","grey","red"),values= c(min(mapped$log2FoldChange), mean(mapped$log2FoldChange), max(mapped$log2FoldChange)))
    payload_id <- string_db()$post_payload(mapped$STRING_id,colors=pal(gene_filtered$log2FoldChange))
    hits <- mapped$STRING_id
    return(list(string_db = string_db(), payload_id = payload_id, hits = hits))
  })

  plotStringdb <- function(){
    if(length(input$string_db_enr_rows_selected)>0){
      selected_genes <- unique(unlist(strsplit(string_db_results()[input$string_db_enr_rows_selected,"gene_symbol"],split = ",")))
    }
    else{
      selected_genes <- stringdb_all()$hits
    }
    return(stringdb_all()$string_db$plot_network(selected_genes, payload_id = stringdb_all()$payload_id) )#maybe set_background(background_vector ) should be used
  }

  output$string_db_enr_title <- renderUI({
    req(gene_results(),input$use_cancer)
    HTML(paste(" ",h3("Protein interaction enrichment results"),sep = '<br/>'))
  })
  output$string_db_enr_text <- renderUI({
    req(gene_results_filtered())
    s1 <- "This table shows a summary of a enrichment analysis in the string database for protein interactions."
    s3 <- "Similarly to the other enrichment analysis, by selecting columns with your curser, you filter on the proteins shown in the Stringdb network below."
    HTML(paste("<p>",paste(s1,s3,sep = '<br/>'),"</p>"))
  })
  string_db_results <- reactive({
    req(stringdb_all())
    temp <- stringdb_all()$string_db$get_enrichment( stringdb_all()$hits )
    temp <- temp[order(temp$fdr),c("p_value","fdr","preferredNames","description","term")]
    colnames(temp) <- c("pvalue","padj","gene_symbol","description","term")
    temp <- temp[,c("description","pvalue","padj","gene_symbol")]
    return(temp)})

  output$string_db_enr <- renderDataTable({
    req(string_db_results())
    datatable(string_db_results())#how to return p-values, and avoid overlap with plot
  })


  output$volcano_title <- renderUI({
    req(gene_results_filtered())
    HTML(paste(" ",h3(paste0("Volcano plot of ", input$phen1)), sep = '<br/>'))
  })

  output$volcano_text <- renderUI({
    req(gene_results_filtered())
    s0 <- "The volcano plot shows a summary of the differential expression regression analysis."
    s1 <- "The plot can be colored by typing a search term for a given column in one of the text boxes to the left."
    s3 <- "You can also export the figure in svg by pressing the camera symbol in the top right corner. If you click on items in the legend they are removed from view."
    s2 <- "You can interact with the plot to zoom in on an area when the zoom option is selected, and you can select datapoints of interest. This will show details on the specific genes, and a boxplot will be plotted showing normalized gene expressions."

    HTML(paste("<p>",paste(s0,s1, s2, s3, sep = '<br/>'),"</p>"))
  })

  output$volcano <- renderPlotly({
    req(gene_results_filtered(), eval(parse(text = input$p))<p_max)
    colnames(gene_results_filtered())
    volcano_plot(input = input, data = gene_results_filtered(), pathway_dic = pathway_dic())})

  output$volcano_title_circ <- renderUI({
    req(gene_results_filtered_circ())
    HTML(paste(" ",h3(paste0("Volcano plot of circRNA for ", input$phen1)), sep = '<br/>'))
  })

  output$volcano_text_circ <- renderUI({
    req(gene_results_filtered_circ(),!is.null(input_data$inp[[paste0("circRNA",1)]]))
    s0 <- "The volcano plot shows a summary of the differential expression regression analysis."
    s1 <- "The plot can be colored by typing a search term for a given column in one of the text boxes to the left."
    s3 <- "You can also export the figure in svg by pressing the camera symbol in the top right corner. If you click on items in the legend they are removed from view."
    s2 <- "You can interact with the plot to zoom in on an area when the zoom option is selected, and you can select datapoints of interest. This will show details on the specific genes, and a boxplot will be plotted showing normalized gene expressions."

    HTML(paste("<p>",paste(s0,s1, s2, s3, sep = '<br/>'),"</p>"))
  })

  output$volcano_c <- renderPlotly({
    req(gene_results_filtered_circ(), eval(parse(text = input$p))<p_max,!is.null(input_data$inp[[paste0("circRNA",1)]]))
    res <- gene_results_filtered_circ()
    volcano_plot(input = input, res, pathway_dic = pathway_dic())
  })


  output$volcano_circ <- renderUI({
    req(gene_results_filtered_circ(), eval(parse(text = input$p))<p_max,!is.null(input_data$inp[[paste0("circRNA",1)]]))
    plotlyOutput("volcano_c",height = "1000px")
  })

  circToLin <- reactive({
    req(gene_results_de())
    a <- gene_results_de()[[paste0("circRNA",1)]][["test"]]
    d <- gene_results_de()[[paste0("circRNA",1)]][["circ_info"]]
    d <- cbind(a,d)
    d <- filter_results(input, d)
    d <- merge(d,ensembl2id(), by = "ensembl_gene_id")
    return(d)
  })

  output$circVsL <- renderPlotly({
    req(circToLin())
    ggplot(circToLin(), aes(x=log2(sum_lin), y=log2(sum_junction),color=gene_biotype)) + geom_point() + geom_abline(slope = 1)
  })

  output$circVsLin <- renderUI({
    req(gene_results_filtered_circ(), eval(parse(text = input$p))<p_max,!is.null(input_data$inp[[paste0("circRNA",1)]]))
    plotlyOutput("circVsL",height = "1000px")
  })

  output$pca_title <- renderUI({
    req(gene_results())
    HTML(paste(" ", h3(paste0("PCA of ",input$phen1)), sep = '<br/>'))
  })
  output$pca_text <- renderUI({
    req(gene_results())
    s1 <- "The PCA plot shows the first three principle components discriminating each sample."
    s2 <- "By using your curser, you can zoom in and out and rotate the view."
    s3 <- "You can export the figure when you have found a fitting angle."
    HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
  })

  output$pca <- renderPlotly({
    req(gene_results_de())
    if(input[["raw_counts1"]]){
      data_norm <- gene_results_de()[["1"]][["dds"]]
      data_norm <- data.frame(assay(varianceStabilizingTransformation(data_norm)) )
      data_norm <- data_norm[rowSums(data_norm)>50,]
    }
    else{
      data_norm <- gene_results_de()[["1"]][["norm_counts"]]
      data_norm <- data_norm[rowSums(data_norm)>50,] #another solution?
    }
    data_norm$ensembl_gene_id <- row.names(data_norm)
    data_norm$ensembl_gene_id <- NULL
    data_norm$Name <- NULL
    data_norm$Description <- NULL
    pheno <- input_data$inp[[paste0("pheno",1)]]
    df_pca  <- prcomp(t(data_norm),scale = T, center = T)
    df_out <- as.data.frame(df_pca$x)
    scores <- df_pca$x
    if(input$pca_pheno > nrow(pheno)){
      showNotification("You need to specify a valid column number",type = "message")
      ann <- pheno[[1]]
    }
    else{
      ann <- pheno[[input$pca_pheno]]
    }
    plot_ly(x=scores[,1], y=scores[,2], z=scores[,3], type = "scatter3d", mode="markers",name = row.names(df_out), color = ann) %>%
      layout(scene = list(xaxis = list(title = "PC1"), yaxis = list(title = "PC2"), zaxis = list(title = "PC3"))) %>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "pca",
          width = 1500,
          height = 1000
        )
      )
  })

  output$gene_results_title <- renderUI({
    req(gene_results())
    HTML(paste(" ", h3("Table of DE results"), sep = '<br/>'))
  })
  output$gene_results_text <- renderUI({
    req(gene_results())
    s1 <- "The table includes analysis results of differential expression along with annotation."
    s2 <- "It can be filtered by setting filter parameters for each column. You can for example change the range of log2FoldChange or search for a pathway of interest."
    s3 <- "The table can also be filtered by inserting a list comma separated gene ids (ensembl or hgnc symbol) in the Post Analysis section. This will search for protein interactions with the given genes, and filter on that given list."
    s4 <- "You can add new datasets by searching the Atlas Database in the Post Analysis. These can be plotted or used as annotation."
    HTML(paste("<p>",paste(s1, s2, s3, s4,sep = '<br/>'),"</p>"))
  })

  output$gene_results_table <- renderDataTable({
    req(gene_results())
    datatable(data = gene_results_2(),caption = "DE results including any extra datasets",filter = list(position = 'top'),escape = F, options = list(autoWidth = TRUE,scrollX = TRUE, columnDefs = list(list(width = '400px', targets = grep(colnames(gene_data$df),pattern = "reactome|kegg|carta|wiki")))))
  })

  output$gene_results_table_circ <- renderDataTable({
    req(gene_results_circ())
    datatable(data = gene_results_circ(),caption = "DE results including any extra datasets",filter = list(position = 'top'),escape = F, options = list(autoWidth = TRUE,scrollX = TRUE, columnDefs = list(list(width = '400px', targets = grep(colnames(gene_data$df),pattern = "reactome|kegg|carta|wiki")))))
  })


  gene_filtering <- reactive({
    if(input$gene_interaction == ""){
      return(gene_data$df)
    }
    else{
      #date_time <- Sys.time()
      #while((as.numeric(Sys.time()) - as.numeric(date_time))<10){} #Sys.sleep does not work
      showNotification("Filtering for protein interactions...",type = "message") #is this working?
      gene_ids <- unlist(strsplit(gsub(" ", "", input$gene_interaction, fixed = TRUE),","))

      if(grepl(gene_ids[1],pattern = "^ENS")){
        gene_ids <- gene_data$df[gene_data$df$ensembl_gene_ids%in%gene_ids,"gene_symbol"]
      }
      gene_network <- read.table(url(paste0("https://string-db.org/api/tsv/interaction_partners?identifiers=",paste(gene_ids,collapse = "%0d"),"&species=",stringdb_id[[input$species]],"&required_score=950&limit=1000")),header = T)$preferredName_B
      gene_network <- c(gene_ids,gene_network)
      return(gene_data$df[gene_data$df$gene_symbol%in%gene_network,])
    }
  })

  output$sign_pheno_exp_title <- renderUI({
    req(gene_results_filtered())
    HTML(paste(" ", h3("Table of matching datasets in the Expression Atlas database according to significant hits"), sep = '<br/>'))
  })

  output$sign_pheno_exp_text <- renderUI({
    req(gene_results())
    s1 <- "This table shows datasets from the Expression Atlas database, which shows a significant correlation with the primary dataset."
    s2 <- "You can include one or more of the datasets by inserting their IDs comma seperated in one of the search fields to the left."
    s3 <- "The added datasets can for example be used to annotate your primary dataset, or you can create a new volcano plot with the input."
    HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
  })
  output$sign_pheno_exp <- renderDataTable({#gene_filtered <- gene_results_filtered()[gene_results_filtered()$padj < eval(parse(text = input$p)) & !is.na(gene_results_filtered()$padj) & abs(gene_results_filtered()$log2FoldChange) > log2(input$fc),]
    req(gene_results(),eval(parse(text = input$p))<p_max)
    if(length(sign_g()$ensembl_gene_id)>100){
      sign_genes <- sign_g()$ensembl_gene_id[1:100] #top 100 to avoid overload
    }
    else{
      sign_genes <- sign_g()$ensembl_gene_id
    }
    datatable(testGenes(sign_genes, ebi_id, input = input),escape = F)
  })

  sign_g <- reactive({ #lowest padj first
    req(gene_results())
    sign_genes <- list()
    temp <- gene_results()[gene_results()$padj < eval(parse(text = input$p)) & !is.na(gene_results()$padj) & abs(gene_results()$log2FoldChange) > log2(input$fc),]
    temp <- temp[order(temp$padj,decreasing = F),]
    sign_genes$gene_symbol <- temp[,"gene_symbol"]
    sign_genes$ensembl_gene_id <- temp[,"ensembl_gene_id"]
    return(sign_genes)
  })

  ##FC cancer comparison heatmap

  output$heatmap_cancer_cor_title <- renderUI({
    req(gene_results(),input$use_cancer)
    HTML(paste(" ", h3("Fold change heatmap between cancer datasets"), sep = '<br/>'))
  })
  output$heatmap_cancer_cor_text <- renderUI({
    req(gene_results(),input$use_cancer)
    s1 <- "This heatmap shows a collective selection of datasets, showing significant fold change in a subset of genes."
    s2 <- "Non-significant and missing values are set to zero. The fold change values are capped at 3 to give a better visual interpretation."
    s3 <- "You can always hover your curser over a datapoint to get information about the x and y axis."
    HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
  })

  output$heatmap_cancer_cor <- renderUI({
    req(gene_results(),input$use_cancer)
    plotlyOutput("heatmap_cancer_c", height = "1000px",width = "1500px")#plotlyOutput("heatmap_cancer_c")
  })

  output$heatmap_cancer_c <- renderPlotly({
    req(gene_data$fc_cancer,gene_results(),input$use_cancer)
    #data <- gene_data$df[,c("ensembl_gene_id","log2FoldChange","gene_symbol")] #should be adjusted if deseq?
    #colnames(data)[2] <- input$phen1
    canc_FC <- gene_data$fc_cancer#merge(data,cancer_gene_FC(),by = "ensembl_gene_id")
    canc_FC$ensembl_gene_id <- NULL
    canc_FC$ensembl_gene_id_human <- NULL
    colnames(canc_FC)[colnames(canc_FC)=="log2FoldChange"] <- input$phen1
    canc_FC$padj <- NULL
    row.names(canc_FC) <-  make.names(canc_FC$gene_symbol,unique = T)
    canc_FC$gene_symbol <- NULL
    #canc_FC[is.na(canc_FC)] <- 0
    canc_FC[canc_FC>3] <- 3
    canc_FC[canc_FC<(-3)] <- -3
    if(input$chain){
      canc_FC <- canc_FC[input$FC_d_cancer_rows_all[1:100],]
    }
    else{
      canc_FC <- canc_FC[order(rowSums(abs(canc_FC),na.rm = T),decreasing = T)[1:100],]
    }
    heatmaply(as.matrix(canc_FC),scale_fill_gradient_fun = scale_fill_gradient( low = input$col_low, high = input$col_high, na.value = "grey")) %>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "heatmap_dataset_comparison",
          width = 2000,
          height = 1500
        ))
  })

  ####NEURO

  output$heatmap_neuro_cor_title <- renderUI({
    req(gene_results(),input$use_neuro)
    HTML(paste(" ", h3("Fold change heatmap from cancer datasets"), sep = '<br/>'))
  })
  output$heatmap_neuro_cor_text <- renderUI({
    req(gene_results(),input$use_neuro)
    s1 <- "This heatmap shows a collective selection of datasets, showing significant fold change in a subset of genes."
    s2 <- "Non-significant and missing values are set to zero. The fold change values are capped at 3 to give a better visual interpretation."
    s3 <- "You can always hover your curser over a datapoint to get information about the x and y axis."
    HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
  })

  output$heatmap_neuro_cor <- renderUI({
    req(gene_results(),input$use_neuro)
    plotlyOutput("heatmap_neuro_c", height = "1000px",width = "1500px")#plotlyOutput("heatmap_neuro_c")
  })

  output$heatmap_neuro_c <- renderPlotly({
    req(gene_data$fc_neuro,gene_results(),input$use_neuro)
    #data <- gene_data$df[,c("ensembl_gene_id","log2FoldChange","gene_symbol")] #should be adjusted if deseq?
    #colnames(data)[2] <- input$phen1
    neuro_FC <- gene_data$fc_neuro#merge(data,neuro_gene_FC(),by = "ensembl_gene_id")
    neuro_FC$ensembl_gene_id <- NULL
    neuro_FC$ensembl_gene_id_human <- NULL
    colnames(neuro_FC)[colnames(neuro_FC)=="log2FoldChange"] <- input$phen1
    neuro_FC$padj <- NULL
    row.names(neuro_FC) <-  make.names(neuro_FC$gene_symbol,unique = T)
    neuro_FC$gene_symbol <- NULL
    #neuro_FC[is.na(neuro_FC)] <- 0
    neuro_FC[neuro_FC>3] <- 3
    neuro_FC[neuro_FC<(-3)] <- -3
    if(input$chain){
      neuro_FC <- neuro_FC[input$FC_d_neuro_rows_all[1:100],]
    }
    else{
      neuro_FC <- neuro_FC[order(rowSums(abs(neuro_FC),na.rm = T),decreasing = T)[1:100],]
    }
    heatmaply(as.matrix(neuro_FC),scale_fill_gradient_fun = scale_fill_gradient( low = input$col_low, high = input$col_high, na.value = "grey")) %>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = "heatmap_dataset_comparison",
          width = 1500,
          height = 2000
        ))
  })

  ##


  observeEvent(event_data("plotly_selected"), {
    req(!is.null(event_data("plotly_selected")))
    d <- event_data("plotly_selected")
    d$curveNumber <- NULL
    d$pointNumber <- NULL
    colnames(d)[1:2] <- c("logFC","-log10P")
    d$key[grep(d$key,pattern = ";$")] <- paste0(d$key[grep(d$key,pattern = ";$")],".")
    ids <- data.frame(sapply(d$key, function(x) strsplit(x, split = ";")))
    ids <- data.frame(matrix(unlist(ids), nrow=length(ids), byrow=T))
    d$key <- NULL
    d$ensembl_gene_id <- as.vector(unlist(ids[,1]))
    d$gene_symbol <- as.vector(unlist(ids[,2]))
    d$gene_biotype <- as.vector(unlist(ids[,3]))
    d <- append_secondary_data(input, d)
    d <- append_annotation(input, d)
    i <- ncol(d)
    i <- ncol(ids) - length(names(pathway_dic()))
    for(p in names(pathway_dic())){
      i <- i + 1
      d[p] <- as.vector(unlist(ids[,i]))
    }
    d <- d[order(d$ensembl_gene_id),]
    gene_set <- d$ensembl_gene_id
    data <- gene_results_de()[["1"]]$norm_counts
    count_sub <- data.frame(data[row.names(data)%in%gene_set,])
    count_sub <- count_sub[order(row.names(count_sub)),]
    row.names(count_sub) <- d$gene_symbol
    d2 <- data.frame(t(count_sub))
    pheno <- input_data$inp[[paste0("pheno",1)]]
    d2$phenotype <- pheno[[1]]##subset count data for gene ids?
    d2 <- melt(data = d2, id = c("phenotype"))

    output$volcano_selected <- renderUI({
      dataTableOutput("volcano_s")
    })

    output$volcano_s <- renderDataTable({
      #req(d)
      datatable(d,escape = F,caption = "Information on genes selected from volcano plot",options = list(autoWidth = TRUE,columnDefs = list(list(width = '400px', targets = grep(colnames(d),pattern = "reactome|kegg|carta|wiki"))),scrollX = TRUE))# %>% formatRound(columns = c(6:ncol(d)), digits = 3)})
      #}, sanitize.text.function = function(x) x, digits = 3)#, sanitize.text.function = function(x) x)
    })

    # output$delta_cor_title <- renderUI({
    #   req(gene_results_de())
    #   HTML(paste(" ",h3("Delta values of correlation between selected gene expressions (cor(case) - cor(control))"),sep = '<br/>'))
    # })
    # output$delta_cor <- renderUI({
    #   plotlyOutput("delta_c")
    # })
    #
    # output$delta_c <- renderPlotly({ #remove
    #   pheno <- gene_results_de()[["1"]][["phenotypes"]] #$phenotype#ifelse(grepl(colnames(our_exp),pattern = "^H"),"Control", "Case") #another way to define phenotype needs to be added
    #   ref <- names(table(pheno))[grepl(names(table(pheno)),pattern = "control|normal|reference",ignore.case = T)]
    #   cases <- count_sub[,pheno != ref]
    #   control <- count_sub[,pheno == ref]
    #   delta_cor <- cor(t(cases),method = "spearman") - cor(t(control),method = "spearman")
    #   heatmaply(as.matrix(delta_cor),Rowv = F,Colv = F,scale_fill_gradient_fun = scale_fill_gradient(low = input$col_low, high = input$col_high)) %>%
    #     config(
    #       toImageButtonOptions = list(
    #         format = "svg",
    #         filename = "delta_cor",
    #         width = 1500,
    #         height = 1000
    #       ))
    # })

    output$boxplot_search <- renderUI({
      plotlyOutput("boxplot_s")
    })

    output$boxplot_search_text <- renderUI({
      req(gene_results_filtered())
      s1 <- "The boxplot shows the expression of the genes selected."
      s2 <- "You can export the figure in svg by pressing the camera symbol."
      s3 <- "By hovering over the plot you can see the values for each gene."
      HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
    })

    output$boxplot_s <- renderPlotly({ #show count data for the selected data (volcano plot)
      #req(input$data_search_rows_current) #create boxplots of other continues variables like methylation

      ggplotly(ggplot(d2, aes(x=variable, y=value,color=phenotype)) + geom_boxplot() + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + xlab("gene symbol") + ylab("scaled counts") + ggtitle("Boxplot of expression values of genes selected from volcano plot")) %>%
        config(
          toImageButtonOptions = list(
            format = "svg",
            filename = "boxplot_expression",
            width = 1500,
            height = 1000
          ))
    })
  })

  #Pathways Enrichment
  ######################################

  enrich_out <- reactive({ #For all imported datasets: run enrichR, merge them
    req(gene_results())
    temp <- gene_results()
    res <- list()
    for(i in 1:input$nfiles){#create exception if not human
      if(i>1){
        padj_id <- paste0("padj_",input[[paste0("pheno",i)]])
        temp2 <- temp[,grepl(colnames(temp),pattern = paste(padj_id,"gene_symbol",sep = "|"))] #,log2FC_id
        temp2 <- temp2[,order(colnames(temp2))]
        colnames(temp2) <- c("gene_symbol","padj") #,"log2FoldChange"
      }
      else{
        temp2 <- temp[,c("gene_symbol","padj")]
      }
      temp2 <- temp2[order(temp2$padj,decreasing = F),]
      temp2 <- temp2[!is.na(temp2$padj),][1:1000,]
      sign_genes <- temp2[temp2$padj < 0.05,"gene_symbol"] #input$p
      df <- do.call("rbind", enrichr(sign_genes, databases = dbs))
      df <- df[!duplicated(df$Term),]#for some reason there are duplicates copies
      df <- df[order(df$Adjusted.P.value,decreasing = F),]
      colnames(df) <- c("Term","Overlap","pvalue","padj","oldp","old adjp","OR","Score","gene_symbol")
      df$gene_symbol <- gsub(df$gene_symbol,pattern = ";",replacement = ",")
      df <- df[,c("Term","Overlap","pvalue","padj","OR","gene_symbol")]
      res[[input[[paste0("phen",i)]]]] <- df #if pheno id is equal it overwrites, maybe check
    }
    for(n in names(res)){
      if(!exists("df2")){
        df2 <- res[[n]]
      }
      else{
        d3 <- res[[n]]
        df3 <- d3[,c("Term","padj","gene_symbol")]
        colnames(df3)[-1] <- paste(c("padj","gene_symbol"),n)
        df2 <- merge(df2,df3,"Term",all = T)
      }
    }
    return(df2)
  })

  output$sign_pathways_title <- renderUI({
    req(gene_results())
    HTML(paste(" ",h3("Pathway enrichment analysis"), sep = '<br/>'))
  })
  output$sign_pathways <- renderDataTable({
    req(gene_results())
    enrich_out <- enrich_out()
    datatable(enrich_out)
  })
  output$sign_pathways_text <- renderUI({
    req(gene_results_filtered())
    s1 <- "This table shows a summary of a enrichment analysis by the R package enrichR including GO molecular function, GO cellular function and Reactome and KEGG pathways."
    s2 <- "You can search for a given GO, pathway or related phenotype in the search field. The table can be downloaded by interacting with the left interface."
    s3 <- "By selecting columns with your curser, you filter on the proteins shown in the Stringdb network below."
    HTML(paste("<p>",paste(s1, s2,s3,sep = '<br/>'),"</p>"))
  })

  ##Gene to gene correlation of significant hits in main dataset

  output$gene_gene_title <- renderUI({
    req(gene_results())
    HTML(paste(" ",h3("Gene to gene correlation of selected genes in main dataset"), sep = '<br/>'))
  })
  output$gene_gene <- renderPlotly({
    req(gene_results())
    data <- data.frame(gene_results_de()[["1"]]$norm_counts)
    sample_ids <- colnames(data)
    if(sum(grepl(row.names(data), pattern = "ENS"))> nrow(data)/2){#check if gene ids are ensembl, otherwise assume hgnc
      data$ensembl_gene_id <- row.names(data)
      data <- merge(data,ensembl2id(),by = "ensembl_gene_id",all.x = T) #also adds gene_biotype
    }
    else {
      data$gene_symbol <- row.names(data)
      data <- merge(data,ensembl2id(),by = "gene_symbol",all.x = T)
    }
    row.names(data) <- make.names(ifelse(is.na(data$gene_symbol)|data$gene_symbol == "",data$ensembl_gene_id,data$gene_symbol),unique = T)

    if(input$chain){
      #data <- data[row.names(data)%in%gene_set,]
      gene_symbol <- gene_results()[input$gene_results_table_rows_current,"gene_symbol"]
      #gene_symbol <- gene_results()[input$gene_results_table_rows_current,"gene_symbol"]
      data <- data[row.names(data)%in%gene_symbol,]#gene_results_table_rows_all is not ordered like gene_results, fix?

    }
    else{
      gene_symbol <- unique(gene_results()[order(gene_results()$padj,decreasing = F)[1:200],"gene_symbol"])
      #gene_symbol <- gene_results()[order(gene_results()$padj,decreasing = F)[1:200],"gene_symbol"]
      data <- data[row.names(data)%in%gene_symbol,]
      #data <- merge(data,ensembl2idl(),by = "ensembl_gene_id",all.x = T)
      #data$gene_symbol <- gene_symbol[ensembl%in%row.names(data)] #ifelse(is.na(data$gene_symbol)|data$gene_symbol == "",data$ensembl_gene_id,data$gene_symbol)
    }

    #data <- merge(data,,by = "ensembl_gene_id",all.x = T)
    #gene_symbol[row.names(data)%in%ensembl]
    data <- data[,sample_ids]
    data <- data[apply(X = data,1, function(x) var(x)!=0),]
    heatmaply(cor(data.frame(t(data)),method = "spearman"),scale_fill_gradient_fun = scale_fill_gradient(low = input$col_low, high = input$col_high)) %>%
      config(toImageButtonOptions = list(format = "svg",
                                         filename = "gene_gene_interaction",
                                         width = 2500,
                                         height = 2500))
  })
  output$gene_gene_text <- renderUI({
    req(gene_results_filtered())
    s1 <- "Gene-gene expression correlation (spearman) of a selection of genes. The default uses the DE tophits. By using the chaining option in the left section the genes shown will be filtered by the current rows shown in the DE table."
    s2 <- "You can export the figure in svg by pressing the camera symbol."
    s3 <- "By hovering over the plot you can see the values for each gene pair."
    HTML(paste("<p>",paste(s1, s2, s3,sep = '<br/>'),"</p>"))
  })

  output$gene_sample_title <- renderUI({
    req(gene_results())
    HTML(paste(" ",h3("Gene to sample correlation of selected genes in main dataset"), sep = '<br/>'))
  })
  output$gene_sample <- renderPlotly({
    req(gene_results())
    data <- data.frame(gene_results_de()[["1"]]$norm_counts)
    pheno <- data.frame(gene_results_de()[["1"]]$phenotypes)
    #names(pheno) <- "phenotype"
    sample_ids <- colnames(data)
    if(sum(grepl(row.names(data), pattern = "ENS"))> nrow(data)/2){#check if gene ids are ensembl, otherwise assume hgnc
      data$ensembl_gene_id <- row.names(data)
      data <- merge(data,ensembl2id(),by = "ensembl_gene_id",all.x = T) #also adds gene_biotype
    }
    else {
      data$gene_symbol <- row.names(data)
      data <- merge(data,ensembl2id(),by = "gene_symbol",all.x = T)
    }
    row.names(data) <- make.names(ifelse(is.na(data$gene_symbol)|data$gene_symbol == "",data$ensembl_gene_id,data$gene_symbol),unique = T)

    if(input$chain){
      #data <- data[row.names(data)%in%gene_set,]
      gene_symbol <- gene_results()[input$gene_results_table_rows_current,"gene_symbol"]
      #gene_symbol <- gene_results()[input$gene_results_table_rows_current,"gene_symbol"]
      data <- data[row.names(data)%in%gene_symbol,]#gene_results_table_rows_all is not ordered like gene_results, fix?

    }
    else{
      gene_symbol <- unique(gene_results()[order(gene_results()$padj,decreasing = F)[1:200],"gene_symbol"])
      #gene_symbol <- gene_results()[order(gene_results()$padj,decreasing = F)[1:200],"gene_symbol"]
      data <- data[row.names(data)%in%gene_symbol,]
      #data <- merge(data,ensembl2idl(),by = "ensembl_gene_id",all.x = T)
      #data$gene_symbol <- gene_symbol[ensembl%in%row.names(data)] #ifelse(is.na(data$gene_symbol)|data$gene_symbol == "",data$ensembl_gene_id,data$gene_symbol)
    }

    #data <- merge(data,,by = "ensembl_gene_id",all.x = T)
    #gene_symbol[row.names(data)%in%ensembl]
    data <- data[,sample_ids]
    names(pheno) <- "phenotype"
    data <- data[apply(X = data,1, function(x) var(x)!=0),]
    heatmaply(data,scale_fill_gradient_fun = scale_fill_gradient(low = input$col_low, high = input$col_high), col_side_colors = pheno,distfun = "spearman") %>%
      config(toImageButtonOptions = list(format = "svg",
                                         filename = "gene_sample_interaction",
                                         width = 2500,
                                         height = 2500))
  })



  #onclick("button", {
  observe({
    req(sign_g(),eval(parse(text = input$p)) < p_max)
    if(length(input$string_db_enr_rows_selected)>0|length(input$sign_pathways_rows_selected)>0){
      selected_genes <- unique(unlist(strsplit(string_db_results()[input$string_db_enr_rows_selected,"gene_symbol"],split = ",")))
      selected_genes <- c(selected_genes, unique(unlist(strsplit(enrich_out()[input$sign_pathways_rows_selected,"gene_symbol"],split = ","))))
      selected_genes <- unique(selected_genes)
      ens_id <- gene_results()[gene_results()$gene_symbol%in%selected_genes,"ensembl_gene_id"] #fix to convert symbol to ENS
      selected_exp <- gene_results_de()[["1"]]$norm_counts[row.names(gene_results_de()[["1"]]$norm_counts)%in%ens_id,]
      row.names(selected_exp) <- selected_genes[selected_genes%in%gene_results()$gene_symbol]
      selected_exp <- data.frame(t(selected_exp))
      selected_exp$phenotype <- gene_results_de()[["1"]][["phenotypes"]] #$phenotype
      d2 <- melt(selected_exp,id = "phenotype")
      output$stringdb_b <-  renderPlotly({ggplotly(ggplot(d2, aes(x=variable, y=value,color=phenotype)) + geom_boxplot() +
                                                     theme(axis.text.x=element_text(angle = -90, hjust = 0)) + xlab("gene symbol") + ylab("scaled counts") +
                                                     ggtitle("Boxplot of genes in stringdb interactions")) %>%
          config(
            toImageButtonOptions = list(
              format = "svg",
              filename = "boxplot_network_expression",
              width = 1500,
              height = 1000
            ))})
      output$stringdb_box <- renderUI({plotlyOutput("stringdb_b")})
      js$loadStringData(selected_genes,stringdb_id[[input$species]])
    }
    else{
      js$loadStringData(sign_g()$gene_symbol,stringdb_id[[input$species]])
    }
  })
  output$stringdb_text <- renderUI({
    req(gene_results())
    HTML(paste(" ",p("String db protein interaction network, as default this will include all the significant hits defined by the user."), sep = '<br/>'))
  })
  output$stringdb_title <- renderUI({
    req(gene_results())
    HTML(paste(" ",h3("Stringdb network"), sep = '<br/>'))
  })

  #output$getPath <- getPathway('WP554')
  #output$reactome <- renderPlot({

  #  fc <- gene_results_de()$log2FoldChange[!duplicated(gene_results_de()$ENTREZID)]
  #  names(fc) <- gene_results_de()$ENTREZID[!duplicated(gene_results_de()$ENTREZID)]
  #  viewPathway("mRNA Splicing - Major Pathway", readable=TRUE, foldChange = fc) #it quickly becomes chaotic, with too many interactions
  #})

}

ui <- fluidPage(

  titlePanel("Omiics NGS Analysis"),
  sidebarLayout(
    sidebarPanel(
      h2("Data Upload"),
      numericInput("nfiles", "Number of paired datasets", value = 1, min = 1, step = 1),
      p("The first dataset will be treated as your primary dataset"),
      numericInput("afiles", "Number of annotation datasets", value = 0, min = 0, step = 1),
      uiOutput("circRNAfiles"),
      uiOutput("fileInputs"),

      tags$hr(),
      #checkboxGroupInput(inputId = "dbs", "Pathways to use in enrichment analysis:",c("GO Molecular Function" = "GO_Molecular_Function_2018","GO Cellular Component"="GO_Cellular_Component_2018","KEGG"="KEGG_2019_Human","Reactome"="Reactome_2016"),selected = "KEGG_2019_Human"),
      radioButtons("species", "Species",choices = c(human = "human",mouse = "mouse"),selected = "human",inline = T),
      checkboxInput(inputId = "batch_correction", label = "Use batch correction", value = F),
      p("This requires a batch annotation in the third column of the phenotype data"),
      checkboxInput(inputId = "gene_filter", label = "Use filter to pre-exclude genes with no variance in expression", value = T),
      textInput(inputId = "p", label = "p cutoff:",value = '10**-14'),
      numericInput(inputId = "fc", label = "FC cutoff:", value = 2),

      tags$hr(),
      h2("Secondary Parameters"),
      checkboxInput(inputId = "use_cancer", label = "Use summary stats from cancers", value = F),
      checkboxInput(inputId = "use_neuro", label = "Use summary stats from neurological diseases", value = F),
      checkboxInput(inputId = "chain", label = "Chain output with filtered table", value = F),
      textInput(inputId = "experiment_id", label = "Use alternative dataset for volcano plot", value = ""),
      textInput(inputId = "volcano_col", label = "A column to annotate color in volcano plot",value = "gene_biotype"),
      textInput(inputId = "col_high", label = "Volcano color of continuous high values", value = "red"),
      textInput(inputId = "col_low", label = "Volcano color of continuous low values", value = "blue"),
      #checkboxInput(inputId = "pca_pheno", label = "Annotate PCA with phenotype", value = F),
      numericInput(inputId = "pca_pheno", label = "PCA annotation column", value = 1),
      checkboxInput(inputId = "log_scale", label = "Scale color values"),
      #,rat = "rat",chicken = "chicken",zebrafish = "zebrafish"
      actionButton(inputId = "start", label = "Start Analysis!"),

      tags$hr(),
      h2("Post Analysis"),
      numericInput(inputId = "k_cluster", label = "Kmeans clustering", value = 6),
      actionButton(inputId = "k_cluster_run", label = "Run kmeans cluster!"),
      textInput("gene_interaction", label = "Comma separated gene ids for Stringdb interaction search",value = ""),
      #actionButton(inputId = "string_gene_filter", label = "Filter on gene ids!"),
      textInput("atlas_query", label = "Enter Atlas phenotype query", value = ""),
      p("Search ebi expression atlas for data related to the search term"),
      actionButton(inputId = "Atlas_search", label = "Search Atlas database!"),
      textInput("Atlas_ids", label = "Enter Atlas ids", value = ""),
      p("Comma separated ebi dataset IDs"),
      actionButton(inputId = "Atlas_run", label = "Add Atlas dataset(s)!"),

      tags$hr(),
      downloadButton(outputId = "significant_genes", label = "Download DE genes"),
      #downloadButton(outputId = "volcano_out",label = "Download volcano plot"),
      downloadButton(outputId = "enrichment_analysis",label = "Download enrichment analysis"),
      downloadButton(outputId = "download_filtered",label = "Download filtered gene table"),
      #downloadButton(outputId = "string_db",label = "Download stringdb network")
    ),
    #####################################
    # Main panel for displaying outputs ----
    mainPanel(
      fluidRow(
        htmlOutput(outputId = "pca_title"),
        plotlyOutput(outputId = "pca",height = "1000px"),
        htmlOutput(outputId = "pca_text"),

        htmlOutput(outputId = "volcano_title"),
        plotlyOutput(outputId = "volcano",height = "1000px"),
        htmlOutput(outputId = "volcano_text"),

        uiOutput("volcano_selected"),
        uiOutput("boxplot_search"),
        htmlOutput(outputId = "boxplot_search_text"),

        # htmlOutput("delta_cor_title"),
        # uiOutput("delta_cor"),

        htmlOutput(outputId = "gene_results_title"),
        dataTableOutput("gene_results_table"),
        htmlOutput(outputId = "gene_results_text"),

        htmlOutput("FC_data_title_cancer"),
        uiOutput(outputId = 'FC_data_cancer'),
        htmlOutput("FC_data_text_cancer"),

        uiOutput('FC_data_pathways_cancer'),

        htmlOutput("heatmap_cancer_cor_title"),
        uiOutput(outputId = "heatmap_cancer_cor"),
        htmlOutput("heatmap_cancer_cor_text"),

        ###
        htmlOutput("FC_data_title_neuro"),
        uiOutput(outputId = 'FC_data_neuro'),
        htmlOutput("FC_data_text_neuro"),

        uiOutput('FC_data_pathways_neuro'),

        htmlOutput("heatmap_neuro_cor_title"),
        uiOutput(outputId = "heatmap_neuro_cor"),
        htmlOutput("heatmap_neuro_cor_text"),
        ##

        #plotOutput(outputId = "heatmap_cancer_cor", height = "1000px"),

        #htmlOutput("pathways_meta_title"),
        #plotlyOutput(outputId = "pathways_meta"),

        htmlOutput(outputId = "gene_gene_title"),
        plotlyOutput(outputId = "gene_gene",height = "1000px",width = "1500px"),
        htmlOutput(outputId = "gene_gene_text"),

        htmlOutput(outputId = "gene_sample_title"),
        plotlyOutput(outputId = "gene_sample",height = "1000px",width = "1500px"),
        htmlOutput(outputId = "gene_sample_text"),

        htmlOutput("sign_pheno_exp_title"),
        dataTableOutput(outputId = "sign_pheno_exp"), #comparison with expression signature in database
        htmlOutput("sign_pheno_exp_text"),

        htmlOutput("sign_pathways_title"),
        dataTableOutput(outputId = "sign_pathways"),
        htmlOutput("sign_pathways_text"),

        htmlOutput("string_db_enr_title"),
        dataTableOutput(outputId = "string_db_enr"), #Enrichment analysis with stringdb
        htmlOutput(outputId = "string_db_enr_text"),
        uiOutput(outputId = "stringdb_box"),

        useShinyjs(),##seems to be no option for expression coloring in the embedded version
        extendShinyjs(functions = "loadStringData",text = "shinyjs.loadStringData = function(params) {
      var defaultParams = {
        genes : null,
        spec: '9606'
      };
      params = shinyjs.getParams(params,defaultParams);
      getSTRING('https://string-db.org', {'species': params.spec,'identifiers': params.genes,'network_flavor':'confidence'})
    };"),
        tags$head(tags$script(src = "http://string-db.org/javascript/combined_embedded_network_v2.0.2.js")),
        htmlOutput(outputId = "stringdb_title"),
        #stringdb_interactive
        tags$div(id = "stringEmbedded"),
        htmlOutput(outputId = "stringdb_text"),
        tags$head(tags$style(HTML("h3 {color: grey;font-size: 20px;font-family: Arial}
                              .shiny-notification {
                                       position:fixed;
                                       top: calc(50%);
                                       left: calc(50%);}

                              ")))

      )))
)


shinyApp(ui, server)


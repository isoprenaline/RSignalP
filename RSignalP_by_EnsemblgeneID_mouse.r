#' SignalP
#'
#' A function to Predicts the presence of signal peptides by SignalP-5.0 Server
#' @param GeneID, should be a file address str
#' @param b, should be a file address str
#' @param c, should be a file address str
#' @param d, should be a file address str
#' @param e, should be a file address str
#' @importFrom biomaRt useMart getBM exportFASTA
#' @importFrom reshape2 melt
#' @importFrom RSelenium remoteDriver
#' @export
#' @examples RsignalP(Gene="~/GeneID.csv");RsignalP(Gene="~/GeneID.csv",b,c,d,e)

#GeneID <- Ensembl gene IDs.csv, should be a csv file address str


RSignalP <- function(GeneID,b="~/RSignalP/AA Seq.fa",c="~/RSignalP/Results Of SignalP.txt",d="~/RSignalP/SP ProteinID.csv",e="~/RSignalP/SP GeneID.csv"){

  library(biomaRt)
  require(reshape2)
  library(RSelenium)

  dir.create("~/RSignalP")
  # input1
  EnsemblGeneID <- read.csv(file=GeneID, header=F, stringsAsFactors = FALSE)
  # define biomart object
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  # query biomart
  results <- apply(EnsemblGeneID, 2, function(x){
    getBM(attributes = "ensembl_peptide_id",
          filters = "ensembl_gene_id", values = x,
          mart = mart)})
  results <- melt(results)
  results <- results[,1]
  protein <- getSequence(id= results,
                         type="ensembl_peptide_id",
                         seqType="peptide",
                         mart= mart)
  #output1
  exportFASTA(protein, file = b)

  remDr <- remoteDriver(
    remoteServerAddr = "localhost",
    port = 4444,
    browserName = "chrome")
  remDr$open()
  remDr$navigate("http://www.cbs.dtu.dk/services/SignalP/index.php/")
  webElem <- remDr$findElement(using ='xpath','/html/body/div/div[2]/div[2]/div[2]/form/div[1]/div/div[2]/input[2]')
  webElem$clickElement()

  webElem <- remDr$findElement(using = 'xpath','//*[@id="fasta"]')
  webElem$clickElement()

  webElem <- remDr$findElement(using = 'id','uploadfile')

  # input2
  path_upload <- normalizePath("~/RSignalP/AA Seq.fa")
  webElem$sendKeysToElement(list(path_upload))
  webElem <- remDr$findElement(using = 'xpath','/html/body/div/div[2]/div[2]/div[2]/form/div[2]/div/input')
  webElem$clickElement()

  webElem <-NULL
  while(is.null(webElem)){
    webElem <- tryCatch({remDr$findElement(using = 'xpath', value = '//*[@id="single-button"]')},
                        error = function(e){NULL})
    print("Waiting fot the Result...")
    Sys.sleep(30)
  }

  webElem <- remDr$findElement(using = 'xpath','//*[@id="single-button"]')
  webElem$clickElement()

  a.elem <- webElem$findChildElement(using = 'xpath', '//*[@id="bs-example-navbar-collapse-1"]/ul[2]/div/ul/li[2]/a')
  path_download <- a.elem$getElementAttribute("href")
  path_download1 <- as.data.frame(path_download)
  path_download1 <- as.vector(path_download1[1,])
  #output2
  download.file(url = path_download1, c)

  my_data <- read.delim("~/RSignalP/Results Of SignalP.txt", header = F)
  colname_sp <- c("ensembl_peptide_id", "Prediction","SP(Sec/SPI)",	"OTHER",	"CS Position")

  names(my_data) <- colname_sp
  results <- grep("OTHER", my_data$Prediction)
  data <- my_data[-results,]
  proteinID <- data$ensembl_peptide_id
  geneID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","ensembl_peptide_id", "description"),
                  filters = "ensembl_peptide_id", values = proteinID,
                  mart = mart)

  co_sp <- merge(data, geneID, by = "ensembl_peptide_id")

  #output3
  write.csv(co_sp,file=d,quote=F)

  geneID1 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                   filters = "ensembl_peptide_id", values = data$ensembl_peptide_id,
                   mart = mart)
  #output4
  write.csv(geneID1, file=e,quote=F)

  print("Done, output 4 files in ~/RSignalP")

}


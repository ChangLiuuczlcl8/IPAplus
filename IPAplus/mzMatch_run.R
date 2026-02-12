install.packages ("remotes")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("https://github.com/andzajan/mzmatch.R.git",
                        build_opts=c("--no-multiarch"), INSTALL_opts=c("--no-test-load"))

library(mzmatch.R)
mzmatch.init(version.1=FALSE, memorysize = 5*1024)
setwd("Dataset/MTBLS2207/input/mzmatch/") 
mzmatch.R.Setup("sample_setup.tsv", projectFolder = getwd())
xseto <- xcmsSet(sampleList$filenames, method='centWave', ppm=5, peakwidth=c(5, 20),
                 snthresh=4, prefilter=c(1,100000), integrate=1, mzdiff=0.01, noise=100000,
                 verbose.columns=TRUE, fitgauss=FALSE,nSlaves=3)
PeakML.xcms.write.SingleMeasurement(xset=xseto, filename=sampleList$outputfilenames,
                                    ionisation="positive", ppm=5, addscans=0,
                                    ApodisationFilter=TRUE, nSlaves=3) 
mzmatch.ipeak.Combine(sampleList=sampleList, v=T, rtwindow=30, ppm=5, combination="set",
                      nSlaves=3)

dir.create(file.path("combined", "csv"), showWarnings = FALSE)
files <- list.files(path="combined", pattern='*.peakml', recursive=FALSE)
for (file in files) {
  mzmatch.ipeak.sort.RelatedPeaks(i=paste("combined/",file,sep=""),o=paste("combined/",file,sep=""),ppm=3,rtwindow=0.1)
  allpeaks <- PeakML.Read(paste("combined/",file,sep=""), ionisation="detect", Rawpath=NULL)
  allpeaksmatrix <- allpeaks[["peakDataMtx"]]
  write.csv(allpeaksmatrix, paste("combined/csv/",substring(file,1, nchar(file)-11),"csv",sep=""), row.names=FALSE)
}

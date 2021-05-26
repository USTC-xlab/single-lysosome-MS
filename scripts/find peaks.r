library(xcms)
library(MassSpecWavelet)

mzdata_path <- "./dir"
mzdata_files <- list.files(mzdata_path, recursive = TRUE, full.names = TRUE)

grp <- rep("sample", length(mzdata_files))
pd <- data.frame(filename = basename(mzdata_files), sample_group = grp)

## Load the data.
ham_raw <- readMSData(files = mzdata_files,
                      pdata = new("NAnnotatedDataFrame", pd),
                      mode = "onDisk")

## Define the parameters for the peak detection
msw <- MSWParam(scales = scales, nearbyPeak = TRUE, 
                winSize.noise = winSize.noise,SNR.method = "data.mean", snthresh = snthresh, ampTh = ampTh)

ham_prep <- findChromPeaks(ham_raw, param = msw)

head(chromPeaks(ham_prep))

mzc_prm <- MzClustParam(sampleGroups = ham_prep$sample_group, ppm = ppm, minFraction = minFraction)
ham_prep_0 <- groupChromPeaks(ham_prep, param = mzc_prm)

data <- cbind(featureDefinitions(ham_prep_0), featureValues(ham_prep_0, value = "into"))
data <- subset(data, select = -peakidx)
write.csv(data, "res_raw_0.csv")

#' buildSpectral Object
#'
#' @param data
#' @param meta
#' @param isAbsorbance
#' @param smoothingAmount
#'
#' @return
#' @export
#'
#' @examples

buildSpectraObject <- function(data, meta, isAbsorbance, smoothingAmount=NULL) {
  #' Takes raw data from text files and builds a MySpectra object
  #' Will convert absortion into reflectance if isAbsorbance = TRUE
  #' Will smooth the spectra using the Savitzky-Golay filter wet data if smoothingAmount is not NULL

  Spectra <- data[,grepl('X',colnames(data))]

  if (isAbsorbance) {
    Spectra <- 1 - Spectra
  }

  Bands <- gsub('X','',colnames(Spectra))

  IDs <- rownames(data)

  Spectra <- as.matrix(Spectra)

  # Smooth the spectra using the Savitzky-Golay filter wet data
  if (!is.null(smoothingAmount)) {
    if (!smoothingAmount %in% c(3, 5, 11)) {
      stop('smoothingAmount must be NULL, 3, 5, OR 11')
    }
    Spectra <- do.call(rbind, lapply(1:nrow(Spectra), function(spectrum){
      sgolayfilt(as.numeric(Spectra[spectrum,]),p=2,n=smoothingAmount,m = 0)
    }))
  }

  dimnames(Spectra) <- NULL

  MyspectraObject <- new('eSpectra')

  MyspectraObject@Meta <- meta
  MyspectraObject@Instrument <- "FieldSpec"
  MyspectraObject@Spectra <- Spectra
  MyspectraObject@Bands <- Bands
  MyspectraObject@Units <- 'nm'
  MyspectraObject@RowsAreSpectra <- TRUE
  MyspectraObject@Type <- 'Reflectance'
  MyspectraObject@ID <- IDs

  MyspectraObject
}

#' SNV
#'Do a standard normal variate transformation for baseline correction
#' @param spectra
#'
#' @return
#' @export
#'
#' @examples
SNV <- function(spectra){

  do.call(rbind,lapply(1:nrow(spectra),function(i){
    spectrum <- as.numeric(spectra[i,])
    meanSpec <- mean(spectrum)
    sdSpec <- sd(spectrum)
    finalSpec <- (spectrum-meanSpec)/sdSpec
  }))
}

#' plotAndCheck
#'plots spectra so you can check it is correct
#' @param spectrum
#' @param title
#'
#' @return
#' @export
#'
#' @examples
plotAndCheck <- function (spectrum, title) {
  plot(
    as.numeric(spectrum@Bands),
    spectrum@Spectra[1,],
    type='l',
    main=title,
    xlab='nm',
    ylab='Reflectance'
  )
}

#' trim_spectra
#'
#' @param Spectra
#' @param k
#'
#' @return
#' @export
#'
#' @examples
trim_spectra <- function(Spectra,k){
  do.call(rbind,lapply(1:nrow(Spectra),function(spec){
    as.numeric(Spectra[spec,seq(from = 1,to = ncol(Spectra),by = k)])}))
}

#' Set gap size

SMOOTHING_AMOUNT <- 11


# Set spectra range (wavelength in nm)
WAVELENGTH_RANGE <- c(400, 2450)

WAVELENGTH_RESOLUTION <- 5

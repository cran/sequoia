#' @title Sequoia report
#'
#' @description Create a summary report of the genotype data and reconstructed
#'   pedigree, utilising a range of functions in the \code{sequoia} package.
#'   For many of the plots and tables, brief help on interpretation is provided.
#'
#' @details  The R markdown (.Rmd) file for the report is included with the
#'   package, and may freely be modified and re-used.
#'
#' @param SeqList list with output from \code{\link{sequoia}}. If not provided,
#'   all other input is ignored and an example report using the griffin data is
#'   generated.
#' @param GenoM  genotype matrix that was used to run \code{\link{sequoia}}.
#' @param MaybeRel list with output from \code{\link{GetMaybeRel}}.
#' @param Conf  list with output from \code{\link{EstConf}}.
#' @param comments  optional comments to print at the top of the report.
#' @param printcode include R code in the report ('yes', default) or not ('no').
#'   The R code includes some brief notes and references to relevant functions.
#' @param output_file  filename (incl. path) for the report. Extension (.html,
#'   .pdf) not needed.
#' @param output_format  format for the report, e.g. \code{'html_document'},
#'   \code{rmarkdown::pdf_document(toc=TRUE, extra_dependencies = 'booktabs')},
#'   or \code{bookdown::gitbook(split_by='none')}. See
#'   \code{\link[rmarkdown]{render}} for details and options. Note that
#'   \code{'bookdown::gitbook'}, and possibly others, appear to only work after
#'   copying the .rmd file to the output folder.
#' @param copy_rmd copy the .rmd file to the output folder? (i.e. same folder as
#'   \code{output_file}).
#'
#' @return A file with the report in the specified output format.
#'
#' @examples
#' \dontrun{
#' sequoia_report(comments = 'This is a test',
#'                file = file.path(getwd(), 'test_griffins'))
#' }
#'
#' @export

sequoia_report <- function(SeqList = NA,
                           GenoM = NA,
                           MaybeRel = NA,
                           Conf = NA,
                           comments = '',
                           printcode = 'yes',
                           output_file = '',
                           output_format = 'html_document',
                           copy_rmd = FALSE) {

  if (all(is.na(SeqList))) {
    ANS <- readline(prompt = paste("No SeqList provided; an example report will be",
                                   "generated based on the griffin data. Continue Y/N? "))
    if (!substr(ANS, 1, 1) %in% c("Y", "y")) stop(call.=FALSE)
  }
  if (output_file == '')  stop('Please specify output_file')

  rmd_file <- system.file('rmd', 'sequoia_report.Rmd', package='sequoia')
  if (copy_rmd) {
    out_dir <- dirname(output_file)
    file.copy(from = rmd_file, to = file.path(out_dir, 'sequoia_report.Rmd'))
    rmd_file <- file.path(out_dir, 'sequoia_report.Rmd')
  }

  rmarkdown::render(input = rmd_file,
                    output_format = output_format,
                    output_file = output_file,
                    params = list(output_sequoia = SeqList,
                                  genotypes = GenoM,
                                  output_GetMaybeRel = MaybeRel,
                                  output_EstConf = Conf,
                                  comments = comments,
                                  printcode = printcode))
}

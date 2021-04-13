##' \pkg{mpmm}
##'
##' fit movement persistence mixed-effects models to animal location data
##'
##' @name mpmm-package
##' @aliases mpmm-package
##' @docType package
##' @author Ian Jonsen
##'
##' @seealso mpmm
##' @references Jonsen ID, McMahon CR, Patterson TA, et al. (2019) Movement responses to environment: fast inference of variation among southern elephant seals with a mixed effects model. Ecology. 100(1):e02566 https://doi.org/10.1002/ecy.2566
##'
##' @keywords mpmm
##' @importFrom lme4 nobars findbars subbars mkReTrms
##' @importFrom glmmTMB getReStruc splitForm
##' @importFrom Matrix t
##' @importFrom dplyr %>% arrange count mutate left_join group_by select
##' @importFrom dplyr bind_cols
##' @importFrom TMB MakeADFun sdreport newtonOption
##' @importFrom tibble tibble as_tibble
##' @importFrom parallel detectCores
##' @importFrom ggplot2 ggplot geom_line aes xlab ylab theme_bw theme ylim xlim
##' @importFrom ggplot2 element_text facet_wrap element_blank
##' @importFrom stats plogis
##' @importFrom tidyr pivot_longer everything
##' @importFrom wesanderson wes_palette
##' @importFrom stats pnorm AIC BIC logLik var getCall pchisq anova
##' @importFrom methods is
##' @importFrom utils globalVariables
NULL

##' @name ellie.ice
##' @docType data
##' @title foieGras-filtered Southern elephant seal Argos satellite data with
##' environmental covariates (11 individuals)
##' @format .RData
##' @keywords data
##' @description Example elephant seal Argos tracking data with environmental
##' covariates. Data were sourced from the Integrated Marine Observing System
##' (IMOS) - IMOS is supported by the Australian Government through the National
##' Collaborative Research Infrastructure Strategy and the Super Science
##' Initiative.
NULL

##' @name ellie.ice.short
##' @docType data
##' @title foieGras-filtered Southern elephant seal Argos satellite data with
##' environmental covariates (4 individuals)
##' @format .RData
##' @keywords data
##' @description Example elephant seal Argos tracking data with environmental
##' covariates. Data were sourced from the Integrated Marine Observing System
##' (IMOS) - IMOS is supported by the Australian Government through the National
##' Collaborative Research Infrastructure Strategy and the Super Science
##' Initiative.
NULL

## stop R CMD check generating NOTES about global variables
globalVariables(c(".", "id", "tid", "model.matrix", "g", "lg.se", "predictor"))

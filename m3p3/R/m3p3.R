#' simulates a single dose-escalation studies with modified 3+3 design
#'
#' @param truerate  a vector of true rates for DLT probabilities at different dose levels
#' @param targetrate the set target rate of DLT probability to determine the MTD
#' @param seed if not empty, the seed to use for random generation
#'
#' @return the output from \code{\link{m.sim3p3}}
#' @export
#' @importFrom UBCRM CreData
#' @importFrom UBCRM updata
#' @importFrom UBCRM troisPtrois

#' @importFrom stats rbinom
#' @importFrom Iso pava

#'
#' @examples m.sim3p3(c(0.25,0.35,0.4,0.45,0.5,0.6), 0.25)
m.sim3p3<-function (truerate, targetrate, seed = NULL)
{
  data <- CreData(length(truerate))
  nextdose <- 1
  while (nextdose %in% data$dose) {
    lastdose <- nextdose
    ndlt <- sum(rbinom(3, 1, truerate[lastdose]))
    data <- updata(data, lastdose, 3, ndlt)
    nextdose <- troisPtrois(data, lastdose)$nextdose
  }
  mtd <- troisPtrois(data, lastdose)$mtd
  if ( is.na(mtd)  ){
    list(data = data, mtd = mtd, lastdose = lastdose)
  } else{
    maxs=length(truerate)*6
    ren=maxs-sum(data[,2])
    ret=rbinom(1,ren,truerate[mtd])
    data[mtd,2]=data[mtd,2]+ren
    data[mtd,3]=data[mtd,3]+ret

    p=data$ndlt[data$npt>0]/data$npt[data$npt>0]
    p.iso=pava(p,data$npt[data$npt>0])
    phi=targetrate
    temp=abs(p.iso-phi)
    new.mtd=which(temp==min(temp)  )[1]
    list(data = data, mtd = new.mtd, lastdose = lastdose)
  }
}







#' simulates n dose-escalation studies with the modified 3+3 design
#'
#' @param truerate a vector of true rates for DLT probabilities at different dose levels
#' @param targetrate the set target rate of DLT probability to determine the MTD
#' @param n number of studies to simulate
#' @param r integer, number of digits for percentages in output
#' @param seed if not empty, the seed to use for random generation
#'
#' @return the output from \code{\link{m.ssim3p3}}
#' @export
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom UBCRM CreData

#'
#' @examples m.ssim3p3(c(0.25,0.35,0.4,0.45,0.5,0.6), 0.25, 10000)
m.ssim3p3 <- function (truerate, targetrate, n, r = 2, seed = NULL)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  lp <- length(truerate)
  fmtd <- rep(0, lp)
  lastdose <- 0
  mnpt <- mndlt <- matrix(rep(0, lp * n), lp, n)
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
  for (i in 1:n) {
    sim <- m.sim3p3(truerate, targetrate)
    mnpt[, i] <- sim$data$npt
    mndlt[, i] <- sim$data$ndlt
    if (sim$mtd %in% sim$data$dose) {
      lastdose <- lastdose + sim$lastdose
      fmtd[sim$mtd] <- fmtd[sim$mtd] + 1
      setTxtProgressBar(pb, i/n)
    }
  }
  close(pb)
  data <- CreData(lp)
  data$npt <- round(apply(mnpt, 1, mean), r)
  data$ndlt <- round(apply(mndlt, 1, mean), r)
  pdlt <- round(apply(mndlt, 1, sum)/apply(mnpt, 1, sum), r)
  exp = apply(mnpt, 1, sum) * 100/sum(mnpt)
  overshoot <- c(100, rep(0, lp))
  for (i in 1:lp) {
    overshoot[i + 1] <- overshoot[i] - exp[i]
  }
  list(data = cbind(data, pdlt, recommendation = fmtd * 100/n,
                    experimentation = round(exp, r), overshoot = round(overshoot[-1],
                                                                       r)), norecommendation = round(100 - sum(fmtd * 100/n),
                                                                                                     r), mean.npt = sum(data$npt), mean.ndlt = sum(data$ndlt),
       mean.lastdose = round(lastdose/n, r))
}




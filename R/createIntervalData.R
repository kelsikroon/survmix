#' Create interval censored data
#'
#' @param patient_info Data frame containing patient information such as Idwoman, Study, Danasc, Age, Darecl, Arm, DaCens, StCens
#' @param hist Data frame containing resutls from histological tests for patients. Must have column names Study, Idwoman, VisitHisto, Idhisto, Datehisto, Typehisto, Istorandom, Classhisto.
#' @param cyt Data frame containing results from cytologyical tests for patients. Must have column names Study, Idwoman, VisitCyto, Idtest, HPVpres, cyto_idHPV, Datecyto, Type, CytMed, Classcyto, Racctest, FUP
#' @param hpv  Data frame containing results from HPV tests for patients
#'
#' @returns A data frame containing the left and right intervals of CIN2+ detection for women who were HPV+ at baseline.
#' @export

create.interval.data <- function(patient_info, hist, cyt, hpv){
  rownames(patient_info) <- patient_info$Idwoman

  hist$Datehisto <- as.Date(hist$Datehisto, "%d-%m-%Y")
  cyt$Datecyto <- as.Date(cyt$Datecyto, "%d-%m-%Y")

  # keep only women that are hrHPV+ at their first visit (A)
  hpv_pos <- hpv[hpv$ResHPV1 ==1 & hpv$VisitHPV=='A',]

  # store baseline date for each woman
  baselines <- patient_info[patient_info$Idwoman %in% hpv_pos$Idwoman,]$Darecl

  hpv_pos_ids <- rep(NA, dim(hpv_pos)[1])
  for (i in 1:dim(hpv_pos)[1]){
    id <- hpv_pos$Idwoman[i]
    baseline <- baselines[i]
    # only keep woman ID if their first visit matches the date of recruitment
    if (hpv_pos$DateHPV[i]==baseline) hpv_pos_ids[i] <- id
  }
  hpv_pos_copy <- hpv_pos_ids

  baseline_pos <- hpv_pos_ids[!is.na(hpv_pos_ids)]

  cyt_HPVpos <- cyt[cyt$Idwoman %in% baseline_pos,]
  baseline_cyt <- cyt_HPVpos[cyt_HPVpos$VisitCyto=='A',]$Classcyto
  names(baseline_cyt) <- cyt_HPVpos[cyt_HPVpos$VisitCyto=='A',]$Idwoman

  hist_HPVpos <- hist[hist$Idwoman %in% names(baseline_cyt),]

  hist_HPVpos$cin2plus <- ifelse(hist_HPVpos$Classhisto %in% c(7, 8, 10, 12), 1, 0)

  ids <- unique(hist_HPVpos$Idwoman)
  right <- rep(NA, length(ids))
  left <- rep(NA, length(ids))

  for (i in 1:length(ids)){
    histo <- hist_HPVpos[hist_HPVpos$Idwoman == ids[i],c(2, 5, 8, 9)]
    cyto <- cyt[cyt$Idwoman==ids[i],c(2, 7, 10)]
    baseline <-  as.Date(patient_info[patient_info$Idwoman ==ids[i],]$Darecl, "%d-%m-%Y")

    # if there are no histology results, then the left interval is the last Pap 1 cytology result
    # and the right interval is infinity
    if (nrow(histo)==0){
      cyt_index <- length(cyto$Classcyto)-match(1,rev(cyto$Classcyto))+1 # index of last Pap1 cytology result
      cytodate <- as.Date(ifelse(is.na(cyt_index), baseline, cyto$Datecyto[cyt_index]), origin='1970-01-01') # date of last Pap1
      left[i] <- as.numeric(difftime(cytodate, baseline)) # time in days since baseline
      right[i] <- Inf # right is infinity because they did not develop CIN2+ in follow-up
      next
    }

    # if all histology results are 0 then CIN2+ did not develop during follow-up
    # so the left interval is the max date between last CIN1 result or last Pap 1
    # result, and the right interval is infinity
    if (sum(histo$cin2plus) ==0){
      histodate <- histo$Datehisto[length(histo$Datehisto)]
      cytodate <- cyto$Datecyto[length(cyto$Datecyto)-match(1,rev(cyto$Classcyto))+1]
      if (is.na(cytodate)) {
        left[i] <- difftime(histodate, baseline)
      }else{
        left[i] <- difftime(as.Date(max(histodate, cytodate)), baseline)
      }
      right[i] <- Inf
      next
    }

    # we find the index of the first time the histology result is CIN2+
    # and then save the save of the first CIN2+ result
    right_index <- match(1, histo$cin2plus)
    first_cin <- histo$Datehisto[right_index]

    # we check if there was a cytology result >Pap 1 within 3 months before the first CIN2+ result
    # as this would be the index smear that led to histology referral
    if (as.numeric(difftime(first_cin, max(cyto[cyto$Datecyto<first_cin & cyto$Classcyto>1 & cyto$Classcyto!=99,]$Datecyto))) < 90){
      r_int <-  max(cyto[cyto$Datecyto<first_cin &  cyto$Classcyto>1 & cyto$Classcyto!=99,]$Datecyto) # index smear date
    }else{
      r_int <- first_cin # histology date
    }
    right[i] <- as.numeric(difftime(r_int, baseline, units='days'))

    # left interval:
    # filter cytology date
    cyto_temp <- cyto[cyto$Datecyto < first_cin,]
    cyto_temp <- cyto_temp[cyto_temp$Classcyto ==1,]
    histo_temp <- histo[histo$Datehisto < histo$Datehisto[right_index],]

    left_index <- ifelse(nrow(cyto_temp)==0, 0, max(cyto_temp$Datecyto)) # find the last date lower than Pap1
    cytodate <- ifelse(left_index==0 |(left_index==-Inf), 0, as.numeric(difftime(as.Date(left_index, origin='1970-01-01'), baseline))) # get the left interval in days since baseline


    left_index <- ifelse(nrow(histo_temp)==0, 0, # if there are no previous histology results
                         ifelse(histo_temp$cin2plus[1]==1, 0, # if the first histology result is CIN2+ then left is zero
                                Position(function(x) 1-x, histo_temp$cin2plus, right=T))) # otherwise find the last time <CIN2+

    histodate <- ifelse(left_index==0, 0, as.numeric(difftime(histo_temp$Datehisto[left_index], baseline)))

    left[i] <- max(histodate, cytodate)
  }
  # - check which women are censored and check the date of censoring:
  #     if the date is after a right interval then ignore the censoring
  #     if the date is before a right interval, then censor becomes the left interval and right=Inf (but check with Hans)
  cens_dates <- unname(sapply(patient_info[patient_info$Idwoman %in% ids,]$DaCens, function(x) as.Date(paste("01-",x,sep=""), "%d-%m-%Y")))
  cens_dates <- as.Date(cens_dates, origin='1970-01-01')
  censoring <- patient_info[patient_info$Idwoman %in%ids,]$StCens
  for (i in 1:length(ids)){
    baseline <- as.Date(patient_info[patient_info$Idwoman ==ids[i],]$Darecl, "%d-%m-%Y")
    censdate <- difftime(cens_dates[i], baseline)
    if (is.na(censdate)){
      remove_id <- ids[i]
      next
    }
    if (censoring[i] %in% c(1, 2)){
      if(right[i]==Inf & left[i]< censdate){
        left[i] <- censdate
        next
      }
    }
  }

  {
    new_data <- data.frame(left=left/365, right= right/365)
    new_data$right[new_data$right < 0.25] <- 0
    new_data$left[new_data$right ==0] <- 0
    new_data$z <- ifelse(new_data$right==0, 1,
                         ifelse(new_data$left==0, NA, 0))
    rownames(new_data) <- ids
    new_data$age <- patient_info[rownames(new_data),]$Age
    new_data$age <- ifelse(new_data$age >39, 1, 0)
    new_data$hpv16 <- hpv[hpv$Idwoman%in% ids & hpv$VisitHPV=='A',]$HPV16_2
    new_data$hpv16[new_data$hpv16 ==9] <- 0
    new_data$hpv18 <- hpv[hpv$Idwoman%in% ids & hpv$VisitHPV=='A',]$HPV18_2
    new_data$hpv18[new_data$hpv18 ==9] <- 0
    new_data$cyt_bmd <- ifelse(baseline_cyt[rownames(new_data)] %in% c("2","2a", "3"), 1, 0)
    new_data$cyt_severe <- ifelse(baseline_cyt[rownames(new_data)] >=4, 1, 0)
    new_data <- new_data[!baseline_cyt[rownames(new_data)] %in% c("99", "", "0"),]
    new_data <- new_data[rownames(new_data)!=remove_id,]

  }
  return(new_data)
}

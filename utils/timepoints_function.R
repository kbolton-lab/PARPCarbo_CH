require(tidyverse)

individual_timepoints <- function(timepoint) {
    # If the variant was called at this time point, then we should use the VAF that the caller spits out
    if (as.logical(timepoint["was_called"])) {
        return("called")
        # If it wasn't called, then we have to check if the pileup found it for this time point
    } else if (as.logical(timepoint["has_pileup"])) {
        # If the variant passes the fisher's exact test, then we can rescue it
        if (as.logical(timepoint["fisher_passed"])) {
            return("rescued")
            # If the variant does not pass. The pileup still found it, but we cannot assume the pileup is right, we have to assume the limit of detection
        } else {
            return("at_limit")
        }
        # If it wasn't called and pileup did not detect it, then this variant for this time point is missing
    } else {
        return("missing")
    }
}

set_vaf <- function(df) {
    # df <- df %>% tidyr::separate(information, into = c("lenient_pass", "average_AF", "was_called", "has_pileup", "fisher_passed", "pileup_vaf"), sep = "\\|")
    vaf <- ifelse(df["lenient_pass"], df["average_AF"], NA)
    vaf_ignore_noise <- ifelse(df["lenient_pass"], df["average_AF"], NA)
    timepoint <- ifelse(is.na(vaf), individual_timepoints(df), "lenient_pass")
    vaf <- ifelse(timepoint == "called", df["average_AF"], vaf)
    vaf_ignore_noise <- ifelse(timepoint == "called", df["average_AF"], vaf_ignore_noise)
    vaf_ignore_noise <- ifelse(timepoint == "at_limit", df["PILEUP_VAF"], vaf_ignore_noise) # If the variant was not called, but the pileup found it, we should use the pileup VAF
    vaf_ignore_noise <- ifelse(!is.na(df["PILEUP_TOTAL_DEPTH"]) & timepoint == "missing", 0.5 / as.numeric(df["PILEUP_TOTAL_DEPTH"]), vaf_ignore_noise)
    return(c(vaf, vaf_ignore_noise, timepoint))
}

time_points <- function(variant, ignore_noise = TRUE, rescue = TRUE) {
    variant <- variant %>%
        mutate(
            lenient_pass = case_when(
                as.logical(pon_FP_pass_XGB) &
                    as.logical(long100_indel_pass_XGB) & as.logical(long_indel_pass_XGB) & as.logical(di_tri_vard_pass_XGB) & as.logical(bcbio_pass_XGB) & as.logical(zscore_pass_XGB) &
                    as.logical(PASS_BY_1) &
                    FP_Filter_DETP20_XGB == 0 & FP_Filter_MMQS100_XGB == 0 & FP_Filter_MMQSD50_XGB == 0 & FP_Filter_NRC_XGB == 0 & FP_Filter_PB10_XGB == 0 & FP_Filter_RLD25_XGB == 0 ~ TRUE,
                TRUE ~ FALSE
            )
        ) %>%
        mutate(lenient_pass = ifelse(FP_Filter_MVC4_XGB == 1 & PILEUP_ALT_DEPTH > 5, FALSE, lenient_pass)) %>%
        mutate(lenient_pass = ifelse(FP_Filter_SB1_XGB == 1 & PILEUP_ALT_DEPTH > 10, FALSE, lenient_pass)) %>%
        mutate(lenient_pass = ifelse(is.na(lenient_pass), FALSE, lenient_pass))

    # Check to see if the new definition creates more "PASS" all time points
    # If the variant was CALLED + PASS in ALL time points
    if (sum(as.integer(as.logical(variant$putative_driver))) == length(variant$putative_driver)) {
        all_time_points <- "true"
        # If the variant was CALLED + Lenient PASS in ALL time points
        # This will handle variants that only failed because of LowVAF or SB1/MVC4
    } else if (sum(as.integer(variant$lenient_pass)) == length(variant$lenient_pass)) {
        all_time_points <- "true"
    } else {
        all_time_points <- "false"
    }
    variant <- variant %>% select(lenient_pass, average_AF, was_called, has_pileup, fisher_passed, PILEUP_VAF, PILEUP_TOTAL_DEPTH, whereincycle)
    
    variant <- data.frame(whereincycle = 1:6) %>%
        left_join(variant, by = "whereincycle") %>%
        mutate(
            lenient_pass = ifelse(is.na(lenient_pass), FALSE, lenient_pass),
            was_called = ifelse(is.na(was_called), FALSE, was_called),
            has_pileup = ifelse(is.na(has_pileup), FALSE, has_pileup),
            fisher_passed = ifelse(is.na(fisher_passed), FALSE, fisher_passed)
        )

    res <- as.data.frame(t(apply(variant, 1, set_vaf)))
    names(res) <- c("vaf", "vaf_ignore_noise", "timepoint")
    variant <- cbind(variant, res)

    if ("rescued" %in% variant$timepoint) {
        # If one of the variants was rescued, we should use the pileup IF it has pileup
        variant <- variant %>% 
            mutate(
                vaf_rescue = ifelse(has_pileup, PILEUP_VAF, vaf),
                vaf_ignore_noise_rescue = ifelse(has_pileup, PILEUP_VAF, vaf_ignore_noise)
            )
        rescued <- "true"
    } else {
        variant <- variant %>% 
            mutate(
                vaf_rescue = vaf,
                vaf_ignore_noise_rescue = vaf_ignore_noise
            )
        rescued <- "false"
    }

    missing <- ifelse(any(variant$timepoint[!is.na(variant$PILEUP_VAF)] == "missing"), "true", "false")
    at_limit <- ifelse(any(variant$timepoint[!is.na(variant$PILEUP_VAF)] == "at_limit"), "true", "false")
    called <- ifelse(any(variant$timepoint[!is.na(variant$PILEUP_VAF)] == "called"), "true", "false")

    res <- as.data.frame(t(variant %>% arrange(whereincycle) %>% select(vaf, vaf_ignore_noise, vaf_rescue, vaf_ignore_noise_rescue)))
    return_string <- paste(res[1, ], collapse = " ")
    if (ignore_noise) {
        return_string <- paste(return_string, paste(res[2, ], collapse = " "), sep = " ")
    }
    if (rescue) {
        return_string <- paste(return_string, paste(res[3, ], collapse = " "), sep = " ")
    }
    res <- paste(all_time_points, return_string, rescued, missing, called, at_limit, sep = " ")

    return(res)
}

# Functions to Process Timepoints AFTER the above
timepoints_to_columns <- function(df) {
    df %>%
      separate(result, c("all_time_points", "prepreVAF_original", "preVAF_original", "duringVAF_original", "postVAF_original", "postpostVAF_original", "extra_tp_original",
                     "prepreVAF_ignore_noise", "preVAF_ignore_noise", "duringVAF_ignore_noise", "postVAF_ignore_noise", "postpostVAF_ignore_noise", "extra_tp_ignore_noise", 
                     "prepreVAF_rescue", "preVAF_rescue", "duringVAF_rescue", "postVAF_rescue", "postpostVAF_rescue", "extra_tp_rescue",
                     "rescued", "missing", "called", "at_limit"), sep = " ", extra = "merge", fill = "right") %>%
  select(-extra_tp_original, -extra_tp_ignore_noise, -extra_tp_rescue) %>%
  mutate_if(is.character, ~na_if(., "NA"))
}

timepoints_set_vaf <- function(timepoints_df, limit_of_detection = FALSE) {
    # 1. Set the VAFs to the original caller VAFs
    timepoints_df <- timepoints_df %>%
      mutate(
        prepreVAF_original = as.numeric(prepreVAF_original),
        preVAF_original = as.numeric(preVAF_original),
        duringVAF_original = as.numeric(duringVAF_original),
        postVAF_original = as.numeric(postVAF_original),
        postpostVAF_original = as.numeric(postpostVAF_original),
        prepreVAF = as.numeric(prepreVAF_original),
        preVAF = as.numeric(preVAF_original),
        duringVAF = as.numeric(duringVAF_original),
        postVAF = as.numeric(postVAF_original),
        postpostVAF = as.numeric(postpostVAF_original),
      )
    # 2. If the variant can be rescued, we will use the rescued VAFs
    # 2a. We need to change the VAFs to the rescued VAFs if the variant was rescued (apples to oranges)
    timepoints_df <- timepoints_df %>%
        mutate(
            prepreVAF = ifelse(rescued == "true", as.numeric(prepreVAF_rescue), prepreVAF),
            preVAF = ifelse(rescued == "true", as.numeric(preVAF_rescue), preVAF),
            duringVAF = ifelse(rescued == "true", as.numeric(duringVAF_rescue), duringVAF),
            postVAF = ifelse(rescued == "true", as.numeric(postVAF_rescue), postVAF),
            postpostVAF = ifelse(rescued == "true", as.numeric(postpostVAF_rescue), postpostVAF)
        )
    
    # 3. If the variant was NOT rescuable, we have an option to explore the limit of detection
    if (limit_of_detection) {
        timepoints_df <- timepoints_df %>%
            mutate(
                prepreVAF = ifelse(is.na(prepreVAF), as.numeric(prepreVAF_ignore_noise), prepreVAF),
                preVAF = ifelse(is.na(preVAF), as.numeric(preVAF_ignore_noise), preVAF),
                duringVAF = ifelse(is.na(duringVAF), as.numeric(duringVAF_ignore_noise), duringVAF),
                postVAF = ifelse(is.na(postVAF), as.numeric(postVAF_ignore_noise), postVAF),
                postpostVAF = ifelse(is.na(postpostVAF), as.numeric(postpostVAF_ignore_noise), postpostVAF)
            )
    }

    # 4. Finally, we set all NA values to 0
    timepoints_df <- timepoints_df %>%
        mutate(
            prepreVAF = ifelse(is.na(prepreVAF), 0, prepreVAF),
            preVAF = ifelse(is.na(preVAF), 0, preVAF),
            duringVAF = ifelse(is.na(duringVAF), 0, duringVAF),
            postVAF = ifelse(is.na(postVAF), 0, postVAF),
            postpostVAF = ifelse(is.na(postpostVAF), 0, postpostVAF)
        )

    return(timepoints_df)
}

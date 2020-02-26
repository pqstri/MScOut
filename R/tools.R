require("tidyverse")


stabilize <- function() {
}


EDSS_relapse <-
  function(relapse_data, edss_data,
           relapse_padding, puid, confirmed_relapse,
           relapse_date, edss_date) {

    # adjust working variables
    relapse_preJoint <- select(
      relapse_data, puid = puid,
      confirmed_relapse  = confirmed_relapse,
      relapse_date       = relapse_date) %>%

      # set tollerance window around relapse date
      mutate(relapse_pre  = relapse_date - relapse_padding,
             relapse_post = relapse_date + relapse_padding) %>%
      filter(!is.na(relapse_date)) %>%

      # group relapse date in one-row-per-patient
      group_by(puid) %>%
      summarise(relapse_pre  = list(relapse_pre),
                relapse_post = list(relapse_post),
                relapse_conf = list(confirmed_relapse)) %>%

      rowwise() %>%
      mutate(relapse_conf_pre  = list(relapse_pre[unlist(relapse_conf) == 1]),
             relapse_conf_post = list(relapse_post[unlist(relapse_conf) == 1]))

    # merge relapse and EDSS information
    edss_rel_data <- edss_data %>%
      rename(puid = puid, edss_date = edss_date) %>%
      left_join(relapse_preJoint) %>%
      rowwise() %>%

      # identify EDSS measurements within the relapse window
      mutate(within_relapse = any(
        unlist(relapse_pre) < edss_date &
          unlist(relapse_post) > edss_date)) %>%

      # identify EDSS measurements within the confirmed relapse window
      mutate(within_conf_relapse = any(
        unlist(relapse_conf_pre) < edss_date &
          unlist(relapse_conf_post) > edss_date))

    # return elaborated data set
    edss_rel_data
  }


cross_confirm <- function(p, w, m) {
  # p <- c(24, 97, 180, 330)
  # w <- 90
  # m <- 15
  if (all(is.na(p))) {return(NA)}
  if(length(p) == 0) {return(NA)}

  l <- length(p)

  step.1 <- data.frame(pr = p, length = l) %>%
    uncount(length) %>%
    mutate(pl = rep(p, l)) %>%
    filter(pl != pr) %>%
    mutate(pl_l = pl+w-m, pl_u = pl+w+m)

  if (nrow(step.1) == 0) {return(NA)}

  step.2 <- rowwise(step.1) %>%
    mutate(within = pr %in% pl_l:pl_u,
           # overlap = pr-pl
    ) %>%
    filter(within) %>%
    ungroup()

  if (nrow(step.2) == 0) {return(NA)}

  min(step.2$pl)
}


confirm_dispach <- function(p, w, m) {
  purr::map_dbl(p, function(x) cross_confirm(x, w, m))
}


EDSS_progressions <- function(edss_rel_data, conf_time_threshold, conf_time_padding) {

  step.1 <- edss_rel_data %>%

    # adjust name references
    rename(edss          = edss,
           edss_baseline = edss_baseline) %>%

    # remove unwanted visits
    # such as unscheduled or screening
    filter(edss_date >= 0) %>%

    # remove EDSS visits within relapses
    filter(!within_conf_relapse) %>% #!!! BUG

    # set baseline

    # compute delta from current visit EDSS and baseline EDSS
    mutate(delta_EDSS = edss - edss_baseline) %>%

    # for each patient
    group_by(puid) %>%

    # check if progression requirement is met at each visit
    mutate(edss_class      = cut(edss_baseline, edss_thresholds),
           min_progression = edss_steps[as.numeric(edss_class)],
           edss_change     = edss - edss_baseline,
           progression     = as.numeric(edss_change >= min_progression)) %>%


    mutate(prog_dates = list(edss_date[progression == 1]))

  step.1$prog_conf_time <-
    confirm_dispach(p = step.1$prog_dates,
                    w = conf_time_threshold,
                    m = conf_time_padding)

  step.1
}


summarizer <- function(EDSS_conf_prog) {

  no_events <- EDSS_conf_prog %>%
    group_by(puid) %>%
    arrange(edss_date) %>%

    filter(!any(!is.na(prog_conf_time))) %>%
    summarise(event = 0,
              time = max(edss_date))

  events <- EDSS_conf_prog %>%
    group_by(puid) %>%
    arrange(edss_date) %>%

    filter(!is.na(prog_conf_time)) %>%
    summarise(event = 1,
              time = first(prog_conf_time))

  bind_rows(events, no_events)

}

merge_back_to <- function(df_end, df_start, puid = "puid", method = first()) {
  df_start %>%
    rename("puid" = puid) %>%
    left_join(df_end, by = c("puid")) %>%
    group_by(puid) %>%
    summarise_all(method) %>%
    rename(puid = "puid") %>%
    mutate(event = ifelse(is.na(event), 0, event),
           time = ifelse(is.na(time), 1, time))
}



compute_progression <- function(edss_df, relapse_df,
                                puid, confirmed_relapse, relapse_date, edss_date,
                                relapse_padding, conf_time_threshold, conf_time_padding, merge_method = first) {

  suppressWarnings(suppressMessages(
  edss_df %>%
    EDSS_relapse(relapse_data      = relapse_df,
                 relapse_padding   = relapse_padding,
                 puid              = puid,
                 confirmed_relapse = confirmed_relapse,
                 relapse_date      = relapse_date,
                 edss_date         = edss_date) %>%
    EDSS_progressions(conf_time_threshold = conf_time_threshold,
                      conf_time_padding   = conf_time_padding) %>%
    summarizer() %>%
    merge_back_to(df_start = edss_df, puid = puid, method = merge_method)
  ))
}



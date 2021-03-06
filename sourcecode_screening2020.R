if (!("pacman" %in% rownames(installed.packages()))) {install.packages("pacman")}
pacman::p_load(pacman,
               tidyverse,
               readxl,
               lubridate,
               mgcv,
               gratia,
               furrr,
               tictoc,
               magrittr,
               beepr,
               RcppRoll,
               scales)

select <- dplyr::select
periods <- function(data,
                    variable,
                    threshold = 3,
                    filter_less_than = NA){
  variable_enquo <- enquo(variable)
  groups_vec <- data %>% group_vars() %>% syms
  out <- data %>%
    dplyr::arrange(!!variable_enquo) %>%
    mutate(period = cumsum(c(1, diff(!!variable_enquo) >= threshold+1))) %>%
    group_by(period, add=T) %>%
    mutate(n_year_period = n_distinct(!!variable_enquo))

  if(!is_empty(groups_vec)){out <- out %>% group_by(!!!groups_vec)}else{out <- out %>% ungroup}

  if(!is.na(filter_less_than)){out <- out %>%
    filter(n_year_period >= filter_less_than) %>%
    mutate(period = cumsum(c(1, diff(!!variable_enquo) >= threshold+1)))
  }

  return(out)
}

trim_tails <- function(range = c(-Inf, Inf)) {
  trans_new("trim_tails",
    transform = function(x) {force(range)
                             desired_breaks <- extended_breaks(n = 7)(x[x >= range[1] & x <= range[2]])
                             break_increment <- diff(desired_breaks)[1]
                             x[x < range[1]] <- range[1] - break_increment
                             x[x > range[2]] <- range[2] + break_increment
                             x},
    inverse = function(x) {x},
    breaks = function(x) {force(range)
                          extended_breaks(n = 7)(x)},
    format = function(x) {force(range)
                          x[1] <- paste("<", range[1])
                          x[length(x)] <- paste(">", range[2])
                          x}
  )
}
## from: https://stackoverflow.com/questions/44628130/ggplot2-dealing-with-extremes-values-by-setting-a-continuous-color-scale
# by: Brian W. Davis

model <- function(x, link = "identity", formula, knots, opt = "nlminb") {
  y <- gamm(
    data = x,
    formula = formula,
    family = gaussian(link = link),
    method = "REML", control = list(opt = opt),
    knots = knots,
    correlation = corCAR1(form = ~decimaldate)
  )$gam
  y$data <- x
  y$autocor <- TRUE
  return(y)
}

model_gam <- function(x, link = "identity", formula, knots) {
  y <- gam(
    data = x,
    formula = formula,
    family = gaussian(link = link),
    knots = knots,
    method = "REML"
  )
  y$autocor <- FALSE
  return(y)
}

model_gam_t <- function(x, link = "identity", formula, knots) {
  y <- gam(
    data = x,
    formula = formula,
    family = scat(link = link),
    knots = knots,
    method = "REML"
  )
  y$autocor <- FALSE
  return(y)
}

modeling <- function(x, link = "identity", autocor = autocor, tdist = tdist) {
    formula <- variable ~ s(decimaldate, k = round(nrow(x) / 2)) + s(decimal_during, bs = "cc", k = 13)
    knots <- list(decimal_during = seq(0, 1, length=13))
  x <- drop_na(x, variable)
  if(autocor == TRUE)
    {out <- try(model(x, link, formula, knots))
  if ("try-error" %in% class(out)) {
    out <- try(model(x, link, formula, knots, opt = "optim"))
  }
  if ("try-error" %in% class(out)) {
    out <- model_gam(x, link, formula, knots)
  }}else{if(tdist == T){out <- model_gam_t(x, link, formula, knots)}
    else{out <- model_gam(x, link, formula, knots)}}

  return(out)
}

inverse.link <- function(x, link) {
  switch(link,
         "identity" = x,
         "logit" = 1 / (1 + exp(-x)),
         "log" = exp(x)
  )
}

link.func <- function(x, link) {
  switch(link,
         "identity" = x,
         "logit" = log(x / (1 - x)),
         "log" = log(x)
  )
}

screeningmodeling <- function(.data,
                     datevar, # Date variable (in date format)
                     values, # value variable
                     link = "identity",
                     autocor = TRUE,
                     conf.type = "conf",
                     tdist = FALSE, # only works with autocor = FALSE
                     beep = FALSE,
                     ...){ # variables to nest under (station id e.g., in some cases variable name)

  nestvars <- enquos(...)
  datevar <- enquo(datevar)
  variable <- enquo(values)
  plan(multiprocess, gc = TRUE, workers = parallel::detectCores(logical = F))
  tictoc::tic()
  .data %>%
    mutate(variable = !!variable,
           date = !!datevar) %>%
    select(date, variable, !!!nestvars) %>%
    group_by(!!!nestvars, date) %>%
    summarise_at("variable", mean) %>%
    ungroup() %>%
    drop_na(variable) %>%
    group_by(!!!nestvars) %>%
    complete(date = full_seq(date, 1)) %>%
    mutate(decimaldate = decimal_date(date),
           month = month(date),
           decimal_during = decimaldate - year(date)) %>%
    nest() %>%
    ungroup() %>%
    mutate(
      fit = future_map(data, possibly(~ modeling(.x,
                                                 link = link,
                                                 autocor = autocor, tdist = tdist),
                                      otherwise = NA_integer_,
                                      quiet = F),
                       .progress = T),
      fderiv = map2(fit, data, possibly(~ fderiv(.x, term = "decimaldate", newdata = .y), otherwise = NA_integer_)),
      fderiv_confint = map2(fderiv, data, possibly(~ .x %>% confint(type = conf.type), otherwise = NA_integer_)),
      predict = map2(fit, data, possibly(~ predict(.x, newdata = .y, type = "terms") %>% as_tibble(), otherwise = NA_integer_)),
      fitted = map2(fit, data, possibly(~ predict(.x, newdata = .y, type = "link"), otherwise = NA_integer_)),
      autocor = map_lgl(fit, possibly(~.x$autocor, otherwise = NA_integer_)),
      intercept = map_dbl(fit, possibly(~ coef(.x) %>% .[1], otherwise = NA_integer_))
    ) %>%
    group_by(!!!nestvars) %>%
    dplyr::select(!!!nestvars, autocor, everything())->
    output
  tictoc::toc()
  if(beep){beepr::beep()}
  return(output)
}



plot_screeningtrends <- function(.output, y_id = NULL, sorting = NULL, wrappingvar = NULL){
  sorting <- enquo(sorting)
  if(rlang::quo_is_null(sorting)){sorting <- enquo(y_id)}
  y_id <- enquo(y_id)
  grouping <- group_vars(.output)
  wrapping <- enquo(wrappingvar)

  .output %>%
    ungroup() %>%
    mutate(!!y_id := reorder(!!y_id, !!sorting %>% desc)) %>%
    group_by(!!!syms(grouping)) %>%
    unnest(cols = c(fderiv_confint, data)) %>%
    select(!!!syms(grouping), date, est, lower, upper) %>%
    mutate(signif = !Vectorize(between)(0, lower, upper), # make a logical for significant change
      sign = sign(est)*signif) %>% # calculate sign and replace insignificant signs of derivatives with 0
    select(-lower, -upper, -est) %>%
    arrange(!!!syms(grouping), date) %>%
    nest() %>%
    group_by(!!!syms(grouping)) %>%
    mutate(
      periods = map(data, ~ .x %>%
                      transmute(id = paste0(signif, sign)) %$% # combine sign of significant estimates and significance
                      id %>% # pass onto
                      rle() %>% # find lengths of periods with the same values
                      .$lengths %>% # extract lengths
                      add(1) %>% # add on to find beginning of each period except first one
                      .[-length(.)] %>% # remove last one due to being the last rownumber+1
                      cumsum() %>% # calculate cumulative sums to find the positions
                      c(1, ., length(.x$sign) - 1) %>% # add the beginning of the first period, and the end of the last one
                      unique()), # filter out duplicates that may occur from last step
      data = map2(data, periods, ~ slice(.x, .y))
    ) %>% # filter out beginnings and ultimate end
    unnest(cols=c(data, periods)) %>% # each date is now a beginning or the ultimate end
    mutate(
      beginning = date,
      end = lead(date)
    ) %>% # this puts the next beginning date as the ending of the previous period
    filter(is.na(end) == F) %>% # filter out NAs, i.e. the last observation
    select(-periods, -date, -signif) %>% # select out only beginning, end, station, and significant sign of period
    ungroup() %>%
    mutate(sign = sign %>% factor(
      levels = c(-1, 0, 1),
      labels = c("Decreasing", "None", "Increasing"), ordered = T
    )) -> data

  data %>%
    ggplot(aes(
      y = !!y_id,
      ymin = as.numeric(!!y_id) - 0.4,
      ymax = as.numeric(!!y_id) + 0.4,
      xmin = beginning,
      xmax = end,
      fill = sign
    )) +
    geom_rect() +
    scale_x_date(
      date_labels = "%Y",
      expand = expansion(mult = 0.01),
      date_minor_breaks = "1 year", minor_breaks = waiver()
    ) +
    scale_fill_manual(values = c(
      "Decreasing" = "#56B4E9",
      "None" = "#F0E442",
      "Increasing" = "#D55E00"
    )) +
    labs(y = NULL, fill = "Significant trend changes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black")) ->
    plotobj

  if (!rlang::quo_is_null(wrapping)) { plotobj <- plotobj + facet_wrap(vars(!!wrapping)) }
  return(plotobj)
}

plot_screeningtrends_pvalues <- function(.output, y_id = NULL, sorting = NULL, wrappingvar = NULL){
  sorting <- enquo(sorting)
  if(rlang::quo_is_null(sorting)){sorting <- enquo(y_id)}
  y_id <- enquo(y_id)
  grouping <- group_vars(.output)
  wrapping <- enquo(wrappingvar)

  .output %>%
    ungroup() %>%
    mutate(!!y_id := reorder(!!y_id, !!sorting %>% desc)) %>%
    group_by(!!!syms(grouping)) %>%
    mutate(fderiv = map(fderiv, ~tibble(deriv = .x$derivatives[[1]][[1]], deriv_se = .x$derivatives[[1]][[2]]))) %>%
    unnest(cols = c(fderiv, fderiv_confint, data)) %>%
    arrange(!!!syms(grouping), date) %>%
    mutate(z = deriv/deriv_se,
           signif = case_when(lower > 0 ~ T,# make a logical for significant change
                              upper < 0 ~ T,
                              lower < 0 & upper > 0 ~ F),
           pvalue = (1-pnorm(abs(z)))*2,
           shannon = -log10(pvalue),
           shannon_sign = shannon * sign(deriv),
           shannon_sign = ifelse(is.infinite(shannon_sign), sign(shannon_sign)*-log10(.Machine$double.eps), shannon_sign),
           sign = sign(est)*signif) %>%
    mutate(
      beginning = date,
      end = lead(date)
    ) %>%
    filter(is.na(end) == F) %>%
    ungroup() -> data

  plotobj <- data %>%
    ggplot(aes(y=!!y_id,
               ymin=!!y_id %>% factor %>% as.numeric %>% subtract(0.4),
               ymax=!!y_id %>% factor %>% as.numeric %>% add(0.4),
               fill=shannon_sign, 
               xmin = beginning,
               xmax = end))+
    geom_rect()+
    scale_fill_gradient2(low = "#56B4E9", mid="#F0E442",
                         high = "#D55E00", midpoint = 0,
                         breaks=seq(-16,16,2),
                         labels=function(x){sign(x)*10^-abs(x)},
                         trans=trim_tails(range=c(log10(0.00001), -log10(0.00001))))+
    labs(fill="p-value *\nsign") +
    theme(axis.text.x = element_text(angle = 40, vjust = 0.9, colour = "black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank()
    )
  if (!rlang::quo_is_null(wrapping)) { plotobj <- plotobj + facet_wrap(vars(!!wrapping)) }
  return(plotobj)
}

plot_screeningtrends_relative <- function(.output, y_id = NULL, sorting = NULL, wrappingvar = NULL){
  sorting <- enquo(sorting)
  if(rlang::quo_is_null(sorting)){sorting <- enquo(y_id)}
  y_id <- enquo(y_id)
  grouping <- group_vars(.output)
  wrapping <- enquo(wrappingvar)

  .output %>%
    ungroup() %>%
    mutate(!!y_id := reorder(!!y_id, !!sorting %>% desc)) %>%
    group_by(!!!syms(grouping)) %>%
    mutate(fderiv = map(fderiv, ~tibble(deriv = .x$derivatives[[1]][[1]], deriv_se = .x$derivatives[[1]][[2]]))) %>%
    unnest(cols = c(predict, fderiv, data)) %>%
    arrange(!!!syms(grouping), date) %>%
    group_by(!!!syms(grouping)) %>%
    mutate(relative = round(deriv, 2)/(intercept+`s(decimaldate)`)) %>%
    mutate(
      beginning = date,
      end = lead(date)
    ) %>%
    ungroup %>%
    filter(is.na(end) == F) %>%
    ungroup() -> data

  plotobj <- data %>%
    ggplot(aes(y=!!y_id,
               ymin=!!y_id %>% factor %>% as.numeric %>% subtract(0.4),
               ymax=!!y_id %>% factor %>% as.numeric %>% add(0.4),
               fill=relative, 
               xmin = beginning,
               xmax = end))+
    geom_rect()+
    scale_fill_gradient2(low = "#56B4E9", mid = "#F0E442", high = "#D55E00", midpoint = 0, 
                         labels=scales::percent, trans = trim_tails(range=c(-1.5, 1.5)))+
    theme(axis.text.x = element_text(angle = 40, vjust = 0.9, colour = "black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank()
    )
  if (!rlang::quo_is_null(wrapping)) { plotobj <- plotobj + facet_wrap(vars(!!wrapping)) }
  return(plotobj)
}

plot_screeningtrends_reference <- function(.output, y_id = NULL, sorting = NULL, wrappingvar = NULL){

  sorting <- enquo(sorting)
  if(rlang::quo_is_null(sorting)){sorting <- enquo(y_id)}
  y_id <- enquo(y_id)
  grouping <- group_vars(.output)
  wrapping <- enquo(wrappingvar)


  .output %>%
    ungroup() %>%
    mutate(!!y_id := reorder(!!y_id, !!sorting %>% desc)) %>%
    group_by(!!!syms(grouping)) %>%
    mutate(fderiv = map(fderiv, ~tibble(deriv = .x$derivatives[[1]][[1]], deriv_se = .x$derivatives[[1]][[2]]))) %>%
    unnest(cols = c(predict, fderiv, data)) %>%
    arrange(!!!syms(grouping), date) %>%
    group_by(!!!syms(grouping)) %>%
    mutate(relative = (intercept+`s(decimaldate)`)/(mean(`s(decimaldate)`[decimaldate<decimaldate[1]+3])+intercept)) %>%
    mutate(
      beginning = date,
      end = lead(date)
    ) %>%
    ungroup %>%
    filter(is.na(end) == F) %>%
    ungroup() -> data

  plotobj <- data %>%
    ggplot(aes(y=!!y_id,
               ymin=!!y_id %>% factor %>% as.numeric %>% subtract(0.4),
               ymax=!!y_id %>% factor %>% as.numeric %>% add(0.4),
               fill=relative, 
               xmin = beginning,
               xmax = end))+
    geom_rect()+
    scale_fill_gradient2(low = "#56B4E9", mid="#F0E442",high = "#D55E00", midpoint = 1, 
                         labels=scales::percent, trans=trim_tails(range=c(0.1, 2)))+
    theme(axis.text.x = element_text(angle = 40, vjust = 0.9, colour = "black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank()
    )
  if (!rlang::quo_is_null(wrapping)) { plotobj <- plotobj + facet_wrap(vars(!!wrapping)) }
  return(plotobj)
}


plot_proportions <- function(.output, adjust = FALSE, 
                             wrappingvar = NULL){
  grouping <- group_vars(.output)
  wrapping <- enquo(wrappingvar)

  .output %>%
    unnest(cols = c(fderiv_confint, data)) %>%
    select(!!!syms(grouping), date, est, lower, upper) %>%
    mutate(signif = !Vectorize(between)(0, lower, upper), # make a logical for significant change
           sign = sign(est)*signif) %>% # calculate sign and replace insignificant signs of derivatives with 0
    select(-lower, -upper, -est) %>%
    arrange(!!!syms(grouping), date) %>%
    nest() %>%
    group_by(!!!syms(grouping)) %>%
    mutate(
      periods = map(data, ~ .x %>%
                      transmute(id = paste0(signif, sign)) %$% # combine sign of significant estimates and significance
                      id %>% # pass onto
                      rle() %>% # find lengths of periods with the same values
                      .$lengths %>% # extract lengths
                      add(1) %>% # add on to find beginning of each period except first one
                      .[-length(.)] %>% # remove last one due to being the last rownumber+1
                      cumsum() %>% # calculate cumulative sums to find the positions
                      c(1, ., length(.x$sign) - 1) %>% # add the beginning of the first period, and the end of the last one
                      unique()), # filter out duplicates that may occur from last step
      data = map2(data, periods, ~ slice(.x, .y))
    ) %>% # filter out beginnings and ultimate end
    unnest(cols=c(data, periods)) %>% # each date is now a beginning or the ultimate end
    mutate(
      beginning = date,
      end = lead(date)
    ) %>% # this puts the next beginning date as the ending of the previous period
    filter(is.na(end) == F) %>% # filter out NAs, i.e. the last observation
    select(-periods, -date, -signif) %>% # select out only beginning, end, station, and significant sign of period
    ungroup() %>%
    mutate(sign = sign %>% factor(
      levels = c(-1, 0, 1) %>% rev,
      labels = c("Decreasing", "None", "Increasing") %>% rev, ordered = T
    )) %>%
    mutate(date = Vectorize(seq.Date)(from = beginning, to = end - 1, 1)) %>%
    unnest(date) %>%
    # group_by(date) %>%
    select(!!!syms(grouping), sign, date) ->
    data1
  if (rlang::quo_is_null(wrapping) == FALSE) {
    data1 %>%
      select(!!!syms(grouping)) %>%
      distinct() %>%
      group_by(!!wrapping) %>% #### deprecation warning
      summarise(n_stations = n()) ->
      number_of_stations

    data1 %>% count(!!wrapping, date) %>% full_join(number_of_stations) %>%
      mutate(p_stations = n/max(n_stations)) %>%
      select(-n) -> stations_date_proportion

    data1  %>%
      group_by(!!wrapping) %>%
      count(date, sign, .drop = F) %>%
      group_by(!!wrapping, date) %>%
      mutate(proportion = n / sum(n)) %>%
      full_join(stations_date_proportion) ->
      data
  }else{
    data1 %>%
      select(!!!syms(grouping)) %>%
      distinct() %>%
      summarise(n_stations = n()) ->
      number_of_stations

    data1 %>% count(date) %>% mutate(n_stations = number_of_stations$n_stations) %>%
      mutate(p_stations = n/max(n_stations)) %>%
      select(-n) -> stations_date_proportion

    data1  %>%
      count(date, sign, .drop = F) %>%
      group_by(date) %>%
      mutate(proportion = n / sum(n)) %>%
      full_join(stations_date_proportion) ->
      data
  }



  if(adjust == TRUE){
    data %>%
      mutate(proportion=proportion*p_stations) ->
      data
  }

  plotobj <- data %>%
    ggplot(aes(x = date, y = proportion, fill = sign)) +
    scale_y_continuous(labels = scales::percent) +
    geom_col(position = "stack", width = 1) +
    scale_fill_manual(values = c(
      "Decreasing" = "#56B4E9",
      "None" = "#F0E442",
      "Increasing" = "#D55E00"
    )) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"),
      panel.ontop = TRUE,
      panel.grid.major.y = element_line(colour = alpha("white", 0.2)),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(colour = alpha("white", 0.3))
    ) +
    labs(x = "Time", y = "Proportion of Stations", fill = "Significant Trend Changes",
         caption = "Dashed line shows proportion of number of stations for that timepoint")+
    geom_line(inherit.aes = F,
              mapping=aes(x=date, y=p_stations), linetype="dashed")

  if (rlang::quo_is_null(wrapping) == FALSE) { plotobj <- plotobj + facet_wrap(vars(!!wrapping)) }
  return(plotobj)
}


plot_data <- function(.output){
  grouping <- group_vars(.output)

  .output %>%
    unnest(cols = c(fderiv_confint, data)) %>%
    select(!!!syms(grouping), date, est, lower, upper) %>%
    mutate(signif = !Vectorize(between)(0, lower, upper), # make a logical for significant change
           sign = sign(est)*signif) %>% # calculate sign and replace insignificant signs of derivatives with 0
    select(-lower, -upper, -est) %>%
    arrange(!!!syms(grouping), date) %>%
    nest() %>%
    group_by(!!!syms(grouping)) %>%
    mutate(
      periods = map(data, ~ .x %>%
                      transmute(id = paste0(signif, sign)) %$% # combine sign of significant estimates and significance
                      id %>% # pass onto
                      rle() %>% # find lengths of periods with the same values
                      .$lengths %>% # extract lengths
                      add(1) %>% # add on to find beginning of each period except first one
                      .[-length(.)] %>% # remove last one due to being the last rownumber+1
                      cumsum() %>% # calculate cumulative sums to find the positions
                      c(1, ., length(.x$sign) - 1) %>% # add the beginning of the first period, and the end of the last one
                      unique()), # filter out duplicates that may occur from last step
      data = map2(data, periods, ~ slice(.x, .y))
    ) %>% # filter out beginnings and ultimate end
    unnest(cols=c(data, periods)) %>% # each date is now a beginning or the ultimate end
    mutate(
      beginning = date,
      end = lead(date)
    ) %>% # this puts the next beginning date as the ending of the previous period
    filter(is.na(end) == F) %>% # filter out NAs, i.e. the last observation
    select(-periods, -date, -signif) %>% # select out only beginning, end, station, and significant sign of period
    ungroup() %>%
    mutate(sign = sign %>% factor(
      levels = c(-1, 0, 1) %>% rev,
      labels = c("Decreasing", "None", "Increasing") %>% rev, ordered = T
    )) -> data

  return(data)
}

splines_and_derivative <- function(gam_object, n_eval = 200, n_sim = 100, eps = 0.0000001) {
  number_of_smooths <- gam_object$smooth %>% length()
  data <- model.frame(gam_object)

  Vc <- vcov(gam_object,
             unconditional = TRUE
  )

  tibble(id = 1:number_of_smooths) %>%
    mutate(
      smooth_obj = map(id, ~ gam_object$smooth[[.x]]),
      dims = map_dbl(smooth_obj, ~ .$dim)
    ) %>%
    filter(dims == 1) %>%
    mutate(
      label = map_chr(smooth_obj, ~ .x$label),
      term = map_chr(smooth_obj, ~ .x$term),
      by = map_chr(smooth_obj, ~ .x$by),
      data = map(term, ~ tibble(x = data[[.x]])),
      grid = ifelse(by != "NA",
                    map(data, ~ tibble(x = seq(floor(min(.x[[1]])), ceiling(max(.x[[1]])), length.out = n_eval)) %>% mutate(by = 1)),
                    map(data, ~ tibble(x = seq(floor(min(.x[[1]])), ceiling(max(.x[[1]])), length.out = n_eval)))),
      grid = pmap(
        list(grid, term, by),
        function(grid, term, by) {
          if (ncol(grid) == 1) {
            names(grid) <- term
          } else {
            names(grid) <- c(term, by)
          }
          return(grid)
        }
      ),
      grid2 = ifelse(by != "NA",
                     map(data, ~ tibble(x = seq(floor(min(.x[[1]])), ceiling(max(.x[[1]])), length.out = n_eval)) %>% mutate(x = x + eps, by = 1)),
                     map(data, ~ tibble(x = seq(floor(min(.x[[1]])), ceiling(max(.x[[1]])), length.out = n_eval)) %>% mutate(x = x + eps))
      ),
      grid2 = pmap(
        list(grid2, term, by),
        function(grid, term, by) {
          if (ncol(grid) == 1) {
            names(grid) <- term
          } else {
            names(grid) <- c(term, by)
          }
          return(grid)
        }
      ),
      first = map_dbl(smooth_obj, ~ .x$first.para),
      last = map_dbl(smooth_obj, ~ .x$last.para),
      coefs = map2(first, last, ~ coef(gam_object)[.x:.y]),
      Vc = map2(first, last, ~ Vc[.x:.y, .x:.y]),
      Bu = map(Vc, ~ MASS::mvrnorm(
        n = n_sim,
        mu = rep(
          0,
          nrow(.x)
        ),
        Sigma = .x
      )),
      knotvalues = map2(smooth_obj, grid, ~ PredictMat(.x, .y)),
      knotvalues_deriv = map2(smooth_obj, grid2, ~ PredictMat(.x, .y)),
      knotvalues_deriv = map2(knotvalues, knotvalues_deriv, ~ (.y - .x) / eps),
      Se = map2(knotvalues, Vc, ~ sqrt(rowSums(.x * (.x %*% .y)))),
      Se_deriv = map2(knotvalues_deriv, Vc, ~ sqrt(rowSums(.x * (.x %*% .y)))),
      crit = pmap_dbl(list(knotvalues, Bu, Se), function(Xp, Bu, Se) {
        quantile(apply(abs((Xp %*% t(Bu)) / Se), 2, max),
                 type = 8,
                 prob = 0.95
        )
      }),
      crit_deriv = pmap_dbl(list(knotvalues_deriv, Bu, Se_deriv), function(Xp, Bu, Se) {
        quantile(apply(abs((Xp %*% t(Bu)) / Se), 2, max),
                 type = 8,
                 prob = 0.95
        )
      }),
      est = map2(knotvalues, coefs, ~ .x %*% .y),
      est_deriv = map2(knotvalues_deriv, coefs, ~ .x %*% .y),
      error_margin = map2(crit, Se, ~ .x * .y),
      error_margin_deriv = map2(crit_deriv, Se_deriv, ~ .x * .y),
      pval = map2(est, error_margin, ~ (1 - pnorm(abs(.x / (.y * 0.5102137))))),
      pval_deriv = map2(est_deriv, error_margin_deriv, ~ (1 - pnorm(abs(.x / (.y * 0.5102137))))),
      grid = map2(.x = grid, .y = term, ~rename_at(.x, .y, ~"x"))
    ) %>%
    unnest(cols=c(grid, Se, Se_deriv, est, est_deriv, error_margin, error_margin_deriv, pval, pval_deriv),
           names_repair = "universal") %>%
    mutate(
      error_margin_point = 1.96 * Se,
      error_margin_point_deriv = 1.96 * Se_deriv,
      pval_point = 1 - pnorm(abs(est / (Se))),
      pval_point_deriv = 1 - pnorm(abs(est_deriv / (Se_deriv)))
    ) %>%
    unnest(c(est, pval_point, pval_point_deriv, pval, pval_deriv, est_deriv)) %>%
    select_if(~!(class(.))%in%"list") ->
    output


  class(output) <- c("splines", class(output))
  return(output)
}

plot.splines <- function(splines_obj) {
  deriv_vars <- c(
    "crit_deriv",
    "Se_deriv",
    "est_deriv",
    "error_margin_deriv",
    "pval_deriv",
    "error_margin_point_deriv",
    "pval_point_deriv"
  )
  spline_vars <- c(
    "crit",
    "Se",
    "est",
    "error_margin",
    "pval",
    "error_margin_point",
    "pval_point"
  )

  splines_obj %>%
    select(-one_of(spline_vars)) %>%
    rename(
      crit = crit_deriv,
      Se = Se_deriv,
      est = est_deriv,
      error_margin = error_margin_deriv,
      pval = pval_deriv,
      error_margin_point = error_margin_point_deriv,
      pval_point = pval_point_deriv
    ) %>%
    mutate(type = "deriv") %>%
    select_if(~!(class(.))%in%"list")->
    deriv_data

  splines_obj %>%
    select(-one_of(deriv_vars)) %>%
    mutate(type = "spline") %>%
    select_if(~!(class(.))%in%"list")->
    spline_data

  spline_data %>%
    full_join(deriv_data) %>%
    mutate(type = type %>% factor(levels = c("spline", "deriv"))) %>%
    ggplot(aes(x = x, y = est)) +
    geom_line() +
    facet_wrap(type ~ label, scales = "free", nrow = 2) +
    geom_ribbon(alpha = 0.4, mapping = aes(ymin = est - error_margin, ymax = est + error_margin)) +
    geom_ribbon(alpha = 0.3, mapping = aes(ymin = est - error_margin_point, ymax = est + error_margin_point))
}


plot_individual_trend <- function(x, y=NULL, title=NULL){
  if(nrow(x) != 1){stop("Filter out the variable (and/or station) you are interested in.")}
  annualterm <- predict(x$fit[[1]], newdata=x$data[[1]], type="terms")[,1]
  intercept <- x$fit[[1]]$coef["(Intercept)"]
  x$fderiv_confint[[1]] %>%
    as_tibble %>%
    rowwise %>%
    mutate(signif = !between(0, lower, upper),
           sign=sign(est),
           signif_sign = signif*sign) %>%
    ungroup %>% bind_cols(x$data[[1]],.) %>%
    mutate(trend = annualterm+intercept) -> dt
    if(x$fit[[1]]$family$link[1]=="log"){dt$trend=exp(dt$trend)}
    dt%>%
    drop_na(variable) %>%
    ggplot(aes(x=date))+geom_line(aes(y=variable))+
    geom_line(aes(y=trend, color=as_factor(signif_sign), group=c(0)), size=1.5)+
    scale_color_manual(values=c("-1" = "#56B4E9",
                                "0" = "#F0E442",
                                "1" = "#D55E00"))+
    theme_bw()+
    labs(x="Date", y=y,title=title)+
    theme(legend.position = "none") -> p
  return(p)
}


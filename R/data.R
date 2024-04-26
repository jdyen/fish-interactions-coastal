# internal function to download fish data from AAEDB
fetch_fish <- function(target = "northern", recompile = FALSE) {
  
  # check target is `southern` or `northern`
  if (!target %in% c("southern", "northern")) 
    stop("target must be one of `southern` or `northern`", call. = FALSE)
  
  # check if data exist
  fish_exists <- any(grepl(paste0("fish-compiled-", target, ".qs"), dir("data/")))
  
  # if data exist and !recompile, load saved version. Re-extract otherwise
  if (fish_exists & !recompile) {
    
    # load data
    cpue <- qread(paste0("data/fish-compiled-", target, ".qs"))
    
  } else {
    
    # set appropriate filters
    .species <- .species_northern
    .waterbodies <- .waterbodies_northern
    if (target == "southern") {
      .species <- .species_southern
      .waterbodies <- .waterbodies_southern
    }

    # grab cpue data from AAEDB, filtered to targets and
    #   for three different size categories
    categories <- c("small", "medium", "large")
    thresholds <- c(0, 13.5, 50.0, Inf)
    if (target == "southern") 
      thresholds <- c(0, 8.0, 16.0, Inf)
    cpue <- vector("list", length = length(thresholds))
    for (i in seq_along(thresholds)) {
      cpue[[i]] <- fetch_cpue(
        c(1, 2, 4, 6, 9, 10:13, 15:50),
        criterion = list(var = "length_cm", lower = thresholds[i], upper = thresholds[i + 1])
      ) |>
        filter(
          scientific_name %in% !!.species,
          waterbody %in% !!.waterbodies,
          survey_year >= 2010
        ) |>
        collect() |>
        mutate(category = categories[i])
    }

    # collapse these into a single data set, with a comparison
    #    (unfiltered by sizes so we can check values are consistent)
    cpue <- bind_rows(
      cpue,
      fetch_cpue(c(1, 2, 4, 6, 9, 10:13, 15:50)) |>
      filter(
        scientific_name %in% !!.species,
        waterbody %in% !!.waterbodies,
        survey_year >= 2010
      ) |>
      collect() |>
      mutate(category = "comparison")
    )    
        
    # group some species together
    cpue <- cpue |>
      mutate(
        species = scientific_name,
        species = ifelse(grepl("Philypnodon", species), "Philypnodon spp.", species),
        species = ifelse(grepl("unidentified eel|Anguilla", species), "Anguilla spp.", species),
        species = ifelse(grepl("Craterocephalus", species), "Craterocephalus spp.", species),
        species = ifelse(grepl("Carassius|Cyprinus", species), "Carp Goldfish", species)
      )
    
    # need to sum up multiple records of species groups at a site
    cpue <- cpue |>
      group_by(id_project, id_survey, species, category) |>
      summarise(catch = sum(catch)) |>
      ungroup() |>
      left_join(
        cpue |> 
          distinct(id_project, id_survey, species, category, waterbody, id_site, survey_date, survey_year, gear_type, effort_h),
        by = c("id_project", "id_survey", "species", "category")
      )
    
    # add some site info
    site_info <- cpue |> fetch_site_info() |> collect()
    st_geometry(site_info) <- st_as_sfc(site_info$geom_pnt, crs = 4283)
    
    # ignored for now, can use VEWH reach table to add reach info    
    vewh_reaches <- fetch_table("eflow_reaches_20171214", "projects") |>
      collect()
    st_geometry(vewh_reaches) <- st_as_sfc(vewh_reaches$geom, crs = 4283)
    site_info <- site_info |>
      st_join(vewh_reaches |> select(vewh_reach), join = st_within) |>
      mutate(reach_no = ifelse(is.na(reach_no), vewh_reach, reach_no)) |>
      select(-vewh_reach)
    
    # grab a few Ovens sites and duplicate for id_project = 9
    ovens_sub <- site_info |>
      filter(id_site %in% c(3160, 3162, 3163, 3183, 3185)) |>
      mutate(id_project = 9)
    site_info <- bind_rows(site_info, ovens_sub)
    
    # add reach info for unknown reaches and then append to CPUE data
    site_info <- site_info |>
      mutate(
        id_site = as.numeric(id_site),
        reach_no = ifelse(id_site %in% c(4464, 2058), 3, reach_no),
        reach_no = ifelse(id_site %in% c(4466, 4467, 4470, 4471), 4, reach_no),
        reach_no = ifelse(id_site %in% c(4468, 4066, 4068, 4069, 4073), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4067, 4070:4072), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3193), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3194, 4266), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4060, 4061), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4382, 4383), 3, reach_no),
        reach_no = ifelse(id_site %in% c(4109, 4110, 4115), 4, reach_no),
        reach_no = ifelse(id_site %in% c(4116), 5, reach_no),
        reach_no = ifelse(id_site %in% c(3133), 0, reach_no),
        reach_no = ifelse(id_site %in% c(3134), 5, reach_no),
        reach_no = ifelse(id_site %in% c(1643:1646, 3164, 4225, 4229, 4231), 0, reach_no),
        reach_no = ifelse(id_site %in% c(3160:3163, 3182, 4172:4185, 4194:4197, 4199:4204, 4208:4212, 4217:4224, 4228, 4232:4241), 4, reach_no),
        reach_no = ifelse(id_site %in% c(3322, 3324, 3325), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3167), 3, reach_no),
        reach_no = ifelse(id_site %in% c(3112, 3177:3180, 4451), 5, reach_no),
        reach_no = ifelse(id_site %in% c(3113, 3181), 4, reach_no),
        reach_no = ifelse(waterbody == "Broken Creek" & is.na(reach_no), 5, reach_no),
        reach_no = ifelse(waterbody == "Ovens River" & is.na(reach_no), 5, reach_no),
        reach_no = ifelse(id_site %in% c(2302, 2306, 3049), 2, reach_no),
        reach_no = ifelse(id_site %in% c(1880, 1783, 3041, 2096, 2398), 3, reach_no),
        reach_no = ifelse(id_site %in% c(4603, 2239, 2926), 4, reach_no),
        reach_no = ifelse(id_site %in% c(2231), 5, reach_no)
      )
    
    # add reach and lat/long info, filter to remove Buffalo sites and some sites
    #   really far upstream and add reach info for Ovens sites in id_project 9
    cpue <- cpue |>
      left_join(
        site_info |>
          select(id_project, id_site, reach_no, latitude, longitude),
        by = c("id_project", "id_site")
      ) |>
      filter(!(id_site %in% c(4226, 4227, 4230, 1864, 2347, 2338, 1799, 2385))) |>
      mutate(waterbody = ifelse(id_site %in% c(3134), "Thomson River", waterbody))
    
    # filter out some upper reaches of the Moorabool and Macalister (0, 2, 0)
    if (target == "southern") {
    cpue <- cpue |>
      filter(
        !(waterbody == "Macalister River" & reach_no == 0),
        !(waterbody == "Moorabool River" & reach_no %in% c(0, 2))
      )
    }
    
    # filter out sites without geom information
    cpue <- cpue |>
      filter(
        id_site %in% !!(site_info |> filter(!is.na(geom_pnt)) |> pull(id_site) |> unique())
      )
    
    # remove sites without category, replace small-bodied species estiamtes 
    #  with total catch
    cpue <- cpue |>
      filter(!is.na(category))
    
    # filter to upper Loddon if southern or lower Loddon otherwise
    if (target == "southern") {
      cpue <- cpue |>
        filter(
          !(waterbody == "Loddon River" & reach_no > 2),
          !(waterbody == "Macalister River" & reach_no != 1),
          !(waterbody == "MacKenzie River" & reach_no != 2),
          !(waterbody == "Moorabool River" & reach_no == 1),
          !(waterbody == "Thomson River" & reach_no > 3)
        )
    } else {
      cpue <- cpue |>
        filter(
          !(waterbody == "Loddon River" & reach_no %in% c(2, 3)),
          !(waterbody == "Broken River" & reach_no %in% c(1, 2)),
          !(waterbody == "Campaspe River" & reach_no %in% c(2, 3)),
          !(waterbody == "Ovens River" & reach_no == 0)
        )
    }
    
    # issues in these SB species
    if (target == "northern") {
      sb_species <- c(
        "Retropinna semoni",
        "Melanotaenia fluviatilis"
      )
      cpue_compare <- cpue |>
        filter(category == "comparison", species %in% sb_species) |>
        mutate(category = "small")
      cpue <- cpue |>
        left_join(
          cpue_compare |> 
            select(id_survey, species, category, catch) |>
            rename(catch_sb = catch),
          by = c("id_survey", "species", "category")
        ) |>
        mutate(catch = ifelse(!is.na(catch_sb), catch_sb, catch))
    }

    # save this
    qsave(cpue, file = paste0("data/fish-compiled-", target, ".qs"))
    
  }
  
  # return
  cpue
  
}

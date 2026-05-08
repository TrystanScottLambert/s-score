# Script to calculate the Stotal metric from Robotham+2011
#
# Three required parquet files:
#  - generated group catalog with id_group, number_members
#  - generated galaxy ctalog with id_galaxy, id_group
#  - mock galaxy catalog with id_galaxy_sky, id_fof (since we are using shark)

gen_gal_path   <- "generated_galaxies.parquet"
gen_group_path <- "generated_groups.parquet"
mock_gal_path  <- "mock_galaxies.parquet"

min_group_size <- 5L


bijcheck <- function(group_ids_1, group_ids_2, min_group_size = 2L) {
  stopifnot(length(group_ids_1) == length(group_ids_2))

  # Frequency tables, ignoring -1 (ungrouped)
  count_table_1 <- table(group_ids_1[group_ids_1 != -1])
  count_table_2 <- table(group_ids_2[group_ids_2 != -1])

  # Groups in catalogue 1 that pass the size cut
  valid_groups_1 <- as.integer(names(count_table_1))[count_table_1 >= min_group_size]

  if (length(valid_groups_1) == 0L) {
    return(list(e_num = 0L, e_den = 0L, q_num = 0, q_den = 0))
  }

  # Iterate over every valid group in catalogue 1
  per_group <- lapply(valid_groups_1, function(group_id) {
    members      <- which(group_ids_1 == group_id)
    overlap_grps <- group_ids_2[members]
    overlap_val  <- overlap_grps[overlap_grps != -1]

    n1 <- as.numeric(count_table_1[as.character(group_id)])

    if (length(overlap_val) > 0L) {
      temptab <- table(overlap_val)
      grp2    <- as.integer(names(temptab))
      counts  <- as.integer(temptab)

      n2 <- as.numeric(count_table_2[as.character(grp2)])

      frac_1 <- counts / n1
      frac_2 <- counts / n2

      # Treat each ungrouped (-1) member as its own singleton match,
      # mirroring the Rust reference implementation.
      n_iso <- sum(overlap_grps == -1)
      if (n_iso > 0L) {
        frac_1 <- c(frac_1, rep(1 / n1, n_iso))
        frac_2 <- c(frac_2, rep(1,      n_iso))
      }

      # Best two-way match = largest purity product 
      # which.max returns the FIRST index of the max, matching the Rust code.
      best <- which.max(frac_1 * frac_2)
      list(q1 = frac_1[best], q2 = frac_2[best], n1 = n1)
    } else {
      # All overlapping members are ungrouped in catalogue 2
      list(q1 = 1 / n1, q2 = 1, n1 = n1)
    }
  })

  q1_vec <- vapply(per_group, `[[`, numeric(1), "q1")
  q2_vec <- vapply(per_group, `[[`, numeric(1), "q2")
  n1_vec <- vapply(per_group, `[[`, numeric(1), "n1")

  list(
    e_num = sum(q1_vec > 0.5 & q2_vec > 0.5),
    e_den = length(n1_vec),
    q_num = sum(q1_vec * n1_vec),
    q_den = sum(n1_vec)
  )
}


# Checking that group ids that only appear once get labeled as -1
# This is to ensure that groups that might only have one member because of cuts
# or because users forget to account for this, are handled.
drop_singletons <- function(group_ids, label = "groups") {
  real <- group_ids != -1
  tab  <- table(group_ids[real])
  singletons <- as.integer(names(tab))[tab == 1L]
  if (length(singletons) > 0L) {
    message(sprintf("  %s: reassigning %d singleton ids to -1",
                    label, length(singletons)))
    group_ids[group_ids %in% singletons] <- -1L
  } else {
    message(sprintf("  %s: no singletons found", label))
  }
  group_ids
}


message("Reading parquet files...")
gen_gal   <- read_parquet(gen_gal_path)
gen_group <- read_parquet(gen_group_path)
mock_gal  <- read_parquet(mock_gal_path)

message(sprintf("  generated galaxies : %d", nrow(gen_gal)))
message(sprintf("  generated groups   : %d", nrow(gen_group)))
message(sprintf("  mock galaxies      : %d", nrow(mock_gal)))


message("Matching galaxies between generated and mock catalogues...")

joined <- gen_gal %>%
  select(id_galaxy, id_group) %>%
  inner_join(
    mock_gal %>% select(id_galaxy_sky, id_fof),
    by = c("id_galaxy" = "id_galaxy_sky")
  )

message(sprintf("  matched galaxies   : %d", nrow(joined)))

if (nrow(joined) == 0L) {
  stop("No galaxies matched between the generated and mock catalogues. ",
       "Check that id_galaxy and id_galaxy_sky use the same convention.")
}

message("Cleaning singleton groups...")
gen_ids  <- drop_singletons(gen_ids,  label = "generated (id_group)")
mock_ids <- drop_singletons(mock_ids, label = "mock      (id_fof)  ")


message("Running bijective comparison...")
fof_vs_mock  <- bijcheck(gen_ids,  mock_ids, min_group_size)
mock_vs_fof  <- bijcheck(mock_ids, gen_ids,  min_group_size)


E_FoF  <- fof_vs_mock$e_num  / fof_vs_mock$e_den
E_mock <- mock_vs_fof$e_num  / mock_vs_fof$e_den
Q_FoF  <- fof_vs_mock$q_num  / fof_vs_mock$q_den
Q_mock <- mock_vs_fof$q_num  / mock_vs_fof$q_den

E_tot <- E_FoF * E_mock
Q_tot <- Q_FoF * Q_mock
S_tot <- E_tot * Q_tot

cat("\n--- Robotham+2011 bijective comparison ---\n")
cat(sprintf("  E_FoF  = %.4f   (%d / %d generated groups bijective)\n",
            E_FoF, fof_vs_mock$e_num, fof_vs_mock$e_den))
cat(sprintf("  E_mock = %.4f   (%d / %d mock groups bijective)\n",
            E_mock, mock_vs_fof$e_num, mock_vs_fof$e_den))
cat(sprintf("  Q_FoF  = %.4f\n", Q_FoF))
cat(sprintf("  Q_mock = %.4f\n", Q_mock))
cat("------------------------------------------\n")
cat(sprintf("  E_tot  = %.4f\n", E_tot))
cat(sprintf("  Q_tot  = %.4f\n", Q_tot))
cat(sprintf("  S_tot  = %.4f\n", S_tot))
cat("------------------------------------------\n")

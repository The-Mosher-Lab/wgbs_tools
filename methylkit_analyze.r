#!/usr/bin/env Rscript

# Author: Jeffrey Grover
# Purpose: Run methylKit on some files output by MethylDackel
# Created: 2019-06-04
# Usage: Load a .methylKit file from MethylDackel and execute methylKit's DMR
# calling on it to output windows of differential methylation

library('argparser', quietly = TRUE)


# This function will analyze methylkit files with the desired parameters
methylkit_analyze <- function(control_files, experimental_files, sample_id_control,
                              sample_id_experimental, c_context, min_coverage,
                              window_size, step_size, cores, q_val, diff,
                              output_dir, pool) {

  # Load the methylkit library here
  # It's bad generally to do this but it takes a while to load.
  suppressMessages(library('methylKit', quietly = TRUE))

  # Read the files
  meth_obj <- methRead(
    location = as.list(c(control_files, experimental_files)),
    sample.id = as.list(c(
      rep(sample_id_control, times = length(control_files)),
      rep(sample_id_experimental, times = length(experimental_files))
      )
    ),
    assembly = 'unimportant_unnecessary_option',
    treatment = c(
      rep(0, times = length(control_files)),
      rep(1, times = length(experimental_files))
    ),
    context = c_context
  )

  # Pool samples if desired
  if (pool == TRUE) {
    meth_obj <- pool(meth_obj)
  }

  # Use a read coverage filter if desired
  if (min_coverage > 0) {
    meth_obj <- filterByCoverage(meth_obj, lo.count = min_coverage)
  }

  # normalizeCoverage for the filtered object
  meth_obj <- normalizeCoverage(meth_obj)

  # Create windows
  meth_obj <- tileMethylCounts(meth_obj,
                               win.size = window_size,
                               step.size = step_size,
                               mc.cores = cores)

  # Unite the tiles
  meth_obj <- unite(meth_obj, mc.cores = cores)

  # Get differentially methylated windows, using multicore support
  meth_diff <- calculateDiffMeth(meth_obj, mc.cores = cores)

  meth_diff_windows <- getMethylDiff(
    meth_diff,
    difference = diff,
    qvalue = q_val
  )

  meth_diff_hyper <- getMethylDiff(
    meth_diff,
    difference = diff,
    qvalue = q_val,
    type = "hyper"
  )

  meth_diff_hypo <- getMethylDiff(
    meth_diff,
    difference = diff,
    qvalue = q_val,
    type = "hypo"
  )

  # Print out the results
  print(paste0(c_context, ' DMWs:', nrow(meth_diff_windows)))
  print(paste0(c_context, ' Hyper-DMWs:', nrow(meth_diff_hyper)))
  print(paste0(c_context, ' Hypo-DMWs:', nrow(meth_diff_hypo)))

  # Export results
  write.csv(meth_diff_windows,
            paste0(output_dir, '/', sample_id_experimental, '_', sample_id_control,
                  '_',c_context, '_norm_window_', window_size, '_', step_size,
                  '_d', diff, '_q', q_val, '.csv'
                  )
            )
  write.csv(meth_diff_hyper,
            paste0(output_dir, '/', sample_id_experimental, '_', sample_id_control,
                  '_',c_context, '_norm_window_', window_size, '_', step_size,
                  '_d', diff, '_q', q_val, '_hyper.csv'
                  )
            )
  write.csv(meth_diff_hypo,
            paste0(output_dir, '/', sample_id_experimental, '_', sample_id_control,
                  '_',c_context, '_norm_window_', window_size, '_', step_size,
                  '_d', diff, '_q', q_val, '_hypo.csv'
                  )
            )
}


# Parse command-line arguments

get_args <- function () {
  parser <- arg_parser('Analyze input files with methylKit and determine differentially methylated windows.')
  parser <- add_argument(parser, '--control',
                        help='Control sample files, comma separated')
  parser <- add_argument(parser, '--experimental',
                        help='Experimental sample files, comma separated')
  parser <- add_argument(parser, '--control_id',
                        help='Sample ID for the control sample')
  parser <- add_argument(parser, '--experimental_id',
                        help='Sample ID for the experimental sample')
  parser <- add_argument(parser, '--context',
                        help='The cytosine context for the comparison to be done',
                        default='CG')
  parser <- add_argument(parser, '--min_cov',
                        help='Minimum number of reads to filter sites by',
                        type='integer',
                        default=0)
  parser <- add_argument(parser, '--window_size',
                        help='Size of the window to use for DMR calling',
                        type='integer',
                        default=300)
  parser <- add_argument(parser, '--step_size',
                        help='How far the windows should tile during DMR calling',
                        type='integer',
                        default=100)
  parser <- add_argument(parser, '--threads',
                        help='Number of CPU threads to use for DMR calling',
                        default=1)
  parser <- add_argument(parser, '--q_val',
                        help='FDR-corrected p-value (q-value) for DMR calls',
                        type='numeric',
                        default=0.05)
  parser <- add_argument(parser, '--diff_meth',
                        help='Difference in methylation (%) for DMR calls',
                        type='numeric',
                        default=25)
  parser <- add_argument(parser, '--pool',
                         help='Combine replicates',
                         default=FALSE)
  return(parse_args(parser))
}


# Main definition

main <- function (args) {
  control_files <- unlist(strsplit(args$control, ','))
  experimental_files <- unlist(strsplit(args$experimental, ','))

  output_dir <- 'methylkit_analyze'
  dir.create(output_dir)

  # Run the analysis and output results

  methylkit_analyze(control_files, experimental_files, args$control_id,
                    args$experimental_id, args$context, args$min_cov,
                    args$window_size, args$step_size, args$threads, args$q_val,
                    args$diff_meth, output_dir, args$pool)
}


# Call main

main(get_args())

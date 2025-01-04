
folder2do <- list.dirs( 'SJARACNe', recursive = FALSE)
cmd <- lapply( c( 'TF', 'SIG'), function( type_1) {
  paste( 'sjaracne', 'local', '-e', shQuote( sapply( folder2do, list.files, '\\.exp\\.txt$', full.names = TRUE)), '-g', shQuote( sapply( file.path( folder2do, type_1), list.files, full.names = TRUE)), '-o', shQuote( file.path( folder2do, type_1)), '-n 100 -pc 1e-2 -pb 1e-5')
})
cmd <- unlist( cmd)
.null <- parallel::mclapply( cmd, system, mc.cores = 2)







.cluster_muts <- function(u,a = NULL, mut_variable = 'mutations') 
{
  tu = descendantSids[[u]]
  mdf.u = amd[ amd$sequence_name %in% tu, ]
  ## find comparator ancestor 
  if ( is.null(a))
    a = .get_comparator_ancestor(u)
  asids = setdiff( descendantSids[[a]] ,  descendantSids[[u]] )
  mdf.a = amd[ amd$sequence_name %in% asids , ]
  if ( nrow( mdf.a ) == 0  |  nrow( mdf.u ) == 0 )
    return( list(defining = NA, all = NA ) )
  vtabu = sort( table( do.call( c, strsplit( mdf.u[[mut_variable]], split='\\|' )  ) ) / nrow( mdf.u ) )
  vtaba = sort( table( do.call( c, strsplit( mdf.a[[mut_variable]], split='\\|' )  ) ) / nrow( mdf.a ) )
  
  umuts = names( vtabu[ vtabu > mutation_cluster_frequency_threshold ] )
  defining_muts = setdiff( names( vtabu[ vtabu > mutation_cluster_frequency_threshold ] )
                           , names(vtaba[ vtaba > mutation_cluster_frequency_threshold ]) )
  list(defining=defining_muts, all=umuts  ) 
}

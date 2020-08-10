#dups between groups function
# dupsBetweenGroups <- function (df, idcol) {
#     # df: the data frame
#     # idcol: the column which identifies the group each row belongs to
# 
#     # Get the data columns to use for finding matches
#     datacols <- setdiff(names(df), idcol)
#     datacols <- setdiff(names(allSamplesDE), "Sample")
# 
#     # Sort by idcol, then datacols. Save order so we can undo the sorting later.
#     sortorder <- do.call(order, df)
#     df <- df[sortorder,]
#     sortorder <- do.call(order, allSamplesDE)
#     allSamplesDE <- allSamplesDE[sortorder,]
# 
#     # Find duplicates within each id group (first copy not marked)
#     dupWithin <- duplicated(df)
#     dupWithin <- duplicated(allSamplesDE)
#     # With duplicates within each group filtered out, find duplicates between groups. 
#     # Need to scan up and down with duplicated() because first copy is not marked.
#     dupBetween = rep(NA, nrow(allSamplesDE))
#     dupBetween[!dupWithin] <- duplicated(allSamplesDE[!dupWithin,])
#     dupBetween[!dupWithin] <- duplicated(allSamplesDE[!dupWithin,], fromLast=TRUE) | dupBetween[!dupWithin]
#     
#     dupBetween = rep(NA, nrow(df))
#     dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols])
#     dupBetween[!dupWithin] <- duplicated(df[!dupWithin,datacols], fromLast=TRUE) | dupBetween[!dupWithin]
#     
#     # ============= Replace NA's with previous non-NA value ==============
#     # This is why we sorted earlier - it was necessary to do this part efficiently
#     # Get indexes of non-NA's
#     goodIdx <- !is.na(dupBetween)
# 
#     # These are the non-NA values from x only
#     # Add a leading NA for later use when we index into this vector
#     goodVals <- c(NA, dupBetween[goodIdx])
# 
#     # Fill the indices of the output vector with the indices pulled from
#     # these offsets of goodVals. Add 1 to avoid indexing to zero.
#     fillIdx <- cumsum(goodIdx)+1
# 
#     # The original vector, now with gaps filled
#     dupBetween <- goodVals[fillIdx]
# 
#     # Undo the original sort
#     dupBetween[sortorder] <- dupBetween
# 
#     # Return the vector of which entries are duplicated across groups
#     return(dupBetween)
# }
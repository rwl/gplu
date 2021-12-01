use crate::internal::off;

/// maxmatch does maximum matching
///
/// maxmatch uses depth-first search to find an augmenting path from
/// each column node to get the maximum matching.
///
/// Alex Pothen and Chin-Ju Fan, Penn State University, 1988
/// last modifed: Alex Pothen July 1990
/// last bcs modifications:  John Lewis, Sept. 1990
///
/// input variables :
///
///    nrows -- number of row nodes in the graph.
///    ncols -- number of column nodes in the graph.
///    colstr, rowind -- adjacency structure of graph, stored by
///                      columns
///
/// output variables :
///
///    rowset -- describe the matching.
///              rowset (row) = col > 0 means column "col" is matched
///                                     to row "row"
///                           = 0       means "row" is an unmatched
///                                     node.
///    colset -- describe the matching.
///              colset (col) = row > 0 means row "row" is matched to
///                             column "col"
///                           = 0       means "col" is an unmatched
///                                     node.
pub fn maxmatch(
    nrows: usize,
    ncols: usize,
    colstr: &[usize],
    rowind: &[usize],
    prevcl: &mut [usize],
    prevrw: &mut [usize],
    marker: &mut [usize],
    tryrow: &mut [usize],
    nxtchp: &mut [isize],
) -> Result<(/*rowset*/ Vec<usize>, /*colset*/ Vec<usize>), String> {
    // Working variables :
    //
    //     prevrw (ncols) -- pointer toward the root of the depth-first
    //                       search from a column to a row.
    //     prevcl (ncols) -- pointer toward the root of the depth-first
    //                       search from a column to a column.
    //                       the pair (prevrw,prevcl) represent a
    //                       matched pair.
    //     marker (nrows) -- marker (row) <= the index of the root of the
    //                       current depth-first search.  row has been
    //                       visited in current pass when equality holds.
    //     tryrow (ncols) -- tryrow (col) is a pointer into rowind to
    //                       the next row to be explored from column col
    //                       in the depth-first search.
    //     nxtchp (ncols) -- nxtchp (col) is a pointer into rowind to the
    //                       next row to be explored from column col for
    //                       the cheap assignment.  set to -1 when
    //                       all rows have been considered for
    //                       cheap assignment
    let mut row: usize;
    let mut prow: usize;
    let mut pcol: usize;
    let mut nextrw: isize;
    let mut lastrw: usize;

    let mut rowset = vec![0; nrows];
    let mut colset = vec![0; ncols];

    for nodec in 1..=ncols {
        // Initialize node 'col' as the root of the path.
        let mut col = nodec;
        prevrw[col - off] = 0;
        prevcl[col - off] = 0;
        nxtchp[col - off] = colstr[col - off] as isize;

        // Main loop begins here. Each time through, try to find a
        // cheap assignment from node col.
        'l100: loop {
            nextrw = nxtchp[col - off];
            lastrw = colstr[col + 1 - off] - 1;

            'l400: loop {
                if nextrw > 0 {
                    for xrow in nextrw as usize..=lastrw {
                        row = rowind[xrow - off];
                        if rowset[row - off] == 0 {
                            break 'l400;
                        }
                    }

                    // Mark column when all adjacent rows have been
                    // considered for cheap assignment.
                    nxtchp[col - off] = -1;
                }

                // Each time through, take a step forward if possible, or
                // backtrack if not .  Quit when backtracking takes us back
                // to the beginning of the search.

                tryrow[col - off] = colstr[col - off];
                nextrw = tryrow[col - off] as isize;
                //lastrw = colstr [col+1-off] - 1;

                if lastrw >= nextrw as usize {
                    for xrow in nextrw as usize..=lastrw {
                        // next line inserted by Alex Pothen, July 1990
                        // ii  = xrow
                        row = rowind[xrow - off];
                        if marker[row - off] < nodec {
                            // Row is unvisited yet for this pass.
                            // Take a forward step.

                            tryrow[col - off] = xrow + 1;
                            marker[row - off] = nodec;
                            let nxtcol = rowset[row - off];

                            if nxtcol < 0 {
                                return Err(
                                    "maxmatch: search reached a forbidden column".to_string()
                                );
                            } else if nxtcol == col {
                                return Err("maxmatch: search followed a matching edge".to_string());
                            } else if nxtcol > 0 {
                                // The forward step led to a matched row
                                // try to extend augmenting path from
                                // the column matched by this row.

                                prevcl[nxtcol - off] = col;
                                prevrw[nxtcol - off] = row;
                                tryrow[nxtcol - off] = colstr[nxtcol - off];
                                col = nxtcol;
                                continue 'l100;
                            } else {
                                // Unmatched row
                                break 'l400;
                            }
                        }
                        //l300: continue
                    }
                }

                // No forward step -- backtrack.
                // If we backtrack all the way, the search is done

                col = prevcl[col - off];
                if col > 0 {
                    continue 'l100;
                } else {
                    break 'l100;
                }
                break;
            }

            // Update the matching by alternating the matching
            // edge backward toward the root.
            rowset[row - off] = col;
            prow = prevrw[col - off];
            pcol = prevcl[col - off];

            while pcol > 0 {
                if rowset[prow - off] != col {
                    return Err(format!("maxmatch: pointer toward root disagrees with matching. prevcl[{}]={} but colset[{}]={}", col, row, row, rowset[row-off]));
                }
                rowset[prow - off] = pcol;
                col = pcol;
                prow = prevrw[pcol - off];
                pcol = prevcl[pcol - off];
            }
            break;
        }
    }

    // Compute the matching from the view of column nodes.
    for row in 1..=nrows {
        let col = rowset[row - off];
        if col > 0 {
            colset[col - off] = row;
        }
    }

    Ok((rowset, colset))
}

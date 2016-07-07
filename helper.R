#'
#' EuPathDB gene information table GO term parser
#'
#' Note: EuPathDB currently includes some GO annotations corresponding to
#' obsolete terms. For example, the L. major gene LmjF.19.1390
#' (http://tritrypdb.org/tritrypdb/showRecord.do?name=GeneRecordClasses.GeneRecordClass&source_id=LmjF.19.1390&project_id=TriTrypDB)
#' includes the term "GO:0003702" on the website and txt data file. The
#' term has been deprecated, and does not have a category associated with it
#' on the website. These will not be included in the final database.
#'
#' @author Keith Hughitt
#'
#' @param filepath Location of TriTrypDB gene information table.
#' @param verbose  Whether or not to enable verbose output.
#' @return Returns a dataframe where each line includes a gene/GO terms pair
#'         along with some addition information about the GO term. Note that
#'         because each gene may have multiple GO terms, a single gene ID may
#'         appear on multiple lines.
#'
parse_go_terms = function (filepath) {
    if (file_ext(filepath) == 'gz') {
        fp = gzfile(filepath, open='rb')
    } else {
        fp = file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    N = 1e5
    gene_ids = c()
    go_rows = data.frame(GO=rep("", N),
                         ONTOLOGY=rep("", N), GO_TERM_NAME=rep("", N),
                         SOURCE=rep("", N), EVIDENCE=rep("", N),
                         stringsAsFactors=FALSE)

    # Counter to keep track of row number
    i = j = 1

    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id = .get_value(x)
            i = i + 1
        }

        # Gene Ontology terms
        else if (grepl("^GO:", x)) {
            gene_ids[j] = gene_id
            go_rows[j,] = c(head(unlist(strsplit(x, '\t')), 5))
            j = j + 1
        }
    }

    # get rid of unallocated rows
    go_rows = go_rows[1:j-1,]

    # drop unneeded columns
    go_rows = go_rows[,c('GO', 'EVIDENCE')]

    # add gene id column
    go_rows = cbind(GID=gene_ids, go_rows)

    # close file pointer
    close(fp)

    # Drop duplicate entries for different evidence codes and return result
    return(unique(go_rows))
}

#'
#' Queries a EuPathDB database for GO term annotations for a species.
#'
#' @author Keith Hughitt
#'
#' @param database URL of database to query (e.g. tritrypdb.org)
#' @organism Full organism name as it appears in the database 
#'  (e.g. "Leishmania major strain Friedlin")
#'
#' @returns A one-to-many dataframe mapping from gene ids to GO terms.
#'
retrieve_go_terms <- function (database, organism) {
    # URL-encode organism name
    organism <- utils::URLencode(organism)

    # Construct query URI
    entry_url <- 'webservices/GeneQuestions/GenesByTaxon.json'
    query_url <- sprintf("%s/%s?organism=%s&o-tables=GoTerms", database,
                         entry_url, organism)

    # Fetch JSON
    message(sprintf("Querying: %s", query_url))
    json <- RCurl::getURL(query_url)

    # Convert to R object
    result <- jsonlite::fromJSON(json)
    records <- result$response$recordset$records

    # dataframe to store result
    gene_go_mapping <- data.frame()

    # Iterate over genes and store GO terms
    for (i in 1:nrow(records)) {
        gene_id    <- records[i,'id']
        gene_table <- records[i,'tables'][[1]]

        # If no annotations found for gene, stop here
        if (length((gene_table$rows)[[1]]) == 0) {
            next 
        }

        # otherwise parse reuslt for gene
        table_fields <- gene_table$rows[[1]]$fields

        # Example entry in table_fields:
        # [[8]]
        #            name              value
        # 1         go_id         GO:0005524
        # 2      ontology Molecular Function
        # 3  go_term_name        ATP binding
        # 4        source           Interpro
        # 5 evidence_code                IEA
        # 6        is_not               <NA>

        # retrieve GO ids and evidence codes
        go_ids <- sapply(table_fields, function(x) { x$value[1] })
        evidence <- sapply(table_fields, function(x) { x$value[5] })

        # add to output dataframe
        gene_go_mapping <- rbind(gene_go_mapping, 
                                 cbind(GID=gene_id, GO=go_ids, EVIDENCE=evidence))
    }

    # Keep the strongest evidence code available for each annotation and
    # drop the rest.
    # http://geneontology.org/page/guide-go-evidence-codes

    # Separate into automatically inferred annotations (those with evidence
    # code IEA, "Inferred from Electronic Annotation", and all others.
    iea_annotations <- gene_go_mapping[gene_go_mapping$EVIDENCE == 'IEA',]                                                                                                 
    other_annotations <- gene_go_mapping[gene_go_mapping$EVIDENCE != 'IEA',]

    # Deduplicate other annotations, and add an IEA entry back in for
    # each gene/GO pair not supported by another type of evidence
    other_annotations <- other_annotations[!duplicated(other_annotations[,1:2]),]

    gene_go_mapping <- rbind(other_annotations, iea_annotations)

    return(gene_go_mapping)
}

#'
#' EuPathDB gene information table InterPro domain parser
#'
#' @author Keith Hughitt
#'
#' @param filepath Location of TriTrypDB gene information table.
#' @param verbose  Whether or not to enable verbose output.
#'
#' @return Returns a dataframe where each line includes a gene/domain pairs.
#'
parse_interpro_domains = function (filepath) {
    if (file_ext(filepath) == 'gz') {
        fp = gzfile(filepath, open='rb')
    } else {
        fp = file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    #N = 1e5
    #gene_ids = c()
    #interpro_rows = data.frame(GO=rep("", N),
    #                     ONTOLOGY=rep("", N), GO_TERM_NAME=rep("", N),
    #                     SOURCE=rep("", N), EVIDENCE=rep("", N),
    #                     stringsAsFactors=FALSE)

    # InterPro table columns
    cols = c('name', 'interpro_id', 'primary_id', 'secondary_id', 'description',
             'start_min', 'end_min', 'evalue')

    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id = .get_value(x)
        }

        # Parse InterPro table
        else if (grepl("TABLE: InterPro Domains", x)) {
            # Skip column header row
            trash = readLines(fp, n=1)

            # Continue reading until end of table
            raw_table = ""

            entry = readLines(fp, n=1)

            while(length(entry) != 0) {
                if (raw_table == "") {
                    raw_table = entry
                } else {
                    raw_table = paste(raw_table, entry, sep='\n')
                }
                entry = readLines(fp, n=1)
            }

            # If table length is greater than 0, read ino
            buffer = textConnection(raw_table)

            interpro_table = read.delim(buffer, header=FALSE, col.names=cols)

        }
    }


    # add gene id column
    go_rows = cbind(GID=gene_ids, go_rows)

    # close file pointer
    close(fp)

    # TODO: Determine source of non-unique rows in the dataframe
    return(unique(go_rows))
}


#'
#' EuPathDB gene information table gene type parser
#'
#' @author Keith Hughitt
#'
#' @param filepath Location of TriTrypDB gene information table.
#' @return Returns a dataframe mapping gene ids to gene types
#'
parse_gene_types = function (filepath) {
    if (file_ext(filepath) == 'gz') {
        fp = gzfile(filepath, open='rb')
    } else {
        fp = file(filepath, open='r')
    }

    # Create empty vector to store dataframe rows
    N = 1e5
    gene_ids   = c()
    gene_types = c()

    # Counter to keep track of row number
    i = 1

    # Iterate through lines in file
    while (length(x <- readLines(fp, n=1, warn=FALSE)) > 0) {
        # Gene ID
        if(grepl("^Gene ID", x)) {
            gene_id = .get_value(x)
        }
        # Gene type
        else if (grepl("^Gene Type:", x)) {
            gene_ids[i]   = gene_id
            gene_types[i] = .get_value(x)
            i = i + 1
        }
    }
    # close file pointer
    close(fp)

    return(data.frame(GID=gene_ids, TYPE=gene_types))
}

#
# kegg_to_genedb
#
# Takes a list of KEGG gene identifiers and returns a list of GeneDB
# ids corresponding to those genes. 
#
kegg_to_genedb = function(kegg_ids, gene_mapping) {
    # query gene ids 10 at a time (max allowed)
    result = c()

    for (x in split(kegg_ids, ceiling(seq_along(kegg_ids) / 10))) {
        query = keggGet(x)
        for (item in query) {
            dblinks = item$DBLINKS
            genedb_id = dblinks[grepl('GeneDB', dblinks)]
            if (length(genedb_id) > 0) {
                # get old-style t. cruzi identifier
                old_id = substring(genedb_id, 9)

                # if possible, map to new id and add to results
                if (!is.null(gene_mapping[[old_id]])) {
                    result = append(result, gene_mapping[[old_id]])
                }
            }
    
        }
    }
    return(result)
}

#
# Parses a key: value string and returns the value
#
.get_value = function(x) {
    return(gsub(" ","", tail(unlist(strsplit(x, ': ')), n=1), fixed=TRUE))
}

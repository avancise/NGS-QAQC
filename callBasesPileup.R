callBasesPileup <- function(plp.fname, min.cov = 7, max.prop.only = 11,
                            prop.thresh = 0.8, binom.pr.thresh = 0.8,
                            min.prop = 0.2) {

  stopifnot(require(dplyr))
  stopifnot(require(tidyr))
  stopifnot(require(readr))
  stopifnot(require(ape))

  # read pileup file
  plp <- read_csv(
    plp.fname,
    col_types = cols("c", "c", "i", "i", "c", "c", "l", "d", "i", "d", "d", "d", "d", "d", "d"),
    progress = FALSE
  ) %>%
    rename(position = ref.pos) %>%
    as.data.frame

  base.cols <- c("A", "C", "G", "T", "N", "-")
  max.pos <- max(plp$position)

  # extract base frequencies from proportions
  base.freq <- round(plp[, base.cols] * plp$coverage, 0)
  base.freq$position <- plp$position
  base.freq$id <- plp$id
  base.freq$coverage <- plp$coverage

  # get total proportions of each base at each position
  position.base.props <- base.freq %>%
    select(-id, -coverage) %>%
    group_by(position) %>%
    summarize_each(funs(sum)) %>%
    mutate(total.reads = rowSums(.[, -1])) %>%
    mutate_each(funs(. / total.reads), -position, -total.reads) %>%
    gather(base, pool.prop, -position, -total.reads)

  # evaluate if frequency and binomial probability are good for base calling
  freq.prob.check <- base.freq %>%
    gather(base, freq, -id, -position, -coverage) %>%
    left_join(position.base.props, by = c("position", "base")) %>%
    filter(freq > 0) %>%
    mutate(
      read.prop = freq / coverage,
      pr.base = pbinom(freq, coverage, pool.prop),
      freq.good = coverage >= min.cov & read.prop == 1 | coverage >= max.prop.only & read.prop >= prop.thresh,
      prob.good = coverage >= max.prop.only &
                  pr.base >= binom.pr.thresh &
                  read.prop >= min.prop
    )

  # call bases at each site for each individual
  baseCallFunc <- function(base, freq.good, pr.base, prob.good) {
    base.call <- base[freq.good]
    if(length(base.call) == 1) return(base.call)
    max.pr <- which.max(pr.base)
    if(prob.good[max.pr]) return(base[max.pr])
    "N"
  }
  base.calls <- freq.prob.check %>%
    group_by(id, position) %>%
    summarize(base = baseCallFunc(base, freq.good, pr.base, prob.good))

  # translate to DNAbin format full sequences
  dna.seqs <- as.matrix(as.DNAbin(do.call(rbind, by(base.calls, base.calls$id, function(x) {
    x.seq <- rep("N", max.pos)
    x.seq[x$position] <- x$base
    x.seq
  }))))

  # identify positions with variable reads for each individual
  plp$variable.bases.in.reads <- apply(plp[, base.cols], 1, function(x) sum(x == 1) == 0)
  plp <- plp[, c("position", "id", "variable.bases.in.reads")]

  # identify positions with variable calls
  vs <- data.frame(position = 1:max.pos)
  vs$variable.calls.at.position <- apply(as.character(dna.seqs), 2, function(x) {
    length(unique(x)) > 1
  })

  # add column of variable info to base frequency data.frame
  base.freq <- base.freq[, c("position", "id", "coverage", base.cols)] %>%
    left_join(base.calls, by = c("position", "id")) %>%
    rename(called.base = base) %>%
    left_join(plp, by = c("position", "id")) %>%
    left_join(vs, by = "position") %>%
    arrange(position, id)

  # write fasta file
  fname <- gsub(".csv", ".fasta", plp.fname)
  write.dna(
    dna.seqs, fname, format = "fasta", nbcol = -1,
    colsep = "", indent = 0, blocksep = 0
  )

  # return data.frames of metrics for positions and ids, fasta filename, and DNAbin sequences
  list(
    base.freq.prob.narrow = freq.prob.check,
    base.freq.wide = base.freq,
    fasta.fname = fname,
    dna.seqs = dna.seqs
  )
}
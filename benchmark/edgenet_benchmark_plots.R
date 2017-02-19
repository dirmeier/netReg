library(stringr)
library(ggplot2)
library(microbenchmark)

bench.fls <- list.files("../../../results/", full.names=T)
bench.df <- c()
ben.l <- list()
for (fl in bench.fls)
{
  p <- as.numeric(stringr::str_match(fl, "benchmark_([[:digit:]]+)_.*")[[2]])
  r <- readRDS(fl)
  bench.df <- rbind(bench.df, cbind(p,r ))
  ben.l[[paste(p)]] <- r
}
bench.df

bench.df$expr <- factor(bench.df$expr)

ggplot2::ggplot(bench.df) +
  ggplot2::geom_boxplot(aes(as.factor(p), time, colour=as.factor(expr))) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))  +
  ggplot2::theme_bw() 
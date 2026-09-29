#' Plot IREA response or polarization results
#'
#' @param result An `irea_result`, a Seurat object containing
#'   `object@tools$IREA`, or a named list of results for dotplots and heatmaps.
#' @param type One of `"compass"`, `"radar"`, `"dotplot"`, or `"heatmap"`.
#' @param fdr_cutoff Significance threshold used for plot colouring.
#' @return An editable ggplot object.
#' @export
IREAPlot <- function(result, type=c("compass","radar","dotplot","heatmap"),
                     fdr_cutoff=0.05) {
  type <- match.arg(type)
  if (inherits(result,"Seurat")) result <- result@tools$IREA
  multi <- is.list(result) && !inherits(result,"irea_result")
  if (multi) {
    if (type %in% c("compass","radar"))
      stop("Compass and radar plots require one IREA result.",call.=FALSE)
    if (!length(result) || !all(vapply(result,inherits,logical(1),"irea_result")))
      stop("Every list element must be an IREA result.",call.=FALSE)
    if (length(unique(vapply(result,function(x)x$parameters$analysis,character(1)))) != 1L)
      stop("Comparison results must use the same analysis.",call.=FALSE)
    groups <- names(result)
    if (is.null(groups) || any(!nzchar(groups))) groups <- paste0("Result ",seq_along(result))
    dat <- do.call(rbind,Map(function(x,name) transform(x$table,group=name),result,groups))
  } else {
    if (!inherits(result,"irea_result")) stop("result must be an IREA result.",call.=FALSE)
    dat <- result$table
    dat$group <- "Result"
  }
  if (!is.numeric(fdr_cutoff) || length(fdr_cutoff)!=1L || is.na(fdr_cutoff) ||
      fdr_cutoff <= 0 || fdr_cutoff > 1) stop("Invalid fdr_cutoff.",call.=FALSE)
  dat$significant <- !is.na(dat$fdr) & dat$fdr < fdr_cutoff
  if (type == "compass") {
    if (result$parameters$analysis != "cytokine_response")
      stop("Compass plot requires a cytokine response result.",call.=FALSE)
    dat$direction <- ifelse(!dat$significant,"Not significant",
                            ifelse(dat$effect >= 0,"Positive","Negative"))
    dat$term <- factor(dat$term,levels=dat$term[order(dat$effect)])
    return(ggplot2::ggplot(dat,ggplot2::aes(x=.data$term,y=abs(.data$effect),fill=.data$direction)) +
      ggplot2::geom_col(width=1) + ggplot2::coord_polar() +
      ggplot2::scale_fill_manual(values=c(Positive="#BC3C29",Negative="#2878B5",
                                          `Not significant`="#B8B8B8")) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(size=5)) +
      ggplot2::labs(x=NULL,y="Absolute enrichment effect",fill="Response"))
  }
  if (type == "radar") {
    if (result$parameters$analysis != "cell_polarization")
      stop("Radar plot requires a cell polarization result.",call.=FALSE)
    theta <- pi/2 - 2*pi*(seq_len(nrow(dat))-1)/nrow(dat)
    dat$x <- dat$radar_score*cos(theta)
    dat$y <- dat$radar_score*sin(theta)
    axes <- data.frame(x=cos(theta),y=sin(theta),term=dat$term)
    circles <- do.call(rbind,lapply(c(0.25,0.5,0.75,1),function(radius) {
      t <- seq(0,2*pi,length.out=181)
      data.frame(x=radius*cos(t),y=radius*sin(t),radius=factor(radius))
    }))
    return(ggplot2::ggplot() +
      ggplot2::geom_path(data=circles,ggplot2::aes(x=.data$x,y=.data$y,group=.data$radius),color="#D9D9D9") +
      ggplot2::geom_segment(data=axes,ggplot2::aes(x=0,y=0,xend=.data$x,yend=.data$y),color="#D9D9D9") +
      ggplot2::geom_polygon(data=dat,ggplot2::aes(x=.data$x,y=.data$y),fill="#2F5597",
                            alpha=0.25,color="#2F5597") +
      ggplot2::geom_point(data=dat,ggplot2::aes(x=.data$x,y=.data$y),color="#2F5597") +
      ggplot2::geom_text(data=axes,ggplot2::aes(x=1.13*.data$x,y=1.13*.data$y,label=.data$term),size=3) +
      ggplot2::coord_fixed(xlim=c(-1.25,1.25),ylim=c(-1.25,1.25)) +
      ggplot2::theme_void() +
      ggplot2::labs(title="Cell polarization",subtitle="Normalized score (0 to 1)"))
  }
  dat$term <- factor(dat$term,levels=unique(dat$term[order(dat$effect)]))
  if (type == "dotplot") {
    return(ggplot2::ggplot(dat,ggplot2::aes(x=.data$group,y=.data$term)) +
      ggplot2::geom_point(ggplot2::aes(size=-log10(pmax(.data$fdr,.Machine$double.xmin)),
                                        color=.data$effect)) +
      ggplot2::scale_color_gradient2() +
      ggplot2::labs(x=NULL,y=NULL,size="-log10(FDR)",color="Effect"))
  }
  ggplot2::ggplot(dat,ggplot2::aes(x=.data$group,y=.data$term,fill=.data$effect)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient2() +
    ggplot2::labs(x=NULL,y=NULL,fill="Enrichment effect")
}

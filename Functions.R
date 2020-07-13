# ---------------------------------------------------------------------------- #
plot_tSNE_Annot = function(
  sce       = NULL,
  var       = 'status',
  col       = NULL,
  title     = NULL,
  proj_name = 'UMAP',
  outpath   = './',
  width     = 4,
  height    = 3,
  size      = 0.5,
  zlim      = NULL,
  midpoint  = NULL
){
  
  meta = colData(sce) %>%
    as.data.frame() %>%
    mutate(score = .[[var]]) %>%
    cbind(
      reducedDim(sce, proj_name) %>% 
        magrittr::set_colnames(str_replace(colnames(.),'sub',''))
    ) 
    
  proj_name = str_replace(proj_name,'_SUB','')
  if(!all(is.na(meta$score))){
    g = meta %>%  {
      ggplot(data = .,aes_string(paste(proj_name, '1',sep='_'), paste(proj_name, '2',sep='_'))) +
        geom_point(data=., size=size, color='lightgray') +
        theme(
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)
        ) +
        xlab(paste(proj_name,'1', sep='_')) +
        ylab(paste(proj_name,'1', sep='_'))
    }
    
    if(!is.null(title))
      g = g + ggtitle(title)
    
    # Sets the colors
    if(is.numeric(meta$score) | is.integer(meta$score)){
      if(is.null(zlim))
        zlim = range(meta$score)
      
      if(is.null(midpoint))
        midpoint = median(meta$score)
      
      g = g +
        geom_point(data = subset(meta, score != 0), aes(color=score), size=size) +
        scale_color_gradient2(
          high     = globals$params$colors$high,
          mid      = globals$params$colors$mid,
          low      = globals$params$colors$low,
          limits   = zlim,
          midpoint = midpoint
        )
    }else{
      if(is.null(col))
        col = ggplotColors(length(unique(meta[[var]])))
      
      g = g +
        geom_point(aes(color=score), size=size) +
        scale_color_manual(values = col)
    }
    
    print(g)
  }
  
}


# ---------------------------------------------------------------------------- #
plot_tSNE_Gene = function(
  sce       = NULL,
  gname     = NULL,
  type      = 'logcounts',
  proj_name = 'UMAP',
  col       = NULL,
  width     = 4,
  height    = 3,
  size      = 0.5,
  plot      = TRUE
){
  
  meta = colData(sce) %>%
    as.data.frame() %>%
    cbind(reducedDim(sce, proj_name))
  
  gid = rowData(reg_sce)$gene_name == gname
  g = meta %>%
    mutate(score = as.vector(assays(reg_sce)[[type]][gid,])) %>%  {
      ggplot(data = .,aes_string(paste(proj_name, '1',sep='_'), paste(proj_name, '2',sep='_'))) +
        geom_point(size=size, color='lightgray') +
        geom_point(data = subset(., score != 0), aes(color=score), size=size) +
        scale_color_gradient2(
          high = 'firebrick',
          mid  = 'lightgray',
          low  = 'cornflowerblue'
        ) +
        theme(
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)
        ) +
        xlab(paste(proj_name, '1')) +
        ylab(paste(proj_name, '2')) +
        ggtitle(gname)
    }
  
  if(plot){
    print(g)
  }else{
    return(g)
  }
}

# ---------------------------------------------------------------------------- #
plot_Violin_Gene = function(
  sce       = NULL,
  gname     = NULL,
  cell_type = NULL,
  type      = 'logcounts',
  col       = NULL,
  width     = 4,
  height    = 3,
  size      = 0.5,
  plot      = TRUE,
  colors    = NULL
){
  
  meta = colData(sce) %>%
    as.data.frame() 
  
  if(is.null(colors)){
    colors = ggplotColors(length(unique(meta[[cell_type]])))
    colors = setNames(colors, unique(meta[[cell_type]]))
  }
  
  gid = rowData(sce)$gene_name == gname
  g = meta %>%
    mutate(cell_type = .[[cell_type]]) %>%
    mutate(index = as.numeric(cell_type)) %>%
    mutate(cell_type = as.character(cell_type)) %>%
    mutate(score = as.vector(assays(sce)[[type]][gid,.$cell_id])) %>%  {
      ggplot(data = .,aes(reorder(cell_type,index, .fun='median'), score, fill=cell_type, color=cell_type)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values  = colors) +
        scale_color_manual(values = colors) +
        theme(
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(size=16),
          axis.text.y  = element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5)
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab('Cell Type') +
        ylab('log2 Normalized Counts') +
        ggtitle(gname)
    }
  
  if(plot){
    print(g)
  }else{
    return(g)
  }
}


# ---------------------------------------------------------------------------- #
# given a number return n colors
ggplotColors =  function(g){
  d <- 360/g
  h <- cumsum(c(15, rep(d,g - 1)))
  hcl(h = h, c = 100, l = 65)
}
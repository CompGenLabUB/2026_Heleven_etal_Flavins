make_pw_plots <- function(DATA, gna, gnb, PRFX='') {

    # Specify column names for the plots as symbols using sym(), to make rlang "!!" work
    GNA <- sym(gna);
    GNB <- sym(gnb);
    #
    DP <- DATA[  (DATA$cluster_label %in% cell.ord), ];
    DN <- DATA[ !(DATA$cluster_label %in% cell.ord), ];
    #
    DD <- rbind(DN,
                DP[ order(DP[ , gna] + DP[ , gnb], DP[ , gna], DP[, gnb], DP$cell_id), ]);
    #
    # BY CELL CLUSTER
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                          color = cluster_label)) +
             geom_point() +
             theme_bw() +
             xlab(paste0(GNA," [",GN.labels[[GNA]],"]")) +
             ylab(paste0(GNB," [",GN.labels[[GNB]],"]")) +
             labs(color="Cell\nCluster") +
             scale_color_manual(breaks=cell.ord,
                                values=cell.colors,
                                na.value=cell.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(plot=p,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_bycellcluster.png"),
           height=10, width=11, unit="in", dpi=600);
    q <- p + facet_wrap( . ~ region_label);
    print(q);
    ggsave(plot=q,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_bycellcluster_facetedbyregion.png"),
           height=15, width=22, unit="in", dpi=600);
    k <- p + facet_wrap( . ~ cluster_shrt);
    print(k);
    ggsave(plot=k,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_bycellcluster_facetedbycells.png"),
           height=15, width=16, unit="in", dpi=600);
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                        color = cluster_label)) +
             geom_jitter(width = 0.15, height=0.15) +
             theme_bw() +
             xlab(paste0(GNA," [",GN.labels[[GNA]],"]")) +
             ylab(paste0(GNB," [",GN.labels[[GNB]],"]")) +
             labs(color="Cell\nCluster") +
             scale_color_manual(breaks=cell.ord,
                                values=cell.colors,
                                na.value=cell.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(plot=p,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_bycellcluster.png"),
           height=10, width=11, unit="in", dpi=600);
    q <- p + facet_wrap( . ~ region_label);
    print(q);
    ggsave(plot=q,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_bycellcluster_facetedbyregion.png"),
           height=15, width=22, unit="in", dpi=600);
    k <- p + facet_wrap( . ~ cluster_shrt);
    print(k);
    ggsave(plot=k,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_bycellcluster_facetedbycells.png"),
           height=15, width=16, unit="in", dpi=600);
    #
    # BY CELL TYPE
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                          color = cluster_shrt)) +
             geom_point() +
             theme_bw() +
             xlab(paste0(GNA," [",GN.labels[[GNA]],"]")) +
             ylab(paste0(GNB," [",GN.labels[[GNB]],"]")) +
             labs(color="Cell\nType") +
             scale_color_manual(breaks=CELL.ord,
                                values=CELL.colors,
                                na.value=CELL.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(plot=p,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_bycelltype.png"),
           height=10, width=11, unit="in", dpi=600);
    q <- p + facet_wrap( . ~ region_label);
    print(q);
    ggsave(plot=q,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_bycelltype_facetedbyregion.png"),
           height=15, width=22, unit="in", dpi=600);
    k <- p + facet_wrap( . ~ cluster_shrt);
    print(k);
    ggsave(plot=k,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_bycelltype_facetedbycells.png"),
           height=15, width=16, unit="in", dpi=600);
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                        color = cluster_shrt)) +
             geom_jitter(width = 0.15, height=0.15) +
             theme_bw() +
             xlab(paste0(GNA," [",GN.labels[[GNA]],"]")) +
             ylab(paste0(GNB," [",GN.labels[[GNB]],"]")) +
             labs(color="Cell\nType") +
             scale_color_manual(breaks=CELL.ord,
                                values=CELL.colors,
                                na.value=CELL.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(plot=p,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_bycelltype.png"),
           height=10, width=11, unit="in", dpi=600);
    q <- p + facet_wrap( . ~ region_label);
    print(q);
    ggsave(plot=q,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_bycelltype_facetedbyregion.png"),
           height=15, width=22, unit="in", dpi=600);
    k <- p + facet_wrap( . ~ cluster_shrt);
    print(k);
    ggsave(plot=k,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_bycelltype_facetedbycells.png"),
           height=15, width=16, unit="in", dpi=600);
    #
    # BY REGION
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                          color = region_label)) +
             geom_point() +
             theme_bw() +
             xlab(paste0(GNA," [",GN.labels[[GNA]],"]")) +
             ylab(paste0(GNB," [",GN.labels[[GNB]],"]")) +
             labs(color="Sample\nRegion") +
             scale_color_manual(breaks=REG.ord,
                                values=REG.colors,
                                na.value=REG.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(plot=p,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_byregion.png"),
           height=10, width=11, unit="in", dpi=600);
    q <- p + facet_wrap( . ~ region_label);
    print(q);
    ggsave(plot=q,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_byregion_facetedbyregion.png"),
           height=15, width=22, unit="in", dpi=600);
    k <- p + facet_wrap( . ~ cluster_shrt);
    print(k);
    ggsave(plot=k,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionplot_byregion_facetedbycells.png"),
           height=15, width=16, unit="in", dpi=600);
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                        color = region_label)) +
             geom_jitter(width = 0.15, height=0.15) +
             theme_bw() +
             xlab(paste0(GNA," [",GN.labels[[GNA]],"]")) +
             ylab(paste0(GNB," [",GN.labels[[GNB]],"]")) +
             labs(color="Sample\nRegion") +
             scale_color_manual(breaks=REG.ord,
                                values=REG.colors,
                                na.value=REG.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(plot=p,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_byregion.png"),
           height=10, width=11, unit="in", dpi=600);
    q <- p + facet_wrap( . ~ region_label);
    print(q);
    ggsave(plot=q,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_byregion_facetedbyregion.png"),
           height=15, width=22, unit="in", dpi=600);
    k <- p + facet_wrap( . ~ cluster_shrt);
    print(k);
    ggsave(plot=k,
           paste0(PRFX,GNA,"-x-",GNB,"_coexpressionjitterplot_byregion_facetedbycells.png"),
           height=15, width=16, unit="in", dpi=600);

} # make_pw_plots

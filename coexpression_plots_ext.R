make_pw_plots2 <- function(DATA, gna, gnb,
                           cell.ord, cell.colors,
                           CELL.ord, CELL.colors,
                           REG.ord, REG.colors,
                           thyset.labels, PFX='', SFX='', DEV="png", RG.flg=TRUE) {

    # Specify column names for the plots as symbols using sym(), to make rlang "!!" work
    GNA <- sym(gna);
    GNB <- sym(gnb);
    GNA.sym <- thyset.labels[ gna, "symbol" ];
    GNB.sym <- thyset.labels[ gnb, "symbol" ]
    GNA.s <- gsub(" / ",":",GNA.sym);
    GNB.s <- gsub(" / ",":",GNB.sym);
    #
    DP <- DATA[  (DATA$cluster_label2 %in% cell.ord), ];
    DN <- DATA[ !(DATA$cluster_label2 %in% cell.ord), ];
    #
    DD <- rbind(DN,
                DP[ order(DP[ , gna] + DP[ , gnb], DP[ , gna], DP[, gnb], DP$cell_id), ]);
    #
    # BY CELL CLUSTER
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                          color = cluster_label2)) +
             geom_point() +
             theme_bw() +
             xlab(paste0(GNA," [",GNA.sym,"]")) +
             ylab(paste0(GNB," [",GNB.sym,"]")) +
             labs(color="Cell\nCluster") +
             guides(color = guide_legend(ncol = 1)) +
             scale_color_manual(breaks=cell.ord,
                                values=cell.colors,
                                na.value=cell.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_bycellcluster",SFX,".",DEV),
           plot=p, device=DEV,
           height=10, width=11, unit="in", dpi=600);
    if (RG.flg) {
        q <- p + facet_wrap( . ~ region_label);
        print(q);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_bycellcluster_facetedbyregion",SFX,".",DEV),
               plot=q, device=DEV,
               height=15, width=22, unit="in", dpi=600);
    };
    k <- p + facet_wrap( . ~ label_short);
    print(k);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_bycellcluster_facetedbycells",SFX,".",DEV),
           plot=k, device=DEV,
           height=15, width=16, unit="in", dpi=600);
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                        color = cluster_label2)) +
             geom_jitter(width = 0.15, height=0.15) +
             theme_bw() +
             xlab(paste0(GNA," [",GNA.sym,"]")) +
             ylab(paste0(GNB," [",GNB.sym,"]")) +
             labs(color="Cell\nCluster") +
             guides(color = guide_legend(ncol = 1)) +
             scale_color_manual(breaks=cell.ord,
                                values=cell.colors,
                                na.value=cell.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_bycellcluster",SFX,".",DEV),
           plot=p, device=DEV,
           height=10, width=11, unit="in", dpi=600);
    if (RG.flg) {
        q <- p + facet_wrap( . ~ region_label);
        print(q);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_bycellcluster_facetedbyregion",SFX,".",DEV),
               plot=q, device=DEV,
               height=15, width=22, unit="in", dpi=600);
    };
    k <- p + facet_wrap( . ~ label_short);
    print(k);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_bycellcluster_facetedbycells",SFX,".",DEV),
           plot=k, device=DEV,
           height=15, width=16, unit="in", dpi=600);
    #
    # BY CELL TYPE
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                          color = label_short)) +
             geom_point() +
             theme_bw() +
             xlab(paste0(GNA," [",GNA.sym,"]")) +
             ylab(paste0(GNB," [",GNB.sym,"]")) +
             labs(color="Cell\nType") +
             guides(color = guide_legend(ncol = 1)) +
             scale_color_manual(breaks=CELL.ord,
                                values=CELL.colors,
                                na.value=CELL.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_bycelltype",SFX,".",DEV),
           plot=p, device=DEV,
           height=10, width=11, unit="in", dpi=600);
    if (RG.flg) {
        q <- p + facet_wrap( . ~ region_label);
        print(q);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_bycelltype_facetedbyregion",SFX,".",DEV),
               plot=q, device=DEV,
               height=15, width=22, unit="in", dpi=600);
    };
    k <- p + facet_wrap( . ~ label_short);
    print(k);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_bycelltype_facetedbycells",SFX,".",DEV),
           plot=k, device=DEV,
           height=15, width=16, unit="in", dpi=600);
    #
    p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                        color = label_short)) +
             geom_jitter(width = 0.15, height=0.15) +
             theme_bw() +
             xlab(paste0(GNA," [",GNA.sym,"]")) +
             ylab(paste0(GNB," [",GNB.sym,"]")) +
             labs(color="Cell\nType") +
             guides(color = guide_legend(ncol = 1)) +
             scale_color_manual(breaks=CELL.ord,
                                values=CELL.colors,
                                na.value=CELL.colors[["Other"]]) +
             theme(legend.title.align=0.5);
    print(p);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_bycelltype",SFX,".",DEV),
           plot=p, device=DEV,
           height=10, width=11, unit="in", dpi=600);
    if (RG.flg) {
        q <- p + facet_wrap( . ~ region_label);
        print(q);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_bycelltype_facetedbyregion",SFX,".",DEV),
               plot=q, device=DEV,
               height=15, width=22, unit="in", dpi=600);
    };
    k <- p + facet_wrap( . ~ label_short);
    print(k);
    ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_bycelltype_facetedbycells",SFX,".",DEV),
           plot=k, device=DEV,
           height=15, width=16, unit="in", dpi=600);
    #
    # BY REGION
    #
    if (RG.flg) {
        p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                              color = region_label)) +
                 geom_point() +
                 theme_bw() +
                 xlab(paste0(GNA," [",GNA.sym,"]")) +
                 ylab(paste0(GNB," [",GNB.sym,"]")) +
                 labs(color="Sample\nRegion") +
                 guides(color = guide_legend(ncol = 1)) +
                 scale_color_manual(breaks=REG.ord,
                                    values=REG.colors,
                                    na.value=REG.colors[["Other"]]) +
                 theme(legend.title.align=0.5);
        print(p);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_byregion",SFX,".",DEV),
               plot=p, device=DEV,
               height=10, width=11, unit="in", dpi=600);
        # if (RG.flg) {
            q <- p + facet_wrap( . ~ region_label);
            print(q);
            ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_byregion_facetedbyregion",SFX,".",DEV),
                   plot=q, device=DEV,
                   height=15, width=22, unit="in", dpi=600);
        # };
        k <- p + facet_wrap( . ~ label_short);
        print(k);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionplot_byregion_facetedbycells",SFX,".",DEV),
               plot=k, device=DEV,
               height=15, width=16, unit="in", dpi=600);
        #
        p <- ggplot(DD, aes(x = !!GNA, y = !!GNB,
                            color = region_label)) +
                 geom_jitter(width = 0.15, height=0.15) +
                 theme_bw() +
                 xlab(paste0(GNA," [",GNA.sym,"]")) +
                 ylab(paste0(GNB," [",GNB.sym,"]")) +
                 labs(color="Sample\nRegion") +
                 guides(color = guide_legend(ncol = 1)) +
                 scale_color_manual(breaks=REG.ord,
                                    values=REG.colors,
                                    na.value=REG.colors[["Other"]]) +
                 theme(legend.title.align=0.5);
        print(p);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_byregion",SFX,".",DEV),
               plot=p, device=DEV,
               height=10, width=11, unit="in", dpi=600);
        # if (RG.flg) {
            q <- p + facet_wrap( . ~ region_label);
            print(q);
            ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_byregion_facetedbyregion",SFX,".",DEV),
                   plot=q, device=DEV,
                   height=15, width=22, unit="in", dpi=600);
        # };
        k <- p + facet_wrap( . ~ label_short);
        print(k);
        ggsave(paste0(PFX,GNA,":",GNA.s,"-x-",GNB,":",GNB.s,"_coexpressionjitterplot_byregion_facetedbycells",SFX,".",DEV),
               plot=k, device=DEV,
               height=15, width=16, unit="in", dpi=600);
    };

} # make_pw_plots

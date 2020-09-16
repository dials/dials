from __future__ import absolute_import, division, print_function

import logging

# modified version of the ab_cluster function so we can access the scipy dendrogram object
from xfel.clustering.cluster import Cluster

logger = logging.getLogger(__name__)


class UnitCellCluster(Cluster):
    def ab_cluster(
        self,
        threshold=10000,
        method="distance",
        linkage_method="single",
        log=False,
        ax=None,
        write_file_lists=True,
        schnell=False,
        doplot=True,
        labels="default",
    ):
        """
        Hierarchical clustering using the unit cell dimentions.

        :param threshold: the threshold to use for prunning the tree into clusters.
        :param method: which clustering method from scipy to use when creating the tree (see scipy.cluster.hierarchy)
        :param linkage_method: which linkage method from scipy to use when creating the linkages. x (see scipy.cluster.hierarchy)
        :param log: if True, use log scale on y axis.
        :param ax: if a matplotlib axes object is provided, plot to this.
                   Otherwise, create a new axes object and display on screen.
        :param write_file_lists: if True, write out the files that make up each cluster.
        :param schnell: if True, use simple euclidian distance, otherwise, use Andrews-Bernstein
                        distance from Andrews & Bernstein J Appl Cryst 47:346 (2014) on the Niggli cells.
        :param doplot: Boolean flag for if the plotting should be done at all.
                       Runs faster if switched off.
        :param labels: 'default' will not display any labels for more than 100 images, but will display
                       file names for fewer. This can be manually overidden with a boolean flag.
        :return: A list of Clusters ordered by largest Cluster to smallest

        .. note::
          Use 'schnell' option with caution, since it can cause strange behaviour
          around symmetry boundaries.
        """

        import numpy as np

        from cctbx.uctbx.determine_unit_cell import NCDist
        from xfel.clustering.singleframe import SingleFrame

        logger.info("Hierarchical clustering of unit cells")
        import scipy.cluster.hierarchy as hcluster
        import scipy.spatial.distance as dist

        # 1. Create a numpy array of G6 cells
        g6_cells = np.array([SingleFrame.make_g6(image.uc) for image in self.members])

        # 2. Do hierarchichal clustering, using the find_distance method above.
        if schnell:
            logger.info("Using Euclidean distance")
            metric = "euclidean"
        else:
            logger.info(
                "Using Andrews-Bernstein distance from Andrews & Bernstein "
                "J Appl Cryst 47:346 (2014)"
            )
            metric = NCDist
        pair_distances = dist.pdist(g6_cells, metric=metric)
        if len(pair_distances) > 0:
            logger.info("Distances have been calculated")
            this_linkage = hcluster.linkage(
                pair_distances, method=linkage_method, metric=metric
            )
            cluster_ids = hcluster.fcluster(this_linkage, threshold, criterion=method)
            logger.debug("Clusters have been calculated")
        else:
            logger.debug("No distances were calculated. Aborting clustering.")
            return [], None

        # 3. Create an array of sub-cluster objects from the clustering
        sub_clusters = []
        for cluster in range(max(cluster_ids)):
            info_string = (
                "Made using ab_cluster with t={}," " {} method, and {} linkage"
            ).format(threshold, method, linkage_method)
            sub_clusters.append(
                self.make_sub_cluster(
                    [
                        self.members[i]
                        for i in range(len(self.members))
                        if cluster_ids[i] == cluster + 1
                    ],
                    "cluster_{}".format(cluster + 1),
                    info_string,
                )
            )

        sub_clusters = sorted(sub_clusters, key=lambda x: len(x.members))
        # Rename to order by size
        for num, cluster in enumerate(sub_clusters):
            cluster.cname = "cluster_{}".format(num + 1)

        # 3.5 optionally write out the clusters to files.
        if write_file_lists:
            for cluster in sub_clusters:
                if len(cluster.members) > 1:
                    cluster.dump_file_list(out_file_name="{}.lst".format(cluster.cname))

        if labels is True:
            labels = [image.name for image in self.members]
        elif labels is False:
            labels = ["" for _ in self.members]
        elif labels == "default":
            if len(self.members) > 100:
                labels = ["" for _ in self.members]
            else:
                labels = [image.name for image in self.members]
        else:
            labels = [getattr(v, labels, "") for v in self.members]

        if doplot:
            import matplotlib.pyplot as plt

            # 4. Plot a dendogram to the axes if no axis is passed, otherwise just
            #    return the axes object
            if ax is None:
                fig = plt.figure("Distance Dendogram")
                ax = fig.gca()
                direct_visualisation = True
            else:
                direct_visualisation = False

        dendrogram = hcluster.dendrogram(
            this_linkage,
            labels=labels,
            p=200,
            truncate_mode="lastp",  # show only the last p merged clusters
            leaf_font_size=8,
            leaf_rotation=90.0,
            color_threshold=threshold,
            ax=ax,
            no_plot=not doplot,
        )

        if doplot:
            if log:
                ax.set_yscale("symlog", linthreshx=(-1, 1))
            else:
                ax.set_ylim(-ax.get_ylim()[1] / 100, ax.get_ylim()[1])

            if direct_visualisation:
                fig.savefig("{}_dendogram.pdf".format(self.cname))
                plt.show()

        return sub_clusters, dendrogram, ax

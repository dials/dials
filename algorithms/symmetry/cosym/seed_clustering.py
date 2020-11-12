"""Seed clustering method for cosym analysis."""
from __future__ import absolute_import, division, print_function

import logging

logger = logging.getLogger(__name__)

import math

import numpy as np
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy
from sklearn import metrics
from sklearn.neighbors import NearestNeighbors

from libtbx import Auto
from libtbx.utils import Sorry
from scitbx.array_family import flex


class seed_clustering(object):
    """Perform seed clustering of coordinates.

    Labels points into clusters such that cluster contains exactly one copy
    of each dataset, then performs silhouettete analysis on the resulting
    clusters to determine the true number of clusters present, under the
    constraint that only equal-sized clusterings are valid, i.e. each
    dataset should appear an equal number of times in each cluster.

    See also:
      https://en.wikipedia.org/wiki/Silhouette_(clustering)
      http://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html

    Attributes:
      cluster_labels (scitbx.array_family.flex.int): A label for each coordinate.
    """

    def __init__(
        self, coordinates, n_datasets, n_sym_ops, min_silhouette_score, n_clusters=Auto
    ):
        """Initialise a seed_clustering object.

        Args:
          coordinates (scitbx.array_family.flex.double): The input array of coordinates
            on which to perform the analysis. The dimensions of the array should
            be (dim, `n_datasets` * `n_sym_ops`), where dim is the number of
            dimensions being used for the analysis.
          n_datasets (int): The number of datasets.
          n_sym_ops (int): The number of symmetry operations.
          min_silhouette_score (float): The minimum silhouette score to be used
            in automatic determination of the number of clusters.
          n_clusters (int): Optionally override the automatic determination of the
            number of clusters.
        """
        self.coords = coordinates

        self.cluster_labels = self._label_clusters_first_pass(n_datasets, n_sym_ops)

        if flex.max(self.cluster_labels) == 0:
            # assume single cluster
            return

        dist_mat, linkage_matrix = self._hierarchical_clustering()
        self.cluster_labels, _ = self._silhouette_analysis(
            self.cluster_labels,
            linkage_matrix,
            n_clusters=n_clusters,
            min_silhouette_score=min_silhouette_score,
        )

    def _label_clusters_first_pass(self, n_datasets, n_sym_ops):
        """First pass labelling of clusters.

        Labels points into clusters such that cluster contains exactly one copy
        of each dataset.

        Args:
          n_datasets (int): The number of datasets.
          n_sym_ops (int): The number of symmetry operations.

        Returns:
          cluster_labels (scitbx.array_family.flex.int): A label for each coordinate, labelled from
          0 .. n_sym_ops.
        """
        # initialise cluster labels: -1 signifies doesn't belong to a cluster
        cluster_labels = flex.int(self.coords.all()[0], -1)
        X_orig = self.coords.as_numpy_array()

        cluster_id = 0
        while cluster_labels.count(-1) > 0:
            dataset_ids = (
                flex.int_range(n_datasets * n_sym_ops) % n_datasets
            ).as_numpy_array()
            coord_ids = flex.int_range(dataset_ids.size).as_numpy_array()

            # select only those points that don't already belong to a cluster
            sel = np.where(cluster_labels == -1)
            X = X_orig[sel]
            dataset_ids = dataset_ids[sel]
            coord_ids = coord_ids[sel]

            # choose a high density point as seed for cluster
            nbrs = NearestNeighbors(
                n_neighbors=min(11, len(X)), algorithm="brute", metric="cosine"
            ).fit(X)
            distances, indices = nbrs.kneighbors(X)
            average_distance = flex.double([dist[1:].mean() for dist in distances])
            i = flex.min_index(average_distance)

            d_id = dataset_ids[i]
            cluster = np.array([coord_ids[i]])
            cluster_dataset_ids = np.array([d_id])
            xis = np.array([X[i]])

            for j in range(n_datasets - 1):
                # select only those rows that don't correspond to a dataset already
                # present in current cluster
                sel = np.where(dataset_ids != d_id)
                X = X[sel]
                dataset_ids = dataset_ids[sel]
                coord_ids = coord_ids[sel]

                assert len(X) > 0

                # Find nearest neighbour in cosine-space to the current cluster centroid
                nbrs = NearestNeighbors(
                    n_neighbors=min(1, len(X)), algorithm="brute", metric="cosine"
                ).fit(X)
                distances, indices = nbrs.kneighbors([xis.mean(axis=0)])
                k = indices[0][0]
                d_id = dataset_ids[k]
                cluster = np.append(cluster, coord_ids[k])
                cluster_dataset_ids = np.append(cluster_dataset_ids, d_id)
                xis = np.append(xis, [X[k]], axis=0)

            # label this cluster
            cluster_labels.set_selected(flex.size_t(cluster.tolist()), cluster_id)
            cluster_id += 1
        return cluster_labels

    def _hierarchical_clustering(self):
        """Perform hierarchical clustering on cluster centroids.

        Returns:
          Tuple[numpy.ndarray, numpy.ndarray]:
            A tuple containing
            the distance matrix as output by :func:`scipy.spatial.distance.pdist` and
            the linkage matrix as output by :func:`scipy.cluster.hierarchy.linkage`.
        """
        cluster_centroids = []
        X = self.coords.as_numpy_array()
        for i in set(self.cluster_labels):
            cluster_centroids.append(
                X[(self.cluster_labels == i).iselection().as_numpy_array()].mean(axis=0)
            )

        # hierarchical clustering of cluster centroids, using cosine metric
        dist_mat = ssd.pdist(cluster_centroids, metric="cosine")
        return dist_mat, hierarchy.linkage(dist_mat, method="average")

    def _silhouette_analysis(
        self, cluster_labels, linkage_matrix, n_clusters, min_silhouette_score
    ):
        """Compare valid equal-sized clustering using silhouette scores.

        Args:
          cluster_labels (scitbx.array_family.flex.int):
          linkage_matrix (numpy.ndarray): The hierarchical clustering of centroids of the
            initial clustering as produced by
            :func:`scipy.cluster.hierarchy.linkage`.
          n_clusters (int): Optionally override the automatic determination of the
            number of clusters.
          min_silhouette_score (float): The minimum silhouette score to be used
            in automatic determination of the number of clusters.

        Returns:
          cluster_labels (scitbx.array_family.flex.int): A label for each coordinate.
        """
        eps = 1e-6
        X = self.coords.as_numpy_array()

        cluster_labels_input = cluster_labels
        distances = linkage_matrix[::, 2]
        distances = np.insert(distances, 0, 0)
        silhouette_scores = flex.double()
        thresholds = flex.double()
        threshold_n_clusters = flex.size_t()
        for threshold in distances[1:]:
            cluster_labels = cluster_labels_input.deep_copy()
            labels = hierarchy.fcluster(
                linkage_matrix, threshold - eps, criterion="distance"
            ).tolist()
            counts = [labels.count(l) for l in set(labels)]
            if len(set(counts)) > 1:
                # only equal-sized clusters are valid
                continue

            n = len(set(labels))
            if n == 1:
                continue
            elif n_clusters is not Auto and n != n_clusters:
                continue
            for i in range(len(labels)):
                cluster_labels.set_selected(
                    cluster_labels_input == i, int(labels[i] - 1)
                )
            if len(set(cluster_labels)) == X.shape[0]:
                # silhouette coefficient not defined if 1 dataset per cluster
                # not sure what the default value should be
                sample_silhouette_values = np.full(cluster_labels.size(), 0)
            else:
                # Compute the silhouette scores for each sample
                sample_silhouette_values = metrics.silhouette_samples(
                    X, cluster_labels.as_numpy_array(), metric="cosine"
                )
            silhouette_avg = sample_silhouette_values.mean()
            silhouette_scores.append(silhouette_avg)
            thresholds.append(threshold)
            threshold_n_clusters.append(n)

            count_negative = (sample_silhouette_values < 0).sum()
            logger.info("Clustering:")
            logger.info("  Number of clusters: %i" % n)
            logger.info(
                "  Threshold score: %.3f (%.1f deg)"
                % (threshold, math.degrees(math.acos(1 - threshold)))
            )
            logger.info("  Silhouette score: %.3f" % silhouette_avg)
            logger.info(
                "  -ve silhouette scores: %.1f%%"
                % (100 * count_negative / sample_silhouette_values.size)
            )

        if n_clusters is Auto:
            idx = flex.max_index(silhouette_scores)
        else:
            idx = flex.first_index(threshold_n_clusters, n_clusters)
            if idx is None:
                raise Sorry("No valid clustering with %i clusters" % n_clusters)

        if n_clusters is Auto and silhouette_scores[idx] < min_silhouette_score:
            # assume single cluster
            cluster_labels = flex.int(cluster_labels.size(), 0)
        else:
            threshold = thresholds[idx] - eps
            labels = hierarchy.fcluster(linkage_matrix, threshold, criterion="distance")
            cluster_labels = flex.double(cluster_labels.size(), -1)
            for i in range(len(labels)):
                cluster_labels.set_selected(
                    cluster_labels_input == i, float(labels[i] - 1)
                )

        return cluster_labels, threshold

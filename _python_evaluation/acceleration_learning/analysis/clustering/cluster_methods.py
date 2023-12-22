"""Python module providing clustering methods"""
# Type hinting
from typing import List
from timeit import default_timer as timer
# Plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
# Color maps for plots
import matplotlib.cm as cm
import numpy as np

# Distance Time Wrapping
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw

# CLustering algorithms
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn_extra.cluster import KMedoids
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from tslearn.clustering import TimeSeriesKMeans
from tslearn.preprocessing import TimeSeriesScalerMeanVariance, TimeSeriesResampler

import pandas as pd

def cluster_with_kmeans(n_dimension_array):
    """ Cluster the data and display the found clusters
    
    Parameters:
    -- n_dimension_array: a numpy array of n dimensions
    """
    # Shape inputs
    reducer = PCA(n_components=2)
    reduced_data = reducer.fit_transform(n_dimension_array)
    for i in range(2, 5):
        kmeans = KMeans(init="k-means++", n_clusters=i, n_init=4)
        kmeans.fit(reduced_data)

        # Step size of the mesh. Decrease to increase the quality of the VQ.
        h = 0.02  # point in the mesh [x_min, x_max]x[y_min, y_max].

        # Plot the decision boundary. For that, we will assign a color to each
        x_min, x_max = reduced_data[:, 0].min() - 1, reduced_data[:, 0].max() + 1
        y_min, y_max = reduced_data[:, 1].min() - 1, reduced_data[:, 1].max() + 1
        xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

        # Obtain labels for each point in mesh. Use last trained model.
        Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

        # Put the result into a color plot
        Z = Z.reshape(xx.shape)
        plt.figure(i-1)
        plt.clf()
        plt.imshow(
            Z,
            interpolation="nearest",
            extent=(xx.min(), xx.max(), yy.min(), yy.max()),
            cmap=plt.cm.Paired,
            aspect="auto",
            origin="lower",
        )
        plt.plot(reduced_data[:, 0], reduced_data[:, 1], "k.", markersize=2)
        # Plot the centroids as a white X
        centroids = kmeans.cluster_centers_
        print("Centroids", centroids)
        plt.scatter(
            centroids[:, 0],
            centroids[:, 1],
            marker="x",
            s=169,
            linewidths=3,
            color="w",
            zorder=10,
        )      
        plt.title(
            "K-means clustering on the digits dataset (PCA-reduced data)\n"
            "Centroids are marked with white cross"
        )
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        plt.xticks(())
        plt.yticks(())
        plt.show()

def cluster_with_gaussian_mixture_model(n_dimension_array):
    X = n_dimension_array

    def get_initial_means(X, init_params, r):
        # Run a GaussianMixture with max_iter=0 to output the initalization means
        gmm = GaussianMixture(
            n_components=4, init_params=init_params, tol=1e-9, max_iter=0, random_state=r
        ).fit(X)
        return gmm.means_
    methods = ["kmeans", "random_from_data", "k-means++", "random"]
    colors = ["navy", "turquoise", "cornflowerblue", "darkorange"]
    times_init = {}
    relative_times = {}

    plt.figure(figsize=(4 * len(methods) // 2, 6))
    plt.subplots_adjust(
        bottom=0.1, top=0.9, hspace=0.15, wspace=0.05, left=0.05, right=0.95
    )

    for n, method in enumerate(methods):
        r = np.random.RandomState(seed=1234)
        plt.subplot(2, len(methods) // 2, n + 1)

        start = timer()
        ini = get_initial_means(X, method, r)
        end = timer()
        init_time = end - start

        gmm = GaussianMixture(
            n_components=4, means_init=ini, tol=1e-9, max_iter=2000, random_state=r
        ).fit(X)

        times_init[method] = init_time
        for i, color in enumerate(colors):
            data = X[gmm.predict(X) == i]
            plt.scatter(data[:, 0], data[:, 1], color=color, marker="x")

        plt.scatter(
            ini[:, 0], ini[:, 1], s=75, marker="D", c="orange", lw=1.5, edgecolors="black"
        )
        relative_times[method] = times_init[method] / times_init[methods[0]]

        plt.xticks(())
        plt.yticks(())
        plt.title(method, loc="left", fontsize=12)
        plt.title(
            "Iter %i | Init Time %.2fx" % (gmm.n_iter_, relative_times[method]),
            loc="right",
            fontsize=10,
        )
    plt.suptitle("GMM iterations and relative time taken to initialize")
    plt.show()

def distance_time_wrapping(first_array, second_array):
    """
    an array: [
        [x1,y1], [x2,y2], ...
    ]
    """
    mpl.rcParams['figure.dpi'] = 300
    savefig_options = dict(format="png", dpi=300, bbox_inches="tight")
    
    x1 = first_array
    x2 = second_array
    distance, warp_path = fastdtw(x1, x2, dist=euclidean)

    print("DTW distance", distance)

    fig, ax = plt.subplots(figsize=(16, 12))

    # Remove the border and axes ticks
    fig.patch.set_visible(False)
    ax.axis('off')

    for [map_x, map_y] in warp_path:
        ax.plot([map_x, map_y], [x1[map_x], x2[map_y]], '-k')

    ax.plot(x1, color='blue', marker='o', markersize=10, linewidth=5)
    ax.plot(x2, color='red', marker='o', markersize=10, linewidth=5)
    ax.tick_params(axis="both", which="major", labelsize=18)

    fig.savefig("ex2_dtw_distance.png", **savefig_options)

def cluster_with_time_series_kmeans(n_dimension_array, n_clusters=3):
    X_train = n_dimension_array
    size_of_output_series = 40

    seed = 0
    np.random.seed(seed)
    X_train = TimeSeriesScalerMeanVariance().fit_transform(X_train)
    # Make time series shorter
    X_train = TimeSeriesResampler(sz=size_of_output_series).fit_transform(X_train)
    sz = X_train.shape[1]
    km = TimeSeriesKMeans(n_clusters=n_clusters, verbose=True, random_state=seed)
    y_pred = km.fit_predict(X_train)
    
    plt.figure()
    for yi in range(n_clusters):
        plt.subplot(3, n_clusters, yi + 1)
        for xx in X_train[y_pred == yi]:
            plt.plot(xx.ravel(), "k-", alpha=.2)
        plt.plot(km.cluster_centers_[yi].ravel(), "r-")
        plt.xlim(0, sz)
        plt.ylim(-4, 4)
        plt.text(0.55, 0.85,'Cluster %d' % (yi + 1),
                transform=plt.gca().transAxes)
        if yi == 1:
            plt.title("Euclidean $k$-means")

    # DBA-k-means
    print("DBA k-means")
    dba_km = TimeSeriesKMeans(n_clusters=n_clusters,
                            n_init=2,
                            metric="dtw",
                            verbose=True,
                            max_iter_barycenter=10,
                            random_state=seed)
    y_pred = dba_km.fit_predict(X_train)

    for yi in range(n_clusters):
        plt.subplot(3, n_clusters, n_clusters + 1 + yi)
        for xx in X_train[y_pred == yi]:
            plt.plot(xx.ravel(), "k-", alpha=.2)
        plt.plot(dba_km.cluster_centers_[yi].ravel(), "r-")
        plt.xlim(0, sz)
        plt.ylim(-4, 4)
        plt.text(0.55, 0.85,'Cluster %d' % (yi + 1),
                transform=plt.gca().transAxes)
        if yi == 1:
            plt.title("DBA $k$-means")

    # Soft-DTW-k-means
    print("Soft-DTW k-means")
    sdtw_km = TimeSeriesKMeans(n_clusters=n_clusters,
                            metric="softdtw",
                            metric_params={"gamma": .01},
                            verbose=True,
                            random_state=seed)
    y_pred = sdtw_km.fit_predict(X_train)

    for yi in range(n_clusters):
        plt.subplot(3, n_clusters, n_clusters*2 + 1 + yi)
        for xx in X_train[y_pred == yi]:
            plt.plot(xx.ravel(), "k-", alpha=.2)
        plt.plot(sdtw_km.cluster_centers_[yi].ravel(), "r-")
        plt.xlim(0, sz)
        plt.ylim(-4, 4)
        plt.text(0.55, 0.85,'Cluster %d' % (yi + 1),
                transform=plt.gca().transAxes)
        if yi == 1:
            plt.title("Soft-DTW $k$-means")

    plt.tight_layout()
    plt.show()

def cluster_with_kmedoid(n_dimension_array, n_clusters=3):
    x = n_dimension_array
    scaler = MinMaxScaler().fit(x)
    x_scaled = x # scaler.transform(x)
    k_medoids = KMedoids(n_clusters=n_clusters, random_state=0)
    k_medoids.fit(x_scaled)
    y_kmed = k_medoids.fit_predict(x_scaled)
    
    rows = []
    for i, sample in enumerate(n_dimension_array):
        color = ""
        if y_kmed[i] == 0:
            color = "red"
        elif y_kmed[i] == 1:
            color = "blue"
        else:
            color = "green"
        rows.append(list(x_scaled[i]) + [color])
    df = pd.DataFrame(rows, columns=['Time', 'Overshoot', 'Acceleration', "Cluster"])
    # Creating figure
    fig = plt.figure(0)
    ax = plt.axes(projection ="3d")
    x = np.array(df['Time'])
    y = np.array(df['Overshoot'])
    z = np.array(df['Acceleration'])
    ax.scatter(x,y,z, marker=".", c=df["Cluster"], s = 40, cmap="RdBu", alpha=0.5)
    ax.scatter(k_medoids.cluster_centers_[:, 0], k_medoids.cluster_centers_[:,1], k_medoids.cluster_centers_[:,2], s = 80, c = 'black', label = 'Centroids', alpha=1)
    ax.set_xlabel('Times to peak')
    ax.set_ylabel('Overshoots')
    ax.set_zlabel('Max acceleration')
    
    # plt.scatter(x_scaled[y_kmed == 0, 0], x_scaled[y_kmed == 0, 1], s = 100, c = 'red', label = 'C1')
    # plt.scatter(x_scaled[y_kmed == 1, 0], x_scaled[y_kmed == 1, 1], s = 100, c = 'blue', label = 'C2')
    # plt.scatter(x_scaled[y_kmed == 2, 0], x_scaled[y_kmed == 2, 1], s = 100, c = 'green', label = 'C3')
    # plt.scatter(k_medoids.cluster_centers_[:, 0], k_medoids.cluster_centers_[:,1], s = 100, c = 'yellow', label = 'Centroids')
    plt.legend()
    plt.savefig(r'C:\git\kdp_hlb_evalframework\_temp\analysis\clustering\cluster.png')

    medoids = list(zip(k_medoids.cluster_centers_[:, 0], k_medoids.cluster_centers_[:, 1], k_medoids.cluster_centers_[:, 2]))
    return medoids

def kmean_cluster(database_df: pd.DataFrame, n_clusters: int, database_parameters: List[str]):
    """Find the group of data"""
    data = database_df[database_parameters]
    model = KMeans(n_clusters=n_clusters).fit(data)
    centers = model.cluster_centers_
    predicted = model.predict(data)
    first_col = database_parameters[0]
    second_col = database_parameters[1]
    plt.title(f"K-means clusters of {first_col} parameter and {second_col} parameter")
    plt.scatter(data[first_col], data[second_col], c=predicted, cmap=cm.rainbow)
    plt.scatter(
        centers[:, data.columns.get_loc(first_col)],
        centers[:, data.columns.get_loc(second_col)],
        c='black',
        s=200,
        alpha=0.5
    )
    plt.xlabel(first_col)
    plt.ylabel(second_col)
    plt.show()

def make_ellipses(gmm: GaussianMixture, ax, n_groups: int):
    """Draw ellipses around the gaussian mixture clusters"""
    colors = ["cyan", "orange", "yellow"]
    for n, color in zip(range(n_groups), colors):
        if gmm.covariance_type == "full":
            covariances = gmm.covariances_[n][:2, :2]
        elif gmm.covariance_type == "tied":
            covariances = gmm.covariances_[:2, :2]
        elif gmm.covariance_type == "diag":
            covariances = np.diag(gmm.covariances_[n][:2])
        elif gmm.covariance_type == "spherical":
            covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_[n]
        v, w = np.linalg.eigh(covariances)
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        ell = mpl.patches.Ellipse(
            gmm.means_[n, :2], v[0], v[1], angle=180 + angle, color=color
        )
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ax.add_artist(ell)
        ax.set_aspect("equal", "datalim")

def gaussian_mixture_cluster(database_df: pd.DataFrame, n_clusters: int, database_parameters: List[str]):
    """Find the group of data with a Gaussian mixture clustering method"""
    data = database_df[database_parameters]
    covariance_type = ["spherical", "diag", "tied", "full"]
    estimator = GaussianMixture(
        n_components=n_clusters, covariance_type="full", max_iter=20, random_state=0
    )
    estimator.fit(data)
    h = plt.subplot(1, 1, 1)
    make_ellipses(estimator, h, n_clusters)
    predicted = estimator.predict(data)
    centers = estimator.means_
    first_col = database_parameters[0]
    second_col = database_parameters[1]
    plt.title(f"Gaussian mixture clusters of {first_col} parameter and {second_col} parameter")
    plt.scatter(data[first_col], data[second_col], c=predicted, cmap=cm.rainbow)
    plt.scatter(
        centers[:, data.columns.get_loc(first_col)],
        centers[:, data.columns.get_loc(second_col)],
        c='black',
        s=50,
        alpha=0.5
    )
    plt.xlabel(first_col)
    plt.ylabel(second_col)
    plt.show()

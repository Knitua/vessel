from types import SimpleNamespace

import numpy as np

import svv.forest.connect.assign as assign_mod


def _tree(points):
    data = np.zeros((len(points), 31), dtype=float)
    data[:, 3:6] = points
    data[:, 15:17] = np.nan
    return SimpleNamespace(data=data)


def test_assign_network_dense_fallback_adds_missing_connection_functions(monkeypatch):
    tree_0_points = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    tree_1_points = np.array([[0.0, 1.0, 0.0], [10.0, 1.0, 0.0]])
    forest = SimpleNamespace(
        networks=[[_tree(tree_0_points), _tree(tree_1_points)]],
        n_trees_per_network=[2],
        convex=True,
    )

    def sparse_matching_fails(_cost):
        raise ValueError("no full matching")

    monkeypatch.setattr(assign_mod, "min_weight_full_bipartite_matching", sparse_matching_fails)
    monkeypatch.setattr(
        assign_mod,
        "linear_sum_assignment",
        lambda _cost: (np.array([0, 1]), np.array([1, 0])),
    )

    assignments, connections = assign_mod.assign_network(forest, neighbors=1)

    assert assignments == [[0, 1], [1, 0]]
    assert len(connections) == 1
    assert len(connections[0]) == 2
    np.testing.assert_allclose(connections[0][0](np.array([0.0, 1.0])), [tree_0_points[0], tree_1_points[1]])
    np.testing.assert_allclose(connections[0][1](np.array([0.0, 1.0])), [tree_0_points[1], tree_1_points[0]])

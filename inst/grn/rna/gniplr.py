"""GNIPLR inference and network export."""
from pathlib import Path

import pandas as pd


def gniplr_sor_g(gene_1, gene_2):
    import numpy as np

    com = []
    for j in range(len(gene_1)):
        com.append((gene_1[j], gene_2[j]))
    sort_1 = sorted(com, key=lambda x: x[0])
    a1 = np.array([x[0] for x in sort_1], dtype=float).reshape(-1, 1)
    b1 = np.array([x[1] for x in sort_1], dtype=float).reshape(-1, 1)
    return a1, b1


def gniplr_lasso_regression(degree, alpha):
    from sklearn.linear_model import Lasso
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import PolynomialFeatures, StandardScaler

    return Pipeline([
        ("poly", PolynomialFeatures(degree=degree)),
        ("std_scaler", StandardScaler()),
        ("lasso_reg", Lasso(alpha=alpha))])


def gniplr_lasso_regress(gene_1, gene_2, degree, alpha):
    import numpy as np

    le_g1 = len(gene_1)
    t = np.linspace(min(gene_1), max(gene_1), num=le_g1)
    a1, b1 = gniplr_sor_g(gene_1, gene_2)
    lasso_reg = gniplr_lasso_regression(degree, alpha)
    lasso_reg.fit(a1, t)
    gene2 = lasso_reg.predict(b1)
    gene1 = lasso_reg.predict(a1)
    return gene1, gene2


def gniplr_granger(gene_1, gene_2, n):
    """GNIPLR granger function with explicit sample size n."""
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import f
    from sklearn.linear_model import LinearRegression

    com = [(gene_1[j], gene_2[j]) for j in range(n)]
    sort_1 = sorted(com, key=lambda x: x[0])
    a1 = [x[0] for x in sort_1]
    b1 = [x[1] for x in sort_1]
    a1_cha = [300 * a1[k] / 90001 for k in range(n)]
    b1_cha = [300 * b1[k] / 90001 for k in range(n)]
    model = LinearRegression()


    a1_t1_1 = np.array(a1_cha[1:n])
    b1_t_1 = np.array(b1_cha[0:n - 1])
    b1_t1_1 = np.array(b1_cha[1:n])
    b1_t1_t2_1 = sm.add_constant(b1_t1_1)
    a1_t1_t2_b1_1 = sm.add_constant(np.vstack((a1_t1_1, b1_t1_1)).T)
    model.fit(a1_t1_t2_b1_1, b1_t_1)
    b1_rssu_1 = np.sum((model.predict(a1_t1_t2_b1_1) - b1_t_1) ** 2)
    model.fit(b1_t1_t2_1, b1_t_1)
    b1_rssr_1 = np.sum((model.predict(b1_t1_t2_1) - b1_t_1) ** 2)
    f1_1 = ((b1_rssr_1 - b1_rssu_1) / 1) / (b1_rssu_1 / (n - 3))
    p1_1 = f.sf(f1_1, 1, n - 3)


    a1_t1 = np.array(a1_cha[1:n - 1])
    a1_t2 = np.array(a1_cha[2:n])
    b1_t = np.array(b1_cha[0:n - 2])
    b1_t1 = np.array(b1_cha[1:n - 1])
    b1_t2 = np.array(b1_cha[2:n])
    b1_t1_t2 = sm.add_constant(np.vstack((b1_t1, b1_t2)).T)
    a1_t1_t2_b1 = sm.add_constant(np.vstack((a1_t1, a1_t2, b1_t1, b1_t2)).T)
    model.fit(a1_t1_t2_b1, b1_t)
    b1_rssu_2 = np.sum((model.predict(a1_t1_t2_b1) - b1_t) ** 2)
    model.fit(b1_t1_t2, b1_t)
    b1_rssr_2 = np.sum((model.predict(b1_t1_t2) - b1_t) ** 2)
    f1_2 = ((b1_rssr_2 - b1_rssu_2) / 2) / (b1_rssu_2 / (n - 2 - 2 - 1))
    p1_2 = f.sf(f1_2, 2, n - 2 - 2 - 1)


    a1_t1_3 = np.array(a1_cha[1:n - 2])
    a1_t1_t3_3 = np.array(a1_cha[2:n - 1])
    a1_t2_3 = np.array(a1_cha[3:n])
    b1_t_3 = np.array(b1_cha[0:n - 3])
    b1_t1_3 = np.array(b1_cha[1:n - 2])
    b1_t1_t3_3 = np.array(b1_cha[2:n - 1])
    b1_t2_3 = np.array(b1_cha[3:n])
    b1_t1_t2_3 = sm.add_constant(np.vstack((b1_t1_3, b1_t2_3, b1_t1_t3_3)).T)
    a1_t1_t2_b1_3 = sm.add_constant(
        np.vstack((a1_t1_3, a1_t2_3, a1_t1_t3_3, b1_t1_3, b1_t2_3, b1_t1_t3_3)).T
    )
    model.fit(a1_t1_t2_b1_3, b1_t_3)
    b1_rssu_3 = np.sum((model.predict(a1_t1_t2_b1_3) - b1_t_3) ** 2)
    model.fit(b1_t1_t2_3, b1_t_3)
    b1_rssr_3 = np.sum((model.predict(b1_t1_t2_3) - b1_t_3) ** 2)
    f1_3 = ((b1_rssr_3 - b1_rssu_3) / 3) / (b1_rssu_3 / (n - 7))
    p1_3 = f.sf(f1_3, 3, n - 7)

    return min([p1_1, p1_2, p1_3])


def infer_gniplr(expression, out_file, lasso_degree=30, lasso_alpha=0.04):
    import numpy as np

    genes = list(expression.index)
    vectors = {
        gene: expression.loc[gene].to_numpy(dtype=float)
        for gene in genes
    }
    rows = []
    for gene_a in genes:
        for gene_b in genes:
            if gene_a == gene_b:
                continue
            gene1, gene2 = gniplr_lasso_regress(
                vectors[gene_a], vectors[gene_b], lasso_degree, lasso_alpha
            )
            p_min = gniplr_granger(gene1, gene2, len(vectors[gene_a]))
            weight = -np.log10(max(p_min, 1e-300))
            rows.append((gene_a, gene_b, weight))
    table = pd.DataFrame(rows, columns=["regulator", "target", "weight"])
    table = table.sort_values("weight", ascending=False)
    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_file, index=False)
